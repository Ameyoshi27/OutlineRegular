#include "MaskVectorizer.h"

#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

struct DatasetCloser {
    void operator()(GDALDataset* dataset) const
    {
        if (dataset) GDALClose(dataset);
    }
};

using DatasetPtr = std::unique_ptr<GDALDataset, DatasetCloser>;

bool RemoveShapefileFamily(const std::filesystem::path& shpPath)
{
    static const char* const extensions[] = {
        ".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix",
        ".sbn", ".sbx", ".fbn", ".fbx", ".ain", ".aih", ".atx"
    };
    const std::filesystem::path stem = shpPath.parent_path() / shpPath.stem();
    for (const char* extension : extensions) {
        std::filesystem::path file = stem;
        file += extension;
        std::error_code ec;
        if (!std::filesystem::remove(file, ec) && ec) {
            std::cerr << "[SHP] Cannot remove " << file.string()
                      << ": " << ec.message() << std::endl;
            return false;
        }
    }
    return true;
}

struct PolygonRecord {
    GIntBig fid = OGRNullFID;
    std::unique_ptr<OGRGeometry> geometry;
    std::unique_ptr<OGRPolygon> exteriorShell;
    OGREnvelope envelope = {};
    double area = 0.0;
    bool remove = false;
    std::size_t containerIndex = std::numeric_limits<std::size_t>::max();
};

long long RemoveContainedPolygons(OGRLayer* layer)
{
    if (!layer) return 0;
    std::vector<PolygonRecord> polygons;
    polygons.reserve(static_cast<std::size_t>(
        std::max<GIntBig>(0, layer->GetFeatureCount(FALSE))));

    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        OGRGeometry* sourceGeometry = feature->GetGeometryRef();
        if (sourceGeometry && !sourceGeometry->IsEmpty()) {
            PolygonRecord record;
            record.fid = feature->GetFID();
            record.geometry.reset(sourceGeometry->clone());
            record.geometry->getEnvelope(&record.envelope);
            if (const auto* surface = dynamic_cast<const OGRSurface*>(record.geometry.get())) {
                record.area = surface->get_Area();
            }
            if (const auto* polygon = dynamic_cast<const OGRPolygon*>(record.geometry.get())) {
                if (const OGRLinearRing* exterior = polygon->getExteriorRing()) {
                    record.exteriorShell = std::make_unique<OGRPolygon>();
                    record.exteriorShell->addRingDirectly(new OGRLinearRing(*exterior));
                }
            }
            if (record.area > 0.0) polygons.push_back(std::move(record));
        }
        OGRFeature::DestroyFeature(feature);
    }

    std::sort(polygons.begin(), polygons.end(),
              [](const PolygonRecord& a, const PolygonRecord& b) {
                  return a.area > b.area;
              });
    for (std::size_t innerIndex = 1; innerIndex < polygons.size(); ++innerIndex) {
        PolygonRecord& inner = polygons[innerIndex];
        // 从大到小遍历，使嵌套结构归属到最近的包含壳。
        // (先遇到的外层即最近容器)。
        for (std::size_t outerIndex = innerIndex; outerIndex-- > 0;) {
            const PolygonRecord& outer = polygons[outerIndex];
            if (!outer.exteriorShell) continue;
            if (inner.area >= outer.area * (1.0 - 1e-9)) continue;
            if (outer.envelope.MinX > inner.envelope.MinX ||
                outer.envelope.MinY > inner.envelope.MinY ||
                outer.envelope.MaxX < inner.envelope.MaxX ||
                outer.envelope.MaxY < inner.envelope.MaxY) {
                continue;
            }
            // 多边形化的标签孤岛既是父的洞又是子多边形；
            // 对外壳判定，父多边形对应的洞
            // 不会掩盖包含关系。
            if (outer.exteriorShell->Contains(inner.geometry.get())) {
                inner.containerIndex = outerIndex;
                break;
            }
        }
    }

    // 决定删除方向。多批次叠加的掩膜可能
    // 用粗外层实例包住另一批次的精细实例：
    // 外环只是内层并集微扩一圈的重复伪影
    // (实测：内层覆盖外层 90-100%)。
    // 外壳应丢弃、保留精细内层。
    // 真实包含(院落建筑等)内层覆盖率低，
    // 维持原行为(删内层)。
    std::vector<double> innerAreaSum(polygons.size(), 0.0);
    std::vector<int> innerCount(polygons.size(), 0);
    for (std::size_t i = 0; i < polygons.size(); ++i) {
        if (polygons[i].containerIndex == std::numeric_limits<std::size_t>::max()) {
            continue;
        }
        innerAreaSum[polygons[i].containerIndex] += polygons[i].area;
        ++innerCount[polygons[i].containerIndex];
    }
    std::vector<bool> isWrapper(polygons.size(), false);
    for (std::size_t i = 0; i < polygons.size(); ++i) {
        if (innerCount[i] == 0) continue;
        const double shellArea = std::max(polygons[i].area, 1e-9);
        const double coverage = innerAreaSum[i] / shellArea;
        const bool wrapper = (innerCount[i] >= 2 && coverage >= 0.60) ||
                             (innerCount[i] == 1 && coverage >= 0.75);
        if (wrapper) {
            isWrapper[i] = true;
            std::cerr << "[Mask] wrapper artifact: outer_area=" << polygons[i].area
                      << " inner_count=" << innerCount[i]
                      << " inner_coverage=" << coverage
                      << " -> drop outer, keep inner" << std::endl;
        }
    }
    for (std::size_t i = 0; i < polygons.size(); ++i) {
        if (polygons[i].containerIndex == std::numeric_limits<std::size_t>::max()) {
            continue;
        }
        polygons[i].remove = !isWrapper[polygons[i].containerIndex];
    }
    for (std::size_t i = 0; i < polygons.size(); ++i) {
        if (isWrapper[i]) polygons[i].remove = true;
    }

    long long removed = 0;
    for (std::size_t outerIndex = 0; outerIndex < polygons.size(); ++outerIndex) {
        PolygonRecord& outer = polygons[outerIndex];
        const auto* polygon = dynamic_cast<const OGRPolygon*>(outer.geometry.get());
        if (!polygon || polygon->getNumInteriorRings() == 0) continue;

        std::vector<const OGRGeometry*> removedChildren;
        for (const PolygonRecord& child : polygons) {
            if (child.remove && child.containerIndex == outerIndex) {
                removedChildren.push_back(child.geometry.get());
            }
        }
        if (removedChildren.empty()) continue;

        OGRPolygon rebuilt;
        rebuilt.assignSpatialReference(polygon->getSpatialReference());
        rebuilt.addRingDirectly(new OGRLinearRing(*polygon->getExteriorRing()));
        for (int ringIndex = 0; ringIndex < polygon->getNumInteriorRings(); ++ringIndex) {
            const OGRLinearRing* ring = polygon->getInteriorRing(ringIndex);
            OGRPolygon hole;
            hole.addRingDirectly(new OGRLinearRing(*ring));
            bool belongsToRemovedChild = false;
            for (const OGRGeometry* child : removedChildren) {
                if (hole.Intersects(child)) {
                    belongsToRemovedChild = true;
                    break;
                }
            }
            if (!belongsToRemovedChild) {
                rebuilt.addRingDirectly(new OGRLinearRing(*ring));
            }
        }

        OGRFeature* feature = layer->GetFeature(outer.fid);
        if (feature) {
            feature->SetGeometry(&rebuilt);
            layer->SetFeature(feature);
            OGRFeature::DestroyFeature(feature);
        }
    }

    for (const PolygonRecord& polygon : polygons) {
        if (polygon.remove && polygon.fid != OGRNullFID &&
            layer->DeleteFeature(polygon.fid) == OGRERR_NONE) {
            ++removed;
        }
    }
    if (removed > 0) layer->SyncToDisk();
    return removed;
}

// 刻意用米表示。像素半径由 GeoTIFF 仿射变换换算，
// 使同一设置适配不同分辨率影像。
constexpr double kSeedErosionDistanceMeters = 1.2;
constexpr double kMinSeedAreaSquareMeters = 5.0;
// 普通 RGB 掩膜保留实例颜色。
// (如渲染叠加图)才走二值兜底。
constexpr std::size_t kMaxReliableColorLabels = 65536;
// 近黑背景阈值：三通道均 <= 此值的像素判为背景。
// 生产掩膜实测：裁剪瓦片含 ~140 万个 RGB=(1,1,1) 像素，视觉上是黑色
// 背景但通过了旧的"任一通道非零"前景测试，形成包裹全部建筑的巨型
// 连通域。真实实例色标都是饱和色(be1401/e91501 等)，不会是近黑。
constexpr int kMaskBackgroundMaxChannel = 2;
constexpr double kNarrowWaistMaxWidthRatio = 0.38;
constexpr double kNarrowWaistMinCoreAreaRatio = 0.10;
constexpr double kNarrowWaistMinSplitGapPixels = 2.0;
// 批次叠加外壳阈值(像素级)。
// 多批次叠加可能用粗外层包住另一批
// 的精细实例；外环只是内层并集的
// 微扩重复(生产掩膜实测：外壳的
// is 90-100% while real containment (courtyard buildings) stays below 10%.
constexpr double kWrapperMultiInnerCoverage = 0.60;   // >= 2 inner components
constexpr double kWrapperSingleInnerCoverage = 0.75;  // single inner component
constexpr double kWrapperContainerVoteFrac = 0.95;    // ring votes must be one-sided
constexpr double kWrapperMinComponentArea = 4.0;      // m2, ignore specks on both sides
// 发丝接缝/桥接分量：多批次叠加在实例间
// 留下 1-2px 抗锯齿缝(实测颜色
// d50400/f70e00 族)。发丝在两实例之间时
// 会被分水岭吸收后焊成连通域；
// 横跨桥接时直接连通。两者都破坏
// 真实实例分隔(案例 id=1949/id=3400)。
// 发丝置为背景使实例断开。
// 真实建筑不可能是发丝状，形状门槛零误伤。
constexpr int kHairlineMaxThicknessPx = 3;   // bbox thin side
constexpr int kHairlineMinLengthPx = 15;     // bbox long side
constexpr double kHairlineFillRatio = 0.55;  // pixels / (thin x long bbox): line-like

long long RemoveHairlineSeamComponents(cv::Mat& colorImage, double pixelAreaMeters)
{
    const int w = colorImage.cols;
    const int h = colorImage.rows;
    std::vector<int> comp(static_cast<std::size_t>(w) * h, -1);
    struct HairInfo {
        int minX = 0, minY = 0, maxX = 0, maxY = 0;
        std::size_t pixels = 0;
    };
    std::vector<HairInfo> comps;
    std::vector<int> stack;
    for (int y = 0; y < h; ++y) {
        const int* row = colorImage.ptr<int>(y);
        for (int x = 0; x < w; ++x) {
            const std::size_t idx = static_cast<std::size_t>(y) * w + x;
            if (row[x] == 0 || comp[idx] >= 0) continue;
            const int id = static_cast<int>(comps.size());
            HairInfo info;
            info.minX = x; info.maxX = x; info.minY = y; info.maxY = y;
            stack.clear();
            stack.push_back(static_cast<int>(idx));
            comp[idx] = id;
            const int color = row[x];
            while (!stack.empty()) {
                const int cur = stack.back();
                stack.pop_back();
                const int cy = cur / w;
                const int cx = cur % w;
                info.pixels++;
                info.minX = std::min(info.minX, cx);
                info.maxX = std::max(info.maxX, cx);
                info.minY = std::min(info.minY, cy);
                info.maxY = std::max(info.maxY, cy);
                const int nb[4] = {cur - 1, cur + 1, cur - w, cur + w};
                const bool ok[4] = {cx > 0, cx < w - 1, cy > 0, cy < h - 1};
                for (int k = 0; k < 4; ++k) {
                    if (!ok[k]) continue;
                    const int n = nb[k];
                    if (comp[n] < 0 && colorImage.ptr<int>(n / w)[n % w] == color) {
                        comp[n] = id;
                        stack.push_back(n);
                    }
                }
            }
            comps.push_back(info);
        }
    }

    std::vector<bool> isHairline(comps.size(), false);
    long long removed = 0;
    for (std::size_t i = 0; i < comps.size(); ++i) {
        const int bw = comps[i].maxX - comps[i].minX + 1;
        const int bh = comps[i].maxY - comps[i].minY + 1;
        const int thin = std::min(bw, bh);
        const int longSide = std::max(bw, bh);
        if (thin > kHairlineMaxThicknessPx) continue;
        if (longSide < kHairlineMinLengthPx) continue;
        const double fill = static_cast<double>(comps[i].pixels) /
                            static_cast<double>(thin * longSide);
        if (fill < kHairlineFillRatio) continue;
        isHairline[i] = true;
        ++removed;
        if (removed <= 40) {
            std::cerr << "[Mask] hairline seam removed: color_id px=" << comps[i].pixels
                      << " bbox=" << bw << "x" << bh
                      << " at px[" << comps[i].minX << "," << comps[i].minY << "]"
                      << " (" << comps[i].pixels * pixelAreaMeters << " m2)" << std::endl;
        }
    }
    if (removed > 0) {
        for (int y = 0; y < h; ++y) {
            int* row = colorImage.ptr<int>(y);
            for (int x = 0; x < w; ++x) {
                const int c = comp[static_cast<std::size_t>(y) * w + x];
                if (c >= 0 && isHairline[c]) row[x] = 0;
            }
        }
    }
    return removed;
}

// 在分水岭标注前于像素级剥除外壳：
// 多数内层实例小于种子面积阈值，
// 会被外壳标签吸收，多边形级
// 修复根本看不到它们。
long long RemoveBatchWrapperRings(cv::Mat& colorImage, double pixelAreaMeters)
{
    const int w = colorImage.cols;
    const int h = colorImage.rows;
    std::vector<int> comp(static_cast<std::size_t>(w) * h, -1);

    struct CompInfo {
        std::size_t pixels = 0;
        int minX = 0, minY = 0, maxX = 0, maxY = 0;
        int container = -1;
        double containerVoteFrac = 0.0;
    };
    std::vector<CompInfo> comps;
    std::vector<int> stack;
    for (int y = 0; y < h; ++y) {
        const int* row = colorImage.ptr<int>(y);
        for (int x = 0; x < w; ++x) {
            const std::size_t idx = static_cast<std::size_t>(y) * w + x;
            if (row[x] == 0 || comp[idx] >= 0) continue;
            const int id = static_cast<int>(comps.size());
            CompInfo info;
            info.minX = x; info.maxX = x; info.minY = y; info.maxY = y;
            stack.clear();
            stack.push_back(static_cast<int>(idx));
            comp[idx] = id;
            const int color = row[x];
            while (!stack.empty()) {
                const int cur = stack.back();
                stack.pop_back();
                const int cy = cur / w;
                const int cx = cur % w;
                info.pixels++;
                info.minX = std::min(info.minX, cx);
                info.maxX = std::max(info.maxX, cx);
                info.minY = std::min(info.minY, cy);
                info.maxY = std::max(info.maxY, cy);
                const int nb[4] = {cur - 1, cur + 1, cur - w, cur + w};
                const bool ok[4] = {cx > 0, cx < w - 1, cy > 0, cy < h - 1};
                for (int k = 0; k < 4; ++k) {
                    if (!ok[k]) continue;
                    const int n = nb[k];
                    if (comp[n] < 0 &&
                        colorImage.ptr<int>(n / w)[n % w] == color) {
                        comp[n] = id;
                        stack.push_back(n);
                    }
                }
            }
            comps.push_back(info);
        }
    }

    // 不同连通域之间的边界投票(双向)。
    std::vector<std::unordered_map<int, std::size_t>> votes(comps.size());
    for (int y = 0; y < h; ++y) {
        const int* row = colorImage.ptr<int>(y);
        for (int x = 0; x < w; ++x) {
            const std::size_t idx = static_cast<std::size_t>(y) * w + x;
            const int a = comp[idx];
            if (a < 0) continue;   // background never needs a container
            if (x + 1 < w) {
                const int b = comp[idx + 1];
                if (b >= 0 && b != a) { ++votes[a][b]; ++votes[b][a]; }
            }
            if (y + 1 < h) {
                const int b = comp[idx + w];
                if (b >= 0 && b != a) { ++votes[a][b]; ++votes[b][a]; }
            }
        }
    }
    for (std::size_t i = 0; i < comps.size(); ++i) {
        // 容器 = 真正包围当前分量的最小连通域：
        // 必须接触(共享边)且包围盒包含我们
        // (方向校验——投票是对称的)。紧贴的
        // 内层实例共享边界，单一投票比例
        // 门槛无效；接触+包含才有效。
        int best = -1;
        std::size_t bestPixels = 0;
        for (const auto& v : votes[i]) {
            const std::size_t j = static_cast<std::size_t>(v.first);
            if (v.second < 5) continue;                       // needs real contact
            if (comps[j].pixels <= comps[i].pixels) continue;  // must be bigger
            if (comps[j].minX > comps[i].minX || comps[j].minY > comps[i].minY ||
                comps[j].maxX < comps[i].maxX || comps[j].maxY < comps[i].maxY) {
                continue;                                     // must contain us
            }
            if (best < 0 || comps[j].pixels < bestPixels) {
                best = v.first;
                bestPixels = comps[j].pixels;
            }
        }
        if (best >= 0) {
            comps[i].container = best;
            comps[i].containerVoteFrac = 1.0;
        }
    }

    // 按容器判定外壳(双侧最小面积门槛
    // 过滤碎屑：那是标签噪声不是批次重复)。
    const std::size_t minPixels = static_cast<std::size_t>(
        std::max(1.0, kWrapperMinComponentArea / std::max(pixelAreaMeters, 1e-9)));
    std::vector<std::size_t> innerPixels(comps.size(), 0);
    std::vector<int> innerCount(comps.size(), 0);
    for (std::size_t i = 0; i < comps.size(); ++i) {
        if (comps[i].container < 0) continue;
        if (comps[i].pixels < minPixels) continue;
        if (comps[static_cast<std::size_t>(comps[i].container)].pixels < minPixels) continue;
        innerPixels[comps[i].container] += comps[i].pixels;
        ++innerCount[comps[i].container];
    }
    std::vector<bool> isWrapper(comps.size(), false);
    long long wrappers = 0;
    for (std::size_t i = 0; i < comps.size(); ++i) {
        if (innerCount[i] == 0) continue;
        const double coverage =
            static_cast<double>(innerPixels[i]) / static_cast<double>(comps[i].pixels);
        const bool coveragePass = (innerCount[i] >= 2 && coverage >= kWrapperMultiInnerCoverage) ||
                                  (innerCount[i] == 1 && coverage >= kWrapperSingleInnerCoverage);
        if (!coveragePass) continue;

        // 环厚度剖面：批次重复伪影外环是
        // 处处 1-3px 薄边；而画幅裁断残留
        // (单栋真楼，细批次被瓦片裁断，
        // 粗批次剩余成为环)一侧有厚实体——
        // 此类必须交给合并阶段
        // 而不是在此删除。
        const bool thinRing = [&]() {
            std::vector<int> innerSet;
            for (std::size_t k = 0; k < comps.size(); ++k) {
                if (comps[k].container == static_cast<int>(i) &&
                    comps[k].pixels >= minPixels) {
                    innerSet.push_back(static_cast<int>(k));
                }
            }
            const int margin = 8;
            const int wx0 = std::max(0, comps[i].minX - margin);
            const int wy0 = std::max(0, comps[i].minY - margin);
            const int wx1 = std::min(w - 1, comps[i].maxX + margin);
            const int wy1 = std::min(h - 1, comps[i].maxY + margin);
            const int ww = wx1 - wx0 + 1;
            const int wh = wy1 - wy0 + 1;
            std::vector<std::int16_t> dist(static_cast<std::size_t>(ww) * wh, -1);
            std::vector<int> queue;
            for (int y = wy0; y <= wy1; ++y) {
                for (int x = wx0; x <= wx1; ++x) {
                    const int c = comp[static_cast<std::size_t>(y) * w + x];
                    bool isInner = false;
                    for (int innerId : innerSet) {
                        if (c == innerId) { isInner = true; break; }
                    }
                    if (isInner) {
                        const std::size_t wi = static_cast<std::size_t>(y - wy0) * ww + (x - wx0);
                        dist[wi] = 0;
                        queue.push_back(static_cast<int>(wi));
                    }
                }
            }
            // BFS from the inner union, capped: everything farther than the
            // cap stays "far" and counts against the thin-ring test.
            const std::int16_t cap = 6;
            for (std::size_t head = 0; head < queue.size(); ++head) {
                const int wi = queue[head];
                const std::int16_t d = dist[static_cast<std::size_t>(wi)];
                if (d >= cap) continue;
                const int cx = wi % ww;
                const int cy = wi / ww;
                const int nb[4] = {wi - 1, wi + 1, wi - ww, wi + ww};
                const bool ok[4] = {cx > 0, cx < ww - 1, cy > 0, cy < wh - 1};
                for (int k = 0; k < 4; ++k) {
                    if (!ok[k]) continue;
                    const std::size_t n = static_cast<std::size_t>(nb[k]);
                    if (dist[n] >= 0) continue;
                    dist[n] = static_cast<std::int16_t>(d + 1);
                    queue.push_back(static_cast<int>(n));
                }
            }
            std::vector<int> outerDists;
            for (int y = comps[i].minY; y <= comps[i].maxY; ++y) {
                for (int x = comps[i].minX; x <= comps[i].maxX; ++x) {
                    if (comp[static_cast<std::size_t>(y) * w + x] != static_cast<int>(i)) continue;
                    outerDists.push_back(dist[static_cast<std::size_t>(y - wy0) * ww + (x - wx0)]);
                }
            }
            if (outerDists.empty()) return false;
            std::sort(outerDists.begin(), outerDists.end());
            const int p90 = outerDists[outerDists.size() / 10 * 9];
            const int p98 = outerDists[outerDists.size() - 1 - outerDists.size() / 50];
            return p90 <= 3 && p98 <= cap;
        }();

        if (!thinRing) {
            std::cerr << "[Mask] contained-but-thick ring kept for merge: outer="
                      << comps[i].pixels * pixelAreaMeters
                      << " m2 at px[" << comps[i].minX << "," << comps[i].minY << "]"
                      << ", inner_coverage=" << coverage
                      << " (likely batch-edge remainder of one real building)" << std::endl;
            continue;
        }
        isWrapper[i] = true;
        ++wrappers;
        std::cerr << "[Mask] wrapper ring removed: outer=" << comps[i].pixels * pixelAreaMeters
                  << " m2 at px[" << comps[i].minX << "," << comps[i].minY << " "
                  << comps[i].maxX - comps[i].minX + 1 << "x" << comps[i].maxY - comps[i].minY + 1
                  << "], inner_count=" << innerCount[i]
                  << ", inner_coverage=" << coverage << std::endl;
    }

    if (wrappers > 0) {
        for (int y = 0; y < h; ++y) {
            int* row = colorImage.ptr<int>(y);
            for (int x = 0; x < w; ++x) {
                const int c = comp[static_cast<std::size_t>(y) * w + x];
                if (c >= 0 && isWrapper[c]) row[x] = 0;
            }
        }
    }
    return wrappers;
}

struct PixelBounds {
    int minX = std::numeric_limits<int>::max();
    int minY = std::numeric_limits<int>::max();
    int maxX = -1;
    int maxY = -1;

    void include(int x, int y)
    {
        minX = std::min(minX, x);
        minY = std::min(minY, y);
        maxX = std::max(maxX, x);
        maxY = std::max(maxY, y);
    }

    cv::Rect rect() const
    {
        return cv::Rect(minX, minY, maxX - minX + 1, maxY - minY + 1);
    }
};

struct SplitCandidate {
    int label = 0;
    int seedIndex = 0;
    int seedRow = -1;
    int seedCol = -1;
    double coreArea = 0.0;
    double coreWidth = 0.0;
};

std::vector<int> SelectColorBands(GDALDataset* dataset)
{
    std::vector<int> selected;
    const int bandCount = dataset ? dataset->GetRasterCount() : 0;
    for (int i = 1; i <= bandCount; ++i) {
        GDALRasterBand* band = dataset->GetRasterBand(i);
        if (!band) continue;
        const GDALColorInterp interp = band->GetColorInterpretation();
        if (interp == GCI_RedBand || interp == GCI_GreenBand ||
            interp == GCI_BlueBand || interp == GCI_GrayIndex ||
            interp == GCI_PaletteIndex) {
            selected.push_back(i);
        }
    }

    if (selected.empty()) {
        for (int i = 1; i <= std::min(3, bandCount); ++i) selected.push_back(i);
    }
    return selected;
}

int ToByte(double value)
{
    if (!std::isfinite(value)) return 0;
    return std::max(0, std::min(255, static_cast<int>(std::lround(value))));
}

int PackPixel(const std::vector<std::vector<double>>& rows,
              std::size_t x,
              const std::vector<bool>& hasNoData,
              const std::vector<double>& noDataValues)
{
    int channels[3] = {0, 0, 0};
    const std::size_t used = std::min<std::size_t>(3, rows.size());
    bool nonBlack = false;
    for (std::size_t b = 0; b < used; ++b) {
        const double value = rows[b][x];
        if ((!hasNoData[b] || value != noDataValues[b]) &&
            std::isfinite(value)) {
            channels[b] = ToByte(value);
            nonBlack = nonBlack || channels[b] != 0;
        }
    }
    if (!nonBlack) return 0;

    // A one-band mask is still a valid label image. Replicate it so its
    // packed label remains stable and can be polygonized like RGB data.
    if (used == 1) channels[1] = channels[2] = channels[0];
    else if (used == 2) channels[2] = channels[1];
    return (channels[0] << 16) | (channels[1] << 8) | channels[2];
}

void SetDatasetSpatialInfo(GDALDataset* dataset,
                           const double geotransform[6],
                           const char* projection)
{
    dataset->SetGeoTransform(const_cast<double*>(geotransform));
    if (projection && projection[0] != '\0') dataset->SetProjection(projection);
}

// 对大连通域，通过检测窄腰判断是否拆分
// whether the distance-transform core contains multiple separated peaks.
bool SplitNarrowWaistComponent(const cv::Mat& componentMask,
                               cv::Mat& outputLabels,
                               int& nextLabel,
                               int parentColor,
                               std::unordered_map<int, int>& labelToParent,
                               long long& narrowWaistSplitCount)
{
    CV_Assert(componentMask.type() == CV_8U);
    if (cv::countNonZero(componentMask) == 0) return false;

    cv::Mat dist;
    cv::distanceTransform(componentMask, dist, cv::DIST_L2, 5);
    double minVal = 0.0, maxVal = 0.0;
    cv::Point minLoc, maxLoc;
    cv::minMaxLoc(dist, &minVal, &maxVal, &minLoc, &maxLoc);
    if (maxVal <= 1.0) return false;

    const double area = static_cast<double>(cv::countNonZero(componentMask));
    const double coreThreshold = std::max(1.5, 0.42 * maxVal);
    cv::Mat coreMask = dist >= coreThreshold;
    cv::Mat coreLabels, coreStats, coreCentroids;
    const int coreCount = cv::connectedComponentsWithStats(
        coreMask, coreLabels, coreStats, coreCentroids, 8, CV_32S);
    if (coreCount < 3) return false; // 0背景 + 1个核心，不够切

    std::vector<SplitCandidate> candidates;
    candidates.reserve(static_cast<std::size_t>(coreCount - 1));
    for (int c = 1; c < coreCount; ++c) {
        const int coreArea = coreStats.at<int>(c, cv::CC_STAT_AREA);
        if (coreArea <= 0) continue;
        const double coreAreaRatio = static_cast<double>(coreArea) / std::max(area, 1.0);
        if (coreAreaRatio < kNarrowWaistMinCoreAreaRatio) continue;
        cv::Mat coreComponent = coreLabels == c;
        double cmin = 0.0, cmax = 0.0;
        cv::Point cminLoc, cmaxLoc;
        cv::minMaxLoc(dist, &cmin, &cmax, &cminLoc, &cmaxLoc, coreComponent);
        if (cmax < kNarrowWaistMinSplitGapPixels) continue;
        SplitCandidate cand;
        cand.seedIndex = c;
        cand.seedRow = cmaxLoc.y;
        cand.seedCol = cmaxLoc.x;
        cand.coreArea = static_cast<double>(coreArea);
        cand.coreWidth = cmax;
        candidates.push_back(cand);
    }

    if (candidates.size() < 2) return false;

    std::sort(candidates.begin(), candidates.end(),
              [](const SplitCandidate& a, const SplitCandidate& b) {
                  if (a.coreWidth != b.coreWidth) return a.coreWidth > b.coreWidth;
                  return a.coreArea > b.coreArea;
              });
    if (candidates.front().coreWidth < kNarrowWaistMinSplitGapPixels) return false;
    if (candidates[1].coreWidth <
        candidates.front().coreWidth * kNarrowWaistMaxWidthRatio) {
        return false;
    }

    cv::Mat marker = cv::Mat::zeros(componentMask.size(), CV_32S);
    int localSeed = 1;
    const double minSeedSeparation = std::max(
        3.0, 0.10 * std::hypot(componentMask.cols, componentMask.rows));
    std::vector<SplitCandidate> selected;
    for (const auto& cand : candidates) {
        bool separated = true;
        for (const auto& chosen : selected) {
            if (std::hypot(cand.seedCol - chosen.seedCol,
                           cand.seedRow - chosen.seedRow) < minSeedSeparation) {
                separated = false;
                break;
            }
        }
        if (!separated) continue;
        selected.push_back(cand);
        marker.at<int>(cand.seedRow, cand.seedCol) = localSeed++;
        if (localSeed > 4) break;
    }
    if (localSeed <= 2) return false;

    cv::Mat seedPixels;
    cv::compare(marker, 0, seedPixels, cv::CMP_GT);
    cv::Mat inverseSeeds;
    cv::bitwise_not(seedPixels, inverseSeeds);
    cv::Mat nearest;
    cv::distanceTransform(inverseSeeds, dist, nearest, cv::DIST_L2, 5, cv::DIST_LABEL_PIXEL);

    std::unordered_map<int, int> nearestToLocalSeed;
    for (int y = 0; y < marker.rows; ++y) {
        const int* markerRow = marker.ptr<int>(y);
        const int* nearestRow = nearest.ptr<int>(y);
        for (int x = 0; x < marker.cols; ++x) {
            if (markerRow[x] > 0 && nearestRow[x] > 0) {
                nearestToLocalSeed[nearestRow[x]] = markerRow[x];
            }
        }
    }

    std::vector<int> outputForSeed(static_cast<std::size_t>(localSeed + 1), 0);
    for (int s = 1; s < localSeed; ++s) {
        outputForSeed[s] = nextLabel++;
        labelToParent[outputForSeed[s]] = parentColor;
    }

    bool changed = false;
    for (int y = 0; y < componentMask.rows; ++y) {
        const uchar* compRow = componentMask.ptr<uchar>(y);
        const int* nearestRow = nearest.ptr<int>(y);
        int* outRow = outputLabels.ptr<int>(y);
        for (int x = 0; x < componentMask.cols; ++x) {
            if (!compRow[x]) continue;
            const auto it = nearestToLocalSeed.find(nearestRow[x]);
            const int s = it == nearestToLocalSeed.end() ? 0 : it->second;
            if (s > 0 && s < localSeed) {
                outRow[x] = outputForSeed[s];
                changed = true;
            }
        }
    }

    if (changed) ++narrowWaistSplitCount;
    return changed;
}

// 构建标签图并保留原始前景。腐蚀用于
// used only to find stable cores; it is never used as the final footprint.
cv::Mat SeparateColorComponents(const cv::Mat& colorImage,
                                int erosionRadiusPixels,
                                double pixelArea,
                                long long& seedCount,
                                long long& splitComponentCount,
                                long long& narrowWaistSplitCount,
                                std::unordered_map<int, int>& labelToParent)
{
    CV_Assert(colorImage.type() == CV_32S);
    cv::Mat kernel = cv::getStructuringElement(
        cv::MORPH_ELLIPSE,
        cv::Size(2 * erosionRadiusPixels + 1, 2 * erosionRadiusPixels + 1));

    std::unordered_map<int, PixelBounds> colorBounds;
    colorBounds.reserve(64);
    for (int y = 0; y < colorImage.rows; ++y) {
        const int* row = colorImage.ptr<int>(y);
        for (int x = 0; x < colorImage.cols; ++x) {
            if (row[x] != 0) colorBounds[row[x]].include(x, y);
        }
    }

    cv::Mat seedMarkers = cv::Mat::zeros(colorImage.size(), CV_32S);
    int nextSeed = 1;
    std::vector<int> seedColorVec;  // seedColorVec[g] = 全局种子 g 所属的颜色(实例)
    seedColorVec.push_back(0);      // 占位，种子编号从 1 开始
    const int minSeedPixels = std::max(
        1, static_cast<int>(std::ceil(kMinSeedAreaSquareMeters /
                                     std::max(pixelArea, 1e-9))));
    for (const auto& colorItem : colorBounds) {
        const int color = colorItem.first;
        const cv::Rect roi = colorItem.second.rect();
        const cv::Mat colorRoi = colorImage(roi);
        cv::Mat markerRoi = seedMarkers(roi);
        cv::Mat colorMask;
        cv::compare(colorRoi, color, colorMask, cv::CMP_EQ);

        cv::Mat eroded;
        cv::erode(colorMask, eroded, kernel, cv::Point(-1, -1), 1,
                  cv::BORDER_CONSTANT, cv::Scalar(0));
        cv::Mat seedComponents;
        cv::Mat seedStats;
        cv::Mat seedCentroids;
        const int seedCountForColor = cv::connectedComponentsWithStats(
            eroded, seedComponents, seedStats, seedCentroids, 8, CV_32S);
        std::vector<int> localToGlobal(static_cast<std::size_t>(seedCountForColor), 0);
        for (int seed = 1; seed < seedCountForColor; ++seed) {
            if (seedStats.at<int>(seed, cv::CC_STAT_AREA) >= minSeedPixels) {
                localToGlobal[seed] = nextSeed;
                seedColorVec.push_back(color);  // 记录该全局种子的颜色(实例)
                ++nextSeed;
            }
        }
        for (int y = 0; y < seedComponents.rows; ++y) {
            const int* sourceRow = seedComponents.ptr<int>(y);
            int* markerRow = markerRoi.ptr<int>(y);
            for (int x = 0; x < seedComponents.cols; ++x) {
                const int localSeed = sourceRow[x];
                if (localSeed > 0 && localToGlobal[localSeed] > 0) {
                    markerRow[x] = localToGlobal[localSeed];
                }
            }
        }
    }

    seedCount = nextSeed - 1;
    cv::Mat binaryMask;
    cv::compare(colorImage, 0, binaryMask, cv::CMP_GT);
    cv::Mat originalComponents;
    cv::Mat originalStats;
    cv::Mat originalCentroids;
    const int originalCount = cv::connectedComponentsWithStats(
        binaryMask, originalComponents, originalStats, originalCentroids, 8, CV_32S);
    cv::Mat result = cv::Mat::zeros(colorImage.size(), CV_32S);
    std::vector<int> fallbackLabel(static_cast<std::size_t>(originalCount), 0);
    int nextOutputLabel = 1;
    for (int component = 1; component < originalCount; ++component) {
        fallbackLabel[component] = nextOutputLabel++;
        // fallback 像素(跨分量归属的边界)不参与同源合并：留空，后续按小面积过滤
    }
    if (seedCount == 0) {
        for (int y = 0; y < result.rows; ++y) {
            const int* originals = originalComponents.ptr<int>(y);
            int* output = result.ptr<int>(y);
            for (int x = 0; x < result.cols; ++x) {
                if (originals[x] > 0) output[x] = fallbackLabel[originals[x]];
            }
        }
        return result;
    }

    cv::Mat seedPixels;
    cv::compare(seedMarkers, 0, seedPixels, cv::CMP_GT);
    cv::Mat inverseSeeds;
    cv::bitwise_not(seedPixels, inverseSeeds);
    cv::Mat distances;
    cv::Mat nearestSeed;
    cv::distanceTransform(inverseSeeds, distances, nearestSeed,
                          cv::DIST_L2, 5, cv::DIST_LABEL_CCOMP);

    // distanceTransform labels the connected zero regions in the same scan
    // order, but map explicitly so this does not depend on that detail.
    std::vector<int> nearestToMarker(static_cast<std::size_t>(seedCount + 1), 0);
    std::vector<int> seedToOriginal(static_cast<std::size_t>(seedCount + 1), 0);
    for (int y = 0; y < result.rows; ++y) {
        const int* markers = seedMarkers.ptr<int>(y);
        const int* nearest = nearestSeed.ptr<int>(y);
        const int* originals = originalComponents.ptr<int>(y);
        for (int x = 0; x < result.cols; ++x) {
            if (markers[x] > 0 && nearest[x] > 0) {
                nearestToMarker[nearest[x]] = markers[x];
                seedToOriginal[markers[x]] = originals[x];
            }
        }
    }

    std::vector<int> outputForSeed(static_cast<std::size_t>(seedCount + 1), 0);
    std::vector<int> seedsPerOriginal(static_cast<std::size_t>(originalCount), 0);
    for (int seed = 1; seed <= seedCount; ++seed) {
        outputForSeed[seed] = nextOutputLabel++;
        const int component = seedToOriginal[seed];
        if (component > 0) ++seedsPerOriginal[component];
        // parent 用"颜色(实例)"而非二值连通分量：同色=同一栋→合并，跨色=不同栋→分开
        if (seed < static_cast<int>(seedColorVec.size())) {
            labelToParent[outputForSeed[seed]] = seedColorVec[seed];
        }
    }
    for (int component = 1; component < originalCount; ++component) {
        if (seedsPerOriginal[component] > 1) ++splitComponentCount;
    }

    for (int y = 0; y < result.rows; ++y) {
        const int* originals = originalComponents.ptr<int>(y);
        const int* nearest = nearestSeed.ptr<int>(y);
        int* output = result.ptr<int>(y);
        for (int x = 0; x < result.cols; ++x) {
            const int component = originals[x];
            if (component <= 0) continue;
            const int nearestLabel = nearest[x];
            const int seed = nearestLabel > 0 && nearestLabel <= seedCount
                ? nearestToMarker[nearestLabel] : 0;
            output[x] = seed > 0 && seedToOriginal[seed] == component
                ? outputForSeed[seed] : fallbackLabel[component];
        }
    }

    const int minWaistSplitPixels = std::max(80, 4 * minSeedPixels);
    for (int component = 1; component < originalCount; ++component) {
        if (seedsPerOriginal[component] > 1) continue;
        const int areaPixels = originalStats.at<int>(component, cv::CC_STAT_AREA);
        if (areaPixels < minWaistSplitPixels) continue;

        const int x = originalStats.at<int>(component, cv::CC_STAT_LEFT);
        const int y = originalStats.at<int>(component, cv::CC_STAT_TOP);
        const int w = originalStats.at<int>(component, cv::CC_STAT_WIDTH);
        const int h = originalStats.at<int>(component, cv::CC_STAT_HEIGHT);
        if (w <= 0 || h <= 0) continue;
        cv::Rect roi(x, y, w, h);
        cv::Mat originalRoi = originalComponents(roi);
        cv::Mat componentMask;
        cv::compare(originalRoi, component, componentMask, cv::CMP_EQ);

        int parentColor = 0;
        const cv::Mat colorRoi = colorImage(roi);
        for (int yy = 0; yy < componentMask.rows && parentColor == 0; ++yy) {
            const uchar* maskRow = componentMask.ptr<uchar>(yy);
            const int* colorRow = colorRoi.ptr<int>(yy);
            for (int xx = 0; xx < componentMask.cols; ++xx) {
                if (maskRow[xx] && colorRow[xx] != 0) {
                    parentColor = colorRow[xx];
                    break;
                }
            }
        }
        if (parentColor == 0) continue;

        cv::Mat resultRoi = result(roi);
        SplitNarrowWaistComponent(
            componentMask, resultRoi, nextOutputLabel, parentColor,
            labelToParent, narrowWaistSplitCount);
    }
    return result;
}

} // namespace

bool VectorizeBuildingMask(const std::string& tifPath,
                           const std::string& outputShpPath,
                           MaskVectorizationStats& stats)
{
    stats = {};
    DatasetPtr source(static_cast<GDALDataset*>(
        GDALOpenEx(tifPath.c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr)));
    if (!source) {
        std::cerr << "[Mask] Cannot open GeoTIFF: " << tifPath << std::endl;
        return false;
    }

    stats.width = source->GetRasterXSize();
    stats.height = source->GetRasterYSize();
    stats.sourceBands = source->GetRasterCount();
    const std::vector<int> colorBands = SelectColorBands(source.get());
    stats.usedBands = static_cast<int>(std::min<std::size_t>(3, colorBands.size()));
    if (stats.width <= 0 || stats.height <= 0 || colorBands.empty()) {
        std::cerr << "[Mask] Raster has no usable color bands." << std::endl;
        return false;
    }

    double geotransform[6] = {};
    if (source->GetGeoTransform(geotransform) != CE_None) {
        std::cerr << "[Mask] GeoTIFF has no valid geotransform." << std::endl;
        return false;
    }
    stats.pixelSizeX = std::hypot(geotransform[1], geotransform[2]);
    stats.pixelSizeY = std::hypot(geotransform[4], geotransform[5]);
    const double pixelSize = std::max(1e-9, (stats.pixelSizeX + stats.pixelSizeY) * 0.5);
    stats.erosionRadiusPixels = std::max(
        1, static_cast<int>(std::ceil(kSeedErosionDistanceMeters / pixelSize)));

    const double cornerX[4] = {
        geotransform[0],
        geotransform[0] + stats.width * geotransform[1],
        geotransform[0] + stats.height * geotransform[2],
        geotransform[0] + stats.width * geotransform[1] + stats.height * geotransform[2]
    };
    const double cornerY[4] = {
        geotransform[3],
        geotransform[3] + stats.width * geotransform[4],
        geotransform[3] + stats.height * geotransform[5],
        geotransform[3] + stats.width * geotransform[4] + stats.height * geotransform[5]
    };
    stats.minX = *std::min_element(cornerX, cornerX + 4);
    stats.maxX = *std::max_element(cornerX, cornerX + 4);
    stats.minY = *std::min_element(cornerY, cornerY + 4);
    stats.maxY = *std::max_element(cornerY, cornerY + 4);

    cv::Mat colorImage(stats.height, stats.width, CV_32S, cv::Scalar(0));
    std::vector<std::vector<double>> sourceRows(
        stats.usedBands, std::vector<double>(static_cast<std::size_t>(stats.width)));
    std::vector<double> noDataValues(stats.usedBands, 0.0);
    std::vector<bool> hasNoData(stats.usedBands, false);
    for (int b = 0; b < stats.usedBands; ++b) {
        int valid = FALSE;
        noDataValues[b] = source->GetRasterBand(colorBands[b])->GetNoDataValue(&valid);
        hasNoData[b] = valid != FALSE;
    }

    // GDAL mask band / validity mask (per color band). Any band explicitly
    // marking a pixel invalid => background, regardless of channel values.
    GDALRasterBand* firstBand = source->GetRasterBand(colorBands[0]);
    GDALRasterBand* maskBand = firstBand ? firstBand->GetMaskBand() : nullptr;
    const bool hasGdalMask =
        maskBand && maskBand->GetMaskFlags() == GMF_PER_DATASET;
    std::vector<uint8_t> maskRow(
        hasGdalMask ? static_cast<std::size_t>(stats.width) : 0);

    long long exactBlack = 0;
    long long nearBlackBg = 0;
    long long noDataBg = 0;
    long long maskBg = 0;
    for (int y = 0; y < stats.height; ++y) {
        for (int b = 0; b < stats.usedBands; ++b) {
            if (source->GetRasterBand(colorBands[b])->RasterIO(
                    GF_Read, 0, y, stats.width, 1, sourceRows[b].data(),
                    stats.width, 1, GDT_Float64, 0, 0, nullptr) != CE_None) {
                std::cerr << "[Mask] Failed reading raster row " << y << "." << std::endl;
                return false;
            }
        }
        if (hasGdalMask &&
            maskBand->RasterIO(GF_Read, 0, y, stats.width, 1,
                               maskRow.data(), stats.width, 1, GDT_Byte,
                               0, 0, nullptr) != CE_None) {
            std::cerr << "[Mask] Failed reading mask band row " << y << "." << std::endl;
            return false;
        }
        int* row = colorImage.ptr<int>(y);
        for (int x = 0; x < stats.width; ++x) {
            if (hasGdalMask && maskRow[static_cast<std::size_t>(x)] == 0) {
                row[x] = 0;
                ++maskBg;
                continue;
            }
            row[x] = PackPixel(sourceRows, static_cast<std::size_t>(x),
                               hasNoData, noDataValues);
            if (row[x] == 0) {
                // Background via NoData or exact black: distinguish for log.
                bool anyNoData = false;
                for (int b = 0; b < stats.usedBands; ++b) {
                    const double v = sourceRows[b][x];
                    if (hasNoData[b] && v == noDataValues[b]) { anyNoData = true; break; }
                }
                if (anyNoData) ++noDataBg;
                else ++exactBlack;
                continue;
            }
            // Near-black foreground suppression: all channels <= threshold.
            const int r = (row[x] >> 16) & 0xFF;
            const int g = (row[x] >> 8) & 0xFF;
            const int bch = row[x] & 0xFF;
            if (r <= kMaskBackgroundMaxChannel &&
                g <= kMaskBackgroundMaxChannel &&
                bch <= kMaskBackgroundMaxChannel) {
                row[x] = 0;
                ++nearBlackBg;
            }
            if (row[x] != 0) ++stats.buildingPixels;
        }
    }
    std::cout << "[Mask background] rule=near_black_max_channel_"
              << kMaskBackgroundMaxChannel
              << " exact_black=" << exactBlack
              << " near_black_bg=" << nearBlackBg
              << " nodata_bg=" << noDataBg
              << " gdal_mask_bg=" << maskBg
              << " foreground=" << stats.buildingPixels << std::endl;
    if (stats.buildingPixels == 0) {
        std::cerr << "[Mask] No non-black building pixels were found." << std::endl;
        return false;
    }

    std::unordered_set<int> colors;
    colors.reserve(64);
    bool tooManyColors = false;
    {
        for (int y = 0; y < colorImage.rows; ++y) {
            int* row = colorImage.ptr<int>(y);
            for (int x = 0; x < colorImage.cols; ++x) {
                if (row[x] == 0) continue;
                if (!tooManyColors) {
                    colors.insert(row[x]);
                    tooManyColors = colors.size() > kMaxReliableColorLabels;
                }
                if (tooManyColors) row[x] = 1;
            }
        }
    }
    if (tooManyColors) {
        // 影像可能含抗锯齿或连续颜色
        // than semantic instance labels. Normalize the prefix scanned before
        // the threshold was reached as well.
        colorImage.setTo(1, colorImage != 0);
    }
    stats.sourceColorCount = static_cast<long long>(colors.size());
    stats.colorLabelsPreserved = !tooManyColors;
    if (tooManyColors) {
        std::cout << "[Mask warn] color count exceeded " << kMaxReliableColorLabels
                  << " — the mask may have been resampled with bilinear/cubic"
                  << " interpolation instead of nearest-neighbour."
                  << " Instance labels degrade; consider re-cropping with"
                  << " gdalwarp -r near." << std::endl;
    }
    if (!tooManyColors) {
        RemoveHairlineSeamComponents(
            colorImage, stats.pixelSizeX * stats.pixelSizeY);
        RemoveBatchWrapperRings(
            colorImage, stats.pixelSizeX * stats.pixelSizeY);
    }
    std::unordered_map<int, int> labelToParent;
    cv::Mat separated = SeparateColorComponents(
        colorImage, stats.erosionRadiusPixels,
        stats.pixelSizeX * stats.pixelSizeY,
        stats.seedCount, stats.splitComponentCount,
        stats.narrowWaistSplitCount, labelToParent);

    GDALDriver* memDriver = GetGDALDriverManager()->GetDriverByName("MEM");
    if (!memDriver) {
        std::cerr << "[Mask] GDAL MEM driver is unavailable." << std::endl;
        return false;
    }
    DatasetPtr labels(memDriver->Create("", stats.width, stats.height, 1,
                                        GDT_Int32, nullptr));
    DatasetPtr valid(memDriver->Create("", stats.width, stats.height, 1,
                                      GDT_Byte, nullptr));
    if (!labels || !valid) {
        std::cerr << "[Mask] Cannot create in-memory label raster." << std::endl;
        return false;
    }
    const char* projection = source->GetProjectionRef();
    SetDatasetSpatialInfo(labels.get(), geotransform, projection);
    SetDatasetSpatialInfo(valid.get(), geotransform, projection);

    cv::Mat validMask;
    cv::compare(colorImage, 0, validMask, cv::CMP_GT);
    if (labels->GetRasterBand(1)->RasterIO(
            GF_Write, 0, 0, stats.width, stats.height, separated.ptr<int>(),
            stats.width, stats.height, GDT_Int32, 0, 0, nullptr) != CE_None ||
        valid->GetRasterBand(1)->RasterIO(
            GF_Write, 0, 0, stats.width, stats.height, validMask.ptr(),
            stats.width, stats.height, GDT_Byte, 0, 0, nullptr) != CE_None) {
        std::cerr << "[Mask] Failed writing separated label raster." << std::endl;
        return false;
    }
    valid->GetRasterBand(1)->SetNoDataValue(0.0);

    GDALDriver* shpDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    if (!shpDriver) {
        std::cerr << "[Mask] ESRI Shapefile driver is unavailable." << std::endl;
        return false;
    }
    const std::filesystem::path outputPath(outputShpPath);
    if (!outputPath.parent_path().empty()) {
        std::filesystem::create_directories(outputPath.parent_path());
    }
    if (!RemoveShapefileFamily(outputPath)) {
        std::cerr << "[Mask] Cannot replace existing Shapefile. Close it in GIS software first: "
                  << outputShpPath << std::endl;
        return false;
    }

    DatasetPtr output(shpDriver->Create(outputShpPath.c_str(), 0, 0, 0,
                                        GDT_Unknown, nullptr));
    if (!output) {
        std::cerr << "[Mask] Cannot create initial outline Shapefile: " << outputShpPath << std::endl;
        return false;
    }
    OGRSpatialReference spatialRef;
    OGRSpatialReference* spatialRefPtr = nullptr;
    if (projection && projection[0] != '\0' &&
        spatialRef.importFromWkt(projection) == OGRERR_NONE) {
        spatialRefPtr = &spatialRef;
    }
    OGRLayer* layer = output->CreateLayer("initial_building_outline", spatialRefPtr,
                                          wkbPolygon, nullptr);
    if (!layer) {
        std::cerr << "[Mask] Cannot create initial outline layer." << std::endl;
        return false;
    }
    OGRFieldDefn valueField("mask", OFTInteger);
    if (layer->CreateField(&valueField) != OGRERR_NONE) {
        std::cerr << "[Mask] Cannot create mask attribute." << std::endl;
        return false;
    }
    OGRFieldDefn parentField("parent", OFTInteger);
    if (layer->CreateField(&parentField) != OGRERR_NONE) {
        std::cerr << "[Mask] Cannot create parent attribute." << std::endl;
        return false;
    }
    // 矢量化时一次性分配的稳定多边形标识。
    // survives the merge/split stages (CopyFieldValues propagates it to parts
    // and merged features) while FIDs churn, so initial outlines, debug
    // hypotheses and the regularized output can be joined on it.
    OGRFieldDefn idField("id", OFTInteger64);
    if (layer->CreateField(&idField) != OGRERR_NONE) {
        std::cerr << "[Mask] Cannot create id attribute." << std::endl;
        return false;
    }
    if (GDALPolygonize(labels->GetRasterBand(1), valid->GetRasterBand(1),
                       layer, 0, nullptr, nullptr, nullptr) != CE_None) {
        std::cerr << "[Mask] GDALPolygonize failed." << std::endl;
        return false;
    }
    // 给每个多边形写入它所属的"原始连通分量(parent)"，供后续按同源合并过分割碎片
    {
        OGRFeatureDefn* defn = layer->GetLayerDefn();
        const int maskIdx = defn->GetFieldIndex("mask");
        const int parentIdx = defn->GetFieldIndex("parent");
        const int idIdx = defn->GetFieldIndex("id");
        if (maskIdx >= 0 && parentIdx >= 0 && idIdx >= 0) {
            layer->ResetReading();
            while (OGRFeature* feature = layer->GetNextFeature()) {
                const int maskValue = feature->GetFieldAsInteger(maskIdx);
                const auto it = labelToParent.find(maskValue);
                if (it != labelToParent.end()) {
                    feature->SetField(parentIdx, it->second);
                }
                feature->SetField(idIdx, feature->GetFID());
                layer->SetFeature(feature);
                OGRFeature::DestroyFeature(feature);
            }
        }
    }
    stats.containedPolygonsRemoved = RemoveContainedPolygons(layer);
    stats.polygonCount = layer->GetFeatureCount(TRUE);
    output->FlushCache();
    return stats.polygonCount > 0;
}

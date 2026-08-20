// =============================================================================
// inspect_mask_nesting.cpp
// 检查掩膜中"分幅叠加"造成的嵌套结构: 一个颜色的大连通域(外批次粗轮廓)
// 内部包含若干其他颜色的小连通域(内批次细轮廓), 外层仅比内层并集略宽。
// 用法: inspect_mask_nesting <mask.tif> [--window x y w h]
// 输出: 连通域统计、包含关系、"薄壳wrapper"判定与案例。
// =============================================================================
#include <gdal_priv.h>
#include <cpl_conv.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <unordered_map>
#include <vector>

namespace {

struct Component {
    int64_t id = -1;          // color id (r<<16|g<<8|b)
    int compId = -1;          // connected component index
    int64_t pixels = 0;
    int minX = 0, minY = 0, maxX = 0, maxY = 0;
    int containerComp = -1;   // surrounding component (-1 background/-2 mixed)
    double containerRingFrac = 0.0;
};

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::printf("usage: %s <mask.tif> [--window x y w h]\n", argv[0]);
        return 1;
    }
    GDALAllRegister();
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(argv[1], GA_ReadOnly));
    if (!ds) {
        std::printf("cannot open %s\n", argv[1]);
        return 1;
    }
    const int width = ds->GetRasterXSize();
    const int height = ds->GetRasterYSize();
    const int bands = ds->GetRasterCount();
    double adf[6] = {0};
    ds->GetGeoTransform(adf);
    std::printf("raster %dx%d bands=%d pixel=%.3fx%.3f type=%s\n",
                width, height, bands, adf[1], adf[5],
                GDALGetDataTypeName(ds->GetRasterBand(1)->GetRasterDataType()));

    // 读取前3个band为uint32颜色 (背景=0)
    std::vector<uint8_t> r(static_cast<size_t>(width) * height);
    std::vector<uint8_t> g(r.size()), b(r.size());
    if (ds->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, width, height,
            r.data(), width, height, GDT_Byte, 0, 0) != CE_None ||
        (bands > 1 && ds->GetRasterBand(2)->RasterIO(GF_Read, 0, 0, width, height,
            g.data(), width, height, GDT_Byte, 0, 0) != CE_None) ||
        (bands > 2 && ds->GetRasterBand(3)->RasterIO(GF_Read, 0, 0, width, height,
            b.data(), width, height, GDT_Byte, 0, 0) != CE_None)) {
        std::printf("read failed\n");
        return 1;
    }
    GDALClose(ds);

    std::vector<uint32_t> color(static_cast<size_t>(width) * height);
    for (size_t i = 0; i < color.size(); ++i) {
        color[i] = bands > 1
            ? (static_cast<uint32_t>(r[i]) << 16) |
              (static_cast<uint32_t>(g[i]) << 8) | static_cast<uint32_t>(b[i])
            : static_cast<uint32_t>(r[i]);
    }

    // 4邻域连通域标记
    std::vector<int32_t> comp(color.size(), -1);
    std::vector<Component> comps;
    std::vector<std::vector<int32_t>> stackRefs;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t idx = static_cast<size_t>(y) * width + x;
            if (color[idx] == 0 || comp[idx] >= 0) continue;
            const int compId = static_cast<int>(comps.size());
            Component c;
            c.compId = compId;
            c.id = color[idx];
            c.minX = x; c.maxX = x; c.minY = y; c.maxY = y;
            std::vector<int32_t> stack;
            stack.push_back(static_cast<int32_t>(idx));
            comp[idx] = compId;
            while (!stack.empty()) {
                const int32_t cur = stack.back();
                stack.pop_back();
                const int cy = cur / width;
                const int cx = cur % width;
                c.pixels++;
                c.minX = std::min(c.minX, cx); c.maxX = std::max(c.maxX, cx);
                c.minY = std::min(c.minY, cy); c.maxY = std::max(c.maxY, cy);
                const int32_t nb[4] = {cur - 1, cur + 1, cur - width, cur + width};
                const bool valid[4] = {cx > 0, cx < width - 1, cy > 0, cy < height - 1};
                for (int k = 0; k < 4; ++k) {
                    if (!valid[k]) continue;
                    const int32_t n = nb[k];
                    if (color[static_cast<size_t>(n)] == color[idx] && comp[n] < 0) {
                        comp[n] = compId;
                        stack.push_back(n);
                    }
                }
            }
            comps.push_back(c);
        }
    }
    std::printf("前景连通域: %zu 个\n", comps.size());

    // 包含关系: 对每个连通域, 检查其边界外一圈像素属于哪个连通域
    std::unordered_map<int, int> containerVotes;
    for (auto& c : comps) {
        std::unordered_map<int, int> votes;
        int ring = 0;
        for (int y = c.minY; y <= c.maxY; ++y) {
            for (int x = c.minX; x <= c.maxX; ++x) {
                const size_t idx = static_cast<size_t>(y) * width + x;
                if (comp[idx] != c.compId) continue;
                // 4邻是否属于自己
                const int nbx[8] = {x-1,x+1,x,x,y-1,y-1,y+1,y+1};
                const int nby[8] = {y,y,y-1,y+1,y-1,y+1,y-1,y+1};
                for (int k = 0; k < 8; ++k) {
                    const int nx = nbx[k], ny = nby[k];
                    if (nx < 0 || ny < 0 || nx >= width || ny >= height) continue;
                    const size_t nidx = static_cast<size_t>(ny) * width + nx;
                    if (comp[nidx] == c.compId) continue;
                    ring++;
                    votes[comp[nidx]]++;   // -1=背景 或 其他连通域
                }
            }
        }
        if (ring == 0) continue;
        int bestId = -2, bestCount = -1;
        for (const auto& v : votes) {
            if (v.second > bestCount) { bestCount = v.second; bestId = v.first; }
        }
        if (bestId >= 0 && static_cast<double>(bestCount) / ring >= 0.90) {
            c.containerComp = bestId;
            c.containerRingFrac = static_cast<double>(bestCount) / ring;
        }
    }

    // wrapper判定: 外层C包含至少1个内层, 且外层像素中未被内层并集覆盖的比例小
    // (薄壳=外层只比内层并集宽1~3像素)。先建 compId -> covered 累计。
    std::vector<int64_t> coveredByInner(comps.size(), 0);
    std::vector<int> innerCount(comps.size(), 0);
    for (const auto& c : comps) {
        if (c.containerComp < 0) continue;
        coveredByInner[c.containerComp] += c.pixels;
        innerCount[c.containerComp]++;
    }
    int wrapperCount = 0, containedCount = 0, bigGapCount = 0;
    std::printf("\n嵌套结构 (外层含>=1其他颜色连通域):\n");
    for (const auto& c : comps) {
        if (innerCount[c.compId] == 0) continue;
        containedCount += innerCount[c.compId];
        const double uncoveredFrac =
            1.0 - static_cast<double>(coveredByInner[c.compId]) /
                      static_cast<double>(c.pixels);
        const bool thinShell = uncoveredFrac <= 0.35;
        if (thinShell) wrapperCount++; else bigGapCount++;
        if (wrapperCount + bigGapCount <= 25) {
            std::printf("  外层comp=%d color=%06llx px=%lld bbox=[%d,%d %dx%d] 内层数=%d "
                        "内层占外层面积比=%.2f uncovered=%.2f %s\n",
                        c.compId,
                        static_cast<unsigned long long>(c.id),
                        static_cast<long long>(c.pixels),
                        c.minX, c.minY, c.maxX - c.minX + 1, c.maxY - c.minY + 1,
                        innerCount[c.compId],
                        static_cast<double>(coveredByInner[c.compId]) /
                            static_cast<double>(c.pixels),
                        uncoveredFrac,
                        thinShell ? "<< 薄壳wrapper(批次叠加伪影)" : "(大间隙,可能是院落等真实包含)");
        }
    }
    std::printf("\n统计: 含内层的外层=%d (薄壳wrapper=%d, 大间隙=%d), 被包含内层总数=%d\n",
                wrapperCount + bigGapCount, wrapperCount, bigGapCount, containedCount);

    // 可选窗口转储
    if (argc >= 7 && std::strcmp(argv[2], "--window") == 0) {
        const int wx = atoi(argv[3]), wy = atoi(argv[4]);
        const int ww = atoi(argv[5]), wh = atoi(argv[6]);
        std::printf("\n窗口[%d,%d %dx%d] 颜色分布:\n", wx, wy, ww, wh);
        std::unordered_map<uint32_t, int> cnt;
        for (int y = wy; y < wy + wh && y < height; ++y) {
            for (int x = wx; x < wx + ww && x < width; ++x) {
                cnt[color[static_cast<size_t>(y) * width + x]]++;
            }
        }
        std::vector<std::pair<uint32_t, int>> v(cnt.begin(), cnt.end());
        std::sort(v.begin(), v.end(),
                  [](auto& a, auto& b) { return a.second > b.second; });
        for (const auto& p : v) {
            if (p.second < 5) break;
            std::printf("  color=%06x px=%d\n", p.first, p.second);
        }
    }
    return 0;
}

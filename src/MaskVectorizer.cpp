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
        // Search from the smallest larger polygon to the largest so nested
        // structures are assigned to their immediate containing shell.
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
            // Polygonized label islands appear as both a child polygon and a
            // hole in the parent. Test against the exterior shell so that the
            // parent's corresponding hole does not hide this containment.
            if (outer.exteriorShell->Contains(inner.geometry.get())) {
                inner.remove = true;
                inner.containerIndex = outerIndex;
                break;
            }
        }
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

// This is deliberately expressed in metres. The pixel radius is derived from
// the GeoTIFF geotransform so the same setting works for different rasters.
constexpr double kSeedErosionDistanceMeters = 1.2;
constexpr double kMinSeedAreaSquareMeters = 5.0;
// Keep instance colors in normal RGB masks. The fallback is only for images
// whose values are effectively continuous (for example, a rendered overlay).
constexpr std::size_t kMaxReliableColorLabels = 65536;
constexpr double kNarrowWaistMaxWidthRatio = 0.38;
constexpr double kNarrowWaistMinCoreAreaRatio = 0.10;
constexpr double kNarrowWaistMinSplitGapPixels = 2.0;

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

// For a large single connected component, detect a narrow waist by checking
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

// Build a label image while retaining the original foreground. Erosion is
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

    for (int y = 0; y < stats.height; ++y) {
        for (int b = 0; b < stats.usedBands; ++b) {
            if (source->GetRasterBand(colorBands[b])->RasterIO(
                    GF_Read, 0, y, stats.width, 1, sourceRows[b].data(),
                    stats.width, 1, GDT_Float64, 0, 0, nullptr) != CE_None) {
                std::cerr << "[Mask] Failed reading raster row " << y << "." << std::endl;
                return false;
            }
        }
        int* row = colorImage.ptr<int>(y);
        for (int x = 0; x < stats.width; ++x) {
            row[x] = PackPixel(sourceRows, static_cast<std::size_t>(x),
                               hasNoData, noDataValues);
            if (row[x] != 0) ++stats.buildingPixels;
        }
    }
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
        // The image likely contains anti-aliased or continuous colors rather
        // than semantic instance labels. Normalize the prefix scanned before
        // the threshold was reached as well.
        colorImage.setTo(1, colorImage != 0);
    }
    stats.sourceColorCount = static_cast<long long>(colors.size());
    stats.colorLabelsPreserved = !tooManyColors;
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
    // Stable per-polygon identifier assigned once at vectorization time. It
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

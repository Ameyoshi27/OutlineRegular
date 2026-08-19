#pragma once

#include <cstddef>
#include <string>

struct MaskVectorizationStats {
    int width = 0;
    int height = 0;
    int sourceBands = 0;
    int usedBands = 0;
    long long sourceColorCount = 0;
    bool colorLabelsPreserved = false;
    long long seedCount = 0;
    long long splitComponentCount = 0;
    long long narrowWaistSplitCount = 0;
    int erosionRadiusPixels = 0;
    double pixelSizeX = 0.0;
    double pixelSizeY = 0.0;
    std::size_t buildingPixels = 0;
    long long polygonCount = 0;
    long long containedPolygonsRemoved = 0;
    double minX = 0.0;
    double minY = 0.0;
    double maxX = 0.0;
    double maxY = 0.0;
}; 

// Converts a color building mask GeoTIFF into separated building polygons.
// Color labels are preserved. Narrow same-label connections are split by
// eroded cores, while the original foreground pixels are retained.
// The output Shapefile uses the raster geotransform and spatial reference.
bool VectorizeBuildingMask(const std::string& tifPath,
                           const std::string& outputShpPath,
                           MaskVectorizationStats& stats);

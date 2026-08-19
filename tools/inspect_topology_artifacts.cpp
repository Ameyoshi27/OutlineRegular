#include <ogrsf_frmts.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

struct Point2 {
    double x = 0.0;
    double y = 0.0;
};

double dist(const Point2& a, const Point2& b)
{
    return std::hypot(a.x - b.x, a.y - b.y);
}

double signedArea(const std::vector<Point2>& pts)
{
    double area = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& a = pts[i];
        const auto& b = pts[(i + 1) % pts.size()];
        area += a.x * b.y - b.x * a.y;
    }
    return 0.5 * area;
}

double angleDiff90(double angle)
{
    angle = std::fmod(std::abs(angle), M_PI / 2.0);
    return std::min(angle, M_PI / 2.0 - angle);
}

double fold90(double angle)
{
    angle = std::fmod(angle, M_PI / 2.0);
    if (angle < 0.0) angle += M_PI / 2.0;
    return angle;
}

double foldedDistance(double a, double b)
{
    double d = std::fmod(std::abs(a - b), M_PI / 2.0);
    return std::min(d, M_PI / 2.0 - d);
}

std::vector<double> dominantAngles(const std::vector<Point2>& pts)
{
    std::vector<std::pair<double, double>> edges;
    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& a = pts[i];
        const auto& b = pts[(i + 1) % pts.size()];
        const double len = dist(a, b);
        if (len < 1e-6) continue;
        edges.push_back({fold90(std::atan2(b.y - a.y, b.x - a.x)), len});
    }
    std::vector<double> result;
    while (!edges.empty() && result.size() < 3) {
        double bestWeight = -1.0;
        double bestAngle = edges.front().first;
        for (const auto& e : edges) {
            double w = 0.0;
            for (const auto& f : edges) {
                if (foldedDistance(e.first, f.first) < 10.0 * M_PI / 180.0) w += f.second;
            }
            if (w > bestWeight) {
                bestWeight = w;
                bestAngle = e.first;
            }
        }
        result.push_back(bestAngle);
        edges.erase(std::remove_if(edges.begin(), edges.end(), [&](const auto& e) {
            return foldedDistance(e.first, bestAngle) < 12.0 * M_PI / 180.0;
        }), edges.end());
    }
    return result;
}

std::vector<Point2> exteriorPoints(const OGRGeometry* geometry)
{
    std::vector<Point2> pts;
    if (!geometry || wkbFlatten(geometry->getGeometryType()) != wkbPolygon) return pts;
    const auto* polygon = geometry->toPolygon();
    const auto* ring = polygon ? polygon->getExteriorRing() : nullptr;
    if (!ring) return pts;
    const int n = ring->getNumPoints();
    for (int i = 0; i < n; ++i) {
        Point2 p{ring->getX(i), ring->getY(i)};
        if (!pts.empty() && dist(pts.front(), p) < 1e-8 && i + 1 == n) continue;
        pts.push_back(p);
    }
    return pts;
}

void inspectRing(const std::vector<Point2>& pts, int fid)
{
    if (pts.size() < 3) return;
    const double area = std::abs(signedArea(pts));
    const double scale = std::sqrt(std::max(area, 1.0));
    const double shortEdge = std::clamp(0.025 * scale, 0.8, 3.0);
    const auto dirs = dominantAngles(pts);

    int shortEdges = 0;
    int offAxisEdges = 0;
    int sharpCorners = 0;
    int bevelCorners = 0;
    double minLen = 1e100;
    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& prev = pts[(i + pts.size() - 1) % pts.size()];
        const auto& cur = pts[i];
        const auto& next = pts[(i + 1) % pts.size()];
        const double len = dist(cur, next);
        minLen = std::min(minLen, len);
        if (len < shortEdge) ++shortEdges;
        const double a = std::atan2(next.y - cur.y, next.x - cur.x);
        double best = 1e100;
        for (double d : dirs) best = std::min(best, angleDiff90(a - d));
        if (best > 12.0 * M_PI / 180.0 && len < 2.0 * shortEdge) ++offAxisEdges;

        const double ax = prev.x - cur.x;
        const double ay = prev.y - cur.y;
        const double bx = next.x - cur.x;
        const double by = next.y - cur.y;
        const double na = std::hypot(ax, ay);
        const double nb = std::hypot(bx, by);
        if (na > 1e-8 && nb > 1e-8) {
            double c = (ax * bx + ay * by) / (na * nb);
            c = std::clamp(c, -1.0, 1.0);
            const double deg = std::acos(c) * 180.0 / M_PI;
            if (deg < 45.0) ++sharpCorners;
            if (std::abs(deg - 135.0) < 18.0 && std::min(na, nb) < 2.0 * shortEdge) ++bevelCorners;
        }
    }

    if (shortEdges || offAxisEdges || sharpCorners || bevelCorners) {
        std::cout << "fid=" << fid
                  << " vertices=" << pts.size()
                  << " area=" << area
                  << " minLen=" << minLen
                  << " shortEdges=" << shortEdges
                  << " offAxisShortEdges=" << offAxisEdges
                  << " sharpCorners=" << sharpCorners
                  << " bevelCorners=" << bevelCorners
                  << " dirs=";
        for (double d : dirs) std::cout << d * 180.0 / M_PI << ",";
        std::cout << "\n";
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: inspect_topology_artifacts <shp>\n";
        return 1;
    }
    GDALAllRegister();
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpenEx(argv[1], GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!ds) {
        std::cerr << "cannot open " << argv[1] << "\n";
        return 1;
    }
    OGRLayer* layer = ds->GetLayer(0);
    layer->ResetReading();
    int featureCount = 0;
    while (OGRFeature* feature = layer->GetNextFeature()) {
        ++featureCount;
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry && wkbFlatten(geometry->getGeometryType()) == wkbPolygon) {
            inspectRing(exteriorPoints(geometry), static_cast<int>(feature->GetFID()));
        }
        OGRFeature::DestroyFeature(feature);
    }
    std::cout << "features=" << featureCount << "\n";
    GDALClose(ds);
    return 0;
}

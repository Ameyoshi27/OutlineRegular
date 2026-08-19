#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr double kNarrowNeckMaxWidth = 4.1;
constexpr double kNarrowNeckMinBoundarySeparation = 4.0;
constexpr double kNarrowNeckMinPartArea = 20.0;

struct Point2D64 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct CandidateInfo {
    std::size_t a = 0;
    std::size_t b = 0;
    double width = 0.0;
    double arcA = 0.0;
    double arcB = 0.0;
    double areaA = 0.0;
    double areaB = 0.0;
    std::string reason;
};

double Distance2D64(const Point2D64& a, const Point2D64& b)
{
    return std::hypot(a.x - b.x, a.y - b.y);
}

std::vector<Point2D64> ExtractExteriorRing64(OGRPolygon* polygon)
{
    std::vector<Point2D64> ring;
    if (!polygon || !polygon->getExteriorRing()) return ring;
    OGRLinearRing* ogrRing = polygon->getExteriorRing();
    for (int i = 0; i < ogrRing->getNumPoints(); ++i) {
        ring.push_back({ogrRing->getX(i), ogrRing->getY(i), ogrRing->getZ(i)});
    }
    if (ring.size() >= 2 && Distance2D64(ring.front(), ring.back()) < 1e-9) {
        ring.pop_back();
    }
    return ring;
}

double GeometryArea(const OGRGeometry* geometry)
{
    return geometry ? OGR_G_Area(OGRGeometry::ToHandle(const_cast<OGRGeometry*>(geometry))) : 0.0;
}

bool PointInPolygon2D64(const Point2D64& p, const std::vector<Point2D64>& poly)
{
    bool inside = false;
    const std::size_t n = poly.size();
    if (n < 3) return false;
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const auto& pi = poly[i];
        const auto& pj = poly[j];
        const bool intersect = ((pi.y > p.y) != (pj.y > p.y)) &&
            (p.x < (pj.x - pi.x) * (p.y - pi.y) / ((pj.y - pi.y) + 1e-12) + pi.x);
        if (intersect) inside = !inside;
    }
    return inside;
}

double Orient2D64(const Point2D64& a, const Point2D64& b, const Point2D64& c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool OnSegment2D64(const Point2D64& a, const Point2D64& b, const Point2D64& p)
{
    constexpr double eps = 1e-8;
    if (std::abs(Orient2D64(a, b, p)) > eps) return false;
    return p.x >= std::min(a.x, b.x) - eps && p.x <= std::max(a.x, b.x) + eps &&
           p.y >= std::min(a.y, b.y) - eps && p.y <= std::max(a.y, b.y) + eps;
}

bool SegmentsIntersect2D64(
    const Point2D64& a,
    const Point2D64& b,
    const Point2D64& c,
    const Point2D64& d)
{
    constexpr double eps = 1e-8;
    const double o1 = Orient2D64(a, b, c);
    const double o2 = Orient2D64(a, b, d);
    const double o3 = Orient2D64(c, d, a);
    const double o4 = Orient2D64(c, d, b);

    if (((o1 > eps && o2 < -eps) || (o1 < -eps && o2 > eps)) &&
        ((o3 > eps && o4 < -eps) || (o3 < -eps && o4 > eps))) {
        return true;
    }
    if (std::abs(o1) <= eps && OnSegment2D64(a, b, c)) return true;
    if (std::abs(o2) <= eps && OnSegment2D64(a, b, d)) return true;
    if (std::abs(o3) <= eps && OnSegment2D64(c, d, a)) return true;
    if (std::abs(o4) <= eps && OnSegment2D64(c, d, b)) return true;
    return false;
}

bool EdgeTouchesVertex(std::size_t edgeIndex, std::size_t vertexIndex, std::size_t vertexCount)
{
    return edgeIndex == vertexIndex || ((edgeIndex + 1) % vertexCount) == vertexIndex;
}

bool NeckCutLineInsidePolygon2D64(
    const std::vector<Point2D64>& ring,
    std::size_t vertexA,
    std::size_t vertexB,
    std::string& reason)
{
    const std::size_t n = ring.size();
    if (n < 4 || vertexA >= n || vertexB >= n || vertexA == vertexB) {
        reason = "bad_indices";
        return false;
    }
    const Point2D64& a = ring[vertexA];
    const Point2D64& b = ring[vertexB];

    constexpr int sampleCount = 7;
    for (int s = 1; s < sampleCount; ++s) {
        const double t = static_cast<double>(s) / sampleCount;
        Point2D64 p;
        p.x = a.x * (1.0 - t) + b.x * t;
        p.y = a.y * (1.0 - t) + b.y * t;
        p.z = a.z * (1.0 - t) + b.z * t;
        if (!PointInPolygon2D64(p, ring)) {
            reason = "sample_outside@" + std::to_string(s) +
                "=(" + std::to_string(p.x) + "," + std::to_string(p.y) + ")";
            return false;
        }
    }

    for (std::size_t edge = 0; edge < n; ++edge) {
        if (EdgeTouchesVertex(edge, vertexA, n) ||
            EdgeTouchesVertex(edge, vertexB, n)) {
            continue;
        }
        if (SegmentsIntersect2D64(a, b, ring[edge], ring[(edge + 1) % n])) {
            reason = "intersects_edge_" + std::to_string(edge) + "_" + std::to_string((edge + 1) % n);
            return false;
        }
    }
    reason = "inside";
    return true;
}

OGRPolygon* MakePolygon64(const std::vector<Point2D64>& ring)
{
    if (ring.size() < 3) return nullptr;
    auto* polygon = new OGRPolygon();
    OGRLinearRing ogrRing;
    for (const auto& p : ring) {
        ogrRing.addPoint(p.x, p.y, p.z);
    }
    ogrRing.addPoint(ring.front().x, ring.front().y, ring.front().z);
    polygon->addRing(&ogrRing);
    return polygon;
}

std::vector<CandidateInfo> FindAllCandidates(const std::vector<Point2D64>& ring);

void DiagnosePair(const std::vector<Point2D64>& ring, int ia, int ib, const std::string& label)
{
    const std::size_t n = ring.size();
    std::cout << "\n[" << label << "] vertexA=" << ia << " vertexB=" << ib << "\n";
    if (ia < 0 || ib < 0 || static_cast<std::size_t>(ia) >= n || static_cast<std::size_t>(ib) >= n) {
        std::cout << "  invalid index for ring size=" << n << "\n";
        return;
    }
    if (ia == ib) {
        std::cout << "  same vertex\n";
        return;
    }
    std::size_t a = static_cast<std::size_t>(std::min(ia, ib));
    std::size_t b = static_cast<std::size_t>(std::max(ia, ib));

    std::vector<double> edgeLengths(n, 0.0);
    std::vector<double> prefix(n + 1, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        edgeLengths[i] = Distance2D64(ring[i], ring[(i + 1) % n]);
        prefix[i + 1] = prefix[i] + edgeLengths[i];
    }
    const double perimeter = prefix.back();
    const std::size_t forwardEdges = b - a;
    const std::size_t backwardEdges = n - forwardEdges;
    const double arcA = prefix[b] - prefix[a];
    const double arcB = perimeter - arcA;
    const double width = Distance2D64(ring[a], ring[b]);
    std::string insideReason;
    const bool inside = NeckCutLineInsidePolygon2D64(ring, a, b, insideReason);

    std::vector<Point2D64> ringA;
    std::vector<Point2D64> ringB;
    for (std::size_t index = a;; index = (index + 1) % n) {
        ringA.push_back(ring[index]);
        if (index == b) break;
    }
    for (std::size_t index = b;; index = (index + 1) % n) {
        ringB.push_back(ring[index]);
        if (index == a) break;
    }
    std::unique_ptr<OGRPolygon> polygonA(MakePolygon64(ringA));
    std::unique_ptr<OGRPolygon> polygonB(MakePolygon64(ringB));
    const double areaA = GeometryArea(polygonA.get());
    const double areaB = GeometryArea(polygonB.get());

    std::cout << "  A=(" << ring[a].x << "," << ring[a].y << ") B=(" << ring[b].x << "," << ring[b].y << ")\n";
    std::cout << "  width=" << width << " passWidth=" << (width > 1e-6 && width <= kNarrowNeckMaxWidth) << "\n";
    std::cout << "  arcA=" << arcA << " arcB=" << arcB
              << " passArc=" << (arcA >= kNarrowNeckMinBoundarySeparation && arcB >= kNarrowNeckMinBoundarySeparation) << "\n";
    std::cout << "  forwardEdges=" << forwardEdges << " backwardEdges=" << backwardEdges
              << " passNonAdjacent=" << (forwardEdges >= 2 && backwardEdges >= 2) << "\n";
    std::cout << "  inside=" << inside << " reason=" << insideReason << "\n";
    std::cout << "  areaA=" << areaA << " areaB=" << areaB
              << " passArea=" << (areaA >= kNarrowNeckMinPartArea && areaB >= kNarrowNeckMinPartArea) << "\n";
    const auto candidates = FindAllCandidates(ring);
    std::cout << "  validCandidates=" << candidates.size() << "\n";
    for (std::size_t k = 0; k < std::min<std::size_t>(candidates.size(), 8); ++k) {
        const auto& c = candidates[k];
        std::cout << "    #" << k
                  << " a=" << c.a
                  << " b=" << c.b
                  << " width=" << c.width
                  << " arcA=" << c.arcA
                  << " arcB=" << c.arcB
                  << " areaA=" << c.areaA
                  << " areaB=" << c.areaB
                  << "\n";
    }
}

std::vector<CandidateInfo> FindAllCandidates(const std::vector<Point2D64>& ring)
{
    std::vector<CandidateInfo> candidates;
    const std::size_t n = ring.size();
    if (n < 4) return candidates;

    std::vector<double> edgeLengths(n, 0.0);
    std::vector<double> prefix(n + 1, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        edgeLengths[i] = Distance2D64(ring[i], ring[(i + 1) % n]);
        prefix[i + 1] = prefix[i] + edgeLengths[i];
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const std::size_t forwardEdges = j - i;
            const std::size_t backwardEdges = n - forwardEdges;
            if (forwardEdges < 2 || backwardEdges < 2) continue;

            const double arcA = prefix[j] - prefix[i];
            const double arcB = prefix.back() - arcA;
            const double width = Distance2D64(ring[i], ring[j]);
            if (width <= 1e-6 || width > kNarrowNeckMaxWidth) continue;
            if (arcA < kNarrowNeckMinBoundarySeparation ||
                arcB < kNarrowNeckMinBoundarySeparation) continue;

            std::string reason;
            if (!NeckCutLineInsidePolygon2D64(ring, i, j, reason)) continue;

            std::vector<Point2D64> ringA;
            std::vector<Point2D64> ringB;
            for (std::size_t index = i;; index = (index + 1) % n) {
                ringA.push_back(ring[index]);
                if (index == j) break;
            }
            for (std::size_t index = j;; index = (index + 1) % n) {
                ringB.push_back(ring[index]);
                if (index == i) break;
            }

            std::unique_ptr<OGRPolygon> polygonA(MakePolygon64(ringA));
            std::unique_ptr<OGRPolygon> polygonB(MakePolygon64(ringB));
            CandidateInfo info;
            info.a = i;
            info.b = j;
            info.width = width;
            info.arcA = arcA;
            info.arcB = arcB;
            info.areaA = GeometryArea(polygonA.get());
            info.areaB = GeometryArea(polygonB.get());
            info.reason = reason;
            candidates.push_back(std::move(info));
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const CandidateInfo& lhs, const CandidateInfo& rhs) {
        if (std::abs(lhs.width - rhs.width) > 1e-9) return lhs.width < rhs.width;
        return std::min(lhs.arcA, lhs.arcB) > std::min(rhs.arcA, rhs.arcB);
    });
    return candidates;
}

OGRFeature* GetFeatureBySequentialIndex(OGRLayer* layer, int index)
{
    if (!layer || index < 0) return nullptr;
    layer->ResetReading();
    int current = 0;
    while (OGRFeature* feature = layer->GetNextFeature()) {
        if (current == index) return feature;
        OGRFeature::DestroyFeature(feature);
        ++current;
    }
    return nullptr;
}

void DiagnoseFeature(OGRFeature* feature, const std::string& label)
{
    if (!feature) {
        std::cout << "\n[" << label << "] feature not found\n";
        return;
    }
    OGRGeometry* geometry = feature->GetGeometryRef();
    OGRPolygon* polygon = geometry ? geometry->toPolygon() : nullptr;
    if (!polygon) {
        std::cout << "\n[" << label << "] not a polygon\n";
        return;
    }
    const auto ring = ExtractExteriorRing64(polygon);
    std::cout << "\n[" << label << "] fid=" << feature->GetFID()
              << " ringVertices=" << ring.size()
              << " area=" << GeometryArea(polygon) << "\n";
    DiagnosePair(ring, 139, 152, "0-based 139/152");
    DiagnosePair(ring, 138, 151, "1-based 139/152");
}

}  // namespace

int main(int argc, char** argv)
{
    std::cout << std::setprecision(15);
    const char* path = argc > 1 ? argv[1] : "D:/outlineRegular/outlineRegular/build_deps_release/Release/initial_building_outline.shp";
    GDALAllRegister();
    GDALDataset* dataset = static_cast<GDALDataset*>(
        GDALOpenEx(path, GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!dataset) {
        std::cerr << "Cannot open " << path << "\n";
        return 1;
    }
    OGRLayer* layer = dataset->GetLayer(0);
    if (!layer) {
        std::cerr << "No layer\n";
        GDALClose(dataset);
        return 1;
    }
    std::cout << "Layer feature count=" << layer->GetFeatureCount() << "\n";

    OGRFeature* fidFeature = layer->GetFeature(1332);
    DiagnoseFeature(fidFeature, "FID 1332");
    OGRFeature::DestroyFeature(fidFeature);

    OGRFeature* seq0Feature = GetFeatureBySequentialIndex(layer, 1332);
    DiagnoseFeature(seq0Feature, "sequential 0-based #1332");
    OGRFeature::DestroyFeature(seq0Feature);

    OGRFeature* seq1Feature = GetFeatureBySequentialIndex(layer, 1331);
    DiagnoseFeature(seq1Feature, "sequential 1-based #1332");
    OGRFeature::DestroyFeature(seq1Feature);

    GDALClose(dataset);
    return 0;
}

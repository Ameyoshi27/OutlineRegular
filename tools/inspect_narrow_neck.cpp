#include <ogrsf_frmts.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr double kMaxWidth = 2.0;
constexpr double kMinBoundarySeparation = 4.0;
constexpr double kBoundarySeparationRatio = 0.02;
constexpr double kMinPartArea = 20.0;
constexpr int kMaxCuts = 12;

struct Point {
    double x = 0.0;
    double y = 0.0;
};

struct Candidate {
    double ax = 0.0;
    double ay = 0.0;
    double bx = 0.0;
    double by = 0.0;
    double width = std::numeric_limits<double>::max();
    double boundarySeparation = 0.0;
    std::size_t i = 0;
    std::size_t j = 0;
};

struct Stats {
    long long edgePairs = 0;
    long long adjacent = 0;
    long long widthFail = 0;
    long long separationFail = 0;
    long long insideFail = 0;
    long long candidates = 0;
};

double dist(const Point& a, const Point& b)
{
    return std::hypot(a.x - b.x, a.y - b.y);
}

double clamp01(double v)
{
    return std::max(0.0, std::min(1.0, v));
}

bool pointInPolygon(const Point& p, const std::vector<Point>& poly)
{
    bool inside = false;
    const std::size_t n = poly.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const auto& pi = poly[i];
        const auto& pj = poly[j];
        if (((pi.y > p.y) != (pj.y > p.y)) &&
            (p.x < (pj.x - pi.x) * (p.y - pi.y) / ((pj.y - pi.y) + 1e-12) + pi.x)) {
            inside = !inside;
        }
    }
    return inside;
}

Candidate closestSegments(
    const Point& a0, const Point& a1,
    const Point& b0, const Point& b1)
{
    const double ux = a1.x - a0.x;
    const double uy = a1.y - a0.y;
    const double vx = b1.x - b0.x;
    const double vy = b1.y - b0.y;
    const double wx = a0.x - b0.x;
    const double aa = ux * ux + uy * uy;
    const double bb = ux * vx + uy * vy;
    const double cc = vx * vx + vy * vy;
    const double dd = ux * wx;
    const double ee = vx * wx;
    const double denom = aa * cc - bb * bb;
    double s = 0.0;
    double t = 0.0;
    if (aa <= 1e-12 && cc <= 1e-12) {
        s = t = 0.0;
    } else if (aa <= 1e-12) {
        s = 0.0;
        t = clamp01(ee / cc);
    } else if (cc <= 1e-12) {
        t = 0.0;
        s = clamp01(-dd / aa);
    } else if (std::abs(denom) > 1e-12) {
        s = clamp01((bb * ee - cc * dd) / denom);
        t = clamp01((aa * ee - bb * dd) / denom);
    } else {
        s = 0.0;
        t = clamp01(ee / cc);
    }
    const double px = a0.x + s * ux;
    const double py = a0.y + s * uy;
    t = cc > 1e-12 ? clamp01(((px - b0.x) * vx + (py - b0.y) * vy) / cc) : 0.0;
    const double qx = b0.x + t * vx;
    const double qy = b0.y + t * vy;
    s = aa > 1e-12 ? clamp01(((qx - a0.x) * ux + (qy - a0.y) * uy) / aa) : 0.0;

    Candidate c;
    c.ax = a0.x + s * ux;
    c.ay = a0.y + s * uy;
    c.bx = b0.x + t * vx;
    c.by = b0.y + t * vy;
    c.width = std::hypot(c.ax - c.bx, c.ay - c.by);
    return c;
}

bool lineMostlyInside(const Candidate& c, const std::vector<Point>& ring)
{
    for (int i = 1; i <= 5; ++i) {
        const double t = static_cast<double>(i) / 6.0;
        Point p;
        p.x = c.ax * (1.0 - t) + c.bx * t;
        p.y = c.ay * (1.0 - t) + c.by * t;
        if (!pointInPolygon(p, ring)) return false;
    }
    return true;
}

std::vector<Point> exteriorRing(OGRPolygon* polygon)
{
    std::vector<Point> ring;
    if (!polygon || !polygon->getExteriorRing()) return ring;
    auto* ogrRing = polygon->getExteriorRing();
    for (int i = 0; i < ogrRing->getNumPoints(); ++i) {
        ring.push_back({ogrRing->getX(i), ogrRing->getY(i)});
    }
    if (ring.size() >= 2 && dist(ring.front(), ring.back()) < 1e-9) {
        ring.pop_back();
    }
    return ring;
}

double area(OGRGeometry* g)
{
    return g ? OGR_G_Area(OGRGeometry::ToHandle(g)) : 0.0;
}

void collectParts(OGRGeometry* geometry, std::vector<std::unique_ptr<OGRGeometry>>& parts)
{
    if (!geometry || geometry->IsEmpty()) return;
    const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    if (type == wkbPolygon) {
        if (area(geometry) >= kMinPartArea) parts.emplace_back(geometry->clone());
    } else if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
        auto* c = geometry->toGeometryCollection();
        for (int i = 0; c && i < c->getNumGeometries(); ++i) {
            collectParts(c->getGeometryRef(i), parts);
        }
    }
}

bool findBest(const std::vector<Point>& ring, Candidate& best, Stats& stats)
{
    if (ring.size() < 8) return false;
    const std::size_t n = ring.size();
    std::vector<double> edgeLengths(n, 0.0);
    std::vector<double> prefix(n + 1, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        edgeLengths[i] = dist(ring[i], ring[(i + 1) % n]);
        prefix[i + 1] = prefix[i] + edgeLengths[i];
    }
    const double perimeter = prefix.back();
    const double minSep = std::max(kMinBoundarySeparation, perimeter * kBoundarySeparationRatio);

    bool found = false;
    double bestScore = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < n; ++i) {
        if (edgeLengths[i] < 1e-6) continue;
        for (std::size_t j = i + 2; j < n; ++j) {
            ++stats.edgePairs;
            if ((i == 0 && j + 1 >= n) || edgeLengths[j] < 1e-6) {
                ++stats.adjacent;
                continue;
            }
            const std::size_t edgeGap = std::min(j - i, n - (j - i));
            if (edgeGap < 3) {
                ++stats.adjacent;
                continue;
            }
            const double width = dist(ring[i], ring[j]);
            if (width <= 1e-6 || width > kMaxWidth) {
                ++stats.widthFail;
                continue;
            }
            const double along = prefix[j] - prefix[i];
            Candidate c;
            c.i = i;
            c.j = j;
            c.ax = ring[i].x;
            c.ay = ring[i].y;
            c.bx = ring[j].x;
            c.by = ring[j].y;
            c.width = width;
            c.boundarySeparation = std::min(along, perimeter - along);
            if (c.boundarySeparation < minSep) {
                ++stats.separationFail;
                continue;
            }
            ++stats.candidates;
            if (!found ||
                c.width < best.width - 1e-6 ||
                (std::abs(c.width - best.width) <= 1e-6 &&
                 c.boundarySeparation > best.boundarySeparation)) {
                bestScore = c.width;
                best = c;
                found = true;
            }
        }
    }
    return found;
}

double polygonArea(const std::vector<Point>& ring)
{
    if (ring.size() < 3) return 0.0;
    double a = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto& p = ring[i];
        const auto& q = ring[(i + 1) % ring.size()];
        a += p.x * q.y - q.x * p.y;
    }
    return std::abs(a) * 0.5;
}

bool splitRing(const std::vector<Point>& ring, std::vector<Point>& a, std::vector<Point>& b, Candidate& c)
{
    Stats stats;
    if (!findBest(ring, c, stats)) return false;
    const std::size_t n = ring.size();
    for (std::size_t index = c.i;; index = (index + 1) % n) {
        a.push_back(ring[index]);
        if (index == c.j) break;
    }
    for (std::size_t index = c.j;; index = (index + 1) % n) {
        b.push_back(ring[index]);
        if (index == c.i) break;
    }
    return a.size() >= 3 && b.size() >= 3 &&
           polygonArea(a) >= kMinPartArea && polygonArea(b) >= kMinPartArea;
}

void simulateRecursiveSplit(std::vector<Point> ring)
{
    std::vector<std::vector<Point>> pending;
    pending.push_back(std::move(ring));
    int cuts = 0;
    for (int pass = 0; pass < kMaxCuts && cuts < kMaxCuts; ++pass) {
        bool split = false;
        std::vector<std::vector<Point>> next;
        for (const auto& item : pending) {
            Candidate c;
            std::vector<Point> a;
            std::vector<Point> b;
            if (cuts < kMaxCuts && splitRing(item, a, b, c)) {
                std::cout << "  cut " << cuts + 1
                          << ": n=" << item.size()
                          << ", vertices=" << c.i << "," << c.j
                          << ", width=" << c.width
                          << ", sep=" << c.boundarySeparation
                          << ", area " << polygonArea(item)
                          << " -> " << polygonArea(a) << " + " << polygonArea(b)
                          << ", points " << a.size() << " + " << b.size()
                          << "\n";
                next.push_back(std::move(a));
                next.push_back(std::move(b));
                ++cuts;
                split = true;
            } else {
                next.push_back(item);
            }
        }
        pending = std::move(next);
        if (!split) break;
    }
    std::cout << "  final parts=" << pending.size() << "\n";
    for (std::size_t i = 0; i < pending.size(); ++i) {
        std::cout << "    part " << i << ": points=" << pending[i].size()
                  << ", area=" << polygonArea(pending[i]) << "\n";
    }
}

void inspectPolygon(OGRPolygon* polygon)
{
    const std::vector<Point> ring = exteriorRing(polygon);
    Stats stats;
    Candidate best;
    const bool has = findBest(ring, best, stats);
    OGREnvelope env;
    polygon->getEnvelope(&env);
    std::cout << "ring points=" << ring.size()
              << ", area=" << area(polygon)
              << ", bbox=(" << env.MinX << "," << env.MinY << ")-("
              << env.MaxX << "," << env.MaxY << ")\n";
    std::cout << "pairs=" << stats.edgePairs
              << ", adjacent=" << stats.adjacent
              << ", widthFail=" << stats.widthFail
              << ", sepFail=" << stats.separationFail
              << ", insideFail=" << stats.insideFail
              << ", candidates=" << stats.candidates << "\n";
    if (!has) {
        std::cout << "best candidate: none\n";
        simulateRecursiveSplit(ring);
        return;
    }
    std::cout << "best candidate: width=" << best.width
              << ", boundarySep=" << best.boundarySeparation
              << ", edges=" << best.i << "," << best.j
              << ", line=(" << best.ax << "," << best.ay << ")-("
              << best.bx << "," << best.by << ")\n";
    simulateRecursiveSplit(ring);
}

} // namespace

int main(int argc, char** argv)
{
    const char* path = argc > 1 ? argv[1] : "test\\1.shp";
    GDALAllRegister();
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpenEx(path, GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!ds) {
        std::cerr << "Cannot open " << path << "\n";
        return 1;
    }
    OGRLayer* layer = ds->GetLayer(0);
    int featureIndex = 0;
    layer->ResetReading();
    while (OGRFeature* f = layer->GetNextFeature()) {
        std::cout << "feature " << featureIndex++ << ", fid=" << f->GetFID() << "\n";
        OGRGeometry* g = f->GetGeometryRef();
        if (g && wkbFlatten(g->getGeometryType()) == wkbPolygon) {
            inspectPolygon(g->toPolygon());
        } else if (g && (wkbFlatten(g->getGeometryType()) == wkbMultiPolygon ||
                         wkbFlatten(g->getGeometryType()) == wkbGeometryCollection)) {
            auto* c = g->toGeometryCollection();
            for (int i = 0; i < c->getNumGeometries(); ++i) {
                if (wkbFlatten(c->getGeometryRef(i)->getGeometryType()) == wkbPolygon) {
                    std::cout << "  part " << i << "\n";
                    inspectPolygon(c->getGeometryRef(i)->toPolygon());
                }
            }
        }
        OGRFeature::DestroyFeature(f);
    }
    GDALClose(ds);
    return 0;
}

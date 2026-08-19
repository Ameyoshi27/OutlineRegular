#include <ogrsf_frmts.h>
#include <laszip/laszip_api.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

using GeometryPtr = std::unique_ptr<OGRGeometry>;

double area(const OGRGeometry* g)
{
    return g ? OGR_G_Area(OGRGeometry::ToHandle(const_cast<OGRGeometry*>(g))) : 0.0;
}

double length(OGRGeometry* g)
{
    return g ? OGR_G_Length(OGRGeometry::ToHandle(g)) : 0.0;
}

double perimeter(OGRGeometry* g)
{
    if (!g) return 0.0;
    GeometryPtr b(g->Boundary());
    return length(b.get());
}

GeometryPtr loadUnionGeometry(const std::string& path)
{
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpenEx(path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!ds) {
        std::cerr << "cannot open " << path << "\n";
        return nullptr;
    }
    OGRLayer* layer = ds->GetLayer(0);
    GeometryPtr result;
    layer->ResetReading();
    while (OGRFeature* f = layer->GetNextFeature()) {
        OGRGeometry* g = f->GetGeometryRef();
        if (g && !g->IsEmpty()) {
            if (!result) {
                result.reset(g->clone());
            } else {
                GeometryPtr u(result->Union(g));
                if (u && !u->IsEmpty()) result = std::move(u);
            }
        }
        OGRFeature::DestroyFeature(f);
    }
    GDALClose(ds);
    return result;
}

int countExteriorPoints(const OGRGeometry* g)
{
    if (!g) return 0;
    const OGRwkbGeometryType type = wkbFlatten(g->getGeometryType());
    if (type == wkbPolygon) {
        const auto* p = g->toPolygon();
        return p && p->getExteriorRing() ? p->getExteriorRing()->getNumPoints() : 0;
    }
    if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
        const auto* c = g->toGeometryCollection();
        int total = 0;
        for (int i = 0; c && i < c->getNumGeometries(); ++i) {
            total += countExteriorPoints(c->getGeometryRef(i));
        }
        return total;
    }
    return 0;
}

int countPolygons(const OGRGeometry* g)
{
    if (!g) return 0;
    const OGRwkbGeometryType type = wkbFlatten(g->getGeometryType());
    if (type == wkbPolygon) return 1;
    if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
        const auto* c = g->toGeometryCollection();
        int total = 0;
        for (int i = 0; c && i < c->getNumGeometries(); ++i) {
            total += countPolygons(c->getGeometryRef(i));
        }
        return total;
    }
    return 0;
}

void printSummary(const std::string& name, OGRGeometry* g)
{
    if (!g) {
        std::cout << name << ": empty\n";
        return;
    }
    OGREnvelope e;
    g->getEnvelope(&e);
    std::cout << name
              << ": type=" << OGRGeometryTypeToName(g->getGeometryType())
              << ", polygons=" << countPolygons(g)
              << ", exterior_points=" << countExteriorPoints(g)
              << ", area=" << area(g)
              << ", perimeter=" << perimeter(g)
              << ", bbox=(" << e.MinX << "," << e.MinY << ")-("
              << e.MaxX << "," << e.MaxY << ")\n";
}

std::vector<std::pair<double, double>> exteriorXY(const OGRGeometry* g)
{
    std::vector<std::pair<double, double>> pts;
    if (!g) return pts;
    const OGRwkbGeometryType type = wkbFlatten(g->getGeometryType());
    const OGRPolygon* p = nullptr;
    if (type == wkbPolygon) {
        p = g->toPolygon();
    } else if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
        const auto* c = g->toGeometryCollection();
        for (int i = 0; c && i < c->getNumGeometries(); ++i) {
            if (wkbFlatten(c->getGeometryRef(i)->getGeometryType()) == wkbPolygon) {
                p = c->getGeometryRef(i)->toPolygon();
                break;
            }
        }
    }
    if (!p || !p->getExteriorRing()) return pts;
    const auto* r = p->getExteriorRing();
    for (int i = 0; i < r->getNumPoints(); ++i) {
        pts.push_back({r->getX(i), r->getY(i)});
    }
    return pts;
}

void printVertices(const std::string& name, OGRGeometry* g)
{
    auto pts = exteriorXY(g);
    if (pts.size() >= 2) {
        const auto& a = pts.front();
        const auto& b = pts.back();
        if (std::hypot(a.first - b.first, a.second - b.second) < 1e-9) {
            pts.pop_back();
        }
    }
    std::cout << name << " vertices:\n";
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        const auto& q = pts[(i + 1) % pts.size()];
        const double len = std::hypot(q.first - p.first, q.second - p.second);
        double angle = std::atan2(q.second - p.second, q.first - p.first) * 180.0 / 3.141592653589793;
        if (angle < 0.0) angle += 180.0;
        if (angle >= 180.0) angle -= 180.0;
        std::cout << "  " << i << ": (" << p.first << "," << p.second
                  << "), edge_len=" << len << ", angle=" << angle << "\n";
    }
}

void compare(const std::string& aName, OGRGeometry* a,
             const std::string& bName, OGRGeometry* b)
{
    if (!a || !b) return;
    GeometryPtr inter(a->Intersection(b));
    GeometryPtr uni(a->Union(b));
    GeometryPtr aMinusB(a->Difference(b));
    GeometryPtr bMinusA(b->Difference(a));
    const double ia = area(inter.get());
    const double ua = std::max(area(uni.get()), 1e-9);
    std::cout << aName << " vs " << bName
              << ": intersection=" << ia
              << ", union=" << ua
              << ", IoU=" << ia / ua
              << ", " << aName << "-only=" << area(aMinusB.get())
              << ", " << bName << "-only=" << area(bMinusA.get()) << "\n";
    if (aMinusB && !aMinusB->IsEmpty()) {
        OGREnvelope e;
        aMinusB->getEnvelope(&e);
        std::cout << "  " << aName << "-only bbox=("
                  << e.MinX << "," << e.MinY << ")-("
                  << e.MaxX << "," << e.MaxY << ")\n";
    }
}

void inspectLas(const std::string& path)
{
    laszip_POINTER reader = nullptr;
    if (laszip_create(&reader) != 0) {
        std::cout << "LAS: cannot create reader\n";
        return;
    }
    laszip_BOOL isCompressed = FALSE;
    if (laszip_open_reader(reader, path.c_str(), &isCompressed) != 0) {
        laszip_CHAR* err = nullptr;
        laszip_get_error(reader, &err);
        std::cout << "LAS: cannot open: " << (err ? err : "?") << "\n";
        laszip_destroy(reader);
        return;
    }
    laszip_header_struct* header = nullptr;
    laszip_get_header_pointer(reader, &header);
    laszip_point_struct* point = nullptr;
    laszip_get_point_pointer(reader, &point);
    laszip_I64 n = header->extended_number_of_point_records > 0
        ? static_cast<laszip_I64>(header->extended_number_of_point_records)
        : static_cast<laszip_I64>(header->number_of_point_records);
    double minX = 1e100, minY = 1e100, minZ = 1e100;
    double maxX = -1e100, maxY = -1e100, maxZ = -1e100;
    for (laszip_I64 i = 0; i < n; ++i) {
        if (laszip_read_point(reader) != 0) break;
        laszip_F64 xyz[3];
        laszip_get_coordinates(reader, xyz);
        minX = std::min(minX, xyz[0]); maxX = std::max(maxX, xyz[0]);
        minY = std::min(minY, xyz[1]); maxY = std::max(maxY, xyz[1]);
        minZ = std::min(minZ, xyz[2]); maxZ = std::max(maxZ, xyz[2]);
    }
    std::cout << "support LAS: points=" << n
              << ", bbox=(" << minX << "," << minY << "," << minZ << ")-("
              << maxX << "," << maxY << "," << maxZ << ")"
              << ", header_offset=(" << header->x_offset << ","
              << header->y_offset << "," << header->z_offset << ")\n";
    laszip_close_reader(reader);
    laszip_destroy(reader);
}

} // namespace

int main()
{
    GDALAllRegister();
    auto initial = loadUnionGeometry("test\\1_initial.shp");
    auto best = loadUnionGeometry("test\\1_besthypothesis.shp");
    auto final = loadUnionGeometry("test\\1.shp");
    printSummary("initial", initial.get());
    printSummary("best", best.get());
    printSummary("final", final.get());
    compare("initial", initial.get(), "best", best.get());
    compare("initial", initial.get(), "final", final.get());
    compare("best", best.get(), "final", final.get());
    printVertices("best", best.get());
    printVertices("final", final.get());
    inspectLas("test\\1_debug_support_points.las");
    return 0;
}

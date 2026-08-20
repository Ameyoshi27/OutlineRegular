/*
 * @Descripttion: LOD底图规则化
 * @version:
 * @Author: JiangTao
 * @Date: 2025/08/24, 16:45:54
 * @LastEditors: JiangTao
 * @LastEditTime: 2025/09/24, 16:46:49
 */
/* ====================================================================
 * 本文件实现 outlineRegular 类，是建筑物轮廓规则化的核心实现。
 *
 * 核心目标：将点云提取出的不规则多边形轮廓，纠正成主要由直角/正交边
 *          组成的规整轮廓（即"规则化"）。
 *
 * regular_Contour() 的整体处理流程：
 *   1. 输入点云分辨率估计，判断轮廓是否近似圆/椭圆（若是则直接返回）；
 *   2. 自适应计算能量权重 Lambda；
 *   3. generatePolygonalHypotheses：基于 VDP 算法生成多组候选假设（从
 *      3 边到 n 边的简化多边形）；
 *   4. computeTotalEnergy：对每个假设计算"数据能量+模型能量+规则性能量
 *      +支撑能量+DLG先验能量"，选出总能量最小的最优假设；
 *   5. removeCollinearPoints：移除接近 180° 的冗余共线顶点；
 *   6. optimizeWithHardConstraints：用 Ceres 做硬约束优化（平行/垂直
 *      边分组 + 中点锚定 + 跨层墙面对齐），重建正交顶点；
 *   7. OptimizeFinal_points：拓扑后处理修复（短边/尖刺/自交/平行阶梯等）；
 *   8. 越界检测与回退：若结果非法或越界则回退到优化前假设。
 *
 * 文件还包含：匿名命名空间内的 2D 几何工具函数、若干 Ceres 残差结构体、
 * 以及 CGAL 轮廓规则化封装等。
 * ==================================================================== */
#include "outlineRegular.h"
#include <pcl/io/vtk_io.h>
#include <pcl/common/angles.h>
#include <pcl/common/distances.h>
#include <pcl/common/centroid.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>  // for pcl::computeCentroid
#include <pcl/surface/convex_hull.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <Eigen/Dense>          // for PCA (assuming Eigen is available)
#include <limits>
#include <iostream>             // for cerr (already used)
#include <chrono>               // 计时
#include <set>
#include <array>
#include <algorithm>
#include <numeric>
#include <ceres/ceres.h>
#include <ceres/rotation.h> // 旋转矩阵函数
#include <limits> // 无限值
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Shape_regularization/regularize_contours.h>
#include <CGAL/linear_least_squares_fitting_2.h> // 线性最小二乘拟合
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2/Stop_above_cost_threshold.h>
//#include <CGAL/Polygon_vertical_orthogonalization_2.h>
//#include "lodgenerate.h"  // standalone tool does not depend on the original Qt LOD application
 //#include <CGAL/PCA_dimension_tags.h> // PCA
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Line_2 Line_2;
typedef std::vector<Point_2> Contour;
typedef CGAL::Shape_regularization::Contours::Multiple_directions_2<Kernel, Contour> Contour_directions;
typedef CGAL::Exact_predicates_exact_constructions_kernel K; 
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;

using namespace std;

namespace {

constexpr std::size_t kSupportDirectionMinPoints = 30;
constexpr std::size_t kSupportDirectionMinPairs = 300;
constexpr double kSupportDirectionPairRadius = 2.5;
constexpr double kSupportDirectionMinPairDistance = 0.35;
constexpr double kSupportDirectionStrongPeakRatio = 0.14; // Uniform +/-5 degree background is 11/90.
constexpr double kSupportDirectionCorrectionDeg = 9.0;
// Original-ring PCA gate: a raster-staircase outline's overall trend is a
// stable direction prior for elongated footprints (validated to ~0.5 deg on
// real cases). A direction selection that deviates beyond these gates needs
// strong wall evidence, otherwise it snaps back to the ring trend. Blocks the
// "VDP hull diagonal" failure mode (18 deg chosen for a ~12 deg building).
constexpr double kPcaDirectionReliableAxisRatio = 1.6; // elongation (std ratio) for a tight 5 deg gate
constexpr double kPcaDirectionUsableAxisRatio = 1.3;   // moderate elongation -> loose 12 deg gate
constexpr double kPcaDirectionGateStrongDeg = 5.0;
constexpr double kPcaDirectionGateWeakDeg = 12.0;

bool isSimplePolygon2D(const std::vector<pcl::PointXYZ>& pts);
double undirectedAngleDifference(double a, double b);
double foldedLineAngle90(double angle);
double foldedAngleDistance90(double a, double b);
void removeNearlyCollinearPoints2D(std::vector<pcl::PointXYZ>& polygon, double threshold_deg);

// ===== clampDouble =====
// 作用：将 value 限制在 [low, high] 区间内
double clampDouble(double value, double low, double high)
{
    return std::max(low, std::min(value, high));
}

// ===== normalizeAnglePi =====
// 作用：将角度归一化到 (-PI, PI] 区间
double normalizeAnglePi(double angle)
{
    while (angle <= -M_PI) angle += 2.0 * M_PI;
    while (angle > M_PI) angle -= 2.0 * M_PI;
    return angle;
}

// ===== polygonSignedArea2D =====
// 作用：计算 2D 多边形的有符号面积（正=逆时针，负=顺时针）
// 返回：带符号面积值；点数<3 时返回 0
double polygonSignedArea2D(const std::vector<pcl::PointXYZ>& pts)
{
    if (pts.size() < 3) return 0.0;

    double area = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        const auto& q = pts[(i + 1) % pts.size()];
        area += static_cast<double>(p.x) * q.y - static_cast<double>(q.x) * p.y;
    }
    return 0.5 * area;
}

// ===== polygonArea2D =====
// 作用：计算 2D 多边形的绝对面积
// 返回：面积的绝对值
double polygonArea2D(const std::vector<pcl::PointXYZ>& pts)
{
    return std::abs(polygonSignedArea2D(pts));
}

Polygon_2 toCgalPolygon2D(const std::vector<pcl::PointXYZ>& pts)
{
    Polygon_2 polygon;
    polygon.container().reserve(pts.size());

    for (size_t i = 0; i < pts.size(); ++i) {
        if (i + 1 == pts.size() &&
            std::abs(static_cast<double>(pts[i].x) - pts.front().x) < 1e-8 &&
            std::abs(static_cast<double>(pts[i].y) - pts.front().y) < 1e-8) {
            continue;
        }
        polygon.push_back(Polygon_2::Point_2(
            static_cast<double>(pts[i].x),
            static_cast<double>(pts[i].y)));
    }

    if (polygon.size() >= 3 &&
        polygon.orientation() == CGAL::CLOCKWISE) {
        polygon.reverse_orientation();
    }

    return polygon;
}

double polygonWithHolesArea2D(const Polygon_with_holes_2& polygon)
{
    double area = std::abs(CGAL::to_double(polygon.outer_boundary().area()));
    for (auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole) {
        area -= std::abs(CGAL::to_double(hole->area()));
    }
    return std::max(0.0, area);
}

double polygonIoU2D(const std::vector<pcl::PointXYZ>& a,
    const std::vector<pcl::PointXYZ>& b)
{
    const double area_a = polygonArea2D(a);
    const double area_b = polygonArea2D(b);
    if (area_a < 1e-6 || area_b < 1e-6) return 0.0;
    if (!isSimplePolygon2D(a) || !isSimplePolygon2D(b)) return 0.0;

    Polygon_2 poly_a = toCgalPolygon2D(a);
    Polygon_2 poly_b = toCgalPolygon2D(b);
    if (poly_a.size() < 3 || poly_b.size() < 3 ||
        !poly_a.is_simple() || !poly_b.is_simple()) {
        return 0.0;
    }

    std::vector<Polygon_with_holes_2> intersections;
    try {
        CGAL::intersection(poly_a, poly_b, std::back_inserter(intersections));
    }
    catch (...) {
        return 0.0;
    }

    double intersection_area = 0.0;
    for (const auto& polygon : intersections) {
        intersection_area += polygonWithHolesArea2D(polygon);
    }

    const double union_area = area_a + area_b - intersection_area;
    if (union_area < 1e-6) return 0.0;
    return clampDouble(intersection_area / union_area, 0.0, 1.0);
}

// ===== cross2D =====
// 作用：计算向量 ab 与 ac 的 2D 叉积（用于判断方向/共线）
// 返回：叉积值
double cross2D(const pcl::PointXYZ& a, const pcl::PointXYZ& b, const pcl::PointXYZ& c)
{
    return (static_cast<double>(b.x) - a.x) * (static_cast<double>(c.y) - a.y) -
        (static_cast<double>(b.y) - a.y) * (static_cast<double>(c.x) - a.x);
}

// ===== onSegment2D =====
// 作用：判断点 p 是否在线段 ab 上（含端点，容差 eps）
// 返回：true 表示在线段上
bool onSegment2D(const pcl::PointXYZ& a, const pcl::PointXYZ& b, const pcl::PointXYZ& p)
{
    const double eps = 1e-8;
    if (std::abs(cross2D(a, b, p)) > eps) return false;
    return p.x >= std::min(a.x, b.x) - eps && p.x <= std::max(a.x, b.x) + eps &&
        p.y >= std::min(a.y, b.y) - eps && p.y <= std::max(a.y, b.y) + eps;
}

// ===== segmentsIntersect2D =====
// 作用：判断线段 ab 与线段 cd 是否在 2D 平面上相交（含端点重合/共线重叠）
// 参数：a,b - 第一条线段端点; c,d - 第二条线段端点
// 返回：true 表示两线段相交
bool segmentsIntersect2D(const pcl::PointXYZ& a, const pcl::PointXYZ& b,
    const pcl::PointXYZ& c, const pcl::PointXYZ& d)
{
    const double eps = 1e-8;
    double c1 = cross2D(a, b, c);
    double c2 = cross2D(a, b, d);
    double c3 = cross2D(c, d, a);
    double c4 = cross2D(c, d, b);

    if (((c1 > eps && c2 < -eps) || (c1 < -eps && c2 > eps)) &&
        ((c3 > eps && c4 < -eps) || (c3 < -eps && c4 > eps))) {
        return true;
    }

    return onSegment2D(a, b, c) || onSegment2D(a, b, d) ||
        onSegment2D(c, d, a) || onSegment2D(c, d, b);
}

// ===== isSimplePolygon2D =====
// 作用：判断多边形是否为简单多边形（无自交边）
// 返回：true 表示无自交（是简单多边形）
bool isSimplePolygon2D(const std::vector<pcl::PointXYZ>& pts)
{
    const size_t n = pts.size();
    if (n < 3) return false;

    for (size_t i = 0; i < n; ++i) {
        size_t i2 = (i + 1) % n;
        for (size_t j = i + 1; j < n; ++j) {
            size_t j2 = (j + 1) % n;
            if (i == j || i2 == j || j2 == i) continue;
            if (i == 0 && j2 == 0) continue;

            if (segmentsIntersect2D(pts[i], pts[i2], pts[j], pts[j2])) {
                return false;
            }
        }
    }
    return true;
}

// ===== removeDuplicatePoints2D =====
// 作用：移除 2D 多边形中距离过近的重复点（包括首尾闭合重复点）
// 参数：epsilon - 判定重复的距离阈值
void removeDuplicatePoints2D(std::vector<pcl::PointXYZ>& pts, float epsilon)
{
    if (pts.empty()) return;

    std::vector<pcl::PointXYZ> out;
    out.reserve(pts.size());
    const double eps2 = static_cast<double>(epsilon) * epsilon;

    auto same2D = [eps2](const pcl::PointXYZ& a, const pcl::PointXYZ& b) {
        const double dx = static_cast<double>(a.x) - b.x;
        const double dy = static_cast<double>(a.y) - b.y;
        return dx * dx + dy * dy <= eps2;
    };

    for (const auto& p : pts) {
        if (out.empty() || !same2D(out.back(), p)) {
            out.push_back(p);
        }
    }

    if (out.size() > 1 && same2D(out.front(), out.back())) {
        out.pop_back();
    }

    pts.swap(out);
}

// ===== circularBinDistance =====
// 作用：计算环形直方图中两个 bin 之间的最短距离（考虑环绕）
int circularBinDistance(int a, int b, int size)
{
    int diff = std::abs(a - b);
    return std::min(diff, size - diff);
}

// ===== dominantLineAngles2D =====
// 作用：基于加权（边长）直方图提取多边形的主要方向角（最多 max_peaks 个）
// 参数：pts - 多边形顶点; max_peaks - 最多提取的主方向数
// 返回：主方向角列表（弧度，映射到 [0, PI/2)）
std::vector<double> dominantLineAngles2D(const std::vector<pcl::PointXYZ>& pts, int max_peaks = 3)
{
    const int bin_size = 180;
    std::vector<double> histogram(bin_size, 0.0);

    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        const auto& q = pts[(i + 1) % pts.size()];
        double dx = q.x - p.x;
        double dy = q.y - p.y;
        double len = std::hypot(dx, dy);
        if (len < 1e-6) continue;

        double angle_deg = std::atan2(dy, dx) * 180.0 / M_PI;
        while (angle_deg < 0.0) angle_deg += 90.0;
        while (angle_deg >= 90.0) angle_deg -= 90.0;

        int bin = static_cast<int>(angle_deg * 2.0);
        if (bin >= 0 && bin < bin_size) histogram[bin] += len;
    }

    std::vector<double> smoothed(bin_size, 0.0);
    for (int i = 0; i < bin_size; ++i) {
        int prev = (i - 1 + bin_size) % bin_size;
        int next = (i + 1) % bin_size;
        smoothed[i] = 0.25 * histogram[prev] + 0.5 * histogram[i] + 0.25 * histogram[next];
    }

    double total_weight = 0.0;
    for (double w : smoothed) total_weight += w;

    std::vector<int> bins(bin_size);
    for (int i = 0; i < bin_size; ++i) bins[i] = i;
    std::sort(bins.begin(), bins.end(), [&](int lhs, int rhs) {
        return smoothed[lhs] > smoothed[rhs];
        });

    std::vector<int> selected_bins;
    const int min_separation_bins = 24; // 12 degrees at 0.5 degree/bin.
    double min_peak_weight = std::max(total_weight * 0.18, 1e-6);

    for (int bin : bins) {
        if (selected_bins.empty()) {
            selected_bins.push_back(bin);
        }
        else {
            if (smoothed[bin] < min_peak_weight) break;

            bool separated = true;
            for (int chosen : selected_bins) {
                if (circularBinDistance(bin, chosen, bin_size) < min_separation_bins) {
                    separated = false;
                    break;
                }
            }
            if (separated) selected_bins.push_back(bin);
        }

        if (static_cast<int>(selected_bins.size()) >= max_peaks) break;
    }

    std::vector<double> peaks;
    for (int bin : selected_bins) {
        peaks.push_back((static_cast<double>(bin) / 2.0) * M_PI / 180.0);
    }
    if (peaks.empty()) {
        peaks.push_back(0.0);
    }
    return peaks;
}

struct HausdorffMbrCandidate {
    double angle = 0.0;
    double area = 0.0;
    double mean_distance = 0.0;
    double q90_distance = 0.0;
    double covered_ratio = 0.0;
    double score = std::numeric_limits<double>::max();
};

std::vector<pcl::PointXYZ> convexHull2DMonotonic(std::vector<pcl::PointXYZ> points)
{
    std::vector<pcl::PointXYZ> hull;
    if (points.size() < 3) return hull;

    std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) {
        if (std::abs(static_cast<double>(a.x) - b.x) > 1e-6) return a.x < b.x;
        return a.y < b.y;
    });
    points.erase(std::unique(points.begin(), points.end(), [](const auto& a, const auto& b) {
        return std::hypot(static_cast<double>(a.x) - b.x, static_cast<double>(a.y) - b.y) < 1e-6;
    }), points.end());
    if (points.size() < 3) return hull;

    auto cross = [](const pcl::PointXYZ& o, const pcl::PointXYZ& a, const pcl::PointXYZ& b) {
        return (static_cast<double>(a.x) - o.x) * (static_cast<double>(b.y) - o.y) -
            (static_cast<double>(a.y) - o.y) * (static_cast<double>(b.x) - o.x);
    };

    std::vector<pcl::PointXYZ> lower;
    for (const auto& p : points) {
        while (lower.size() >= 2 && cross(lower[lower.size() - 2], lower.back(), p) <= 0.0) {
            lower.pop_back();
        }
        lower.push_back(p);
    }

    std::vector<pcl::PointXYZ> upper;
    for (auto it = points.rbegin(); it != points.rend(); ++it) {
        while (upper.size() >= 2 && cross(upper[upper.size() - 2], upper.back(), *it) <= 0.0) {
            upper.pop_back();
        }
        upper.push_back(*it);
    }

    lower.pop_back();
    upper.pop_back();
    hull = lower;
    hull.insert(hull.end(), upper.begin(), upper.end());
    return hull;
}

double distancePointToMbrBoundary(
    const pcl::PointXYZ& p,
    double angle,
    double min_u,
    double max_u,
    double min_v,
    double max_v)
{
    const double ux = std::cos(angle);
    const double uy = std::sin(angle);
    const double vx = -uy;
    const double vy = ux;
    const double u = static_cast<double>(p.x) * ux + static_cast<double>(p.y) * uy;
    const double v = static_cast<double>(p.x) * vx + static_cast<double>(p.y) * vy;
    const double du = std::max({ min_u - u, 0.0, u - max_u });
    const double dv = std::max({ min_v - v, 0.0, v - max_v });
    if (du > 0.0 || dv > 0.0) return std::hypot(du, dv);
    return std::min({ u - min_u, max_u - u, v - min_v, max_v - v });
}

bool computeHausdorffMbrCandidate(
    const std::vector<pcl::PointXYZ>& polygon,
    double angle,
    double distance_threshold,
    HausdorffMbrCandidate& candidate)
{
    if (polygon.size() < 3) return false;
    angle = foldedLineAngle90(angle);
    const double ux = std::cos(angle);
    const double uy = std::sin(angle);
    const double vx = -uy;
    const double vy = ux;

    double min_u = std::numeric_limits<double>::max();
    double max_u = -std::numeric_limits<double>::max();
    double min_v = std::numeric_limits<double>::max();
    double max_v = -std::numeric_limits<double>::max();
    for (const auto& p : polygon) {
        const double u = static_cast<double>(p.x) * ux + static_cast<double>(p.y) * uy;
        const double v = static_cast<double>(p.x) * vx + static_cast<double>(p.y) * vy;
        min_u = std::min(min_u, u);
        max_u = std::max(max_u, u);
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
    }
    const double width = max_u - min_u;
    const double height = max_v - min_v;
    if (width < 1e-6 || height < 1e-6) return false;

    std::vector<double> distances;
    distances.reserve(polygon.size());
    double weighted_distance = 0.0;
    double total_length = 0.0;
    double covered_length = 0.0;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double length = std::hypot(static_cast<double>(b.x) - a.x,
            static_cast<double>(b.y) - a.y);
        if (length < 1e-6) continue;

        pcl::PointXYZ mid;
        mid.x = static_cast<float>(0.5 * (static_cast<double>(a.x) + b.x));
        mid.y = static_cast<float>(0.5 * (static_cast<double>(a.y) + b.y));
        mid.z = a.z;
        const double d0 = distancePointToMbrBoundary(a, angle, min_u, max_u, min_v, max_v);
        const double d1 = distancePointToMbrBoundary(b, angle, min_u, max_u, min_v, max_v);
        const double dm = distancePointToMbrBoundary(mid, angle, min_u, max_u, min_v, max_v);
        const double d = std::max(dm, 0.5 * (d0 + d1));
        distances.push_back(d);
        weighted_distance += d * length;
        total_length += length;
        if (d <= distance_threshold) covered_length += length;
    }
    if (distances.empty() || total_length < 1e-6) return false;
    std::sort(distances.begin(), distances.end());
    const size_t q90_index = std::min(distances.size() - 1,
        static_cast<size_t>(std::floor(0.90 * static_cast<double>(distances.size() - 1))));

    candidate.angle = angle;
    candidate.area = width * height;
    candidate.mean_distance = weighted_distance / total_length;
    candidate.q90_distance = distances[q90_index];
    candidate.covered_ratio = covered_length / total_length;
    candidate.score = candidate.mean_distance + 0.65 * candidate.q90_distance +
        distance_threshold * std::max(0.0, 0.70 - candidate.covered_ratio);
    return std::isfinite(candidate.score);
}

std::vector<double> hausdorffMbrLineAngles2D(
    const std::vector<pcl::PointXYZ>& polygon,
    double resolution,
    int max_angles = 1)
{
    std::vector<double> result;
    if (polygon.size() < 3 || max_angles <= 0) return result;

    const double area = std::max(polygonArea2D(polygon), 1.0);
    const double scale = std::sqrt(area);
    const double distance_threshold = clampDouble(
        std::max(2.5 * resolution, 0.018 * scale), 0.35, 2.0);

    std::vector<double> candidate_angles;
    auto append_angle = [&](double angle, double tolerance_deg) {
        angle = foldedLineAngle90(angle);
        for (double existing : candidate_angles) {
            if (foldedAngleDistance90(existing, angle) < tolerance_deg * M_PI / 180.0) return;
        }
        candidate_angles.push_back(angle);
    };

    const std::vector<pcl::PointXYZ> hull = convexHull2DMonotonic(polygon);
    for (size_t i = 0; i < hull.size(); ++i) {
        const auto& a = hull[i];
        const auto& b = hull[(i + 1) % hull.size()];
        const double length = std::hypot(static_cast<double>(b.x) - a.x,
            static_cast<double>(b.y) - a.y);
        if (length < std::max(1.0, 0.02 * scale)) continue;
        append_angle(std::atan2(static_cast<double>(b.y) - a.y,
            static_cast<double>(b.x) - a.x), 2.0);
    }
    for (double angle : dominantLineAngles2D(polygon, 3)) {
        append_angle(angle, 2.0);
    }
    pcl::PointCloud<pcl::PointXYZ>::Ptr polygonCloud(new pcl::PointCloud<pcl::PointXYZ>);
    polygonCloud->points.assign(polygon.begin(), polygon.end());
    double pcaAngle = 0.0;
    double pcaAxisRatio = 1.0;
    if (outlineRegular::estimatePcaDirection2D(polygonCloud, pcaAngle, pcaAxisRatio)) {
        append_angle(pcaAngle, 2.0);
    }
    if (candidate_angles.empty()) return dominantLineAngles2D(polygon, max_angles);

    std::vector<HausdorffMbrCandidate> scored;
    scored.reserve(candidate_angles.size());
    double min_area = std::numeric_limits<double>::max();
    for (double angle : candidate_angles) {
        HausdorffMbrCandidate candidate;
        if (!computeHausdorffMbrCandidate(polygon, angle, distance_threshold, candidate)) continue;
        min_area = std::min(min_area, candidate.area);
        scored.push_back(candidate);
    }
    if (scored.empty()) return dominantLineAngles2D(polygon, max_angles);

    for (auto& candidate : scored) {
        const double area_ratio = candidate.area / std::max(min_area, 1e-6);
        candidate.score += distance_threshold * 0.20 * std::max(0.0, area_ratio - 1.0);
    }
    std::sort(scored.begin(), scored.end(), [](const auto& a, const auto& b) {
        return a.score < b.score;
    });

    for (const auto& candidate : scored) {
        bool duplicate = false;
        for (double kept : result) {
            if (foldedAngleDistance90(kept, candidate.angle) < 10.0 * M_PI / 180.0) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) continue;
        result.push_back(candidate.angle);
        if (static_cast<int>(result.size()) >= max_angles) break;
    }
    if (result.empty()) result.push_back(scored.front().angle);

    const auto& best = scored.front();
    std::cerr << "[HausdorffMBR] candidates=" << scored.size()
        << " best_angle_deg=" << best.angle * 180.0 / M_PI
        << " mean=" << best.mean_distance
        << " q90=" << best.q90_distance
        << " covered=" << best.covered_ratio
        << " score=" << best.score << std::endl;
    return result;
}

// ===== distancePointToSegmentWithProjection =====
// 作用：计算点到线段的距离，同时输出投影参数 t（0=起点, 1=终点）
// 参数：point - 查询点; seg_start,seg_end - 线段两端; t - 投影比例输出
// 返回：点到线段的最短欧氏距离
double distancePointToSegmentWithProjection(const pcl::PointXYZ& point,
    const pcl::PointXYZ& seg_start,
    const pcl::PointXYZ& seg_end,
    double& t)
{
    double dx = seg_end.x - seg_start.x;
    double dy = seg_end.y - seg_start.y;
    double len_sq = dx * dx + dy * dy;
    if (len_sq < 1e-10) {
        t = 0.0;
        return std::hypot(point.x - seg_start.x, point.y - seg_start.y);
    }

    t = ((point.x - seg_start.x) * dx + (point.y - seg_start.y) * dy) / len_sq;
    double tc = clampDouble(t, 0.0, 1.0);
    double px = seg_start.x + tc * dx;
    double py = seg_start.y + tc * dy;
    return std::hypot(point.x - px, point.y - py);
}

double distanceXYToSegmentWithProjection(
    double px,
    double py,
    const pcl::PointXYZ& seg_start,
    const pcl::PointXYZ& seg_end,
    double& t)
{
    const double ax = static_cast<double>(seg_start.x);
    const double ay = static_cast<double>(seg_start.y);
    const double bx = static_cast<double>(seg_end.x);
    const double by = static_cast<double>(seg_end.y);
    const double dx = bx - ax;
    const double dy = by - ay;
    const double len_sq = dx * dx + dy * dy;
    if (len_sq < 1e-10) {
        t = 0.0;
        return std::hypot(px - ax, py - ay);
    }

    t = ((px - ax) * dx + (py - ay) * dy) / len_sq;
    const double tc = clampDouble(t, 0.0, 1.0);
    const double qx = ax + tc * dx;
    const double qy = ay + tc * dy;
    return std::hypot(px - qx, py - qy);
}

double distanceXYToSegment2D(
    double px,
    double py,
    const pcl::PointXYZ& seg_start,
    const pcl::PointXYZ& seg_end)
{
    double t = 0.0;
    return distanceXYToSegmentWithProjection(px, py, seg_start, seg_end, t);
}

// ===== edgeSupportPenalty2D =====
// 作用：评估多边形各边对支撑点云的覆盖程度，对缺乏点云支撑的边施加惩罚
// 参数：pts - 多边形顶点; support_cloud - 支撑点云; support_distance - 关联距离阈值; area_hint - 面积提示（用于归一化）
// 返回：惩罚能量值，越大表示多边形与点云的偏离越严重
double edgeSupportPenalty2D(const std::vector<pcl::PointXYZ>& pts,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support_cloud,
    double support_distance,
    double area_hint)
{
    if (!support_cloud || support_cloud->empty() || pts.size() < 3) return 0.0;

    double scale = std::sqrt(std::max(area_hint, 1.0));
    double penalty = 0.0;
    const double min_coverage = 0.35;

    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& a = pts[i];
        const auto& b = pts[(i + 1) % pts.size()];
        double len = std::hypot(b.x - a.x, b.y - a.y);
        if (len < support_distance * 2.0) continue;

        int bin_count = static_cast<int>(clampDouble(len / std::max(support_distance, 0.1), 6.0, 30.0));
        std::vector<bool> covered(bin_count, false);
        int support_count = 0;

        for (const auto& p : support_cloud->points) {
            double t = 0.0;
            double dist = distancePointToSegmentWithProjection(p, a, b, t);
            if (dist <= support_distance && t >= -0.05 && t <= 1.05) {
                int bin = static_cast<int>(clampDouble(t, 0.0, 0.999999) * bin_count);
                covered[bin] = true;
                ++support_count;
            }
        }

        int covered_count = 0;
        for (bool has_support : covered) {
            if (has_support) ++covered_count;
        }

        double coverage = static_cast<double>(covered_count) / std::max(1, bin_count);
        double expected_support = std::max(2.0, len / std::max(support_distance * 1.5, 0.1));
        double support_ratio = clampDouble(support_count / expected_support, 0.0, 1.0);
        double confidence = 0.5 * coverage + 0.5 * support_ratio;

        if (confidence < min_coverage) {
            double normalized = (min_coverage - confidence) / min_coverage;
            penalty += (len / scale) * normalized * normalized;
        }
    }

    return penalty;
}

// ===== angleDistanceToOrthogonalSystem =====
// 作用：计算线角度与基准角度在正交系统（0/90度）下的角度偏差
// 参数：line_angle - 待评估线角度; base_line_angle - 基准线角度
// 返回：到最近的 0° 或 90° 倍数的偏差（弧度）
double angleDistanceToOrthogonalSystem(double line_angle, double base_line_angle)
{
    double diff = std::abs(normalizeAnglePi(line_angle - base_line_angle));
    diff = std::fmod(diff, M_PI);
    if (diff > M_PI / 2.0) diff = M_PI - diff;
    return std::min(diff, std::abs(M_PI / 2.0 - diff));
}

// ===== vertexAngleRad =====
// 作用：计算顶点 cur 处的内角（弧度），由 prev-cur-next 三点确定
// 返回：内角值 [0, PI]；退化情况返回 0
double vertexAngleRad(const pcl::PointXYZ& prev, const pcl::PointXYZ& cur, const pcl::PointXYZ& next)
{
    Eigen::Vector2d a(prev.x - cur.x, prev.y - cur.y);
    Eigen::Vector2d b(next.x - cur.x, next.y - cur.y);
    double na = a.norm();
    double nb = b.norm();
    if (na < 1e-8 || nb < 1e-8) return 0.0;

    double dot = a.dot(b) / (na * nb);
    dot = clampDouble(dot, -1.0, 1.0);
    return std::acos(dot);
}

// ===== regularityPenalty2D =====
// 作用：综合评估多边形的规则性惩罚（角度偏差 + 短边 + 极端转角）
// 参数：pts - 多边形顶点; area_hint - 面积提示（用于尺度归一化）
// 返回：规则性惩罚能量，越规整越接近 0
double regularityPenalty2D(const std::vector<pcl::PointXYZ>& pts, double area_hint)
{
    if (pts.size() < 3) return 0.0;

    double scale = std::sqrt(std::max(area_hint, 1.0));
    double min_edge = clampDouble(0.015 * scale, 0.3, 3.0);
    std::vector<double> base_angles = dominantLineAngles2D(pts, 3);
    const double free_angle = 5.0 * M_PI / 180.0;
    const double soft_angle = 15.0 * M_PI / 180.0;

    double angle_penalty = 0.0;
    double short_edge_penalty = 0.0;
    double corner_penalty = 0.0;

    for (size_t i = 0; i < pts.size(); ++i) {
        const auto& p = pts[i];
        const auto& q = pts[(i + 1) % pts.size()];
        double dx = q.x - p.x;
        double dy = q.y - p.y;
        double len = std::hypot(dx, dy);
        if (len < 1e-6) continue;

        double angle_diff = std::numeric_limits<double>::max();
        for (double base_angle : base_angles) {
            angle_diff = std::min(angle_diff, angleDistanceToOrthogonalSystem(std::atan2(dy, dx), base_angle));
        }
        if (angle_diff > free_angle) {
            double normalized = (angle_diff - free_angle) / soft_angle;
            angle_penalty += (len / scale) * normalized * normalized;
        }

        if (len < min_edge) {
            double normalized = (min_edge - len) / min_edge;
            short_edge_penalty += normalized * normalized;
        }

        double angle = vertexAngleRad(pts[(i + pts.size() - 1) % pts.size()], pts[i], pts[(i + 1) % pts.size()]);
        const double min_corner = 45.0 * M_PI / 180.0;
        const double collinear_corner = 175.0 * M_PI / 180.0;
        if (angle < min_corner) {
            double normalized = (min_corner - angle) / min_corner;
            corner_penalty += normalized * normalized;
        }
        else if (angle > collinear_corner) {
            double normalized = (angle - collinear_corner) / (M_PI - collinear_corner);
            corner_penalty += 0.5 * normalized * normalized;
        }
    }

    return angle_penalty + 2.0 * short_edge_penalty + 0.5 * corner_penalty;
}

// 规则化参数集合：根据分辨率和面积自适应生成的一组调优阈值
struct OutlineTuning {
    double resolution;          // 点云分辨率（平均点间距）
    double area;                // 参考面积
    double association_distance; // 点到边的数据关联距离阈值
    double prune_distance;      // 短边剔除阈值
    double fine_prune_distance; // 精细短边剔除阈值（迭代后处理）
    double repair_distance;     // 拓扑修复距离阈值
    double huber_delta;         // Ceres Huber 鲁棒核参数
    double mid_anchor_weight;   // 中点锚定约束权重
    double angle_tolerance;     // 边分类角度容差（弧度）
};

// ===== makeOutlineTuning =====
// 作用：根据点云分辨率和参考面积，自适应计算一组规则化调优参数
// 返回：填充好的 OutlineTuning 结构体
OutlineTuning makeOutlineTuning(double resolution, double area)
{
    double safe_area = std::max(area, 1.0);
    double scale = std::sqrt(safe_area);
    if (!std::isfinite(resolution) || resolution <= 1e-6) {
        resolution = clampDouble(0.005 * scale, 0.05, 0.5);
    }

    OutlineTuning tuning;
    tuning.resolution = resolution;
    tuning.area = safe_area;
    tuning.association_distance = clampDouble(std::max(2.5 * resolution, 0.012 * scale), 0.2, 2.0);
    tuning.prune_distance = clampDouble(std::max(3.0 * resolution, 0.008 * scale), 0.2, 2.0);
    tuning.fine_prune_distance = clampDouble(std::max(1.5 * resolution, 0.004 * scale), 0.1, 1.0);
    // Topology artifacts are governed more by the footprint scale than by the
    // sampled point spacing.  The previous 0.028 * scale term left almost all
    // medium-sized buildings pinned to the 1 m lower bound.
    tuning.repair_distance = clampDouble(std::max(8.0 * resolution, 0.10 * scale), 0.8, 4.0);
    tuning.huber_delta = clampDouble(std::max(1.5 * resolution, 0.05), 0.05, 0.5);
    tuning.mid_anchor_weight = clampDouble(0.5 / std::max(tuning.huber_delta, 0.05), 2.0, 8.0);
    tuning.angle_tolerance = 25.0 * M_PI / 180.0;
    return tuning;
}

double undirectedAngleDifference(double a, double b)
{
    double diff = std::fmod(std::abs(a - b), M_PI);
    return std::min(diff, M_PI - diff);
}

std::vector<pcl::PointXYZ> robustEdgeInliers(
    const std::vector<pcl::PointXYZ>& candidates,
    const pcl::PointXYZ& edge_start,
    const pcl::PointXYZ& edge_end,
    const OutlineTuning& tuning)
{
    if (candidates.size() < 8) return candidates;

    std::vector<pcl::PointXYZ> active = candidates;
    const double edge_angle = std::atan2(edge_end.y - edge_start.y, edge_end.x - edge_start.x);
    const double max_angle_error = 20.0 * M_PI / 180.0;

    // Iteratively fit a 2D TLS line and reject points outside a robust MAD band.
    // The direction gate prevents a nearby facade from replacing the mask edge.
    for (int iteration = 0; iteration < 2; ++iteration) {
        double cx = 0.0, cy = 0.0;
        for (const auto& p : active) {
            cx += p.x;
            cy += p.y;
        }
        cx /= active.size();
        cy /= active.size();

        double xx = 0.0, xy = 0.0, yy = 0.0;
        for (const auto& p : active) {
            const double dx = p.x - cx;
            const double dy = p.y - cy;
            xx += dx * dx;
            xy += dx * dy;
            yy += dy * dy;
        }
        const double fitted_angle = 0.5 * std::atan2(2.0 * xy, xx - yy);
        if (undirectedAngleDifference(fitted_angle, edge_angle) > max_angle_error) {
            return candidates;
        }

        const double nx = -std::sin(fitted_angle);
        const double ny = std::cos(fitted_angle);
        std::vector<double> distances;
        distances.reserve(active.size());
        for (const auto& p : active) {
            distances.push_back(std::abs((p.x - cx) * nx + (p.y - cy) * ny));
        }
        std::vector<double> sorted = distances;
        std::sort(sorted.begin(), sorted.end());
        const double median = sorted[sorted.size() / 2];
        for (double& value : sorted) value = std::abs(value - median);
        std::sort(sorted.begin(), sorted.end());
        const double mad = sorted[sorted.size() / 2];
        const double threshold = clampDouble(
            median + 2.5 * std::max(1.4826 * mad, 0.03),
            0.08, tuning.association_distance);

        std::vector<pcl::PointXYZ> filtered;
        filtered.reserve(active.size());
        for (size_t i = 0; i < active.size(); ++i) {
            if (distances[i] <= threshold) filtered.push_back(active[i]);
        }
        if (filtered.size() < 6 || filtered.size() * 2 < candidates.size()) return candidates;
        active.swap(filtered);
    }

    // Keep comparable residual counts per edge and sample across its full span.
    constexpr size_t kMaxResidualsPerEdge = 120;
    if (active.size() > kMaxResidualsPerEdge) {
        const double dx = edge_end.x - edge_start.x;
        const double dy = edge_end.y - edge_start.y;
        std::sort(active.begin(), active.end(), [&](const auto& lhs, const auto& rhs) {
            const double lt = (lhs.x - edge_start.x) * dx + (lhs.y - edge_start.y) * dy;
            const double rt = (rhs.x - edge_start.x) * dx + (rhs.y - edge_start.y) * dy;
            return lt < rt;
        });
        std::vector<pcl::PointXYZ> sampled;
        sampled.reserve(kMaxResidualsPerEdge);
        for (size_t i = 0; i < kMaxResidualsPerEdge; ++i) {
            const size_t index = i * (active.size() - 1) / (kMaxResidualsPerEdge - 1);
            sampled.push_back(active[index]);
        }
        active.swap(sampled);
    }
    return active;
}

void robustEdgeInliersWithWeights(
    std::vector<pcl::PointXYZ>& candidates,
    std::vector<double>& weights,
    const pcl::PointXYZ& edge_start,
    const pcl::PointXYZ& edge_end,
    const OutlineTuning& tuning)
{
    if (weights.size() != candidates.size()) {
        candidates = robustEdgeInliers(candidates, edge_start, edge_end, tuning);
        weights.assign(candidates.size(), 1.0);
        return;
    }
    if (candidates.size() < 8) return;

    std::vector<pcl::PointXYZ> active = candidates;
    std::vector<double> active_weights = weights;
    const double edge_angle = std::atan2(edge_end.y - edge_start.y, edge_end.x - edge_start.x);
    const double max_angle_error = 20.0 * M_PI / 180.0;

    for (int iteration = 0; iteration < 2; ++iteration) {
        double cx = 0.0, cy = 0.0;
        for (const auto& p : active) {
            cx += p.x;
            cy += p.y;
        }
        cx /= active.size();
        cy /= active.size();

        double xx = 0.0, xy = 0.0, yy = 0.0;
        for (const auto& p : active) {
            const double dx = p.x - cx;
            const double dy = p.y - cy;
            xx += dx * dx;
            xy += dx * dy;
            yy += dy * dy;
        }
        const double fitted_angle = 0.5 * std::atan2(2.0 * xy, xx - yy);
        if (undirectedAngleDifference(fitted_angle, edge_angle) > max_angle_error) {
            return;
        }

        const double nx = -std::sin(fitted_angle);
        const double ny = std::cos(fitted_angle);
        std::vector<double> distances;
        distances.reserve(active.size());
        for (const auto& p : active) {
            distances.push_back(std::abs((p.x - cx) * nx + (p.y - cy) * ny));
        }
        std::vector<double> sorted = distances;
        std::sort(sorted.begin(), sorted.end());
        const double median = sorted[sorted.size() / 2];
        for (double& value : sorted) value = std::abs(value - median);
        std::sort(sorted.begin(), sorted.end());
        const double mad = sorted[sorted.size() / 2];
        const double threshold = clampDouble(
            median + 2.5 * std::max(1.4826 * mad, 0.03),
            0.08, tuning.association_distance);

        std::vector<pcl::PointXYZ> filtered;
        std::vector<double> filtered_weights;
        filtered.reserve(active.size());
        filtered_weights.reserve(active_weights.size());
        for (size_t i = 0; i < active.size(); ++i) {
            if (distances[i] <= threshold) {
                filtered.push_back(active[i]);
                filtered_weights.push_back(active_weights[i]);
            }
        }
        if (filtered.size() < 6 || filtered.size() * 2 < candidates.size()) return;
        active.swap(filtered);
        active_weights.swap(filtered_weights);
    }

    constexpr size_t kMaxResidualsPerEdge = 120;
    if (active.size() > kMaxResidualsPerEdge) {
        const double dx = edge_end.x - edge_start.x;
        const double dy = edge_end.y - edge_start.y;
        std::vector<size_t> order(active.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t lhs, size_t rhs) {
            const auto& lp = active[lhs];
            const auto& rp = active[rhs];
            const double lt = (lp.x - edge_start.x) * dx + (lp.y - edge_start.y) * dy;
            const double rt = (rp.x - edge_start.x) * dx + (rp.y - edge_start.y) * dy;
            return lt < rt;
        });

        std::vector<pcl::PointXYZ> sampled;
        std::vector<double> sampled_weights;
        sampled.reserve(kMaxResidualsPerEdge);
        sampled_weights.reserve(kMaxResidualsPerEdge);
        for (size_t i = 0; i < kMaxResidualsPerEdge; ++i) {
            const size_t index = order[i * (active.size() - 1) / (kMaxResidualsPerEdge - 1)];
            sampled.push_back(active[index]);
            sampled_weights.push_back(active_weights[index]);
        }
        active.swap(sampled);
        active_weights.swap(sampled_weights);
    }

    candidates.swap(active);
    weights.swap(active_weights);
}

struct DirectionSystem {
    double angle = 0.0;         // Line angle folded into [0, PI/2).
    double weight = 0.0;        // Geometry + support weighted evidence.
    double length = 0.0;        // Total contributing edge length.
    double support_score = 0.0; // 0..1, averaged over contributing edges.
    size_t edge_count = 0;
    double longest_edge = 0.0;
};

struct LocalChainEdge {
    size_t index = 0;
    int direction_id = -1;
    bool constrained = false;
    double line_angle = 0.0;
    double snapped_line_angle = 0.0;
    double length = 0.0;
    double angle_error = 0.0;
};

double foldedLineAngle90(double angle)
{
    angle = std::fmod(angle, M_PI / 2.0);
    if (angle < 0.0) angle += M_PI / 2.0;
    return angle;
}

double foldedAngleDistance90(double a, double b)
{
    const double period = M_PI / 2.0;
    double diff = std::fmod(std::abs(a - b), period);
    return std::min(diff, period - diff);
}

double lineAngleToNearestOrthogonalAxis(double line_angle, double base_angle)
{
    const double normal_line = normalizeAnglePi(line_angle);
    const double axis0 = normalizeAnglePi(base_angle);
    const double axis1 = normalizeAnglePi(base_angle + M_PI / 2.0);
    return undirectedAngleDifference(normal_line, axis0) <=
        undirectedAngleDifference(normal_line, axis1) ? axis0 : axis1;
}

Eigen::Vector3d lineFromPointAndAngle(const pcl::PointXYZ& point, double line_angle)
{
    const double normal_angle = line_angle - M_PI / 2.0;
    const double a = std::cos(normal_angle);
    const double b = std::sin(normal_angle);
    const double d = point.x * a + point.y * b;
    return Eigen::Vector3d(a, b, -d);
}

pcl::PointXYZ intersectionOfLines2D(
    const Eigen::Vector3d& lhs,
    const Eigen::Vector3d& rhs,
    float z_value)
{
    const double det = lhs.x() * rhs.y() - rhs.x() * lhs.y();
    pcl::PointXYZ point;
    point.z = z_value;
    if (std::abs(det) < 1e-8) {
        point.x = std::numeric_limits<float>::quiet_NaN();
        point.y = std::numeric_limits<float>::quiet_NaN();
        return point;
    }

    point.x = static_cast<float>((lhs.y() * (-rhs.z()) - rhs.y() * (-lhs.z())) / det);
    point.y = static_cast<float>((rhs.x() * (-lhs.z()) - lhs.x() * (-rhs.z())) / det);
    return point;
}

double circularMeanAngle90(const std::vector<DirectionSystem>& systems)
{
    double x = 0.0;
    double y = 0.0;
    double total = 0.0;
    for (const auto& item : systems) {
        x += item.weight * std::cos(4.0 * item.angle);
        y += item.weight * std::sin(4.0 * item.angle);
        total += item.weight;
    }
    if (total <= 1e-9) return systems.empty() ? 0.0 : systems.front().angle;
    double mean = 0.25 * std::atan2(y, x);
    return foldedLineAngle90(mean);
}

std::vector<DirectionSystem> detectDirectionSystems(
    const std::vector<pcl::PointXYZ>& polygon,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support,
    const OutlineTuning& tuning)
{
    std::vector<DirectionSystem> evidence;
    if (polygon.size() < 3) return evidence;

    const double area = std::max(polygonArea2D(polygon), 1.0);
    const double scale = std::sqrt(area);
    const double min_length = clampDouble(0.025 * scale, 0.8, 2.5);
    const double max_median = std::min(0.45, tuning.association_distance * 0.9);

    double total_length = 0.0;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double length = std::hypot(dx, dy);
        if (length < min_length) continue;
        total_length += length;

        double support_score = 0.0;
        if (support && !support->empty()) {
            std::vector<bool> covered(8, false);
            size_t inliers = 0;
            std::vector<double> distances;
            distances.reserve(32);
            for (const auto& p : support->points) {
                const double t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / (length * length);
                if (t <= 0.05 || t >= 0.95) continue;
                const double distance = std::abs((p.x - a.x) * dy - (p.y - a.y) * dx) / length;
                if (distance <= tuning.association_distance) {
                    ++inliers;
                    covered[std::min<size_t>(7, static_cast<size_t>(t * 8.0))] = true;
                    distances.push_back(distance);
                }
            }

            const size_t covered_bins = static_cast<size_t>(
                std::count(covered.begin(), covered.end(), true));
            if (!distances.empty()) {
                std::sort(distances.begin(), distances.end());
                const double median_distance = distances[distances.size() / 2];
                const double coverage_score = static_cast<double>(covered_bins) / covered.size();
                const double inlier_score = clampDouble(
                    static_cast<double>(inliers) / std::max(10.0, length / std::max(tuning.association_distance, 0.1)),
                    0.0, 1.0);
                const double fit_score = clampDouble(
                    1.0 - median_distance / std::max(max_median, 1e-6),
                    0.0, 1.0);
                support_score = 0.45 * coverage_score + 0.30 * inlier_score + 0.25 * fit_score;
            }
        }

        DirectionSystem item;
        item.angle = foldedLineAngle90(std::atan2(dy, dx));
        item.length = length;
        item.weight = length * (0.65 + 0.35 * support_score);
        item.support_score = support_score;
        item.edge_count = 1;
        item.longest_edge = length;
        evidence.push_back(item);
    }

    if (evidence.empty()) return evidence;

    std::sort(evidence.begin(), evidence.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.weight > rhs.weight;
    });

    const double cluster_tolerance = 8.0 * M_PI / 180.0;
    std::vector<DirectionSystem> clusters;
    for (const auto& item : evidence) {
        int best_index = -1;
        double best_distance = cluster_tolerance;
        for (size_t i = 0; i < clusters.size(); ++i) {
            const double distance = foldedAngleDistance90(item.angle, clusters[i].angle);
            if (distance < best_distance) {
                best_distance = distance;
                best_index = static_cast<int>(i);
            }
        }

        if (best_index < 0) {
            clusters.push_back(item);
        }
        else {
            auto& cluster = clusters[best_index];
            std::vector<DirectionSystem> for_mean;
            for_mean.push_back(cluster);
            for_mean.push_back(item);
            const double old_weight = cluster.weight;
            cluster.angle = circularMeanAngle90(for_mean);
            cluster.weight += item.weight;
            cluster.length += item.length;
            cluster.support_score =
                (cluster.support_score * old_weight + item.support_score * item.weight) /
                std::max(cluster.weight, 1e-9);
            cluster.edge_count += item.edge_count;
            cluster.longest_edge = std::max(cluster.longest_edge, item.longest_edge);
        }
    }

    std::sort(clusters.begin(), clusters.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.weight > rhs.weight;
    });

    const double total_weight = [&]() {
        double sum = 0.0;
        for (const auto& item : clusters) sum += item.weight;
        return std::max(sum, 1e-9);
    }();

    std::vector<DirectionSystem> accepted;
    accepted.reserve(std::min<size_t>(clusters.size(), 4));
    for (size_t i = 0; i < clusters.size(); ++i) {
        const double ratio = clusters[i].weight / total_weight;
        const bool is_primary = accepted.empty();
        const double primary_ratio = accepted.empty()
            ? ratio
            : accepted.front().weight / total_weight;
        const double long_edge_threshold =
            clampDouble(0.12 * total_length, 8.0, 25.0);
        const bool strong_secondary =
            ratio >= 0.22 &&
            clusters[i].edge_count >= 4 &&
            clusters[i].support_score >= 0.40;
        const bool long_secondary =
            ratio >= 0.20 &&
            clusters[i].edge_count >= 3 &&
            clusters[i].longest_edge >= long_edge_threshold &&
            clusters[i].support_score >= 0.45;
        const bool very_long_single =
            ratio >= 0.28 &&
            clusters[i].longest_edge >= 0.24 * total_length &&
            clusters[i].support_score >= 0.55;
        const bool dominant_primary_guard =
            primary_ratio >= 0.62 &&
            ratio < 0.22 &&
            clusters[i].longest_edge < clampDouble(0.18 * total_length, 14.0, 35.0);
        const bool accept_secondary =
            !dominant_primary_guard &&
            (strong_secondary || long_secondary || very_long_single);

        if (is_primary || accept_secondary) {
            bool duplicate = false;
            for (const auto& kept : accepted) {
                if (foldedAngleDistance90(clusters[i].angle, kept.angle) < 12.0 * M_PI / 180.0) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) accepted.push_back(clusters[i]);
        }
        if (accepted.size() >= 4) break;
    }

    std::cerr << "[DirectionSystems] candidates=" << clusters.size()
        << " accepted=" << accepted.size();
    for (const auto& item : accepted) {
        std::cerr << " {angle_deg=" << item.angle * 180.0 / M_PI
            << ", ratio=" << item.weight / total_weight
            << ", edges=" << item.edge_count
            << ", longest=" << item.longest_edge
            << ", support=" << item.support_score << "}";
    }
    std::cerr << std::endl;

    return accepted;
}

std::vector<DirectionSystem> buildMultiDirectionCandidates(
    const std::vector<pcl::PointXYZ>& polygon,
    double resolution)
{
    std::vector<DirectionSystem> result;
    if (polygon.size() < 4) return result;

    const double area = std::max(polygonArea2D(polygon), 1.0);
    const double scale = std::sqrt(area);
    const double min_length = clampDouble(std::max(2.5 * resolution, 0.018 * scale), 1.0, 4.0);
    const double cluster_tolerance = 10.0 * M_PI / 180.0;

    std::vector<DirectionSystem> clusters;
    double total_length = 0.0;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double length = std::hypot(static_cast<double>(b.x) - a.x,
            static_cast<double>(b.y) - a.y);
        if (length < min_length) continue;
        total_length += length;
        const double angle = foldedLineAngle90(std::atan2(
            static_cast<double>(b.y) - a.y,
            static_cast<double>(b.x) - a.x));

        int best_index = -1;
        double best_distance = cluster_tolerance;
        for (size_t k = 0; k < clusters.size(); ++k) {
            const double distance = foldedAngleDistance90(angle, clusters[k].angle);
            if (distance < best_distance) {
                best_distance = distance;
                best_index = static_cast<int>(k);
            }
        }

        DirectionSystem item;
        item.angle = angle;
        item.weight = length;
        item.length = length;
        item.support_score = 0.65;
        item.edge_count = 1;
        item.longest_edge = length;
        if (best_index < 0) {
            clusters.push_back(item);
        } else {
            auto& cluster = clusters[static_cast<size_t>(best_index)];
            std::vector<DirectionSystem> mean_items = { cluster, item };
            cluster.angle = circularMeanAngle90(mean_items);
            cluster.weight += item.weight;
            cluster.length += item.length;
            cluster.edge_count += 1;
            cluster.longest_edge = std::max(cluster.longest_edge, length);
        }
    }
    if (clusters.empty()) return result;

    std::sort(clusters.begin(), clusters.end(), [](const auto& a, const auto& b) {
        return a.length > b.length;
    });

    auto append_system = [&](DirectionSystem system) {
        system.angle = foldedLineAngle90(system.angle);
        for (const auto& kept : result) {
            if (foldedAngleDistance90(kept.angle, system.angle) < 12.0 * M_PI / 180.0) {
                return;
            }
        }
        result.push_back(system);
    };

    std::vector<double> mbr_angles = hausdorffMbrLineAngles2D(polygon, resolution, 1);
    if (!mbr_angles.empty()) {
        DirectionSystem primary = clusters.front();
        primary.angle = mbr_angles.front();
        primary.support_score = std::max(primary.support_score, 0.65);
        append_system(primary);
    }

    for (const auto& cluster : clusters) {
        const double ratio = cluster.length / std::max(total_length, 1e-9);
        if (result.empty() || ratio >= 0.08 || cluster.longest_edge >= clampDouble(0.10 * total_length, 6.0, 25.0)) {
            DirectionSystem system = cluster;
            system.support_score = std::max(system.support_score, 0.65);
            append_system(system);
        }
        if (result.size() >= 4) break;
    }

    std::cerr << "[MultiDirectionCandidates] count=" << result.size();
    for (const auto& system : result) {
        std::cerr << " {angle_deg=" << system.angle * 180.0 / M_PI
            << ", len=" << system.length
            << ", edges=" << system.edge_count
            << ", longest=" << system.longest_edge << "}";
    }
    std::cerr << std::endl;
    return result;
}

std::vector<LocalChainEdge> classifyLocalChainEdges(
    const std::vector<pcl::PointXYZ>& polygon,
    const std::vector<DirectionSystem>& systems,
    const OutlineTuning& tuning)
{
    std::vector<LocalChainEdge> edges;
    if (polygon.size() < 3 || systems.empty()) return edges;

    const double scale = std::sqrt(std::max(polygonArea2D(polygon), 1.0));
    const double min_constrained_length = clampDouble(0.018 * scale, 0.9, 4.0);
    const double angle_tolerance = 14.0 * M_PI / 180.0;
    edges.reserve(polygon.size());

    for (size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        LocalChainEdge edge;
        edge.index = i;
        edge.length = std::hypot(dx, dy);
        edge.line_angle = std::atan2(dy, dx);

        double best_error = std::numeric_limits<double>::max();
        int best_id = -1;
        double best_snapped = edge.line_angle;
        for (size_t k = 0; k < systems.size(); ++k) {
            const double snapped = lineAngleToNearestOrthogonalAxis(edge.line_angle, systems[k].angle);
            const double error = undirectedAngleDifference(edge.line_angle, snapped);
            if (error < best_error) {
                best_error = error;
                best_id = static_cast<int>(k);
                best_snapped = snapped;
            }
        }

        edge.direction_id = best_id;
        edge.angle_error = best_error;
        edge.snapped_line_angle = best_snapped;
        edge.constrained =
            best_id >= 0 &&
            edge.length >= min_constrained_length &&
            best_error <= angle_tolerance;
        edges.push_back(edge);
    }

    return edges;
}

bool hasCredibleMultiDirectionChains(
    const std::vector<pcl::PointXYZ>& polygon,
    const std::vector<DirectionSystem>& systems,
    const OutlineTuning& tuning)
{
    if (polygon.size() < 5 || systems.size() < 2) return false;

    const std::vector<LocalChainEdge> edges =
        classifyLocalChainEdges(polygon, systems, tuning);
    if (edges.size() != polygon.size()) return false;

    const double perimeter = [&]() {
        double sum = 0.0;
        for (const auto& edge : edges) sum += edge.length;
        return std::max(sum, 1e-9);
    }();
    const double scale = std::sqrt(std::max(polygonArea2D(polygon), 1.0));
    const double min_chain_length = clampDouble(0.10 * perimeter, 6.0, 25.0);

    std::vector<double> direction_length(systems.size(), 0.0);
    std::vector<size_t> direction_edges(systems.size(), 0);
    for (const auto& edge : edges) {
        if (!edge.constrained || edge.direction_id < 0 ||
            static_cast<size_t>(edge.direction_id) >= systems.size()) {
            continue;
        }
        direction_length[static_cast<size_t>(edge.direction_id)] += edge.length;
        ++direction_edges[static_cast<size_t>(edge.direction_id)];
    }

    std::vector<double> longest_chain(systems.size(), 0.0);
    std::vector<size_t> chain_edges(systems.size(), 0);
    std::vector<size_t> max_chain_edges(systems.size(), 0);
    const size_t n = edges.size();
    size_t start = 0;
    while (start < n && !edges[start].constrained) ++start;
    if (start == n) {
        std::cerr << "[MultiChain] rejected: no constrained edges" << std::endl;
        return false;
    }

    // Linearize the ring from the first constrained edge and merge same-direction runs.
    for (size_t offset = 0; offset < n; ) {
        const size_t idx = (start + offset) % n;
        const auto& edge = edges[idx];
        if (!edge.constrained || edge.direction_id < 0 ||
            static_cast<size_t>(edge.direction_id) >= systems.size()) {
            ++offset;
            continue;
        }
        const int id = edge.direction_id;
        double length = 0.0;
        size_t count = 0;
        while (offset < n) {
            const size_t j = (start + offset) % n;
            const auto& e = edges[j];
            if (!e.constrained || e.direction_id != id) break;
            length += e.length;
            ++count;
            ++offset;
        }
        const size_t sid = static_cast<size_t>(id);
        if (length > longest_chain[sid]) {
            longest_chain[sid] = length;
            max_chain_edges[sid] = count;
        }
        chain_edges[sid] += count;
    }

    bool has_secondary = false;
    for (size_t i = 1; i < systems.size(); ++i) {
        const double ratio = direction_length[i] / perimeter;
        const double chain_ratio = longest_chain[i] / perimeter;
        const bool enough_total = ratio >= 0.16 && direction_edges[i] >= 3;
        // 连续同向链通道：阶梯状斜翼表现为若干条同向短边连成一串。
        const bool contiguous_chain =
            longest_chain[i] >= min_chain_length &&
            chain_ratio >= 0.10 &&
            max_chain_edges[i] >= 2;
        // 单长边翼通道：VDP 简化后的斜翼长边常是一条完整长边(如 40m+)，
        // 与相邻边之间隔着正交连接边，永远凑不出"连续同向链"。此时以
        // "系统最长边自身达到链长门槛"作为等价证据。防阶梯伪影由外层
        // enough_total(总量+边数)与 enough_support(方向支撑)把关，伪影
        // 边短(通常 <3m)且支撑弱，不会从此通道漏入。
        const bool single_long_edge_chain =
            systems[i].longest_edge >= min_chain_length &&
            systems[i].longest_edge >= 0.07 * perimeter;
        const bool enough_chain = contiguous_chain || single_long_edge_chain;
        const bool enough_support = systems[i].support_score >= 0.45;
        std::cerr << "[MultiChain] dir=" << i
            << " angle_deg=" << systems[i].angle * 180.0 / M_PI
            << " len_ratio=" << ratio
            << " longest_chain=" << longest_chain[i]
            << " chain_ratio=" << chain_ratio
            << " edges=" << direction_edges[i]
            << " chain_edges=" << max_chain_edges[i]
            << " longest_edge=" << systems[i].longest_edge
            << " support=" << systems[i].support_score
            << " channel=" << (contiguous_chain ? "contiguous"
                : (single_long_edge_chain ? "single_long_edge" : "none"))
            << std::endl;
        if (enough_total && enough_chain && enough_support) {
            has_secondary = true;
        }
    }

    const bool primary_ok =
        direction_length[0] / perimeter >= 0.20 &&
        longest_chain[0] >= clampDouble(0.08 * perimeter, 5.0, 25.0);
    std::cerr << "[MultiChain] primary_len_ratio=" << direction_length[0] / perimeter
        << " primary_longest=" << longest_chain[0]
        << " scale=" << scale
        << " accepted=" << ((primary_ok && has_secondary) ? 1 : 0)
        << std::endl;
    return primary_ok && has_secondary;
}

bool regularizeByLocalChains(
    const std::vector<pcl::PointXYZ>& input,
    const std::vector<DirectionSystem>& systems,
    const OutlineTuning& tuning,
    std::vector<pcl::PointXYZ>& output)
{
    output.clear();
    if (input.size() < 4 || systems.size() < 2 || !isSimplePolygon2D(input)) {
        return false;
    }

    const std::vector<LocalChainEdge> edges =
        classifyLocalChainEdges(input, systems, tuning);
    if (edges.size() != input.size()) return false;

    size_t constrained_edges = 0;
    std::vector<Eigen::Vector3d> snapped_lines(input.size());
    std::vector<bool> has_line(input.size(), false);

    for (const auto& edge : edges) {
        if (!edge.constrained) continue;
        ++constrained_edges;
        const size_t i = edge.index;
        pcl::PointXYZ mid;
        mid.x = static_cast<float>(0.5 * (input[i].x + input[(i + 1) % input.size()].x));
        mid.y = static_cast<float>(0.5 * (input[i].y + input[(i + 1) % input.size()].y));
        mid.z = input[i].z;
        snapped_lines[i] = lineFromPointAndAngle(mid, edge.snapped_line_angle);
        has_line[i] = true;
    }

    if (constrained_edges < std::max<size_t>(4, input.size() / 2)) {
        std::cerr << "[LocalChain] skip: constrained_edges=" << constrained_edges
            << "/" << input.size() << std::endl;
        return false;
    }

    output = input;
    size_t updated_vertices = 0;
    size_t skipped_parallel = 0;
    size_t skipped_move = 0;
    const double max_move = clampDouble(0.10 * std::sqrt(std::max(polygonArea2D(input), 1.0)), 2.0, 8.0);
    for (size_t i = 0; i < input.size(); ++i) {
        const size_t previous_edge = (i + input.size() - 1) % input.size();
        const size_t next_edge = i;
        if (!has_line[previous_edge] || !has_line[next_edge]) continue;

        pcl::PointXYZ point = intersectionOfLines2D(
            snapped_lines[previous_edge], snapped_lines[next_edge], input[i].z);
        if (std::isnan(point.x) || std::isnan(point.y)) {
            ++skipped_parallel;
            continue;
        }
        const double move = std::hypot(point.x - input[i].x, point.y - input[i].y);
        if (move > max_move) {
            ++skipped_move;
            continue;
        }
        output[i] = point;
        ++updated_vertices;
    }

    removeDuplicatePoints2D(output, static_cast<float>(tuning.fine_prune_distance));
    removeNearlyCollinearPoints2D(output, 12.0);
    removeDuplicatePoints2D(output, static_cast<float>(tuning.fine_prune_distance));

    if (updated_vertices < 2 || output.size() < 3 || !isSimplePolygon2D(output)) {
        std::cerr << "[LocalChain] rejected: updated_vertices=" << updated_vertices
            << " vertices=" << output.size()
            << " constrained_edges=" << constrained_edges << "/" << input.size()
            << " skipped_parallel=" << skipped_parallel
            << " skipped_move=" << skipped_move
            << " max_move=" << max_move
            << std::endl;
        return false;
    }

    const double iou = polygonIoU2D(output, input);
    const double area_ratio = polygonArea2D(output) / std::max(polygonArea2D(input), 1e-6);
    const bool accepted = iou >= 0.82 && area_ratio >= 0.80 && area_ratio <= 1.25;
    std::cerr << "[LocalChain] constrained_edges=" << constrained_edges
        << "/" << input.size()
        << " updated_vertices=" << updated_vertices
        << " skipped_parallel=" << skipped_parallel
        << " skipped_move=" << skipped_move
        << " iou=" << iou
        << " area_ratio=" << area_ratio
        << " accepted=" << (accepted ? 1 : 0)
        << std::endl;
    return accepted;
}

bool hasMultipleDirectionSystems(const std::vector<pcl::PointXYZ>& polygon,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support,
    const OutlineTuning& tuning)
{
    const std::vector<DirectionSystem> systems =
        detectDirectionSystems(polygon, support, tuning);
    if (systems.empty()) {
        std::cerr << "[DirectionSystems] reliable_edges=0 mode=StrictOrthogonal\n";
        return false;
    }
    const bool multiple_systems = systems.size() >= 2;
    std::cerr << "[DirectionSystems] mode="
        << (multiple_systems ? "AllowDiagonal" : "StrictOrthogonal")
        << std::endl;
    return multiple_systems;
#if 0
    if (polygon.size() < 3 || !support || support->empty()) return false;
    const double area = std::max(polygonArea2D(polygon), 1.0);
    const double min_length = clampDouble(0.04 * std::sqrt(area), 1.0, 3.0);

    struct DirectionEvidence {
        double angle = 0.0;  // Folded into [0, PI/2).
        double weight = 0.0;
        double length = 0.0;
    };
    std::vector<DirectionEvidence> evidence;
    double reliable_perimeter = 0.0;

    for (size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double length = std::hypot(dx, dy);
        if (length < min_length) continue;

        std::vector<bool> covered(8, false);
        size_t inliers = 0;
        std::vector<double> distances;
        for (const auto& p : support->points) {
            const double t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / (length * length);
            if (t <= 0.05 || t >= 0.95) continue;
            const double distance = std::abs((p.x - a.x) * dy - (p.y - a.y) * dx) / length;
            if (distance <= tuning.association_distance) {
                ++inliers;
                covered[std::min<size_t>(7, static_cast<size_t>(t * 8.0))] = true;
                distances.push_back(distance);
            }
        }
        const size_t covered_bins = static_cast<size_t>(std::count(covered.begin(), covered.end(), true));
        if (inliers < 10 || covered_bins < 5) continue;
        std::sort(distances.begin(), distances.end());
        const double median_distance = distances[distances.size() / 2];
        const double max_median = std::min(0.35, tuning.association_distance * 0.75);
        if (median_distance > max_median) continue;

        double angle = std::fmod(std::atan2(dy, dx), M_PI / 2.0);
        if (angle < 0.0) angle += M_PI / 2.0;
        const double coverage = static_cast<double>(covered_bins) / covered.size();
        const double fit_score = clampDouble(1.0 - median_distance / std::max(max_median, 1e-6), 0.2, 1.0);
        evidence.push_back({ angle, length * coverage * fit_score, length });
        reliable_perimeter += length;
    }

    if (evidence.empty()) {
        std::cerr << "[DirectionSystems] reliable_edges=0 mode=StrictOrthogonal\n";
        return false;
    }

    auto foldedDistance = [](double a, double b) {
        const double period = M_PI / 2.0;
        double diff = std::fmod(std::abs(a - b), period);
        return std::min(diff, period - diff);
    };
    const double cluster_tolerance = 10.0 * M_PI / 180.0;
    auto clusterWeight = [&](double center, const std::vector<bool>* eligible,
        size_t* edge_count, double* longest_edge) {
        double weight = 0.0;
        *edge_count = 0;
        *longest_edge = 0.0;
        for (size_t i = 0; i < evidence.size(); ++i) {
            if (eligible && !(*eligible)[i]) continue;
            if (foldedDistance(evidence[i].angle, center) <= cluster_tolerance) {
                weight += evidence[i].weight;
                ++(*edge_count);
                *longest_edge = std::max(*longest_edge, evidence[i].length);
            }
        }
        return weight;
    };

    double total_weight = 0.0;
    for (const auto& item : evidence) total_weight += item.weight;
    double primary_weight = -1.0;
    double primary_angle = evidence.front().angle;
    for (const auto& candidate : evidence) {
        size_t count = 0;
        double longest = 0.0;
        const double weight = clusterWeight(candidate.angle, nullptr, &count, &longest);
        if (weight > primary_weight) {
            primary_weight = weight;
            primary_angle = candidate.angle;
        }
    }

    std::vector<bool> secondary_eligible(evidence.size(), false);
    for (size_t i = 0; i < evidence.size(); ++i) {
        secondary_eligible[i] = foldedDistance(evidence[i].angle, primary_angle) >=
            15.0 * M_PI / 180.0;
    }
    double secondary_weight = 0.0;
    size_t secondary_edges = 0;
    double secondary_longest = 0.0;
    for (size_t i = 0; i < evidence.size(); ++i) {
        if (!secondary_eligible[i]) continue;
        size_t count = 0;
        double longest = 0.0;
        const double weight = clusterWeight(
            evidence[i].angle, &secondary_eligible, &count, &longest);
        if (weight > secondary_weight) {
            secondary_weight = weight;
            secondary_edges = count;
            secondary_longest = longest;
        }
    }

    const double primary_ratio = primary_weight / std::max(total_weight, 1e-9);
    const double secondary_ratio = secondary_weight / std::max(total_weight, 1e-9);
    const bool secondary_spatial_support = secondary_edges >= 2 ||
        secondary_longest >= 0.22 * std::max(reliable_perimeter, 1e-6);
    const bool multiple_systems = secondary_ratio >= 0.20 && secondary_spatial_support;
    std::cerr << "[DirectionSystems] reliable_edges=" << evidence.size()
        << " primary_ratio=" << primary_ratio
        << " secondary_ratio=" << secondary_ratio
        << " secondary_edges=" << secondary_edges
        << " mode=" << (multiple_systems ? "AllowDiagonal" : "StrictOrthogonal")
        << std::endl;
    return multiple_systems;
#endif
}

struct PreservedArcSegment {
    size_t start = 0;
    size_t end = 0;
    pcl::PointXYZ start_point;
    pcl::PointXYZ end_point;
    bool use_bezier = false;
    double cx = 0.0;
    double cy = 0.0;
    double radius = 0.0;
    double start_angle = 0.0;
    double sweep_angle = 0.0;
    double rmse = 0.0;
    Eigen::Vector2d bezier_p0 = Eigen::Vector2d::Zero();
    Eigen::Vector2d bezier_p1 = Eigen::Vector2d::Zero();
    Eigen::Vector2d bezier_p2 = Eigen::Vector2d::Zero();
    Eigen::Vector2d bezier_p3 = Eigen::Vector2d::Zero();
    double curve_length = 0.0;
};

bool forceOrthogonalPolygonToAngle(
    const std::vector<pcl::PointXYZ>& input,
    double main_angle,
    std::vector<pcl::PointXYZ>& result)
{
    result.clear();
    if (input.size() < 4) return false;

    std::vector<pcl::PointXYZ> working = input;
    auto axisOfEdge = [main_angle](const pcl::PointXYZ& a, const pcl::PointXYZ& b) {
        const double angle = std::atan2(b.y - a.y, b.x - a.x);
        return undirectedAngleDifference(angle, main_angle) <=
            undirectedAngleDifference(angle, main_angle + M_PI / 2.0) ? 0 : 1;
    };

    // A valid orthogonal ring alternates its two axes. Remove vertices whose
    // incoming and outgoing edges would snap to the same axis before solving
    // line intersections; otherwise the two lines are parallel and no corner
    // exists. Odd vertex counts are resolved by the same cyclic pass.
    bool changed = true;
    while (changed && working.size() >= 4) {
        changed = false;
        for (size_t i = 0; i < working.size(); ++i) {
            const size_t previous = (i + working.size() - 1) % working.size();
            const size_t next = (i + 1) % working.size();
            if (axisOfEdge(working[previous], working[i]) ==
                axisOfEdge(working[i], working[next])) {
                working.erase(working.begin() + i);
                changed = true;
                break;
            }
        }
    }
    if (working.size() < 4 || working.size() % 2 != 0) return false;

    std::vector<Eigen::Vector3d> lines;
    lines.reserve(working.size());
    for (size_t i = 0; i < working.size(); ++i) {
        const auto& a = working[i];
        const auto& b = working[(i + 1) % working.size()];
        const double theta = axisOfEdge(a, b) == 0 ? main_angle - M_PI / 2.0 : main_angle;
        const double d = 0.5 * ((a.x + b.x) * std::cos(theta) +
            (a.y + b.y) * std::sin(theta));
        lines.emplace_back(std::cos(theta), std::sin(theta), -d);
    }

    result.reserve(lines.size());
    for (size_t i = 0; i < lines.size(); ++i) {
        const size_t previous = (i + lines.size() - 1) % lines.size();
        const double det = lines[previous].x() * lines[i].y() -
            lines[i].x() * lines[previous].y();
        if (std::abs(det) < 1e-8) {
            result.clear();
            return false;
        }
        const double previous_d = -lines[previous].z();
        const double current_d = -lines[i].z();
        const double x = (previous_d * lines[i].y() -
            current_d * lines[previous].y()) / det;
        const double y = (lines[previous].x() * current_d -
            lines[i].x() * previous_d) / det;
        result.emplace_back(static_cast<float>(x), static_cast<float>(y), working[0].z);
    }
    removeDuplicatePoints2D(result, 0.05f);
    if (result.size() < 4 || result.size() % 2 != 0 || !isSimplePolygon2D(result)) {
        result.clear();
        return false;
    }
    for (size_t i = 0; i < result.size(); ++i) {
        const double angle = std::atan2(result[(i + 1) % result.size()].y - result[i].y,
            result[(i + 1) % result.size()].x - result[i].x);
        const double error = std::min(undirectedAngleDifference(angle, main_angle),
            undirectedAngleDifference(angle, main_angle + M_PI / 2.0));
        if (error > 0.25 * M_PI / 180.0) {
            result.clear();
            return false;
        }
    }
    return true;
}

bool forceOrthogonalPolygon(const std::vector<pcl::PointXYZ>& input,
    std::vector<pcl::PointXYZ>& result)
{
    const auto dominant = dominantLineAngles2D(input, 1);
    if (dominant.empty()) {
        result.clear();
        return false;
    }
    return forceOrthogonalPolygonToAngle(input, dominant.front(), result);
}

bool isStrictOrthogonalToMainAngle(
    const std::vector<pcl::PointXYZ>& polygon,
    double main_angle,
    double tolerance_deg,
    double min_edge_length)
{
    if (polygon.size() < 4) return false;
    const double tolerance = tolerance_deg * M_PI / 180.0;
    double checked_length = 0.0;
    double off_axis_length = 0.0;
    size_t checked_edges = 0;

    for (size_t i = 0; i < polygon.size(); ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % polygon.size()];
        const double length = std::hypot(b.x - a.x, b.y - a.y);
        if (length < min_edge_length) continue;
        const double angle = std::atan2(b.y - a.y, b.x - a.x);
        const double error = std::min(
            undirectedAngleDifference(angle, main_angle),
            undirectedAngleDifference(angle, main_angle + M_PI / 2.0));
        checked_length += length;
        ++checked_edges;
        if (error > tolerance) off_axis_length += length;
    }

    if (checked_edges < 4 || checked_length < 1e-6) return false;
    return off_axis_length / checked_length <= 0.02;
}

bool fitCircle2D(const std::vector<pcl::PointXYZ>& pts,
    size_t start, size_t count,
    double& cx, double& cy, double& radius, double& rmse)
{
    if (pts.size() < 3 || count < 6) return false;

    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_yy = 0.0;
    double sum_xy = 0.0, sum_xxx = 0.0, sum_yyy = 0.0;
    double sum_xyy = 0.0, sum_xxy = 0.0;
    for (size_t k = 0; k < count; ++k) {
        const auto& p = pts[(start + k) % pts.size()];
        const double x = p.x;
        const double y = p.y;
        const double xx = x * x;
        const double yy = y * y;
        sum_x += x;
        sum_y += y;
        sum_xx += xx;
        sum_yy += yy;
        sum_xy += x * y;
        sum_xxx += xx * x;
        sum_yyy += yy * y;
        sum_xyy += x * yy;
        sum_xxy += xx * y;
    }

    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    A << sum_xx, sum_xy, sum_x,
         sum_xy, sum_yy, sum_y,
         sum_x,  sum_y,  static_cast<double>(count);
    b << -(sum_xxx + sum_xyy),
         -(sum_xxy + sum_yyy),
         -(sum_xx + sum_yy);
    if (std::abs(A.determinant()) < 1e-9) return false;

    const Eigen::Vector3d sol = A.ldlt().solve(b);
    if (!sol.allFinite()) return false;
    cx = -0.5 * sol.x();
    cy = -0.5 * sol.y();
    const double r2 = cx * cx + cy * cy - sol.z();
    if (!std::isfinite(r2) || r2 <= 0.0) return false;
    radius = std::sqrt(r2);
    if (!std::isfinite(radius) || radius <= 0.0) return false;

    double err2 = 0.0;
    for (size_t k = 0; k < count; ++k) {
        const auto& p = pts[(start + k) % pts.size()];
        const double d = std::hypot(p.x - cx, p.y - cy);
        const double e = d - radius;
        err2 += e * e;
    }
    rmse = std::sqrt(err2 / static_cast<double>(count));
    return std::isfinite(rmse);
}

double unwrapArcSweep(const std::vector<pcl::PointXYZ>& pts,
    size_t start, size_t count, double cx, double cy, double& start_angle)
{
    if (count < 2) return 0.0;
    start_angle = std::atan2(pts[start % pts.size()].y - cy, pts[start % pts.size()].x - cx);
    double previous = start_angle;
    double unwrapped = start_angle;
    for (size_t k = 1; k < count; ++k) {
        double angle = std::atan2(
            pts[(start + k) % pts.size()].y - cy,
            pts[(start + k) % pts.size()].x - cx);
        while (angle - previous > M_PI) angle -= 2.0 * M_PI;
        while (angle - previous < -M_PI) angle += 2.0 * M_PI;
        unwrapped += angle - previous;
        previous = angle;
    }
    return unwrapped - start_angle;
}

double polylineLength(const std::vector<pcl::PointXYZ>& pts, size_t start, size_t count)
{
    if (pts.empty() || count < 2) return 0.0;
    double length = 0.0;
    for (size_t k = 1; k < count; ++k) {
        length += std::hypot(
            pts[(start + k) % pts.size()].x - pts[(start + k - 1) % pts.size()].x,
            pts[(start + k) % pts.size()].y - pts[(start + k - 1) % pts.size()].y);
    }
    return length;
}

double lineFitRmse(const std::vector<pcl::PointXYZ>& pts, size_t start, size_t count)
{
    if (pts.empty() || count < 3) return 0.0;
    const auto& a = pts[start % pts.size()];
    const auto& b = pts[(start + count - 1) % pts.size()];
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double len = std::hypot(dx, dy);
    if (len < 1e-6) return 0.0;
    double err2 = 0.0;
    for (size_t k = 1; k + 1 < count; ++k) {
        const auto& p = pts[(start + k) % pts.size()];
        const double distance = std::abs((p.x - a.x) * dy - (p.y - a.y) * dx) / len;
        err2 += distance * distance;
    }
    return std::sqrt(err2 / static_cast<double>(std::max<size_t>(1, count - 2)));
}

double accumulatedTurnAngle(const std::vector<pcl::PointXYZ>& pts, size_t start, size_t count)
{
    if (pts.empty() || count < 4) return 0.0;
    double total = 0.0;
    for (size_t k = 1; k + 1 < count; ++k) {
        const auto& a = pts[(start + k - 1) % pts.size()];
        const auto& b = pts[(start + k) % pts.size()];
        const auto& c = pts[(start + k + 1) % pts.size()];
        const double a1 = std::atan2(b.y - a.y, b.x - a.x);
        const double a2 = std::atan2(c.y - b.y, c.x - b.x);
        total += std::abs(normalizeAnglePi(a2 - a1));
    }
    return total;
}

bool fitCubicBezier2D(const std::vector<pcl::PointXYZ>& pts,
    size_t start,
    size_t count,
    Eigen::Vector2d& p0,
    Eigen::Vector2d& p1,
    Eigen::Vector2d& p2,
    Eigen::Vector2d& p3,
    double& rmse)
{
    if (pts.empty() || count < 6) return false;
    p0 = Eigen::Vector2d(pts[start % pts.size()].x, pts[start % pts.size()].y);
    p3 = Eigen::Vector2d(pts[(start + count - 1) % pts.size()].x,
                         pts[(start + count - 1) % pts.size()].y);

    std::vector<double> t_values(count, 0.0);
    double total_length = 0.0;
    for (size_t k = 1; k < count; ++k) {
        total_length += std::hypot(
            pts[(start + k) % pts.size()].x - pts[(start + k - 1) % pts.size()].x,
            pts[(start + k) % pts.size()].y - pts[(start + k - 1) % pts.size()].y);
        t_values[k] = total_length;
    }
    if (total_length < 1e-6) return false;
    for (double& t : t_values) t /= total_length;

    Eigen::Matrix2d normal = Eigen::Matrix2d::Zero();
    Eigen::Vector2d rhs_x = Eigen::Vector2d::Zero();
    Eigen::Vector2d rhs_y = Eigen::Vector2d::Zero();
    for (size_t k = 1; k + 1 < count; ++k) {
        const double t = t_values[k];
        const double u = 1.0 - t;
        const double b0 = u * u * u;
        const double b1 = 3.0 * u * u * t;
        const double b2 = 3.0 * u * t * t;
        const double b3 = t * t * t;
        const auto& p = pts[(start + k) % pts.size()];
        const double qx = p.x - b0 * p0.x() - b3 * p3.x();
        const double qy = p.y - b0 * p0.y() - b3 * p3.y();
        Eigen::Vector2d row(b1, b2);
        normal += row * row.transpose();
        rhs_x += row * qx;
        rhs_y += row * qy;
    }
    if (std::abs(normal.determinant()) < 1e-10) return false;
    const Eigen::Vector2d control_x = normal.ldlt().solve(rhs_x);
    const Eigen::Vector2d control_y = normal.ldlt().solve(rhs_y);
    if (!control_x.allFinite() || !control_y.allFinite()) return false;
    p1 = Eigen::Vector2d(control_x.x(), control_y.x());
    p2 = Eigen::Vector2d(control_x.y(), control_y.y());

    double err2 = 0.0;
    for (size_t k = 0; k < count; ++k) {
        const double t = t_values[k];
        const double u = 1.0 - t;
        const Eigen::Vector2d q =
            (u * u * u) * p0 +
            (3.0 * u * u * t) * p1 +
            (3.0 * u * t * t) * p2 +
            (t * t * t) * p3;
        const auto& p = pts[(start + k) % pts.size()];
        const double dx = q.x() - p.x;
        const double dy = q.y() - p.y;
        err2 += dx * dx + dy * dy;
    }
    rmse = std::sqrt(err2 / static_cast<double>(count));
    return std::isfinite(rmse);
}

std::vector<pcl::PointXYZ> sampleArc(const PreservedArcSegment& arc, double step)
{
    std::vector<pcl::PointXYZ> result;
    if (arc.use_bezier) {
        const double chord_length = std::hypot(
            arc.end_point.x - arc.start_point.x,
            arc.end_point.y - arc.start_point.y);
        const int segments = std::max(4, static_cast<int>(
            std::ceil(std::max(arc.curve_length, chord_length) /
                      std::max(step, 0.15))));
        result.reserve(static_cast<size_t>(segments + 1));
        for (int i = 0; i <= segments; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(segments);
            const double u = 1.0 - t;
            const Eigen::Vector2d q =
                (u * u * u) * arc.bezier_p0 +
                (3.0 * u * u * t) * arc.bezier_p1 +
                (3.0 * u * t * t) * arc.bezier_p2 +
                (t * t * t) * arc.bezier_p3;
            pcl::PointXYZ p;
            p.x = static_cast<float>(q.x());
            p.y = static_cast<float>(q.y());
            p.z = static_cast<float>(arc.start_point.z * (1.0 - t) + arc.end_point.z * t);
            result.push_back(p);
        }
        result.front() = arc.start_point;
        result.back() = arc.end_point;
        return result;
    }

    const double length = std::abs(arc.sweep_angle) * arc.radius;
    const int segments = std::max(4, static_cast<int>(std::ceil(length / std::max(step, 0.15))));
    result.reserve(static_cast<size_t>(segments + 1));
    for (int i = 0; i <= segments; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(segments);
        const double a = arc.start_angle + arc.sweep_angle * t;
        pcl::PointXYZ p;
        p.x = static_cast<float>(arc.cx + arc.radius * std::cos(a));
        p.y = static_cast<float>(arc.cy + arc.radius * std::sin(a));
        p.z = static_cast<float>(arc.start_point.z * (1.0 - t) + arc.end_point.z * t);
        result.push_back(p);
    }
    result.front() = arc.start_point;
    result.back() = arc.end_point;
    return result;
}

std::vector<PreservedArcSegment> detectPreservedArcs(
    const std::vector<pcl::PointXYZ>& ring,
    double resolution)
{
    std::vector<PreservedArcSegment> arcs;
    if (ring.size() < 10) return arcs;

    const double area_scale = std::sqrt(std::max(polygonArea2D(ring), 1.0));
    const size_t min_count = 7;
    const size_t max_count = std::min<size_t>(ring.size() / 2, 48);
    if (max_count < min_count) return arcs;
    const double max_rmse = clampDouble(
        std::max(2.5 * resolution, 0.006 * area_scale), 0.12, 0.45);
    const double min_radius = clampDouble(0.02 * area_scale, 1.0, 8.0);
    const double max_radius = std::max(6.0, 3.0 * area_scale);
    const double min_sweep = 20.0 * M_PI / 180.0;
    const double max_sweep = 170.0 * M_PI / 180.0;
    const double min_arc_length = clampDouble(0.04 * area_scale, 2.0, 8.0);
    const double max_curve_rmse = clampDouble(
        std::max(3.0 * resolution, 0.008 * area_scale), 0.15, 0.6);
    const double min_line_curve_ratio = 1.8;
    const double min_turn = 12.0 * M_PI / 180.0;

    struct Candidate {
        PreservedArcSegment arc;
        double score = 0.0;
        size_t count = 0;
    };
    std::vector<Candidate> candidates;
    for (size_t start = 0; start < ring.size(); ++start) {
        for (size_t count = min_count; count <= max_count; ++count) {
            double cx = 0.0, cy = 0.0, radius = 0.0, rmse = 0.0;
            if (!fitCircle2D(ring, start, count, cx, cy, radius, rmse)) continue;
            if (radius < min_radius || radius > max_radius || rmse > max_rmse) continue;
            double start_angle = 0.0;
            const double sweep = unwrapArcSweep(ring, start, count, cx, cy, start_angle);
            const double abs_sweep = std::abs(sweep);
            const double arc_length = abs_sweep * radius;
            if (abs_sweep < min_sweep || abs_sweep > max_sweep || arc_length < min_arc_length) continue;

            PreservedArcSegment arc;
            arc.start = start;
            arc.end = (start + count - 1) % ring.size();
            arc.start_point = ring[arc.start];
            arc.end_point = ring[arc.end];
            arc.cx = cx;
            arc.cy = cy;
            arc.radius = radius;
            arc.start_angle = start_angle;
            arc.sweep_angle = sweep;
            arc.rmse = rmse;
            arc.curve_length = arc_length;
            const double score = arc_length / std::max(rmse, 0.03);
            candidates.push_back({ arc, score, count });
        }
    }

    for (size_t start = 0; start < ring.size(); ++start) {
        for (size_t count = min_count; count <= max_count; ++count) {
            const double curve_length = polylineLength(ring, start, count);
            if (curve_length < min_arc_length) continue;
            const double straight_rmse = lineFitRmse(ring, start, count);
            if (straight_rmse < max_curve_rmse * min_line_curve_ratio) continue;
            if (accumulatedTurnAngle(ring, start, count) < min_turn) continue;

            Eigen::Vector2d p0, p1, p2, p3;
            double bezier_rmse = 0.0;
            if (!fitCubicBezier2D(ring, start, count, p0, p1, p2, p3, bezier_rmse)) continue;
            if (bezier_rmse > max_curve_rmse) continue;
            if (straight_rmse / std::max(bezier_rmse, 0.03) < min_line_curve_ratio) continue;

            PreservedArcSegment curve;
            curve.start = start;
            curve.end = (start + count - 1) % ring.size();
            curve.start_point = ring[curve.start];
            curve.end_point = ring[curve.end];
            curve.use_bezier = true;
            curve.rmse = bezier_rmse;
            curve.curve_length = curve_length;
            curve.bezier_p0 = p0;
            curve.bezier_p1 = p1;
            curve.bezier_p2 = p2;
            curve.bezier_p3 = p3;
            const double score = curve_length * straight_rmse / std::max(bezier_rmse, 0.03);
            candidates.push_back({ curve, score, count });
        }
    }

    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b) { return a.score > b.score; });

    std::vector<bool> used(ring.size(), false);
    for (const auto& candidate : candidates) {
        size_t occupied = 0;
        for (size_t k = 0; k < candidate.count; ++k) {
            if (used[(candidate.arc.start + k) % ring.size()]) ++occupied;
        }
        if (occupied > candidate.count / 5) continue;
        for (size_t k = 0; k < candidate.count; ++k) {
            used[(candidate.arc.start + k) % ring.size()] = true;
        }
        arcs.push_back(candidate.arc);
        if (arcs.size() >= 4) break;
    }

    std::sort(arcs.begin(), arcs.end(),
        [](const PreservedArcSegment& a, const PreservedArcSegment& b) { return a.start < b.start; });
    if (!arcs.empty()) {
        const size_t bezier_count = static_cast<size_t>(std::count_if(
            arcs.begin(), arcs.end(), [](const auto& arc) { return arc.use_bezier; }));
        std::cerr << "[Curve] detected preserved curves=" << arcs.size()
            << " bezier=" << bezier_count << std::endl;
    }
    return arcs;
}

bool restorePreservedArcs(
    std::vector<pcl::PointXYZ>& polygon,
    const std::vector<PreservedArcSegment>& arcs,
    double resolution)
{
    if (polygon.size() < 3 || arcs.empty()) return false;
    bool changed = false;
    const double step = clampDouble(std::max(2.0 * resolution, 0.25), 0.25, 0.8);
    const double area_scale = std::sqrt(std::max(polygonArea2D(polygon), 1.0));
    const double endpoint_tolerance = clampDouble(0.03 * area_scale, 1.0, 4.0);

    for (const auto& arc : arcs) {
        if (polygon.size() < 3) break;
        int best_edge = -1;
        double best_score = std::numeric_limits<double>::max();
        for (size_t i = 0; i < polygon.size(); ++i) {
            const auto& a = polygon[i];
            const auto& b = polygon[(i + 1) % polygon.size()];
            const double direct = std::hypot(a.x - arc.start_point.x, a.y - arc.start_point.y) +
                std::hypot(b.x - arc.end_point.x, b.y - arc.end_point.y);
            const double reverse = std::hypot(a.x - arc.end_point.x, a.y - arc.end_point.y) +
                std::hypot(b.x - arc.start_point.x, b.y - arc.start_point.y);
            const double score = std::min(direct, reverse);
            if (score < best_score) {
                best_score = score;
                best_edge = static_cast<int>(i);
            }
        }
        if (best_edge < 0 || best_score > 2.0 * endpoint_tolerance) continue;

        const size_t edge = static_cast<size_t>(best_edge);
        const auto& a = polygon[edge];
        const auto& b = polygon[(edge + 1) % polygon.size()];
        const double direct = std::hypot(a.x - arc.start_point.x, a.y - arc.start_point.y) +
            std::hypot(b.x - arc.end_point.x, b.y - arc.end_point.y);
        const double reverse = std::hypot(a.x - arc.end_point.x, a.y - arc.end_point.y) +
            std::hypot(b.x - arc.start_point.x, b.y - arc.start_point.y);
        PreservedArcSegment oriented = arc;
        if (reverse < direct) {
            oriented.start_point = arc.end_point;
            oriented.end_point = arc.start_point;
            oriented.start_angle = arc.start_angle + arc.sweep_angle;
            oriented.sweep_angle = -arc.sweep_angle;
        }

        std::vector<pcl::PointXYZ> samples = sampleArc(oriented, step);
        if (samples.size() < 3) continue;
        samples.front() = polygon[edge];
        samples.back() = polygon[(edge + 1) % polygon.size()];

        std::vector<pcl::PointXYZ> next;
        next.reserve(polygon.size() + samples.size());
        for (size_t i = 0; i < polygon.size(); ++i) {
            if (i == edge) {
                next.insert(next.end(), samples.begin(), samples.end() - 1);
            } else {
                next.push_back(polygon[i]);
            }
        }
        removeDuplicatePoints2D(next, 0.05f);
        if (next.size() >= 3 && isSimplePolygon2D(next)) {
            polygon.swap(next);
            changed = true;
            if (arc.use_bezier) {
                std::cerr << "[Curve] restored bezier length=" << arc.curve_length
                    << " rmse=" << arc.rmse << std::endl;
            } else {
                std::cerr << "[Curve] restored arc radius=" << arc.radius
                    << " sweep_deg=" << std::abs(arc.sweep_angle) * 180.0 / M_PI
                    << " rmse=" << arc.rmse << std::endl;
            }
        }
    }
    return changed;
}

void simplifySameOrthogonalAxisVertices(std::vector<pcl::PointXYZ>& polygon)
{
    if (polygon.size() < 4) return;
    bool changed = true;
    while (changed && polygon.size() >= 4) {
        changed = false;
        const auto dominant = dominantLineAngles2D(polygon, 1);
        if (dominant.empty()) return;
        const double main_angle = dominant.front();
        auto axis = [main_angle](const pcl::PointXYZ& a, const pcl::PointXYZ& b) {
            const double angle = std::atan2(b.y - a.y, b.x - a.x);
            return undirectedAngleDifference(angle, main_angle) <=
                undirectedAngleDifference(angle, main_angle + M_PI / 2.0) ? 0 : 1;
        };
        for (size_t i = 0; i < polygon.size(); ++i) {
            const size_t previous = (i + polygon.size() - 1) % polygon.size();
            const size_t next = (i + 1) % polygon.size();
            if (axis(polygon[previous], polygon[i]) == axis(polygon[i], polygon[next])) {
                polygon.erase(polygon.begin() + i);
                changed = true;
                break;
            }
        }
    }
}

void removeNearlyCollinearPoints2D(std::vector<pcl::PointXYZ>& polygon, double threshold_deg)
{
    if (polygon.size() < 4) return;
    bool changed = true;
    const double threshold_rad = threshold_deg * M_PI / 180.0;
    while (changed && polygon.size() >= 4) {
        changed = false;
        for (size_t i = 0; i < polygon.size(); ++i) {
            const auto& prev = polygon[(i + polygon.size() - 1) % polygon.size()];
            const auto& cur = polygon[i];
            const auto& next = polygon[(i + 1) % polygon.size()];
            Eigen::Vector2d a(prev.x - cur.x, prev.y - cur.y);
            Eigen::Vector2d b(next.x - cur.x, next.y - cur.y);
            if (a.norm() < 1e-8 || b.norm() < 1e-8) continue;
            const double cos_angle = clampDouble(a.dot(b) / (a.norm() * b.norm()), -1.0, 1.0);
            const double angle = std::acos(cos_angle);
            if (std::abs(M_PI - angle) <= threshold_rad) {
                std::vector<pcl::PointXYZ> candidate = polygon;
                candidate.erase(candidate.begin() + i);
                if (candidate.size() >= 3 && isSimplePolygon2D(candidate)) {
                    polygon.swap(candidate);
                    changed = true;
                    break;
                }
            }
        }
    }
}

// ===== isRegularizedPolygonAcceptable =====
// 作用：检查规则化后的候选多边形是否可接受（简单、非退化、面积比/IoU合理）
// 参数：candidate - 候选多边形; reference - 参考多边形（用于面积比校验）
// 返回：true 表示结果可接受
bool isRegularizedPolygonAcceptable(const std::vector<pcl::PointXYZ>& candidate,
    const std::vector<pcl::PointXYZ>& reference)
{
    if (candidate.size() < 3) {
        std::cerr << "[AcceptCheck] reject final polygon: vertex_count="
            << candidate.size() << std::endl;
        return false;
    }
    if (!isSimplePolygon2D(candidate)) {
        std::cerr << "[AcceptCheck] reject final polygon: self-intersection or invalid ring"
            << std::endl;
        return false;
    }

    double candidate_area = polygonArea2D(candidate);
    if (candidate_area < 1e-6) {
        std::cerr << "[AcceptCheck] reject final polygon: near-zero area"
            << std::endl;
        return false;
    }

    double reference_area = polygonArea2D(reference);
    if (reference_area > 1e-6) {
        double ratio = candidate_area / reference_area;
        if (ratio < 0.75 || ratio > 2.8) {
            std::cerr << "[AcceptCheck] reject final polygon: area_ratio="
                << ratio << std::endl;
            return false;
        }

        const double iou = polygonIoU2D(candidate, reference);
        if (iou < 0.75) {
            std::cerr << "[AcceptCheck] reject final polygon: IoU="
                << iou << ", area_ratio=" << ratio << std::endl;
            return false;
        }
    }

    return true;
}

struct BoundaryDistanceStats {
    double mean = std::numeric_limits<double>::max();
    double q90 = std::numeric_limits<double>::max();
};

void logPolygonBBox2D(
    const char* label,
    const std::vector<pcl::PointXYZ>& polygon)
{
    if (!label || polygon.empty()) return;
    double min_x = polygon.front().x;
    double max_x = polygon.front().x;
    double min_y = polygon.front().y;
    double max_y = polygon.front().y;
    for (const auto& p : polygon) {
        min_x = std::min(min_x, static_cast<double>(p.x));
        max_x = std::max(max_x, static_cast<double>(p.x));
        min_y = std::min(min_y, static_cast<double>(p.y));
        max_y = std::max(max_y, static_cast<double>(p.y));
    }
    std::cerr << label << " bbox X=" << min_x << ".." << max_x
        << " Y=" << min_y << ".." << max_y
        << " n=" << polygon.size() << std::endl;
}

BoundaryDistanceStats boundaryDistanceStats2D(
    const std::vector<pcl::PointXYZ>& source,
    const std::vector<pcl::PointXYZ>& target)
{
    BoundaryDistanceStats stats;
    if (source.size() < 3 || target.size() < 3) return stats;

    std::vector<double> distances;
    distances.reserve(source.size() * 2);
    for (size_t i = 0; i < source.size(); ++i) {
        const auto& a = source[i];
        const auto& b = source[(i + 1) % source.size()];
        const double samples[2][2] = {
            { static_cast<double>(a.x), static_cast<double>(a.y) },
            {
                static_cast<double>(a.x) +
                    0.5 * (static_cast<double>(b.x) - static_cast<double>(a.x)),
                static_cast<double>(a.y) +
                    0.5 * (static_cast<double>(b.y) - static_cast<double>(a.y))
            }
        };
        for (const auto& sample : samples) {
            double best = std::numeric_limits<double>::max();
            for (size_t j = 0; j < target.size(); ++j) {
                best = std::min(best, distanceXYToSegment2D(
                    sample[0], sample[1],
                    target[j], target[(j + 1) % target.size()]));
            }
            if (std::isfinite(best)) distances.push_back(best);
        }
    }
    if (distances.empty()) return stats;

    const double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
    std::sort(distances.begin(), distances.end());
    const size_t q90_index = std::min(distances.size() - 1,
        static_cast<size_t>(std::floor(0.90 * static_cast<double>(distances.size() - 1))));
    stats.mean = sum / static_cast<double>(distances.size());
    stats.q90 = distances[q90_index];
    return stats;
}

bool isSingleDirectionCandidateAcceptable(
    const std::vector<pcl::PointXYZ>& candidate,
    const std::vector<pcl::PointXYZ>& reference,
    const OutlineTuning& tuning)
{
    if (isRegularizedPolygonAcceptable(candidate, reference)) return true;
    if (candidate.size() < 3 || reference.size() < 3 || !isSimplePolygon2D(candidate)) {
        return false;
    }

    const double reference_area = std::max(polygonArea2D(reference), 1e-6);
    const double candidate_area = polygonArea2D(candidate);
    const double area_ratio = candidate_area / reference_area;
    if (area_ratio < 0.65 || area_ratio > 1.45) {
        std::cerr << "[SingleQuality] reject area_ratio=" << area_ratio << std::endl;
        return false;
    }

    const double scale = std::sqrt(reference_area);
    const double mean_limit = clampDouble(
        std::max(2.5 * tuning.resolution, 0.035 * scale), 0.8, 3.0);
    const double q90_limit = clampDouble(
        std::max(4.0 * tuning.resolution, 0.075 * scale), 1.5, 5.0);

    const BoundaryDistanceStats cand_to_ref =
        boundaryDistanceStats2D(candidate, reference);
    const BoundaryDistanceStats ref_to_cand =
        boundaryDistanceStats2D(reference, candidate);
    const double mean = std::max(cand_to_ref.mean, ref_to_cand.mean);
    const double q90 = std::max(cand_to_ref.q90, ref_to_cand.q90);
    bool ok = mean <= mean_limit && q90 <= q90_limit;
    if (!ok) {
        const double iou = polygonIoU2D(candidate, reference);
        // 宽松通道的面积放大必须设上限：大建筑(scale 可达 100m+)原本会
        // 得到 8m/20m 的容忍度，足以放行"整个斜翼被单方向化"的系统性
        // 形变(实测 IoU=0.69 / mean=4.08 / q90=11.7 被放行)。
        const double relaxed_mean_limit =
            std::min(std::max(2.5, 0.08 * scale), 4.5);
        const double relaxed_q90_limit =
            std::min(std::max(4.5, 0.20 * scale), 9.0);
        const bool relaxed_ok =
            iou >= 0.72 &&
            area_ratio >= 0.74 && area_ratio <= 1.25 &&
            mean <= relaxed_mean_limit &&
            q90 <= relaxed_q90_limit;
        if (relaxed_ok) {
            std::cerr << "[SingleQuality] relaxed accept"
                << " IoU=" << iou
                << " area_ratio=" << area_ratio
                << " mean=" << mean
                << " q90=" << q90 << std::endl;
            ok = true;
        }
    }
    std::cerr << "[SingleQuality] fallback_check area_ratio=" << area_ratio
        << " mean=" << mean << "/" << mean_limit
        << " q90=" << q90 << "/" << q90_limit
        << " accepted=" << (ok ? 1 : 0) << std::endl;
    if (!ok && std::isfinite(mean) && mean > 50.0) {
        logPolygonBBox2D("[SingleQuality] candidate", candidate);
        logPolygonBBox2D("[SingleQuality] reference", reference);
    }
    return ok;
}
}

bool outlineRegular::estimateSupportDirection2D(
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support,
    double& angle,
    double& peakRatio,
    std::size_t& pairCount)
{
    angle = 0.0;
    peakRatio = 0.0;
    pairCount = 0;
    constexpr std::size_t kMaxSamples = 5000;
    if (!support || support->size() < kSupportDirectionMinPoints) return false;

    pcl::PointCloud<pcl::PointXYZ>::Ptr support2D(new pcl::PointCloud<pcl::PointXYZ>);
    support2D->resize(support->size());
    for (std::size_t i = 0; i < support->size(); ++i) {
        support2D->points[i].x = support->points[i].x;
        support2D->points[i].y = support->points[i].y;
        support2D->points[i].z = 0.0f;
    }
    pcl::KdTreeFLANN<pcl::PointXYZ> tree;
    tree.setInputCloud(support2D);
    std::array<double, 90> bins{};
    const std::size_t stride = std::max<std::size_t>(1, support->size() / kMaxSamples);
    std::vector<int> indices;
    std::vector<float> distances;
    double totalWeight = 0.0;
    for (std::size_t i = 0; i < support->size(); i += stride) {
        indices.clear();
        distances.clear();
        tree.radiusSearch(support2D->points[i], static_cast<float>(kSupportDirectionPairRadius), indices, distances, 96);
        for (std::size_t k = 0; k < indices.size(); ++k) {
            const int j = indices[k];
            if (j <= static_cast<int>(i) ||
                distances[k] < kSupportDirectionMinPairDistance * kSupportDirectionMinPairDistance) continue;
            const auto& a = support->points[i];
            const auto& b = support->points[static_cast<std::size_t>(j)];
            double degrees = foldedLineAngle90(std::atan2(b.y - a.y, b.x - a.x)) * 180.0 / M_PI;
            const int bin = std::min(89, std::max(0, static_cast<int>(std::floor(degrees))));
            const double weight = std::sqrt(std::max(0.0f, distances[k]));
            bins[static_cast<std::size_t>(bin)] += weight;
            totalWeight += weight;
            ++pairCount;
        }
    }
    if (pairCount < kSupportDirectionMinPairs || totalWeight <= 1e-9) return false;

    int bestBin = 0;
    double bestWindow = -1.0;
    for (int center = 0; center < 90; ++center) {
        double window = 0.0;
        for (int d = -5; d <= 5; ++d) {
            window += bins[static_cast<std::size_t>((center + d + 90) % 90)];
        }
        if (window > bestWindow) {
            bestWindow = window;
            bestBin = center;
        }
    }
    peakRatio = bestWindow / totalWeight;
    angle = (static_cast<double>(bestBin) + 0.5) * M_PI / 180.0;
    return peakRatio >= kSupportDirectionStrongPeakRatio;
}

bool outlineRegular::estimatePcaDirection2D(
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& points,
    double& angle,
    double& axisRatio)
{
    angle = 0.0;
    axisRatio = 1.0;
    if (!points || points->size() < 3) return false;
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    for (const auto& p : points->points) center += Eigen::Vector2d(p.x, p.y);
    center /= static_cast<double>(points->size());
    Eigen::Matrix2d covariance = Eigen::Matrix2d::Zero();
    for (const auto& p : points->points) {
        const Eigen::Vector2d d = Eigen::Vector2d(p.x, p.y) - center;
        covariance += d * d.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(covariance);
    if (solver.info() != Eigen::Success) return false;
    const auto values = solver.eigenvalues();
    const auto axis = solver.eigenvectors().col(1);
    angle = foldedLineAngle90(std::atan2(axis.y(), axis.x()));
    axisRatio = std::sqrt(std::max(values[1], 1e-9) / std::max(values[0], 1e-9));
    return std::isfinite(angle) && std::isfinite(axisRatio);
}

void outlineRegular::setSupportDirectionHint(
    double angle, double peakRatio, std::size_t pairCount)
{
    has_support_direction_hint_ = true;
    support_direction_hint_ = foldedLineAngle90(angle);
    support_direction_peak_ratio_ = peakRatio;
    support_direction_pair_count_ = pairCount;
}

// ===== PointToLineDistanceCost =====
// 作用：Ceres 残差结构体——约束受约束边（平行/垂直）的点到直线距离
// 公式：x * cos(theta) + y * sin(theta) - d = 0
struct PointToLineDistanceCost {
    // 构造函数：初始化时传入点的坐标 (x, y) 和 角度偏移量 angle_offset
    PointToLineDistanceCost(double x, double y, double angle_offset, double sqrt_weight = 1.0)
        : x_(x), y_(y), angle_offset_(angle_offset), sqrt_weight_(sqrt_weight) {}

    // Ceres 要求的重载运算符 ()，用于计算残差
    // T 是模板类型，Ceres 会传入自动求导的特殊类型
    template <typename T>
    bool operator()(const T* const base_theta, const T* const d, T* residual) const {

        // 【核心逻辑】：当前这条边的实际角度 = 全局主朝向 base_theta + 固定的偏移量
        // 这里的 angle_offset_ 通常是 0, 90, 180 或 270 度
        T theta = base_theta[0] + T(angle_offset_);

        // 计算点到直线的距离（使用 Hessian Normal Form：海塞正规式）
        // 公式：r = x * cos(theta) + y * sin(theta) - d
        // 优化器的目标就是让这个 residual[0] 趋近于 0
        residual[0] = T(sqrt_weight_) *
            (T(x_) * ceres::cos(theta) + T(y_) * ceres::sin(theta) - d[0]);
        return true;
    }

    // 静态工厂函数：用于创建 Ceres 能够识别的自动求导代价函数
    static ceres::CostFunction* Create(double x, double y, double angle_offset, double sqrt_weight = 1.0) {
        // <结构体, 残差维度(1), 参数块1维度(base_theta是1维), 参数块2维度(d是1维)>
        return (new ceres::AutoDiffCostFunction<PointToLineDistanceCost, 1, 1, 1>(
            new PointToLineDistanceCost(x, y, angle_offset, sqrt_weight)));
    }

    double x_, y_, angle_offset_, sqrt_weight_; // 存储点的坐标、偏移量和残差权重
};

// ===== FreeLineCost =====
// 作用：Ceres 残差结构体——对于不属于平行/垂直系统的自由边，使用独立的 theta
struct FreeLineCost {
    // 构造函数：只传入点的坐标
    FreeLineCost(double x, double y, double sqrt_weight = 1.0)
        : x_(x), y_(y), sqrt_weight_(sqrt_weight) {}

    template <typename T>
    bool operator()(const T* const theta, const T* const d, T* residual) const {
        // 【不同点】：这里的 theta 是独立的，不依赖于任何全局变量
        // 每条自由边都有自己独立的旋转角度参数
        residual[0] = T(sqrt_weight_) *
            (T(x_) * ceres::cos(theta[0]) + T(y_) * ceres::sin(theta[0]) - d[0]);
        return true;
    }

    static ceres::CostFunction* Create(double x, double y, double sqrt_weight = 1.0) {
        // <结构体, 残差维度(1), 参数块1维度(独立的theta), 参数块2维度(d)>
        return (new ceres::AutoDiffCostFunction<FreeLineCost, 1, 1, 1>(
            new FreeLineCost(x, y, sqrt_weight)));
    }
    double x_, y_, sqrt_weight_;
};

// ===== VertexPreservationResidual =====
// 作用：Ceres 残差结构体——顶点保形约束（两条相邻直线的交点应接近原始顶点）
struct VertexPreservationResidual {
    // 增加 offset1, offset2 参数
    VertexPreservationResidual(double x0, double y0, double offset1, double offset2)
        : x0_(x0), y0_(y0), offset1_(offset1), offset2_(offset2) {}

    template <typename T>
    bool operator()(const T* const theta1_param, const T* const d1,
        const T* const theta2_param, const T* const d2,
        T* residual) const {
        // 1. 恢复真实角度：参数值 + 偏移量
        // 对于自由边，offset 为 0，theta1_param 就是真实角度
        // 对于约束边，offset 为 0/90/180...，theta1_param 是 base_theta
        T theta1 = theta1_param[0] + T(offset1_);
        T theta2 = theta2_param[0] + T(offset2_);

        // 2. 直线1: x*cos(t1) + y*sin(t1) = d1
        //    直线2: x*cos(t2) + y*sin(t2) = d2
        T cos1 = ceres::cos(theta1);
        T sin1 = ceres::sin(theta1);
        T cos2 = ceres::cos(theta2);
        T sin2 = ceres::sin(theta2);

        // 3. 计算行列式 det = sin(t2 - t1)
        T det = cos1 * sin2 - cos2 * sin1;

        // 4. 平行保护
        if (ceres::abs(det) < T(1e-5)) {
            // 平行无交点，给予较大惩罚让其分开或忽略
            residual[0] = T(0.0);
            residual[1] = T(0.0);
            return true;
        }

        // 5. Cramer法则求交点
        T x = (d1[0] * sin2 - d2[0] * sin1) / det;
        T y = (d2[0] * cos1 - d1[0] * cos2) / det;

        // 6. 计算与原始顶点的偏差
        residual[0] = x - T(x0_);
        residual[1] = y - T(y0_);
        return true;
    }

    static ceres::CostFunction* Create(double x0, double y0, double offset1, double offset2) {
        return new ceres::AutoDiffCostFunction<VertexPreservationResidual, 2, 1, 1, 1, 1>(
            new VertexPreservationResidual(x0, y0, offset1, offset2));
    }

private:
    double x0_, y0_;
    double offset1_, offset2_;
};

// ===== ParameterRegularizer =====
// 作用：Ceres 残差结构体——正则化约束，防止参数偏离初始值太远
struct ParameterRegularizer {
    ParameterRegularizer(double initial_value, double weight)
        : initial_value_(initial_value), weight_(weight) {}

    template <typename T>
    bool operator()(const T* const current_value, T* residual) const {
        // 残差 = (当前值 - 初始值) * 权重
        residual[0] = (current_value[0] - T(initial_value_)) * T(weight_);
        return true;
    }

    static ceres::CostFunction* Create(double initial_value, double weight) {
        return new ceres::AutoDiffCostFunction<ParameterRegularizer, 1, 1>(
            new ParameterRegularizer(initial_value, weight));
    }

    double initial_value_, weight_;
};

// ===== LineAnchorCost =====
// 作用：Ceres 残差结构体——锚点约束，防止直线沿轴向滑动
// 原理：约束直线的方程必须经过原始线段的中点
struct LineAnchorCost {
    LineAnchorCost(double mid_x, double mid_y, double offset, double weight)
        : mid_x_(mid_x), mid_y_(mid_y), offset_(offset), weight_(weight) {}

    template <typename T>
    bool operator()(const T* const base_theta, const T* const d, T* residual) const {
        T theta = base_theta[0] + T(offset_);
        // 直线方程：x*cos(theta) + y*sin(theta) - d = 0
        residual[0] = (mid_x_ * ceres::cos(theta) + mid_y_ * ceres::sin(theta) - d[0]) * T(weight_);
        return true;
    }

    static ceres::CostFunction* Create(double mid_x, double mid_y, double offset, double weight) {
        return new ceres::AutoDiffCostFunction<LineAnchorCost, 1, 1, 1>(
            new LineAnchorCost(mid_x, mid_y, offset, weight));
    }

    double mid_x_, mid_y_, offset_, weight_;
};

// ===== FreeLineAnchorCost =====
// 作用：Ceres 残差结构体——自由边的锚点约束（针对没有 base_theta 的边）
struct FreeLineAnchorCost {
    FreeLineAnchorCost(double mid_x, double mid_y, double weight)
        : mid_x_(mid_x), mid_y_(mid_y), weight_(weight) {}

    template <typename T>
    bool operator()(const T* const theta, const T* const d, T* residual) const {
        residual[0] = (mid_x_ * ceres::cos(theta[0]) + mid_y_ * ceres::sin(theta[0]) - d[0]) * T(weight_);
        return true;
    }

    static ceres::CostFunction* Create(double mid_x, double mid_y, double weight) {
        return new ceres::AutoDiffCostFunction<FreeLineAnchorCost, 1, 1, 1>(
            new FreeLineAnchorCost(mid_x, mid_y, weight));
    }

    double mid_x_, mid_y_, weight_;
};


// ===== outlineRegular (默认构造) =====
// 作用：默认构造函数（空对象，成员需后续设置）
outlineRegular::outlineRegular() {};

// ===== outlineRegular (两参数构造) =====
// 作用：用原始轮廓点和 RANSAC 内点点云初始化规则化对象
// 参数：p - 原始轮廓顶点; ransac_inner_cloud - RANSAC 拟合内点点云
outlineRegular::outlineRegular(vector<pcl::PointXYZ> p, pcl::PointCloud<pcl::PointXYZ>::Ptr ransac_inner_cloud)
{
    final_points.reset(new pcl::PointCloud<pcl::PointXYZ>);
    this->original_points = p;
    this->ransac_inner_cloud = ransac_inner_cloud;
    if (ransac_inner_cloud) {
        support_weights_.assign(ransac_inner_cloud->size(), 1.0);
    }
    if(!original_points.empty())GraphZ=original_points[0].z;

    // 计算ransac_inner_cloud的二维包围盒并扩展20.0，用于最终结果越界检测
    ransac_min_x = 1e10; ransac_max_x = -1e10;
    ransac_min_y = 1e10; ransac_max_y = -1e10;
    if (ransac_inner_cloud && !ransac_inner_cloud->empty()) {
        for (const auto& pt : ransac_inner_cloud->points) {
            if (pt.x < ransac_min_x) ransac_min_x = pt.x;
            if (pt.x > ransac_max_x) ransac_max_x = pt.x;
            if (pt.y < ransac_min_y) ransac_min_y = pt.y;
            if (pt.y > ransac_max_y) ransac_max_y = pt.y;
        }
    }
    const double bbox_margin = 20.0;
    ransac_min_x -= bbox_margin; ransac_max_x += bbox_margin;
    ransac_min_y -= bbox_margin; ransac_max_y += bbox_margin;

    pcl::PointCloud<pcl::PointXYZ> p_cloud;
    for (auto& op : original_points) {p_cloud.points.push_back(op);}
    //pcl::io::savePCDFileBinaryCompressed("original_points.pcd", p_cloud);
    //pcl::io::savePCDFileBinaryCompressed("ransac_inner_cloud.pcd", *ransac_inner_cloud);

}

outlineRegular::outlineRegular(vector<pcl::PointXYZ> p,
    pcl::PointCloud<pcl::PointXYZ>::Ptr ransac_inner_cloud,
    const std::vector<double>& support_weights)
    : outlineRegular(std::move(p), ransac_inner_cloud)
{
    if (ransac_inner_cloud && support_weights.size() == ransac_inner_cloud->size()) {
        support_weights_ = support_weights;
    }
    else if (ransac_inner_cloud) {
        support_weights_.assign(ransac_inner_cloud->size(), 1.0);
    }
}

// ===== outlineRegular (三参数构造) =====
// 作用：在两参数构造基础上，额外传入 DLG（底图）先验多边形用于跨层引导
// 参数：p - 原始轮廓顶点; ransac_inner_cloud - RANSAC 内点点云; dlg_polygon - DLG 底图多边形
outlineRegular::outlineRegular(vector<pcl::PointXYZ> p,
    pcl::PointCloud<pcl::PointXYZ>::Ptr ransac_inner_cloud,
    const std::vector<pcl::PointXYZ>& dlg_polygon)
    : outlineRegular(std::move(p), ransac_inner_cloud)
{
    dlg_polygon_ = dlg_polygon;
    initializeDLGPrior();
}

// ===== initializeDLGPrior =====
// 作用：初始化 DLG（底图）先验——根据 DLG 多边形与原始轮廓的重心/面积比，
//       计算置信度并决定是否启用方向先验和位置先验
void outlineRegular::initializeDLGPrior()
{
    use_dlg_direction_ = false;
    use_dlg_position_ = false;
    dlg_confidence_ = 0.0;
    if (dlg_polygon_.size() < 3 || original_points.size() < 3) return;

    const double dlg_area = polygonArea2D(dlg_polygon_);
    const double contour_area = polygonArea2D(original_points);
    if (dlg_area < 1e-6 || contour_area < 1e-6) return;

    // lambda：计算多边形顶点的 2D 重心
    auto centroid = [](const std::vector<pcl::PointXYZ>& polygon) {
        Eigen::Vector2d center(0.0, 0.0);
        for (const auto& p : polygon) center += Eigen::Vector2d(p.x, p.y);
        return center / static_cast<double>(polygon.size());
    };

    const double scale = std::sqrt(dlg_area);
    const double center_ratio = (centroid(original_points) - centroid(dlg_polygon_)).norm() /
        std::max(scale, 1.0);
    const double area_ratio = contour_area / dlg_area;
    if (center_ratio > 0.45 || area_ratio < 0.03 || area_ratio > 1.50) return;

    const double center_score = clampDouble(1.0 - center_ratio / 0.45, 0.0, 1.0);
    const double area_score = clampDouble(1.0 - std::abs(std::min(area_ratio, 1.0) - 1.0), 0.25, 1.0);
    dlg_confidence_ = center_score * area_score;
    use_dlg_direction_ = dlg_confidence_ >= 0.25;

    // DLG usually represents the maximum footprint. Recessed upper floors only inherit direction.
    use_dlg_position_ = use_dlg_direction_ && area_ratio >= 0.72 &&
        area_ratio <= 1.25 && center_ratio <= 0.18;

    std::cerr << "[DLG prior] area_ratio=" << area_ratio
        << " center_ratio=" << center_ratio
        << " confidence=" << dlg_confidence_
        << " position=" << use_dlg_position_ << std::endl;
}

// ===== computeDLGPriorEnergy =====
// 作用：计算假设多边形相对于 DLG 底图的先验能量（点到 DLG 边的截断距离平方均值）
// 参数：hypothesis - 候选假设多边形
// 返回：DLG 先验能量（未启用时返回 0）
double outlineRegular::computeDLGPriorEnergy(const std::vector<pcl::PointXYZ>& hypothesis)
{
    if (!use_dlg_position_ || hypothesis.size() < 3 || dlg_polygon_.size() < 3) return 0.0;

    double sum = 0.0;
    size_t sample_count = 0;
    const double truncation = clampDouble(0.04 * std::sqrt(polygonArea2D(dlg_polygon_)), 0.5, 2.0);
    for (size_t i = 0; i < hypothesis.size(); ++i) {
        const auto& p = hypothesis[i];
        const auto& q = hypothesis[(i + 1) % hypothesis.size()];
        pcl::PointXYZ samples[2] = {
            p,
            pcl::PointXYZ(0.5f * (p.x + q.x), 0.5f * (p.y + q.y), p.z)
        };
        for (const auto& sample : samples) {
            double min_distance = std::numeric_limits<double>::max();
            for (size_t j = 0; j < dlg_polygon_.size(); ++j) {
                min_distance = std::min(min_distance, computePointToSegmentDistance(
                    sample, dlg_polygon_[j], dlg_polygon_[(j + 1) % dlg_polygon_.size()]));
            }
            min_distance = std::min(min_distance, truncation);
            sum += min_distance * min_distance;
            ++sample_count;
        }
    }
    return sample_count == 0 ? 0.0 : dlg_confidence_ * sum / sample_count;
}

// ===== ~outlineRegular (析构) =====
// 作用：析构函数（空实现）
outlineRegular::~outlineRegular()
{

}

// ===== setInterFloorRegularizationContext =====
// 作用：设置跨楼层规则化上下文——提供建筑物主线角度和参考墙面，用于后续墙面对齐优化
// 参数：building_line_angles - 建筑物已知主线角度列表; reference_walls - 参考墙面列表;
//       wall_snap_weight - 墙面吸附权重; max_wall_snap_distance - 最大吸附距离;
//       max_wall_snap_angle_deg - 最大吸附角度偏差（度）
void outlineRegular::setInterFloorRegularizationContext(
    const std::vector<double>& building_line_angles,
    const std::vector<ReferenceWallLine>& reference_walls,
    double wall_snap_weight,
    double max_wall_snap_distance,
    double max_wall_snap_angle_deg)
{
    building_line_angles_ = building_line_angles;
    reference_walls_ = reference_walls;
    inter_floor_wall_snap_weight_ = std::max(0.0, wall_snap_weight);
    inter_floor_max_wall_snap_distance_ = std::max(0.0, max_wall_snap_distance);
    inter_floor_max_wall_snap_angle_ = std::max(0.0, max_wall_snap_angle_deg) * M_PI / 180.0;
}

// ===== regular_Contour =====
// 作用：轮廓规则化的主入口。完整流程：分辨率估计 -> 圆/椭圆检测 ->
//       生成假设 -> 能量选优 -> 共线点移除 -> Ceres 硬约束优化 ->
//       拓扑后处理 -> 越界检测与回退，最终将结果写入 final_points
void outlineRegular::regular_Contour()
{
    // 确保输入点数为足够的点来生成闭合多边形
    if (original_points.empty() || original_points.size() < 3) {
        return;
    }

    final_points->clear();
    best_hypothesis.clear();
    best_energy_hypothesis_.clear();

    pcl::PointCloud<pcl::PointXYZ>::Ptr fitting_cloud = ransac_inner_cloud;
    pcl::PointCloud<pcl::PointXYZ>::Ptr fallback_cloud;
    if (!fitting_cloud || fitting_cloud->empty()) {
        fallback_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
        fallback_cloud->points.assign(original_points.begin(), original_points.end());
        fitting_cloud = fallback_cloud;
    }

    // ---- 计时(定位规则化各阶段耗时) ----
    auto _t_total = std::chrono::steady_clock::now();
    auto _t0 = _t_total;
    std::cerr << "[RegTime] === polygon: original_points=" << original_points.size()
              << " support=" << fitting_cloud->size() << " ===" << std::endl;

    // 计算点云分辨率 (平均点间距)
    XG::modularFunction mfn;
    double resolution = mfn.computeModelResolution(*fitting_cloud);
    if (!std::isfinite(resolution) || resolution <= 1e-6) {
        resolution = makeOutlineTuning(resolution, computeOBBArea(original_points)).resolution;
    }
    const std::vector<PreservedArcSegment> preserved_arcs =
        detectPreservedArcs(original_points, resolution);

    //判断当前轮廓是否近似圆或椭圆
    if (fitting_cloud->size() >= 6) {
        ShapeResult res = ShapeDetector::analyzeCloud(fitting_cloud, 2 * resolution, 0.15);
        if (res.type == ShapeResult::CIRCLE || res.type == ShapeResult::ELLIPSE)
        {
            best_energy_hypothesis_ = original_points;
            for (auto& op : original_points)
            {
                final_points->points.push_back(op);
            }
            return;
        }
    }
    {
        auto _t1 = std::chrono::steady_clock::now();
        std::cerr << "[RegTime]   resolution+shapedetect: "
                  << std::chrono::duration<double, std::milli>(_t1 - _t0).count() << " ms" << std::endl;
        _t0 = _t1;
    }

    // 自动计算 Lambda
    // 这里传入 original_points 用来算面积 A
    double adaptive_lambda = computeAdaptiveLambda(resolution, original_points);

    // 生成多边形假设
    std::vector<std::vector<pcl::PointXYZ>> hypotheses;
    _t0 = std::chrono::steady_clock::now();
    generatePolygonalHypotheses(original_points, hypotheses);
    //generateHypothesesByAreaSubtraction(original_points, hypotheses);
    {
        auto _t1 = std::chrono::steady_clock::now();
        std::cerr << "[RegTime]   generateHypotheses: "
                  << std::chrono::duration<double, std::milli>(_t1 - _t0).count()
                  << " ms (hyps=" << hypotheses.size() << ")" << std::endl;
        _t0 = _t1;
    }
    double min_energy = std::numeric_limits<double>::max();
    
    //saveAllHypotheses(hypotheses);//保存所有假设
    // 遍历假设，找到总能量最小的
    int best_index = -1; // 用于记录最优索引
    size_t skipped_invalid_hypotheses = 0;
    _t0 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < hypotheses.size(); ++i)
    {
        removeDuplicatePoints2D(hypotheses[i], 1e-4f);
        if (hypotheses[i].size() < 3 || !isSimplePolygon2D(hypotheses[i])) {
            ++skipped_invalid_hypotheses;
            continue;
        }

        // 使用 hypotheses[i] 获取当前假设
        double total_energy = computeTotalEnergy(hypotheses[i], original_points, adaptive_lambda);

        if (total_energy < min_energy) {
            min_energy = total_energy;
            best_hypothesis = hypotheses[i];
            best_index = static_cast<int>(i); // 记录当前索引
        }
    }
    {
        auto _t1 = std::chrono::steady_clock::now();
        std::cerr << "[RegTime]   energy loop: "
                  << std::chrono::duration<double, std::milli>(_t1 - _t0).count() << " ms" << std::endl;
        if (skipped_invalid_hypotheses > 0) {
            std::cerr << "[Hypothesis] skipped invalid/self-intersecting candidates="
                << skipped_invalid_hypotheses << std::endl;
        }
        _t0 = _t1;
    }
    // 在循环结束后输出
    if (best_index != -1) {
        std::cout << "最优假设的索引为: " << best_index << "，最小能量为: " << min_energy << std::endl;
        if (refineThinFeatureHypothesis(best_hypothesis, original_points, adaptive_lambda, resolution)) {
            min_energy = computeTotalEnergy(best_hypothesis, original_points, adaptive_lambda, false);
        }
        best_energy_hypothesis_ = best_hypothesis;
    }
    else {
        std::cout << "未找到有效假设。" << std::endl;
    }

    if (best_index == -1 || best_hypothesis.empty()) {
        best_energy_hypothesis_ = original_points;
        for (const auto& op : original_points) {
            final_points->points.push_back(op);
        }
        return;
    }

    if (best_hypothesis.size() == 3)
    {
        std::cout << "best_hypothesis点数为3，将直接用于规则化结果" << std::endl;
        final_points->clear();
        for (auto& op : best_hypothesis) { final_points->points.push_back(op); }
        return;
    }

    // 优化假设转换为轮廓
    if (!best_hypothesis.empty())
    {
        //保存最优假设
        const std::vector<pcl::PointXYZ> pre_cleanup_hypothesis = best_hypothesis;
        pcl::PointCloud<pcl::PointXYZ>::Ptr in_cloud(new pcl::PointCloud<pcl::PointXYZ>);
        for (auto& op : best_hypothesis) { in_cloud->points.push_back(op); }
        //pcl::io::savePCDFileBinaryCompressed("best_hypothesis.pcd", *in_cloud);
        
        // 移除接近 180° 的点
        removeCollinearPoints(best_hypothesis, 20.0);
        {
            const size_t before_micro = best_hypothesis.size();
            OutlineTuning pre_tuning = makeOutlineTuning(
                resolution, polygonArea2D(best_hypothesis));
            removeDuplicatePoints2D(
                best_hypothesis, static_cast<float>(pre_tuning.fine_prune_distance));
            repairTopologyBatch(best_hypothesis, pre_tuning.repair_distance);
            resolveShortEdgeCluster(best_hypothesis,
                clampDouble(pre_tuning.repair_distance * 2.0,
                            pre_tuning.repair_distance, 7.0));
            removeDuplicatePoints2D(
                best_hypothesis, static_cast<float>(pre_tuning.fine_prune_distance));
            removeCollinearPoints(best_hypothesis, 20.0);
            if (best_hypothesis.size() != before_micro) {
                std::cerr << "[TopoPre] before=" << before_micro
                          << " after=" << best_hypothesis.size() << std::endl;
            }
        }
        if (best_hypothesis.size() < 3 || !isSimplePolygon2D(best_hypothesis)) {
            std::cerr << "[Hypothesis] pre-Ceres cleanup created invalid ring; "
                "rollback to VDP hypothesis" << std::endl;
            best_hypothesis = pre_cleanup_hypothesis;
        }

        if (best_hypothesis.size() == 3)
        {
            std::cout << "best_hypothesis点数为3，将直接用于规则化结果" << std::endl;
            final_points->clear();
            for (auto& op : best_hypothesis) { final_points->points.push_back(op); }
            return;
        }

        pcl::PointCloud<pcl::PointXYZ>::Ptr in_cloud1(new pcl::PointCloud<pcl::PointXYZ>);
        for (auto& op : best_hypothesis) { in_cloud1->points.push_back(op); }
        //pcl::io::savePCDFileBinaryCompressed("best_hypothesis_afterremove180.pcd", *in_cloud1);
        /*for (auto& op : best_hypothesis) {
            final_points->points.push_back(op);
        }*/

        /*for (size_t i = 0; i < best_hypothesis.size(); ++i)
        {
            final_points->points.push_back(best_hypothesis[i]);
        }
        return;*/

        // Ceres 优化强制平行和垂直约束
        //std::vector<std::pair<int, int>> perp_pairs, parallel_pairs;
        //computeConstraintPairs(best_hypothesis, perp_pairs, parallel_pairs);
        //ceresOptimize(best_hypothesis, ransac_inner_cloud, perp_pairs, parallel_pairs);
        //优化final_points     
        pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud(new pcl::PointCloud<pcl::PointXYZ>);
        for (auto& op : best_hypothesis) { input_cloud->points.push_back(op); }
        
        std::vector<pcl::PointXYZ> fallback_hypothesis = best_hypothesis;

        const OutlineTuning model_tuning = makeOutlineTuning(
            resolution, polygonArea2D(best_hypothesis));

        // ===== 多方向证据前置 =====
        // 在尝试单方向正则化之前，先在最优假设上检测方向系统证据：
        // 若存在可信的第二方向系统（斜翼、L 形翼等），直接进入多方向
        // 分支，避免"单方向宽松放行"掩盖斜翼被强行正交化的系统性形变
        // (实测案例：IoU=0.69 / mean=4.1m / q90=11.7m 被 relaxed 通道放行)。
        // 检测不到多方向证据的建筑仍走原 SingleFirst 路径，行为不变。
        std::vector<DirectionSystem> direction_systems =
            detectDirectionSystems(best_hypothesis, fitting_cloud, model_tuning);
        bool credible_multi_direction =
            direction_systems.size() >= 2 &&
            hasCredibleMultiDirectionChains(best_hypothesis, direction_systems, model_tuning);
        if (!credible_multi_direction) {
            std::vector<DirectionSystem> relaxed_systems =
                buildMultiDirectionCandidates(best_hypothesis, resolution);
            if (relaxed_systems.size() >= 2 &&
                hasCredibleMultiDirectionChains(best_hypothesis, relaxed_systems, model_tuning)) {
                direction_systems = std::move(relaxed_systems);
                credible_multi_direction = true;
                std::cerr << "[BuildingMode] relaxed multi-direction evidence accepted" << std::endl;
            }
        }

        // Prefer a single orthogonal direction first.  Many mask-derived
        // footprints contain local stair-step or short diagonal artifacts that
        // can make secondary direction detection too eager.  If one-direction
        // regularization still explains the best hypothesis well, keep it and
        // avoid the less stable multi-direction branch.  Buildings with
        // credible multi-direction evidence (checked above) skip this block.
        std::vector<double> single_first_line_angles;
        if (!credible_multi_direction)
        {
            auto _single_t0 = std::chrono::steady_clock::now();
            std::vector<double> single_line_angles = hausdorffMbrLineAngles2D(
                fallback_hypothesis, resolution, 1);
            if (single_line_angles.empty()) {
                single_line_angles = dominantLineAngles2D(fallback_hypothesis, 1);
            }
            if (!single_line_angles.empty()) {
                double wallAngle = 0.0;
                double wallPeakRatio = 0.0;
                std::size_t wallPairCount = 0;
                bool strongWallDirection = false;
                if (has_support_direction_hint_) {
                    wallAngle = support_direction_hint_;
                    wallPeakRatio = support_direction_peak_ratio_;
                    wallPairCount = support_direction_pair_count_;
                    strongWallDirection = wallPairCount >= kSupportDirectionMinPairs &&
                        wallPeakRatio >= kSupportDirectionStrongPeakRatio;
                } else {
                    strongWallDirection = estimateSupportDirection2D(
                        fitting_cloud, wallAngle, wallPeakRatio, wallPairCount);
                }
                const double originalAngle = single_line_angles.front();
                const double difference = foldedAngleDistance90(originalAngle, wallAngle);
                constexpr double kWallDirectionCorrectionThreshold =
                    kSupportDirectionCorrectionDeg * M_PI / 180.0;
                const bool corrected = strongWallDirection &&
                    difference > kWallDirectionCorrectionThreshold;
                if (corrected) single_line_angles.front() = wallAngle;

                // PCA gate on the ORIGINAL outline ring (not the VDP
                // hypothesis): without strong wall evidence a selection far
                // from the outline's overall trend is a hull-diagonal
                // artifact, so snap it back. Wall-corrected selections are
                // exempt (geometry evidence outranks the outline prior).
                bool pcaFallback = false;
                double pcaAngleDeg = 0.0;
                double pcaAxisRatio = 1.0;
                double pcaDeviationDeg = 0.0;
                {
                    pcl::PointCloud<pcl::PointXYZ>::Ptr ringCloud(
                        new pcl::PointCloud<pcl::PointXYZ>);
                    ringCloud->points.assign(
                        original_points.begin(), original_points.end());
                    double ringPcaAngle = 0.0;
                    if (!corrected &&
                        estimatePcaDirection2D(ringCloud, ringPcaAngle, pcaAxisRatio)) {
                        const double gateDeg =
                            pcaAxisRatio >= kPcaDirectionReliableAxisRatio
                                ? kPcaDirectionGateStrongDeg
                                : (pcaAxisRatio >= kPcaDirectionUsableAxisRatio
                                       ? kPcaDirectionGateWeakDeg
                                       : 1e9);
                        pcaDeviationDeg =
                            foldedAngleDistance90(single_line_angles.front(), ringPcaAngle) *
                            180.0 / M_PI;
                        if (pcaDeviationDeg > gateDeg && !strongWallDirection) {
                            single_line_angles.front() = ringPcaAngle;
                            pcaFallback = true;
                        }
                        pcaAngleDeg = ringPcaAngle * 180.0 / M_PI;
                    }
                }
                std::cerr << "[DirectionDecision] mbr_deg=" << originalAngle * 180.0 / M_PI
                          << " wall_deg=" << wallAngle * 180.0 / M_PI
                          << " diff_deg=" << difference * 180.0 / M_PI
                          << " peak_ratio=" << wallPeakRatio
                          << " pairs=" << wallPairCount
                          << " strong=" << (strongWallDirection ? 1 : 0)
                          << " pca_deg=" << pcaAngleDeg
                          << " pca_axis=" << pcaAxisRatio
                          << " pca_dev=" << pcaDeviationDeg
                          << " pca_fallback=" << (pcaFallback ? 1 : 0)
                          << " selected_deg=" << single_line_angles.front() * 180.0 / M_PI
                          << " corrected=" << (corrected ? 1 : 0) << std::endl;
            }
            single_first_line_angles = single_line_angles;
            if (!single_line_angles.empty() && fallback_hypothesis.size() >= 4) {
                std::vector<pcl::PointXYZ> single_ceres = fallback_hypothesis;
                optimizeWithHardConstraints(single_ceres, fitting_cloud, false, single_line_angles);

                pcl::PointCloud<pcl::PointXYZ>::Ptr single_cloud(new pcl::PointCloud<pcl::PointXYZ>);
                for (const auto& p : single_ceres) single_cloud->push_back(p);
                OptimizeFinal_points(single_cloud);

                std::vector<pcl::PointXYZ> single_topology(
                    single_cloud->points.begin(), single_cloud->points.end());
                std::vector<std::vector<pcl::PointXYZ>> single_sources = {
                    single_topology, single_ceres, fallback_hypothesis
                };

                bool single_accepted = false;
                for (const auto& source : single_sources) {
                    std::vector<pcl::PointXYZ> candidate;
                    if (!forceOrthogonalPolygonToAngle(
                            source, single_line_angles.front(), candidate)) {
                        continue;
                    }
                    if (!isStrictOrthogonalToMainAngle(
                            candidate, single_line_angles.front(), 1.0,
                            model_tuning.fine_prune_distance)) {
                        continue;
                    }
                    if (!isSingleDirectionCandidateAcceptable(
                            candidate, fallback_hypothesis, model_tuning)) {
                        continue;
                    }

                    final_points->clear();
                    for (const auto& p : candidate) final_points->push_back(p);
                    best_hypothesis = candidate;
                    single_accepted = true;
                    break;
                }

                auto _single_t1 = std::chrono::steady_clock::now();
                std::cerr << "[SingleFirst] angle_deg="
                    << single_line_angles.front() * 180.0 / M_PI
                    << " accepted=" << (single_accepted ? 1 : 0)
                    << " time="
                    << std::chrono::duration<double, std::milli>(_single_t1 - _single_t0).count()
                    << " ms" << std::endl;

                if (single_accepted) {
                    if (!preserved_arcs.empty() && final_points && final_points->size() >= 3) {
                        std::vector<pcl::PointXYZ> arc_candidate(
                            final_points->points.begin(), final_points->points.end());
                        if (restorePreservedArcs(
                                arc_candidate, preserved_arcs, resolution) &&
                            isRegularizedPolygonAcceptable(arc_candidate, fallback_hypothesis)) {
                            final_points->clear();
                            for (const auto& p : arc_candidate) final_points->push_back(p);
                        }
                    }
                    std::cerr << "[RegTime] TOTAL: "
                        << std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - _t_total).count()
                        << " ms" << std::endl;
                    return;
                }
            }
        }

        // direction_systems / credible_multi_direction 已在 SingleFirst 之前
        // 完成证据前置检测（见上方"多方向证据前置"块），此处直接使用。
        std::vector<double> preferred_line_angles;
        const bool allow_diagonal_edges = credible_multi_direction;
        if (!allow_diagonal_edges && !single_first_line_angles.empty()) {
            preferred_line_angles.push_back(single_first_line_angles.front());
            std::cerr << "[SingleDirection] locked angle_deg="
                << single_first_line_angles.front() * 180.0 / M_PI << std::endl;
        }
        else {
            preferred_line_angles.reserve(direction_systems.size());
            for (const auto& system : direction_systems) {
                preferred_line_angles.push_back(system.angle);
            }
        }
        std::cerr << "[BuildingMode] "
            << (allow_diagonal_edges ? "AllowDiagonal" : "StrictOrthogonal")
            << std::endl;

        if (allow_diagonal_edges) {
            std::vector<pcl::PointXYZ> local_chain_candidate;
            if (regularizeByLocalChains(
                    best_hypothesis, direction_systems, model_tuning, local_chain_candidate) &&
                isRegularizedPolygonAcceptable(local_chain_candidate, fallback_hypothesis)) {
                std::cerr << "[LocalChain] accepted as final candidate before global Ceres"
                    << std::endl;
                best_hypothesis = local_chain_candidate;
                final_points->clear();
                for (const auto& p : best_hypothesis) final_points->push_back(p);
                auto _t_end = std::chrono::steady_clock::now();
                std::cerr << "[RegTime] TOTAL: "
                    << std::chrono::duration<double, std::milli>(_t_end - _t_total).count()
                    << " ms" << std::endl;
                return;
            }
        }

        _t0 = std::chrono::steady_clock::now();
        optimizeWithHardConstraints(
            best_hypothesis, fitting_cloud, allow_diagonal_edges, preferred_line_angles);
        {
            auto _t1 = std::chrono::steady_clock::now();
            std::cerr << "[RegTime]   optimizeWithHardConstraints(ceres): "
                      << std::chrono::duration<double, std::milli>(_t1 - _t0).count() << " ms" << std::endl;
            _t0 = _t1;
        }

        pcl::PointCloud<pcl::PointXYZ>::Ptr in_cloud2(new pcl::PointCloud<pcl::PointXYZ>);
        for (auto& op : best_hypothesis) { in_cloud2->points.push_back(op); }
        //pcl::io::savePCDFileBinaryCompressed("best_hypothesis_afterceres.pcd", *in_cloud2);

        // optimizeWithHardConstraints already rebuilds best_hypothesis through line intersections.
        final_points->clear();
        for (const auto& op : best_hypothesis) {
            final_points->points.push_back(op);
        }
        std::vector<pcl::PointXYZ> ceres_candidate(final_points->points.begin(), final_points->points.end());

        _t0 = std::chrono::steady_clock::now();
        OptimizeFinal_points(final_points);
        {
            auto _t1 = std::chrono::steady_clock::now();
            std::cerr << "[RegTime]   OptimizeFinal_points(topo): "
                      << std::chrono::duration<double, std::milli>(_t1 - _t0).count() << " ms" << std::endl;
            _t0 = _t1;
        }

        std::vector<pcl::PointXYZ> topology_candidate(final_points->points.begin(), final_points->points.end());
        if (!allow_diagonal_edges && final_points->size() >= 4) {
            std::vector<pcl::PointXYZ> topology_result(
                final_points->points.begin(), final_points->points.end());
            optimizeWithHardConstraints(topology_result, fitting_cloud, false, preferred_line_angles);
            std::vector<pcl::PointXYZ> strict_result = topology_result;
            final_points->clear();
            for (const auto& p : strict_result) final_points->push_back(p);
            OptimizeFinal_points(final_points);
            topology_candidate.assign(final_points->points.begin(), final_points->points.end());
        }

        if (!allow_diagonal_edges && final_points->size() >= 4) {
            std::vector<pcl::PointXYZ> current_result(
                final_points->points.begin(), final_points->points.end());
            std::vector<pcl::PointXYZ> orthogonal_result;
            const double strict_angle = preferred_line_angles.empty() ?
                dominantLineAngles2D(current_result, 1).front() : preferred_line_angles.front();
            if (forceOrthogonalPolygonToAngle(
                    current_result, strict_angle, orthogonal_result) &&
                isStrictOrthogonalToMainAngle(
                    orthogonal_result, strict_angle, 1.0,
                    model_tuning.fine_prune_distance) &&
                isSingleDirectionCandidateAcceptable(
                    orthogonal_result, fallback_hypothesis, model_tuning)) {
                final_points->clear();
                for (const auto& p : orthogonal_result) final_points->push_back(p);
                topology_candidate = orthogonal_result;
                std::cerr << "[Orthogonal] final strict orthogonal snap accepted" << std::endl;
            }
        }

        if (!allow_diagonal_edges && final_points->size() >= 4) {
            const double strict_main_angle = preferred_line_angles.empty()
                ? dominantLineAngles2D(fallback_hypothesis, 1).front()
                : preferred_line_angles.front();
            std::vector<pcl::PointXYZ> strict_final(
                final_points->points.begin(), final_points->points.end());
            if (!isStrictOrthogonalToMainAngle(
                    strict_final, strict_main_angle, 3.0,
                    model_tuning.fine_prune_distance)) {
                std::vector<std::vector<pcl::PointXYZ>> strict_sources = {
                    strict_final, topology_candidate, ceres_candidate, fallback_hypothesis
                };
                bool strict_recovered = false;
                for (const auto& source : strict_sources) {
                    std::vector<pcl::PointXYZ> candidate;
                    if (!forceOrthogonalPolygonToAngle(
                            source, strict_main_angle, candidate)) {
                        continue;
                    }
                    if (!isStrictOrthogonalToMainAngle(
                            candidate, strict_main_angle, 1.0,
                            model_tuning.fine_prune_distance)) {
                        continue;
                    }
                    if (!isSingleDirectionCandidateAcceptable(
                            candidate, fallback_hypothesis, model_tuning)) {
                        continue;
                    }

                    final_points->clear();
                    for (const auto& p : candidate) final_points->push_back(p);
                    topology_candidate = candidate;
                    strict_recovered = true;
                    std::cerr << "[Orthogonal] recovered strict orthogonal candidate" << std::endl;
                    break;
                }
                if (!strict_recovered) {
                    std::cerr << "[Orthogonal] strict candidate unavailable; rollback to pre-Ceres hypothesis"
                        << std::endl;
                    final_points->clear();
                    best_hypothesis = fallback_hypothesis;
                    for (const auto& p : fallback_hypothesis) final_points->push_back(p);
                }
            }
        }

        auto regularizedCandidateOk = [&](const std::vector<pcl::PointXYZ>& candidate) {
            if (isRegularizedPolygonAcceptable(candidate, fallback_hypothesis)) return true;
            if (!allow_diagonal_edges &&
                isStrictOrthogonalToMainAngle(
                    candidate,
                    preferred_line_angles.empty() ? dominantLineAngles2D(fallback_hypothesis, 1).front()
                                                  : preferred_line_angles.front(),
                    3.0, model_tuning.fine_prune_distance)) {
                return isSingleDirectionCandidateAcceptable(
                    candidate, fallback_hypothesis, model_tuning);
            }
            return false;
        };

        std::vector<pcl::PointXYZ> final_candidate(final_points->points.begin(), final_points->points.end());
        if (!regularizedCandidateOk(final_candidate)) {
            bool accepted = false;
            const double strict_main_angle = (!allow_diagonal_edges && !preferred_line_angles.empty())
                ? preferred_line_angles.front()
                : 0.0;
            const auto strict_ok = [&](const std::vector<pcl::PointXYZ>& candidate) {
                return allow_diagonal_edges || isStrictOrthogonalToMainAngle(
                    candidate, strict_main_angle, 3.0, model_tuning.fine_prune_distance);
            };
            if (strict_ok(topology_candidate) &&
                regularizedCandidateOk(topology_candidate)) {
                std::cerr << "[outlineRegular] fallback to topology candidate" << std::endl;
                final_points->clear();
                for (const auto& p : topology_candidate) final_points->push_back(p);
                accepted = true;
            }
            else if (strict_ok(ceres_candidate) &&
                regularizedCandidateOk(ceres_candidate)) {
                std::cerr << "[outlineRegular] fallback to ceres candidate" << std::endl;
                final_points->clear();
                for (const auto& p : ceres_candidate) final_points->push_back(p);
                accepted = true;
            }
            if (!accepted) {
                std::cerr << "[outlineRegular] final polygon invalid, rollback to pre-Ceres hypothesis" << std::endl;
                final_points->clear();
                best_hypothesis = fallback_hypothesis;
                for (const auto& op : best_hypothesis) {
                    final_points->points.push_back(op);
                }
            }
        }

        // 越界检测：final_points中有点超出ransac内点包围盒+20范围，则回退到best_hypothesis
        if (final_points && !final_points->empty() && !best_hypothesis.empty()) {
            bool out_of_bounds = false;
            for (const auto& pt : final_points->points) {
                if (pt.x < ransac_min_x || pt.x > ransac_max_x ||
                    pt.y < ransac_min_y || pt.y > ransac_max_y) {
                    out_of_bounds = true;
                    std::cerr << "[outlineRegular] final_points越界: ("
                        << pt.x << "," << pt.y << "), bounds=["
                        << ransac_min_x << "," << ransac_max_x << ","
                        << ransac_min_y << "," << ransac_max_y << "]" << std::endl;
                    break;
                }
            }
            if (out_of_bounds) {
                std::cerr << "[outlineRegular] 检测到越界点，回退到pre-Ceres hypothesis" << std::endl;
                final_points->clear();
                best_hypothesis = fallback_hypothesis;
                for (auto& op : best_hypothesis) { final_points->points.push_back(op); }
            }
        }

        if (!preserved_arcs.empty() && final_points && final_points->size() >= 3) {
            std::vector<pcl::PointXYZ> arc_candidate(
                final_points->points.begin(), final_points->points.end());
            if (restorePreservedArcs(arc_candidate, preserved_arcs, resolution) &&
                isRegularizedPolygonAcceptable(arc_candidate, fallback_hypothesis)) {
                final_points->clear();
                for (const auto& p : arc_candidate) final_points->push_back(p);
            }
        }

        //pcl::io::savePCDFileBinaryCompressed("final_points.pcd", *final_points);
        //std::vector<pcl::PointXYZ> output_vertices;
       // regularize_cgal(best_hypothesis, ransac_inner_cloud, output_vertices);
        // 将优化结果添加到 final_points 中
        //if (!output_vertices.empty())
        //{
        // for (const auto& point : output_vertices)
        // final_points->push_back(point);
      // }
       // else
        //{
       // final_points->push_back(point);
       // }
       

    }

    std::cerr << "[RegTime] TOTAL: "
              << std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - _t_total).count()
              << " ms" << std::endl;
}

// ===== distSq =====
// 作用：计算两点间的欧几里得距离平方（避免开方，用于性能敏感的比较）
// 返回：距离平方
double distSq(const pcl::PointXYZ& p1, const pcl::PointXYZ& p2) {
    return std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2);
}

// ===== computePerpendicularDistance =====
// 作用：计算点到线段的垂直距离（面积法），退化时回退为点到点距离
// 参数：p - 查询点; line_start, line_end - 线段两端
// 返回：垂直距离
double computePerpendicularDistance(const pcl::PointXYZ& p, const pcl::PointXYZ& line_start, const pcl::PointXYZ& line_end) {
    double dx = line_end.x - line_start.x;
    double dy = line_end.y - line_start.y;
    double mag = std::sqrt(dx * dx + dy * dy);
    if (mag < 1e-6) return std::sqrt(distSq(p, line_start)); // 退化为点到点

    // 面积法计算距离: Area = 0.5 * base * height => height = 2 * Area / base
    // 叉乘 (x1*y2 - x2*y1)
    double area_x2 = std::abs(p.x * dy - p.y * dx + line_end.x * line_start.y - line_end.y * line_start.x);
    return area_x2 / mag;
}

// ===== generatePolygonalHypotheses =====
// 作用：核心修正——基于论文的 Vertex-Driven Douglas-Peucker (VDP) 实现。
//       从连接最远两点出发，迭代地插入到当前边距离最远的原始点，
//       逐步生成从 3 边到最多 30 边的候选多边形假设序列。
// 参数：original_points - 原始轮廓顶点; hypotheses - 输出的候选假设列表
void outlineRegular::generatePolygonalHypotheses(
    const std::vector<pcl::PointXYZ>& original_points,
    std::vector<std::vector<pcl::PointXYZ>>& hypotheses)
{
    hypotheses.clear();
    size_t total_points = original_points.size();
    if (total_points < 3) return;

    // --- 步骤 1: 初始化 (n=2) ---
    // 论文提到初始形状是连接最远两点的线段 
    // 对于闭合轮廓，这实际上将点云分成了两段半圆

    size_t idx_a = 0;
    size_t idx_b = 0;
    double max_dist_sq = -1.0;

    // 简单寻找最远的两个点 O(N^2)，数据量小可以接受，数据量大可用凸包优化
    for (size_t i = 0; i < total_points; ++i) {
        for (size_t j = i + 1; j < total_points; ++j) {
            double d = distSq(original_points[i], original_points[j]);
            if (d > max_dist_sq) {
                max_dist_sq = d;
                idx_a = i;
                idx_b = j;
            }
        }
    }

    // current_indices 存储当前多边形的顶点在 original_points 中的索引
    // 必须保持有序，以便正确映射回原始点云片段
    // 初始多边形由两个点组成，形成闭合回路: A -> B -> A
    // 为了方便逻辑，我们存成有序列表。注意 VDP 插入时需要保持顺序。
    // 这里我们先构建一个以 idx_a 起始的链表概念

    std::list<size_t> active_indices;
    // 确保顺序：总是从小索引到大索引，处理环绕时特殊处理
    if (idx_a > idx_b) std::swap(idx_a, idx_b);
    active_indices.push_back(idx_a);
    active_indices.push_back(idx_b);

    // --- 步骤 2: 迭代增加顶点 (n=3 到 n=N) ---
    // 论文: "iteratively adding the point with the farthest point-to-edge distance" 

    // 我们至少需要生成到 3 个点才开始保存假设
    // 设置最大迭代次数，防止死循环
    size_t max_vertices = std::min((size_t)30, total_points); // 通常建筑轮廓不需要超过30个点

    while (active_indices.size() < max_vertices) {
        double max_error = -1.0;
        size_t best_point_idx = 0;
        auto best_insert_pos = active_indices.end(); // 插入位置

        // 遍历当前多边形的每一条边
        auto it = active_indices.begin();
        while (it != active_indices.end()) {
            size_t p_start_idx = *it;

            // 获取下一节点，如果是最后一个，则连回第一个 (闭合)
            auto next_it = std::next(it);
            size_t p_end_idx = (next_it == active_indices.end()) ? active_indices.front() : *next_it;

            pcl::PointXYZ p_start = original_points[p_start_idx];
            pcl::PointXYZ p_end = original_points[p_end_idx];

            // 扫描这两个顶点之间的所有原始点
            // 处理环绕情况：如果 p_end_idx < p_start_idx，说明跨越了数组尾部
            size_t curr = (p_start_idx + 1) % total_points;

            while (curr != p_end_idx) {
                double dist = computePerpendicularDistance(original_points[curr], p_start, p_end);

                // 找到这一轮中误差最大的点
                if (dist > max_error) {
                    max_error = dist;
                    best_point_idx = curr;
                    // 记录插入位置：应该插在 it 之后
                    best_insert_pos = (next_it == active_indices.end()) ? active_indices.end() : next_it;
                }
                curr = (curr + 1) % total_points;
            }

            it++;
        }

        // 如果误差已经非常小，或者找不到更远的点（max_error未更新），提前结束
        if (max_error < 0.01) break;

        // 插入最远点
        active_indices.insert(best_insert_pos, best_point_idx);

        // --- 保存假设 (仅当 n >= 3 时) ---
        if (active_indices.size() >= 3) {
            std::vector<pcl::PointXYZ> poly;
            for (auto idx : active_indices) {
                poly.push_back(original_points[idx]);
            }
            removeDuplicatePoints2D(poly, 1e-4f);
            if (poly.size() >= 3 && isSimplePolygon2D(poly)) {
                hypotheses.push_back(poly);
            }
            else {
                std::cerr << "[Hypothesis] skip self-intersecting VDP candidate vertices="
                    << poly.size() << std::endl;
            }
        }
    }
}

// 生成多边形假设（3 边到 n 边）
//void outlineRegular::generatePolygonalHypotheses(const std::vector<pcl::PointXYZ>& origial_points, std::vector<std::vector<pcl::PointXYZ>>& hypotheses)
//{
//    size_t n = origial_points.size();
//    if (n < 3) {
//        std::cerr << "Not enough unique points to form a polygon" << std::endl;
//        return;
//    }
//    // 转换为 Points_dp 格式
//    std::vector<XG::modularFunction::Points_dp> points_dp;
//    for (const auto& p : origial_points) {
//        points_dp.emplace_back(p.x, p.y);
//    }
//    XG::modularFunction mfn;
//    std::vector<XG::modularFunction::Points_dp> p1, p2;
//    std::vector<XG::modularFunction::Points_dp> p11, p22;
//    mfn.splitClosedPolylineByMaxSegment(points_dp, p1, p2); // 将闭合点集分成两个子折线 p1 和 p2
//    // 转换为 PCL 点云
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_p1(new pcl::PointCloud<pcl::PointXYZ>);
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_p11(new pcl::PointCloud<pcl::PointXYZ>);
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_p2(new pcl::PointCloud<pcl::PointXYZ>);
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_p22(new pcl::PointCloud<pcl::PointXYZ>);
//    for (const auto& pt : p1)
//    {
//        pcl::PointXYZ ptt;
//        ptt.x = pt.x;
//        ptt.y = pt.y;
//        ptt.z = 0;
//        cloud_p1->push_back(ptt);
//    }
//    for (const auto& pt : p2)
//    {
//        pcl::PointXYZ ptt;
//        ptt.x = pt.x;
//        ptt.y = pt.y;
//        ptt.z = 0;
//        cloud_p2->push_back(ptt);
//    }
//    cloud_p11 = mfn.dedupPCLCloud2D(cloud_p1, 0.1); // 去重
//    cloud_p22 = mfn.dedupPCLCloud2D(cloud_p2, 0.1);
//    mfn.sortByCoordinate2D(cloud_p11); // 排序
//    mfn.sortByCoordinate2D(cloud_p22);
//    for (const auto& pt : *cloud_p11)
//    {
//        XG::modularFunction::Points_dp p;
//        p.x = pt.x;
//        p.y = pt.y;
//        p11.emplace_back(p);
//    }
//    for (const auto& pt : *cloud_p22)
//    {
//        XG::modularFunction::Points_dp p;
//        p.x = pt.x;
//        p.y = pt.y;
//        p22.emplace_back(p);
//    }
//    cloud_p1->clear();
//    cloud_p2->clear();
//    cloud_p11->clear();
//    cloud_p22->clear();
//    p1.clear();
//    std::vector<XG::modularFunction::Points_dp>().swap(p1); // 清空p1
//    p2.clear();
//    std::vector<XG::modularFunction::Points_dp>().swap(p2);
//    p1 = p11;
//    p2 = p22;
//    p11.clear();
//    std::vector<XG::modularFunction::Points_dp>().swap(p11);
//    p22.clear();
//    std::vector<XG::modularFunction::Points_dp>().swap(p22);
//
//    // 检查 p1 和 p2 的点数，决定哪些分支生成假设
//    bool use_p1 = (p1.size() >= 3);
//    bool use_p2 = (p2.size() >= 3);
//    if (!use_p1 && !use_p2) {
//        std::cerr << "Both p1 and p2 have fewer than 3 points, no hypotheses generated." << std::endl;
//        hypotheses.clear();
//        return;
//    }
//    if (!use_p1) {
//        std::cerr << "p1 has fewer than 3 points, skipping p1 branch." << std::endl;
//    }
//    if (!use_p2) {
//        std::cerr << "p2 has fewer than 3 points, skipping p2 branch." << std::endl;
//    }
//
//    vector<std::vector<XG::modularFunction::Points_dp>> vec_simplified_points, vec2_simplified_points;
//    vector<double> epsilon_1, epsilon_2;
//    double MIN_EPSILON = 2;
//    // 生成多边形假设 3 边到 n 边
//    hypotheses.clear();
//
//    // 只对合格的分支生成简化假设
//    if (use_p1) {
//        for (size_t k = 3; k <= p1.size(); ++k) {
//            std::vector<XG::modularFunction::Points_dp> simplified_points;
//            double epsilon = 0.01; // 初始 epsilon值根据数据的粗糙度调整
//            vertexDrivenDouglasPeucker(p1, epsilon, simplified_points, k);
//            vec_simplified_points.emplace_back(simplified_points);
//            epsilon_1.emplace_back(epsilon);
//            if (epsilon < MIN_EPSILON)
//                break;
//        }
//    }
//    if (use_p2) {
//        for (size_t k = 3; k <= p2.size(); ++k) {
//            std::vector<XG::modularFunction::Points_dp> simplified_points_;
//            double epsilon = 0.01; // 初始 epsilon值根据数据的粗糙度调整
//            vertexDrivenDouglasPeucker(p2, epsilon, simplified_points_, k);
//            vec2_simplified_points.emplace_back(simplified_points_);
//            epsilon_2.emplace_back(epsilon);
//            if (epsilon < MIN_EPSILON)
//                break;
//        }
//    }
//
//    // 标准化简化结果（只对使用的分支）
//    if (use_p1) {
//        juged_sytle(vec_simplified_points);
//    }
//    if (use_p2) {
//        juged_sytle(vec2_simplified_points);
//    }
//
//    // 提取中间顶点与 epsilon，并构建 hypotheses
//    if (use_p1 && use_p2) {
//        // 两者都使用：正常合并
//        std::vector<std::pair<XG::modularFunction::Points_dp, double>> sp_p1;
//        if (!vec_simplified_points.empty()) {
//            size_t last_idx = vec_simplified_points.size() - 1;
//            if (last_idx >= 0 && !vec_simplified_points[last_idx].empty() && vec_simplified_points[last_idx].size() >= 3) {  // last_idx >0 改为 >=0
//                int n1 = static_cast<int>(last_idx);
//                std::cerr << "p1: epsilon_1.size()=" << epsilon_1.size() << ", vec_simplified_points[" << n1 << "].size()=" << vec_simplified_points[n1].size() << std::endl;
//                size_t num_mid_points = vec_simplified_points[n1].size() - 2; // 排除首尾
//                if (epsilon_1.size() >= num_mid_points) {
//                    for (size_t i = 1; i < vec_simplified_points[n1].size() - 1; ++i) {
//                        std::pair<XG::modularFunction::Points_dp, double> p;
//                        p.first = vec_simplified_points[n1][i];
//                        p.second = epsilon_1[i - 1];
//                        sp_p1.emplace_back(p);
//                    }
//                }
//                else {
//                    std::cerr << "Warning: epsilon_1 size mismatch for p1, skipping intermediates." << std::endl;
//                }
//            }
//            else {
//                std::cerr << "Warning: vec_simplified_points invalid (empty or <3 points in last), skipping p1." << std::endl;
//            }
//        }
//        else {
//            std::cerr << "Warning: vec_simplified_points empty, skipping p1." << std::endl;
//        }
//
//        std::vector<std::pair<XG::modularFunction::Points_dp, double>> sp_p2;
//        if (!vec2_simplified_points.empty()) {
//            size_t last_idx = vec2_simplified_points.size() - 1;
//            if (last_idx >= 0 && !vec2_simplified_points[last_idx].empty() && vec2_simplified_points[last_idx].size() >= 3) {
//                int n2 = static_cast<int>(last_idx);
//                std::cerr << "p2: epsilon_2.size()=" << epsilon_2.size() << ", vec2_simplified_points[" << n2 << "].size()=" << vec2_simplified_points[n2].size() << std::endl;
//                size_t num_mid_points = vec2_simplified_points[n2].size() - 2;
//                if (epsilon_2.size() >= num_mid_points) {
//                    for (size_t i = 1; i < vec2_simplified_points[n2].size() - 1; ++i) {
//                        std::pair<XG::modularFunction::Points_dp, double> p;
//                        p.first = vec2_simplified_points[n2][i];
//                        p.second = epsilon_2[i - 1];
//                        sp_p2.emplace_back(p);
//                    }
//                }
//                else {
//                    std::cerr << "Warning: epsilon_2 size mismatch for p2, skipping intermediates." << std::endl;
//                }
//            }
//            else {
//                std::cerr << "Warning: vec2_simplified_points invalid, skipping p2." << std::endl;
//            }
//        }
//        else {
//            std::cerr << "Warning: vec2_simplified_points empty, skipping p2." << std::endl;
//        }
//
//        // 合并并构建 hypotheses
//        std::vector<XG::modularFunction::Points_dp> hp = mergeAndSort(sp_p1, sp_p2);
//        if (!hp.empty()) {
//            // star/end 从 p1 的第一个简化（假设 p1 优先；可调整为 vec_simplified_points[0]）
//            pcl::PointXYZ star_point, end_point;
//            if (!vec_simplified_points.empty() && !vec_simplified_points[0].empty() && vec_simplified_points[0].size() >= 2) {
//                star_point.x = vec_simplified_points[0][0].x;
//                star_point.y = vec_simplified_points[0][0].y;
//                star_point.z = 0;
//                end_point.x = vec_simplified_points[0].back().x;
//                end_point.y = vec_simplified_points[0].back().y;
//                end_point.z = 0;
//            }
//            else {
//                std::cerr << "Warning: Cannot get star/end from p1 simplified, using defaults (0,0)." << std::endl;
//                star_point = pcl::PointXYZ{ 0, 0, 0 };
//                end_point = pcl::PointXYZ{ 0, 0, 0 };
//            }
//            hypotheses.reserve(hp.size());
//            std::vector<pcl::PointXYZ> temp_points;
//            temp_points.reserve(hp.size());
//            for (size_t i = 0; i < hp.size(); ++i) {
//                pcl::PointXYZ pcl_point;
//                pcl_point.x = hp[i].x;
//                pcl_point.y = hp[i].y;
//                pcl_point.z = GraphZ; // 假设 GraphZ 是全局 z 值
//                temp_points.push_back(pcl_point);
//                std::vector<pcl::PointXYZ> temp_points_;
//                temp_points_.emplace_back(star_point);
//                for (int j = 0; j < static_cast<int>(temp_points.size()); ++j) {  // 显式 cast
//                    temp_points_.emplace_back(temp_points[j]);
//                }
//                temp_points_.emplace_back(end_point);
//                hypotheses.push_back(temp_points_);
//            }
//        }
//    }
//    else {
//        // 只有一个分支使用：提取活跃 sp，模拟合并，构建 hypotheses
//        bool is_p1_active = use_p1;
//        auto& active_simplified = is_p1_active ? vec_simplified_points : vec2_simplified_points;
//        auto& active_epsilon = is_p1_active ? epsilon_1 : epsilon_2;
//        std::vector<std::pair<XG::modularFunction::Points_dp, double>> sp_active;
//        if (!active_simplified.empty()) {
//            size_t last_idx = active_simplified.size() - 1;
//            if (last_idx >= 0 && !active_simplified[last_idx].empty() && active_simplified[last_idx].size() >= 3) {
//                int n_active = static_cast<int>(last_idx);
//                std::cerr << (is_p1_active ? "p1" : "p2") << ": active_epsilon.size()=" << active_epsilon.size() << ", active_simplified[" << n_active << "].size()=" << active_simplified[n_active].size() << std::endl;
//                size_t num_mid_points = active_simplified[n_active].size() - 2;
//                if (active_epsilon.size() >= num_mid_points) {
//                    for (size_t i = 1; i < active_simplified[n_active].size() - 1; ++i) {
//                        std::pair<XG::modularFunction::Points_dp, double> p;
//                        p.first = active_simplified[n_active][i];
//                        p.second = active_epsilon[i - 1];
//                        sp_active.emplace_back(p);
//                    }
//                }
//                else {
//                    std::cerr << "Warning: active_epsilon size mismatch, skipping intermediates." << std::endl;
//                }
//            }
//            else {
//                std::cerr << "Warning: active_simplified invalid (empty or <3 points in last), skipping." << std::endl;
//            }
//        }
//        else {
//            std::cerr << "Warning: active_simplified empty, skipping." << std::endl;
//        }
//
//        // 模拟合并（单个 source）
//        std::vector<XG::modularFunction::Points_dp> hp;
//        if (!sp_active.empty()) {
//            std::vector<SortableElement> temp;
//            temp.reserve(sp_active.size());
//            int source = is_p1_active ? 0 : 1;
//            for (size_t i = 0; i < sp_active.size(); ++i) {
//                temp.emplace_back(SortableElement{
//                    sp_active[i].first, // Points_dp 坐标信息
//                    sp_active[i].second, // score
//                    source, // source
//                    static_cast<size_t>(i) // 原始索引
//                    });
//            }
//            std::sort(temp.begin(), temp.end());
//            hp.reserve(temp.size());
//            for (const auto& elem : temp) {
//                hp.push_back(elem.point);
//            }
//        }
//
//        // 定义 star_point 和 end_point 从活跃分支
//        pcl::PointXYZ star_point, end_point;
//        if (!active_simplified.empty() && !active_simplified[0].empty() && active_simplified[0].size() >= 2) {
//            star_point.x = active_simplified[0][0].x;
//            star_point.y = active_simplified[0][0].y;
//            star_point.z = 0;
//            end_point.x = active_simplified[0].back().x;
//            end_point.y = active_simplified[0].back().y;
//            end_point.z = 0;
//        }
//        else {
//            std::cerr << "Warning: Cannot get star/end from active simplified, using defaults (0,0)." << std::endl;
//            star_point = pcl::PointXYZ{ 0, 0, 0 };
//            end_point = pcl::PointXYZ{ 0, 0, 0 };
//        }
//
//        // 构建 hypotheses
//        if (!hp.empty()) {
//            hypotheses.reserve(hp.size());
//            std::vector<pcl::PointXYZ> temp_points;
//            temp_points.reserve(hp.size());
//            for (size_t i = 0; i < hp.size(); ++i) {
//                pcl::PointXYZ pcl_point;
//                pcl_point.x = hp[i].x;
//                pcl_point.y = hp[i].y;
//                pcl_point.z = GraphZ; // 假设 GraphZ 是全局 z 值
//                temp_points.push_back(pcl_point);
//                std::vector<pcl::PointXYZ> temp_points_;
//                temp_points_.emplace_back(star_point);
//                for (int j = 0; j < static_cast<int>(temp_points.size()); ++j) {
//                    temp_points_.emplace_back(temp_points[j]);
//                }
//                temp_points_.emplace_back(end_point);
//                hypotheses.push_back(temp_points_);
//            }
//        }
//    }
//}

// ===== mergeAndSort =====
// 作用：将两段子折线的简化中间顶点按自定义规则（score 降序 > 来源优先 > 索引）
//       合并并排序，用于构建候选假设的顶点序列
// 参数：sp_p1, sp_p2 - 两段子折线的(顶点, 分数)对
// 返回：排序后的顶点列表
std::vector<XG::modularFunction::Points_dp> outlineRegular::mergeAndSort(
    const std::vector<std::pair<XG::modularFunction::Points_dp, double>>& sp_p1,
    const std::vector<std::pair<XG::modularFunction::Points_dp, double>>& sp_p2) {
    // 预处理存储容器：动态调整权重和有效性的关键
    std::vector<SortableElement> temp;
    temp.reserve(sp_p1.size() + sp_p2.size());
    // 添加sp_p1元素，source=0 标记来源
    for (size_t i = 0; i < sp_p1.size(); ++i) {
        temp.emplace_back(SortableElement{
            sp_p1[i].first, // Points_dp 坐标信息
            sp_p1[i].second, // score
            0, // source=sp_p1
            i // 原始索引
            });
    }
    // 添加sp_p2元素，source=1 标记来源
    for (size_t i = 0; i < sp_p2.size(); ++i) {
        temp.emplace_back(SortableElement{
            sp_p2[i].first,
            sp_p2[i].second,
            1, // source=sp_p2
            i
            });
    }
    // 一次排序std::sort O(N log N) 时间复杂度
    std::sort(temp.begin(), temp.end());
    // 取出排序后的预处理存储容器
    std::vector<XG::modularFunction::Points_dp> hp;
    hp.reserve(temp.size());
    for (const auto& elem : temp) {
        hp.push_back(elem.point);
    }
    return hp;
}

// ===== vertexDrivenDouglasPeucker =====
// 作用：修改后的 VDP 算法——通过二分法调整 epsilon，使 Douglas-Peucker 简化后的
//       顶点数恰好等于目标顶点数 target_vertices
// 参数：points - 输入折线; epsilon - 输出使用的距离阈值; out - 简化结果; target_vertices - 目标顶点数
void outlineRegular::vertexDrivenDouglasPeucker(const std::vector<XG::modularFunction::Points_dp>& points, double& epsilon,
    std::vector<XG::modularFunction::Points_dp>& out, size_t target_vertices)
{
    if (points.size() < 2) {
        throw std::invalid_argument("Not enough points to simplify");
    }
    XG::modularFunction mfn;
    std::vector<XG::modularFunction::Points_dp> simplified_points;
    // 二分法调整 epsilon 值
    double min_epsilon = 0.0, max_epsilon = 1000.0;
    size_t max_iterations = 100;
    size_t iteration = 0;
    std::vector<XG::modularFunction::Points_dp> points_;
    points_ = points;

    while (iteration < max_iterations) {
        simplified_points.clear();
        mfn.douglasPeucker(points_, epsilon, simplified_points);
        size_t current_vertices = simplified_points.size();
        if (current_vertices == target_vertices) {
            out = simplified_points;
            //cerr << epsilon << endl;
            return;
        }
        else if (current_vertices > target_vertices)//简化边数太多，增大 epsilon
        {
            min_epsilon = epsilon;
            epsilon = (epsilon + max_epsilon) / 2.0;
        }
        else//简化过度，边数太少，减小 epsilon
        {
            max_epsilon = epsilon;
            epsilon = (epsilon + min_epsilon) / 2.0;
        }
        iteration++;
    }

    out = simplified_points;
}

// ===== calculateTriangleArea =====
// 作用：计算三个 2D 点构成的三角形面积（绝对值）
// 返回：三角形面积
double calculateTriangleArea(const pcl::PointXYZ& p1, const pcl::PointXYZ& p2, const pcl::PointXYZ& p3) {
    return 0.5 * std::abs(p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y));
}

// ===== generateHypothesesByAreaSubtraction =====
// 作用：基于最小堆的点删除法生成多边形假设序列。
//       用双向链表维护拓扑，每次删除面积贡献最小的顶点，逐步从 n 边简化到 3 边。
// 参数：original_points - 原始轮廓顶点; hypotheses - 输出的假设列表（从 3 边到 n 边）
void outlineRegular::generateHypothesesByAreaSubtraction(
    const std::vector<pcl::PointXYZ>& original_points,
    std::vector<std::vector<pcl::PointXYZ>>& hypotheses)
{
    int n = static_cast<int>(original_points.size());
    if (n < 3) return;

    hypotheses.clear();

    // 1. 初始化双向链表和顶点数据
    std::vector<VertexNode> nodes(n);
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;

    for (int i = 0; i < n; ++i) {
        nodes[i].id = i;
        nodes[i].prev = (i - 1 + n) % n;
        nodes[i].next = (i + 1) % n;
        nodes[i].is_deleted = false;
    }

    // 2. 计算初始面积并放入最小堆
    for (int i = 0; i < n; ++i) {
        nodes[i].area = calculateTriangleArea(
            original_points[nodes[i].prev],
            original_points[i],
            original_points[nodes[i].next]);
        pq.push({ nodes[i].area, i });
    }

    int current_vertex_count = n;

    // 3. 循环删除点，直到剩下3个顶点
    while (current_vertex_count >= 3) {
        // 保存当前状态作为一个多边形假设
        std::vector<pcl::PointXYZ> current_poly;
        int start_node = 0;
        while (nodes[start_node].is_deleted) start_node++; // 找到第一个没被删的点

        int curr = start_node;
        do {
            current_poly.push_back(original_points[curr]);
            curr = nodes[curr].next;
        } while (curr != start_node);

        hypotheses.push_back(current_poly);

        if (current_vertex_count == 3) break; // 已经是最简多边形

        // 查找堆顶，执行删除
        int to_delete = -1;
        while (!pq.empty()) {
            auto top = pq.top();
            pq.pop();
            // 懒惰删除检查：如果堆中的面积与当前节点面积不符，说明是旧数据
            if (!nodes[top.second].is_deleted && std::abs(nodes[top.second].area - top.first) < 1e-9) {
                to_delete = top.second;
                break;
            }
        }

        if (to_delete != -1) {
            // 执行拓扑删除
            nodes[to_delete].is_deleted = true;
            int p = nodes[to_delete].prev;
            int n_node = nodes[to_delete].next;
            nodes[p].next = n_node;
            nodes[n_node].prev = p;

            // 更新邻居的面积（因为它们的邻居变了）
            // lambda：重新计算指定顶点（邻居变化后）的三角形面积并推入堆
            auto update_area = [&](int idx) {
                nodes[idx].area = calculateTriangleArea(
                    original_points[nodes[idx].prev],
                    original_points[idx],
                    original_points[nodes[idx].next]);
                pq.push({ nodes[idx].area, idx });
            };

            update_area(p);
            update_area(n_node);
            current_vertex_count--;
        }
    }

    // 提示：hypotheses[0] 是 n 点，hypotheses.back() 是 3 点
    // 如果需要按 3 到 n 排序，反转一下 vector
    std::reverse(hypotheses.begin(), hypotheses.end());
}

// ===== computePointToSegmentDistance =====
// 作用：计算点到线段的最短距离（含投影截断：垂足落在线段外时取端点距离）
// 参数：point - 查询点; seg_start, seg_end - 线段两端
// 返回：最短距离
double outlineRegular::computePointToSegmentDistance(const pcl::PointXYZ& point, const pcl::PointXYZ& seg_start, const pcl::PointXYZ& seg_end)
{
    double dx = seg_end.x - seg_start.x;
    double dy = seg_end.y - seg_start.y;
    double mag = std::sqrt(dx * dx + dy * dy);
    if (mag < 1e-6) {
        return std::sqrt(std::pow(point.x - seg_start.x, 2) + std::pow(point.y - seg_start.y, 2));
    }
    double u = ((point.x - seg_start.x) * dx + (point.y - seg_start.y) * dy) / (mag * mag);
    if (u < 0.0) {
        return std::sqrt(std::pow(point.x - seg_start.x, 2) + std::pow(point.y - seg_start.y, 2));
    }
    if (u > 1.0) {
        return std::sqrt(std::pow(point.x - seg_end.x, 2) + std::pow(point.y - seg_end.y, 2));
    }
    double px = seg_start.x + u * dx;
    double py = seg_start.y + u * dy;
    return std::sqrt(std::pow(point.x - px, 2) + std::pow(point.y - py, 2));
}

// ===== computeDataEnergy =====
// 作用：计算数据能量——所有原始点到假设多边形最近边的距离平方和
// 参数：hypothesis - 候选多边形; origial_points - 原始点
// 返回：数据能量（越小表示假设对原始点的拟合越好）
double outlineRegular::computeDataEnergy(const std::vector<pcl::PointXYZ>& hypothesis, const std::vector<pcl::PointXYZ>& origial_points)
{
    double data_energy = 0.0;
    for (const auto& point : origial_points) {
        double min_dist = std::numeric_limits<double>::max();
        for (size_t i = 0; i < hypothesis.size(); ++i) {
            const pcl::PointXYZ& seg_start = hypothesis[i];
            const pcl::PointXYZ& seg_end = hypothesis[(i + 1) % hypothesis.size()];
            double dist = computePointToSegmentDistance(point, seg_start, seg_end);
            min_dist = std::min(min_dist, dist);
        }
        data_energy += min_dist * min_dist;
    }
    return data_energy;
}

// ===== computeTotalEnergy =====
// 作用：计算假设多边形的总能量 = 数据能量 + 模型能量 + 规则性能量 + 支撑能量 + DLG先验能量
// 参数：hypothesis - 候选多边形; origial_points - 原始点; lambda - 能量权重
// 返回：总能量（用于假设选优，越小越好）
double outlineRegular::computeTotalEnergy(const std::vector<pcl::PointXYZ>& hypothesis, const std::vector<pcl::PointXYZ>& origial_points, const double lambda, bool verbose)
{
    double data_energy = computeDataEnergy(hypothesis, origial_points)/(2*log(2.0));

    // 计算模型能量
    // 计算边界框面积 A (x/y 范围，忽略 z)
    double A = computeOBBArea(origial_points);
    if (A < 2.0) A = 2.0;  // 防零/退化


    double model_energy = static_cast<double>(hypothesis.size())*log2(A)* lambda;
    double regularity_energy = 2.0 * lambda * log2(A) * regularityPenalty2D(hypothesis, A);
    pcl::PointCloud<pcl::PointXYZ>::Ptr support_cloud = ransac_inner_cloud;
    if (!support_cloud || support_cloud->empty()) {
        support_cloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
        support_cloud->points.assign(origial_points.begin(), origial_points.end());
    }
    OutlineTuning tuning = makeOutlineTuning(0.0, A);
    double support_energy = 3.0 * lambda * log2(A) *
        edgeSupportPenalty2D(hypothesis, support_cloud, tuning.association_distance, A);
	const double dlg_energy = computeDLGPriorEnergy(hypothesis);
    if (verbose) {
        cerr << "data_energy = " << data_energy
            << " model_energy = " << model_energy
            << " regularity_energy = " << regularity_energy
            << " support_energy = " << support_energy
            << " dlg_energy = " << dlg_energy
            << " total = " << data_energy + model_energy + regularity_energy + support_energy + dlg_energy << endl;
    }
    return data_energy + model_energy + regularity_energy + support_energy + dlg_energy;
}

// ===== computeAdaptiveLambda =====
// 作用：根据点云分辨率和包围盒面积，自适应计算能量平衡权重 Lambda
// 参数：resolution - 点云分辨率; points - 轮廓顶点（用于计算面积）
// 返回：钳位到 [1e-5, 0.5] 的 Lambda 值
double outlineRegular::computeAdaptiveLambda(double resolution, const std::vector<pcl::PointXYZ>& points)
{
    // 1. 计算包围盒面积 A (利用你之前写好的函数)
    double A = computeOBBArea(points);
    // 防止 A 过小导致 log2(A) 变成负数或 0
    if (A < 2.0) A = 2.0;
    if (!std::isfinite(resolution) || resolution <= 1e-6) {
        resolution = makeOutlineTuning(resolution, A).resolution;
    }

    // 2. 设定最小支撑点数 k
    // 含义：一个几何特征至少需要多少个激光点支撑才算有效？
    // 机载点云通常较稀疏，建议 k = 3 到 5
    double min_support_points = 4.0;

    // 3. 常数因子 (1 / (2 * ln 2)) ≈ 0.721
    double energy_factor = 1.0 / (2.0 * std::log(2.0));

    // 4. 计算 Lambda
    // 公式： lambda = (Energy_Factor * k * resolution^2) / log2(A)
    double lambda = (energy_factor * min_support_points * std::pow(resolution, 2)) / std::log2(A);

    // 5. 安全边界（防止 lambda 过小或过大）
    // 经验值：lambda 通常在 1e-4 到 0.1 之间（取决于单位是米）
    lambda = clampDouble(lambda, 1e-5, 0.5);

    std::cerr << "[Auto-Lambda] Resolution: " << resolution
        << "m, Area: " << A
        << "m2, Calculated Lambda: " << lambda << std::endl;

    return lambda;
}

// ===== computeOBBArea =====
// 作用：计算点集的最小外接矩形（OBB）面积
// 参数：points - 输入顶点
// 返回：最小外接矩形面积；点数<3 时返回 0
double outlineRegular::computeOBBArea(const std::vector<pcl::PointXYZ>& points) {
    if (points.size() < 3) return 0.0;

    // Define local 2D point struct
    struct Point2D {
        float x, y;
        Point2D(float a = 0.0f, float b = 0.0f) : x(a), y(b) {}
        Point2D operator-(const Point2D& p) const { return { x - p.x, y - p.y }; }
    };

    // Define local helper functions as lambdas
    auto dot = [](const Point2D& a, const Point2D& b) -> float {
        return a.x * b.x + a.y * b.y;
    };

    auto perp = [](const Point2D& a) -> Point2D {
        return { -a.y, a.x };
    };

    auto norm = [&](const Point2D& a) -> float {
        return std::sqrt(dot(a, a));
    };

    // Directly use x, y as 2D coordinates (assuming z=0 or ignored)
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2d(new pcl::PointCloud<pcl::PointXYZ>);
    for (const auto& pt : points) {
        cloud2d->push_back(pcl::PointXYZ(pt.x, pt.y, 0.0f));
    }

    // Compute 2D convex hull
    pcl::PointCloud<pcl::PointXYZ>::Ptr hull_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::ConvexHull<pcl::PointXYZ> chull;
    chull.setInputCloud(cloud2d);
    chull.setDimension(2);  // Enforce 2D computation
    chull.reconstruct(*hull_cloud);

    std::vector<Point2D> hull2d;
    for (const auto& pt : *hull_cloud) {
        hull2d.emplace_back(pt.x, pt.y);
    }

    size_t n = hull2d.size();
    if (n < 3) return 0.0;

    // Ensure counterclockwise order
    float signed_area = 0.0f;
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        signed_area += hull2d[i].x * hull2d[j].y - hull2d[j].x * hull2d[i].y;
    }
    if (signed_area < 0.0f) {
        std::reverse(hull2d.begin(), hull2d.end());
    }

    // Compute minimum area bounding rectangle using O(n^2) method
    double min_area = std::numeric_limits<double>::max();
    for (size_t i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
        Point2D orig = hull2d[i0];
        Point2D diff = hull2d[i1] - orig;
        float len = norm(diff);
        if (len < 1e-6f) continue;

        Point2D U0 = { diff.x / len, diff.y / len };
        Point2D U1 = perp(U0);

        float min0 = 0.0f, max0 = 0.0f, max1 = 0.0f;
        for (size_t j = 0; j < n; ++j) {
            Point2D D = hull2d[j] - orig;
            float d0 = dot(U0, D);
            if (d0 < min0) min0 = d0;
            if (d0 > max0) max0 = d0;
            float d1 = dot(U1, D);
            if (d1 > max1) max1 = d1;
        }
        double area = static_cast<double>(max0 - min0) * static_cast<double>(max1);
        if (area < min_area) min_area = area;
    }

    return min_area;
}

bool outlineRegular::refineThinFeatureHypothesis(
    std::vector<pcl::PointXYZ>& hypothesis,
    const std::vector<pcl::PointXYZ>& original_points,
    double lambda,
    double resolution)
{
    if (hypothesis.size() < 5 || original_points.size() < 5) return false;

    const double area = std::max(polygonArea2D(hypothesis), 1.0);
    const double scale = std::sqrt(area);
    const double short_limit = clampDouble(std::max(4.0 * resolution, 0.035 * scale), 0.8, 4.0);
    const double max_area_delta = 0.08;
    const double min_iou = 0.88;
    const double current_energy = computeTotalEnergy(hypothesis, original_points, lambda, false);

    std::vector<pcl::PointXYZ> best = hypothesis;
    double best_energy = current_energy;
    std::string best_action;

    auto tryCandidate = [&](std::vector<pcl::PointXYZ> candidate, const std::string& action) {
        removeDuplicatePoints2D(candidate, 1e-4f);
        if (candidate.size() < 4 || !isSimplePolygon2D(candidate)) return;

        const double candidate_area = polygonArea2D(candidate);
        const double area_delta = std::abs(candidate_area - area) / std::max(area, 1e-6);
        if (area_delta > max_area_delta) return;

        const double iou = polygonIoU2D(candidate, hypothesis);
        if (iou < min_iou) return;

        // Bias slightly toward simpler repaired rings, but never accept a large data loss.
        const double energy = computeTotalEnergy(candidate, original_points, lambda, false);
        const double complexity_bonus = 0.15 * std::max<int>(0,
            static_cast<int>(hypothesis.size()) - static_cast<int>(candidate.size()));
        if (energy - complexity_bonus < best_energy + 0.05) {
            best = std::move(candidate);
            best_energy = energy;
            best_action = action;
        }
    };

    const size_t n = hypothesis.size();
    for (size_t i = 0; i < n; ++i) {
        const size_t prev = (i + n - 1) % n;
        const size_t next = (i + 1) % n;
        const auto& a = hypothesis[prev];
        const auto& b = hypothesis[i];
        const auto& c = hypothesis[next];
        const double len_ab = std::hypot(b.x - a.x, b.y - a.y);
        const double len_bc = std::hypot(c.x - b.x, c.y - b.y);
        const double chord = std::hypot(c.x - a.x, c.y - a.y);
        if (len_ab > short_limit && len_bc > short_limit) continue;
        if (chord < 1e-6 || chord > short_limit * 5.0) continue;

        const double detour_ratio = (len_ab + len_bc) / chord;
        const double angle = vertexAngleRad(a, b, c);
        if (detour_ratio < 1.12 && angle > 35.0 * M_PI / 180.0) continue;

        std::vector<pcl::PointXYZ> deleted;
        deleted.reserve(n - 1);
        for (size_t k = 0; k < n; ++k) {
            if (k != i) deleted.push_back(hypothesis[k]);
        }
        tryCandidate(std::move(deleted), "delete-apex");

        // Rectangularize a one-apex thin feature by projecting the apex to the
        // orthogonal system of the chord. This gives Ceres a right-angle corner
        // instead of a single spike when the geometry supports keeping it.
        Eigen::Vector2d u(c.x - a.x, c.y - a.y);
        const double u_norm = u.norm();
        if (u_norm < 1e-6) continue;
        u /= u_norm;
        Eigen::Vector2d v(-u.y(), u.x());
        Eigen::Vector2d av(a.x, a.y);
        Eigen::Vector2d bv(b.x, b.y);
        const double along = (bv - av).dot(u);
        const double depth = (bv - av).dot(v);
        if (std::abs(depth) < std::max(0.5, 2.0 * resolution) ||
            std::abs(depth) > short_limit * 3.0 ||
            along < short_limit * 0.25 || along > chord - short_limit * 0.25) {
            continue;
        }

        Eigen::Vector2d p1 = av + along * u;
        Eigen::Vector2d p0 = p1 + depth * v;
        std::vector<pcl::PointXYZ> rectangularized;
        rectangularized.reserve(n + 1);
        for (size_t k = 0; k < n; ++k) {
            if (k == i) {
                pcl::PointXYZ q0 = b;
                q0.x = static_cast<float>(p0.x());
                q0.y = static_cast<float>(p0.y());
                pcl::PointXYZ q1 = b;
                q1.x = static_cast<float>(p1.x());
                q1.y = static_cast<float>(p1.y());
                rectangularized.push_back(q0);
                rectangularized.push_back(q1);
            }
            else {
                rectangularized.push_back(hypothesis[k]);
            }
        }
        tryCandidate(std::move(rectangularized), "rectangularize-apex");
    }

    if (!best_action.empty()) {
        std::cerr << "[ThinFeature] accepted " << best_action
                  << " vertices " << hypothesis.size() << " -> " << best.size()
                  << " energy " << current_energy << " -> " << best_energy << std::endl;
        hypothesis = std::move(best);
        return true;
    }
    return false;
}

// ===== saveAllHypotheses =====
// 作用：将所有候选假设多边形逐个保存为 PCD 文件（调试用）
// 参数：hypotheses - 候选假设列表
void outlineRegular::saveAllHypotheses(const std::vector<std::vector<pcl::PointXYZ>>& hypotheses)
{
    if (hypotheses.empty()) {
        std::cout << "No hypotheses to save." << std::endl;
        return;
    }

    for (size_t i = 0; i < hypotheses.size(); ++i) {
        const auto& current_poly = hypotheses[i];

        // 1. 构造文件名：例如 hypothesis_0_v12.pcd (表示索引0，有12个顶点)
        std::string filename = "hypothesis_" + std::to_string(i) +
            "_v" + std::to_string(current_poly.size()) + ".pcd";

        // 2. 转换为 PCL 点云对象
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
        cloud->width = static_cast<uint32_t>(current_poly.size());
        cloud->height = 1;
        cloud->is_dense = true;
        cloud->points.assign(current_poly.begin(), current_poly.end());

        // 3. 保存
        if (pcl::io::savePCDFileBinaryCompressed(filename, *cloud) == -1) {
            std::cerr << "Failed to save: " << filename << std::endl;
        }
        else {
            std::cout << "Saved: " << filename << std::endl;
        }
    }

    std::cout << "Total " << hypotheses.size() << " files saved." << std::endl;
}

// ===== computeAngle2D =====
// 作用：辅助函数——计算三个点的二维夹角 (0 - 180度)，p2 为顶点
// 返回：夹角（度）；退化情况返回 0
double computeAngle2D(const pcl::PointXYZ& p1, const pcl::PointXYZ& p2, const pcl::PointXYZ& p3) {
    double v1_x = p1.x - p2.x;
    double v1_y = p1.y - p2.y;
    double v2_x = p3.x - p2.x;
    double v2_y = p3.y - p2.y;

    double norm1 = std::hypot(v1_x, v1_y);
    double norm2 = std::hypot(v2_x, v2_y);

    if (norm1 < 1e-6 || norm2 < 1e-6) return 0.0; // 退化情况

    // 点积公式: a . b = |a||b|cos(theta)
    double dot = v1_x * v2_x + v1_y * v2_y;
    double cos_theta = dot / (norm1 * norm2);

    // 截断防止数值误差导致 acos 越界
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;

    return std::acos(cos_theta) * 180.0 / M_PI;
}

// ===== removeCollinearPoints =====
// 作用：核心函数——移除多边形中接近 180° 的共线冗余顶点
// 参数：polygon - 多边形顶点（就地修改）; threshold_deg - 容差（度），如 10.0 表示 170~190° 区间内被移除
void outlineRegular::removeCollinearPoints(std::vector<pcl::PointXYZ>& polygon, double threshold_deg) {
    if (polygon.size() < 3) return;

    bool changed = true;
    int max_iterations = 10; // 防止极限情况下的死循环
    int iter = 0;

    while (changed && iter < max_iterations) {
        changed = false;
        if (polygon.size() < 3) break;

        std::vector<pcl::PointXYZ> next_polygon;
        int n = polygon.size();
        std::vector<bool> to_remove(n, false);

        for (int i = 0; i < n; ++i) {
            int prev = (i + n - 1) % n;
            int next = (i + 1) % n;

            // 计算夹角
            double angle = computeAngle2D(polygon[prev], polygon[i], polygon[next]);

            // 如果接近 180 度，标记删除
            if (std::abs(angle - 180.0) < threshold_deg) {
                to_remove[i] = true;
                changed = true;
                // 一旦发现有可以删除的点，这一轮就标记了变化
            }
        }

        // 构造新多边形
        for (int i = 0; i < n; ++i) {
            if (!to_remove[i]) {
                next_polygon.push_back(polygon[i]);
            }
        }

        // 检查剩余点数，至少保留三角形
        if (next_polygon.size() >= 3) {
            polygon = next_polygon;
        }
        else {
            // 如果删除太狠导致不足3个点，停止删除
            changed = false;
        }

        iter++;
    }

}

// ===== regularizeBuildingFootprint =====
// 作用：通过角度检查移除接近 180° 的冗余顶点来简化建筑足迹多边形
// 参数：input_points - 输入顶点; angle_remove_threshold - 角度移除阈值(度);
//       angle_adjust_threshold - 角度调整阈值(度); distance_threshold - 距离阈值
// 返回：简化后的多边形顶点
std::vector<pcl::PointXYZ> outlineRegular::regularizeBuildingFootprint(std::vector<pcl::PointXYZ> input_points, float angle_remove_threshold,
    float angle_adjust_threshold, float distance_threshold)
{
    // 第一步：转换为 2D 点
    std::vector<Eigen::Vector2f> poly2d;
    for (const auto& point : input_points) {
        poly2d.push_back(Eigen::Vector2f(point.x, point.y));
    }
    // 移除闭合重复点
    if (!poly2d.empty() && poly2d.front() == poly2d.back()) {
        poly2d.pop_back();
    }
    int n = poly2d.size();
    if (n < 3) return input_points; // 点数少于 3 个无法简化轮廓
    // 第二步：通过角度检查移除 2D 点
    std::vector<bool> to_keep(n, true);
    for (int i = 0; i < n; i++) {
        int prev = (i - 1 + n) % n;
        int curr = i;
        int next = (i + 1) % n;
        Eigen::Vector2f P_prev = poly2d[prev];
        Eigen::Vector2f P_curr = poly2d[curr];
        Eigen::Vector2f P_next = poly2d[next];
        // 计算内部角度
        Eigen::Vector2f A = P_prev - P_curr; // 从当前点到前一点
        Eigen::Vector2f B = P_next - P_curr; // 从当前点到后一点
        float dot = A.dot(B);
        float normA = A.norm();
        float normB = B.norm();
        if (normA == 0 || normB == 0) continue; // 退化点
        float cos_theta = dot / (normA * normB);
        cos_theta = std::max(-1.0f, std::min(1.0f, cos_theta)); // 限制到 [-1, 1]
        float theta = acos(cos_theta) * 180.0f / M_PI; // [0, 180]
        float cross = A.x() * B.y() - A.y() * B.x();
        float internal_angle;
        if (cross < 0) { // 逆时针外角 >180°，内角 <180°
            internal_angle = theta;
        }
        else if (cross > 0) { // 逆时针外角 <180°，内角 >180°
            internal_angle = 360.0f - theta;
        }
        else {
            if (dot > 0) {
                internal_angle = 0.0f; // 同方向（直线）
            }
            else {
                internal_angle = 180.0f; // 反方向（折返）
            }
        }
        // 判断是否需要移除该点（接近 180°）
        if (std::abs(internal_angle - 180.0f) < angle_remove_threshold /*||
            std::abs(internal_angle) < angle_remove_threshold*/) {
            // 距离检查
            Eigen::Vector2f V = P_next - P_prev;
            Eigen::Vector2f W = P_curr - P_prev;
            float cross_VW = V.x() * W.y() - V.y() * W.x();
            float d = std::abs(cross_VW) / V.norm();
            if (d < distance_threshold) {
                to_keep[i] = false;
            }
        }
    }
    // 构建简化后的多边形
    std::vector<Eigen::Vector2f> simplified_poly;
    for (int i = 0; i < n; i++) {
        if (to_keep[i]) {
            simplified_poly.push_back(poly2d[i]);
        }
    }
    // 构建输出点云
    std::vector<pcl::PointXYZ> output_points;
    for (const auto& p : simplified_poly) {
        output_points.push_back(pcl::PointXYZ(p.x(), p.y(), 0.0f));
    }
    // 确保闭合（可选）
    /*if (!output_points.empty()) {
        output_points.push_back(output_points.front());
    }*/
    return output_points;
}

//void outlineRegular::computeConstraintPairs(const std::vector<pcl::PointXYZ>& poly, std::vector<std::pair<int, int>>& perp, std::vector<std::pair<int, int>>& parallelPairs)
//{
//    size_t n = poly.size();
//    if (n < 2) return;
//
//    perp.clear();
//    parallelPairs.clear();
//    auto isDuplicate = [](const std::vector<std::pair<int, int>>& v, int a, int b) {
//        for (const auto& p : v) {
//            if ((p.first == a && p.second == b) || (p.first == b && p.second == a)) return true;
//        }
//        return false;
//    };
//    for (size_t i = 0; i < n; ++i) {
//        size_t ni = (i + 1) % n;
//        Eigen::Vector2d dir_i(poly[ni].x - poly[i].x, poly[ni].y - poly[i].y);
//        if (dir_i.norm() < 1e-9) continue;
//        dir_i.normalize();
//        // 检查 i 边与其他边 j (j != i)
//        for (size_t j = 0; j < n; ++j) {
//            if (j == i) continue; // 自身
//            // 跳过相邻边（前一个和后一个）
//            size_t prev = (i + n - 1) % n;
//            size_t next = (i + 1) % n;
//            if (j == prev || j == next) {
//                // 相邻边垂直/平行约束在交点处处理，避免循环
//                // 考虑循环中平行/垂直约束
//            }
//            size_t nj = (j + 1) % n;
//            Eigen::Vector2d dir_j(poly[nj].x - poly[j].x, poly[nj].y - poly[j].y);
//            double norm_j = dir_j.norm();
//            if (norm_j < 1e-9) {
//                std::cerr << "Warning: Degenerate edge " << j << " with length " << norm_j << std::endl;
//                continue;
//            }
//            dir_j.normalize();
//            double cosang = std::abs(dir_i.dot(dir_j));
//            cosang = std::max(-1.0, std::min(1.0, cosang));
//            double ang = std::acos(cosang);
//            ang = std::min(ang, M_PI - ang);
//            // 垂直边（接近 90°）
//            if (std::abs(ang - M_PI / 2) < t_threshold) {
//                if (!isDuplicate(perp, i, j)) {
//                    perp.emplace_back(i, j);
//                    std::cout << "Added perpendicular pair: (" << i << ", " << j << ")" << std::endl;
//                }
//            }
//            // 平行边（小角度）
//            if (ang < p_threshold && j != prev && j != next) {
//                if (!isDuplicate(parallelPairs, i, j)) {
//                    parallelPairs.emplace_back(i, j);
//                    std::cout << "Added parallel pair: (" << i << ", " << j << ")" << std::endl;
//                }
//            }
//        }
//    }
//}

//void outlineRegular::ceresOptimize(const std::vector<pcl::PointXYZ>& poly,
//    const pcl::PointCloud<pcl::PointXYZ>::Ptr& points,
//    const std::vector<std::pair<int, int>>& perp,
//    const std::vector<std::pair<int, int>>& parallelPairs)
//{
//    size_t n = poly.size();
//    if (n < 3) return;
//
//    // 1. 初始化直线参数 (保持你原有的逻辑)
//    std::vector<double> params(3 * n);
//    for (size_t i = 0; i < n; ++i) {
//        auto p1 = poly[i], p2 = poly[(i + 1) % n];
//        double a = p1.y - p2.y, b = p2.x - p1.x, c = p1.x * p2.y - p2.x * p1.y;
//        double norm = std::sqrt(a * a + b * b);
//        params[3 * i + 0] = a / norm;
//        params[3 * i + 1] = b / norm;
//        params[3 * i + 2] = c / norm;
//    }
//
//    // 2. 【关键】进行数据关联
//    auto associations = classifyPointsToLines(points, poly);
//
//    ceres::Problem problem;
//
//    // 3. 添加参数块与 Manifold (处理 a^2+b^2=1 约束)
//    for (size_t i = 0; i < n; ++i) {
//        // 如果你的 Ceres 版本较新 (>2.0)，建议使用 SphereManifold 
//        // 这里的实现依然使用原来的 NormResidual 软约束，但在参数块设置上要小心
//        problem.AddParameterBlock(&params[i * 3], 3);
//
//        // 添加归一化约束 (强约束)
//        problem.AddResidualBlock(
//            new ceres::AutoDiffCostFunction<NormResidual, 1, 3>(new NormResidual()),
//            new ceres::ScaledLoss(nullptr, 10000.0, ceres::TAKE_OWNERSHIP), // 极大权重
//            &params[i * 3]
//        );
//    }
//
//    // 4. 添加“点到直线”的数据拟合残差
//    double cauchy_scale = 0.5;
//    for (size_t i = 0; i < n; ++i) {
//        // 只添加属于第 i 条直线的点
//        for (int pt_idx : associations[i].point_indices) {
//            const auto& pt = points->points[pt_idx];
//            problem.AddResidualBlock(
//                new ceres::AutoDiffCostFunction<PointToLineResidual, 1, 3>(
//                    new PointToLineResidual(pt.x, pt.y)),
//                new ceres::CauchyLoss(cauchy_scale),
//                &params[i * 3] // 只关联第 i 条直线
//            );
//        }
//    }
//
//    // 5. 添加垂直约束
//    for (const auto& pp : perp) {
//        problem.AddResidualBlock(
//            new ceres::AutoDiffCostFunction<PerpResidual, 1, 3, 3>(new PerpResidual()),
//            new ceres::ScaledLoss(nullptr, 500.0, ceres::TAKE_OWNERSHIP), // 权重根据需要调整
//            &params[pp.first * 3], &params[pp.second * 3]);
//    }
//
//    // 6. 添加平行约束
//    for (const auto& pp : parallelPairs) {
//        problem.AddResidualBlock(
//            new ceres::AutoDiffCostFunction<ParallelResidual, 1, 3, 3>(new ParallelResidual()),
//            new ceres::ScaledLoss(nullptr, 500.0, ceres::TAKE_OWNERSHIP),
//            &params[pp.first * 3], &params[pp.second * 3]);
//    }
//
//    // 7. 【关键新增】添加“锚点约束” (Shape Prior)
//    // 防止直线为了满足直角约束而漂移太远
//    for (size_t i = 0; i < n; ++i) {
//        // 让直线尽量穿过初始多边形的两个端点 (p1, p2)
//        // 这样可以限制直线的位置，不仅仅是方向
//        auto p1 = poly[i];
//        auto p2 = poly[(i + 1) % n];
//
//        // 约束：初始端点到优化后直线的距离不应过大
//        problem.AddResidualBlock(
//            new ceres::AutoDiffCostFunction<PointToLineResidual, 1, 3>(new PointToLineResidual(p1.x, p1.y)),
//            new ceres::ScaledLoss(nullptr, 10.0, ceres::TAKE_OWNERSHIP), // 较小权重，允许微调
//            &params[i * 3]
//        );
//        problem.AddResidualBlock(
//            new ceres::AutoDiffCostFunction<PointToLineResidual, 1, 3>(new PointToLineResidual(p2.x, p2.y)),
//            new ceres::ScaledLoss(nullptr, 10.0, ceres::TAKE_OWNERSHIP),
//            &params[i * 3]
//        );
//    }
//
//    // 求解
//    ceres::Solver::Options options;
//    options.linear_solver_type = ceres::DENSE_QR;
//    options.max_num_iterations = 100;
//    ceres::Solver::Summary summary;
//    ceres::Solve(options, &problem, &summary);
//
//    std::cout << summary.BriefReport() << std::endl;
//
//    // 存储结果
//    optimized_params.resize(n);
//    for (size_t i = 0; i < n; ++i) {
//        optimized_params[i] = Eigen::Vector3d(params[3 * i], params[3 * i + 1], params[3 * i + 2]);
//    }
//}

// ===== updateFinalPointsFromLines =====
// 作用：基于优化后的直线参数重建多边形顶点——对相邻两条直线求交得到新顶点。
//       包含平行保护（跳过近平行线）和空间合理性检查（越界拦截），
//       若拓扑坍塌则回退到 best_hypothesis
void outlineRegular::updateFinalPointsFromLines()
{
    final_points->clear();
    size_t n = optimized_params.size();
    if (n < 3) return;

    // --- 步骤 1: 动态计算合法空间范围（防止点“飞走”） ---
    // 基于原始点云的包围盒，并给定一个合理的裕量（例如 20 米）
    // 这样即使 det 稍微大于阈值，如果产生的点在几公里外，也会被拦截
    double min_x = 1e10, max_x = -1e10, min_y = 1e10, max_y = -1e10;
    bool has_valid_bounds = false;
    if (!original_points.empty()) {
        for (const auto& p : original_points) {
            if (p.x < min_x) min_x = p.x; if (p.x > max_x) max_x = p.x;
            if (p.y < min_y) min_y = p.y; if (p.y > max_y) max_y = p.y;
        }
        double margin = 10.0; // 裕量，可根据建筑尺度调整
        min_x -= margin; max_x += margin;
        min_y -= margin; max_y += margin;
        has_valid_bounds = true;
    }

    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;

        // 直线方程参数: ax + by + c = 0
        double a1 = optimized_params[i][0], b1 = optimized_params[i][1], c1 = optimized_params[i][2];
        double a2 = optimized_params[j][0], b2 = optimized_params[j][1], c2 = optimized_params[j][2];

        // 行列式 det = a1*b2 - a2*b1 = sin(theta2 - theta1)
        double det = a1 * b2 - a2 * b1;

        // --- 步骤 2: 增强平行保护 ---
        // 将阈值从 1e-5 提高到 1e-3 (对应夹角约 0.05度)。
        // 当两条线几乎平行时，det 趋于 0，分母爆炸会导致坐标飞出。
        if (std::abs(det) < 1e-3) {
            std::cerr << "[UpdatePoints] Warning: Skipping near-parallel lines " << i << " and " << j
                << " (det=" << det << ")" << std::endl;
            // 跳过此交点，相当于在拓扑上合并了这两条边
            continue;
        }

        pcl::PointXYZ intersection_point;
        intersection_point.x = static_cast<float>((b1 * c2 - b2 * c1) / det);
        intersection_point.y = static_cast<float>((a2 * c1 - a1 * c2) / det);
        intersection_point.z = 0.0f;

        // --- 步骤 3: 空间合理性检查 ---
        // 核心修复：检查计算出的交点是否在原始建筑区域的合理范围内
        if (has_valid_bounds) {
            if (intersection_point.x < min_x || intersection_point.x > max_x ||
                intersection_point.y < min_y || intersection_point.y > max_y) {
                std::cerr << "[UpdatePoints] Critical: Intersection (" << i << "," << j
                    << ") at (" << intersection_point.x << "," << intersection_point.y
                    << ") is out of bounds. Skipping." << std::endl;
                continue;
            }
        }

        final_points->push_back(intersection_point);
    }

    // --- 步骤 4: 结果一致性检查 ---
    if (final_points->size() < 3) {
        std::cerr << "[UpdatePoints] Error: Polygon collapsed to " << final_points->size()
            << " points. Using original hypothesis as fallback." << std::endl;
        // 兜底方案：如果优化后的拓扑完全坍塌，可以考虑恢复为 best_hypothesis
        final_points->clear();
        if (final_points->empty() && !best_hypothesis.empty()) {
            for (const auto& p : best_hypothesis) final_points->push_back(p);
        }
    }
}

// ===== VDPEnergySimplify =====
// 作用：基于 VDP 假设 + 能量最小化的独立简化入口（不经过 Ceres 优化）
// 参数：input_vertices - 输入顶点; output_vertices - 输出最优假设; resolution - 点云分辨率
void outlineRegular::VDPEnergySimplify(const std::vector<pcl::PointXYZ>& input_vertices,
    std::vector<pcl::PointXYZ>& output_vertices, double& resolution)
{
    // 自动计算 Lambda
    // 这里传入 original_points 用来算面积 A
    double adaptive_lambda = computeAdaptiveLambda(resolution, input_vertices);

    // 生成多边形假设
    std::vector<std::vector<pcl::PointXYZ>> hypotheses;
    std::vector<pcl::PointXYZ> best_hypothesis_;
    generatePolygonalHypotheses(input_vertices, hypotheses);
    double min_energy = std::numeric_limits<double>::max();

    // 遍历假设，找到总能量最小的
    int best_index = -1; // 用于记录最优索引
    for (size_t i = 0; i < hypotheses.size(); ++i)
    {
        // 使用 hypotheses[i] 获取当前假设
        double total_energy = computeTotalEnergy(hypotheses[i], input_vertices, adaptive_lambda);

        if (total_energy < min_energy) {
            min_energy = total_energy;
            best_hypothesis_ = hypotheses[i];
            best_index = static_cast<int>(i); // 记录当前索引
        }
    }

    output_vertices= best_hypothesis_;
}


// ===== regularize_cgal =====
// 作用：静态方法——使用 CGAL Polyline_simplification_2 对输入顶点做距离容差简化
// 参数：input_vertices - 输入顶点; output_vertices - 简化后顶点;
//       min_length - 最小边长; max_angle_deg - 最大角度偏差; max_offset - 最大偏移
void outlineRegular::regularize_cgal(const std::vector<pcl::PointXYZ>& input_vertices,
    std::vector<pcl::PointXYZ>& output_vertices,
    double min_length,
    double max_angle_deg,
    double max_offset
    // 可以考虑在函数参数中增加一个 double simplify_tolerance
)
{
    
    if (input_vertices.size() < 3) {
        std::cerr << "Input vertices too few for contour!" << std::endl;
        return;
    }
    double current_z = input_vertices[0].z;

    // 假设 Kernel 和 Point_2 已经在你的类中被 typedef 过了
    typedef CGAL::Polygon_2<Kernel> Polygon_2;
    namespace PS = CGAL::Polyline_simplification_2;

    // 步骤1: 转换为 CGAL Polygon_2 点集，用于简化
    Polygon_2 poly;
    for (const auto& pt : input_vertices) {
        poly.push_back(Kernel::Point_2(pt.x, pt.y));
    }

    // ==========================================================
    // 步骤1.5: 轮廓简化 (Polyline Simplification)
    // ==========================================================
    // 设定简化的容差（例如 0.1 米），距离在这个范围内的冗余点会被剔除
    // 这里作为示例，取 max_offset 的一半，你可以根据实际需求调整
    double simplify_tolerance = max_offset * 0.5;

    // CGAL 的距离代价是基于平方距离的
    PS::Stop_above_cost_threshold stop(simplify_tolerance * simplify_tolerance);
    PS::Squared_distance_cost cost;

    // 执行简化（simplify 会直接返回一个精简后的 Polygon_2）
    poly = PS::simplify(poly, cost, stop);

    std::cout << "Vertices after simplification: " << poly.size() << std::endl;

    if (poly.size() < 3) {
        std::cerr << "Polygon collapsed after simplification!" << std::endl;
        output_vertices.clear();
        output_vertices=input_vertices;
        return;
    }

    //// 将简化后的多边形提取到 Contour 中，供规则化使用
    //Contour contour;
    //for (auto it = poly.vertices_begin(); it != poly.vertices_end(); ++it) {
    //    contour.emplace_back(it->x(), it->y());
    //}
    //// ==========================================================

    //// 步骤2: 方向估计，使用 Multiple_directions_2 支持多方向
    //using Contour_directions = CGAL::Shape_regularization::Contours::Multiple_directions_2<Kernel, Contour>;
    //const bool is_closed = true;
    //Contour_directions directions(
    //    contour, is_closed,
    //    CGAL::parameters::minimum_length(min_length).maximum_angle(max_angle_deg)
    //);
    //std::cout << "Estimated number of principal directions: " << directions.number_of_directions() << std::endl;

    //// 步骤3: 轮廓正则化
    //Contour regularized;
    //CGAL::Shape_regularization::Contours::regularize_closed_contour(
    //    contour, directions, std::back_inserter(regularized),
    //    CGAL::parameters::maximum_offset(max_offset)
    //);


    // 步骤5: 转换为 PCL 格式
    output_vertices.clear();
    for (const auto& pt : poly) {
        output_vertices.emplace_back(static_cast<float>(pt.x()), static_cast<float>(pt.y()), current_z);
    }
    std::cout << "Regularized vertices count: " << output_vertices.size() << std::endl;
}

//// 【核心函数】优化求解
//void outlineRegular::optimizeWithHardConstraints(
//    const std::vector<pcl::PointXYZ>& hypothesis_raw,
//    const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud)
//{
//    // 1. 预处理：剔除短边
//    std::vector<pcl::PointXYZ> hypothesis= hypothesis_raw;
//    pruneShortEdges(hypothesis_raw, hypothesis, 1.0); // 阈值 0.5m，可根据点云密度调整
//
//    size_t n = hypothesis.size();
//    if (n < 3) return;
//
//    std::vector<EdgeRegularInfo> edges(n);
//    std::vector<double> initial_ds(n);
//    std::vector<double> initial_thetas(n);
//    std::vector<pcl::PointXYZ> mid_points(n);
//
//    // 2. 初始化边参数
//    for (size_t i = 0; i < n; ++i) {
//        pcl::PointXYZ p1 = hypothesis[i];
//        pcl::PointXYZ p2 = hypothesis[(i + 1) % n];
//
//        double dx = p2.x - p1.x;
//        double dy = p2.y - p1.y;
//        double line_angle = std::atan2(dy, dx);
//        double normal_angle = line_angle - M_PI / 2.0;
//        double d = p1.x * std::cos(normal_angle) + p1.y * std::sin(normal_angle);
//
//        edges[i] = { (int)i, -1, 0.0, d, normal_angle, {} };
//        initial_ds[i] = d;
//        initial_thetas[i] = normal_angle;
//        mid_points[i].x = (p1.x + p2.x) / 2.0;
//        mid_points[i].y = (p1.y + p2.y) / 2.0;
//    }
//
//    // 3. 计算全局主方向 (使用直方图优化)
//    double base_theta = calculateGlobalDominantAngle(hypothesis);
//
//    // 4. 边分类 (Grouping)
//    auto normalize_ang = [](double a) {
//        while (a <= -M_PI) a += 2 * M_PI;
//        while (a > M_PI) a -= 2 * M_PI; 
//        return a;
//    };
//    double tolerance = 40.0 * M_PI / 180.0;//认为45度加减10度的范围是不用纠正的
//
//    for (size_t i = 0; i < n; ++i) {
//        double diff = normalize_ang(edges[i].current_theta - base_theta);
//
//        // 平行检测
//        if (std::abs(diff) < tolerance || std::abs(std::abs(diff) - M_PI) < tolerance) {
//            edges[i].type = 0;
//            edges[i].theta_offset = (std::abs(diff) < tolerance) ? 0.0 : M_PI;//保证直线参数的符号稳定性，防止优化过程中 $d$ 发生正负号“跳变”，提高收敛速度和精度
//        }
//        // 垂直检测
//        else if (std::abs(std::abs(diff) - M_PI / 2.0) < tolerance) {
//            edges[i].type = 1;
//            edges[i].theta_offset = (diff > 0) ? M_PI / 2.0 : -M_PI / 2.0;
//        }
//        else {
//            edges[i].type = -1; // 自由边
//        }
//    }
//
//    // 5. 数据关联 (带拐角遮蔽优化)
//    for (const auto& pt : cloud->points) {
//        double min_dist = 1e9;
//        int best_edge = -1;
//
//        // 寻找最近的边
//        for (size_t i = 0; i < n; ++i) {
//            double d = computePointToSegmentDistance(pt, hypothesis[i], hypothesis[(i + 1) % n]);
//            if (d < min_dist) { min_dist = d; best_edge = i; }
//        }
//
//        // 阈值过滤 + 拐角遮蔽
//        if (best_edge != -1 && min_dist < 1.0) {
//            pcl::PointXYZ p1 = hypothesis[best_edge];
//            pcl::PointXYZ p2 = hypothesis[(best_edge + 1) % n];
//
//            // 计算投影比例 t
//            double lx = p2.x - p1.x;
//            double ly = p2.y - p1.y;
//            double len_sq = lx * lx + ly * ly;
//
//            if (len_sq > 1e-6) {
//                double t = ((pt.x - p1.x) * lx + (pt.y - p1.y) * ly) / len_sq;
//                // 【关键优化】只保留线段中间 80% 的点，避开拐角噪声
//                if (t > 0.05 && t < 0.95) 
//                {
//                    edges[best_edge].associated_points.push_back(pt);
//                }
//            }
//        }
//    }
//
//    // --- Step 4: 构建 Ceres 问题 ---
//    ceres::Problem problem;
//
//    // 参数块 1: 全局基准角度
//    problem.AddParameterBlock(&base_theta, 1);
//
//    for (int i = 0; i < n; ++i) {
//        // 参数块 2: 每条边的 d
//        problem.AddParameterBlock(&edges[i].d, 1);
//
//        // 如果是自由边，它有自己独立的 theta
//        if (edges[i].type == -1) {
//            problem.AddParameterBlock(&edges[i].current_theta, 1);
//        }
//
//        // 添加残差
//        for (const auto& pt : edges[i].associated_points) {
//            if (edges[i].type != -1) {
//                // 受约束边：共享 base_theta
//                ceres::CostFunction* cost = PointToLineDistanceCost::Create(
//                    pt.x, pt.y, edges[i].theta_offset);
//                problem.AddResidualBlock(cost, new ceres::HuberLoss(0.1), &base_theta, &edges[i].d);
//            }
//            else {
//                // 自由边
//                ceres::CostFunction* cost = FreeLineCost::Create(pt.x, pt.y);
//                problem.AddResidualBlock(cost, new ceres::HuberLoss(0.1), &edges[i].current_theta, &edges[i].d);
//            }
//        }
//    }
//
//    // --- Step 7: 添加顶点锚定约束 (修复版) ---
//    // 目的：防止规则化后的轮廓“飞”离原始轮廓太远
//    // 原理：将原始假设的顶点(hypothesis)视为"强约束点"，要求相邻的两条直线都必须尽量穿过这个顶点
//
//    // 权重设置：建议比普通点云的权重(0.1)大很多，比如 10.0 到 50.0
//    // 这样优化器会优先保证形状不散架，再去微调贴合点云
//    double vertex_anchor_weight = 10.0;
//    ceres::LossFunction* vertex_loss = new ceres::ScaledLoss(
//        new ceres::HuberLoss(0.1), // 使用 Huber 防止个别极端坏点拉坏整体
//        vertex_anchor_weight,
//        ceres::TAKE_OWNERSHIP
//    );
//
//    for (size_t i = 0; i < n; ++i) {
//        // 获取原始顶点 P (连接 edge[i] 和 edge[i+1] 的点)
//        // 注意：根据你的拓扑定义，hypothesis[i] 是 edge[i] 的起点，hypothesis[(i+1)%n] 是 edge[i] 的终点
//        // 而这个终点同时也是 edge[(i+1)%n] 的起点。
//        // 所以我们需要约束：edge[i] 经过 p2，且 edge[(i+1)%n] 也经过 p2
//
//        pcl::PointXYZ p_vertex = hypothesis[(i + 1) % n];
//
//        size_t j = (i + 1) % n; // 下一条边
//
//        // --- 约束 1: 当前边 edge[i] 必须靠近顶点 p_vertex ---
//        if (edges[i].type != -1) {
//            // 受约束边
//            ceres::CostFunction* cost = PointToLineDistanceCost::Create(
//                p_vertex.x, p_vertex.y, edges[i].theta_offset);
//            // 注意：这里每个 AddResidualBlock 只用了一次 base_theta，不会崩溃
//            problem.AddResidualBlock(cost, vertex_loss, &base_theta, &edges[i].d);
//        }
//        else {
//            // 自由边
//            ceres::CostFunction* cost = FreeLineCost::Create(p_vertex.x, p_vertex.y);
//            problem.AddResidualBlock(cost, vertex_loss, &edges[i].current_theta, &edges[i].d);
//        }
//
//        // --- 约束 2: 下一条边 edge[j] 必须靠近同一个顶点 p_vertex ---
//        if (edges[j].type != -1) {
//            // 受约束边
//            ceres::CostFunction* cost = PointToLineDistanceCost::Create(
//                p_vertex.x, p_vertex.y, edges[j].theta_offset);
//            problem.AddResidualBlock(cost, vertex_loss, &base_theta, &edges[j].d);
//        }
//        else {
//            // 自由边
//            ceres::CostFunction* cost = FreeLineCost::Create(p_vertex.x, p_vertex.y);
//            problem.AddResidualBlock(cost, vertex_loss, &edges[j].current_theta, &edges[j].d);
//        }
//    }
//
//    // 7. 求解
//    ceres::Solver::Options options;
//    options.linear_solver_type = ceres::DENSE_QR;
//    options.max_num_iterations = 100;
//    options.minimizer_progress_to_stdout = true;
//    ceres::Solver::Summary summary;
//    ceres::Solve(options, &problem, &summary);
//
//    // --- Step 6: 结果写回 optimized_params 以便重建 ---
//    optimized_params.resize(n);
//    for (int i = 0; i < n; ++i) {
//        double theta;
//        if (edges[i].type == -1) theta = edges[i].current_theta;
//        else theta = base_theta + edges[i].theta_offset;
//
//        // 转换回直线一般式 ax + by + c = 0
//        // x*cos + y*sin - d = 0  => a=cos, b=sin, c=-d
//        optimized_params[i] = Eigen::Vector3d(std::cos(theta), std::sin(theta), -edges[i].d);
//
//    }
//}

// ===== optimizeWithHardConstraints (带外层迭代) =====
// 作用：核心函数——带外层 EM 迭代的硬约束优化求解。
//       每轮迭代执行：短边剔除 -> 边参数初始化 -> 多主方向分类 ->
//       数据关联 -> 构建 Ceres 问题(点拟合+中点锚定+正则化+DLG/跨层约束) ->
//       求解 -> 直线求交重建多边形，直到收敛或达最大迭代次数。
// 参数：hypothesis_raw - 初始假设多边形; cloud - 支撑点云
void outlineRegular::optimizeWithHardConstraints(
    const std::vector<pcl::PointXYZ>& hypothesis_raw,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
    bool allow_diagonal_edges,
    const std::vector<double>& preferred_line_angles)
{
    // 初始化当前多边形
    std::vector<pcl::PointXYZ> current_polygon = hypothesis_raw;
    std::vector<pcl::PointXYZ> last_valid_polygon = current_polygon;
    if (!allow_diagonal_edges) {
        simplifySameOrthogonalAxisVertices(current_polygon);
        if (current_polygon.size() < 4) current_polygon = hypothesis_raw;
    }
    if (current_polygon.size() >= 3 && isSimplePolygon2D(current_polygon)) {
        last_valid_polygon = current_polygon;
    }
    double tuning_area = computeOBBArea(hypothesis_raw);
    if (tuning_area < 1e-6) tuning_area = polygonArea2D(hypothesis_raw);

    double tuning_resolution = 0.0;
    if (cloud && !cloud->empty()) {
        XG::modularFunction mfn;
        tuning_resolution = mfn.computeModelResolution(*cloud);
    }
    OutlineTuning tuning = makeOutlineTuning(tuning_resolution, tuning_area);

    // 引入外层 EM 迭代：解决直线平移后产生新短边的问题
    int max_outer_iters = 10; // 迭代 2-3 次通常即可收敛

    for (int iter = 0; iter < max_outer_iters; ++iter) {
        // 【修改】先备份当前完整的多边形
        std::vector<pcl::PointXYZ> backup_polygon = current_polygon;

        // 1. 预处理：剔除短边 (经过上一轮平移后可能挤压出新的短边)
        pruneShortEdges(current_polygon, current_polygon, tuning.prune_distance);

        size_t n = current_polygon.size();
        if (n < 3)
        { 
            // 【修改】如果修剪后退化了，恢复成修剪前的有效多边形，然后退出
            current_polygon = backup_polygon;
            break; // 如果修剪后无法构成多边形，提前终止
        }

        std::vector<EdgeRegularInfo> edges(n);
        std::vector<pcl::PointXYZ> mid_points(n);

        // 2. 初始化边参数
        for (size_t i = 0; i < n; ++i) {
            pcl::PointXYZ p1 = current_polygon[i];
            pcl::PointXYZ p2 = current_polygon[(i + 1) % n];

            double dx = p2.x - p1.x;
            double dy = p2.y - p1.y;
            double line_angle = std::atan2(dy, dx);
            double normal_angle = line_angle - M_PI / 2.0;
            // Hessian Normal Form: x*cos + y*sin - d = 0
            double d = p1.x * std::cos(normal_angle) + p1.y * std::sin(normal_angle);

            edges[i] = { (int)i, -1, 0.0, d, normal_angle, {}, {} };

            // 记录中点，用于后续的“中点锚定约束”
            mid_points[i].x = (p1.x + p2.x) / 2.0;
            mid_points[i].y = (p1.y + p2.y) / 2.0;
            mid_points[i].z = p1.z;
        }

        // 3. 计算全局主方向。复杂建筑可能有多个非正交主方向，分别给一组 base_theta。
        std::vector<double> base_line_angles;
        auto append_angle = [&](double angle, double tolerance_deg) {
            angle = foldedLineAngle90(angle);
            const double tolerance = tolerance_deg * M_PI / 180.0;
            for (double existing : base_line_angles) {
                if (foldedAngleDistance90(existing, angle) < tolerance) {
                    return;
                }
            }
            base_line_angles.push_back(angle);
        };

        for (double angle : preferred_line_angles) {
            append_angle(angle, 8.0);
        }

        const std::vector<double> polygon_angles = dominantLineAngles2D(
            current_polygon, allow_diagonal_edges ? 4 : 1);
        for (double angle : polygon_angles) {
            append_angle(angle, allow_diagonal_edges ? 10.0 : 20.0);
        }

        if (!allow_diagonal_edges && base_line_angles.size() > 1) {
            base_line_angles.resize(1);
        }
        if (allow_diagonal_edges && base_line_angles.size() > 4) {
            base_line_angles.resize(4);
        }

        if (use_dlg_direction_) {
            const std::vector<double> dlg_angles = dominantLineAngles2D(dlg_polygon_, 3);
            std::vector<double> merged_angles;
            const double merge_tolerance = 20.0 * M_PI / 180.0;
            for (double dlg_angle : dlg_angles) {
                merged_angles.push_back(dlg_angle);
            }
            for (double mesh_angle : base_line_angles) {
                bool represented = false;
                for (double dlg_angle : dlg_angles) {
                    if (angleDistanceToOrthogonalSystem(mesh_angle, dlg_angle) < merge_tolerance) {
                        represented = true;
                        break;
                    }
                }
                if (!represented) merged_angles.push_back(mesh_angle);
            }
            if (!merged_angles.empty()) {
                if (merged_angles.size() > (allow_diagonal_edges ? 4 : 1)) {
                    merged_angles.resize(allow_diagonal_edges ? 4 : 1);
                }
                base_line_angles = std::move(merged_angles);
            }
        }
        for (double building_angle : building_line_angles_) {
            bool represented = false;
            for (double angle : base_line_angles) {
                if (angleDistanceToOrthogonalSystem(angle, building_angle) < 5.0 * M_PI / 180.0) {
                    represented = true;
                    break;
                }
            }
            if (!represented && base_line_angles.size() < (allow_diagonal_edges ? 4 : 1)) {
                base_line_angles.push_back(building_angle);
            }
        }
        std::vector<double> base_thetas;
        base_thetas.reserve(base_line_angles.size());
        for (double line_angle : base_line_angles) {
            base_thetas.push_back(line_angle - M_PI / 2.0);
        }
        if (base_thetas.empty()) {
            base_thetas.push_back(calculateGlobalDominantAngle(current_polygon));
        }
        if (iter == 0) {
            std::cerr << "[DirectionSystems] base_count=" << base_line_angles.size();
            for (double angle : base_line_angles) {
                std::cerr << " angle_deg=" << angle * 180.0 / M_PI;
            }
            std::cerr << std::endl;
        }

       
        //// 4. 边分类与强制吸附 (核心修改：消除盲区，强制 0, 45, 90, 135, 180 度)
        //auto normalize_ang = [](double a) {
        //    while (a <= -M_PI) a += 2 * M_PI;
        //    while (a > M_PI) a -= 2 * M_PI;
        //    return a;
        //};

        //double step = M_PI / 4.0; // 45度间隔的弧度值

        //for (size_t i = 0; i < n; ++i) {
        //    double diff = normalize_ang(edges[i].current_theta - base_theta);

        //    // 四舍五入到最近的 45 度倍数 (-4, -3, -2, -1, 0, 1, 2, 3, 4)
        //    double factor = std::round(diff / step);
        //    edges[i].theta_offset = factor * step;

        //    int abs_factor = std::abs(static_cast<int>(factor));

        //    // 根据绝对倍数分配类型
        //    if (abs_factor == 0 || abs_factor == 4) {
        //        edges[i].type = 0; // 平行 (0 或 180度)
        //    }
        //    else if (abs_factor == 2) {
        //        edges[i].type = 1; // 垂直 (90度)
        //    }
        //    else if (abs_factor == 1) {
        //        edges[i].type = 2; // 45度方向
        //    }
        //    else if (abs_factor == 3) {
        //        edges[i].type = 3; // 135度方向
        //    }
        //    // 此时所有边都被强制约束，不存在自由边 (type == -1) 了
        //}
        // 4. 边分类 (Grouping)
        // lambda：将角度归一化到 (-PI, PI]
        auto normalize_ang = [](double a) {
            while (a <= -M_PI) a += 2 * M_PI;
            while (a > M_PI) a -= 2 * M_PI;
            return a;
        };
        // Strict buildings assign every edge to the nearest axis. Buildings
        // with credible diagonal evidence keep edges outside the normal
        // orthogonal tolerance as independent free directions.
        double tolerance = allow_diagonal_edges
            ? std::min(tuning.angle_tolerance, 15.0 * M_PI / 180.0)
            : M_PI;
        std::vector<int> edge_base_index(n, -1);

        for (size_t i = 0; i < n; ++i) {
            double best_error = tolerance;
            int best_type = -1;
            int best_base = -1;
            double best_offset = 0.0;

            for (size_t k = 0; k < base_thetas.size(); ++k) {
                double diff = normalize_ang(edges[i].current_theta - base_thetas[k]);
                double parallel_error = std::min(std::abs(diff), std::abs(std::abs(diff) - M_PI));
                if (parallel_error < best_error) {
                    best_error = parallel_error;
                    best_type = 0;
                    best_base = static_cast<int>(k);
                    best_offset = (std::abs(diff) < M_PI / 2.0) ? 0.0 : (diff >= 0.0 ? M_PI : -M_PI);
                }

                double perpendicular_error = std::abs(std::abs(diff) - M_PI / 2.0);
                if (perpendicular_error < best_error) {
                    best_error = perpendicular_error;
                    best_type = 1;
                    best_base = static_cast<int>(k);
                    best_offset = (diff > 0.0) ? M_PI / 2.0 : -M_PI / 2.0;
                }
            }

            if (best_type != -1) {
                edges[i].type = best_type;
                edges[i].theta_offset = best_offset;
                edge_base_index[i] = best_base;
            }
            else {
                edges[i].type = -1; // 自由边
            }
        }


        // 5. 数据关联 (带拐角遮蔽优化)
        if (cloud && !cloud->empty()) {
            auto supportBaseWeight = [&](size_t pt_index) {
                if (pt_index < support_weights_.size()) {
                    return clampDouble(support_weights_[pt_index], 0.10, 1.0);
                }
                return 1.0;
            };
            for (size_t pt_index = 0; pt_index < cloud->points.size(); ++pt_index) {
                const auto& pt = cloud->points[pt_index];
                double min_dist = 1e9;
                int best_edge = -1;

                // 寻找最近的边
                for (size_t i = 0; i < n; ++i) {
                    double dist = computePointToSegmentDistance(pt, current_polygon[i], current_polygon[(i + 1) % n]);
                    if (dist < min_dist) { min_dist = dist; best_edge = i; }
                }

                // 阈值过滤 + 拐角遮蔽
                if (best_edge != -1 && min_dist < tuning.association_distance) {
                    pcl::PointXYZ p1 = current_polygon[best_edge];
                    pcl::PointXYZ p2 = current_polygon[(best_edge + 1) % n];
                    double lx = p2.x - p1.x, ly = p2.y - p1.y;
                    double len_sq = lx * lx + ly * ly;

                    if (len_sq > 1e-6) {
                        double t = ((pt.x - p1.x) * lx + (pt.y - p1.y) * ly) / len_sq;
                        // 【关键】只保留线段中间 90% 的点，避开拐角噪声
                        if (t > 0.05 && t < 0.95) {
                            edges[best_edge].associated_points.push_back(pt);
                            edges[best_edge].associated_weights.push_back(supportBaseWeight(pt_index));
                        }
                    }
                }
            }

            size_t associated_count = 0;
            size_t inlier_count = 0;
            for (size_t i = 0; i < n; ++i) {
                associated_count += edges[i].associated_points.size();
                robustEdgeInliersWithWeights(
                    edges[i].associated_points,
                    edges[i].associated_weights,
                    current_polygon[i], current_polygon[(i + 1) % n], tuning);
                inlier_count += edges[i].associated_points.size();
            }
            if (iter == 0) {
                std::cerr << "[Support] associated=" << associated_count
                    << " robust_inliers=" << inlier_count << std::endl;
            }
        }

        // --- Step 6: 构建 Ceres 问题 ---
        ceres::Problem problem;

        // 参数块 1: 每个方向组的基准角度
        std::vector<double> initial_base_thetas = base_thetas;
        for (size_t k = 0; k < base_thetas.size(); ++k) {
            problem.AddParameterBlock(&base_thetas[k], 1);
            problem.AddResidualBlock(ParameterRegularizer::Create(initial_base_thetas[k], 0.5), nullptr, &base_thetas[k]);
        }

        std::vector<double> initial_thetas(n);
        for (size_t i = 0; i < n; ++i) {
            initial_thetas[i] = edges[i].current_theta;
        }

        double support_weight_sum = 0.0;
        double support_weight_min = std::numeric_limits<double>::max();
        double support_weight_max = 0.0;
        double normal_weight_sum = 0.0;
        double normal_weight_min = std::numeric_limits<double>::max();
        double normal_weight_max = 0.0;
        size_t support_weight_count = 0;
        const double weight_sigma = std::max(0.35, tuning.association_distance * 0.45);

        for (int i = 0; i < n; ++i) {
            // 参数块 2: 每条边的平移距离 d
            problem.AddParameterBlock(&edges[i].d, 1);

            // 兜底逻辑：如果有自由边则添加独立 theta 参数（在此策略下不会触发，但保留以防扩充）
            if (edges[i].type == -1) {
                problem.AddParameterBlock(&edges[i].current_theta, 1);
                problem.AddResidualBlock(ParameterRegularizer::Create(initial_thetas[i], 0.5), nullptr, &edges[i].current_theta);
            }

            // 添加点到直线的拟合残差
            for (size_t pt_index = 0; pt_index < edges[i].associated_points.size(); ++pt_index) {
                const auto& pt = edges[i].associated_points[pt_index];
                const double boundary_distance = computePointToSegmentDistance(
                    pt, current_polygon[static_cast<size_t>(i)],
                    current_polygon[(static_cast<size_t>(i) + 1) % n]);
                const double distance_weight = std::exp(
                    -(boundary_distance * boundary_distance) /
                    (2.0 * weight_sigma * weight_sigma));
                const double boundary_weight = clampDouble(0.20 + 0.80 * distance_weight, 0.20, 1.0);
                const double normal_weight = pt_index < edges[i].associated_weights.size()
                    ? clampDouble(edges[i].associated_weights[pt_index], 0.10, 1.0)
                    : 1.0;
                const double support_weight = clampDouble(boundary_weight * normal_weight, 0.10, 1.0);
                const double sqrt_weight = std::sqrt(support_weight);
                support_weight_sum += support_weight;
                support_weight_min = std::min(support_weight_min, support_weight);
                support_weight_max = std::max(support_weight_max, support_weight);
                normal_weight_sum += normal_weight;
                normal_weight_min = std::min(normal_weight_min, normal_weight);
                normal_weight_max = std::max(normal_weight_max, normal_weight);
                ++support_weight_count;

                if (edges[i].type != -1) {
                    // 受约束边：共享 base_theta
                    ceres::CostFunction* cost = PointToLineDistanceCost::Create(
                        pt.x, pt.y, edges[i].theta_offset, sqrt_weight);
                    int base_idx = std::max(0, edge_base_index[i]);
                    problem.AddResidualBlock(cost, new ceres::HuberLoss(tuning.huber_delta), &base_thetas[base_idx], &edges[i].d);
                }
                else {
                    // 自由边
                    ceres::CostFunction* cost = FreeLineCost::Create(pt.x, pt.y, sqrt_weight);
                    problem.AddResidualBlock(cost, new ceres::HuberLoss(tuning.huber_delta), &edges[i].current_theta, &edges[i].d);
                }
            }

            if (use_dlg_position_) {
                const double theta = edges[i].current_theta;
                const double current_line_angle = theta + M_PI / 2.0;
                const double max_angle_error = 15.0 * M_PI / 180.0;
                const double max_offset = clampDouble(0.04 * std::sqrt(tuning.area), 0.5, 2.0);
                double best_offset = max_offset;
                double target_d = edges[i].d;
                bool matched = false;

                for (size_t j = 0; j < dlg_polygon_.size(); ++j) {
                    const auto& a = dlg_polygon_[j];
                    const auto& b = dlg_polygon_[(j + 1) % dlg_polygon_.size()];
                    const double dlg_angle = std::atan2(b.y - a.y, b.x - a.x);
                    double angle_error = std::abs(normalizeAnglePi(current_line_angle - dlg_angle));
                    angle_error = std::min(angle_error, std::abs(M_PI - angle_error));
                    if (angle_error > max_angle_error) continue;

                    const double candidate_d = 0.5 * (a.x + b.x) * std::cos(theta) +
                        0.5 * (a.y + b.y) * std::sin(theta);
                    const double offset = std::abs(candidate_d - edges[i].d);
                    if (offset < best_offset) {
                        best_offset = offset;
                        target_d = candidate_d;
                        matched = true;
                    }
                }

                if (matched) {
                    const double position_weight = 0.35 + 0.45 * dlg_confidence_;
                    problem.AddResidualBlock(
                        ParameterRegularizer::Create(target_d, position_weight),
                        new ceres::HuberLoss(tuning.huber_delta), &edges[i].d);
                }
            }
        }
        if (iter == 0 && support_weight_count > 0) {
            std::cerr << "[SupportWeight] count=" << support_weight_count
                << " avg=" << support_weight_sum / support_weight_count
                << " min=" << support_weight_min
                << " max=" << support_weight_max
                << " normal_avg=" << normal_weight_sum / support_weight_count
                << " normal_min=" << normal_weight_min
                << " normal_max=" << normal_weight_max
                << " sigma=" << weight_sigma
                << std::endl;
        }

        if (inter_floor_wall_snap_weight_ > 0.0 && !reference_walls_.empty()) {
            for (size_t i = 0; i < n; ++i) {
                int base_idx = std::max(0, edge_base_index[i]);
                double theta = (edges[i].type == -1)
                    ? edges[i].current_theta
                    : (base_thetas[base_idx] + edges[i].theta_offset);

                double best_d = 0.0;
                double best_distance = inter_floor_max_wall_snap_distance_;
                bool matched = false;

                for (const auto& ref_wall : reference_walls_) {
                    double ref_theta = ref_wall.normal_theta;
                    double ref_d = ref_wall.d;
                    double angle_error = std::abs(normalize_ang(theta - ref_theta));

                    if (angle_error > M_PI / 2.0) {
                        ref_theta = normalize_ang(ref_theta + M_PI);
                        ref_d = -ref_d;
                        angle_error = std::abs(normalize_ang(theta - ref_theta));
                    }

                    if (angle_error > inter_floor_max_wall_snap_angle_) {
                        continue;
                    }

                    double d_distance = std::abs(edges[i].d - ref_d);
                    if (d_distance < best_distance) {
                        best_distance = d_distance;
                        best_d = ref_d;
                        matched = true;
                    }
                }

                if (matched) {
                    double edge_len = std::hypot(
                        current_polygon[(i + 1) % n].x - current_polygon[i].x,
                        current_polygon[(i + 1) % n].y - current_polygon[i].y);
                    double length_weight = std::max(0.5, std::min(2.0, edge_len / 10.0));
                    problem.AddResidualBlock(
                        ParameterRegularizer::Create(best_d, inter_floor_wall_snap_weight_ * length_weight),
                        new ceres::HuberLoss(tuning.huber_delta),
                        &edges[i].d);
                }
            }
        }

        // --- Step 7: 【优化】添加中点锚定约束 (取代原来的顶点锚定) ---
        // 防止规则化后的直线偏离原始位置太远，但允许两端自由延伸以修正夹角
        double mid_anchor_weight = tuning.mid_anchor_weight; // 适中权重即可
        ceres::LossFunction* mid_loss = new ceres::ScaledLoss(
            new ceres::HuberLoss(tuning.huber_delta), mid_anchor_weight, ceres::TAKE_OWNERSHIP);

        for (size_t i = 0; i < n; ++i) {
            pcl::PointXYZ p_mid = mid_points[i];

            if (edges[i].type != -1) {
                ceres::CostFunction* cost = PointToLineDistanceCost::Create(p_mid.x, p_mid.y, edges[i].theta_offset);
                int base_idx = std::max(0, edge_base_index[i]);
                problem.AddResidualBlock(cost, mid_loss, &base_thetas[base_idx], &edges[i].d);
            }
            else {
                ceres::CostFunction* cost = FreeLineCost::Create(p_mid.x, p_mid.y);
                problem.AddResidualBlock(cost, mid_loss, &edges[i].current_theta, &edges[i].d);
            }
        }

        // 8. 求解
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.max_num_iterations = 100;
        options.minimizer_progress_to_stdout = false; // 放在循环里建议关闭输出，避免刷屏
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        // 9. 结果写回 optimized_params，准备重建
        optimized_params.resize(n);
        for (int i = 0; i < n; ++i) {
            int base_idx = std::max(0, edge_base_index[i]);
            double theta = (edges[i].type == -1) ? edges[i].current_theta : (base_thetas[base_idx] + edges[i].theta_offset);
            // 转换回直线一般式 ax + by + c = 0  => a=cos, b=sin, c=-d
            optimized_params[i] = Eigen::Vector3d(std::cos(theta), std::sin(theta), -edges[i].d);
        }

        // 10. 【核心迭代补充】根据优化后的直线参数，重新计算交点作为下一轮的多边形
        std::vector<pcl::PointXYZ> new_polygon;
        for (size_t i = 0; i < n; ++i) {
            size_t j = (i + 1) % n;
            pcl::PointXYZ intersect_pt = computeIntersection(optimized_params[i], optimized_params[j], current_polygon[0].z);

            // 防御性编程：如果计算交点失败（比如由于误差出现了绝对平行线），防止存入非法点
            if (!std::isnan(intersect_pt.x) && !std::isnan(intersect_pt.y)) {
                new_polygon.push_back(intersect_pt);
            }
        }

        // 传递前做一次短边剔除，防止新产生的退化短边进入下一轮迭代
        if (new_polygon.size() >= 3) {
            std::vector<pcl::PointXYZ> pruned;
            pruneShortEdges(new_polygon, pruned, tuning.fine_prune_distance);
            if (pruned.size() >= 3 && isSimplePolygon2D(pruned)) {
                new_polygon = pruned;
            }
        }

        // 传递给下一次迭代
        const bool new_polygon_valid =
            new_polygon.size() >= 3 &&
            isSimplePolygon2D(new_polygon) &&
            polygonIoU2D(new_polygon, backup_polygon) >= (allow_diagonal_edges ? 0.55 : 0.65);

        if (!new_polygon_valid) {
            std::cerr << "[Ceres] reject invalid iteration polygon; keep last valid polygon"
                << std::endl;
            current_polygon = last_valid_polygon;
            break;
        }

        if (new_polygon.size() == n && new_polygon.size() >= 3) {
            current_polygon = new_polygon;
            last_valid_polygon = current_polygon;
        }
        else if (new_polygon.size() >= 3) {
            // 顶点数变化但多边形仍有效，继续迭代
            current_polygon = new_polygon;
            last_valid_polygon = current_polygon;
        }
        else {
            // 如果求交失败导致拓扑破坏，终止迭代
            break;
        }
    }

    // 所有的迭代结束后，更新最终假设
    if (current_polygon.size() < 3 || !isSimplePolygon2D(current_polygon)) {
        std::cerr << "[Ceres] final optimized polygon invalid; rollback to last valid polygon"
            << std::endl;
        current_polygon = last_valid_polygon;
    }
    best_hypothesis = current_polygon;
}

// ===== computeIntersection =====
// 作用：辅助计算两条一般式直线(ax+by+c=0)的交点
// 参数：l1,l2 - 两条直线的 [a,b,c] 参数; z_value - 输出点的 z 值
// 返回：交点（近平行时返回 NaN）
pcl::PointXYZ outlineRegular::computeIntersection(const Eigen::Vector3d& l1, const Eigen::Vector3d& l2, float z_value) {
    double det = l1[0] * l2[1] - l2[0] * l1[1];
    pcl::PointXYZ pt;
    if (std::abs(det) < 1e-5) {
        // 平行线兜底
        pt.x = std::numeric_limits<float>::quiet_NaN();
        pt.y = std::numeric_limits<float>::quiet_NaN();
    }
    else {
        // 克莱姆法则
        pt.x = static_cast<float>((l1[1] * l2[2] - l2[1] * l1[2]) / det);
        pt.y = static_cast<float>((l2[0] * l1[2] - l1[0] * l2[2]) / det);
    }
    pt.z = z_value;
    return pt;
}

// ===== calculateGlobalDominantAngle =====
// 作用：基于边长加权的角度直方图计算全局主方向（法向量角度）
// 参数：hypothesis - 多边形顶点
// 返回：主方向对应的法向量角度（弧度），即线段角度 - PI/2
double outlineRegular::calculateGlobalDominantAngle(const std::vector<pcl::PointXYZ>& hypothesis) {
    const int bin_size = 180; // 0.5度精度，映射到 90度区间 -> 180个bin
    std::vector<double> histogram(bin_size, 0.0);
    size_t n = hypothesis.size();

    for (size_t i = 0; i < n; ++i) {
        pcl::PointXYZ p1 = hypothesis[i];
        pcl::PointXYZ p2 = hypothesis[(i + 1) % n];
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        double len = std::hypot(dx, dy);

        if (len < 1e-3) continue;

        // 计算角度并映射到 [0, 90)
        double angle_deg = std::atan2(dy, dx) * 180.0 / M_PI;
        // 归一化到 [0, 90)
        while (angle_deg < 0) angle_deg += 90.0;
        while (angle_deg >= 90.0) angle_deg -= 90.0;

        // 映射到 bin 索引 (0-179)
        int bin = static_cast<int>(angle_deg * 2.0);
        if (bin >= 0 && bin < bin_size) {
            histogram[bin] += len; // 长度加权
        }
    }

    // 简单平滑 (3点滑动平均)
    std::vector<double> smoothed = histogram;
    for (int i = 1; i < bin_size - 1; ++i) {
        smoothed[i] = 0.25 * histogram[i - 1] + 0.5 * histogram[i] + 0.25 * histogram[i + 1];
    }
    // 处理边界
    smoothed[0] = 0.5 * histogram[0] + 0.25 * histogram[1] + 0.25 * histogram[bin_size - 1];
    smoothed[bin_size - 1] = 0.5 * histogram[bin_size - 1] + 0.25 * histogram[bin_size - 2] + 0.25 * histogram[0];

    // 找峰值
    int best_bin = 0;
    double max_weight = -1.0;
    for (int i = 0; i < bin_size; ++i) {
        if (smoothed[i] > max_weight) {
            max_weight = smoothed[i];
            best_bin = i;
        }
    }

    // 转回法向量角度 (line_angle - 90度)
    // 注意：直方图给的是线段方向 [0,90)，我们需要返回对应的法向量角度
    // 法向量 = 线段角度 - 90度
    double dom_angle = (static_cast<double>(best_bin) / 2.0) * M_PI / 180.0;
    return dom_angle - M_PI / 2.0;
}

// ===== pruneShortEdges =====
// 作用：短边剔除（增强版）——移除多边形中短于阈值的边对应的冗余顶点，
//       加入共线判断智能合并（优先删除共线中间点）
// 参数：input - 输入顶点; output - 输出顶点; min_edge_length - 最小边长阈值
void outlineRegular::pruneShortEdges(const std::vector<pcl::PointXYZ>& input,
    std::vector<pcl::PointXYZ>& output,
    double min_edge_length)
{
    if (input.empty()) {
        output.clear();
        return;
    }

    output = input;
    bool total_changed = true;
    double collinear_threshold_deg = 15.0; // 夹角 > 165° 视为共线

    while (total_changed) {
        total_changed = false;

        // 1. 检查索引相邻的点对 (i 和 i+1)
        for (size_t i = 0; i + 1 < output.size(); ) {
            double dx = output[i].x - output[i + 1].x;
            double dy = output[i].y - output[i + 1].y;
            double dz = output[i].z - output[i + 1].z;
            double dist = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (dist < min_edge_length) {
                // 【增强】检查是否三个连续点共线，如果是则智能删除中间点
                if (output.size() >= 3) {
                    size_t prev_idx = (i > 0) ? (i - 1) : (output.size() - 1);
                    size_t next_idx = (i + 2 < output.size()) ? (i + 2) : 0;

                    pcl::PointXYZ a = output[prev_idx];
                    pcl::PointXYZ b = output[i];
                    pcl::PointXYZ c = output[i + 1];
                    pcl::PointXYZ d = output[next_idx];

                    // 情况A：a-b-c 共线，且 ab 或 bc 中有短边，删除 b
                    double angle_abc = computeAngle2D(a, b, c);
                    if (std::abs(angle_abc - 180.0) < collinear_threshold_deg) {
                        output.erase(output.begin() + i);
                        total_changed = true;
                        continue; // 不i++，重新检查新i位置
                    }

                    // 情况B：b-c-d 共线，且 bc 或 cd 中有短边，删除 c
                    double angle_bcd = computeAngle2D(b, c, d);
                    if (std::abs(angle_bcd - 180.0) < collinear_threshold_deg) {
                        output.erase(output.begin() + i + 1);
                        total_changed = true;
                        continue;
                    }
                }

                // 默认行为：删除后一个点 (i + 1)
                output.erase(output.begin() + i + 1);
                total_changed = true;
            }
            else {
                i++;
            }
        }

        // 2. 检查首尾闭合
        if (output.size() >= 3) {
            double dx = output.back().x - output.front().x;
            double dy = output.back().y - output.front().y;
            double dz = output.back().z - output.front().z;
            double dist_loop = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (dist_loop < min_edge_length) {
                // 检查首尾是否共线：前两点和倒数第二点
                if (output.size() >= 4) {
                    double angle_last_first_second = computeAngle2D(
                        output[output.size() - 2], output.back(), output[0]);
                    if (std::abs(angle_last_first_second - 180.0) < collinear_threshold_deg) {
                        // 末尾点是共线冗余，删末尾
                        output.pop_back();
                        total_changed = true;
                        continue;
                    }

                    double angle_last_first_front = computeAngle2D(
                        output.back(), output[0], output[1]);
                    if (std::abs(angle_last_first_front - 180.0) < collinear_threshold_deg) {
                        // 首点是共线冗余，删首
                        output.erase(output.begin());
                        total_changed = true;
                        continue;
                    }
                }
                output.erase(output.begin());
                total_changed = true;
            }
        }
        else {
            break;
        }
    }
}

// ===== pointsEqual =====
// 作用：判断两个三维点是否在给定容差内相等
// 参数：a,b - 待比较点; epsilon - 坐标容差
// 返回：true 表示三点坐标差均小于 epsilon
bool pointsEqual(const pcl::PointXYZ& a, const pcl::PointXYZ& b, float epsilon = 0.1f) {
    return std::fabs(a.x - b.x) < epsilon &&
        std::fabs(a.y - b.y) < epsilon &&
        std::fabs(a.z - b.z) < epsilon;
}

// ===== areCollinearRobust =====
// 作用：健壮的共线判断——同时处理 0度(同向直行) 和 180度(反向折返Spike) 两种情况
// 参数：a,b,c - 三个连续点; angle_thresh_deg - 角度容差（度）
// 返回：true 表示三点共线（同向或反向）
bool areCollinearRobust(const pcl::PointXYZ& a, const pcl::PointXYZ& b, const pcl::PointXYZ& c, double angle_thresh_deg = 10.0) {
    Eigen::Vector2d ab(b.x - a.x, b.y - a.y);
    Eigen::Vector2d bc(c.x - b.x, c.y - b.y);

    double norm_ab = ab.norm();
    double norm_bc = bc.norm();

    // 极短边直接视为共线交点（兜底）
    if (norm_ab < 1e-4 || norm_bc < 1e-4) return true;

    ab /= norm_ab;
    bc /= norm_bc;

    // 计算点积，对应 cos(theta)
    double dot = ab.dot(bc);
    dot = std::max(-1.0, std::min(1.0, dot)); // 截断，防止精度越界导致 acos 崩溃

    // 理想直线前进夹角为 0度，完全折返(Spike)夹角为 180度
    double angle_deg = std::acos(dot) * 180.0 / M_PI;

    // 核心修改：同时剔除同向共线点和反向重叠毛刺
    return (std::abs(angle_deg) < angle_thresh_deg) || (std::abs(angle_deg - 180.0) < angle_thresh_deg);
}

// ===== OptimizeFinal_points =====
// 作用：闭环优化——对最终结果点云做去重、移除共线点、批量拓扑修复、再清理共线点。
//       是规则化结果的最后一道后处理关卡。
// 参数：input - 最终结果点云（就地修改）
void outlineRegular::OptimizeFinal_points(pcl::PointCloud<pcl::PointXYZ>::Ptr& input)
{
    if (input->points.size() < 3) return;

    std::vector<pcl::PointXYZ> pts;
    for (size_t i = 0; i < input->points.size(); ++i) {
        pts.push_back(input->points[i]);
    }

    double area = computeOBBArea(pts);
    if (area < 1e-6) area = polygonArea2D(pts);
    double resolution = 0.0;
    if (ransac_inner_cloud && !ransac_inner_cloud->empty()) {
        XG::modularFunction mfn;
        resolution = mfn.computeModelResolution(*ransac_inner_cloud);
    }
    OutlineTuning tuning = makeOutlineTuning(resolution, area);
    const size_t before_count = pts.size();
    const double cluster_threshold = clampDouble(tuning.repair_distance * 2.0, tuning.repair_distance, 7.0);
    std::cerr << "[Topo] before=" << before_count
              << " resolution=" << tuning.resolution
              << " area=" << area
              << " repair_distance=" << tuning.repair_distance
              << " cluster_threshold=" << cluster_threshold << std::endl;

    // 1. 去重
    removeDuplicatePoints2D(pts, static_cast<float>(tuning.fine_prune_distance));

    // 2. 闭环优化，移除共线点
    removeCollinearPoints(pts, 20.0f);

    // 3. 批量拓扑修复。先修单短边模式，再修连续短边簇，最后再跑一轮轻量兜底。
    for (int pass = 0; pass < 2 && pts.size() >= 3; ++pass) {
        const size_t pass_before = pts.size();
        repairTopologyBatch(pts, tuning.repair_distance);
        resolveShortEdgeCluster(pts, cluster_threshold);
        resolveAngleAnomalies(pts, cluster_threshold);
        removeDuplicatePoints2D(pts, static_cast<float>(tuning.fine_prune_distance));
        removeCollinearPoints(pts, 20.0f);
        if (pts.size() == pass_before) break;
    }

    std::cerr << "[Topo] after=" << pts.size()
              << " removed=" << (before_count >= pts.size() ? before_count - pts.size() : 0)
              << std::endl;

    if (pts.size() >= 3) {
        input->points.clear();
        for (const auto& p : pts) {
            input->points.push_back(p);
        }
    }

    //if (!input || input->points.empty()) return;

    //// 1. 将输入的有序边界点云提取到 vector 容器中
    //std::vector<pcl::PointXYZ> boundary_pts(input->points.begin(), input->points.end());

    //// 2. 调用图片中对应的官方 CGAL 规整算法（传入 0.5 米的距离偏差容差）
    //double tolerance = 0.5;
    //cgalContourRegularize(boundary_pts, tolerance);

    //// 3. 将规整后的新几何拓扑顶点刷新回 input 点云指针中
    //input->clear();
    //for (const auto& pt : boundary_pts) {
    //    input->points.push_back(pt);
    //}
}

// ===== resolveShortEdgeIntersections =====
// 作用：建筑轮廓”倒角(Chamfer)”消除——仅在法线不相交（双钝角）时触发，
//       延长相邻两条长边求交，替换短边连接
// 参数：pts - 多边形顶点（就地修改）; short_edge_threshold - 短边判定阈值
void outlineRegular::resolveShortEdgeIntersections(std::vector<pcl::PointXYZ>& pts, double short_edge_threshold)
{
    if (pts.size() < 4) return;

    bool changed = true;
    int max_iters = 20;
    int iter = 0;

    while (changed && iter < max_iters) {
        changed = false;
        int n = pts.size();
        if (n < 4) break;

        for (int i = 0; i < n; ++i) {
            int i_prev = (i - 1 + n) % n; // a
            int i_curr = i;               // b
            int i_next = (i + 1) % n;     // c
            int i_nnext = (i + 2) % n;    // d

            pcl::PointXYZ a = pts[i_prev];
            pcl::PointXYZ b = pts[i_curr];
            pcl::PointXYZ c = pts[i_next];
            pcl::PointXYZ d = pts[i_nnext];

            double bc_len = std::hypot(c.x - b.x, c.y - b.y);
            double ab_len = std::hypot(b.x - a.x, b.y - a.y);
            double cd_len = std::hypot(d.x - c.x, d.y - c.y);

            // 基础约束：当前边是短边
            if (bc_len > 0.01 && bc_len < short_edge_threshold&&ab_len>bc_len&&cd_len>bc_len) {

                // 构造用于判断拐角性质的向量
                Eigen::Vector2d ba(a.x - b.x, a.y - b.y);
                Eigen::Vector2d bc(c.x - b.x, c.y - b.y);
                Eigen::Vector2d cb(b.x - c.x, b.y - c.y);
                Eigen::Vector2d cd(d.x - c.x, d.y - c.y);

                // 【核心改进：你的法线不相交判据】
                // 只有当两个拐角都是钝角时（点积均小于0），才证明这是个“倒角切边”
                if (ba.dot(bc) < -1e-3 && cd.dot(cb) < -1e-3) {

                    double len_ab = ba.norm();
                    double len_cd = cd.norm();

                    // 防过度矫正：相邻墙体必须显著长于倒角本身
                    if (len_ab > bc_len * 1.0 && len_cd > bc_len * 1.0) {

                        // 求直线 ab 和 cd 的交点
                        double A1 = b.y - a.y; double B1 = a.x - b.x;
                        double C1 = A1 * a.x + B1 * a.y;

                        double A2 = d.y - c.y; double B2 = c.x - d.x;
                        double C2 = A2 * c.x + B2 * c.y;

                        double det = A1 * B2 - A2 * B1;

                        if (std::abs(det) > 0.01) {
                            pcl::PointXYZ intersect_pt;
                            intersect_pt.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
                            intersect_pt.y = static_cast<float>((A1 * C2 - A2 * C1) / det);
                            intersect_pt.z = b.z;

                            // 安全检查：交点不能飞得太远
                            if (std::hypot(intersect_pt.x - b.x, intersect_pt.y - b.y) < bc_len * 1.5 &&
                                std::hypot(intersect_pt.x - c.x, intersect_pt.y - c.y) < bc_len * 1.5) {

                                pts[i_curr] = intersect_pt;
                                pts.erase(pts.begin() + i_next);
                                changed = true;
                                break; // 拓扑改变，重新迭代
                            }
                        }
                    }
                }
            }
        }
        iter++;
    }
}

// ===== resolveOrthogonalSpikes =====
// 作用：消除正交边界之间的”直角尖刺/阶梯”——仅在法线相交（存在直角/锐角）时触发，
//       要求 ab 和 cd 近似垂直时延长求交消除尖刺
// 参数：pts - 多边形顶点（就地修改）; max_spike_length - 尖刺最大长度
void outlineRegular::resolveOrthogonalSpikes(std::vector<pcl::PointXYZ>& pts, double max_spike_length)
{
    if (pts.size() < 4) return;

    bool changed = true;
    int max_iters = 20;
    int iter = 0;

    while (changed && iter < max_iters) {
        changed = false;
        int n = pts.size();
        if (n < 4) break;

        for (int i = 0; i < n; ++i) {
            int i_a = (i - 1 + n) % n;
            int i_b = i;
            int i_c = (i + 1) % n;
            int i_d = (i + 2) % n;

            pcl::PointXYZ a = pts[i_a];
            pcl::PointXYZ b = pts[i_b];
            pcl::PointXYZ c = pts[i_c];
            pcl::PointXYZ d = pts[i_d];

            double bc_len = std::hypot(c.x - b.x, c.y - b.y);
            double ab_len = std::hypot(b.x - a.x, b.y - a.y);
            double cd_len = std::hypot(d.x - c.x, d.y - c.y);

            if (bc_len > 0.01 && bc_len < max_spike_length && ab_len>bc_len && cd_len > bc_len) {

                Eigen::Vector2d ba(a.x - b.x, a.y - b.y);
                Eigen::Vector2d bc(c.x - b.x, c.y - b.y);
                Eigen::Vector2d cb(b.x - c.x, b.y - c.y);
                Eigen::Vector2d cd(d.x - c.x, d.y - c.y);

                // 【核心改进：你的法线相交判据】
                // 只要拐角中存在直角或锐角（点积 >= 0），说明 bc 的法线必然与 ab 或 cd 相交，符合锯齿/阶梯特征
                if (ba.dot(bc) >= -1e-3 || cd.dot(cb) >= -1e-3) {

                    // 并且要求 ab 和 cd 近似垂直（正交错位特征）
                    Eigen::Vector2d ab_dir = -ba / ba.norm();
                    Eigen::Vector2d cd_dir = cd / cd.norm();
                    double dot_orthogonal = std::abs(ab_dir.dot(cd_dir));

                    if (dot_orthogonal < 0.5) { // 约 60° ~ 120° 之间

                        double A1 = b.y - a.y; double B1 = a.x - b.x;
                        double C1 = A1 * a.x + B1 * a.y;

                        double A2 = d.y - c.y; double B2 = c.x - d.x;
                        double C2 = A2 * c.x + B2 * c.y;

                        double det = A1 * B2 - A2 * B1;

                        if (std::abs(det) > 0.01) {
                            pcl::PointXYZ p;
                            p.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
                            p.y = static_cast<float>(((A1 * C2 - A2 * C1) / det));
                            p.z = b.z;

                            if (std::hypot(p.x - b.x, p.y - b.y) < max_spike_length * 2.5 &&
                                std::hypot(p.x - c.x, p.y - c.y) < max_spike_length * 2.5) {

                                pts[i_b] = p;
                                pts.erase(pts.begin() + i_c);
                                changed = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
        iter++;
    }
}

// ===== cgalContourRegularize =====
// 作用：调用 CGAL Shape_regularization 官方算法对闭合轮廓做正交规则化
// 参数：pts - 多边形顶点（就地修改）; distance_tolerance - 距离容差（最大偏移）
void outlineRegular::cgalContourRegularize(std::vector<pcl::PointXYZ>& pts, double distance_tolerance) {
    if (pts.size() < 3) {
        std::cerr << "[Warning] 点数少于3，无法构成多边形闭合轮廓！" << std::endl;
        return;
    }

    // 步骤 1: 数据类型转换（将 PCL PointXYZ 转换为 CGAL 二维点平面轮廓）
    Contour contour;
    contour.reserve(pts.size());
    for (const auto& pt : pts) {
        contour.push_back(Point_2(pt.x, pt.y));
    }

    // 步骤 2: 估计建筑物主方向（CGAL 内部的高级几何逻辑）
    // 通过最长边或直方图自适应拟合出多边形的主墙面方向，从而支持多角度正交
    const FT min_length = FT(0.3);   // 边长阈值：小于0.3米的短边不参与基准方向计算，避免锯齿噪声干扰
    const FT max_angle = FT(25.0);  // 角度阈值（角度制）：在25度以内的斜边均会被强行向主方向靠拢
    const bool is_closed = true;     // 建筑物属于闭合多边形

    Contour_directions directions(
        contour,
        is_closed,
        CGAL::parameters::minimum_length(min_length).maximum_angle(max_angle)
    );

    // 步骤 3: 调用 CGAL 官方核心函数进行闭合轮廓规整
    Contour regularized_contour;
    CGAL::Shape_regularization::Contours::regularize_closed_contour(
        contour,
        directions,
        std::back_inserter(regularized_contour),
        CGAL::parameters::maximum_offset(FT(distance_tolerance)) // 对应图片中的距离容差参数
    );

    // 步骤 4: 结果写回 PCL 格式
    if (!regularized_contour.empty()) {
        // 计算原本点云的平均 Z 轴高度，保证规整后房顶高度或地面高度不丢失
        float avg_z = 0.0f;
        for (const auto& pt : pts) avg_z += pt.z;
        avg_z /= pts.size();

        pts.clear();
        pts.reserve(regularized_contour.size());
        for (const auto& cgal_pt : regularized_contour) {
            pcl::PointXYZ pcl_pt;
            pcl_pt.x = static_cast<float>(cgal_pt.x());
            pcl_pt.y = static_cast<float>(cgal_pt.y());
            pcl_pt.z = avg_z;
            pts.push_back(pcl_pt);
        }
    }
    else {
        std::cerr << "[Error] CGAL 轮廓规整失败，返回空结果！" << std::endl;
    }
}

// ===== regularizeBuildingOutline =====
// 作用：调用 CGAL 正交化算法对建筑轮廓做规则化（封装接口，当前核心调用已注释）
// 参数：pts - 多边形顶点（就地修改）; tolerance - 容差
void outlineRegular::regularizeBuildingOutline(std::vector<pcl::PointXYZ>& pts, double tolerance) {
    if (pts.size() < 3) {
        std::cerr << "[Warning] 点数少于3，无法构成闭合多边形轮廓！" << std::endl;
        return;
    }
    double z = pts[0].z;

    // 1. 将 PCL 点云转换为 CGAL 多边形
    Polygon_2 poly;
    for (const auto& pt : pts) {
        poly.push_back(Polygon_2::Point_2(pt.x, pt.y));
    }

    // 2. 调用 CGAL 官方正交化算法
    Polygon_2 orthogonal_poly;
    /*CGAL::polygon_vertical_orthogonalization_2(
        poly,
        std::back_inserter(orthogonal_poly),
        tolerance
    );*/

    // 3. 将规整后的结果写回 PCL 点云
    pts.clear();
    for (auto it = orthogonal_poly.vertices_begin(); it != orthogonal_poly.vertices_end(); ++it) {
        pcl::PointXYZ pt;
        pt.x = CGAL::to_double(it->x());
        pt.y = CGAL::to_double(it->y());
        pt.z = z; // 保持原有高度（若有需要可从原点云获取）
        pts.push_back(pt);
    }
}

// ===== resolveParallelStep =====
// 作用：消除平行边之间的斜向错位——将其强行转化为标准的正交阶梯。
//       当 ab 和 cd 近似平行同向且 bc 为短边时，沿墙面方向投影 b、c
// 参数：pts - 多边形顶点（就地修改）; step_threshold - 阶梯阈值
void outlineRegular::resolveParallelStep(std::vector<pcl::PointXYZ>& pts, double step_threshold)
{
    if (pts.size() < 4) return;
    bool changed = true;
    int max_iters = 20;
    int iter = 0;

    while (changed && iter < max_iters) {
        changed = false;
        int n = pts.size();
        if (n < 4) break;

        for (int i = 0; i < n; ++i) {
            int i_a = (i - 1 + n) % n;
            int i_b = i;
            int i_c = (i + 1) % n;
            int i_d = (i + 2) % n;

            pcl::PointXYZ a = pts[i_a]; pcl::PointXYZ b = pts[i_b];
            pcl::PointXYZ c = pts[i_c]; pcl::PointXYZ d = pts[i_d];

            double bc_len = std::hypot(c.x - b.x, c.y - b.y);

            // 条件1：bc 是一条短边
            if (bc_len > 0.01 && bc_len < step_threshold) {
                Eigen::Vector2d v_ab(b.x - a.x, b.y - a.y);
                Eigen::Vector2d v_cd(d.x - c.x, d.y - c.y);
                double len_ab = v_ab.norm();
                double len_cd = v_cd.norm();

                // 条件2：相邻的 ab 和 cd 必须相对较长，确保它们是主墙面
                if (len_ab > bc_len * 1.0 && len_cd > bc_len * 1.0) {
                    Eigen::Vector2d dir_ab = v_ab / len_ab;
                    Eigen::Vector2d dir_cd = v_cd / len_cd;

                    // 条件3：【核心判断】检查 ab 和 cd 是否近似平行且同向
                    // 修改：将 0.95 (18度) 放宽至 0.866 (30度)，与 Ceres 30度阈值同步
                    if (dir_ab.dot(dir_cd) > 0.866) {

                        // 提取两面墙的平均方向
                        Eigen::Vector2d dir = (dir_ab + dir_cd).normalized();

                        // 找到斜边 bc 的中点 M
                        pcl::PointXYZ M;
                        M.x = (b.x + c.x) / 2.0f;
                        M.y = (b.y + c.y) / 2.0f;
                        M.z = b.z;

                        // 核心几何变换：将 b 和 c 沿着原墙面方向投影，使其截断截面垂直于 dir
                        // 这样新的 b'c' 会完美垂直于主墙面
                        double t_b = (Eigen::Vector2d(M.x - b.x, M.y - b.y).dot(dir)) / (dir_ab.dot(dir));
                        pcl::PointXYZ b_new;
                        b_new.x = static_cast<float>(b.x + t_b * dir_ab.x());
                        b_new.y = static_cast<float>(b.y + t_b * dir_ab.y());
                        b_new.z = b.z;

                        double t_c = (Eigen::Vector2d(M.x - c.x, M.y - c.y).dot(dir)) / (dir_cd.dot(dir));
                        pcl::PointXYZ c_new;
                        c_new.x = static_cast<float>(c.x + t_c * dir_cd.x());
                        c_new.y = static_cast<float>(c.y + t_c * dir_cd.y());
                        c_new.z = c.z;

                        // 防护验证：避免因微小平行误差导致投影点飞走
                        if (std::hypot(b_new.x - b.x, b_new.y - b.y) < bc_len * 2.0 &&
                            std::hypot(c_new.x - c.x, c_new.y - c.y) < bc_len * 2.0) {

                            // 替换原有的 b 和 c，生成完美的正交阶梯
                            pts[i_b] = b_new;
                            pts[i_c] = c_new;
                            changed = true;
                            break; // 拓扑改变，跳出当前 for 循环进入下一轮检查
                        }
                    }
                }
            }
        }
        iter++;
    }
}

// ===== resolveDoubleShortEdgeSpike =====
// 作用：处理由两条连续短边造成的直角拐角破缺——将 b,c,d 替换为 ab 和 de 的交点
// 参数：pts - 多边形顶点（就地修改）; threshold - 短边判定阈值
void outlineRegular::resolveDoubleShortEdgeSpike(std::vector<pcl::PointXYZ>& pts, double threshold)
{
    if (pts.size() < 5) return;

    bool changed = true;
    int max_iters = 20;
    int iter = 0;

    while (changed && iter < max_iters) {
        changed = false;
        int n = pts.size();
        if (n < 5) break;

        for (int i = 0; i < n; ++i) {
            // 提取连续的 5 个点：a, b, c, d, e
            int i_a = (i - 1 + n) % n;
            int i_b = i;
            int i_c = (i + 1) % n;
            int i_d = (i + 2) % n;
            int i_e = (i + 3) % n;

            pcl::PointXYZ a = pts[i_a];
            pcl::PointXYZ b = pts[i_b];
            pcl::PointXYZ c = pts[i_c];
            pcl::PointXYZ d = pts[i_d];
            pcl::PointXYZ e = pts[i_e];

            // 计算两条中间可疑短边的长度
            double bc_len = std::hypot(c.x - b.x, c.y - b.y);
            double cd_len = std::hypot(d.x - c.x, d.y - c.y);

            // 条件1：bc 和 cd 都是短边
            if (bc_len > 0.01 && bc_len < threshold &&
                cd_len > 0.01 && cd_len < threshold) {

                Eigen::Vector2d v_ab(b.x - a.x, b.y - a.y);
                Eigen::Vector2d v_de(e.x - d.x, e.y - d.y);
                double len_ab = v_ab.norm();
                double len_de = v_de.norm();

                // 条件2：两端的 ab 和 de 必须相对较长，确保它们是主墙面
                if (/*len_ab > threshold * 1.2 && len_de > threshold * 1.2*/len_ab>bc_len && len_de>cd_len && len_ab>cd_len && len_de > bc_len) {

                    Eigen::Vector2d dir_ab = v_ab / len_ab;
                    Eigen::Vector2d dir_de = v_de / len_de;

                    // 条件3：两端墙面近似垂直（点积 < 0.5 约等于夹角 60°~120°）
                    double dot = std::abs(dir_ab.dot(dir_de));
                    if (dot < 0.5) {

                        // 计算直线 ab 的方程 A1*x + B1*y = C1
                        double A1 = b.y - a.y; double B1 = a.x - b.x;
                        double C1 = A1 * a.x + B1 * a.y;

                        // 计算直线 de 的方程 A2*x + B2*y = C2
                        double A2 = e.y - d.y; double B2 = d.x - e.x;
                        double C2 = A2 * d.x + B2 * d.y;

                        double det = A1 * B2 - A2 * B1;

                        if (std::abs(det) > 0.01) {
                            pcl::PointXYZ p;
                            p.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
                            p.y = static_cast<float>((A1 * C2 - A2 * C1) / det);
                            p.z = b.z;

                            // 条件4：防过度矫正（交点不能离原本的拐角区域太远）
                            double max_dist = threshold * 2.5;
                            if (std::hypot(p.x - b.x, p.y - b.y) < max_dist &&
                                std::hypot(p.x - c.x, p.y - c.y) < max_dist &&
                                std::hypot(p.x - d.x, p.y - d.y) < max_dist) {

                                // 几何拓扑修改：用交点 p 替换 b，然后删除 c 和 d
                                pts[i_b] = p;

                                // 极其重要：删除两个元素时，必须先删除索引大的，再删除索引小的，防止越界或错位
                                int max_idx = std::max(i_c, i_d);
                                int min_idx = std::min(i_c, i_d);
                                pts.erase(pts.begin() + max_idx);
                                pts.erase(pts.begin() + min_idx);

                                changed = true;
                                break; // 拓扑改变，跳出当前循环，重新计算 n 并遍历
                            }
                        }
                    }
                }
            }
        }
        iter++;
    }
}

// ===== resolveGenericShortEdge =====
// 作用：通用短边折叠——短边夹在两条长边（3倍以上）之间时直接延长求交，
//       不判断角度类型，作为所有模式匹配的兜底
// 参数：pts - 多边形顶点（就地修改）; threshold - 短边判定阈值
void outlineRegular::resolveGenericShortEdge(std::vector<pcl::PointXYZ>& pts, double threshold)
{
    if (pts.size() < 4) return;

    bool changed = true;
    int max_iters = 20;
    int iter = 0;

    while (changed && iter < max_iters) {
        changed = false;
        int n = pts.size();
        if (n < 4) break;

        for (int i = 0; i < n; ++i) {
            int i_a = (i - 1 + n) % n;
            int i_b = i;
            int i_c = (i + 1) % n;
            int i_d = (i + 2) % n;

            pcl::PointXYZ a = pts[i_a];
            pcl::PointXYZ b = pts[i_b];
            pcl::PointXYZ c = pts[i_c];
            pcl::PointXYZ d = pts[i_d];

            double bc_len = std::hypot(c.x - b.x, c.y - b.y);
            double ab_len = std::hypot(b.x - a.x, b.y - a.y);
            double cd_len = std::hypot(d.x - c.x, d.y - c.y);

            // 条件：bc是短边，且两侧边都明显长于它（3倍以上）
            if (bc_len > 0.01 && bc_len < threshold &&
                ab_len > bc_len * 3.0 && cd_len > bc_len * 3.0) {

                double A1 = b.y - a.y; double B1 = a.x - b.x;
                double C1 = A1 * a.x + B1 * a.y;
                double A2 = d.y - c.y; double B2 = c.x - d.x;
                double C2 = A2 * c.x + B2 * c.y;
                double det = A1 * B2 - A2 * B1;

                if (std::abs(det) > 0.01) {
                    pcl::PointXYZ p;
                    p.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
                    p.y = static_cast<float>((A1 * C2 - A2 * C1) / det);
                    p.z = b.z;

                    if (std::hypot(p.x - b.x, p.y - b.y) < threshold * 2.0 &&
                        std::hypot(p.x - c.x, p.y - c.y) < threshold * 2.0) {
                        pts[i_b] = p;
                        pts.erase(pts.begin() + i_c);
                        changed = true;
                        break;
                    }
                }
            }
        }
        iter++;
    }
}

// ===== resolveShortEdgeCluster =====
// 作用：处理由2~4条连续短边组成的尖刺/倒角/小凹槽。相比单短边修复，
//       这里把整段短边簇替换为两侧主边的交点，专门覆盖连续锯齿漏检场景。
void outlineRegular::resolveShortEdgeCluster(std::vector<pcl::PointXYZ>& pts, double threshold)
{
    if (pts.size() < 5) return;

    bool changed = true;
    int rounds = 0;
    const int max_rounds = 12;

    auto dist2d = [](const pcl::PointXYZ& p, const pcl::PointXYZ& q) {
        return std::hypot(static_cast<double>(p.x) - q.x, static_cast<double>(p.y) - q.y);
    };

    auto lineIntersection = [](const pcl::PointXYZ& a,
                               const pcl::PointXYZ& b,
                               const pcl::PointXYZ& c,
                               const pcl::PointXYZ& d,
                               pcl::PointXYZ& out) {
        const double A1 = b.y - a.y;
        const double B1 = a.x - b.x;
        const double C1 = A1 * a.x + B1 * a.y;
        const double A2 = d.y - c.y;
        const double B2 = c.x - d.x;
        const double C2 = A2 * c.x + B2 * c.y;
        const double det = A1 * B2 - A2 * B1;
        if (std::abs(det) < 0.01) return false;

        out.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
        out.y = static_cast<float>((A1 * C2 - A2 * C1) / det);
        out.z = b.z;
        return true;
    };

    auto pointLineDistance = [](const pcl::PointXYZ& p,
                                const pcl::PointXYZ& a,
                                const pcl::PointXYZ& b) {
        const double dx = static_cast<double>(b.x) - a.x;
        const double dy = static_cast<double>(b.y) - a.y;
        const double length = std::hypot(dx, dy);
        if (length < 1e-9) return std::numeric_limits<double>::infinity();
        return std::abs(dx * (a.y - p.y) - (a.x - p.x) * dy) / length;
    };

    while (changed && rounds < max_rounds) {
        changed = false;
        const int n = static_cast<int>(pts.size());
        if (n < 5) break;

        // A spike can consist of one short return edge and one considerably
        // longer edge.  Such a pair is intentionally not an "all short"
        // cluster, so detect its apex from the detour ratio and sharp reversal.
        for (int i = 0; i < n && !changed; ++i) {
            const int i_prev = (i - 1 + n) % n;
            const int i_next = (i + 1) % n;
            const pcl::PointXYZ& prev = pts[i_prev];
            const pcl::PointXYZ& apex = pts[i];
            const pcl::PointXYZ& next = pts[i_next];

            const double incoming = dist2d(prev, apex);
            const double outgoing = dist2d(apex, next);
            const double chord = dist2d(prev, next);
            if (incoming < 0.01 || outgoing < 0.01 || chord < 0.01) continue;
            if (std::min(incoming, outgoing) > threshold) continue;
            if (std::max(incoming, outgoing) > threshold * 2.0) continue;
            if (chord > threshold * 1.5) continue;

            const double ux = (apex.x - prev.x) / incoming;
            const double uy = (apex.y - prev.y) / incoming;
            const double vx = (next.x - apex.x) / outgoing;
            const double vy = (next.y - apex.y) / outgoing;
            const double turn = std::abs(std::atan2(ux * vy - uy * vx,
                                                    ux * vx + uy * vy));
            const double detour_ratio = (incoming + outgoing) / chord;
            if (turn < 135.0 * M_PI / 180.0 || detour_ratio < 1.45) continue;

            std::vector<pcl::PointXYZ> candidate = pts;
            candidate.erase(candidate.begin() + i);
            if (candidate.size() < 3 || !isSimplePolygon2D(candidate)) continue;

            const double old_area = polygonArea2D(pts);
            const double new_area = polygonArea2D(candidate);
            if (old_area > 1e-6) {
                const double ratio = new_area / old_area;
                if (ratio < 0.92 || ratio > 1.08) continue;
            }

            pts.swap(candidate);
            changed = true;
            std::cerr << "[Topo] removed reflex spike apex turn="
                      << turn * 180.0 / M_PI
                      << " detour_ratio=" << detour_ratio << std::endl;
        }
        if (changed) {
            ++rounds;
            continue;
        }

        for (int start = 0; start < n && !changed; ++start) {
            std::vector<pcl::PointXYZ> rotated;
            rotated.reserve(pts.size());
            for (int k = 0; k < n; ++k) {
                rotated.push_back(pts[(start + k) % n]);
            }

            for (int edge_count = 4; edge_count >= 2 && !changed; --edge_count) {
                if (static_cast<int>(rotated.size()) < edge_count + 3) continue;

                const pcl::PointXYZ& a = rotated[0];
                const pcl::PointXYZ& b = rotated[1];
                const pcl::PointXYZ& last = rotated[1 + edge_count];
                const pcl::PointXYZ& d = rotated[2 + edge_count];

                const double len_ab = dist2d(a, b);
                const double len_tail = dist2d(last, d);
                if (len_ab < 0.01 || len_tail < 0.01) continue;

                double chain_len = 0.0;
                double max_edge = 0.0;
                bool all_short = true;
                for (int e = 1; e <= edge_count; ++e) {
                    const double len = dist2d(rotated[e], rotated[e + 1]);
                    chain_len += len;
                    max_edge = std::max(max_edge, len);
                    if (len > threshold || len < 0.01) {
                        all_short = false;
                        break;
                    }
                }
                if (!all_short) continue;

                const double side_min = std::min(len_ab, len_tail);
                if (side_min < std::max(max_edge * 1.2, threshold * 0.6)) continue;
                if (chain_len > threshold * 3.0) continue;

                pcl::PointXYZ p;
                const double ux = (b.x - a.x) / len_ab;
                const double uy = (b.y - a.y) / len_ab;
                const double vx = (d.x - last.x) / len_tail;
                const double vy = (d.y - last.y) / len_tail;
                const double side_cross = std::abs(ux * vy - uy * vx);

                // A short staircase may lie between two portions of the same
                // wall.  Their supporting lines are parallel, so replacing the
                // chain by an intersection is impossible; remove the interior
                // bump and reconnect its two endpoints instead.
                if (side_cross < std::sin(15.0 * M_PI / 180.0)) {
                    double max_offset = 0.0;
                    for (int v = 2; v <= edge_count; ++v) {
                        max_offset = std::max(max_offset,
                            pointLineDistance(rotated[v], b, last));
                    }
                    if (max_offset > threshold * 0.75) continue;

                    std::vector<pcl::PointXYZ> candidate = rotated;
                    candidate.erase(candidate.begin() + 2,
                                    candidate.begin() + 1 + edge_count);
                    if (candidate.size() < 3 || !isSimplePolygon2D(candidate)) continue;

                    const double old_area = polygonArea2D(rotated);
                    const double new_area = polygonArea2D(candidate);
                    if (old_area > 1e-6) {
                        const double ratio = new_area / old_area;
                        if (ratio < 0.92 || ratio > 1.08) continue;
                    }

                    pts.swap(candidate);
                    changed = true;
                    std::cerr << "[Topo] collapsed parallel micro-feature edges="
                              << edge_count << " offset=" << max_offset << std::endl;
                    continue;
                }

                if (!lineIntersection(a, b, last, d, p)) continue;

                double max_drift = 0.0;
                for (int v = 1; v <= edge_count + 1; ++v) {
                    max_drift = std::max(max_drift, dist2d(p, rotated[v]));
                }
                if (max_drift > threshold * 2.4 || max_drift > 8.0) continue;

                std::vector<pcl::PointXYZ> candidate = rotated;
                candidate[1] = p;
                candidate.erase(candidate.begin() + 2, candidate.begin() + 2 + edge_count);
                if (candidate.size() < 3) continue;
                if (!isSimplePolygon2D(candidate)) continue;

                const double old_area = polygonArea2D(rotated);
                const double new_area = polygonArea2D(candidate);
                if (old_area > 1e-6) {
                    const double ratio = new_area / old_area;
                    if (ratio < 0.85 || ratio > 1.15) continue;
                }

                pts.swap(candidate);
                changed = true;
                std::cerr << "[Topo] collapsed corner micro-feature edges="
                          << edge_count << " drift=" << max_drift << std::endl;
            }
        }

        ++rounds;
    }
}

// ===== resolveAngleAnomalies =====
// 作用：处理不一定很短、但角度明显异常的尖角/小缺口。
//       现有短边修复主要依赖边长阈值，容易漏掉 3~6m 的中等长度异常折角。
void outlineRegular::resolveAngleAnomalies(std::vector<pcl::PointXYZ>& pts, double threshold)
{
    if (pts.size() < 5) return;

    auto dist2d = [](const pcl::PointXYZ& a, const pcl::PointXYZ& b) {
        return std::hypot(static_cast<double>(a.x) - b.x,
                          static_cast<double>(a.y) - b.y);
    };

    const double original_area = polygonArea2D(pts);
    const double scale = std::sqrt(std::max(original_area, 1.0));
    const double max_leg = clampDouble(std::max(threshold * 1.6, 0.06 * scale), threshold, 12.0);
    const double max_chord = clampDouble(std::max(threshold * 2.0, 0.08 * scale), threshold * 1.2, 14.0);
    const double max_area_ratio_delta = 0.12;

    bool changed = true;
    int rounds = 0;
    while (changed && rounds < 12 && pts.size() >= 5) {
        changed = false;
        const int n = static_cast<int>(pts.size());
        for (int i = 0; i < n; ++i) {
            const int i_prev = (i - 1 + n) % n;
            const int i_next = (i + 1) % n;
            const auto& prev = pts[i_prev];
            const auto& apex = pts[i];
            const auto& next = pts[i_next];

            const double in_len = dist2d(prev, apex);
            const double out_len = dist2d(apex, next);
            const double chord = dist2d(prev, next);
            if (in_len < 1e-6 || out_len < 1e-6 || chord < 1e-6) continue;
            if (std::min(in_len, out_len) > max_leg) continue;
            if (chord > max_chord) continue;

            Eigen::Vector2d a(prev.x - apex.x, prev.y - apex.y);
            Eigen::Vector2d b(next.x - apex.x, next.y - apex.y);
            const double cos_angle = clampDouble(a.dot(b) / (a.norm() * b.norm()), -1.0, 1.0);
            const double angle_deg = std::acos(cos_angle) * 180.0 / M_PI;
            const double detour_ratio = (in_len + out_len) / chord;

            const bool sharp_spike = angle_deg < 55.0 && detour_ratio > 1.10;
            const bool narrow_notch = angle_deg < 75.0 && detour_ratio > 1.35 &&
                std::min(in_len, out_len) < threshold * 1.8;
            const bool small_bevel = std::abs(angle_deg - 135.0) < 18.0 &&
                std::min(in_len, out_len) < threshold * 1.25 &&
                std::max(in_len, out_len) < max_leg;
            if (!sharp_spike && !narrow_notch && !small_bevel) continue;

            std::vector<pcl::PointXYZ> candidate = pts;
            candidate.erase(candidate.begin() + i);
            if (candidate.size() < 3 || !isSimplePolygon2D(candidate)) continue;

            const double new_area = polygonArea2D(candidate);
            if (original_area > 1e-6) {
                const double delta = std::abs(new_area - polygonArea2D(pts)) / original_area;
                if (delta > max_area_ratio_delta) continue;
            }

            pts.swap(candidate);
            changed = true;
            std::cerr << "[Topo] removed angle anomaly angle=" << angle_deg
                << " detour_ratio=" << detour_ratio
                << " in=" << in_len
                << " out=" << out_len
                << std::endl;
            break;
        }
        ++rounds;
    }
}

// ===== repairTopologyBatch =====
// 作用：批量评分排序+冲突解决框架，将所有拓扑修复模式（倒角/尖刺/双短边/
//       平行阶梯/平行凹槽/通用折叠6种）统一为：候选收集 -> 安全过滤 ->
//       评分排序(越短越优先) -> 冲突解决 -> 从后往前批量执行
// 参数：pts - 多边形顶点（就地修改）; threshold - 修复距离阈值
void outlineRegular::repairTopologyBatch(std::vector<pcl::PointXYZ>& pts, double threshold)
{
    if (pts.size() < 4) return;

    // 内部结构体：记录一个候选拓扑修复（类型+涉及的顶点索引+新顶点+评分）
    struct RepairCandidate {
        enum Type { CHAMFER, SPIKE, DOUBLE_SPIKE, PARALLEL_STEP, PARALLEL_NOTCH, GENERIC };
        Type type;
        int idx_b;
        int idx_c;
        int idx_e;  // 额外索引：DOUBLE_SPIKE时为idx_d，PARALLEL_STEP时=-1
        pcl::PointXYZ new_vertex;  // 替换后的新顶点
        pcl::PointXYZ new_vertex2; // PARALLEL_STEP时c_new
        double score;
    };

    bool changed = true;
    int max_rounds = 20;
    int round = 0;

    while (changed && round < max_rounds) {
        changed = false;
        int n = pts.size();
        if (n < 4) break;

        std::vector<RepairCandidate> candidates;

        // === 扫描：收集所有候选修复 ===
        for (int i = 0; i < n; ++i) {
            int i_a = (i - 1 + n) % n;
            int i_b = i;
            int i_c = (i + 1) % n;
            int i_d = (i + 2) % n;
            int i_e = (i + 3) % n;

            pcl::PointXYZ a = pts[i_a];
            pcl::PointXYZ b = pts[i_b];
            pcl::PointXYZ c = pts[i_c];
            pcl::PointXYZ d = pts[i_d];

            double bc_len = std::hypot(c.x - b.x, c.y - b.y);
            double ab_len = std::hypot(b.x - a.x, b.y - a.y);
            double cd_len = std::hypot(d.x - c.x, d.y - c.y);

            if (bc_len < 0.01) continue;

            Eigen::Vector2d ba(a.x - b.x, a.y - b.y);
            Eigen::Vector2d bc_vec(c.x - b.x, c.y - b.y);
            Eigen::Vector2d cb(b.x - c.x, b.y - c.y);
            Eigen::Vector2d cd_vec(d.x - c.x, d.y - c.y);

            // --- 模式1：倒角(Chamfer) — 至少一侧钝角 ---
            if (bc_len < threshold && ab_len > threshold * 1.5 && cd_len > threshold * 1.5) {
                bool oneSideObtuse = (ba.dot(bc_vec) < -1e-3) || (cd_vec.dot(cb) < -1e-3);
                bool notNormalConvex = !(ba.dot(bc_vec) > 0.3 && cd_vec.dot(cb) > 0.3);

                if (oneSideObtuse && notNormalConvex) {
                    double A1 = b.y - a.y; double B1 = a.x - b.x;
                    double C1 = A1 * a.x + B1 * a.y;
                    double A2 = d.y - c.y; double B2 = c.x - d.x;
                    double C2 = A2 * c.x + B2 * c.y;
                    double det = A1 * B2 - A2 * B1;

                    if (std::abs(det) > 0.01) {
                        pcl::PointXYZ p;
                        p.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
                        p.y = static_cast<float>((A1 * C2 - A2 * C1) / det);
                        p.z = b.z;

                        if (std::hypot(p.x - b.x, p.y - b.y) < bc_len * 3.0 &&
                            std::hypot(p.x - c.x, p.y - c.y) < bc_len * 3.0) {
                            candidates.push_back({ RepairCandidate::CHAMFER, i_b, i_c, -1,
                                                  p, pcl::PointXYZ(), bc_len });
                        }
                    }
                }
            }

            // --- 模式2：尖刺(Spike) — 不再要求正交 ---
            if (bc_len < threshold && ab_len > bc_len && cd_len > bc_len) {
                bool hasSharpCorner = (ba.dot(bc_vec) >= -1e-3) || (cd_vec.dot(cb) >= -1e-3);

                if (hasSharpCorner) {
                    double A1 = b.y - a.y; double B1 = a.x - b.x;
                    double C1 = A1 * a.x + B1 * a.y;
                    double A2 = d.y - c.y; double B2 = c.x - d.x;
                    double C2 = A2 * c.x + B2 * c.y;
                    double det = A1 * B2 - A2 * B1;

                    if (std::abs(det) > 0.01) {
                        pcl::PointXYZ p;
                        p.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
                        p.y = static_cast<float>((A1 * C2 - A2 * C1) / det);
                        p.z = b.z;

                        if (std::hypot(p.x - b.x, p.y - b.y) < threshold * 3.0 &&
                            std::hypot(p.x - c.x, p.y - c.y) < threshold * 3.0) {
                            candidates.push_back({ RepairCandidate::SPIKE, i_b, i_c, -1,
                                                  p, pcl::PointXYZ(), bc_len });
                        }
                    }
                }
            }

            // --- 模式3：双连续短边(Double Spike) ---
            if (n >= 5) {
                pcl::PointXYZ e = pts[i_e];
                double cd_len2 = std::hypot(d.x - c.x, d.y - c.y);

                if (bc_len > 0.01 && bc_len < threshold &&
                    cd_len2 > 0.01 && cd_len2 < threshold) {

                    Eigen::Vector2d v_ab(b.x - a.x, b.y - a.y);
                    Eigen::Vector2d v_de(e.x - d.x, e.y - d.y);
                    double len_ab2 = v_ab.norm();
                    double len_de2 = v_de.norm();

                    if (len_ab2 > bc_len && len_de2 > cd_len2 &&
                        len_ab2 > cd_len2 && len_de2 > bc_len) {

                        Eigen::Vector2d dir_ab = v_ab / len_ab2;
                        Eigen::Vector2d dir_de = v_de / len_de2;
                        double dot = std::abs(dir_ab.dot(dir_de));

                        if (dot < 0.5) {
                            double A1 = b.y - a.y; double B1 = a.x - b.x;
                            double C1 = A1 * a.x + B1 * a.y;
                            double A2 = e.y - d.y; double B2 = d.x - e.x;
                            double C2 = A2 * d.x + B2 * d.y;
                            double det = A1 * B2 - A2 * B1;

                            if (std::abs(det) > 0.01) {
                                pcl::PointXYZ p;
                                p.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
                                p.y = static_cast<float>((A1 * C2 - A2 * C1) / det);
                                p.z = b.z;

                                double max_dist = threshold * 2.5;
                                if (std::hypot(p.x - b.x, p.y - b.y) < max_dist &&
                                    std::hypot(p.x - c.x, p.y - c.y) < max_dist &&
                                    std::hypot(p.x - d.x, p.y - d.y) < max_dist) {
                                    candidates.push_back({ RepairCandidate::DOUBLE_SPIKE,
                                        i_b, i_c, i_d, p, pcl::PointXYZ(), bc_len + cd_len2 });
                                }
                            }
                        }
                    }
                }
            }

            // --- 模式4：平行错位阶梯(Parallel Step) ---
            if (bc_len > 0.01 && bc_len < threshold) {
                Eigen::Vector2d v_ab(b.x - a.x, b.y - a.y);
                Eigen::Vector2d v_cd(d.x - c.x, d.y - c.y);
                double len_ab2 = v_ab.norm();
                double len_cd2 = v_cd.norm();

                if (len_ab2 > bc_len * 1.0 && len_cd2 > bc_len * 1.0) {
                    Eigen::Vector2d dir_ab = v_ab / len_ab2;
                    Eigen::Vector2d dir_cd = v_cd / len_cd2;

                    if (dir_ab.dot(dir_cd) > 0.866) {
                        // 计算a和d到直线bc的有符号距离，判断同侧/异侧
                        Eigen::Vector2d bc_vec(c.x - b.x, c.y - b.y);
                        double bc_n = bc_vec.norm();
                        if (bc_n > 1e-6) {
                            Eigen::Vector2d bc_dir = bc_vec / bc_n;
                            Eigen::Vector2d bc_normal(-bc_dir.y(), bc_dir.x());
                            double signed_a = (a.x - b.x) * bc_normal.x() + (a.y - b.y) * bc_normal.y();
                            double signed_d = (d.x - b.x) * bc_normal.x() + (d.y - b.y) * bc_normal.y();
                            bool same_side = (signed_a * signed_d > 0);

                            if (!same_side) {
                                // === 模式4(异侧): 沿墙面方向投影 ===
                                Eigen::Vector2d dir = (dir_ab + dir_cd).normalized();
                                pcl::PointXYZ M;
                                M.x = (b.x + c.x) / 2.0f;
                                M.y = (b.y + c.y) / 2.0f;
                                M.z = b.z;

                                double t_b = (Eigen::Vector2d(M.x - b.x, M.y - b.y).dot(dir)) / (dir_ab.dot(dir));
                                pcl::PointXYZ b_new;
                                b_new.x = static_cast<float>(b.x + t_b * dir_ab.x());
                                b_new.y = static_cast<float>(b.y + t_b * dir_ab.y());
                                b_new.z = b.z;

                                double t_c = (Eigen::Vector2d(M.x - c.x, M.y - c.y).dot(dir)) / (dir_cd.dot(dir));
                                pcl::PointXYZ c_new;
                                c_new.x = static_cast<float>(c.x + t_c * dir_cd.x());
                                c_new.y = static_cast<float>(c.y + t_c * dir_cd.y());
                                c_new.z = c.z;

                                if (std::hypot(b_new.x - b.x, b_new.y - b.y) < bc_len * 2.0 &&
                                    std::hypot(c_new.x - c.x, c_new.y - c.y) < bc_len * 2.0) {
                                    candidates.push_back({ RepairCandidate::PARALLEL_STEP,
                                        i_b, i_c, -1, b_new, c_new, bc_len });
                                }
                            }
                            else {
                                // === 模式6(同侧): 沿垂线方向投影 ===
                                Eigen::Vector2d wall_dir = (dir_ab + dir_cd).normalized();
                                Eigen::Vector2d perp_dir(-wall_dir.y(), wall_dir.x());

                                // b投影到cd所在直线
                                double A_cd = d.y - c.y, B_cd = c.x - d.x;
                                double C_cd = A_cd * c.x + B_cd * c.y;
                                double denom_b = A_cd * perp_dir.x() + B_cd * perp_dir.y();
                                if (std::abs(denom_b) > 0.01) {
                                    double t_bp = (C_cd - A_cd * b.x - B_cd * b.y) / denom_b;
                                    pcl::PointXYZ b_proj;
                                    b_proj.x = static_cast<float>(b.x + t_bp * perp_dir.x());
                                    b_proj.y = static_cast<float>(b.y + t_bp * perp_dir.y());
                                    b_proj.z = b.z;

                                    // c投影到ab所在直线
                                    double A_ab = b.y - a.y, B_ab = a.x - b.x;
                                    double C_ab = A_ab * a.x + B_ab * a.y;
                                    double denom_c = A_ab * perp_dir.x() + B_ab * perp_dir.y();
                                    if (std::abs(denom_c) > 0.01) {
                                        double t_cp = (C_ab - A_ab * c.x - B_ab * c.y) / denom_c;
                                        pcl::PointXYZ c_proj;
                                        c_proj.x = static_cast<float>(c.x + t_cp * perp_dir.x());
                                        c_proj.y = static_cast<float>(c.y + t_cp * perp_dir.y());
                                        c_proj.z = c.z;

                                        if (std::hypot(b_proj.x - b.x, b_proj.y - b.y) < bc_len * 2.0 &&
                                            std::hypot(c_proj.x - c.x, c_proj.y - c.y) < bc_len * 2.0) {
                                            candidates.push_back({ RepairCandidate::PARALLEL_NOTCH,
                                                i_b, i_c, -1, b_proj, c_proj, bc_len });
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // --- 模式5：通用短边折叠(Generic) — 兜底 ---
            if (bc_len > 0.01 && bc_len < threshold &&
                ab_len > bc_len * 3.0 && cd_len > bc_len * 3.0) {

                double A1 = b.y - a.y; double B1 = a.x - b.x;
                double C1 = A1 * a.x + B1 * a.y;
                double A2 = d.y - c.y; double B2 = c.x - d.x;
                double C2 = A2 * c.x + B2 * c.y;
                double det = A1 * B2 - A2 * B1;

                if (std::abs(det) > 0.01) {
                    pcl::PointXYZ p;
                    p.x = static_cast<float>((C1 * B2 - C2 * B1) / det);
                    p.y = static_cast<float>((A1 * C2 - A2 * C1) / det);
                    p.z = b.z;

                    if (std::hypot(p.x - b.x, p.y - b.y) < threshold * 2.0 &&
                        std::hypot(p.x - c.x, p.y - c.y) < threshold * 2.0) {
                        candidates.push_back({ RepairCandidate::GENERIC, i_b, i_c, -1,
                                              p, pcl::PointXYZ(), bc_len * 2.0 });
                    }
                }
            }
        }

        if (candidates.empty()) break;

        // === 安全性过滤：排除新顶点离原始位置过远的候选 ===
        // 防止近平行线求交产生的飞点，绝对上限5.0m
        const double MAX_ABS_DRIFT = 5.0;
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
            [&](const RepairCandidate& c) {
                pcl::PointXYZ b_orig = pts[c.idx_b];
                pcl::PointXYZ c_orig = pts[c.idx_c];
                double drift_b = std::hypot(c.new_vertex.x - b_orig.x, c.new_vertex.y - b_orig.y);
                if (c.type == RepairCandidate::PARALLEL_STEP || c.type == RepairCandidate::PARALLEL_NOTCH) {
                    double drift_c = std::hypot(c.new_vertex2.x - c_orig.x, c.new_vertex2.y - c_orig.y);
                    return drift_b > MAX_ABS_DRIFT || drift_c > MAX_ABS_DRIFT;
                }
                double drift_c = std::hypot(c.new_vertex.x - c_orig.x, c.new_vertex.y - c_orig.y);
                return drift_b > MAX_ABS_DRIFT || drift_c > MAX_ABS_DRIFT;
            }), candidates.end());

        if (candidates.empty()) break;

        // === 按分数排序（越短越优先） ===
        std::sort(candidates.begin(), candidates.end(),
            [](const RepairCandidate& lhs, const RepairCandidate& rhs) {
                return lhs.score < rhs.score;
            });

        // === 冲突解决 + 批量执行 ===
        std::set<int> occupied;
        // 收集需执行的修复（按顶点冲突过滤）
        std::vector<const RepairCandidate*> accepted;
        for (const auto& cand : candidates) {
            bool conflict = false;
            if (cand.type == RepairCandidate::DOUBLE_SPIKE) {
                if (occupied.count(cand.idx_b) || occupied.count(cand.idx_c) || occupied.count(cand.idx_e))
                    conflict = true;
            }
            else if (cand.type == RepairCandidate::PARALLEL_STEP || cand.type == RepairCandidate::PARALLEL_NOTCH) {
                if (occupied.count(cand.idx_b) || occupied.count(cand.idx_c))
                    conflict = true;
            }
            else {
                if (occupied.count(cand.idx_b) || occupied.count(cand.idx_c))
                    conflict = true;
            }
            if (conflict) continue;

            if (cand.type == RepairCandidate::DOUBLE_SPIKE) {
                occupied.insert(cand.idx_b); occupied.insert(cand.idx_c); occupied.insert(cand.idx_e);
            }
            else if (cand.type == RepairCandidate::PARALLEL_STEP || cand.type == RepairCandidate::PARALLEL_NOTCH) {
                occupied.insert(cand.idx_b); occupied.insert(cand.idx_c);
            }
            else {
                occupied.insert(cand.idx_b); occupied.insert(cand.idx_c);
            }
            accepted.push_back(&cand);
        }

        if (accepted.empty()) break;

        // === 从后往前执行（按delete索引降序，防止错位） ===
        std::sort(accepted.begin(), accepted.end(),
            [](const RepairCandidate* lhs, const RepairCandidate* rhs) {
                int max_l = std::max(lhs->idx_c, lhs->idx_e);
                int max_r = std::max(rhs->idx_c, rhs->idx_e);
                return max_l > max_r;
            });

        for (const auto* cand : accepted) {
            if (cand->type == RepairCandidate::DOUBLE_SPIKE) {
                int max_del = std::max(cand->idx_c, cand->idx_e);
                int min_del = std::min(cand->idx_c, cand->idx_e);
                int sz = (int)pts.size();
                if (max_del < sz) pts.erase(pts.begin() + max_del);
                if (min_del < (int)pts.size()) pts.erase(pts.begin() + min_del);
                if (cand->idx_b < (int)pts.size()) pts[cand->idx_b] = cand->new_vertex;
            }
            else if (cand->type == RepairCandidate::PARALLEL_STEP || cand->type == RepairCandidate::PARALLEL_NOTCH) {
                if (cand->idx_b < (int)pts.size()) pts[cand->idx_b] = cand->new_vertex;
                if (cand->idx_c < (int)pts.size()) pts[cand->idx_c] = cand->new_vertex2;
            }
            else {
                if (cand->idx_b < (int)pts.size()) pts[cand->idx_b] = cand->new_vertex;
                if (cand->idx_c < (int)pts.size() && cand->idx_c != cand->idx_b)
                    pts.erase(pts.begin() + cand->idx_c);
            }
        }

        changed = true;
        round++;
    }
}

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
#include <filesystem>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
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
// 超过此幅度的墙面纠偏视为可疑：墙面证据稀疏时，
// 直方图峰可能来自一小段异常墙碎片
// (实测：仅 400-8000 配对造成 15-45° 的"纠偏"，把建筑
// 转离轮廓/假设/墙面主体)。真正的纠偏
// (已验证案例：32.1°→16.5°)不会超过此上限。
constexpr double kWallDirectionCorrectionMaxDeg = 25.0;
// 原始轮廓 PCA 门控：栅格楼梯轮廓的整体趋势
// 是拉长轮廓的稳定方向先验(真实数据上验证到 ~0.5°)。
// 方向选择偏离先验超过门槛时必须有
// 强墙面证据，否则回退到轮廓趋势。用于拦截
// "VDP 凸包对角线"失效模式(12°建筑被选成18°)。
constexpr double kPcaDirectionReliableAxisRatio = 1.6; // elongation (std ratio) for the PCA fallback gate
constexpr double kPcaDirectionGateStrongDeg = 5.0;
// 小建筑走正常管线(main.cpp 的方向矩形快通道
// 已移除)，但其最优假设必须按单一主方向
// 规则化：多方向规则化产生斜边，
// 在小轮廓上明显不对。数值沿用原
// main.cpp 矩形快通道阈值。
constexpr double kSmallBuildingSingleDirectionArea = 60.0;
// 证据支撑的缺口保护：真实矩形缺口(退台/院湾)
// 在能量假设中保留，却被正交清理链抹掉
// (repair_distance 尺度的拓扑修复、forceOrthogonal-
// PolygonToAngle 的同轴删点)。只有当缺口足够深
// (是建筑特征而非掩膜噪声)且网格
// 支撑点云确认缺口壁上有支撑时才保护。
constexpr double kNotchMinDepthFactor = 0.6;    // x tuning.repair_distance
constexpr double kNotchMinDepthFloor = 1.0;     // meters
constexpr double kNotchMaxWidth = 15.0;         // meters; larger concavities are wings, not notches
constexpr double kNotchSupportRadius = 0.8;     // meters around the notch's 3 edges
constexpr int kNotchMinSupportPoints = 8;
constexpr double kNotchPreserveMidpointTol = 1.2;  // meters, bottom-edge midpoint match
constexpr double kNotchPreserveAngleDeg = 20.0;
constexpr double kNotchPreserveLengthRatio = 0.5;
// ---- Direction chain extraction & multi-peak detection ----
// Split-and-merge continuous edge chain extraction on the closed ring,
// followed by weighted KDE multi-peak detection in [0°, 90°).
constexpr double kDirectionMinChainLength = 0.8;        // m, chains shorter than this don't vote
constexpr double kDirectionMaxAngularResidualDeg = 6.0; // deg, max point-line angular residual for chain fit
constexpr double kDirectionPeakSeparationDeg = 12.0;    // deg, minimum separation between direction peaks
constexpr double kDirectionKdeBandwidthDeg = 4.0;       // deg, KDE smoothing bandwidth
constexpr double kDirectionMinSecondaryRatio = 0.18;    // secondary peak must have >=18% of total weight
constexpr int kDirectionMinIndependentChains = 2;       // secondary peak needs >=2 spatially independent chains
// ---- Structure-aware hypothesis repair (single-spike -> rectangle) ----
// VDP encodes a real rectangular notch/protrusion as ONE extreme vertex (a
// spike); the energy's vertex-count model term actively prefers that
// degenerate encoding for small structures. This stage inspects each spike
// against the ORIGINAL outline points in a local chord frame and rebuilds a
// 3-corner rectangular structure when the raw evidence supports it, BEFORE
// direction regularization (so Ceres optimizes the correct edge set).
constexpr double kSpikeMinDepth = 1.2;        // m, spike height over the P->N chord
constexpr double kSpikeMaxDepth = 8.0;
constexpr double kSpikeChordMin = 1.5;        // m, chord length window
constexpr double kSpikeChordMax = 30.0;
constexpr double kRepairMinWidth = 0.8;       // m, structure opening width
constexpr double kRepairMaxWidth = 12.0;
constexpr double kRepairBottomBand = 0.5;     // m, v-band collecting bottom points
constexpr double kRepairWallBand = 0.6;       // m, u-band collecting wall points
constexpr double kRepairMinWallSpan = 1.0;    // m, wall extent along v
constexpr int kRepairMinRawWallPoints = 4;    // original points per wall
constexpr double kRepairSupportRadius = 0.8;  // m, mesh support near the 3 segments
constexpr int kRepairMinSupportPoints = 6;
constexpr double kRepairDataWeight = 1.0;     // gain weights for acceptance
constexpr double kRepairSupportWeight = 0.5;
constexpr double kRepairPerVertexCost = 0.12; // m of gain demanded per added vertex
constexpr int kRepairMaxPerPolygon = 4;       // accepted repairs per hypothesis
constexpr double kRepairMinNewEdge = 0.35;    // m, shortest edge after insertion
constexpr double kRepairMaxAreaChange = 0.06; // relative polygon area change

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

// 作用：含内环多边形的总面积。
double polygonWithHolesArea2D(const Polygon_with_holes_2& polygon)
{
    double area = std::abs(CGAL::to_double(polygon.outer_boundary().area()));
    for (auto hole = polygon.holes_begin(); hole != polygon.holes_end(); ++hole) {
        area -= std::abs(CGAL::to_double(hole->area()));
    }
    return std::max(0.0, area);
}

// 作用：两多边形的 IoU(交并比)。
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

// 作用：单调链法求二维凸包。
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

// 作用：点到方向外接矩形边界的距离。
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

// 作用：计算某方向的 Hausdorff-MBR 候选评分。
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

// 作用：用 Hausdorff-MBR 评分从候选角中选出主方向(最多 max_angles 个)。
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
    // 拓扑伪影更多由轮廓尺度而非采样点距决定。
    // 之前的 0.028*scale 使多数中等建筑
    // 都被钉在 1m 下限上。
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

    // 迭代拟合 2D TLS 直线，剔除鲁棒 MAD 带外的点。
    // 方向门防止邻近立面替换掩膜边。
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

    // 让每条边的残差数量可比，并沿边全长采样。
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

// 折叠角空间([0°,90°), 周期90°)的加权圆均值(四倍角法)。
// 普通算术平均在 1°/89° 这类实际相差 2° 的角对上会得到错误的 45°。
// angles 与 weights 等长; 退化情形(向量和不显著)返回首个角度。
double circularMeanAngle90(const std::vector<double>& angles,
                           const std::vector<double>& weights)
{
    if (angles.empty()) return 0.0;
    double x = 0.0;
    double y = 0.0;
    double total = 0.0;
    for (std::size_t i = 0; i < angles.size() && i < weights.size(); ++i) {
        x += weights[i] * std::cos(4.0 * angles[i]);
        y += weights[i] * std::sin(4.0 * angles[i]);
        total += weights[i];
    }
    if (total <= 1e-9 ||
        (x * x + y * y) <= 1e-12 * total * total) {
        return angles.front();
    }
    const double mean = 0.25 * std::atan2(y, x);
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

    // 从第一条受约束边开始线性化环，合并同向连续段。
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

// ===== Evidence-backed notch protection =====
// A notch = two nearby reflex vertices joined by a short bottom path. Real
// notches (setback bays) carry mesh support points along their walls; mask
// noise teeth do not. Only evidence-backed notches are protected from the
// orthogonal clean-up chain.
struct NotchFeature {
    std::size_t entry = 0;      // reflex vertex index in the source polygon
    std::size_t exit = 0;       // next nearby reflex vertex
    double width = 0.0;         // bottom path length entry->exit
    double depth = 0.0;         // min of the two adjacent wall lengths
    double midX = 0.0, midY = 0.0;
    double dirX = 1.0, dirY = 0.0;   // unit direction along the bottom
};

std::vector<NotchFeature> DetectEvidenceBackedNotches(
    const std::vector<pcl::PointXYZ>& polygon,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support,
    double minDepth,
    double maxWidth)
{
    std::vector<NotchFeature> notches;
    const std::size_t n = polygon.size();
    if (n < 5 || !support || support->empty()) return notches;

    // 按多数转向符号判定凹顶点。
    double turnSum = 0.0;
    std::vector<double> turns(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = polygon[(i + n - 1) % n];
        const auto& b = polygon[i];
        const auto& c = polygon[(i + 1) % n];
        turns[i] = static_cast<double>(b.x - a.x) * (c.y - b.y) -
                   static_cast<double>(b.y - a.y) * (c.x - b.x);
        turnSum += turns[i];
    }
    const double majority = turnSum >= 0.0 ? 1.0 : -1.0;
    std::vector<std::size_t> reflex;
    for (std::size_t i = 0; i < n; ++i) {
        if (turns[i] * majority < 0.0) reflex.push_back(i);
    }
    if (reflex.size() < 2) return notches;

    // 相邻近距凹顶点构成候选段，然后合并
    // 共享顶点的链：VDP 假设常把一个
    // 矩形缺口近似为 2-3 个凹角的链，而每个子缺口
    // 单独看都太小，保护不住。
    struct NotchSpan {
        std::size_t entry;
        std::size_t exit;
    };
    std::vector<NotchSpan> spans;
    for (std::size_t k = 0; k < reflex.size(); ++k) {
        const std::size_t entry = reflex[k];
        const std::size_t exit = reflex[(k + 1) % reflex.size()];
        const std::size_t gap = (exit + n - entry) % n;
        if (gap == 0 || gap > 2) continue;
        if (!spans.empty() && spans.back().exit == entry) {
            spans.back().exit = exit;
        } else {
            spans.push_back({entry, exit});
        }
    }

    for (const auto& span : spans) {
        const std::size_t entry = span.entry;
        const std::size_t exit = span.exit;

        double width = 0.0;
        double midX = 0.0, midY = 0.0;
        int edgeCount = 0;
        for (std::size_t s = entry; s != exit; s = (s + 1) % n) {
            const auto& a = polygon[s];
            const auto& b = polygon[(s + 1) % n];
            width += std::hypot(b.x - a.x, b.y - a.y);
            midX += 0.5 * (a.x + b.x);
            midY += 0.5 * (a.y + b.y);
            ++edgeCount;
        }
        if (edgeCount == 0) continue;
        midX /= edgeCount;
        midY /= edgeCount;
        const double wallIn = std::hypot(
            polygon[entry].x - polygon[(entry + n - 1) % n].x,
            polygon[entry].y - polygon[(entry + n - 1) % n].y);
        const double wallOut = std::hypot(
            polygon[(exit + 1) % n].x - polygon[exit].x,
            polygon[(exit + 1) % n].y - polygon[exit].y);
        const double depth = std::min(wallIn, wallOut);
        if (depth < minDepth) continue;
        if (width < 0.5 || width > maxWidth) continue;

        // 证据：缺口壁和底边附近的支撑点。
        int supportCount = 0;
        for (const auto& p : support->points) {
            double best = std::numeric_limits<double>::max();
            for (std::size_t s = (entry + n - 1) % n; ; s = (s + 1) % n) {
                double t = 0.0;
                best = std::min(best, distancePointToSegmentWithProjection(
                    p, polygon[s], polygon[(s + 1) % n], t));
                if (s == (exit + 1) % n) break;
            }
            if (best <= kNotchSupportRadius) ++supportCount;
        }
        if (supportCount < kNotchMinSupportPoints) continue;

        NotchFeature notch;
        notch.entry = entry;
        notch.exit = exit;
        notch.width = width;
        notch.depth = depth;
        notch.midX = midX;
        notch.midY = midY;
        const double dirLen = std::max(
            std::hypot(static_cast<double>(polygon[exit].x - polygon[entry].x),
                       static_cast<double>(polygon[exit].y - polygon[entry].y)), 1e-9);
        notch.dirX = (polygon[exit].x - polygon[entry].x) / dirLen;
        notch.dirY = (polygon[exit].y - polygon[entry].y) / dirLen;
        notches.push_back(notch);
    }
    return notches;
}

// 候选多边形仍保留与缺口底边对应的边
// (中点+方向+长度匹配)时为真。被抹掉的缺口不满足。
bool NotchPreservedInCandidate(
    const std::vector<pcl::PointXYZ>& candidate,
    const NotchFeature& notch)
{
    const std::size_t n = candidate.size();
    if (n < 4) return false;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = candidate[i];
        const auto& b = candidate[(i + 1) % n];
        const double len = std::hypot(b.x - a.x, b.y - a.y);
        if (len < kNotchPreserveLengthRatio * notch.width) continue;
        const double midX = 0.5 * (a.x + b.x);
        const double midY = 0.5 * (a.y + b.y);
        if (std::hypot(midX - notch.midX, midY - notch.midY) >
                kNotchPreserveMidpointTol) {
            continue;
        }
        const double ux = (b.x - a.x) / len;
        const double uy = (b.y - a.y) / len;
        const double dot = std::abs(ux * notch.dirX + uy * notch.dirY);
        if (dot < std::cos(kNotchPreserveAngleDeg * M_PI / 180.0)) continue;
        return true;
    }
    return false;
}

std::vector<bool> BuildNotchVertexMask(
    const std::vector<pcl::PointXYZ>& polygon,
    const std::vector<NotchFeature>& notches)
{
    std::vector<bool> mask(polygon.size(), false);
    const std::size_t n = polygon.size();
    for (const auto& notch : notches) {
        std::size_t s = (notch.entry + n - 1) % n;   // include the incoming wall edge start
        while (true) {
            mask[s % n] = true;
            if (s % n == (notch.exit + 1) % n) break;
            s = (s + 1) % n;
        }
    }
    return mask;
}

// ===== Structure-aware hypothesis repair =====
// A spike vertex V with hypothesis neighbours P,N defines a local frame:
// u along the chord P->N, v toward V. Original outline points inside the
// window should show a 3-segment rectangular structure (entry wall, bottom,
// exit wall) when the spike actually represents a real notch/protrusion.
struct SpikeStructure {
    std::size_t spikeIndex = 0;
    bool isNotch = true;
    double depth = 0.0;         // structure depth along v
    double width = 0.0;         // opening width along u
    int rawBottom = 0;          // original points on the bottom band
    int rawWallLeft = 0;        // original points on each wall
    int rawWallRight = 0;
    int cloudSupport = 0;       // mesh support points near the 3 segments
    pcl::PointXYZ c1, c2, c3, c4;   // corners, v=0 entry line to v=depth
};

double MeanDistanceToPolygonBoundary(
    const std::vector<pcl::PointXYZ>& points,
    const std::vector<pcl::PointXYZ>& polygon,
    const std::vector<std::size_t>& indices)
{
    if (indices.empty() || polygon.size() < 3) return 0.0;
    double sum = 0.0;
    for (std::size_t idx : indices) {
        const auto& p = points[idx];
        double best = std::numeric_limits<double>::max();
        for (std::size_t e = 0; e < polygon.size(); ++e) {
            double t = 0.0;
            best = std::min(best, distancePointToSegmentWithProjection(
                p, polygon[e], polygon[(e + 1) % polygon.size()], t));
        }
        sum += best;
    }
    return sum / static_cast<double>(indices.size());
}

bool DetectSpikeStructure(
    const std::vector<pcl::PointXYZ>& polygon,
    std::size_t spikeIndex,
    const std::vector<pcl::PointXYZ>& original,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support,
    SpikeStructure& out)
{
    const std::size_t n = polygon.size();
    if (n < 4 || spikeIndex >= n || original.size() < 6) return false;
    const auto& P = polygon[(spikeIndex + n - 1) % n];
    const auto& V = polygon[spikeIndex];
    const auto& N = polygon[(spikeIndex + 1) % n];
    const double chordX = N.x - P.x;
    const double chordY = N.y - P.y;
    const double chordLen = std::hypot(chordX, chordY);
    if (chordLen < kSpikeChordMin || chordLen > kSpikeChordMax) return false;

    const double ux = chordX / chordLen;
    const double uy = chordY / chordLen;
    // 尖刺高度：V 到弦线的垂距。
    const double crossPN = chordX * (V.y - P.y) - chordY * (V.x - P.x);
    const double height = std::abs(crossPN) / chordLen;
    if (height < kSpikeMinDepth || height > kSpikeMaxDepth) return false;

    // 朝向：多边形多数转向符号；V 在外侧 => 凹口。
    double turnSum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = polygon[i];
        const auto& b = polygon[(i + 1) % n];
        const auto& c = polygon[(i + 2) % n];
        turnSum += static_cast<double>(b.x - a.x) * (c.y - b.y) -
                   static_cast<double>(b.y - a.y) * (c.x - b.x);
    }
    const double majority = turnSum >= 0.0 ? 1.0 : -1.0;
    const double sideSign = crossPN >= 0.0 ? 1.0 : -1.0;
    const bool isNotch = sideSign * majority < 0.0;

    // 局部坐标系：v 轴指向 V。
    const double vx = -uy * sideSign;
    const double vy = ux * sideSign;

    // 弦线附近的原始点窗口。
    std::vector<std::size_t> window;
    for (std::size_t i = 0; i < original.size(); ++i) {
        const double dx = original[i].x - P.x;
        const double dy = original[i].y - P.y;
        const double u = dx * ux + dy * uy;
        const double v = dx * vx + dy * vy;
        if (u < -1.2 || u > chordLen + 1.2) continue;
        if (v < -0.6 || v > height + 2.0) continue;
        window.push_back(i);
    }
    if (window.size() < 6) return false;

    auto uvOf = [&](std::size_t i, double& u, double& v) {
        const double dx = original[i].x - P.x;
        const double dy = original[i].y - P.y;
        u = dx * ux + dy * uy;
        v = dx * vx + dy * vy;
    };
    double vMax = 0.0;
    for (std::size_t i : window) {
        double u = 0.0, v = 0.0;
        uvOf(i, u, v);
        vMax = std::max(vMax, v);
    }
    if (vMax < 0.75 * height) return false;   // spike not backed by raw data

    // 底边带：最深的点群确定开口 [uLo, uHi]。
    std::vector<double> bottomU;
    for (std::size_t i : window) {
        double u = 0.0, v = 0.0;
        uvOf(i, u, v);
        if (v >= vMax - kRepairBottomBand) bottomU.push_back(u);
    }
    if (bottomU.size() < 3) return false;
    std::sort(bottomU.begin(), bottomU.end());
    const double uLo = bottomU[bottomU.size() / 10];
    const double uHi = bottomU[bottomU.size() - 1 - bottomU.size() / 10];
    const double width = uHi - uLo;
    if (width < kRepairMinWidth || width > kRepairMaxWidth) return false;
    if (uLo < -0.8 || uHi > chordLen + 0.8) return false;

    // 墙壁：贴着 u 边界、沿 v 有跨度的点。
    int rawLeft = 0, rawRight = 0, rawBottom = 0;
    double leftVMin = 1e9, leftVMax = -1e9;
    double rightVMin = 1e9, rightVMax = -1e9;
    for (std::size_t i : window) {
        double u = 0.0, v = 0.0;
        uvOf(i, u, v);
        if (v >= vMax - kRepairBottomBand) ++rawBottom;
        if (u <= uLo + kRepairWallBand && v >= 0.35 * vMax) {
            ++rawLeft;
            leftVMin = std::min(leftVMin, v);
            leftVMax = std::max(leftVMax, v);
        }
        if (u >= uHi - kRepairWallBand && v >= 0.35 * vMax) {
            ++rawRight;
            rightVMin = std::min(rightVMin, v);
            rightVMax = std::max(rightVMax, v);
        }
    }
    if (rawLeft < kRepairMinRawWallPoints || rawRight < kRepairMinRawWallPoints) return false;
    if (leftVMax - leftVMin < kRepairMinWallSpan || rightVMax - rightVMin < kRepairMinWallSpan) {
        return false;
    }

    // 角点换回世界坐标。
    auto toXY = [&](double u, double v, pcl::PointXYZ& p) {
        p.x = static_cast<float>(P.x + u * ux + v * vx);
        p.y = static_cast<float>(P.y + u * uy + v * vy);
        p.z = V.z;
    };
    SpikeStructure s;
    s.spikeIndex = spikeIndex;
    s.isNotch = isNotch;
    s.depth = vMax;
    s.width = width;
    s.rawBottom = rawBottom;
    s.rawWallLeft = rawLeft;
    s.rawWallRight = rawRight;
    toXY(uLo, 0.0, s.c1);
    toXY(uLo, vMax, s.c2);
    toXY(uHi, vMax, s.c3);
    toXY(uHi, 0.0, s.c4);

    // 三条新边附近的网格支撑(该处点云被
    // 遮挡时走下方的弱分支)。
    if (support && !support->empty()) {
        const pcl::PointXYZ seg[3][2] = {{s.c1, s.c2}, {s.c2, s.c3}, {s.c3, s.c4}};
        for (const auto& p : support->points) {
            double best = std::numeric_limits<double>::max();
            for (int k = 0; k < 3; ++k) {
                double t = 0.0;
                best = std::min(best, distancePointToSegmentWithProjection(
                    p, seg[k][0], seg[k][1], t));
            }
            if (best <= kRepairSupportRadius) ++s.cloudSupport;
        }
    }
    const bool strongSupport = s.cloudSupport >= kRepairMinSupportPoints;
    const bool weakButConvincing =
        s.cloudSupport >= 2 &&
        s.rawWallLeft >= 2 * kRepairMinRawWallPoints &&
        s.rawWallRight >= 2 * kRepairMinRawWallPoints &&
        s.depth >= 1.5 * kSpikeMinDepth;
    if (!strongSupport && !weakButConvincing) return false;

    out = s;
    return true;
}

// 全局修复统计，每次运行打印一次。
long long g_repairFeaturesExamined = 0;
long long g_repairFeaturesWithCandidates = 0;
long long g_repairFeaturesRepaired = 0;
long long g_repairCandidatesDetected = 0;
long long g_repairAccepted = 0;
long long g_repairVerticesAdded = 0;
long long g_repairMaxVerticesAdded = 0;

// 贪心修复：每轮选出净收益最高、
// 数据/支撑收益超过新增顶点代价的结构，
// 重新检测，最多插入 kRepairMaxPerPolygon 处。
void RepairHypothesisStructures(
    std::vector<pcl::PointXYZ>& polygon,
    const std::vector<pcl::PointXYZ>& original,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support,
    long long sourceFid)
{
    ++g_repairFeaturesExamined;
    if (polygon.size() < 4) return;
    const double baseArea = polygonArea2D(polygon);

    bool anyCandidate = false;
    int accepted = 0;
    int verticesAdded = 0;

    for (int round = 0; round < kRepairMaxPerPolygon; ++round) {
        const std::size_t n = polygon.size();
        double bestScore = 0.0;   // must beat zero (cost included)
        std::size_t bestSpike = n;
        SpikeStructure bestStruct;
        std::vector<pcl::PointXYZ> bestCandidate;

        for (std::size_t i = 0; i < n; ++i) {
            SpikeStructure s;
            if (!DetectSpikeStructure(polygon, i, original, support, s)) continue;
            anyCandidate = true;
            ++g_repairCandidatesDetected;

            // 构建修复后的多边形：P ... V ... N -> P c1 c2 c3 c4 N。
            std::vector<pcl::PointXYZ> candidate;
            candidate.reserve(n + 3);
            for (std::size_t k = 0; k < n; ++k) {
                if (k == s.spikeIndex) {
                    candidate.push_back(s.c1);
                    candidate.push_back(s.c2);
                    candidate.push_back(s.c3);
                    candidate.push_back(s.c4);
                } else {
                    candidate.push_back(polygon[k]);
                }
            }
            removeDuplicatePoints2D(candidate, 1e-4f);
            if (candidate.size() < 4 || !isSimplePolygon2D(candidate)) continue;
            if (std::abs(polygonArea2D(candidate) - baseArea) >
                    kRepairMaxAreaChange * std::max(baseArea, 1.0)) {
                continue;
            }
            double minEdge = std::numeric_limits<double>::max();
            for (std::size_t e = 0; e < candidate.size(); ++e) {
                minEdge = std::min(minEdge, std::hypot(
                    static_cast<double>(candidate[(e + 1) % candidate.size()].x - candidate[e].x),
                    static_cast<double>(candidate[(e + 1) % candidate.size()].y - candidate[e].y)));
            }
            if (minEdge < kRepairMinNewEdge) continue;

            // 收益：尖刺附近的原始点和支撑点
            // 到边界的距离下降须超过顶点代价。
            double dataGain = 0.0;
            double supportGain = 0.0;
            std::vector<std::size_t> rawIdx;
            for (std::size_t q = 0; q < original.size(); ++q) {
                if (std::hypot(original[q].x - polygon[s.spikeIndex].x,
                               original[q].y - polygon[s.spikeIndex].y) <=
                    s.depth + 1.5) {
                    rawIdx.push_back(q);
                }
            }
            std::vector<pcl::PointXYZ> rawStorage;
            for (std::size_t q : rawIdx) rawStorage.push_back(original[q]);
            std::vector<std::size_t> rawAll(rawIdx.size());
            for (std::size_t q = 0; q < rawIdx.size(); ++q) rawAll[q] = q;
            if (!rawAll.empty()) {
                dataGain = MeanDistanceToPolygonBoundary(original, polygon, rawIdx) -
                          MeanDistanceToPolygonBoundary(rawStorage, candidate, rawAll);
            }
            if (support && !support->empty()) {
                const auto& V0 = polygon[s.spikeIndex];
                std::vector<std::size_t> supIdx;
                for (std::size_t q = 0; q < support->size(); ++q) {
                    if (std::hypot(support->points[q].x - V0.x,
                                   support->points[q].y - V0.y) <= s.depth + 1.5) {
                        supIdx.push_back(q);
                    }
                }
                if (!supIdx.empty()) {
                    // support->points uses an aligned allocator; copy into a
                    // plain vector for the shared distance helper.
                    std::vector<pcl::PointXYZ> supStorage;
                    for (std::size_t q : supIdx) supStorage.push_back(support->points[q]);
                    std::vector<std::size_t> supAll(supIdx.size());
                    for (std::size_t q = 0; q < supIdx.size(); ++q) supAll[q] = q;
                    supportGain = MeanDistanceToPolygonBoundary(supStorage, polygon, supAll) -
                                  MeanDistanceToPolygonBoundary(supStorage, candidate, supAll);
                }
            }

            const double cost = 3.0 * kRepairPerVertexCost;
            const double score = kRepairDataWeight * dataGain +
                                 kRepairSupportWeight * supportGain - cost;
            std::cerr << "[HypothesisRepair] source_fid=" << sourceFid
                      << " feature=" << (s.isNotch ? "notch" : "protrusion")
                      << " depth=" << s.depth << " width=" << s.width
                      << " raw_bottom=" << s.rawBottom
                      << " raw_wall=" << s.rawWallLeft << "/" << s.rawWallRight
                      << " cloud_support=" << s.cloudSupport
                      << " data_gain=" << dataGain
                      << " support_gain=" << supportGain
                      << " complexity_cost=" << cost
                      << " score=" << score << std::endl;

            if (score > bestScore) {
                bestScore = score;
                bestSpike = s.spikeIndex;
                bestStruct = s;
                bestCandidate = std::move(candidate);
            }
        }

        if (bestSpike >= polygon.size() || bestCandidate.empty()) break;

        std::cerr << "[HypothesisRepair] source_fid=" << sourceFid
                  << " accepted=1 round=" << round
                  << " vertices=" << polygon.size() << "->" << bestCandidate.size()
                  << " score=" << bestScore << std::endl;
        verticesAdded += static_cast<int>(bestCandidate.size() - polygon.size());
        polygon = std::move(bestCandidate);
        ++accepted;
        ++g_repairAccepted;
    }

    if (anyCandidate) ++g_repairFeaturesWithCandidates;
    if (accepted > 0) {
        ++g_repairFeaturesRepaired;
        g_repairVerticesAdded += verticesAdded;
        g_repairMaxVerticesAdded = std::max(g_repairMaxVerticesAdded,
                                            static_cast<long long>(verticesAdded));
    }
}

void PrintHypothesisRepairStats()
{
    std::cerr << "[HypothesisRepair] summary examined=" << g_repairFeaturesExamined
              << ", with_candidates=" << g_repairFeaturesWithCandidates
              << ", repaired_features=" << g_repairFeaturesRepaired
              << ", accepted_repairs=" << g_repairAccepted
              << ", vertices_added_total=" << g_repairVerticesAdded
              << ", avg_per_repaired="
              << (g_repairFeaturesRepaired > 0
                      ? static_cast<double>(g_repairVerticesAdded) /
                            static_cast<double>(g_repairFeaturesRepaired)
                      : 0.0)
              << ", max_per_feature=" << g_repairMaxVerticesAdded << std::endl;
}

bool forceOrthogonalPolygonToAngle(
    const std::vector<pcl::PointXYZ>& input,
    double main_angle,
    std::vector<pcl::PointXYZ>& result,
    const std::vector<bool>* protectedVertices = nullptr)
{
    result.clear();
    if (input.size() < 4) return false;

    std::vector<pcl::PointXYZ> working = input;
    std::vector<char> protection;
    if (protectedVertices && protectedVertices->size() == input.size()) {
        protection.assign(protectedVertices->begin(), protectedVertices->end());
    }
    auto eraseWorking = [&working, &protection](std::size_t index) {
        working.erase(working.begin() + index);
        if (!protection.empty()) protection.erase(protection.begin() + index);
    };
    auto axisOfEdge = [main_angle](const pcl::PointXYZ& a, const pcl::PointXYZ& b) {
        const double angle = std::atan2(b.y - a.y, b.x - a.x);
        return undirectedAngleDifference(angle, main_angle) <=
            undirectedAngleDifference(angle, main_angle + M_PI / 2.0) ? 0 : 1;
    };

    // A valid orthogonal ring alternates its two axes. Remove vertices whose
    // incoming and outgoing edges would snap to the same axis before solving
    // line intersections; otherwise the two lines are parallel and no corner
    // exists. Odd vertex counts are resolved by the same cyclic pass.
    // 受保护的(证据支撑缺口)顶点绝不在此时删除；
    // 若保护顶点无法满足交替轴模式，下方重建
    // 失败，调用方换下一个候选源。
    bool changed = true;
    while (changed && working.size() >= 4) {
        changed = false;
        for (size_t i = 0; i < working.size(); ++i) {
            if (!protection.empty() && protection[i]) continue;
            const size_t previous = (i + working.size() - 1) % working.size();
            const size_t next = (i + 1) % working.size();
            if (axisOfEdge(working[previous], working[i]) ==
                axisOfEdge(working[i], working[next])) {
                eraseWorking(i);
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

void RemoveClosingDuplicate(std::vector<pcl::PointXYZ>& points);

// ===== Mask-only 局部圆弧/椭圆弧检测与恢复(独立通道) =====
// 不复用 OSGB 的 detectPreservedArcs(Bezier); Mask-only 只允许 CircleArc。
// 检测在平滑环上定区间, rawRing 提供拟合支撑(raw==smooth 时用环采样+严格条件)。
// 恢复在直线候选通过质量检查之后执行。

static void CircleDistStats(const std::vector<pcl::PointXYZ>& pts,
    std::size_t start, std::size_t count, double cx, double cy, double radius,
    double& rmse, double& q90)
{
    std::vector<double> dists;
    dists.reserve(count);
    double err2 = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
        const auto& p = pts[(start + k) % pts.size()];
        const double d = std::abs(std::hypot(p.x - cx, p.y - cy) - radius);
        dists.push_back(d);
        err2 += d * d;
    }
    rmse = std::sqrt(err2 / count);
    std::sort(dists.begin(), dists.end());
    q90 = dists[static_cast<std::size_t>(count * 0.9)];
}

static double LineFitRmseOf(const std::vector<pcl::PointXYZ>& pts,
    std::size_t start, std::size_t count)
{
    if (count < 2) return 0.0;
    double mx = 0.0, my = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
        const auto& p = pts[(start + k) % pts.size()];
        mx += p.x; my += p.y;
    }
    mx /= count; my /= count;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
        const auto& p = pts[(start + k) % pts.size()];
        const double dx = p.x - mx, dy = p.y - my;
        sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
    }
    const double th = 0.5 * std::atan2(2 * sxy, sxx - syy);
    const double dX = std::cos(th), dY = std::sin(th);
    double e2 = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
        const auto& p = pts[(start + k) % pts.size()];
        const double dx = p.x - mx, dy = p.y - my;
        const double perp = -dx * dY + dy * dX;
        e2 += perp * perp;
    }
    return std::sqrt(e2 / count);
}

// v2: 开放曲线转角(仅内部顶点,无闭合)
struct OpenTurnStats {
    double signedTurnDeg = 0.0;
    double absoluteTurnDeg = 0.0;
    double consistency = 0.0;
    int reversalCount = 0;
    int totalTurns = 0;
};
static OpenTurnStats OpenCurveTurnStats(const std::vector<pcl::PointXYZ>& pts)
{
    OpenTurnStats s;
    const int n = (int)pts.size();
    if (n < 3) return s;
    double prevCross = 0.0;
    for (int i = 1; i < n - 1; ++i) {
        const double v1x = pts[i].x - pts[i-1].x;
        const double v1y = pts[i].y - pts[i-1].y;
        const double v2x = pts[i+1].x - pts[i].x;
        const double v2y = pts[i+1].y - pts[i].y;
        const double cross = v1x * v2y - v1y * v2x;
        const double ang = std::atan2(cross, v1x * v2x + v1y * v2y);
        s.signedTurnDeg += ang * 180.0 / M_PI;
        s.absoluteTurnDeg += std::abs(ang) * 180.0 / M_PI;
        ++s.totalTurns;
        if (i > 1 && cross * prevCross < 0) ++s.reversalCount;
        prevCross = cross;
    }
    s.consistency = s.absoluteTurnDeg > 1e-9
        ? std::abs(s.signedTurnDeg) / s.absoluteTurnDeg : 0.0;
    return s;
}

// v2: smooth→raw全局单调映射
struct SmoothRawMap {
    std::vector<std::size_t> rawUnwrapped;
    std::size_t rawStart = 0;
    bool rawReversed = false;
    bool valid = false;
    std::vector<double> rawCumArc;
    double totalRawArc = 0.0;
    std::size_t rawN = 0;
};

static SmoothRawMap BuildSmoothRawMap(
    const std::vector<pcl::PointXYZ>& smoothRing,
    const std::vector<pcl::PointXYZ>& rawRing)
{
    SmoothRawMap m;
    const std::size_t sn = smoothRing.size();
    const std::size_t rn = rawRing.size();
    if (sn < 4 || rn < 4) return m;
    m.rawN = rn;
    auto ringArea = [](const std::vector<pcl::PointXYZ>& r) {
        double a = 0;
        for (std::size_t i = 0; i < r.size(); ++i) {
            const std::size_t j = (i + 1) % r.size();
            a += r[i].x * r[j].y - r[j].x * r[i].y;
        }
        return a;
    };
    m.rawReversed = ((ringArea(smoothRing) > 0) != (ringArea(rawRing) > 0));
    auto nearestRaw = [&](const pcl::PointXYZ& q) {
        std::size_t best = 0; double bd = 1e18;
        for (std::size_t i = 0; i < rn; ++i) {
            const double d = std::hypot(rawRing[i].x - q.x, rawRing[i].y - q.y);
            if (d < bd) { bd = d; best = i; }
        }
        return best;
    };
    m.rawStart = nearestRaw(smoothRing[0]);
    m.rawCumArc.resize(rn);
    m.rawCumArc[0] = 0.0;
    double cum = 0.0;
    for (std::size_t k = 0; k < rn; ++k) {
        const std::size_t cur = m.rawReversed ? (m.rawStart + rn - k) % rn : (m.rawStart + k) % rn;
        const std::size_t nxt = m.rawReversed ? (m.rawStart + rn - k - 1) % rn : (m.rawStart + k + 1) % rn;
        cum += std::hypot(rawRing[nxt].x - rawRing[cur].x, rawRing[nxt].y - rawRing[cur].y);
        if (k + 1 < rn) m.rawCumArc[k + 1] = cum;
    }
    m.totalRawArc = cum;
    m.rawUnwrapped.resize(sn);
    std::size_t rawPtr = 0;
    for (std::size_t si = 0; si < sn; ++si) {
        std::size_t bestOff = 0;
        double bestD = 1e18;
        for (std::size_t off = 0; off < rn; ++off) {
            const std::size_t ri = (rawPtr + off) % rn;
            const double d = std::hypot(rawRing[ri].x - smoothRing[si].x,
                                        rawRing[ri].y - smoothRing[si].y);
            if (d < bestD) { bestD = d; bestOff = off; }
            if (off > 10 && d > bestD * 3) break;
        }
        rawPtr = (rawPtr + bestOff) % rn;
        m.rawUnwrapped[si] = rawPtr;
    }
    m.valid = true;
    return m;
}

struct RawInterval {
    std::size_t rawStartIdx = 0;
    std::size_t rawEndIdx = 0;
    int pointCount = 0;
    double arcLength = 0.0;
    bool valid = false;
};

static RawInterval ExtractRawInterval(
    const SmoothRawMap& map,
    std::size_t smoothStart, std::size_t smoothEnd)
{
    RawInterval ri;
    if (!map.valid || map.rawUnwrapped.empty()) return ri;
    const std::size_t rn = map.rawN;
    const std::size_t uS = map.rawUnwrapped[smoothStart];
    const std::size_t uE = map.rawUnwrapped[smoothEnd];
    std::size_t span = (uE >= uS) ? (uE - uS) : (rn - uS + uE);
    if (span >= rn) return ri;
    ri.rawStartIdx = uS;
    ri.rawEndIdx = uE;
    ri.pointCount = (int)span + 1;
    if (uE < uS) {
        ri.arcLength = map.totalRawArc - map.rawCumArc[uS] + map.rawCumArc[uE];
    } else {
        ri.arcLength = map.rawCumArc[uE] - map.rawCumArc[uS];
    }
    ri.valid = ri.pointCount >= 2;
    return ri;
}

static std::vector<pcl::PointXYZ> CollectRawPoints(
    const std::vector<pcl::PointXYZ>& rawRing,
    const RawInterval& ri)
{
    std::vector<pcl::PointXYZ> pts;
    if (!ri.valid) return pts;
    const std::size_t rn = rawRing.size();
    pts.reserve((std::size_t)ri.pointCount);
    for (std::size_t k = 0; k < (std::size_t)ri.pointCount; ++k) {
        pts.push_back(rawRing[(ri.rawStartIdx + k) % rn]);
    }
    return pts;
}

static std::vector<pcl::PointXYZ> WeightedSample(
    const std::vector<pcl::PointXYZ>& pts, int targetCount)
{
    const int n = (int)pts.size();
    if (n < 2 || targetCount < 2 || n >= targetCount) return pts;
    std::vector<double> cum(n, 0.0);
    for (int i = 1; i < n; ++i)
        cum[i] = cum[i-1] + std::hypot(pts[i].x - pts[i-1].x, pts[i].y - pts[i-1].y);
    const double total = cum[n-1];
    if (total < 1e-9) return pts;
    std::vector<pcl::PointXYZ> out;
    out.reserve(targetCount);
    for (int k = 0; k < targetCount; ++k) {
        const double target = (double)k / (targetCount - 1) * total;
        int lo = 0, hi = n - 1;
        while (lo < hi - 1) {
            const int mid = (lo + hi) / 2;
            if (cum[mid] <= target) lo = mid; else hi = mid;
        }
        const double t = (target - cum[lo]) / std::max(cum[hi] - cum[lo], 1e-12);
        pcl::PointXYZ p;
        p.x = (float)(pts[lo].x + t * (pts[hi].x - pts[lo].x));
        p.y = (float)(pts[lo].y + t * (pts[hi].y - pts[lo].y));
        p.z = pts[lo].z;
        out.push_back(p);
    }
    return out;
}

// v2: 圆弧评估(转角用smoothPts, 拟合用fitPts)
static std::string EvalCircleV2(
    const std::vector<pcl::PointXYZ>& fitPts,
    const std::vector<pcl::PointXYZ>& smoothPts,
    int sourceVertexCount,
    double pixelSize,
    MaskConicArc& out)
{
    const int n = (int)fitPts.size();
    if (n < 5) return "too_few_points";
    double mx = 0, my = 0;
    for (const auto& p : fitPts) { mx += p.x; my += p.y; }
    mx /= n; my /= n;
    std::vector<pcl::PointXYZ> c(n);
    for (int i = 0; i < n; ++i) {
        c[i].x = (float)(fitPts[i].x - mx);
        c[i].y = (float)(fitPts[i].y - my);
        c[i].z = fitPts[i].z;
    }
    double cx = 0, cy = 0, r = 0, rmse = 0;
    if (!fitCircle2D(c, 0, n, cx, cy, r, rmse)) return "circle_fit_failed";
    double q90 = 0;
    CircleDistStats(c, 0, n, cx, cy, r, rmse, q90);
    const double lineRmse = LineFitRmseOf(c, 0, n);
    const double chord = std::hypot(c[n-1].x - c[0].x, c[n-1].y - c[0].y);
    double arcLen = 0;
    for (int i = 0; i < n - 1; ++i)
        arcLen += std::hypot(c[i+1].x - c[i].x, c[i+1].y - c[i].y);
    // 转角: 用smoothPts(去中心化+采样)
    OpenTurnStats ts;
    if (smoothPts.size() >= 3) {
        double smx = 0, smy = 0;
        for (const auto& p : smoothPts) { smx += p.x; smy += p.y; }
        smx /= smoothPts.size(); smy /= smoothPts.size();
        std::vector<pcl::PointXYZ> sc(smoothPts.size());
        for (std::size_t i = 0; i < smoothPts.size(); ++i) {
            sc[i].x = (float)(smoothPts[i].x - smx);
            sc[i].y = (float)(smoothPts[i].y - smy);
            sc[i].z = smoothPts[i].z;
        }
        auto sampled = WeightedSample(sc, std::max(12, (int)sc.size()));
        ts = OpenCurveTurnStats(sampled);
    } else {
        ts = OpenCurveTurnStats(c);
    }
    double sa = 0;
    const double sweep = unwrapArcSweep(c, 0, n, cx, cy, sa);
    const double sweepDeg = std::abs(sweep) * 180.0 / M_PI;

    // 门限
    if (arcLen < 10.0) return "too_short";
    if (chord < 6.0) return "chord_too_short";
    if (ts.absoluteTurnDeg < 20.0) return "turn_too_small";
    if (ts.consistency < 0.80) return "turn_inconsistent";
    if (ts.totalTurns > 0 && (double)ts.reversalCount / ts.totalTurns > 0.10)
        return "reversal_too_many";
    if (rmse > std::max(2.0 * pixelSize, 0.60)) return "circle_rmse";
    if (q90 > std::max(3.0 * pixelSize, 1.00)) return "circle_q90";
    if (lineRmse / std::max(rmse, 0.01) < 2.0) return "line_explains_better";
    if (sweepDeg < 20.0 || sweepDeg > 240.0) return "sweep_out_of_range";

    const double improve = std::max(0.0, 1.0 - rmse / std::max(lineRmse, 0.01));
    const double fitTol = std::max(0.5 * pixelSize, 0.1);
    const double fitQ = std::exp(-(rmse / fitTol) * (rmse / fitTol));
    out.type = MaskCurveType::CircleArc;
    out.cx = cx + mx; out.cy = cy + my; out.radius = r;
    out.startAngle = sa; out.sweepAngle = sweep;
    out.rmse = rmse; out.q90 = q90; out.lineRmse = lineRmse;
    out.arcLength = arcLen; out.chordLength = chord;
    out.sweepDeg = sweepDeg;
    out.supportCount = sourceVertexCount;
    out.score = improve * fitQ * ts.consistency * std::min(1.0, arcLen / 15.0);
    if (out.score <= 0.0) return "zero_score";
    return "";
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

// ===== Direction chain extraction & multi-peak detection =====
// Split-and-merge continuous edge chain extraction on a closed ring,
// followed by weighted KDE multi-peak direction detection in [0°, 90°).
// Chains replace individual edges as direction evidence: a staircase edge
// is one noisy chain, not a "diagonal direction".  Curves decompose into
// short high-RMSE chains that get downweighted, not misread as multi-dir.

struct DirectionChain {
    // 索引：基于原环(旋转前)的索引，闭环安全
    std::size_t startIndexOriginal = 0;
    std::size_t endIndexOriginal = 0;
    // 端点坐标(直接存储，不通过索引间接访问)
    pcl::PointXYZ startPoint;
    pcl::PointXYZ endPoint;
    double length = 0.0;             // meters along the chain
    double angleRad = 0.0;           // principal direction, folded to [0, π/2)
    double rmse = 0.0;               // fit RMSE (m)
    double maxDeviation = 0.0;       // max point-to-line distance (m)
    double centerX = 0.0;            // chain midpoint (for spatial independence)
    double centerY = 0.0;
    double weight = 0.0;             // voting weight for direction estimation
    // 拟合直线参数: normalX*x + normalY*y = lineOffset
    double normalX = 0.0;
    double normalY = 0.0;
    double lineOffset = 0.0;
    // 拓扑 vs 方向：短链保留在拓扑中，方向投票时过滤
    bool isShort = false;            // length < kDirectionMinChainLength
};

struct DirectionPeak {
    double angleRad = 0.0;           // folded to [0, π/2)
    double ratio = 0.0;              // weightedLength / totalWeight
    double weightedLength = 0.0;
    int chainCount = 0;
    double meanRmse = 0.0;
    bool strong = false;
};

struct DirectionEvidence2D {
    bool valid = false;
    DirectionPeak primary;
    DirectionPeak secondary;
    int totalChains = 0;
    double totalWeightedLength = 0.0;
    double fitTolerance = 0.0;        // the tolerance used for chain fitting
    std::vector<DirectionChain> chains;
    std::vector<DirectionPeak> peaks; // all peaks (primary + secondary + minor)
    // evidence source: 0=mask chains, 1=OSGB support, 2=fused
    int source = 0;
};

// Fit a line to ring[start..end] (inclusive, no wraparound in this helper).
// Returns direction, normal, centroid-based line offset, RMSE and max deviation.
bool FitChainLine(
    const std::vector<pcl::PointXYZ>& ring,
    std::size_t start, std::size_t end,
    double& dirX, double& dirY,
    double& rmse, double& maxDev,
    double* outCentroidX = nullptr,
    double* outCentroidY = nullptr,
    double* outNormalX = nullptr,
    double* outNormalY = nullptr,
    double* outLineOffset = nullptr)
{
    if (end <= start) return false;
    const std::size_t n = end - start + 1;
    double mx = 0.0, my = 0.0;
    for (std::size_t i = start; i <= end; ++i) {
        mx += ring[i].x; my += ring[i].y;
    }
    mx /= n; my /= n;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (std::size_t i = start; i <= end; ++i) {
        const double dx = ring[i].x - mx;
        const double dy = ring[i].y - my;
        sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
    }
    double theta = 0.5 * std::atan2(2.0 * sxy, sxx - syy);
    dirX = std::cos(theta);
    dirY = std::sin(theta);
    // normal = perpendicular to direction
    const double nx = -dirY;
    const double ny = dirX;
    // line offset from centroid (not from start point!)
    const double d = nx * mx + ny * my;
    rmse = 0.0; maxDev = 0.0;
    for (std::size_t i = start; i <= end; ++i) {
        const double dx = ring[i].x - mx;
        const double dy = ring[i].y - my;
        const double perp = -dx * dirY + dy * dirX;
        rmse += perp * perp;
        maxDev = std::max(maxDev, std::abs(perp));
    }
    rmse = std::sqrt(rmse / n);
    if (outCentroidX) *outCentroidX = mx;
    if (outCentroidY) *outCentroidY = my;
    if (outNormalX) *outNormalX = nx;
    if (outNormalY) *outNormalY = ny;
    if (outLineOffset) *outLineOffset = d;
    return std::isfinite(rmse) && std::isfinite(maxDev);
}

// Recursive split: if the chain from start to end has max deviation > tolerance,
// split at the max-deviation vertex and recurse both halves.
void SplitChainRecursive(
    const std::vector<pcl::PointXYZ>& ring,
    std::size_t start, std::size_t end,
    double tolerance,
    std::vector<std::pair<std::size_t, std::size_t>>& segments)
{
    if (end <= start + 1) {
        segments.push_back({start, end});
        return;
    }
    double dirX, dirY, rmse, maxDev;
    if (!FitChainLine(ring, start, end, dirX, dirY, rmse, maxDev)) {
        segments.push_back({start, end});
        return;
    }
    if (maxDev <= tolerance) {
        segments.push_back({start, end});
        return;
    }
    // find max deviation vertex — only interior vertices are valid split
    // points. The old version scanned endpoints too and gave up when the
    // argmax landed on one, emitting badly curved spans as single chains
    // (fatal for chain-line intersection topology).
    std::size_t maxIdx = start;      // sentinel: no interior vertex found
    double bestDev = 0.0;
    double mx = 0.0, my = 0.0;
    const std::size_t n = end - start + 1;
    for (std::size_t i = start; i <= end; ++i) { mx += ring[i].x; my += ring[i].y; }
    mx /= n; my /= n;
    for (std::size_t i = start + 1; i < end; ++i) {
        const double dx = ring[i].x - mx;
        const double dy = ring[i].y - my;
        const double perp = -dx * dirY + dy * dirX;
        if (std::abs(perp) > bestDev) {
            bestDev = std::abs(perp);
            maxIdx = i;
        }
    }
    if (maxIdx == start) {
        // 无内部候选分裂点(共线或仅端点偏离)，整段作为叶子
        segments.push_back({start, end});
        return;
    }
    SplitChainRecursive(ring, start, maxIdx, tolerance, segments);
    SplitChainRecursive(ring, maxIdx, end, tolerance, segments);
}

// Merge adjacent segments whose directions are close and merged fit still
// satisfies the tolerance.
void MergeAdjacentSegments(
    const std::vector<pcl::PointXYZ>& ring,
    std::vector<std::pair<std::size_t, std::size_t>>& segments,
    double tolerance,
    double maxAngleDeg)
{
    if (segments.size() < 2) return;
    bool merged = true;
    while (merged) {
        merged = false;
        for (std::size_t i = 0; i < segments.size(); ++i) {
            std::size_t next = (i + 1) % segments.size();
            if (next == 0 && segments.size() == 1) break;
            // check angular proximity
            double dx1, dy1, rmse1, dev1, dx2, dy2, rmse2, dev2;
            if (!FitChainLine(ring, segments[i].first, segments[i].second, dx1, dy1, rmse1, dev1)) continue;
            if (!FitChainLine(ring, segments[next].first, segments[next].second, dx2, dy2, rmse2, dev2)) continue;
            double dot = std::abs(dx1 * dx2 + dy1 * dy2);
            if (dot < std::cos(maxAngleDeg * M_PI / 180.0)) continue;
            // try merged fit
            std::size_t s = segments[i].first;
            std::size_t e = segments[next].second;
            double mDx, mDy, mRmse, mDev;
            // handle wraparound: if next wraps to index 0, merged goes from s to ring.size()-1 + next.second
            // for simplicity, only merge non-wrapping pairs in this pass
            if (next == 0) continue;  // skip wraparound in this simple version
            if (!FitChainLine(ring, s, e, mDx, mDy, mRmse, mDev)) continue;
            if (mDev > tolerance) continue;
            segments[i].second = e;
            segments.erase(segments.begin() + next);
            merged = true;
            break;
        }
    }
}

// Extract continuous edge chains from a closed ring via split-and-merge.
std::vector<DirectionChain> ExtractDirectionChains(
    const std::vector<pcl::PointXYZ>& ring,
    double fitTolerance)
{
    std::vector<DirectionChain> chains;
    if (ring.size() < 4) return chains;

    // ---- 闭环旋转(修复: 用 ring.size() 取模, 不用 extended.size()) ----
    // 从质心最远点起始，避免在真实边中间分裂
    double cx = 0.0, cy = 0.0;
    for (const auto& p : ring) { cx += p.x; cy += p.y; }
    cx /= ring.size(); cy /= ring.size();
    std::size_t startIdx = 0;
    double maxDist = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const double d = std::hypot(ring[i].x - cx, ring[i].y - cy);
        if (d > maxDist) { maxDist = d; startIdx = i; }
    }

    // 正确闭环旋转: rotated[i] = ring[(startIdx+i) % ring.size()]，最后加闭合点
    std::vector<pcl::PointXYZ> rotated;
    rotated.reserve(ring.size() + 1);
    for (std::size_t i = 0; i < ring.size(); ++i) {
        rotated.push_back(ring[(startIdx + i) % ring.size()]);
    }
    rotated.push_back(rotated.front());  // 闭合: 最后一条边 = 最后顶点→首顶点

    std::vector<std::pair<std::size_t, std::size_t>> segments;
    SplitChainRecursive(rotated, 0, rotated.size() - 1, fitTolerance, segments);
    MergeAdjacentSegments(rotated, segments, fitTolerance,
        kDirectionMaxAngularResidualDeg);

    for (const auto& seg : segments) {
        if (seg.second <= seg.first) continue;
        double dirX, dirY, rmse, maxDev;
        double centroidX = 0.0, centroidY = 0.0, nx = 0.0, ny = 0.0, lineOff = 0.0;
        if (!FitChainLine(rotated, seg.first, seg.second, dirX, dirY, rmse, maxDev,
                          &centroidX, &centroidY, &nx, &ny, &lineOff)) continue;
        double length = 0.0;
        for (std::size_t i = seg.first; i < seg.second; ++i) {
            length += std::hypot(rotated[i + 1].x - rotated[i].x,
                                 rotated[i + 1].y - rotated[i].y);
        }
        // 不再按长度过滤——短链保留在拓扑中，标记 isShort 供方向投票过滤

        DirectionChain chain;
        // 原环索引 = (startIdx + rotated_index) % ring.size()
        chain.startIndexOriginal = (startIdx + seg.first) % ring.size();
        chain.endIndexOriginal = (startIdx + seg.second) % ring.size();
        // 直接存储端点坐标(不通过索引间接访问)
        chain.startPoint = rotated[seg.first];
        chain.endPoint = rotated[seg.second];
        chain.length = length;
        chain.angleRad = foldedLineAngle90(std::atan2(dirY, dirX));
        chain.rmse = rmse;
        chain.maxDeviation = maxDev;
        chain.centerX = 0.5 * (rotated[seg.first].x + rotated[seg.second].x);
        chain.centerY = 0.5 * (rotated[seg.first].y + rotated[seg.second].y);
        // 拟合直线参数(基于质心，非起点)
        chain.normalX = nx;
        chain.normalY = ny;
        chain.lineOffset = lineOff;
        chain.isShort = length < kDirectionMinChainLength;
        // voting weight: length × fit quality (短链权重小)
        const double fitScore = std::exp(
            -(rmse / std::max(fitTolerance, 0.01)) * (rmse / std::max(fitTolerance, 0.01)));
        chain.weight = length * fitScore;
        chains.push_back(chain);
    }
    return chains;
}

// Weighted KDE multi-peak detection on chain angles in [0°, 90°).
std::vector<DirectionPeak> DetectDirectionPeaks(
    const std::vector<DirectionChain>& chains,
    double fitTolerance)
{
    std::vector<DirectionPeak> peaks;
    if (chains.empty()) return peaks;

    double totalWeight = 0.0;
    for (const auto& c : chains) totalWeight += c.weight;
    if (totalWeight < 1e-9) return peaks;

    // sample the KDE at 0.5° intervals
    constexpr int kBins = 180;  // 90° / 0.5°
    std::vector<double> kde(kBins, 0.0);
    const double bandwidthRad = kDirectionKdeBandwidthDeg * M_PI / 180.0;
    for (int b = 0; b < kBins; ++b) {
        const double angle = (b + 0.5) * 0.5 * M_PI / 180.0;
        double density = 0.0;
        for (const auto& c : chains) {
            double diff = std::abs(angle - c.angleRad);
            // circular distance in [0, π/2)
            if (diff > M_PI / 4.0) diff = M_PI / 2.0 - diff;
            density += c.weight * std::exp(-0.5 * (diff / bandwidthRad) * (diff / bandwidthRad));
        }
        kde[b] = density;
    }

    // find local maxima — 折叠角空间是循环的(0°=90°), 首尾 bin 必须
    // 按循环数组比较, 否则正好落在 0°/90° 附近的峰会漏检
    std::vector<int> maxima;
    for (int b = 0; b < kBins; ++b) {
        const int prevBin = (b + kBins - 1) % kBins;
        const int nextBin = (b + 1) % kBins;
        if (kde[b] > kde[prevBin] && kde[b] >= kde[nextBin] &&
            kde[b] > totalWeight * 0.02) {
            maxima.push_back(b);
        }
    }
    if (maxima.empty()) return peaks;

    // sort by density descending
    std::sort(maxima.begin(), maxima.end(),
        [&](int a, int b) { return kde[a] > kde[b]; });

    // enforce separation
    const double minSepBins = kDirectionPeakSeparationDeg / 0.5;
    std::vector<int> selected;
    for (int m : maxima) {
        bool separated = true;
        for (int s : selected) {
            double diff = std::abs(m - s);
            if (diff > kBins / 2.0) diff = kBins - diff;
            if (diff < minSepBins) { separated = false; break; }
        }
        if (separated) selected.push_back(m);
        if (selected.size() >= 4) break;
    }

    // build peaks with chain membership
    for (int bin : selected) {
        const double peakAngle = (bin + 0.5) * 0.5 * M_PI / 180.0;
        DirectionPeak peak;
        peak.angleRad = peakAngle;
        double halfWidth = kDirectionPeakSeparationDeg * 0.5 * M_PI / 180.0;
        double rmseSum = 0.0;
        for (const auto& c : chains) {
            double diff = std::abs(c.angleRad - peakAngle);
            if (diff > M_PI / 4.0) diff = M_PI / 2.0 - diff;
            if (diff <= halfWidth) {
                peak.weightedLength += c.weight;
                peak.chainCount++;
                rmseSum += c.rmse;
            }
        }
        if (peak.chainCount > 0) peak.meanRmse = rmseSum / peak.chainCount;
        peak.ratio = peak.weightedLength / totalWeight;
        peak.strong = peak.ratio >= 0.30 && peak.chainCount >= 2;
        peaks.push_back(peak);
    }
    return peaks;
}

// Full direction evidence estimation from a polygon ring.
DirectionEvidence2D EstimateDirectionEvidence(
    const std::vector<pcl::PointXYZ>& ring,
    double pixelSize,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support)
{
    DirectionEvidence2D evidence;
    const double fitTolerance = std::clamp(
        0.5 * std::max(pixelSize, 0.3), 0.15, 0.60);
    evidence.fitTolerance = fitTolerance;

    evidence.chains = ExtractDirectionChains(ring, fitTolerance);
    evidence.totalChains = static_cast<int>(evidence.chains.size());
    evidence.totalWeightedLength = 0.0;
    for (const auto& c : evidence.chains) evidence.totalWeightedLength += c.weight;
    if (evidence.chains.empty()) return evidence;

    evidence.peaks = DetectDirectionPeaks(evidence.chains, fitTolerance);
    if (evidence.peaks.empty()) return evidence;

    evidence.primary = evidence.peaks[0];
    if (evidence.peaks.size() >= 2) {
        evidence.secondary = evidence.peaks[1];
        // multi-direction requires: secondary ratio >= threshold AND >= 2 independent chains
        // single long chain as secondary is not enough
        if (evidence.secondary.ratio >= kDirectionMinSecondaryRatio &&
            evidence.secondary.chainCount >= kDirectionMinIndependentChains) {
            evidence.valid = true;
        }
    }
    evidence.source = support && !support->empty() ? 2 : 0;
    return evidence;
}

// ===== 方向系统诊断(只读通道) =====
// 目的: 验证"加权 KDE + 稳定边链贪心聚类"的混合方向检测是否比现有
// 全局 KDE 更稳。只输出 [DirectionDiag] 日志，不改变
// DirectionDecision / TopologyPreservingRegularize / Ceres 的任何行为。
//
// 设计要点:
// - 检测对象是 ExtractDirectionChains 的稳定边链(非原始像素短边)，
//   链自带长度/拟合RMSE/端点坐标，楼梯残差已在链层面被吸收;
// - 方向系统 = theta / theta+90° 一组: 折叠角空间([0°,90°))天然满足
//   平行边与垂直边同系统(垂直边折叠后与平行边重合);
// - 贪心聚类种子从权重前 N 的稳定链中依次考察，不无条件取最长链，
//   并要求系统成立门槛(成员数>=2 或种子足够长);
// - 短链不作种子也不丢弃: 能归组则归组，位于显著转折(凹凸/窄颈)的
//   标记为受保护结构，其余记为自由边/噪声。
namespace {

constexpr double kDirDiagStableMinLength = 1.5;    // m, 稳定链长度下限(高于投票下限 0.8)
constexpr double kDirDiagSeedMinLength = 3.0;      // m, 单链成系统的最小长度(配合连续性/占比约束)
constexpr double kDirDiagSystemMinSeparationDeg = 20.0; // 候选方向间最小角距: 低于此值视为同一方向的栅格抖动(楼梯链角误差±5~8°, 同向簇可差10~15°)
constexpr double kDirDiagAssignDeg = 16.0;         // 稳定链归组阈值(统一分配, 就近候选)
constexpr double kDirDiagShortAssignDeg = 16.0;    // 短链归组阈值
constexpr double kDirDiagLongRescueDeg = 20.0;     // 未归组长链(>=3m)抢救归组上限
constexpr double kDirDiagStructuralTurnDeg = 45.0; // 相邻链转角超过此值=显著凹凸/窄颈
constexpr int kDirDiagMaxCandidates = 7;           // 候选方向上限(防御异常; 90°/12°分离理论约7)
constexpr int kDirDiagMaxFinalSystems = 5;         // 最终可信系统上限(超过标记 complex 不静默丢弃)
constexpr double kDirDiagSystemMinChainFrac = 0.15; // 多链系统最小支持长度占稳定链总长比(8%时楼梯碎片凑数成系统)
constexpr double kDirDiagUnassignedWeakRatio = 0.05;  // 未归组长链弱支持占比阈值(≤此值不回退)
constexpr double kDirDiagUnassignedHighRatio = 0.10;  // 未归组长链占比>此值且无法成新系统→不确定
constexpr double kDirDiagGrayZoneScoreMin = 0.60;     // 灰区(5%~10%)时单/多方向模型评分下限
constexpr double kDirDiagSpikeAngleDeg = 32.0;     // 尖刺内角阈值(<此值且低面积证据=尖刺)
constexpr double kDirDiagSpikeMaxArea = 0.30;      // m², 尖刺三角面积上限
constexpr double kDirDiagZigzagChordDev = 0.12;    // m, 近共线锯齿顶点偏离弦的上限
constexpr double kDirDiagShortDiagMinLen = 0.8;    // m, 短斜边方向检查的长度下限
constexpr double kDirDiagSupportRadius = 0.6;      // m, 支撑点计数半径
constexpr double kDirDiagMultiGainMin = 0.20;      // 多方向判定: multi比single评分高多少
constexpr double kDirDiagMultiFracMin = 0.22;      // 多方向判定: 次系统权重占比下限
// 多方向必须由一个独立且足够大的次方向系统证明。仅凭评分增益或
// 一簇短链会把栅格楼梯误判成多方向；这两个门槛比旧的 0.18/15m
// 更保守，优先保证单方向建筑不会进入 AllowDiagonal。
constexpr double kDirDiagRobustSecondaryFrac = 0.22;
constexpr double kDirDiagRobustSecondaryConfidence = 0.22;
constexpr int kDirDiagRobustSecondaryChains = 2;

// ---- 双残差: 平滑轮廓定拓扑, 原始像素轮廓提供 Ceres 几何观测 ----
// 总权比 raw:smooth = 70:30；逐点权重按实际关联点数归一化。
constexpr bool kUseRawMaskBoundaryResiduals = true;
constexpr double kRawResidualWeight = 0.70;     // raw 残差总权重占比
constexpr double kSmoothResidualWeight = 0.30;  // smooth 残差总权重占比
constexpr int kRawResidualMinPoints = 24;       // raw 采样点下限(低于则回退)
constexpr int kRawResidualMaxPoints = 5000;     // raw 采样点上限(按弧长占比分配)
constexpr double kRawResidualSpacingMin = 0.3;  // m, raw 采样间距下限
constexpr double kRawResidualSpacingMax = 0.8;  // m, raw 采样间距上限
constexpr bool kDumpRawResidualPoints = true;   // 调试点输出开关
constexpr std::size_t kRawResidualDebugMaxPoints = 250000;

struct DirectionSystemDiag {
    double angleRad = 0.0;        // 系统代表方向(折叠[0°,90°))
    int chainCount = 0;           // 归入的稳定链数
    int shortRegrouped = 0;       // 归入的短链数
    double totalLength = 0.0;     // 稳定链总长度
    double weight = 0.0;          // 加权支持(长度×拟合质量)
    double meanRmse = 0.0;        // 平均拟合残差
    double meanContinuity = 0.0;  // 平均连续性(弦长/弧长, 1=直线)
    double concentration = 0.0;   // 圆集中度(四倍角合向量长度R, 1=成员角完全集中)
    double extent = 0.0;          // 空间覆盖: 成员端点 bbox 对角线
    double kdePeak = 0.0;         // 系统角度处的 KDE 密度
    long long supportPoints = 0;  // 成员链关联的支撑点数
    double confidence = 0.0;      // 方向置信度 [0,1]
};

double DirDiagPointToSegmentDistance(
    double px, double py,
    const pcl::PointXYZ& a, const pcl::PointXYZ& b)
{
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double lenSq = dx * dx + dy * dy;
    if (lenSq < 1e-12) return std::hypot(px - a.x, py - a.y);
    double t = ((px - a.x) * dx + (py - a.y) * dy) / lenSq;
    t = std::clamp(t, 0.0, 1.0);
    return std::hypot(px - a.x - t * dx, py - a.y - t * dy);
}

constexpr double kDirDiagCertaintyConf = 0.15;   // 主系统置信度下限, 低于则方向不确定
constexpr double kDirDiagCertaintyScore = 0.30;  // 单方向评分下限
constexpr double kDirDiagMultiSecondConf = 0.18; // 多方向判定: 次系统置信度下限
constexpr double kDirDiagMultiSecondLength = 15.0; // 多方向判定: 次系统绝对支持长度旁路(m)

struct DirectionChainAssignment {
    bool stable = false;
    double continuity = 0.0;   // 端点直线距/弧长, 1=完美直线
    long long supportPts = 0;
    int system = -1;           // 归组系统, -1=未归组(自由)
    bool protectedStructure = false;
};

struct DirectionSystemBuild {
    std::vector<DirectionSystemDiag> systems;
    std::vector<DirectionChainAssignment> chainInfo;  // 与输入链平行
    int shortRegrouped = 0;
    int shortProtected = 0;
    int shortFree = 0;
    int rescuedLongChains = 0;   // 抢救阈值内归入系统的 >=3m 长链
    int unassignedLongCount = 0;     // 最终仍未归组的 >=3m 稳定长链数
    double unassignedLongLength = 0.0;  // 未归组长链总长(m)
    double unassignedLongWeight = 0.0;  // 未归组长链总权重
    double unassignedLengthRatio = 0.0; // 未归组长链总长 / 稳定链总长
    double unassignedWeightRatio = 0.0; // 未归组长链权重 / 稳定链总权重
    double maxUnassignedLength = 0.0;   // 最长未归组长链(m)
    bool complexDirectionEvidence = false; // 可信系统数超过上限(证据复杂, 未静默丢弃)
    double totalStableWeight = 0.0;
    double totalStableLength = 0.0;    // 稳定链总长
    double totalChainLength = 0.0;     // 全部链总长(≈周长, 供占比判定)
    double singleScore = 0.0;
    double multiScore = 0.0;
    bool multiDirection = false;
    bool directionCertain = false;
    std::string uncertainReason;  // directionCertain=false 的原因
};

// 稳定边链贪心聚类 + 加权KDE 的方向系统构建。
// 方向诊断(只读)与拓扑保持通道(方向接入)共用同一实现。
// 系统判定综合考虑: 支持链总长度(权重占比)、拟合残差(置信度内的
// fitFactor)、方向模型改善幅度(评分增益)、空间范围(systems.extent);
// supportCloud 为可用的支撑/墙面点(mask-only 无 OSGB 时传空)。
DirectionSystemBuild BuildDirectionSystems(
    const std::vector<DirectionChain>& chains,
    double fitTolerance,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& supportCloud)
{
    DirectionSystemBuild build;
    const std::size_t n = chains.size();

    // 逐链指标: 稳定性/连续性/支撑点数 + 稳定链总量立即累计
    // (systemCredible 依赖 totalStableLength, 必须在候选验证前就绪;
    //   此前该值在未归组核算阶段才计算, 初筛时≈0 导致多链系统
    //   的长度门槛失效, 弱系统先被保留再靠 DirectionPrune 补救)
    build.chainInfo.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& c = chains[i];
        auto& info = build.chainInfo[i];
        build.totalChainLength += c.length;
        const double chord = std::hypot(c.endPoint.x - c.startPoint.x,
                                        c.endPoint.y - c.startPoint.y);
        info.continuity = chord / std::max(c.length, 1e-9);
        info.stable = c.length >= kDirDiagStableMinLength &&
                      c.rmse <= fitTolerance;
        if (info.stable) {
            build.totalStableLength += c.length;
            build.totalStableWeight += c.weight;
        }
        if (supportCloud && !supportCloud->empty()) {
            for (const auto& p : supportCloud->points) {
                if (DirDiagPointToSegmentDistance(
                        p.x, p.y, c.startPoint, c.endPoint) <=
                        kDirDiagSupportRadius) {
                    ++info.supportPts;
                }
            }
        }
    }

    // 稳定链加权 KDE 峰
    std::vector<DirectionChain> stableChains;
    for (std::size_t i = 0; i < n; ++i) {
        if (build.chainInfo[i].stable) stableChains.push_back(chains[i]);
    }
    const auto peaks = DetectDirectionPeaks(stableChains, fitTolerance);

    // 多阶段方向系统聚类: 候选发现与最终保留彻底分离。
    // 阶段一(候选发现): 全部稳定链按权重降序, 12°分离, 候选上限
    //   kDirDiagMaxCandidates(防御异常, 不在生成第3个时停止)。
    // 阶段二(初始分配): 就近候选 9° 归组, 圆均值精化。
    // 阶段三(候选验证): 弱候选删除后对全部稳定链重新就近分配,
    //   再精化再分配, 迭代至稳定(≤3轮)——避免删除后留下的
    //   归组真空和漂移角。
    // 阶段四(合并): 精化后角距<8°的系统合并(索引统一重映射)。
    // 阶段五(未归组核算): 长链能成新可信系统则新增重分配;
    //   按支持占比(而非数量)决定是否方向不确定。
    std::vector<DirectionSystemDiag>& systems = build.systems;
    const double minSepRad = kDirDiagSystemMinSeparationDeg * M_PI / 180.0;
    const double assignRad = kDirDiagAssignDeg * M_PI / 180.0;
    const double rescueRad = kDirDiagLongRescueDeg * M_PI / 180.0;

    // 就近分配全部稳定链(9°归组, >=3m 长链 18° 抢救)
    auto assignAllStable = [&]() {
        build.rescuedLongChains = 0;
        for (auto& info : build.chainInfo) {
            if (info.stable) info.system = -1;
        }
        if (systems.empty()) return;
        for (std::size_t i = 0; i < n; ++i) {
            if (!build.chainInfo[i].stable) continue;
            const auto& c = chains[i];
            std::size_t nearest = 0;
            double nearestDist = foldedAngleDistance90(
                c.angleRad, systems.front().angleRad);
            for (std::size_t s = 1; s < systems.size(); ++s) {
                const double d = foldedAngleDistance90(
                    c.angleRad, systems[s].angleRad);
                if (d < nearestDist) { nearestDist = d; nearest = s; }
            }
            if (nearestDist <= assignRad) {
                build.chainInfo[i].system = static_cast<int>(nearest);
            } else if (c.length >= 3.0 && nearestDist <= rescueRad) {
                build.chainInfo[i].system = static_cast<int>(nearest);
                ++build.rescuedLongChains;
            }
        }
    };
    // 汇总系统成员 + 折叠角圆均值精化(算术平均在 1°/89° 相邻角对
    // 上会得到错误的 45°)
    auto aggregateSystems = [&]() {
        std::vector<double> angles(systems.size(), 0.0);
        for (std::size_t s = 0; s < systems.size(); ++s) {
            DirectionSystemDiag fresh;
            fresh.angleRad = systems[s].angleRad;  // 空系统保留原角待清理
            systems[s] = fresh;
            DirectionSystemDiag& sys = systems[s];
            std::vector<double> memberAngles;
            std::vector<double> memberWeights;
            double continuitySum = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                if (build.chainInfo[i].system != static_cast<int>(s)) continue;
                if (!build.chainInfo[i].stable) continue;
                const auto& c = chains[i];
                sys.totalLength += c.length;
                sys.weight += c.weight;
                sys.supportPoints += build.chainInfo[i].supportPts;
                continuitySum += build.chainInfo[i].continuity;
                memberAngles.push_back(c.angleRad);
                memberWeights.push_back(c.weight);
            }
            sys.chainCount = static_cast<int>(memberAngles.size());
            if (!memberAngles.empty()) {
                sys.angleRad = circularMeanAngle90(memberAngles, memberWeights);
                sys.meanContinuity = continuitySum / memberAngles.size();
                // 圆集中度 R = |Σw·e^{i4θ}| / Σw: 真方向系统的成员角
                // 高度集中(平行墙折叠后 ±3°, R≈1); 楼梯/噪声碎片的
                // 角度簇分散(R 低)——用 R 区分真方向与凑数簇
                double cx4 = 0.0, cy4 = 0.0, wSum = 0.0;
                for (std::size_t k = 0; k < memberAngles.size(); ++k) {
                    cx4 += memberWeights[k] * std::cos(4.0 * memberAngles[k]);
                    cy4 += memberWeights[k] * std::sin(4.0 * memberAngles[k]);
                    wSum += memberWeights[k];
                }
                sys.concentration =
                    wSum > 1e-9 ? std::hypot(cx4, cy4) / wSum : 0.0;
            }
            angles[s] = sys.angleRad;
        }
    };
    // 删除指定系统并统一重映射 chainInfo.system 索引
    // (donor<keeper 时 erase 会使 keeper 及其后继索引前移)
    auto removeSystem = [&](std::size_t donor, std::size_t keeper) {
        // keeper 无效(无幸存系统/donor 自身)时成员置为未分配,
        // 否则全灭路径会留下悬空 system 索引, 后续 remap 越界崩溃
        const bool hasKeeper =
            keeper < systems.size() && keeper != donor;
        const int keeperAfter = hasKeeper
            ? static_cast<int>(keeper > donor ? keeper - 1 : keeper)
            : -1;
        for (auto& info : build.chainInfo) {
            int& s = info.system;
            if (s < 0) continue;
            if (s == static_cast<int>(donor)) s = keeperAfter;
            else if (s > static_cast<int>(donor)) --s;
        }
        systems.erase(systems.begin() + donor);
    };
    // 候选可信度: ≥2 条稳定链且支持长度占比达标;
    // 或单条可信长链(连续性≥0.85, 长度≥max(15m, 周长10%))。
    // 所有路线都要求圆集中度 R≥0.80(成员角集中才是真方向系统,
    // 楼梯/噪声碎片簇 R 低, 否则单方向建筑会被拆成 4~6 个
    // 垃圾系统——实测 65/94 的质量投诉源于此)
    auto systemCredible = [&](const DirectionSystemDiag& sys) {
        if (sys.chainCount <= 0) return false;
        if (sys.concentration < 0.80) return false;
        if (sys.chainCount >= 2) {
            return sys.totalLength >=
                kDirDiagSystemMinChainFrac *
                std::max(build.totalStableLength, 1e-9);
        }
        const double minSingle = std::max(
            kDirDiagSeedMinLength * 5.0,
            0.10 * std::max(build.totalChainLength, 1e-9));
        return sys.meanContinuity >= 0.85 && sys.totalLength >= minSingle;
    };

    {
        // ---- 阶段一: 候选发现 ----
        std::vector<std::size_t> sorted;
        for (std::size_t i = 0; i < n; ++i) {
            if (build.chainInfo[i].stable) sorted.push_back(i);
        }
        std::sort(sorted.begin(), sorted.end(),
                  [&](std::size_t a, std::size_t b) {
                      return chains[a].weight > chains[b].weight;
                  });
        std::vector<double> seedAngles;
        for (std::size_t seedIdx : sorted) {
            if (static_cast<int>(seedAngles.size()) >= kDirDiagMaxCandidates) break;
            const double seedAngle = chains[seedIdx].angleRad;
            bool nearExisting = false;
            for (double a : seedAngles) {
                if (foldedAngleDistance90(seedAngle, a) < minSepRad) {
                    nearExisting = true;
                    break;
                }
            }
            if (nearExisting) continue;
            // KDE 峰近旁用峰角校正(校正后复查分离角, 防形成近重复)
            double sysAngle = seedAngle;
            for (const auto& pk : peaks) {
                if (foldedAngleDistance90(seedAngle, pk.angleRad) < minSepRad) {
                    bool peakSeparated = true;
                    for (double a : seedAngles) {
                        if (foldedAngleDistance90(pk.angleRad, a) < minSepRad) {
                            peakSeparated = false;
                            break;
                        }
                    }
                    if (peakSeparated) sysAngle = pk.angleRad;
                    break;
                }
            }
            seedAngles.push_back(sysAngle);
        }
        systems.clear();
        systems.resize(seedAngles.size());
        for (std::size_t s = 0; s < seedAngles.size(); ++s) {
            systems[s].angleRad = seedAngles[s];
        }

        // ---- 阶段二+三: 分配 → 精化 → 验证删弱 → 重分配, 迭代至稳定 ----
        for (int iter = 0; iter < 3; ++iter) {
            assignAllStable();
            aggregateSystems();
            // 删除不可信候选(从后往前, 索引安全)
            bool removedAny = false;
            for (std::size_t s = systems.size(); s-- > 0;) {
                if (!systemCredible(systems[s])) {
                    // 并入最近的幸存系统(无幸存则全部置未分配)
                    std::size_t nearest = systems.size();
                    double nearestDist = 1e9;
                    for (std::size_t t = 0; t < systems.size(); ++t) {
                        if (t == s) continue;
                        const double d = foldedAngleDistance90(
                            systems[s].angleRad, systems[t].angleRad);
                        if (d < nearestDist) { nearestDist = d; nearest = t; }
                    }
                    removeSystem(s, nearest < systems.size() ? nearest : 0);
                    removedAny = true;
                }
            }
            if (!removedAny) break;
        }
        // 空系统清理(种子立了但无人归入)
        {
            std::vector<DirectionSystemDiag> kept;
            std::vector<int> remap(systems.size(), -1);
            for (std::size_t s = 0; s < systems.size(); ++s) {
                if (systems[s].chainCount > 0) {
                    remap[s] = static_cast<int>(kept.size());
                    kept.push_back(systems[s]);
                }
            }
            for (auto& info : build.chainInfo) {
                if (info.system >= 0) {
                    // 越界防御: 悬空索引(不应发生)按未分配处理
                    if (static_cast<std::size_t>(info.system) >= remap.size()) {
                        info.system = -1;
                    } else {
                        info.system = remap[static_cast<std::size_t>(info.system)];
                    }
                }
            }
            systems.swap(kept);
        }

        // ---- 阶段四: 精化后角距<8°的系统合并 ----
        {
            const double mergeRad = 14.0 * M_PI / 180.0;
            int mergeGuard = 0;
            bool merged = true;
            while (merged && systems.size() >= 2) {
                merged = false;
                std::size_t bi = 0, bj = 1;
                double bestDist = 1e9;
                for (std::size_t a = 0; a < systems.size(); ++a) {
                    for (std::size_t b = a + 1; b < systems.size(); ++b) {
                        const double d = foldedAngleDistance90(
                            systems[a].angleRad, systems[b].angleRad);
                        if (d < bestDist) { bestDist = d; bi = a; bj = b; }
                    }
                }
                if (bestDist < mergeRad) {
                    const std::size_t donor =
                        systems[bj].weight > systems[bi].weight ? bi : bj;
                    const std::size_t keeper = donor == bi ? bj : bi;
                    removeSystem(donor, keeper);
                    aggregateSystems();
                    merged = true;
                }
            }
        }

        // ---- 阶段五: 未归组长链核算 ----
        // 5.1 未归组长链尝试组成新的可信系统(≤12° 内聚, 满足可信度)
        //     且最终系统数未超上限 → 新增后全量重分配
        {
            std::vector<std::size_t> unassigned;
            for (std::size_t i = 0; i < n; ++i) {
                if (build.chainInfo[i].stable &&
                    build.chainInfo[i].system < 0 &&
                    chains[i].length >= 3.0) {
                    unassigned.push_back(i);
                }
            }
            if (!unassigned.empty() &&
                static_cast<int>(systems.size()) < kDirDiagMaxFinalSystems) {
                // 未归组长链内部再聚类(权重降序, 12°分离)
                std::sort(unassigned.begin(), unassigned.end(),
                          [&](std::size_t a, std::size_t b) {
                              return chains[a].weight > chains[b].weight;
                          });
                std::vector<std::size_t> group;
                double groupAngle = 0.0;
                for (std::size_t cand : unassigned) {
                    if (group.empty()) {
                        group.push_back(cand);
                        groupAngle = chains[cand].angleRad;
                        continue;
                    }
                    if (foldedAngleDistance90(chains[cand].angleRad, groupAngle) <
                        minSepRad) {
                        group.push_back(cand);
                    }
                }
                // 评估这一组的可信度
                double lenSum = 0.0, wSum = 0.0, contSum = 0.0;
                std::vector<double> ga, gw;
                for (std::size_t m : group) {
                    lenSum += chains[m].length;
                    wSum += chains[m].weight;
                    contSum += build.chainInfo[m].continuity;
                    ga.push_back(chains[m].angleRad);
                    gw.push_back(chains[m].weight);
                }
                const double groupMeanAngle = circularMeanAngle90(ga, gw);
                const double minSingle = std::max(
                    15.0, 0.10 * std::max(build.totalChainLength, 1e-9));
                const bool credible = lenSum >= std::max(minSingle,
                    kDirDiagSystemMinChainFrac *
                    std::max(build.totalStableLength, 1e-9)) &&
                    (group.size() >= 2 ||
                     (contSum / group.size() >= 0.85 && lenSum >= minSingle));
                if (credible && group.size() >= 2) {
                    DirectionSystemDiag newSys;
                    newSys.angleRad = groupMeanAngle;
                    systems.push_back(newSys);
                    assignAllStable();
                    aggregateSystems();
                }
            }
        }
        // 5.2 统计未归组长链占比
        build.unassignedLongCount = 0;
        build.unassignedLongLength = 0.0;
        build.unassignedLongWeight = 0.0;
        build.maxUnassignedLength = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (!build.chainInfo[i].stable || build.chainInfo[i].system >= 0) {
                continue;
            }
            const auto& c = chains[i];
            if (c.length < 3.0) continue;
            ++build.unassignedLongCount;
            build.unassignedLongLength += c.length;
            build.unassignedLongWeight += c.weight;
            build.maxUnassignedLength =
                std::max(build.maxUnassignedLength, c.length);
        }
        // 稳定链总量已在逐链指标阶段累计, 不重复
        build.unassignedLengthRatio =
            build.totalStableLength > 1e-9
                ? build.unassignedLongLength / build.totalStableLength : 0.0;
        build.unassignedWeightRatio =
            build.totalStableWeight > 1e-9
                ? build.unassignedLongWeight / build.totalStableWeight : 0.0;
        // 5.3 可信系统超上限标记(不静默丢弃)
        build.complexDirectionEvidence =
            static_cast<int>(systems.size()) > kDirDiagMaxFinalSystems;
    }

    // 系统统计: rmse/范围/KDE密度/置信度
    // (totalStableWeight/Length 已在未归组核算阶段计算, 不重复累加)
    const double bandwidthRad = kDirectionKdeBandwidthDeg * M_PI / 180.0;
    for (std::size_t s = 0; s < systems.size(); ++s) {
        auto& sys = systems[s];
        double rmseSum = 0.0;
        double minX = 1e18, minY = 1e18, maxX = -1e18, maxY = -1e18;
        for (std::size_t i = 0; i < n; ++i) {
            if (build.chainInfo[i].system != static_cast<int>(s)) continue;
            const auto& c = chains[i];
            if (build.chainInfo[i].stable) rmseSum += c.rmse;
            minX = std::min(minX, (double)c.startPoint.x);
            maxX = std::max(maxX, (double)c.startPoint.x);
            minX = std::min(minX, (double)c.endPoint.x);
            maxX = std::max(maxX, (double)c.endPoint.x);
            minY = std::min(minY, (double)c.startPoint.y);
            maxY = std::max(maxY, (double)c.startPoint.y);
            minY = std::min(minY, (double)c.endPoint.y);
            maxY = std::max(maxY, (double)c.endPoint.y);
        }
        double kdeDensity = 0.0;
        for (const auto& c : chains) {
            kdeDensity += c.weight * std::exp(
                -0.5 * std::pow(foldedAngleDistance90(c.angleRad, sys.angleRad) /
                                bandwidthRad, 2));
        }
        sys.meanRmse = sys.chainCount > 0 ? rmseSum / sys.chainCount : 0.0;
        sys.extent = (minX <= maxX && minY <= maxY)
            ? std::hypot(maxX - minX, maxY - minY) : 0.0;
        sys.kdePeak = kdeDensity;
        const double supportFrac =
            build.totalStableWeight > 1e-9
                ? sys.weight / build.totalStableWeight : 0.0;
        const double fitFactor = std::exp(-sys.meanRmse / std::max(fitTolerance, 0.01));
        sys.confidence = std::clamp(supportFrac * fitFactor, 0.0, 1.0);
    }

    // 系统置信度计算完成后再做一次“方向系统净化”。弱系统不能继续
    // 参与拓扑候选，否则它们虽然各自方向合法，却会给单方向建筑制造
    // 一组没有建筑意义的斜边。保留主系统和真正达到独立支持门槛的次系统，
    // 同时重映射 chainInfo，保证后续格网吸附与 Ceres 使用同一套系统。
    if (systems.size() >= 2 && build.totalStableLength > 1e-9) {
        std::size_t primary = 0;
        for (std::size_t s = 1; s < systems.size(); ++s) {
            if (systems[s].weight > systems[primary].weight) primary = s;
        }

        std::vector<bool> keep(systems.size(), false);
        keep[primary] = true;
        for (std::size_t s = 0; s < systems.size(); ++s) {
            if (s == primary) continue;
            const double lengthFrac = systems[s].totalLength /
                std::max(build.totalStableLength, 1e-9);
            if (systems[s].chainCount >= kDirDiagRobustSecondaryChains &&
                lengthFrac >= kDirDiagRobustSecondaryFrac &&
                systems[s].confidence >= kDirDiagRobustSecondaryConfidence &&
                systems[s].concentration >= 0.80) {
                keep[s] = true;
            }
        }

        std::vector<DirectionSystemDiag> filtered;
        filtered.reserve(systems.size());
        std::vector<int> remap(systems.size(), -1);
        // 主系统固定放在 0 号位；后续单方向代码约定 systems.front()
        // 就是唯一主方向。
        remap[primary] = 0;
        filtered.push_back(systems[primary]);
        for (std::size_t s = 0; s < systems.size(); ++s) {
            if (s == primary || !keep[s]) continue;
            remap[s] = static_cast<int>(filtered.size());
            filtered.push_back(systems[s]);
        }
        for (auto& info : build.chainInfo) {
            if (info.system < 0) continue;
            const auto old = static_cast<std::size_t>(info.system);
            info.system = old < remap.size() ? remap[old] : -1;
            if (info.system < 0 && !filtered.empty()) {
                // 被删除的弱系统后续会按主方向处理；不要留下自由斜边。
                info.system = 0;
            }
        }
        if (filtered.size() != systems.size()) {
            std::cerr << "[DirectionPrune] systems=" << systems.size()
                      << " kept=" << filtered.size()
                      << " primary_angle_deg="
                      << filtered.front().angleRad * 180.0 / M_PI << std::endl;
            systems.swap(filtered);
        }
    }

    // 短链处置: 归组 / 受保护结构 / 自由边
    {
        const double shortAssignRad = kDirDiagShortAssignDeg * M_PI / 180.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (build.chainInfo[i].stable || build.chainInfo[i].system >= 0) continue;
            const auto& c = chains[i];
            int bestSys = -1;
            double bestDist = shortAssignRad;
            for (std::size_t s = 0; s < systems.size(); ++s) {
                const double d = foldedAngleDistance90(
                    c.angleRad, systems[s].angleRad);
                if (d < bestDist) { bestDist = d; bestSys = static_cast<int>(s); }
            }
            if (bestSys >= 0) {
                build.chainInfo[i].system = bestSys;
                ++systems[static_cast<std::size_t>(bestSys)].shortRegrouped;
                ++build.shortRegrouped;
                continue;
            }
            const auto turnDeg = [&](std::size_t a, std::size_t b) {
                return foldedAngleDistance90(
                    chains[a].angleRad, chains[b].angleRad) * 180.0 / M_PI;
            };
            const double turnPrev = turnDeg((i + n - 1) % n, i);
            const double turnNext = turnDeg(i, (i + 1) % n);
            if (std::max(turnPrev, turnNext) >= kDirDiagStructuralTurnDeg) {
                build.chainInfo[i].protectedStructure = true;
                ++build.shortProtected;
            } else {
                ++build.shortFree;
            }
        }
    }

    // 单/多方向模型评分 + 判定
    {
        double totalW = 0.0;
        for (const auto& c : chains) totalW += c.weight;
        if (totalW > 1e-9 && !systems.empty()) {
            const double norm = std::pow(45.0 * M_PI / 180.0, 2);
            auto residualTo = [&](const std::vector<double>& angles) {
                double r = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    double best = foldedAngleDistance90(
                        chains[i].angleRad, angles.front());
                    for (double a : angles) {
                        best = std::min(best, foldedAngleDistance90(
                            chains[i].angleRad, a));
                    }
                    r += chains[i].weight * best * best;
                }
                return r;
            };
            std::vector<double> singleAngle = { systems.front().angleRad };
            std::vector<double> multiAngles;
            for (const auto& s : systems) multiAngles.push_back(s.angleRad);
            build.singleScore = 1.0 - std::sqrt(
                residualTo(singleAngle) / (totalW * norm));
            build.multiScore = 1.0 - std::sqrt(
                residualTo(multiAngles) / (totalW * norm));
            if (systems.size() >= 2) {
                // 找主系统之外支持最强的系统(不一定是列表第 2 个:
                // 中间可能夹着弱系统, 强方向排在第 3 位的建筑真实存在)
                std::size_t best = 1;
                for (std::size_t s = 2; s < systems.size(); ++s) {
                    if (systems[s].weight > systems[best].weight) best = s;
                }
                const DirectionSystemDiag& second = systems[best];
                const double secondFrac =
                    build.totalStableWeight > 1e-9
                        ? second.weight / build.totalStableWeight : 0.0;
                // 多方向判据(两种之一成立即可):
                // (a) 模型改善幅度大: 多系统解释率显著优于单系统;
                // (b) 次系统绝对支持强: 两系统角差小时(如 32°/48° 双翼
                //     建筑)评分增益必然有限, 但大量级的次翼墙长是
                //     不可忽视的真实方向证据。
                const bool byGain =
                    (build.multiScore - build.singleScore) >= kDirDiagMultiGainMin;
                // bySupport 单长链旁路的附加约束(方案A): 单链必须
                // 足够直(连续性)、且总长占轮廓周长一定比例——
                // 真实方向翼的墙连续挺直, 楼梯残差链连续性低
                const double minSupportLength = std::max(
                    kDirDiagMultiSecondLength,
                    0.10 * std::max(build.totalChainLength, 1e-9));
                const bool singleChainCredible =
                    second.chainCount >= 2 ||
                    (second.meanContinuity >= 0.85 &&
                     second.totalLength >= minSupportLength);
                const bool bySupport =
                    second.totalLength >= minSupportLength && singleChainCredible;
                // byGain 路线沿用"≥2 独立链"防孤链噪声;
                // bySupport 路线在附加约束下允许单条可信长墙
                build.multiDirection =
                    (byGain || bySupport) &&
                    secondFrac >= kDirDiagMultiFracMin &&
                    second.chainCount >= (bySupport ? 1 : kDirectionMinIndependentChains) &&
                    second.confidence >= kDirDiagMultiSecondConf;
            }
        }
        // 方向确定性: 支持比例 + 方向组证据, 不做单条边二元否决。
        // - 未归组占比 >10%(长度或权重) 且无法形成新可信系统 → 不确定;
        // - 灰区(5%~10%): 用单/多方向模型评分裁决;
        // - ≤5% 弱支持: 不回退(低支持异常链交由格网吸附兜底)。
        if (systems.empty()) {
            build.uncertainReason = "no_systems";
        } else if (build.unassignedLengthRatio > kDirDiagUnassignedHighRatio ||
                   build.unassignedWeightRatio > kDirDiagUnassignedHighRatio) {
            build.uncertainReason = "unassigned_long_chain";
        } else if (build.unassignedLengthRatio > kDirDiagUnassignedWeakRatio ||
                   build.unassignedWeightRatio > kDirDiagUnassignedWeakRatio) {
            // 灰区: 模型解释率不足才判不确定
            const double bestScore =
                std::max(build.singleScore,
                         build.multiDirection ? build.multiScore : build.singleScore);
            if (bestScore < kDirDiagGrayZoneScoreMin) {
                build.uncertainReason = "unassigned_gray_zone_low_score";
            }
        } else if (systems.front().chainCount < 2) {
            build.uncertainReason = "primary_single_chain";
        } else if (systems.front().confidence < kDirDiagCertaintyConf) {
            build.uncertainReason = "low_confidence";
        } else if (build.singleScore < kDirDiagCertaintyScore) {
            build.uncertainReason = "low_score";
        }
        build.directionCertain =
            build.uncertainReason.empty();
    }
    return build;
}

// 多边形质量检查(拓扑候选/Ceres 结果/VDP 备用结果共用同一标准)。
// initialRing 为参考环。hasDirection=false 时跳过方向检查(无方向
// 上下文的备用结果); multiDirection=true 时 >3m 长边须贴近
// allowedAngles(系统角+自由链角, 折叠角) 12°内, false 时用输出自身
// 加权主导族 10°。返回空串=通过, 否则返回原因。
std::string CheckPolygonQualityVsRing(
    const std::vector<pcl::PointXYZ>& poly,
    const std::vector<pcl::PointXYZ>& initialRing,
    bool hasDirection,
    bool multiDirection,
    const std::vector<double>& systemAngles,
    double maxVertexDisp,
    long long fid = -1,
    int partIndex = 0)
{
    if (poly.size() < 3) return "too_few_vertices";
    if (!isSimplePolygon2D(poly)) return "self_intersecting";
    // 零长边与异常短边链(凹凸台阶的短边是合法拓扑, 只拦退化)
    int shortRun = 0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const double len = std::hypot(
            poly[(i + 1) % poly.size()].x - poly[i].x,
            poly[(i + 1) % poly.size()].y - poly[i].y);
        if (len < 0.02) return "zero_length_edge";
        if (len < 0.15) {
            if (++shortRun >= 3) return "short_edge_run";
        } else {
            shortRun = 0;
        }
    }
    // 参考环统计: 面积/面积加权威心/包围盒对角线/边界距离
    double ringArea = std::abs(polygonArea2D(initialRing));
    double ra = 0.0, rsx = 0.0, rsy = 0.0;
    double ringMinX = 1e18, ringMinY = 1e18, ringMaxX = -1e18, ringMaxY = -1e18;
    for (std::size_t i = 0; i < initialRing.size(); ++i) {
        const auto& p = initialRing[i];
        const auto& q = initialRing[(i + 1) % initialRing.size()];
        const double cross = p.x * q.y - q.x * p.y;
        ra += cross;
        rsx += (p.x + q.x) * cross;
        rsy += (p.y + q.y) * cross;
        ringMinX = std::min(ringMinX, (double)p.x);
        ringMaxX = std::max(ringMaxX, (double)p.x);
        ringMinY = std::min(ringMinY, (double)p.y);
        ringMaxY = std::max(ringMaxY, (double)p.y);
    }
    double ringCx, ringCy;
    if (std::abs(ra) > 1e-12) {
        ringCx = rsx / (3.0 * ra);
        ringCy = rsy / (3.0 * ra);
    } else {
        ringCx = 0.0; ringCy = 0.0;
        for (const auto& p : initialRing) { ringCx += p.x; ringCy += p.y; }
        ringCx /= static_cast<double>(initialRing.size());
        ringCy /= static_cast<double>(initialRing.size());
    }
    const double ringDiagonal = std::hypot(ringMaxX - ringMinX, ringMaxY - ringMinY);
    auto distToRing = [&](const pcl::PointXYZ& pt) {
        double best = std::numeric_limits<double>::max();
        for (std::size_t i = 0; i < initialRing.size(); ++i) {
            const auto& a = initialRing[i];
            const auto& b = initialRing[(i + 1) % initialRing.size()];
            const double dx = b.x - a.x;
            const double dy = b.y - a.y;
            const double lenSq = dx * dx + dy * dy;
            if (lenSq < 1e-12) continue;
            double t = ((pt.x - a.x) * dx + (pt.y - a.y) * dy) / lenSq;
            t = std::clamp(t, 0.0, 1.0);
            const double d = std::hypot(pt.x - a.x - t * dx, pt.y - a.y - t * dy);
            best = std::min(best, d);
        }
        return best;
    };

    const double polyArea = std::abs(polygonArea2D(poly));
    const double ratio = polyArea / std::max(ringArea, 1e-6);
    if (ratio < 0.60 || ratio > 1.60) return "area_ratio=" + std::to_string(ratio);
    // 异常长边(飞点特征)
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const double len = std::hypot(
            poly[(i + 1) % poly.size()].x - poly[i].x,
            poly[(i + 1) % poly.size()].y - poly[i].y);
        if (len > 1.1 * ringDiagonal) return "abnormal_long_edge";
    }
    // 质心位移(面积加权): 规则化不应整体搬动建筑
    {
        double pa2 = 0.0, psx = 0.0, psy = 0.0;
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const auto& p = poly[i];
            const auto& q = poly[(i + 1) % poly.size()];
            const double cross = p.x * q.y - q.x * p.y;
            pa2 += cross;
            psx += (p.x + q.x) * cross;
            psy += (p.y + q.y) * cross;
        }
        double cx, cy;
        if (std::abs(pa2) > 1e-12) {
            cx = psx / (3.0 * pa2);
            cy = psy / (3.0 * pa2);
        } else {
            cx = 0.0; cy = 0.0;
            for (const auto& p : poly) { cx += p.x; cy += p.y; }
            cx /= static_cast<double>(poly.size());
            cy /= static_cast<double>(poly.size());
        }
        if (std::hypot(cx - ringCx, cy - ringCy) > 1.0) return "centroid_shift";
    }
    // 墙面贴合(边中点到环边界<=1.2m)与顶点位移(飞点)
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const auto& p1 = poly[i];
        const auto& p2 = poly[(i + 1) % poly.size()];
        pcl::PointXYZ mid;
        mid.x = 0.5f * (p1.x + p2.x);
        mid.y = 0.5f * (p1.y + p2.y);
        mid.z = p1.z;
        if (distToRing(mid) > 1.20) return "wall_deviation";
    }
    for (const auto& p : poly) {
        if (distToRing(p) > maxVertexDisp) return "flying_vertex";
    }
    // 方向误差: 所有可见直线边(>=0.8m)都须贴近检测到的方向系统。
    // 0.8m 以下通常是像素噪声或受保护的微小结构，不让更长的短斜边
    // 通过质量闸门。
    // 多方向: 容差 12°(优化基准角允许在系统角附近小幅精化),
    //   合法角只有 systemAngles——未归组链不得自我合法化。
    // 单方向: 基准必须是检测到的主系统角 systemAngles.front()(8°),
    //   只有方向上下文缺失时才退回"从结果自身估计主导族"的旧行为。
    if (!hasDirection) return "";
    struct LongEdge { double len; double ang; };
    std::vector<LongEdge> longEdges;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const auto& p1 = poly[i];
        const auto& p2 = poly[(i + 1) % poly.size()];
        const double len = std::hypot(p2.x - p1.x, p2.y - p1.y);
        if (len >= kDirDiagShortDiagMinLen) {
            longEdges.push_back({len, std::atan2(p2.y - p1.y, p2.x - p1.x)});
        }
    }
    if (longEdges.empty()) return "";
    double maxErrDeg = 0.0;
    if (multiDirection) {
        if (systemAngles.empty()) return "";
        for (const auto& e : longEdges) {
            double best = 1e9;
            for (double a : systemAngles) {
                best = std::min(best,
                    foldedAngleDistance90(e.ang, a) * 180.0 / M_PI);
            }
            maxErrDeg = std::max(maxErrDeg, best);
            if (best > 12.0) return "direction_violation";
        }
        if (fid >= 0) {
            std::cerr << "[DirectionQuality] fid=" << fid << " part=" << partIndex
                      << " mode=multi expected_deg=";
            for (std::size_t k = 0; k < systemAngles.size(); ++k) {
                std::cerr << (k ? "/" : "")
                          << systemAngles[k] * 180.0 / M_PI;
            }
            std::cerr << " max_error_deg=" << maxErrDeg << " accepted=1" << std::endl;
        }
        return "";
    }
    // 单方向
    double expectedDeg = 0.0;
    bool usedDetected = !systemAngles.empty();
    if (usedDetected) {
        expectedDeg = systemAngles.front() * 180.0 / M_PI;
        for (const auto& e : longEdges) {
            const double err =
                foldedAngleDistance90(e.ang, systemAngles.front()) * 180.0 / M_PI;
            maxErrDeg = std::max(maxErrDeg, err);
            if (err > 8.0) return "direction_violation";
        }
    } else {
        // 上下文缺失: 从结果自身估计主导族(旧行为兜底)
        double bestRef = 0.0, bestSupport = -1.0;
        for (const auto& cand : longEdges) {
            double support = 0.0;
            for (const auto& e : longEdges) {
                if (foldedAngleDistance90(e.ang, cand.ang) * 180.0 / M_PI <= 10.0) {
                    support += e.len;
                }
            }
            if (support > bestSupport) { bestSupport = support; bestRef = cand.ang; }
        }
        expectedDeg = bestRef * 180.0 / M_PI;
        for (const auto& e : longEdges) {
            const double err =
                foldedAngleDistance90(e.ang, bestRef) * 180.0 / M_PI;
            maxErrDeg = std::max(maxErrDeg, err);
            if (err > 10.0) return "direction_violation";
        }
    }
    if (fid >= 0) {
        std::cerr << "[DirectionQuality] fid=" << fid << " part=" << partIndex
                  << " mode=single expected_deg=" << expectedDeg
                  << (usedDetected ? "" : "(self)")
                  << " max_error_deg=" << maxErrDeg << " accepted=1" << std::endl;
    }
    return "";
}

// ---- 双残差调试点收集器(可选输出 debug_mask_raw_residual_points.shp) ----
struct RawResidualDebugPoint {
    long long fid = 0;
    int part = 0;
    pcl::PointXYZ pt{};
    int src = 0;         // 0=smooth, 1=raw
    double weight = 1.0;
    int assocEdge = -1;  // 关联的候选边索引(-1=未关联)
};
static std::vector<RawResidualDebugPoint>& RawResidualDebugPoints() {
    static std::vector<RawResidualDebugPoint> instance;
    return instance;
}


// A second geometric direction is only credible when the support cloud has an
// independent direction peak as well.  Polygon edges alone are unreliable on
// raster staircases and short repaired notches.
struct SupportDirectionPeaks2D {
    bool valid = false;
    double primaryAngle = 0.0;
    double primaryRatio = 0.0;
    double secondaryAngle = 0.0;
    double secondaryRatio = 0.0;
    std::size_t pairCount = 0;
};

SupportDirectionPeaks2D estimateSupportDirectionPeaks2D(
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& support)
{
    SupportDirectionPeaks2D result;
    constexpr std::size_t kMaxSamples = 5000;
    if (!support || support->size() < kSupportDirectionMinPoints) return result;

    pcl::PointCloud<pcl::PointXYZ>::Ptr support2D(
        new pcl::PointCloud<pcl::PointXYZ>);
    support2D->resize(support->size());
    for (std::size_t i = 0; i < support->size(); ++i) {
        support2D->points[i].x = support->points[i].x;
        support2D->points[i].y = support->points[i].y;
        support2D->points[i].z = 0.0f;
    }

    pcl::KdTreeFLANN<pcl::PointXYZ> tree;
    tree.setInputCloud(support2D);
    std::array<double, 90> bins{};
    const std::size_t stride = std::max<std::size_t>(
        1, support->size() / kMaxSamples);
    std::vector<int> indices;
    std::vector<float> distances;
    double totalWeight = 0.0;

    for (std::size_t i = 0; i < support->size(); i += stride) {
        indices.clear();
        distances.clear();
        tree.radiusSearch(support2D->points[i],
            static_cast<float>(kSupportDirectionPairRadius),
            indices, distances, 96);
        for (std::size_t k = 0; k < indices.size(); ++k) {
            const int j = indices[k];
            if (j <= static_cast<int>(i) ||
                distances[k] < kSupportDirectionMinPairDistance *
                    kSupportDirectionMinPairDistance) {
                continue;
            }
            const auto& a = support->points[i];
            const auto& b = support->points[static_cast<std::size_t>(j)];
            const double degrees = foldedLineAngle90(
                std::atan2(b.y - a.y, b.x - a.x)) * 180.0 / M_PI;
            const int bin = std::min(89, std::max(0,
                static_cast<int>(std::floor(degrees))));
            const double weight = std::sqrt(std::max(0.0f, distances[k]));
            bins[static_cast<std::size_t>(bin)] += weight;
            totalWeight += weight;
            ++result.pairCount;
        }
    }
    if (result.pairCount < kSupportDirectionMinPairs ||
        totalWeight <= 1e-9) {
        return result;
    }

    std::array<double, 90> windows{};
    for (int center = 0; center < 90; ++center) {
        for (int d = -5; d <= 5; ++d) {
            windows[static_cast<std::size_t>(center)] +=
                bins[static_cast<std::size_t>((center + d + 90) % 90)];
        }
    }

    auto bestPeak = [&windows](int excludedCenter) {
        int best = 0;
        double bestValue = -1.0;
        for (int center = 0; center < 90; ++center) {
            const int distance = std::abs(center - excludedCenter);
            const int circularDistance = std::min(distance, 90 - distance);
            if (excludedCenter >= 0 && circularDistance < 15) continue;
            if (windows[static_cast<std::size_t>(center)] > bestValue) {
                bestValue = windows[static_cast<std::size_t>(center)];
                best = center;
            }
        }
        return std::pair<int, double>(best, bestValue);
    };

    const auto primary = bestPeak(-1);
    const auto secondary = bestPeak(primary.first);
    result.primaryAngle = (primary.first + 0.5) * M_PI / 180.0;
    result.secondaryAngle = (secondary.first + 0.5) * M_PI / 180.0;
    result.primaryRatio = primary.second / totalWeight;
    result.secondaryRatio = secondary.second / totalWeight;
    result.valid = result.primaryRatio >= kSupportDirectionStrongPeakRatio &&
        result.secondaryRatio >= kSupportDirectionStrongPeakRatio;
    return result;
}
}

// Close the outer file-local helper namespace.  The direction diagnostic
// helpers above use a nested anonymous namespace, while the class member
// implementations below must remain in the outlineRegular namespace.


}

bool outlineRegular::BuildStrictDirectionalFallback(
    const std::vector<pcl::PointXYZ>& input,
    double mainAngle,
    std::vector<pcl::PointXYZ>& result) const
{
    result.clear();
    if (input.size() < 3 || !std::isfinite(mainAngle)) return false;
    return forceOrthogonalPolygonToAngle(input, mainAngle, result);
}

void outlineRegular::RunDirectionSystemDiagnostic(
    long long fid,
    const std::vector<pcl::PointXYZ>& ring,
    double pixelSize,
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& supportCloud)
{
    // 1. 稳定边链提取(与拓扑通道同参数, 保证诊断与实跑几何一致)
    const double fitTolerance = std::clamp(
        std::max(0.50, 1.2 * pixelSize), 0.15, 0.90);
    const auto chains = ExtractDirectionChains(ring, fitTolerance);
    if (chains.size() < 3) return;

    // 2. 方向系统构建(与拓扑通道共用同一实现)
    const auto build = BuildDirectionSystems(chains, fitTolerance, supportCloud);

    double chordSum = 0.0;
    int shortTotal = 0;
    for (const auto& info : build.chainInfo) {
        chordSum += info.continuity;
        if (!info.stable) ++shortTotal;
    }
    int stableCount = 0;
    for (const auto& info : build.chainInfo) {
        if (info.stable) ++stableCount;
    }

    // 3. 日志输出(格式统一, 带 fid)
    std::cerr << "[DirectionDiag] fid=" << fid
              << " chains=" << chains.size()
              << " stable=" << stableCount
              << " short=" << shortTotal
              << " systems=" << build.systems.size()
              << " kde_peaks=" << build.systems.size()
              << " mean_continuity=" << (chains.empty() ? 0.0 : chordSum / chains.size())
              << " certain=" << (build.directionCertain ? 1 : 0)
              << std::endl;
    for (std::size_t s = 0; s < build.systems.size(); ++s) {
        const auto& sys = build.systems[s];
        std::cerr << "[DirectionDiag] fid=" << fid
                  << " system=" << s
                  << " angle_deg=" << sys.angleRad * 180.0 / M_PI
                  << " chains=" << sys.chainCount
                  << " short_regrouped=" << sys.shortRegrouped
                  << " total_length=" << sys.totalLength
                  << " mean_rmse=" << sys.meanRmse
                  << " extent=" << sys.extent
                  << " kde_peak=" << sys.kdePeak
                  << " support_pts=" << sys.supportPoints
                  << " conf=" << sys.confidence
                  << std::endl;
    }
    if (shortTotal > 0) {
        std::cerr << "[DirectionDiag] fid=" << fid
                  << " short_regrouped=" << build.shortRegrouped
                  << " short_protected=" << build.shortProtected
                  << " short_free=" << build.shortFree
                  << std::endl;
    }
    const char* verdict = !build.directionCertain ? "uncertain"
        : (build.multiDirection ? "multi" : "single");
    std::cerr << "[DirectionDiag] fid=" << fid
              << " single_score=" << build.singleScore
              << " multi_score=" << build.multiScore
              << " verdict=" << verdict
              << std::endl;
    // 多方向案例补链级明细, 便于核对短斜边归属
    if (build.multiDirection) {
        for (std::size_t i = 0; i < chains.size(); ++i) {
            const auto& c = chains[i];
            std::cerr << "[DirectionDiag] fid=" << fid
                      << " chain=" << i
                      << " len=" << c.length
                      << " rmse=" << c.rmse
                      << " cont=" << build.chainInfo[i].continuity
                      << " ang_deg=" << c.angleRad * 180.0 / M_PI
                      << " sys=" << build.chainInfo[i].system
                      << (build.chainInfo[i].protectedStructure ? " PROTECTED" : "")
                      << std::endl;
        }
    }
}


// 结构感知假设修复的运行级汇总(实现在上面的
// 匿名命名空间里)。
void outlineRegular::PrintHypothesisRepairSummary()
{
    PrintHypothesisRepairStats();
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
    std::cerr << "[RegTime] === polygon: fid=" << source_feature_id_
              << " original_points=" << original_points.size()
              << " support=" << fitting_cloud->size() << " ===" << std::endl;

    // 计算点云分辨率 (平均点间距)
    XG::modularFunction mfn;
    double resolution = mfn.computeModelResolution(*fitting_cloud);
    if (!std::isfinite(resolution) || resolution <= 1e-6) {
        resolution = makeOutlineTuning(resolution, computeOBBArea(original_points)).resolution;
    }
    const std::vector<PreservedArcSegment> preserved_arcs =
        curve_restoration_enabled_
            ? detectPreservedArcs(original_points, resolution)
            : std::vector<PreservedArcSegment>{};

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

        // 结构感知假设修复(R1)：VDP 把真实矩形缺口/
        // rectangular notches/protrusions into single spike vertices, which
        // the downstream orthogonal chain cannot recover. Runs AFTER the
        // fallback_hypothesis 捕获之前执行
        // (拓扑修复不会吃掉新建结构，
        // 验收参考与缺口保护都能看到它们)。
        {
            const std::size_t beforeRepair = best_hypothesis.size();
            RepairHypothesisStructures(
                best_hypothesis, original_points, fitting_cloud, source_feature_id_);
            if (best_hypothesis.size() != beforeRepair) {
                best_energy_hypothesis_ = best_hypothesis;
                input_cloud->points.assign(
                    best_hypothesis.begin(), best_hypothesis.end());
            }
        }

        std::vector<pcl::PointXYZ> fallback_hypothesis = best_hypothesis;

        const OutlineTuning model_tuning = makeOutlineTuning(
            resolution, polygonArea2D(best_hypothesis));

        // ===== 多方向证据前置 =====
        // 在尝试单方向正则化之前，先在最优假设上检测方向系统证据：
        // 若存在可信的第二方向系统（斜翼、L 形翼等），直接进入多方向
        // 分支，避免"单方向宽松放行"掩盖斜翼被强行正交化的系统性形变
        // (实测案例：IoU=0.69 / mean=4.1m / q90=11.7m 被 relaxed 通道放行)。
        // 检测不到多方向证据的建筑仍走原 SingleFirst 路径，行为不变。

        // ---- 边链方向证据(split-and-merge + KDE 多峰) ----
        const DirectionEvidence2D chainEvidence =
            EstimateDirectionEvidence(original_points, resolution, fitting_cloud);
        std::cerr << "[DirectionEvidence] fid=" << source_feature_id_
            << " chains=" << chainEvidence.totalChains
            << " total_weight=" << chainEvidence.totalWeightedLength
            << " fit_tol=" << chainEvidence.fitTolerance
            << " source=" << chainEvidence.source << std::endl;
        for (std::size_t pk = 0; pk < chainEvidence.peaks.size(); ++pk) {
            const auto& peak = chainEvidence.peaks[pk];
            std::cerr << "[DirectionPeak] fid=" << source_feature_id_
                << " peak=" << pk
                << " angle=" << peak.angleRad * 180.0 / M_PI
                << " ratio=" << peak.ratio
                << " chains=" << peak.chainCount
                << " rmse=" << peak.meanRmse
                << " strong=" << (peak.strong ? 1 : 0) << std::endl;
        }

        std::vector<DirectionSystem> direction_systems =
            detectDirectionSystems(best_hypothesis, fitting_cloud, model_tuning);
        const SupportDirectionPeaks2D supportPeaks =
            estimateSupportDirectionPeaks2D(fitting_cloud);
        const bool independentWallSecondary =
            supportPeaks.valid &&
            foldedAngleDistance90(
                supportPeaks.primaryAngle, supportPeaks.secondaryAngle) >=
                15.0 * M_PI / 180.0;
        std::cerr << "[MultiDirectionWallEvidence] valid="
            << (supportPeaks.valid ? 1 : 0)
            << " primary_deg=" << supportPeaks.primaryAngle * 180.0 / M_PI
            << " primary_ratio=" << supportPeaks.primaryRatio
            << " secondary_deg=" << supportPeaks.secondaryAngle * 180.0 / M_PI
            << " secondary_ratio=" << supportPeaks.secondaryRatio
            << " pairs=" << supportPeaks.pairCount
            << " independent=" << (independentWallSecondary ? 1 : 0)
            << std::endl;

        bool geometric_multi_direction =
            direction_systems.size() >= 2 &&
            hasCredibleMultiDirectionChains(best_hypothesis, direction_systems, model_tuning);
        // 边链证据门控：次方向峰必须有权重占比和独立链数支撑
        const bool chainMultiDirection = chainEvidence.valid;
        std::cerr << "[MultiChainGate] fid=" << source_feature_id_
            << " chain_multi=" << (chainMultiDirection ? 1 : 0)
            << " chain_secondary_ratio=" << chainEvidence.secondary.ratio
            << " chain_secondary_chains=" << chainEvidence.secondary.chainCount
            << " geometric_multi=" << (geometric_multi_direction ? 1 : 0)
            << " wall_secondary=" << (independentWallSecondary ? 1 : 0) << std::endl;

        bool credible_multi_direction =
            (geometric_multi_direction || chainMultiDirection) &&
            independentWallSecondary;
        if ((geometric_multi_direction || chainMultiDirection) && !independentWallSecondary) {
            std::cerr << "[BuildingMode] strict multi-direction rejected: "
                         "no independent wall secondary peak" << std::endl;
        }
        if (!credible_multi_direction) {
            std::vector<DirectionSystem> relaxed_systems =
                buildMultiDirectionCandidates(best_hypothesis, resolution);
            const bool relaxedGeometricEvidence = relaxed_systems.size() >= 2 &&
                hasCredibleMultiDirectionChains(
                    best_hypothesis, relaxed_systems, model_tuning);
            const bool relaxedChainEvidence = chainMultiDirection;
            if ((relaxedGeometricEvidence || relaxedChainEvidence) && independentWallSecondary) {
                if (relaxed_systems.size() >= 2) {
                    direction_systems = std::move(relaxed_systems);
                    credible_multi_direction = true;
                    std::cerr << "[BuildingMode] relaxed multi-direction evidence accepted" << std::endl;
                } else {
                    std::cerr << "[BuildingModeGuard] rejected_multi reason=less_than_two_systems"
                              << " (relaxed systems=" << relaxed_systems.size() << ")"
                              << std::endl;
                }
            } else if ((relaxedGeometricEvidence || relaxedChainEvidence) && !independentWallSecondary) {
                std::cerr << "[BuildingMode] relaxed multi-direction rejected: "
                             "geometry-only secondary direction" << std::endl;
            }
        }
        // 最终守卫: 多方向必须建立在 >=2 个方向系统之上,
        // 单系统 AllowDiagonal 会输出无系统约束的自由斜边
        if (credible_multi_direction && direction_systems.size() < 2) {
            credible_multi_direction = false;
            std::cerr << "[BuildingModeGuard] rejected_multi reason=less_than_two_systems"
                      << " (final systems=" << direction_systems.size() << ")"
                      << std::endl;
        }

        // 小轮廓(按最优假设面积判定，取代原先
        // main.cpp 的矩形快通道)：强制单一主方向
        // 规则化，避免多方向分支产生
        // 斜边。
        const double smallBuildingHypothesisArea =
            polygonArea2D(fallback_hypothesis);
        const bool forceSingleDirection =
            smallBuildingHypothesisArea > 0.0 &&
            smallBuildingHypothesisArea < kSmallBuildingSingleDirectionArea;
        if (forceSingleDirection && credible_multi_direction) {
            std::cerr << "[SmallBuilding] hypothesis_area=" << smallBuildingHypothesisArea
                      << " < " << kSmallBuildingSingleDirectionArea
                      << " -> force single direction" << std::endl;
            credible_multi_direction = false;
        }

        // 优先尝试单一正交方向。掩膜轮廓常含
        // 楼梯或短斜边伪影，容易让次方向
        // 检测过于激进。若单方向规则化
        // 仍能很好解释最优假设，则保留，
        // 避开稳定性较差的多方向分支。
        // 有可信多方向证据的建筑(上面已查)跳过本块。
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
                constexpr double kWallDirectionCorrectionMax =
                    kWallDirectionCorrectionMaxDeg * M_PI / 180.0;
                bool corrected = strongWallDirection &&
                    difference > kWallDirectionCorrectionThreshold;
                if (corrected && difference > kWallDirectionCorrectionMax) {
                    // 稀疏证据防护：与选择角相差那么远的墙峰
                    // 通常是一小段异常墙碎片，
                    // 而不是建筑的主方向。
                    std::cerr << "[DirectionDecision] reject correction: diff_deg="
                              << difference * 180.0 / M_PI
                              << " exceeds cap " << kWallDirectionCorrectionMaxDeg
                              << " (pairs=" << wallPairCount
                              << " peak_ratio=" << wallPeakRatio << ")" << std::endl;
                    corrected = false;
                }
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
                        // 只有足够拉长的轮廓才允许 PCA 回退：
                        // 近各向同性轮廓的 PCA 轴不稳定(±20°噪声)，
                        // 曾把与墙面一致的 MBR 角
                        // 扳偏 27°(拉长比仅 1.37)，毁掉了
                        // 小建筑的规则化。
                        const double gateDeg =
                            pcaAxisRatio >= kPcaDirectionReliableAxisRatio
                                ? kPcaDirectionGateStrongDeg
                                : 1e9;
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
                const double notchMinDepth = std::max(
                    kNotchMinDepthFloor,
                    kNotchMinDepthFactor * model_tuning.repair_distance);
                // 参考缺口来自原始假设而非候选源：
                // 第一个源(修复后的拓扑候选)可能
                // 已经把缺口抹掉，源内检测
                // 将一无所获、保护失效。
                const auto referenceNotches = DetectEvidenceBackedNotches(
                    fallback_hypothesis, fitting_cloud, notchMinDepth, kNotchMaxWidth);
                if (!referenceNotches.empty()) {
                    std::cerr << "[NotchProtect] reference_notches="
                              << referenceNotches.size() << " (widths:";
                    for (const auto& notch : referenceNotches) {
                        std::cerr << " " << notch.width << "x" << notch.depth;
                    }
                    std::cerr << " m)" << std::endl;
                }
                for (const auto& source : single_sources) {
                    // 源自身的缺口驱动正交吸附内部的
                    // 顶点保护掩码。
                    const auto protectedNotches = DetectEvidenceBackedNotches(
                        source, fitting_cloud, notchMinDepth, kNotchMaxWidth);
                    std::vector<bool> notchMask;
                    const std::vector<bool>* notchMaskPtr = nullptr;
                    if (!protectedNotches.empty()) {
                        notchMask = BuildNotchVertexMask(source, protectedNotches);
                        notchMaskPtr = &notchMask;
                    }
                    std::vector<pcl::PointXYZ> candidate;
                    if (!forceOrthogonalPolygonToAngle(
                            source, single_line_angles.front(), candidate,
                            notchMaskPtr)) {
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
                    // 每个参考缺口都必须在候选中保留。
                    if (!referenceNotches.empty()) {
                        int preservedCount = 0;
                        for (const auto& notch : referenceNotches) {
                            if (NotchPreservedInCandidate(candidate, notch)) ++preservedCount;
                        }
                        if (preservedCount < static_cast<int>(referenceNotches.size())) {
                            std::cerr << "[NotchProtect] reject source: kept "
                                      << preservedCount << "/" << referenceNotches.size()
                                      << " reference notches" << std::endl;
                            continue;
                        }
                    }

                    final_points->clear();
                    for (const auto& p : candidate) final_points->push_back(p);
                    best_hypothesis = candidate;
                    single_accepted = true;
                    break;
                }

                auto _single_t1 = std::chrono::steady_clock::now();
                std::cerr << "[SingleFirst] fid=" << source_feature_id_ << " angle_deg="
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
        std::cerr << "[BuildingMode] fid=" << source_feature_id_
            << " " << (allow_diagonal_edges ? "AllowDiagonal" : "StrictOrthogonal")
            << " angles=";
        for (double a : preferred_line_angles) {
            std::cerr << a * 180.0 / M_PI << ",";
        }
        std::cerr << std::endl;

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
                    // 不用原始假设覆盖当前结果——
                    // pre-Ceres ring here.  That ring may contain diagonal
                    // edges even though this building is in StrictOrthogonal
                    // mode.  The common acceptance block below will try a
                    // conservative orthogonal fallback first.
                    std::cerr << "[Orthogonal] strict candidate unavailable; "
                                 "defer fallback decision" << std::endl;
                }
            }
        }

        auto regularizedCandidateOk = [&](const std::vector<pcl::PointXYZ>& candidate) {
            // 小建筑必须严格正交：下面的宽松
            // acceptable branch below would otherwise pass Ceres results
            // Ceres 结果(真实案例：id=2168/
            // 2154 被强制单方向后仍是混合角度)。
            if (forceSingleDirection &&
                !isStrictOrthogonalToMainAngle(
                    candidate,
                    preferred_line_angles.empty() ? dominantLineAngles2D(fallback_hypothesis, 1).front()
                                                  : preferred_line_angles.front(),
                    3.0, model_tuning.fine_prune_distance)) {
                return false;
            }
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
                // 严格模式绝不能以输出斜边原始
                // hypothesis.  A quality-rejected orthogonal candidate is
                // still preferable to violating the mode contract; it is
                // bounded by the same validity/angle checks as normal output.
                bool snappedFallback = false;
                if (!allow_diagonal_edges && !preferred_line_angles.empty()) {
                    std::vector<pcl::PointXYZ> snapped;
                    const double snapAngle = preferred_line_angles.front();
                    if (forceOrthogonalPolygonToAngle(fallback_hypothesis, snapAngle, snapped) &&
                        isSimplePolygon2D(snapped) &&
                        isStrictOrthogonalToMainAngle(
                            snapped, snapAngle, 3.0,
                            model_tuning.fine_prune_distance)) {
                        std::cerr << "[outlineRegular] strict orthogonal fallback" << std::endl;
                        final_points->clear();
                        for (const auto& p : snapped) final_points->push_back(p);
                        snappedFallback = true;
                    }
                }
                if (!snappedFallback) {
                    std::cerr << "[outlineRegular] final polygon invalid, rollback to pre-Ceres hypothesis" << std::endl;
                    final_points->clear();
                    best_hypothesis = fallback_hypothesis;
                    for (const auto& op : best_hypothesis) {
                        final_points->points.push_back(op);
                    }
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
                bool snappedFallback = false;
                if (!allow_diagonal_edges && !preferred_line_angles.empty()) {
                    std::vector<pcl::PointXYZ> snapped;
                    if (forceOrthogonalPolygonToAngle(
                            fallback_hypothesis, preferred_line_angles.front(), snapped) &&
                        isSimplePolygon2D(snapped) &&
                        isStrictOrthogonalToMainAngle(
                            snapped, preferred_line_angles.front(), 3.0,
                            model_tuning.fine_prune_distance)) {
                        std::cerr << "[outlineRegular] out-of-bounds result; "
                                     "using strict orthogonal fallback" << std::endl;
                        final_points->clear();
                        for (const auto& p : snapped) final_points->points.push_back(p);
                        snappedFallback = true;
                    }
                }
                if (!snappedFallback) {
                    std::cerr << "[outlineRegular] 检测到越界点，回退到pre-Ceres hypothesis" << std::endl;
                    final_points->clear();
                    best_hypothesis = fallback_hypothesis;
                    for (auto& op : best_hypothesis) { final_points->points.push_back(op); }
                }
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

    // 定义局部 2D 点结构
    struct Point2D {
        float x, y;
        Point2D(float a = 0.0f, float b = 0.0f) : x(a), y(b) {}
        Point2D operator-(const Point2D& p) const { return { x - p.x, y - p.y }; }
    };

    // 以 lambda 定义局部辅助函数
    auto dot = [](const Point2D& a, const Point2D& b) -> float {
        return a.x * b.x + a.y * b.y;
    };

    auto perp = [](const Point2D& a) -> Point2D {
        return { -a.y, a.x };
    };

    auto norm = [&](const Point2D& a) -> float {
        return std::sqrt(dot(a, a));
    };

    // 直接用 x, y 作 2D 坐标(假设 z=0 或忽略)
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2d(new pcl::PointCloud<pcl::PointXYZ>);
    for (const auto& pt : points) {
        cloud2d->push_back(pcl::PointXYZ(pt.x, pt.y, 0.0f));
    }

    // 计算二维凸包
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

    // 保证逆时针序
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
    const std::vector<double>& preferred_line_angles,
    const std::vector<int>* fixed_edge_assignments,
    bool preserve_topology)
{
    // 初始化当前多边形
    std::vector<pcl::PointXYZ> current_polygon = hypothesis_raw;
    std::vector<pcl::PointXYZ> last_valid_polygon = current_polygon;
    // 拓扑保持模式：顶点数由边链拓扑决定，禁止删点式预处理
    if (!preserve_topology && !allow_diagonal_edges) {
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
        // 拓扑保持模式下短边可能是链拓扑的一部分(窄颈/凹凸)，不剔除
        if (!preserve_topology) {
            pruneShortEdges(current_polygon, current_polygon, tuning.prune_distance);
        }

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

        // When the caller supplies a direction model, it is the complete
        // legal direction set. Re-discovering angles from the already altered
        // polygon reintroduces diagonal edges during Ceres.
        if (base_line_angles.empty()) {
            const std::vector<double> polygon_angles = dominantLineAngles2D(
                current_polygon, allow_diagonal_edges ? 4 : 1);
            for (double angle : polygon_angles) {
                append_angle(angle, allow_diagonal_edges ? 10.0 : 20.0);
            }
        }

        if (!allow_diagonal_edges && base_line_angles.size() > 1) {
            base_line_angles.resize(1);
        }
        if (allow_diagonal_edges && base_line_angles.size() > 4) {
            base_line_angles.resize(4);
        }

        if (use_dlg_direction_ && preferred_line_angles.empty()) {
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
            if (!preferred_line_angles.empty()) break;
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
        // Every edge must belong to one of the supplied direction systems.
        // A multi-direction building may use several axes, but it may not
        // create an unmodelled free-angle edge.
        std::vector<int> edge_base_index(n, -1);

        for (size_t i = 0; i < n; ++i) {
            double best_error = std::numeric_limits<double>::max();
            int best_type = -1;
            int best_base = -1;
            double best_offset = 0.0;

            // 固定边归属：拓扑通道指定的边只允许吸附到指定基方向，
            // 防止重新分类破坏链拓扑的方向结构
            const int fixed_base = (fixed_edge_assignments != nullptr &&
                                    i < fixed_edge_assignments->size())
                ? (*fixed_edge_assignments)[i] : -1;
            const bool has_fixed_base = fixed_base >= 0 &&
                static_cast<size_t>(fixed_base) < base_thetas.size();

            for (size_t k = 0; k < base_thetas.size(); ++k) {
                if (has_fixed_base && k != static_cast<size_t>(fixed_base)) continue;
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
                // This is only possible for a non-finite edge angle or an
                // empty direction set. Keep the edge out of the solver's
                // free-angle path; the caller will reject the result.
                edges[i].type = 0;
                edges[i].theta_offset = 0.0;
            }
        }


        // 5. 数据关联 (带拐角遮蔽优化)
        const bool hasResidualWeightOverride = cloud &&
            residual_weights_override_.size() == cloud->points.size();
        const double minResidualWeight = hasResidualWeightOverride ? 1e-6 : 0.10;
        if (cloud && !cloud->empty()) {
            auto supportBaseWeight = [&](size_t pt_index) {
                // 双残差模式: 逐点权重覆盖(70:30 总权比), 否则构造权重
                const std::vector<double>& w = !hasResidualWeightOverride
                    ? support_weights_ : residual_weights_override_;
                if (pt_index < w.size()) {
                    return clampDouble(w[pt_index], minResidualWeight, 1.0);
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
                    ? clampDouble(edges[i].associated_weights[pt_index], minResidualWeight, 1.0)
                    : 1.0;
                const double support_weight = clampDouble(
                    boundary_weight * normal_weight, minResidualWeight, 1.0);
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
        // 拓扑保持模式下顶点数固定，有效性由下方的简单多边形+IoU检查兜底
        if (!preserve_topology && new_polygon.size() >= 3) {
            std::vector<pcl::PointXYZ> pruned;
            pruneShortEdges(new_polygon, pruned, tuning.fine_prune_distance);
            if (pruned.size() >= 3 && isSimplePolygon2D(pruned)) {
                new_polygon = pruned;
            }
        }

        // 传递给下一次迭代
        // 注: 平行相邻边(台阶 jog)求交 NaN 丢顶点、退化合并是优化器清理
        // 台阶伪影的自然机制, 不强行锁顶点数(锁住会保留中间态斜边);
        // 输出层面的方向一致性由拓扑通道的 rogue 检查兜底。
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
        pcl::PointXYZ new_vertex2; // PARALLEL_STEPʱc_new
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

// ===== TopologyPreservingRegularize =====
// 拓扑保持规则化实验通道：从初始轮廓的连续边链出发（非 VDP 假设），
// 保留显著凹凸/窄颈/方向转折作为受保护顶点，用等弧长轮廓采样点
// 作为稠密几何残差，Ceres 只优化每条链的线参数(θ,d)。
// v2: 修复闭环旋转/原环索引/拟合直线offset/角点锚定/多方向。
enum class TopologyRegularizationStatus {
    SuccessCeres,
    SuccessTopologyCandidate,
    FailedUseNormalPipeline
};

// 备用结果质量检查的公开入口: 组装方向合法角集合后转发共享实现
std::string outlineRegular::CheckRingQuality(
    const std::vector<pcl::PointXYZ>& poly,
    const std::vector<pcl::PointXYZ>& initialRing,
    const DirectionContextOut& dirContext,
    double maxVertexDisp,
    long long fid,
    int partIndex)
{
    // 合法方向只有 systemAngles; 未归组链角度不再进入验收集合。
    // 三级证据: 完整(严格检查全部受约束边) / 部分(只检查主导长边
    // 与 primaryAngle 一致) / 无证据(几何检查+自估方向兜底)。
    const bool multi = dirContext.completeEvidence && dirContext.multiDirection;
    const std::string base = CheckPolygonQualityVsRing(
        poly, initialRing, /*hasDirection=*/dirContext.completeEvidence, multi,
        dirContext.systemAngles, maxVertexDisp, fid, partIndex);
    if (!base.empty()) return base;
    if (!dirContext.completeEvidence && dirContext.valid &&
        !dirContext.systemAngles.empty()) {
        // 部分证据: 最长边(前3条)必须贴近某个已确认系统角(12°)
        struct LongEdge { double len; double ang; };
        std::vector<LongEdge> edges;
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const auto& p1 = poly[i];
            const auto& p2 = poly[(i + 1) % poly.size()];
            edges.push_back({std::hypot(p2.x - p1.x, p2.y - p1.y),
                             std::atan2(p2.y - p1.y, p2.x - p1.x)});
        }
        std::sort(edges.begin(), edges.end(),
                  [](const LongEdge& a, const LongEdge& b) { return a.len > b.len; });
        const std::size_t checkCount = std::min<std::size_t>(3, edges.size());
        for (std::size_t k = 0; k < checkCount; ++k) {
            if (edges[k].len < 3.0) break;
            bool ok = false;
            for (double a : dirContext.systemAngles) {
                if (foldedAngleDistance90(edges[k].ang, a) * 180.0 / M_PI <= 12.0) {
                    ok = true;
                    break;
                }
            }
            if (!ok) return "partial_direction_violation";
        }
    }
    return "";
}

// 兜底局部规则性检查: 方向检查覆盖 0.8~3m 的短斜边, 加上尖刺与
// 近共线锯齿。仅用于 VDP 结果/最优假设/直接输出的拓扑候选——
// 这些结果没有受保护短边的显式标记, 不做任何豁免。
// 保存双残差调试点(局部坐标 + originOffset)到点 Shapefile
bool outlineRegular::SaveRawResidualDebugDump(
    const std::string& shpPath,
    const Eigen::Vector3d& originOffset,
    OGRSpatialReference* spatialRef)
{
    auto& points = RawResidualDebugPoints();
    // 清理旧文件族
    for (const char* ext : {".shp", ".shx", ".dbf", ".prj", ".cpg"}) {
        std::filesystem::path p = shpPath;
        p.replace_extension(ext);
        std::error_code ec;
        std::filesystem::remove(p, ec);
    }
    if (points.empty()) return false;
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    if (!driver) return false;
    GDALDataset* ds = driver->Create(
        shpPath.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!ds) return false;
    OGRLayer* layer = ds->CreateLayer(
        "debug_mask_raw_residual_points", spatialRef, wkbPoint25D, nullptr);
    if (!layer) { GDALClose(ds); return false; }
    OGRFieldDefn fFid("fid", OFTInteger64);
    OGRFieldDefn fPart("part", OFTInteger);
    OGRFieldDefn fSrc("source", OFTInteger);  // 0=smooth 1=raw
    OGRFieldDefn fW("weight", OFTReal);
    OGRFieldDefn fAssoc("assoc_edge", OFTInteger);     // -1=未关联
    layer->CreateField(&fFid);
    layer->CreateField(&fPart);
    layer->CreateField(&fSrc);
    layer->CreateField(&fW);
    layer->CreateField(&fAssoc);
    for (const auto& p : points) {
        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField(0, p.fid);
        f->SetField(1, p.part);
        f->SetField(2, p.src);
        f->SetField(3, p.weight);
        f->SetField(4, p.assocEdge);
        OGRPoint pt(p.pt.x + originOffset.x(),
                    p.pt.y + originOffset.y(),
                    p.pt.z + originOffset.z());
        f->SetGeometry(&pt);
        layer->CreateFeature(f);
        OGRFeature::DestroyFeature(f);
    }
    ds->FlushCache();
    GDALClose(ds);
    points.clear();
    return true;
}

std::string outlineRegular::CheckFallbackLocalRegularity(
    const std::vector<pcl::PointXYZ>& poly,
    const DirectionContextOut& dirContext,
    long long fid,
    int partIndex,
    const char* stage)
{
    if (poly.size() < 3) return "too_few_vertices";
    const std::size_t m = poly.size();
    int shortDiagonal = 0, spike = 0, zigzag = 0;
    double irregularLength = 0.0, perimeter = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        const auto& p1 = poly[i];
        const auto& p2 = poly[(i + 1) % m];
        perimeter += std::hypot(p2.x - p1.x, p2.y - p1.y);
    }
    // 1) 短斜边: 0.8~3m 的边做方向检查(有方向证据时)
    if (dirContext.valid && !dirContext.systemAngles.empty()) {
        for (std::size_t i = 0; i < m; ++i) {
            const auto& p1 = poly[i];
            const auto& p2 = poly[(i + 1) % m];
            const double len = std::hypot(p2.x - p1.x, p2.y - p1.y);
            if (len < kDirDiagShortDiagMinLen || len >= 3.0) continue;
            const double ang = std::atan2(p2.y - p1.y, p2.x - p1.x);
            bool ok = false;
            for (double a : dirContext.systemAngles) {
                if (foldedAngleDistance90(ang, a) * 180.0 / M_PI <= 12.0) {
                    ok = true;
                    break;
                }
            }
            if (!ok) {
                ++shortDiagonal;
                irregularLength += len;
            }
        }
    }
    // 2) 尖刺与 3) 近共线锯齿: 逐顶点检查
    for (std::size_t i = 0; i < m; ++i) {
        const auto& prev = poly[(i + m - 1) % m];
        const auto& cur = poly[i];
        const auto& next = poly[(i + 1) % m];
        const double v1x = cur.x - prev.x, v1y = cur.y - prev.y;
        const double v2x = next.x - cur.x, v2y = next.y - cur.y;
        const double l1 = std::hypot(v1x, v1y);
        const double l2 = std::hypot(v2x, v2y);
        if (l1 < 1e-9 || l2 < 1e-9) continue;
        const double cosA = (v1x * v2x + v1y * v2y) / (l1 * l2);
        const double angleDeg =
            std::acos(std::clamp(cosA, -1.0, 1.0)) * 180.0 / M_PI;
        // 尖刺: 内角 <32° 且 三角面积小(删除该点几何损失可忽略)
        // 或两侧都是短边(细尖残留)
        if (angleDeg < kDirDiagSpikeAngleDeg) {
            const double triArea = 0.5 * std::abs(
                (cur.x - prev.x) * (next.y - prev.y) -
                (next.x - prev.x) * (cur.y - prev.y));
            if (triArea < kDirDiagSpikeMaxArea ||
                (l1 < 2.0 && l2 < 2.0)) {
                ++spike;
                irregularLength += std::min(l1, l2);
            }
        }
        // 近共线锯齿: 顶点紧贴前后连线但转角明显(>5°)
        {
            const double chordDx = next.x - prev.x;
            const double chordDy = next.y - prev.y;
            const double chordLen = std::hypot(chordDx, chordDy);
            if (chordLen > 1e-9) {
                const double dev = std::abs(
                    (cur.x - prev.x) * chordDy - (cur.y - prev.y) * chordDx) /
                    chordLen;
                const double turnDeg = 180.0 - angleDeg;
                if (dev < kDirDiagZigzagChordDev && turnDeg > 5.0 && turnDeg < 45.0) {
                    ++zigzag;
                    irregularLength += std::min(l1, l2);
                }
            }
        }
    }
    const double irregularRatio =
        perimeter > 1e-9 ? irregularLength / perimeter : 0.0;
    // Direction evidence is a hard contract: no visible straight edge may
    // remain outside the accepted direction systems. Curves are disabled in
    // mask-only mode, so there is no curve exemption here.
    const bool accepted =
        (shortDiagonal == 0 && irregularRatio <= 0.10) &&
        spike == 0 && zigzag < 3;
    std::string reason;
    if (!accepted) {
        if (shortDiagonal > 0) reason = "short_diagonal";
        else if (spike > 0) reason = "spike";
        else reason = "zigzag";
    }
    if (fid >= 0) {
        std::cerr << "[FallbackQuality] fid=" << fid << " part=" << partIndex
                  << " stage=" << (stage ? stage : "?")
                  << " short_diagonal=" << shortDiagonal
                  << " spike=" << spike
                  << " zigzag=" << zigzag
                  << " irregular_ratio=" << irregularRatio
                  << " accepted=" << (accepted ? 1 : 0)
                  << " reason=" << (accepted ? "-" : reason)
                  << std::endl;
    }
    return accepted ? "" : ("fallback_" + reason);
}

// 保守清理: 只删除有明确低面积证据的尖刺和近共线锯齿顶点。
// 真实矩形凹凸的角是 90°(不会被尖刺判据命中), 凹凸底部顶点
// 距弦 0.3~0.5m(不会被锯齿判据命中)——只清确定性噪声。
std::vector<pcl::PointXYZ> outlineRegular::CleanLowEvidenceIrregularities(
    const std::vector<pcl::PointXYZ>& poly)
{
    std::vector<pcl::PointXYZ> current = poly;
    for (int pass = 0; pass < 3; ++pass) {
        if (current.size() < 5) break;
        std::vector<pcl::PointXYZ> next;
        next.reserve(current.size());
        const std::size_t m = current.size();
        std::vector<bool> drop(m, false);
        for (std::size_t i = 0; i < m; ++i) {
            const auto& prev = current[(i + m - 1) % m];
            const auto& cur = current[i];
            const auto& nx = current[(i + 1) % m];
            const double v1x = cur.x - prev.x, v1y = cur.y - prev.y;
            const double v2x = nx.x - cur.x, v2y = nx.y - cur.y;
            const double l1 = std::hypot(v1x, v1y);
            const double l2 = std::hypot(v2x, v2y);
            if (l1 < 1e-9 || l2 < 1e-9) { drop[i] = true; continue; }
            const double cosA = (v1x * v2x + v1y * v2y) / (l1 * l2);
            const double angleDeg =
                std::acos(std::clamp(cosA, -1.0, 1.0)) * 180.0 / M_PI;
            if (angleDeg < kDirDiagSpikeAngleDeg) {
                const double triArea = 0.5 * std::abs(
                    (cur.x - prev.x) * (nx.y - prev.y) -
                    (nx.x - prev.x) * (cur.y - prev.y));
                if (triArea < kDirDiagSpikeMaxArea ||
                    (l1 < 2.0 && l2 < 2.0)) {
                    drop[i] = true;
                    continue;
                }
            }
            const double chordDx = nx.x - prev.x;
            const double chordDy = nx.y - prev.y;
            const double chordLen = std::hypot(chordDx, chordDy);
            if (chordLen > 1e-9) {
                const double dev = std::abs(
                    (cur.x - prev.x) * chordDy - (cur.y - prev.y) * chordDx) /
                    chordLen;
                const double turnDeg = 180.0 - angleDeg;
                if (dev < kDirDiagZigzagChordDev && turnDeg > 5.0 && turnDeg < 45.0) {
                    drop[i] = true;
                }
            }
        }
        bool changed = false;
        for (std::size_t i = 0; i < m; ++i) {
            if (drop[i]) changed = true;
            else next.push_back(current[i]);
        }
        if (!changed || next.size() < 4) break;
        // 清理后自交/退化检查, 失败则放弃本轮结果
        if (!isSimplePolygon2D(next)) break;
        current.swap(next);
    }
    return current;
}

std::vector<pcl::PointXYZ> outlineRegular::TopologyPreservingRegularize(
    const std::vector<pcl::PointXYZ>& initialRing,
    double pixelSize,
    bool& usedFallback,
    DirectionContextOut* dirContext,
    int partIndex,
    const std::vector<pcl::PointXYZ>* rawRing)
{
    usedFallback = false;
    if (initialRing.size() < 6) {
        usedFallback = true;
        return {};
    }

    const double area = polygonArea2D(initialRing);
    if (area < 5.0) {
        usedFallback = true;
        return {};
    }

    // ---- 0.5 圆形/曲线轮廓预检 ----
    // 等周比接近 1 的圆形轮廓在正交格网吸附下会产生劣质多边形;
    // regular_Contour 内置圆/椭圆检测, 此类轮廓直接交给 VDP。
    // 方形建筑等周比≈0.785, 一般建筑 <=0.85, 阈值 0.90 只拦真圆。
    {
        double perim = 0.0;
        for (std::size_t i = 0; i < initialRing.size(); ++i) {
            perim += std::hypot(
                initialRing[(i + 1) % initialRing.size()].x - initialRing[i].x,
                initialRing[(i + 1) % initialRing.size()].y - initialRing[i].y);
        }
        const double roundness =
            4.0 * M_PI * std::abs(area) / std::max(perim * perim, 1e-9);
        if (roundness >= 0.90) {
            usedFallback = true;
            std::cerr << "[TopologyFallbackToVDP] fid=" << source_feature_id_
                      << " reason=contour_circular roundness=" << roundness
                      << std::endl;
            return {};
        }
    }

    // ---- 1. 边链提取(split-and-merge, 修复版) ----
    // 拓扑通道用较宽的拟合容差: 长墙上的栅格锯齿(±0.3~0.5m)应直接被
    // 吸收进长链，只有显著凹凸(>容差)才成为独立拓扑边
    const double fitTolerance = std::clamp(
        std::max(0.50, 1.2 * pixelSize), 0.15, 0.90);
    auto allChains = ExtractDirectionChains(initialRing, fitTolerance);
    if (allChains.size() < 3) {
        usedFallback = true;
        return {};
    }

    // ---- 1.5 楼梯run吸收 ----
    // 栅格化斜边在平滑后会残留"连续短链"楼梯(交替 H/V 或小幅折转)。
    // 极大短链run若整体接近一条直线(maxDev 小)，则整段拟合成单条链——
    // 既恢复真实斜边，也压缩长直边上的锯齿。真实凹凸(直线上大偏差)不会触发。
    std::size_t absorbedRunCount = 0;
    {
        const double stepMaxLen = std::max(1.2, 2.5 * pixelSize);
        const double runMinSpan = std::max(2.5, 5.0 * pixelSize);
        const double runFitTol = std::max(0.35, 1.2 * pixelSize);
        const std::size_t N = initialRing.size();
        const std::size_t M = allChains.size();
        std::vector<bool> isStep(M);
        for (std::size_t c = 0; c < M; ++c) {
            isStep[c] = allChains[c].length <= stepMaxLen;
        }
        std::vector<DirectionChain> workChains;
        std::size_t i = 0;
        while (i < M) {
            if (!isStep[i]) { workChains.push_back(allChains[i]); ++i; continue; }
            std::size_t j = i;
            while (j < M && isStep[j]) ++j;
            double span = 0.0;
            for (std::size_t c = i; c < j; ++c) span += allChains[c].length;
            bool absorbed = false;
            if (j - i >= 3 && span >= runMinSpan) {
                // 收集run覆盖的环顶点(原环索引，含环回)
                const std::size_t s = allChains[i].startIndexOriginal;
                const std::size_t e = allChains[j - 1].endIndexOriginal;
                std::vector<pcl::PointXYZ> pts;
                for (std::size_t k = s;; k = (k + 1) % N) {
                    pts.push_back(initialRing[k]);
                    if (k == e) break;
                }
                if (pts.size() >= 3 && pts.size() <= N) {
                    double dirX, dirY, rmse, maxDev, cx2, cy2, nx, ny, off;
                    if (FitChainLine(pts, 0, pts.size() - 1, dirX, dirY, rmse, maxDev,
                                     &cx2, &cy2, &nx, &ny, &off) &&
                        maxDev <= runFitTol) {
                        DirectionChain m = allChains[i];
                        m.endIndexOriginal = e;
                        m.endPoint = allChains[j - 1].endPoint;
                        m.length = span;
                        m.centerX = 0.5 * (m.startPoint.x + m.endPoint.x);
                        m.centerY = 0.5 * (m.startPoint.y + m.endPoint.y);
                        m.angleRad = foldedLineAngle90(std::atan2(dirY, dirX));
                        m.rmse = rmse;
                        m.maxDeviation = maxDev;
                        m.normalX = nx;
                        m.normalY = ny;
                        m.lineOffset = off;
                        m.isShort = false;
                        const double fitScore = std::exp(
                            -(rmse / std::max(fitTolerance, 0.01)) *
                            (rmse / std::max(fitTolerance, 0.01)));
                        m.weight = span * fitScore;
                        workChains.push_back(m);
                        absorbed = true;
                        ++absorbedRunCount;
                    }
                }
            }
            if (!absorbed) {
                for (std::size_t c = i; c < j; ++c) workChains.push_back(allChains[c]);
            }
            i = j;
        }
        allChains = std::move(workChains);
    }
    if (allChains.size() < 3) {
        usedFallback = true;
        return {};
    }

    // 方向投票只用长链(>=kDirectionMinChainLength)
    std::vector<DirectionChain> directionChains;
    for (const auto& c : allChains) {
        if (!c.isShort) directionChains.push_back(c);
    }

    // ---- 2. 方向系统判定(两阶段聚类+加权KDE, 单次调用统一用于
    //      方向决策/格网吸附/DirectionContextOut/日志) ----
    // mask-only 无 OSGB 墙面证据(supportCloud 传空); 判定综合:
    // 支持链总长度(权重占比)/拟合残差(conf 内 fitFactor)/
    // 模型改善幅度(评分增益)/空间范围(systems.extent)。
    const auto dirBuild = BuildDirectionSystems(allChains, fitTolerance, nullptr);
    const bool isMultiDirection = dirBuild.directionCertain && dirBuild.multiDirection;

    // 统一诊断/应用日志(与实际使用的 DirectionSystemBuild 同源,
    // 同一环的 verdict/分数/系统角必然一致; part 为环序号)
    {
        const char* mode = !dirBuild.directionCertain
            ? "uncertain" : (isMultiDirection ? "multi" : "single");
        std::cerr << "[DirectionApply] fid=" << source_feature_id_
                  << " part=" << partIndex << " mode=" << mode << std::endl;
        std::cerr << "[DirectionApply] fid=" << source_feature_id_
                  << " part=" << partIndex
                  << " systems=" << dirBuild.systems.size() << std::endl;
        for (std::size_t s = 0; s < dirBuild.systems.size(); ++s) {
            const auto& sys = dirBuild.systems[s];
            std::cerr << "[DirectionSystem] fid=" << source_feature_id_
                      << " part=" << partIndex
                      << " index=" << s
                      << " angle_deg=" << sys.angleRad * 180.0 / M_PI
                      << " chains=" << sys.chainCount
                      << " total_length=" << sys.totalLength
                      << " concentration=" << sys.concentration
                      << " confidence=" << sys.confidence << std::endl;
        }
        int stableCount = 0, assignedCount = 0;
        for (const auto& info : dirBuild.chainInfo) {
            if (info.stable) ++stableCount;
            if (info.stable && info.system >= 0) ++assignedCount;
        }
        std::cerr << "[DirectionAssignment] fid=" << source_feature_id_
                  << " part=" << partIndex
                  << " stable=" << stableCount
                  << " assigned=" << assignedCount
                  << " unassigned_long=" << dirBuild.unassignedLongCount
                  << " unassigned_len=" << dirBuild.unassignedLongLength
                  << " unassigned_ratio=" << dirBuild.unassignedLengthRatio
                  << " rescued=" << dirBuild.rescuedLongChains
                  << " protected_short=" << dirBuild.shortProtected << std::endl;
        std::cerr << "[DirectionTotals] fid=" << source_feature_id_
                  << " part=" << partIndex
                  << " stable_length=" << dirBuild.totalStableLength
                  << " stable_weight=" << dirBuild.totalStableWeight << std::endl;
        std::cerr << "[DirectionModel] fid=" << source_feature_id_
                  << " part=" << partIndex
                  << " single_score=" << dirBuild.singleScore
                  << " multi_score=" << dirBuild.multiScore
                  << " gain=" << (dirBuild.multiScore - dirBuild.singleScore)
                  << std::endl;
    }

    // 方向上下文传出: 无论最终是否回退都先填写——部分可信方向
    // 也约束 VDP/最优假设的主导长边, 完全无证据才允许自估兜底
    if (dirContext) {
        dirContext->valid = !dirBuild.systems.empty();
        dirContext->completeEvidence = dirBuild.directionCertain;
        dirContext->multiDirection = isMultiDirection;
        dirContext->primaryAngle = dirBuild.systems.empty()
            ? 0.0 : dirBuild.systems.front().angleRad;
        dirContext->systemAngles.clear();
        for (const auto& s : dirBuild.systems) {
            dirContext->systemAngles.push_back(s.angleRad);
        }
        dirContext->unassignedLengthRatio = dirBuild.unassignedLengthRatio;
    }

    if (!dirBuild.directionCertain) {
        // 方向不确定: 交给 VDP 备用流程, 不强行套用方向系统。
        // 区分原因: 未归组占比超限(方向证据冲突) vs 置信度/评分不足
        usedFallback = true;
        std::cerr << "[TopologyFallbackToVDP] fid=" << source_feature_id_
                  << " part=" << partIndex
                  << " reason="
                  << (dirBuild.uncertainReason.empty()
                          ? "direction_uncertain" : dirBuild.uncertainReason)
                  << std::endl;
        return {};
    }

    // ---- 2.5 按方向系统格网吸附 + 相邻共线链合并 ----
    // (1) 格网吸附——每条链吸附到其归属系统的格网(θ 或 θ+90°):
    //     单方向建筑: 所有链(含短链/未归组链)强制主系统格网, 禁止斜边;
    //     多方向建筑: 各系统链保持各自格网角度; 未归组的稳定链保持
    //     自由方向(不创建新系统); 受保护结构短链保留自身方向;
    // (2) 共线合并——相邻链格网法向平行且偏移差小于阈值时合并为一条物理边。
    // 真实转角(不同系统或偏移差大)不会被合并，凹凸/窄颈拓扑得以保留。
    std::vector<DirectionChain> topoChains = allChains;
    {
        const double offsetMergeTol = std::max(0.45, 1.0 * pixelSize);
        std::vector<int> chainSystem(topoChains.size(), -1);
        for (std::size_t i = 0; i < topoChains.size(); ++i) {
            chainSystem[i] = dirBuild.chainInfo[i].system;
        }
        // (1) 格网吸附: 重写每条吸附链的直线参数
        for (std::size_t i = 0; i < topoChains.size(); ++i) {
            int sys = chainSystem[i];
            if (!isMultiDirection) {
                // 单方向: 一律覆盖到主系统格网。已归组到次系统的链
                // (次系统未达多方向门槛但仍聚了类)若保留其格网,
                // 候选会同时出现两个相差十几度的格网
                sys = 0;
                chainSystem[i] = 0;
            } else if (sys < 0) {
                // A protected/unassigned chain is still a building boundary;
                // assign it to the nearest credible system instead of
                // carrying a free-angle edge into the candidate.
                double bestDist = std::numeric_limits<double>::max();
                int bestSys = -1;
                for (std::size_t s = 0; s < dirBuild.systems.size(); ++s) {
                    const double d = foldedAngleDistance90(
                        topoChains[i].angleRad,
                        dirBuild.systems[s].angleRad);
                    if (d < bestDist) {
                        bestDist = d;
                        bestSys = static_cast<int>(s);
                    }
                }
                if (bestSys >= 0) {
                    sys = bestSys;
                    chainSystem[i] = bestSys;
                }
            }
            if (sys < 0) continue;
            const double pa = dirBuild.systems[static_cast<std::size_t>(sys)].angleRad;
            const double n1x = -std::sin(pa), n1y = std::cos(pa);
            const double n2x = std::cos(pa), n2y = std::sin(pa);
            DirectionChain& c = topoChains[i];
            const double d1 = std::abs(n1x * c.normalX + n1y * c.normalY);
            const double d2 = std::abs(n2x * c.normalX + n2y * c.normalY);
            double gx = (d1 >= d2) ? n1x : n2x;
            double gy = (d1 >= d2) ? n1y : n2y;
            if (gx * c.normalX + gy * c.normalY < 0.0) { gx = -gx; gy = -gy; }
            c.normalX = gx;
            c.normalY = gy;
            c.lineOffset = gx * c.centerX + gy * c.centerY;
        }
        // (2) 共线合并: 仅当两链都已吸附(同格网)且法向平行、偏移差小
        bool mergedAny = true;
        while (mergedAny && topoChains.size() > 3) {
            mergedAny = false;
            for (std::size_t i = 0; i < topoChains.size(); ++i) {
                const std::size_t j = (i + 1) % topoChains.size();
                const DirectionChain& ci = topoChains[i];
                const DirectionChain& cj = topoChains[j];
                if (chainSystem[i] < 0 || chainSystem[j] < 0) continue;
                // 格网法向必须平行(同一物理方向)，垂直的不合并
                if (std::abs(ci.normalX * cj.normalX + ci.normalY * cj.normalY) < 0.9) continue;
                // j 的偏移按 i 的法向符号对齐后比较
                double onx = ci.normalX, ony = ci.normalY;
                if (onx * cj.normalX + ony * cj.normalY < 0.0) { onx = -onx; ony = -ony; }
                const double offJ = onx * cj.centerX + ony * cj.centerY;
                if (std::abs(ci.lineOffset - offJ) > offsetMergeTol) continue;
                // 合并: 保留 i 的区间端点，几何参数取长度加权
                DirectionChain m = ci;
                const double w1 = ci.length;
                const double w2 = cj.length;
                const double wSum = std::max(w1 + w2, 1e-9);
                m.endPoint = cj.endPoint;
                m.length = w1 + w2;
                m.centerX = 0.5 * (m.startPoint.x + m.endPoint.x);
                m.centerY = 0.5 * (m.startPoint.y + m.endPoint.y);
                m.lineOffset = (ci.lineOffset * w1 + offJ * w2) / wSum;
                m.isShort = m.length < kDirectionMinChainLength;
                topoChains[i] = m;
                topoChains.erase(topoChains.begin() + j);
                chainSystem.erase(chainSystem.begin() + j);
                mergedAny = true;
                break;
            }
        }
    }
    if (topoChains.size() < 3) {
        usedFallback = true;
        std::cerr << "[TopologyFallbackToVDP] fid=" << source_feature_id_
                  << " reason=merged_chains_too_few" << std::endl;
        return {};
    }

    std::cerr << "[TopologyChain] fid=" << source_feature_id_
              << " raw_vertices=" << initialRing.size()
              << " all_chains=" << allChains.size()
              << " absorbed_runs=" << absorbedRunCount
              << " merged_topo_chains=" << topoChains.size()
              << " direction_chains=" << directionChains.size()
              << std::endl;

    // ---- 3. 构造候选拓扑 ----
    // 顶点来源优先级:
    //   (1) 有限边段交点: 两链格网线交点, 须同时满足——距原始角点<=1m,
    //       且交点在两链线段方向上的投影参数落在有限范围内
    //       (端点外延不超过 max(1m, 链长50%));
    //   (2) 原始角点投影对: 近平行链/交点过远时, 不无限延长求交,
    //       原始角点分别投影并夹取到两链线段上, 插入一对顶点——
    //       台阶 jog/凹凸/窄颈的拓扑结构因此保留。
    std::vector<pcl::PointXYZ> topologyPolygon;
    topologyPolygon.reserve(topoChains.size());
    const float zVal = initialRing[0].z;
    // 点到直线(法向式 n·x = off)的投影
    auto projectOntoLine = [](const pcl::PointXYZ& p,
                              double nx, double ny, double off) {
        const double t = nx * p.x + ny * p.y - off;
        pcl::PointXYZ q;
        q.x = static_cast<float>(p.x - t * nx);
        q.y = static_cast<float>(p.y - t * ny);
        return q;
    };
    // 投影并夹取到链的有限范围: 沿吸附线的方向参数化
    // (用原始弦方向夹取会把投影点拉离格网线——吸附旋转后的线方向
    // 与弦方向可差 25°, 连接边随之偏斜; 必须沿吸附线本身夹取)
    auto projectClampedToChain = [&](const pcl::PointXYZ& p,
                                     const DirectionChain& c) {
        // 吸附线: n·x = off, 方向 d = (-ny, nx)
        const double dx = -c.normalY, dy = c.normalX;
        // 线上坐标: t = d·x (线本身过 n*off 点)
        auto lineCoord = [&](const pcl::PointXYZ& q) {
            return dx * q.x + dy * q.y;
        };
        const double tSp = lineCoord(c.startPoint);
        const double tEp = lineCoord(c.endPoint);
        const double tMin = std::min(tSp, tEp);
        const double tMax = std::max(tSp, tEp);
        pcl::PointXYZ q = projectOntoLine(p, c.normalX, c.normalY, c.lineOffset);
        const double t = std::clamp(lineCoord(q), tMin, tMax);
        // 沿吸附线重建: 基点 n*off + 方向 d * t
        q.x = static_cast<float>(c.normalX * c.lineOffset + dx * t);
        q.y = static_cast<float>(c.normalY * c.lineOffset + dy * t);
        return q;
    };
    // 点在链线段方向上的参数(0=起点, len=终点)
    auto segParam = [&](const pcl::PointXYZ& p, const DirectionChain& c) {
        const double vx = c.endPoint.x - c.startPoint.x;
        const double vy = c.endPoint.y - c.startPoint.y;
        const double len = std::hypot(vx, vy);
        if (len < 1e-9) return 0.0;
        return ((p.x - c.startPoint.x) * vx + (p.y - c.startPoint.y) * vy) / len;
    };

    for (std::size_t c = 0; c < topoChains.size(); ++c) {
        const auto& chainA = topoChains[c];
        const auto& chainB = topoChains[(c + 1) % topoChains.size()];
        // 原始角点 = 链A终点 = 链B起点(共享环顶点)
        const pcl::PointXYZ rawCorner = chainA.endPoint;

        pcl::PointXYZ vertex = rawCorner;
        vertex.z = zVal;
        bool useIntersection = false;
        {
            const double det = chainA.normalX * chainB.normalY -
                               chainB.normalX * chainA.normalY;
            if (std::abs(det) > 1e-9) {
                pcl::PointXYZ isect;
                isect.x = static_cast<float>(
                    (chainA.lineOffset * chainB.normalY -
                     chainA.normalY * chainB.lineOffset) / det);
                isect.y = static_cast<float>(
                    (chainA.normalX * chainB.lineOffset -
                     chainA.lineOffset * chainB.normalX) / det);
                isect.z = zVal;
                const double shift = std::hypot(
                    isect.x - rawCorner.x, isect.y - rawCorner.y);
                const double lenA = std::hypot(
                    chainA.endPoint.x - chainA.startPoint.x,
                    chainA.endPoint.y - chainA.startPoint.y);
                const double lenB = std::hypot(
                    chainB.endPoint.x - chainB.startPoint.x,
                    chainB.endPoint.y - chainB.startPoint.y);
                const double marginA = std::max(1.0, 0.5 * lenA);
                const double marginB = std::max(1.0, 0.5 * lenB);
                const double tA = segParam(isect, chainA);
                const double tB = segParam(isect, chainB);
                // 有限边段交点的全部条件: 数值有效 + 近角点 + 参数在
                // 两链有限范围内(含适度外延, 允许墙线自然延伸到角点)
                if (std::isfinite(isect.x) && std::isfinite(isect.y) &&
                    shift <= 1.0 &&
                    tA >= -marginA && tA <= lenA + marginA &&
                    tB >= -marginB && tB <= lenB + marginB) {
                    vertex = isect;
                    useIntersection = true;
                }
            }
        }
        if (useIntersection) {
            topologyPolygon.push_back(vertex);
        } else {
            // 近平行链或交点越界: 原始角点投影对(夹取到链线段),
            // 连接边短且贴近格网, 不产生远距离飞点
            pcl::PointXYZ pa = projectClampedToChain(rawCorner, chainA);
            pcl::PointXYZ pb = projectClampedToChain(rawCorner, chainB);
            pa.z = zVal;
            pb.z = zVal;
            topologyPolygon.push_back(pa);
            topologyPolygon.push_back(pb);
        }
    }

    std::vector<double> systemAngles;
    for (const auto& s : dirBuild.systems) systemAngles.push_back(s.angleRad);

    // Remove duplicate vertices without losing the direction context. The
    // candidate is small, so rebuild the assignment from the resulting edge
    // angles; every edge is then fixed to one supplied system before Ceres.
    removeDuplicatePoints2D(topologyPolygon, 0.05f);
    std::vector<int> fixedEdgeAssignments(topologyPolygon.size(), -1);
    for (std::size_t i = 0; i < topologyPolygon.size(); ++i) {
        const auto& a = topologyPolygon[i];
        const auto& b = topologyPolygon[(i + 1) % topologyPolygon.size()];
        const double edgeAngle = std::atan2(b.y - a.y, b.x - a.x);
        double best = std::numeric_limits<double>::max();
        int bestSystem = -1;
        for (std::size_t s = 0; s < systemAngles.size(); ++s) {
            const double d = foldedAngleDistance90(
                edgeAngle, systemAngles[s]);
            if (d < best) {
                best = d;
                bestSystem = static_cast<int>(s);
            }
        }
        fixedEdgeAssignments[i] = bestSystem;
    }

    // ---- 3.5 质量检查器(转发到共享实现, 候选与 Ceres 结果共用; ----
    //      VDP 备用结果经 CheckRingQuality 走同一标准)
    // 合法方向只有系统角: 单方向验收基准=主系统角(8°),
    // 多方向=全部系统角(12°); 未归组链角度不再进入验收集合
    auto qualityCheck = [&](const std::vector<pcl::PointXYZ>& poly,
                            double maxVertexDisp) -> std::string {
        return CheckPolygonQualityVsRing(
            poly, initialRing, /*hasDirection=*/true, isMultiDirection,
            systemAngles, maxVertexDisp,
            source_feature_id_, partIndex);
    };

    // ---- 4. 候选拓扑质量验证(不合格→VDP, 不直接输出) ----
    const std::string candidateReason = qualityCheck(topologyPolygon, 2.5);
    if (!candidateReason.empty()) {
        usedFallback = true;
        std::cerr << "[TopologyCandidateReject] fid=" << source_feature_id_
                  << " reason=" << candidateReason << std::endl;
        std::cerr << "[TopologyFallbackToVDP] fid=" << source_feature_id_
                  << " reason=candidate_quality" << std::endl;
        return {};
    }

    // raw 采样点到 smooth 环边界的距离(门控用)
    auto distToRingPt = [&](const pcl::PointXYZ& pt) {
        double best = std::numeric_limits<double>::max();
        for (std::size_t i = 0; i < initialRing.size(); ++i) {
            const auto& a = initialRing[i];
            const auto& b = initialRing[(i + 1) % initialRing.size()];
            const double dx = b.x - a.x;
            const double dy = b.y - a.y;
            const double lenSq = dx * dx + dy * dy;
            if (lenSq < 1e-12) continue;
            double t2 = ((pt.x - a.x) * dx + (pt.y - a.y) * dy) / lenSq;
            t2 = std::clamp(t2, 0.0, 1.0);
            const double d = std::hypot(pt.x - a.x - t2 * dx, pt.y - a.y - t2 * dy);
            best = std::min(best, d);
        }
        return best;
    };
    // ---- 5. 残差采样: 双轮廓(smooth 定拓扑 + raw 像素轮廓提供几何观测) ----
    // smooth 采样: 拓扑/方向来源自身的等弧长重采样(拐角噪声由求解器
    //   数据关联的中段遮蔽处理);
    // raw 采样(可选): 平滑前原始像素轮廓的等弧长采样——真实边界位置,
    //   通过可靠性门控后才参与; raw 不进入方向检测, 不影响候选顶点。
    pcl::PointCloud<pcl::PointXYZ>::Ptr residualCloud(
        new pcl::PointCloud<pcl::PointXYZ>);
    std::vector<double> residualWeights;
    std::vector<pcl::PointXYZ> smoothSamples;
    std::vector<pcl::PointXYZ> rawSamples;
    std::vector<int> smoothAssociations;
    std::vector<int> rawAssociations;
    bool rawAccepted = false;
    double rawAreaCached = 0.0;
    std::string rawRejectReason;

    auto sampleBoundary = [&](const std::vector<pcl::PointXYZ>& ring,
                              double step,
                              std::size_t maxPoints) {
        std::vector<pcl::PointXYZ> samples;
        for (std::size_t i = 0; i < ring.size(); ++i) {
            const auto& a = ring[i];
            const auto& b = ring[(i + 1) % ring.size()];
            const double len = std::hypot(b.x - a.x, b.y - a.y);
            const int count = std::max(1, static_cast<int>(std::ceil(len / step)));
            for (int k = 0; k < count; ++k) {
                const double t = static_cast<double>(k) / count;
                pcl::PointXYZ pt;
                pt.x = static_cast<float>(a.x + t * (b.x - a.x));
                pt.y = static_cast<float>(a.y + t * (b.y - a.y));
                pt.z = zVal;
                samples.push_back(pt);
            }
        }
        if (maxPoints > 0 && samples.size() > maxPoints) {
            std::vector<pcl::PointXYZ> reduced;
            reduced.reserve(maxPoints);
            for (std::size_t i = 0; i < maxPoints; ++i) {
                reduced.push_back(samples[
                    i * (samples.size() - 1) / (maxPoints - 1)]);
            }
            samples.swap(reduced);
        }
        return samples;
    };

    const double smoothStep = std::clamp(pixelSize, 0.3, 0.8);
    smoothSamples = sampleBoundary(initialRing, smoothStep, 0);
    const double assocMaxDist = std::max(3.0 * pixelSize, 1.2);
    auto associateEdge = [&](const pcl::PointXYZ& pt) {
        int bestEdge = -1;
        double bestDist = assocMaxDist;
        for (std::size_t e = 0; e < topologyPolygon.size(); ++e) {
            const auto& p1 = topologyPolygon[e];
            const auto& p2 = topologyPolygon[(e + 1) % topologyPolygon.size()];
            const double dx = p2.x - p1.x;
            const double dy = p2.y - p1.y;
            const double lenSq = dx * dx + dy * dy;
            if (lenSq < 1e-12) continue;
            const double t = ((pt.x - p1.x) * dx +
                              (pt.y - p1.y) * dy) / lenSq;
            if (t <= 0.05 || t >= 0.95) continue;
            const double d = std::hypot(
                pt.x - p1.x - t * dx, pt.y - p1.y - t * dy);
            if (d < bestDist) {
                bestDist = d;
                bestEdge = static_cast<int>(e);
            }
        }
        return bestEdge;
    };

    // raw 必须与当前 smooth 部件面积和位置一致。raw 点先完成有限边段
    // 关联，未关联点不会进入 Ceres。
    if (kUseRawMaskBoundaryResiduals && rawRing && rawRing->size() >= 3) {
        rawAreaCached = std::abs(polygonArea2D(*rawRing));
        const double smoothArea = std::abs(polygonArea2D(initialRing));
        const double areaRatio = rawAreaCached / std::max(smoothArea, 1e-6);
        if (areaRatio < 0.60 || areaRatio > 1.60) {
            rawRejectReason = "area_ratio=" + std::to_string(areaRatio);
        } else {
            const double rawStep = std::clamp(
                std::max(pixelSize, kRawResidualSpacingMin),
                kRawResidualSpacingMin, kRawResidualSpacingMax);
            rawSamples = sampleBoundary(
                *rawRing, rawStep,
                static_cast<std::size_t>(kRawResidualMaxPoints));
            if (rawSamples.size() < static_cast<std::size_t>(kRawResidualMinPoints)) {
                rawRejectReason = "too_few_points=" +
                    std::to_string(rawSamples.size());
            } else {
                std::vector<double> dists;
                dists.reserve(rawSamples.size());
                for (const auto& pt : rawSamples) {
                    dists.push_back(distToRingPt(pt));
                }
                std::sort(dists.begin(), dists.end());
                const double median = dists[dists.size() / 2];
                const double q90 = dists[static_cast<std::size_t>(
                    0.9 * static_cast<double>(dists.size() - 1))];
                const double medianGate = std::max(1.5 * pixelSize, 0.6);
                const double q90Gate = std::max(3.0 * pixelSize, 1.5);
                if (median > medianGate || q90 > q90Gate) {
                    rawRejectReason = "distance median=" +
                        std::to_string(median) + " q90=" + std::to_string(q90);
                } else {
                    smoothAssociations.reserve(smoothSamples.size());
                    rawAssociations.reserve(rawSamples.size());
                    for (const auto& pt : smoothSamples) {
                        smoothAssociations.push_back(associateEdge(pt));
                    }
                    for (const auto& pt : rawSamples) {
                        rawAssociations.push_back(associateEdge(pt));
                    }
                    const std::size_t associatedRaw = static_cast<std::size_t>(
                        std::count_if(rawAssociations.begin(), rawAssociations.end(),
                                      [](int edge) { return edge >= 0; }));
                    const std::size_t associatedSmooth = static_cast<std::size_t>(
                        std::count_if(smoothAssociations.begin(), smoothAssociations.end(),
                                      [](int edge) { return edge >= 0; }));
                    if (associatedRaw < static_cast<std::size_t>(kRawResidualMinPoints)) {
                        rawRejectReason = "too_few_associated=" +
                            std::to_string(associatedRaw);
                    } else if (associatedSmooth < 6) {
                        rawRejectReason = "too_few_smooth_associated=" +
                            std::to_string(associatedSmooth);
                    } else {
                        rawAccepted = true;
                        const double rawUnit = kRawResidualWeight / associatedRaw;
                        const double smoothUnit =
                            kSmoothResidualWeight / associatedSmooth;
                        const double unitScale = 1.0 / std::max(rawUnit, smoothUnit);
                        const double rawPointWeight = rawUnit * unitScale;
                        const double smoothPointWeight = smoothUnit * unitScale;
                        for (std::size_t i = 0; i < smoothSamples.size(); ++i) {
                            if (smoothAssociations[i] < 0) continue;
                            residualCloud->push_back(smoothSamples[i]);
                            residualWeights.push_back(smoothPointWeight);
                        }
                        for (std::size_t i = 0; i < rawSamples.size(); ++i) {
                            if (rawAssociations[i] < 0) continue;
                            residualCloud->push_back(rawSamples[i]);
                            residualWeights.push_back(rawPointWeight);
                        }
                        std::cerr << "[RawResidual] fid=" << source_feature_id_
                                  << " part=" << partIndex
                                  << " raw_samples=" << rawSamples.size()
                                  << " raw_associated=" << associatedRaw
                                  << " smooth_associated=" << associatedSmooth
                                  << " source=fid"
                                  << " median_to_smooth=" << median
                                  << " q90_to_smooth=" << q90
                                  << " accepted=1" << std::endl;
                        std::cerr << "[DualResidual] fid=" << source_feature_id_
                                  << " part=" << partIndex
                                  << " raw_target=" << kRawResidualWeight
                                  << " smooth_target=" << kSmoothResidualWeight
                                  << " raw_total_weight="
                                  << rawPointWeight * associatedRaw
                                  << " smooth_total_weight="
                                  << smoothPointWeight * associatedSmooth
                                  << " raw_share=" << kRawResidualWeight
                                  << " rejected_raw="
                                  << (rawSamples.size() - associatedRaw)
                                  << std::endl;
                    }
                }
            }
        }
        if (!rawAccepted) {
            std::cerr << "[RawResidual] fid=" << source_feature_id_
                      << " part=" << partIndex
                      << " raw_samples=" << rawSamples.size()
                      << " smooth_samples=" << smoothSamples.size()
                      << " source=fid accepted=0 reason="
                      << rawRejectReason << std::endl;
        }
    }

    if (!rawAccepted) {
        for (const auto& pt : smoothSamples) residualCloud->push_back(pt);
        residualWeights.clear();
    }
    setResidualWeights(residualWeights);
    if (residualCloud->size() < 6) {
        usedFallback = true;
        std::cerr << "[TopologyFallbackToVDP] fid=" << source_feature_id_
                  << " reason=insufficient_residual" << std::endl;
        return {};
    }
    std::cerr << "[DenseResidual] fid=" << source_feature_id_
              << " residual_points=" << residualCloud->size()
              << " raw=" << (rawAccepted ? 1 : 0) << std::endl;

    if (rawAccepted && kDumpRawResidualPoints) {
        auto& debugPoints = RawResidualDebugPoints();
        auto appendDebug = [&](const pcl::PointXYZ& pt, int src,
                               double weight, int edge) {
            if (debugPoints.size() >= kRawResidualDebugMaxPoints) return;
            RawResidualDebugPoint dbg;
            dbg.fid = source_feature_id_;
            dbg.part = partIndex;
            dbg.pt = pt;
            dbg.src = src;
            dbg.weight = edge >= 0 ? weight : 0.0;
            dbg.assocEdge = edge;
            debugPoints.push_back(dbg);
        };
        const std::size_t smoothCount = static_cast<std::size_t>(
            std::count_if(smoothAssociations.begin(), smoothAssociations.end(),
                          [](int edge) { return edge >= 0; }));
        const std::size_t rawCount = static_cast<std::size_t>(
            std::count_if(rawAssociations.begin(), rawAssociations.end(),
                          [](int edge) { return edge >= 0; }));
        const double rawUnit = kRawResidualWeight / std::max<std::size_t>(rawCount, 1);
        const double smoothUnit =
            kSmoothResidualWeight / std::max<std::size_t>(smoothCount, 1);
        const double unitScale = 1.0 / std::max(rawUnit, smoothUnit);
        for (std::size_t i = 0; i < smoothSamples.size(); ++i) {
            appendDebug(smoothSamples[i], 0, smoothUnit * unitScale,
                        smoothAssociations[i]);
        }
        for (std::size_t i = 0; i < rawSamples.size(); ++i) {
            appendDebug(rawSamples[i], 1, rawUnit * unitScale,
                        rawAssociations[i]);
        }
    }

    // ---- 6. Ceres 平差(方向系统角度作为 preferred) ----
    std::vector<double> preferredAngles;
    if (!isMultiDirection) {
        // 单方向: 只给主系统角, 配合求解器单轴强制
        preferredAngles.push_back(dirBuild.systems.front().angleRad);
    } else {
        for (const auto& s : dirBuild.systems) {
            preferredAngles.push_back(s.angleRad);
        }
    }

    // 优化器将结果写入成员 best_hypothesis，输入多边形保持不变；
    // preserve_topology=true：顶点数固定为链数，禁止删点式预处理
    optimizeWithHardConstraints(
        topologyPolygon, residualCloud, isMultiDirection, preferredAngles,
        &fixedEdgeAssignments, true);
    std::vector<pcl::PointXYZ> ceresResult = best_hypothesis;

    // ---- 7. Ceres 结果验收 ----
    // 三级回退: Ceres 结果不合格 → 退回已验证合格的候选拓扑
    // (直接输出的候选同样要过兜底局部规则性检查);
    // 候选不合格已在第 4 步直接转 VDP; VDP 失败由调用方回退最优假设。
    removeDuplicatePoints2D(ceresResult, 0.05f);
    // 候选直接输出前的局部规则性闸门(短斜边/尖刺/锯齿)
    auto candidateUsable = [&]() {
        outlineRegular::DirectionContextOut ctx;
        ctx.valid = true;
        ctx.completeEvidence = true;
        ctx.multiDirection = isMultiDirection;
        ctx.primaryAngle = dirBuild.systems.front().angleRad;
        ctx.systemAngles = systemAngles;
        return CheckFallbackLocalRegularity(
            topologyPolygon, ctx, source_feature_id_, partIndex,
            "topology_candidate");
    };
    if (ceresResult.size() < 3) {
        std::cerr << "[TopologyCeresReject] fid=" << source_feature_id_
                  << " reason=empty_result" << std::endl;
        const std::string localReason = candidateUsable();
        if (!localReason.empty()) {
            usedFallback = true;
            std::cerr << "[TopologyFallbackToVDP] fid=" << source_feature_id_
                      << " part=" << partIndex
                      << " reason=candidate_local_" << localReason << std::endl;
            return {};
        }
        std::cerr << "[TopologyAccept] fid=" << source_feature_id_
                  << " source=candidate vertices=" << topologyPolygon.size() << std::endl;
        return topologyPolygon;
    }
    const std::string ceresReason = qualityCheck(ceresResult, 2.5);
    if (!ceresReason.empty()) {
        std::cerr << "[TopologyCeresReject] fid=" << source_feature_id_
                  << " reason=" << ceresReason << std::endl;
        const std::string localReason = candidateUsable();
        if (!localReason.empty()) {
            usedFallback = true;
            std::cerr << "[TopologyFallbackToVDP] fid=" << source_feature_id_
                      << " part=" << partIndex
                      << " reason=candidate_local_" << localReason << std::endl;
            return {};
        }
        std::cerr << "[TopologyAccept] fid=" << source_feature_id_
                  << " source=candidate vertices=" << topologyPolygon.size() << std::endl;
        return topologyPolygon;
    }
    // 双轮廓质量对比: 结果边界采样到 raw/smooth 边界的 mean/q90
    if (rawAccepted && rawRing) {
        auto distToSegments = [&](const pcl::PointXYZ& pt,
                                  const std::vector<const std::vector<pcl::PointXYZ>*>& rings) {
            double best = std::numeric_limits<double>::max();
            for (const auto* ringPtr : rings) {
                const auto& seg = *ringPtr;
                for (std::size_t i = 0; i < seg.size(); ++i) {
                    const auto& a = seg[i];
                    const auto& b = seg[(i + 1) % seg.size()];
                    const double dx = b.x - a.x;
                    const double dy = b.y - a.y;
                    const double lenSq = dx * dx + dy * dy;
                    if (lenSq < 1e-12) continue;
                    double t = ((pt.x - a.x) * dx + (pt.y - a.y) * dy) / lenSq;
                    t = std::clamp(t, 0.0, 1.0);
                    const double d = std::hypot(
                        pt.x - a.x - t * dx, pt.y - a.y - t * dy);
                    best = std::min(best, d);
                }
            }
            return best;
        };
        std::vector<const std::vector<pcl::PointXYZ>*> refRings = {rawRing};
        std::vector<const std::vector<pcl::PointXYZ>*> smoothRef = {&initialRing};
        std::vector<double> rawDists, smoothDists;
        for (std::size_t i = 0; i < ceresResult.size(); ++i) {
            const auto& p1 = ceresResult[i];
            const auto& p2 = ceresResult[(i + 1) % ceresResult.size()];
            const double len = std::hypot(p2.x - p1.x, p2.y - p1.y);
            const int steps = std::max(1, static_cast<int>(len / 0.5));
            for (int k = 0; k < steps; ++k) {
                pcl::PointXYZ pt;
                const double t = static_cast<double>(k) / steps;
                pt.x = static_cast<float>(p1.x + t * (p2.x - p1.x));
                pt.y = static_cast<float>(p1.y + t * (p2.y - p1.y));
                rawDists.push_back(distToSegments(pt, refRings));
                smoothDists.push_back(distToSegments(pt, smoothRef));
            }
        }
        auto statsOf = [](std::vector<double>& v) {
            std::sort(v.begin(), v.end());
            return std::make_pair(
                v.empty() ? 0.0 : v[v.size() / 2],
                v.empty() ? 0.0 : v[static_cast<std::size_t>(v.size() * 0.9)]);
        };
        const auto rs = statsOf(rawDists);
        const auto ss = statsOf(smoothDists);
        const double resultArea = std::abs(polygonArea2D(ceresResult));
        std::cerr << "[DualQuality] fid=" << source_feature_id_
                  << " part=" << partIndex
                  << " raw_median=" << rs.first << " raw_q90=" << rs.second
                  << " smooth_median=" << ss.first << " smooth_q90=" << ss.second
                  << " area_ratio_raw=" << (resultArea / std::max(rawAreaCached, 1e-6))
                  << " area_ratio_smooth="
                  << (resultArea / std::max(std::abs(polygonArea2D(initialRing)), 1e-6))
                  << std::endl;
    }
    std::cerr << "[TopologyAccept] fid=" << source_feature_id_
              << " source=ceres vertices=" << ceresResult.size()
              << " chains=" << topoChains.size()
              << " mode=" << (isMultiDirection ? "multi" : "single") << std::endl;
    return ceresResult;
}
// ===== Mask-only 圆弧检测/恢复的公共实现(全局作用域) =====
std::vector<MaskConicArc> DetectMaskConicArcs(
    const std::vector<pcl::PointXYZ>& smoothRing,
    const std::vector<pcl::PointXYZ>& rawRing,
    double pixelSize, long long fid, int partIdx)
{
    std::vector<MaskConicArc> candidates;
    const std::size_t n = smoothRing.size();
    if (n < 12) return candidates;
    const SmoothRawMap rawMap = BuildSmoothRawMap(smoothRing, rawRing);
    std::string bestReject = "none";
    double bestRejectScore = -1;

    struct Cand { std::size_t s, e; bool w; MaskConicArc a; };
    std::vector<Cand> all;
    const std::size_t minC = 5, maxC = n * 3 / 4;

    std::vector<double> cumArc(n, 0.0);
    for (std::size_t i = 1; i < n; ++i)
        cumArc[i] = cumArc[i-1] + std::hypot(
            smoothRing[i].x - smoothRing[i-1].x,
            smoothRing[i].y - smoothRing[i-1].y);
    const double totalArc = cumArc[n-1] + std::hypot(
        smoothRing[0].x - smoothRing[n-1].x,
        smoothRing[0].y - smoothRing[n-1].y);

    for (std::size_t start = 0; start < n; ++start) {
        for (std::size_t cnt = minC; cnt <= maxC; ++cnt) {
            if (cnt >= n) break;
            const std::size_t end = (start + cnt - 1) % n;
            const bool w = (start + cnt - 1) >= n;
            double iArc = w ? (totalArc - cumArc[start] + cumArc[end])
                            : (cumArc[end] - cumArc[start]);
            if (iArc < 10.0) continue;

            std::vector<pcl::PointXYZ> fitPts;
            int sourceN = 0;
            if (rawMap.valid) {
                RawInterval ri = ExtractRawInterval(rawMap, start, end);
                if (ri.valid && ri.pointCount >= 3) {
                    fitPts = CollectRawPoints(rawRing, ri);
                    sourceN = ri.pointCount;
                }
            }
            if (fitPts.empty()) {
                const std::size_t c2 = w ? (n - start + end + 1) : (end - start + 1);
                fitPts.reserve(c2);
                for (std::size_t k = 0; k < c2; ++k)
                    fitPts.push_back(smoothRing[(start + k) % n]);
                sourceN = (int)c2;
            }
            bool ok = sourceN >= 10;
            if (!ok && sourceN >= 5 && iArc >= 10.0) {
                fitPts = WeightedSample(fitPts, 12);
                ok = true;
            }
            if (!ok) continue;

            std::vector<pcl::PointXYZ> smoothIvl;
            {
                const std::size_t c3 = w ? (n - start + end + 1) : (end - start + 1);
                smoothIvl.reserve(c3);
                for (std::size_t k = 0; k < c3; ++k)
                    smoothIvl.push_back(smoothRing[(start + k) % n]);
            }
            MaskConicArc a;
            const std::string reason = EvalCircleV2(fitPts, smoothIvl, sourceN, pixelSize, a);
            if (reason.empty()) {
                a.startIdx = start; a.endIdx = end; a.wrapsZero = w;
                a.startPoint = smoothRing[start];
                a.endPoint = smoothRing[end];
                all.push_back({start, end, w, a});
            } else {
                if (a.arcLength > bestRejectScore) {
                    bestRejectScore = a.arcLength;
                    bestReject = reason;
                }
            }
        }
    }
    {
        std::sort(all.begin(), all.end(),
            [](const Cand& a, const Cand& b) { return a.a.score > b.a.score; });
        std::vector<bool> used(n, false);
        for (const auto& c : all) {
            bool ov = false;
            const std::size_t cnt = c.w ? (n - c.s + c.e + 1) : (c.e - c.s + 1);
            for (std::size_t k = 0; k < cnt; ++k)
                if (used[(c.s + k) % n]) { ov = true; break; }
            if (ov) continue;
            for (std::size_t k = 0; k < cnt; ++k) used[(c.s + k) % n] = true;
            candidates.push_back(c.a);
        }
    }
    std::sort(candidates.begin(), candidates.end(),
        [](const MaskConicArc& a, const MaskConicArc& b) {
            return a.startIdx < b.startIdx; });
    std::cerr << "[MaskCurveDetect] fid=" << fid << " part=" << partIdx
              << " raw_candidates=" << all.size()
              << " accepted=" << candidates.size()
              << " best_reject=" << bestReject
              << " pixel_size=" << pixelSize << std::endl;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        const auto& a = candidates[i];
        std::cerr << "[MaskCurve] fid=" << fid << " part=" << partIdx
                  << " index=" << i << " type=circle"
                  << " start=" << a.startIdx << " end=" << a.endIdx
                  << " wrap=" << (a.wrapsZero ? 1 : 0)
                  << " source_n=" << a.supportCount
                  << " length=" << a.arcLength
                  << " chord=" << a.chordLength
                  << " sweep_deg=" << a.sweepDeg
                  << " rmse=" << a.rmse
                  << " line_rmse=" << a.lineRmse
                  << " score=" << a.score << std::endl;
    }
    return candidates;
}
std::vector<pcl::PointXYZ> RestoreMaskConicArcs(
    const std::vector<pcl::PointXYZ>& candidate,
    const std::vector<MaskConicArc>& arcs,
    double pixelSize, long long fid, int partIdx)
{
    if (arcs.empty() || candidate.size() < 4) return candidate;
    std::vector<pcl::PointXYZ> result = candidate;
    for (std::size_t ai = 0; ai < arcs.size(); ++ai) {
        const auto& arc = arcs[ai];
        auto proj = [&](const pcl::PointXYZ& pt) {
            double bd = 1e18; std::size_t be = 0; double bt = 0.0;
            for (std::size_t e = 0; e < result.size(); ++e) {
                const auto& a = result[e];
                const auto& b = result[(e + 1) % result.size()];
                const double dx = b.x - a.x, dy = b.y - a.y;
                const double ls = dx * dx + dy * dy;
                if (ls < 1e-12) continue;
                double t = ((pt.x - a.x) * dx + (pt.y - a.y) * dy) / ls;
                t = std::clamp(t, 0.0, 1.0);
                const double d = std::hypot(pt.x - a.x - t * dx, pt.y - a.y - t * dy);
                if (d < bd) { bd = d; be = e; bt = t; }
            }
            return std::make_tuple(be, bt, bd);
        };
        auto [sE, sT, sD] = proj(arc.startPoint);
        auto [eE, eT, eD] = proj(arc.endPoint);
        const double lim = std::max(3.0 * pixelSize, 1.5);
        if (sD > lim || eD > lim) {
            std::cerr << "[MaskCurveRestore] fid=" << fid << " part=" << partIdx
                      << " index=" << ai << " restored=0 reason=endpoint_dist"
                      << " s=" << sD << " e=" << eD << std::endl;
            continue;
        }
        const std::size_t rn = result.size();
        double fwdL = 0.0, bwdL = 0.0;
        std::size_t fwdE = 0, bwdE = 0;
        { std::size_t e = sE;
          while (e != eE) {
              const auto& a = result[e]; const auto& b = result[(e + 1) % rn];
              fwdL += std::hypot(b.x - a.x, b.y - a.y);
              e = (e + 1) % rn; ++fwdE; if (fwdE > rn) break; } }
        { std::size_t e = sE;
          while (e != eE) {
              const auto& a = result[(e + rn - 1) % rn]; const auto& b = result[e];
              bwdL += std::hypot(b.x - a.x, b.y - a.y);
              e = (e + rn - 1) % rn; ++bwdE; if (bwdE > rn) break; } }
        const bool fwd = std::abs(fwdL - arc.arcLength) <= std::abs(bwdL - arc.arcLength);
        pcl::PointXYZ pS; { const auto& a = result[sE]; const auto& b = result[(sE + 1) % rn];
          pS.x = (float)(a.x + sT * (b.x - a.x)); pS.y = (float)(a.y + sT * (b.y - a.y));
          pS.z = arc.startPoint.z; }
        pcl::PointXYZ pE; { const auto& a = result[eE]; const auto& b = result[(eE + 1) % rn];
          pE.x = (float)(a.x + eT * (b.x - a.x)); pE.y = (float)(a.y + eT * (b.y - a.y));
          pE.z = arc.endPoint.z; }
        const double srcLen = arc.chordLength;
        const double dstLen = std::hypot(pE.x - pS.x, pE.y - pS.y);
        if (srcLen < 1e-6 || dstLen < 1e-6) continue;
        const double scale = dstLen / srcLen;
        if (scale < 0.85 || scale > 1.20) continue;
        const double sa = std::atan2(arc.endPoint.y - arc.startPoint.y,
                                      arc.endPoint.x - arc.startPoint.x);
        const double da = std::atan2(pE.y - pS.y, pE.x - pS.x);
        const double rot = da - sa;
        const double cR = std::cos(rot), sR = std::sin(rot);
        auto xf = [&](const pcl::PointXYZ& p) {
            const double dx = p.x - arc.startPoint.x;
            const double dy = p.y - arc.startPoint.y;
            pcl::PointXYZ q;
            q.x = (float)(pS.x + scale * (dx * cR - dy * sR));
            q.y = (float)(pS.y + scale * (dx * sR + dy * cR));
            q.z = p.z; return q; };
        const double step = std::clamp(pixelSize, 0.3, 0.8);
        const double aLen = std::abs(arc.sweepAngle) * arc.radius;
        const int steps = std::max(8, (int)std::ceil(aLen / step));
        std::vector<pcl::PointXYZ> cs; cs.reserve(steps);
        for (int k = 0; k < steps; ++k) {
            const double t = (double)k / (steps - 1);
            const double ang = arc.startAngle + t * arc.sweepAngle;
            pcl::PointXYZ p;
            p.x = (float)(arc.cx + arc.radius * std::cos(ang));
            p.y = (float)(arc.cy + arc.radius * std::sin(ang));
            p.z = arc.startPoint.z; cs.push_back(xf(p)); }
        std::vector<bool> inSpan(rn, false);
        if (fwd) { std::size_t e = (sE + 1) % rn;
            for (std::size_t g = 0; g < rn && e != (eE + 1) % rn; e = (e + 1) % rn, ++g)
                inSpan[e] = true; }
        else { std::size_t e = eE;
            for (std::size_t g = 0; g < rn && e != sE; e = (e + 1) % rn, ++g)
                inSpan[e] = true; }
        std::vector<pcl::PointXYZ> nr;
        if (fwd) { for (std::size_t k = 0; k < rn; ++k) {
                const std::size_t idx = (eE + 1 + k) % rn;
                if (!inSpan[idx]) nr.push_back(result[idx]); } }
        else { for (std::size_t k = 0; k < rn; ++k) {
                const std::size_t idx = (sE + 1 + k) % rn;
                if (!inSpan[idx]) nr.push_back(result[idx]); } }
        for (const auto& p : cs) nr.push_back(p);
        if (nr.size() > 1) { const auto& f = nr.front(); const auto& l = nr.back();
            if (std::hypot(f.x - l.x, f.y - l.y) < 0.01) nr.pop_back(); }
        if (nr.size() < 4) continue;
        if (!isSimplePolygon2D(nr)) continue;
        const double oldA = std::abs(polygonArea2D(result));
        const double newA = std::abs(polygonArea2D(nr));
        const double ar = newA / std::max(oldA, 1e-6);
        if (ar < 0.90 || ar > 1.10) continue;
        result = nr;
        std::cerr << "[MaskCurveRestore] fid=" << fid << " part=" << partIdx
                  << " index=" << ai << " restored=1 scale=" << scale
                  << " samples=" << cs.size() << std::endl;
    }
    return result;
}
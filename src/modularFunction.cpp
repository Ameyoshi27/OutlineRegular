// =============================================================================
// modularFunction.cpp
// 作用：实现 modularFunction 类的 Douglas-Peucker 抽稀 与 点到线段距离 计算。
// =============================================================================

#include "modularFunction.h"

#include <algorithm>

namespace XG {

// ===== douglasPeucker =====
// 作用：Douglas-Peucker 递归抽稀算法——在保留轮廓形状的前提下，按阈值 epsilon
//       去掉偏离直线很小的中间点，从而减少多边形顶点数。
// 参数：points  - 待简化的有序点列；
//       epsilon - 距离阈值，点到首尾连线的距离小于它则被丢弃；
//       out     - 输出简化后的点列。
void modularFunction::douglasPeucker(std::vector<Points_dp>& points, double epsilon, std::vector<Points_dp>& out)
{
    out.clear();
    if (points.size() < 2) {
        throw std::invalid_argument("Not enough points to simplify");
    }

    // 找出离首尾连线最远的点(即对形状贡献最大的点)
    double max_distance = 0.0;
    std::size_t index = 0;
    const Points_dp start_point = points.front();
    const Points_dp end_point = points.back();

    for (std::size_t i = 1; i + 1 < points.size(); ++i) {
        const double dist = perpendicularDistance(points[i], start_point, end_point);
        if (dist > max_distance) {
            index = i;
            max_distance = dist;
        }
    }

    if (max_distance > epsilon) {
        // 最大距离仍超阈值：以该点为切分点，对左右两段分别递归
        std::vector<Points_dp> rec_results1;
        std::vector<Points_dp> rec_results2;
        std::vector<Points_dp> first_line(points.begin(), points.begin() + index + 1);
        std::vector<Points_dp> last_line(points.begin() + index, points.end());

        douglasPeucker(first_line, epsilon, rec_results1);
        douglasPeucker(last_line, epsilon, rec_results2);

        // 拼接两段结果(去掉重复的切分点)
        out = rec_results1;
        if (!rec_results2.empty()) {
            out.insert(out.end(), rec_results2.begin() + 1, rec_results2.end());
        }
    } else {
        // 最大距离足够小：用首尾两点代替整段
        out.push_back(points.front());
        out.push_back(points.back());
    }
}

// ===== perpendicularDistance =====
// 作用：计算点 point 到线段 lineStart-lineEnd 的最短距离。
//       (把点到直线的距离 clamp 到线段范围内，退化为端点距离时直接返回。)
// 参数：point     - 待计算的点；
//       lineStart - 线段起点；
//       lineEnd   - 线段终点。
// 返回：最短距离。
double modularFunction::perpendicularDistance(const Points_dp& point, const Points_dp& lineStart, const Points_dp& lineEnd) const
{
    const double a = point.x - lineStart.x;
    const double b = point.y - lineStart.y;
    const double c = lineEnd.x - lineStart.x;
    const double d = lineEnd.y - lineStart.y;
    const double len_sq = c * c + d * d;
    if (len_sq == 0.0) {
        // 线段退化为点：直接返回两点距离
        return std::sqrt(a * a + b * b);
    }

    // 投影参数 param 限制在 [0,1]，即垂足落在线段上
    double param = (a * c + b * d) / len_sq;
    param = std::max(0.0, std::min(1.0, param));
    const double xx = lineStart.x + param * c;
    const double yy = lineStart.y + param * d;

    const double dx = point.x - xx;
    const double dy = point.y - yy;
    return std::sqrt(dx * dx + dy * dy);
}

} // namespace XG

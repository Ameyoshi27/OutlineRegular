// =============================================================================
// modularFunction.h
// 作用：通用的几何小工具(命名空间 XG)。
//       提供 Douglas-Peucker 多边形抽稀(简化)算法 和 点云模型分辨率估计，
//       被 outlineRegular(轮廓规则化) 调用，用于简化多边形顶点、估计点云密度。
// =============================================================================

#pragma once

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace XG {

// 几何工具类
class modularFunction {
public:
    // 二维点(双精度)，供 Douglas-Peucker 算法使用
    struct Points_dp {
        double x = 0.0;
        double y = 0.0;

        Points_dp() = default;
        Points_dp(double x_, double y_) : x(x_), y(y_) {}

        // 相等判断(带微小容差)，用于去重
        bool operator==(const Points_dp& other) const {
            const double eps = 1e-9;
            return std::fabs(x - other.x) < eps && std::fabs(y - other.y) < eps;
        }
    };

    // Douglas-Peucker 抽稀：按阈值 epsilon 简化 points，结果写入 out
    void douglasPeucker(std::vector<Points_dp>& points, double epsilon, std::vector<Points_dp>& out);

    // 计算点云的模型分辨率(平均近邻距离)：只统计距离 <0.5 的近邻并取平均；
    // 用于衡量点云密度，进而决定规则化时的参数(如自适应 lambda、最小边长)。
    template <typename PointT>
    double computeModelResolution(pcl::PointCloud<PointT>& cloud)
    {
        if (cloud.size() < 2) {
            return 0.0;
        }

        double sum = 0.0;
        int count = 0;
        const double max_neighbor_distance = 0.5; // 只考虑 0.5 以内的近邻，避免远点干扰

        // 暴力 O(n^2) 搜索每个点的最近邻距离
        for (std::size_t i = 0; i < cloud.size(); ++i) {
            double best = std::numeric_limits<double>::max();
            for (std::size_t j = 0; j < cloud.size(); ++j) {
                if (i == j) continue;
                const double dx = static_cast<double>(cloud[i].x) - cloud[j].x;
                const double dy = static_cast<double>(cloud[i].y) - cloud[j].y;
                const double dz = static_cast<double>(cloud[i].z) - cloud[j].z;
                const double dist2 = dx * dx + dy * dy + dz * dz;
                if (dist2 > 0.0 && dist2 < best) {
                    best = dist2;
                }
            }
            if (std::isfinite(best)) {
                const double dist = std::sqrt(best);
                if (dist > 0.0 && dist < max_neighbor_distance) {
                    sum += dist;
                    ++count;
                }
            }
        }

        return count > 0 ? sum / static_cast<double>(count) : 0.0;
    }

private:
    // 计算点 point 到线段 lineStart-lineEnd 的最短距离(供 Douglas-Peucker 使用)
    double perpendicularDistance(const Points_dp& point, const Points_dp& lineStart, const Points_dp& lineEnd) const;
};

} // namespace XG

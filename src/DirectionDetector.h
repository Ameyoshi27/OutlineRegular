// =============================================================
// DirectionDetector.h
// 统一主方向检测模块: 在平滑初始轮廓上一次性检测方向系统,
// 结果作为"真值"传给后续所有简化/平差路径。
//
// 设计原则:
// - 每个建筑环只检测一次, 输出包含所有方向系统
// - 后续 Ceres 中 θ_base 固定为检测结果, 不再优化
// - 不设自由边: 每条边必须归入平行或垂直于某个系统
// =============================================================
#pragma once

#include <pcl/point_types.h>
#include <vector>
#include <string>
#include <cmath>

// 检测结果: 一个方向系统
// 候选(active=false)只表示"方向迹象"; active=true 表示具有足够
// 物理墙长或统计权重, 允许约束规则化。主峰恒 active。
struct DetectedDirectionSystem {
    double angleRad = 0.0;          // 折叠角 [0, π/2)
    int chainCount = 0;             // 归入的稳定链数
    double totalLength = 0.0;       // 稳定链总长(m)
    double weight = 0.0;            // 加权支持(长度^指数)
    double meanRmse = 0.0;          // 平均拟合残差
    double concentration = 0.0;     // 圆集中度(四倍角合向量长度)
    double extent = 0.0;            // 空间覆盖(bbox对角线)
    double confidence = 0.0;        // 综合置信度 [0,1]
    double prominence = 0.0;        // KDE 峰显著度(相对主峰高度 [0,1])
    double height = 0.0;            // KDE 峰高(相对主峰归一, 主峰=1.0)
    double avgEdgeLen = 0.0;        // 系统内平均边长(m), 诊断锯齿伪系统
    bool active = false;            // 生效方向(下游只消费 active 系统)
    // 直墙证据救回: prominence 不达标但 ±5° 内有独立直墙支撑的峰。
    // 救回峰的角度用直链 4θ 加权圆均值精化, active 走直墙独立通道
    bool evidenceBacked = false;
    int straightSupportCount = 0;
    double straightSupportLength = 0.0;
    // 输入证据(来自原始边链统计, 不从规则化结果反推; 供 A/B 评估)
    double inputLongestChain = 0.0;      // 最长支持直链(m)
    int inputSupportRuns = 0;            // 支持链沿环序的连续 run 数
    double inputLongestRunRatio = 0.0;   // 最长 run / 支持总长
    double inputPerimeterRatio = 0.0;    // 支持总长 / 周长
};

// 检测结果: 整体
struct DetectedDirectionResult {
    bool valid = false;              // 是否有可信方向
    bool multiDirection = false;    // 是否多方向(active系统数>=2)
    double primaryAngle = 0.0;      // 主方向折叠角
    std::vector<DetectedDirectionSystem> systems;
    // PCA 纠偏诊断: rawPrimaryAngle 是 KDE 原始峰角,
    // primaryAngle 可能已被单方向 PCA 纠偏覆盖
    double rawPrimaryAngle = 0.0;
    bool primaryRefined = false;
    std::string refinementReason;
    // 全局 PCA 主轴诊断(供下游 alternate_direction 重试)
    bool pcaValid = false;
    double pcaAngleRad = 0.0;
    double pcaAxisRatio = 0.0;
    double pcaAnisotropy = 0.0;
    // 诊断
    int totalChains = 0;
    int stableChains = 0;
    double totalStableLength = 0.0;
    double concentration = 0.0;     // 主系统集中度
    std::string rejectReason;       // 无效时的原因

    // 便捷: 所有系统角度
    std::vector<double> systemAngles() const {
        std::vector<double> angles;
        for (const auto& s : systems) angles.push_back(s.angleRad);
        return angles;
    }
};

// =============================================================
// 主检测函数
// 输入:
//   smoothRing - 平滑后的初始轮廓环(局部坐标)
//   rawRing    - 原始像素轮廓(可选, 为空则只用平滑环)
//   pixelSize  - TIF 真实像素尺寸(m)
//   fid, partIdx - 日志标识
// 输出:
//   DetectedDirectionResult
// =============================================================
DetectedDirectionResult DetectBuildingDirection(
    const std::vector<pcl::PointXYZ>& smoothRing,
    const std::vector<pcl::PointXYZ>& rawRing,
    double pixelSize,
    long long fid = -1,
    int partIdx = 0);

// =============================================================
// 辅助: 边到最近方向系统的分配
// 返回: 系统索引(平行或垂直都算), -1 表示无法分配
// =============================================================
int AssignEdgeToDirectionSystem(
    double edgeAngle,
    const std::vector<double>& systemAngles,
    double maxAssignDeg = 20.0);

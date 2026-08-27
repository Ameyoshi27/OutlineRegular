// =============================================================
// DirectionDetector.cpp
// Unified building direction detector for topology, VDP and fallback paths.
//
// 核心估计器: 边链化 + 长度加权圆环 KDE + 峰检测(众数)。
// 设计动机(替代早期贪心种子聚类 + 圆均值):
// - 链化恢复空间连续性: 平滑环上斜墙是方向渐变的短边串,
//   逐边投票会在角度空间被摊薄, 而精确对齐栅格轴的楼梯边
//   高度集中, 密度峰被伪方向占据。以链为单位投票(链角=端点
//   向量)恢复墙体的真实长度话语权;
// - 峰位置=众数, 不被归组进来的异质边拉偏(均值污染);
// - 软密度贡献替代硬归组, 弥散证据不再触发整体拒绝;
// - 短链以 len^1.5 降权参与而非被排除, 第二方向的结构
//   证据不再丢失;
// - 锯齿伪方向在密度上是低而宽的隆起, 被 prominence 准入
//   自然过滤, 不依赖链数/份额组合阈值。
// =============================================================
#include "DirectionDetector.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <utility>

namespace {

// ---- 核心参数(物理意义明确, 不随数据集漂移) ----
constexpr double kKdeBandwidthDeg = 5.0;       // 墙体方向自然散布量级
constexpr double kLenExponent = 1.5;           // 长边主导, 短边保留贡献
constexpr double kPeakProminenceRatio = 0.12;  // 次峰至少高出鞍部 12% 主峰高度
constexpr double kMinPeakSeparationDeg = 10.0; // 更近的方向在像素尺度不可分
constexpr int kMaxFinalSystems = 5;
// ---- 候选峰 → 生效方向的判据 ----
// weightShare: len^指数权重对长边超线性放大, 短而真实的墙体在权重上吃亏
// lengthShare: 候选系统之间的实际墙长占比, 物理支撑的直接度量
// minPhysicalLength: 低于该墙长的方向不约束规则化(绝对量 + 周长比例)
// 不使用 chainCount(累计的是链内原始边数, 语义不是独立链数)
constexpr double kActiveSystemMinWeightShare = 0.12;
constexpr double kActiveSystemMinLengthShare = 0.18;
constexpr double kActiveSystemMinLengthMeters = 6.0;
constexpr double kActiveSystemMinPerimeterShare = 0.08;
constexpr double kChainMergeTurnDeg = 25.0;    // 链内相邻边最大转向

constexpr int kKdeBins = 360;                  // 折叠空间 [0°,90°), 0.25°/格

double foldedAngle90(double angle) {
    angle = std::fmod(angle, M_PI / 2.0);
    if (angle < 0.0) angle += M_PI / 2.0;
    return angle;
}

double foldedDistance90(double a, double b) {
    const double period = M_PI / 2.0;
    double difference = std::fmod(std::abs(a - b), period);
    return std::min(difference, period - difference);
}

double wrappedTurn(double a, double b) {
    double turn = std::abs(a - b);
    if (turn > M_PI) turn = 2.0 * M_PI - turn;
    return turn;
}

struct DetChain {
    double angleRad = 0.0;  // 折叠角 [0, π/2)
    double length = 0.0;
    double weight = 0.0;    // length^kLenExponent
    int edgeCount = 0;
};

// ---- 边链化: 空间连续且转向平缓的边合并成链, 链角=端点向量 ----
std::vector<DetChain> buildChainsFromRing(
    const std::vector<pcl::PointXYZ>& ring, double pixelSize) {
    struct RingEdge {
        double dirRad = 0.0;   // 有向角 [-π, π)
        double length = 0.0;
        std::size_t from = 0;
        std::size_t to = 0;
    };
    std::vector<RingEdge> edgeList;
    edgeList.reserve(ring.size());
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto& first = ring[i];
        const auto& last = ring[(i + 1) % ring.size()];
        const double length = std::hypot(last.x - first.x, last.y - first.y);
        if (length < 1e-6) continue;
        RingEdge edge;
        edge.dirRad = std::atan2(last.y - first.y, last.x - first.x);
        edge.length = length;
        edge.from = i;
        edge.to = (i + 1) % ring.size();
        edgeList.push_back(edge);
    }
    if (edgeList.empty()) return {};

    // 环遍历起点选在转向最大处, 避免一条链被首尾切断
    std::size_t startEdge = 0;
    double maxTurn = -1.0;
    for (std::size_t i = 0; i < edgeList.size(); ++i) {
        const std::size_t j = (i + 1) % edgeList.size();
        const double turn = wrappedTurn(edgeList[j].dirRad, edgeList[i].dirRad);
        if (turn > maxTurn) {
            maxTurn = turn;
            startEdge = j;
        }
    }

    const double mergeTurn = kChainMergeTurnDeg * M_PI / 180.0;
    const double minChainLength = std::max(1.0, 2.0 * pixelSize);
    std::vector<DetChain> chains;
    std::size_t chainFromVertex = 0;
    double chainDir = 0.0;
    double chainLen = 0.0;
    int chainEdges = 0;
    const auto closeChain = [&](std::size_t endVertex) {
        if (chainEdges == 0) return;
        const auto& a = ring[chainFromVertex];
        const auto& b = ring[endVertex];
        const double spanX = b.x - a.x;
        const double spanY = b.y - a.y;
        if (std::hypot(spanX, spanY) > 1e-6 && chainLen >= minChainLength) {
            DetChain chain;
            chain.angleRad = foldedAngle90(std::atan2(spanY, spanX));
            chain.length = chainLen;
            chain.weight = std::pow(chainLen, kLenExponent);
            chain.edgeCount = chainEdges;
            chains.push_back(chain);
        }
        chainEdges = 0;
    };

    for (std::size_t step = 0; step < edgeList.size(); ++step) {
        const std::size_t idx = (startEdge + step) % edgeList.size();
        const auto& edge = edgeList[idx];
        if (chainEdges == 0) {
            chainFromVertex = edge.from;
            chainDir = edge.dirRad;
            chainLen = edge.length;
            chainEdges = 1;
            continue;
        }
        if (wrappedTurn(edge.dirRad, chainDir) < mergeTurn) {
            chainDir = std::atan2(
                std::sin(chainDir) * chainLen + std::sin(edge.dirRad) * edge.length,
                std::cos(chainDir) * chainLen + std::cos(edge.dirRad) * edge.length);
            chainLen += edge.length;
            ++chainEdges;
        } else {
            closeChain(edge.from);  // 链终点 = 上一边终点 = 新边起点
            chainFromVertex = edge.from;
            chainDir = edge.dirRad;
            chainLen = edge.length;
            chainEdges = 1;
        }
    }
    // 收尾: 最后一条链的终点 = 第一条链的起点
    if (chainEdges > 0) {
        closeChain(edgeList[startEdge].from);
    }
    return chains;
}

// 环形域上格点 idx 的 prominence: 两侧各走到第一个更高格点,
// 每侧取路径最低点, prominence = h - max(两侧最低点)。
// 全局最高峰两侧绕一整圈, 退化为 h - 全局最小。
double circularProminence(const std::vector<double>& density,
                          std::size_t idx) {
    const std::size_t n = density.size();
    const double height = density[idx];
    double pathMin[2] = {height, height};
    int side = 0;
    for (int direction : {-1, 1}) {
        for (std::size_t step = 1; step < n; ++step) {
            const std::size_t j =
                (idx + n + direction * static_cast<int>(step)) % n;
            pathMin[side] = std::min(pathMin[side], density[j]);
            if (density[j] > height) break;
        }
        ++side;
    }
    return height - std::max(pathMin[0], pathMin[1]);
}

} // namespace

DetectedDirectionResult DetectBuildingDirection(
    const std::vector<pcl::PointXYZ>& smoothRing,
    const std::vector<pcl::PointXYZ>& rawRing,
    double pixelSize,
    long long fid,
    int partIdx) {
    (void)rawRing;
    DetectedDirectionResult result;
    if (smoothRing.size() < 3) {
        result.rejectReason = "too_few_vertices";
        return result;
    }

    // ---- 链化: 恢复墙体的空间连续性 ----
    const std::vector<DetChain> chains =
        buildChainsFromRing(smoothRing, pixelSize);
    if (chains.empty()) {
        result.rejectReason = "no_valid_edges";
        return result;
    }
    double totalPerimeter = 0.0;
    double totalWeight = 0.0;
    int totalEdges = 0;
    for (const auto& chain : chains) {
        totalPerimeter += chain.length;
        totalWeight += chain.weight;
        totalEdges += chain.edgeCount;
    }
    result.totalChains = totalEdges;
    result.stableChains = static_cast<int>(chains.size());
    result.totalStableLength = totalPerimeter;

    // ---- 变带宽长度加权圆环 KDE(样本=链) ----
    // 每条链的方向不确定度由像素分辨率与链长决定:
    // σ_i = 2·atan(pixelSize/len), 长边是窄而可靠的方向证据,
    // 短边是宽而模糊的证据。统一带宽会把长边的窄峰摊薄,
    // 让精确对齐栅格轴的边(测量偏差系统性地小)虚假占优。
    const double binStep = (M_PI / 2.0) / kKdeBins;
    const double sigmaMin = 2.0 * M_PI / 180.0;
    const double sigmaMax = 12.0 * M_PI / 180.0;
    const int maxCutoffBins =
        static_cast<int>(4.0 * sigmaMax / binStep) + 1;

    std::vector<double> density(kKdeBins, 0.0);
    for (const auto& chain : chains) {
        const double sigma = std::clamp(
            2.0 * std::atan2(pixelSize, chain.length), sigmaMin, sigmaMax);
        const double invTwoSigmaSq = 1.0 / (2.0 * sigma * sigma);
        const int cutoffBins = std::min(
            maxCutoffBins, static_cast<int>(4.0 * sigma / binStep) + 1);
        const int center = static_cast<int>(chain.angleRad / binStep);
        for (int offset = -cutoffBins; offset <= cutoffBins; ++offset) {
            const int bin = ((center + offset) % kKdeBins + kKdeBins) % kKdeBins;
            const double delta = offset * binStep;
            density[bin] += chain.weight * std::exp(-delta * delta * invTwoSigmaSq);
        }
    }
    if (totalWeight > 1e-9) {
        for (auto& value : density) value /= totalWeight;
    }

    // ---- 峰检测: 局部极大 → prominence 准入 → 间隔过滤 ----
    struct Peak {
        std::size_t bin = 0;
        double height = 0.0;
        double prominence = 0.0;
        double angleRad = 0.0;  // 抛物线插值后的峰位置
    };
    std::vector<Peak> peaks;
    for (std::size_t b = 0; b < kKdeBins; ++b) {
        const double prev = density[(b + kKdeBins - 1) % kKdeBins];
        const double next = density[(b + 1) % kKdeBins];
        if (density[b] > prev && density[b] >= next && density[b] > 1e-12) {
            Peak peak;
            peak.bin = b;
            peak.height = density[b];
            peaks.push_back(peak);
        }
    }
    if (peaks.empty()) {
        const auto maxIt = std::max_element(density.begin(), density.end());
        Peak peak;
        peak.bin = static_cast<std::size_t>(maxIt - density.begin());
        peak.height = *maxIt;
        peaks.push_back(peak);
    }
    std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) {
        return a.height > b.height;
    });

    if (fid >= 0) {  // 密度诊断: 前 6 个局部极大(角度/高度)
        std::cerr << "[DirectionKde] fid=" << fid << " part=" << partIdx
                  << " chains=" << chains.size();
        for (std::size_t i = 0; i < peaks.size() && i < 6; ++i) {
            std::cerr << " max" << i << "="
                      << peaks[i].bin * binStep * 180.0 / M_PI
                      << "deg/h=" << peaks[i].height;
        }
        std::cerr << std::endl;
    }

    const double mainHeight = peaks.front().height;
    const double minSeparation = kMinPeakSeparationDeg * M_PI / 180.0;
    std::vector<Peak> accepted;
    for (auto& peak : peaks) {
        if (static_cast<int>(accepted.size()) >= kMaxFinalSystems) break;
        peak.prominence = circularProminence(density, peak.bin);
        if (!accepted.empty() &&
            peak.prominence < kPeakProminenceRatio * mainHeight) {
            continue;
        }
        // 抛物线插值细化峰位置
        const double left = density[(peak.bin + kKdeBins - 1) % kKdeBins];
        const double right = density[(peak.bin + 1) % kKdeBins];
        const double center = density[peak.bin];
        const double denom = left - 2.0 * center + right;
        const double shift = std::abs(denom) > 1e-15
            ? 0.5 * (left - right) / denom : 0.0;
        peak.angleRad = foldedAngle90(
            (static_cast<double>(peak.bin) + std::clamp(shift, -1.0, 1.0))
            * binStep);
        bool tooClose = false;
        for (const auto& kept : accepted) {
            if (foldedDistance90(peak.angleRad, kept.angleRad) < minSeparation) {
                tooClose = true;
                break;
            }
        }
        if (!tooClose) accepted.push_back(peak);
    }

    // ---- 每峰统计: ±1 统计窗口内的链(软归属的统计近似) ----
    const double statWindow = kKdeBandwidthDeg * M_PI / 180.0;
    for (const auto& peak : accepted) {
        DetectedDirectionSystem system;
        system.angleRad = peak.angleRad;
        system.prominence = mainHeight > 1e-12
            ? peak.prominence / mainHeight : 0.0;
        double x4 = 0.0, y4 = 0.0;
        for (const auto& chain : chains) {
            if (foldedDistance90(chain.angleRad, peak.angleRad) > statWindow) continue;
            system.totalLength += chain.length;
            system.weight += chain.weight;
            system.avgEdgeLen += chain.length;
            system.chainCount += chain.edgeCount;
            x4 += chain.weight * std::cos(4.0 * chain.angleRad);
            y4 += chain.weight * std::sin(4.0 * chain.angleRad);
        }
        if (system.chainCount > 0) {
            system.avgEdgeLen /= system.chainCount;
        }
        system.concentration = system.weight > 1e-9
            ? std::hypot(x4, y4) / system.weight : 0.0;
        const double share = totalWeight > 1e-9
            ? system.weight / totalWeight : 0.0;
        system.confidence = std::min(1.0, share * system.concentration);
        result.systems.push_back(system);
    }

    // ---- 候选峰 → 生效方向: 逐峰独立裁决 ----
    // 主峰恒生效; 其余峰需足够物理墙长 且(权重份额或墙长份额达标)。
    // multiDirection 是派生状态(active 数量), 不再由固定 systems[1]
    // 的权重份额决定, 也不因 multi=true 放行全部弱候选。
    double candidateLengthSum = 0.0;
    for (const auto& system : result.systems) {
        candidateLengthSum += system.totalLength;
    }
    const double minPhysicalLength = std::max(
        kActiveSystemMinLengthMeters,
        kActiveSystemMinPerimeterShare * totalPerimeter);
    int activeCount = 0;
    for (std::size_t i = 0; i < result.systems.size(); ++i) {
        DetectedDirectionSystem& system = result.systems[i];
        const double weightShare = totalWeight > 1e-9
            ? system.weight / totalWeight : 0.0;
        const double lengthShare = candidateLengthSum > 1e-9
            ? system.totalLength / candidateLengthSum : 0.0;
        const char* reason = "insufficient_support";
        if (i == 0) {
            system.active = true;
            reason = "primary";
        } else {
            const bool enoughLength =
                system.totalLength >= minPhysicalLength;
            const bool enoughSupport =
                weightShare >= kActiveSystemMinWeightShare ||
                lengthShare >= kActiveSystemMinLengthShare;
            system.active = enoughLength && enoughSupport;
            if (system.active) {
                reason = weightShare >= kActiveSystemMinWeightShare
                    ? "accepted_weight" : "accepted_length";
            }
        }
        if (system.active) ++activeCount;
        if (fid >= 0) {
            std::cerr << "[DirectionActivation] fid=" << fid
                      << " part=" << partIdx
                      << " rank=" << (i + 1)
                      << " angle_deg=" << system.angleRad * 180.0 / M_PI
                      << " weight_share=" << weightShare
                      << " length_share=" << lengthShare
                      << " total_length=" << system.totalLength
                      << " min_length=" << minPhysicalLength
                      << " active=" << (system.active ? 1 : 0)
                      << " reason=" << reason << std::endl;
        }
    }

    result.valid = true;
    result.primaryAngle = result.systems.front().angleRad;
    result.concentration = result.systems.front().concentration;
    result.multiDirection = activeCount >= 2;

    if (fid >= 0) {
        std::cerr << "[DirectionDetect] fid=" << fid << " part=" << partIdx
                  << " valid=" << (result.valid ? 1 : 0)
                  << " candidate_systems=" << result.systems.size()
                  << " active_systems=" << activeCount
                  << " multi=" << (result.multiDirection ? 1 : 0)
                  << " primary_deg=" << result.primaryAngle * 180.0 / M_PI
                  << " chains=" << chains.size()
                  << " perimeter=" << totalPerimeter
                  << " vertices=" << smoothRing.size()
                  << " bandwidth_deg=" << kKdeBandwidthDeg
                  << " prom_cut=" << kPeakProminenceRatio;
        if (!result.valid) std::cerr << " reject=" << result.rejectReason;
        std::cerr << std::endl;
        for (std::size_t i = 0; i < result.systems.size(); ++i) {
            const auto& system = result.systems[i];
            std::cerr << "[DirectionDetect] fid=" << fid
                      << " system=" << i
                      << " angle_deg=" << system.angleRad * 180.0 / M_PI
                      << " edges=" << system.chainCount
                      << " length=" << system.totalLength
                      << " avg_edge=" << system.avgEdgeLen
                      << " concentration=" << system.concentration
                      << " prom=" << system.prominence
                      << " confidence=" << system.confidence << std::endl;
        }
    }
    return result;
}

int AssignEdgeToDirectionSystem(
    double edgeAngle,
    const std::vector<double>& systemAngles,
    double maxAssignDeg) {
    if (systemAngles.empty()) return -1;
    const double tolerance = maxAssignDeg * M_PI / 180.0;
    int best = -1;
    double bestDistance = tolerance;
    for (std::size_t i = 0; i < systemAngles.size(); ++i) {
        const double distance = foldedDistance90(edgeAngle, systemAngles[i]);
        if (distance < bestDistance) {
            best = static_cast<int>(i);
            bestDistance = distance;
        }
    }
    return best;
}

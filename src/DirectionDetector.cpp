// =============================================================
// DirectionDetector.cpp
// 统一主方向检测: 两阶段贪心聚类 + 加权KDE
// 从 outlineRegular.cpp 的 BuildDirectionSystems 提取,
// 作为独立模块供所有规则化路径使用。
// =============================================================
#include "DirectionDetector.h"
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <limits>

// ---- 内部常量(与原 BuildDirectionSystems 保持一致) ----
namespace {

constexpr double kStableMinLength = 1.5;        // m, 稳定链长度下限
constexpr double kSystemMinSeparationDeg = 20.0; // 系统间最小角距(低于=同方向栅格抖动)
constexpr double kAssignDeg = 16.0;              // 稳定链归组阈值
constexpr double kShortAssignDeg = 16.0;         // 短链归组阈值
constexpr double kLongRescueDeg = 20.0;          // 未归组长链抢救上限
constexpr double kMergeDeg = 14.0;               // 精化后角距过近的系统合并
constexpr int kMaxCandidates = 7;                // 候选方向上限
constexpr int kMaxFinalSystems = 5;              // 最终系统上限
constexpr double kSystemMinChainFrac = 0.15;     // 多链系统最小支持长度占比
constexpr double kSeedMinLength = 3.0;           // 单链成系统的最小长度
constexpr double kKdeBandwidthDeg = 4.0;         // KDE 平滑带宽
constexpr double kCertaintyConf = 0.15;          // 主系统置信度下限
constexpr double kCertaintyScore = 0.30;         // 单方向评分下限
constexpr double kUnassignedHighRatio = 0.10;    // 未归组占比>此值→不确定
constexpr double kUnassignedWeakRatio = 0.05;    // 弱未归组阈值

} // namespace

// ---- 折叠角工具函数 ----
static double foldedAngle90(double angle) {
    angle = std::fmod(angle, M_PI / 2.0);
    if (angle < 0) angle += M_PI / 2.0;
    return angle;
}

static double foldedDistance90(double a, double b) {
    const double period = M_PI / 2.0;
    double diff = std::fmod(std::abs(a - b), period);
    return std::min(diff, period - diff);
}

static double circularMean90(const std::vector<double>& angles,
                              const std::vector<double>& weights) {
    if (angles.empty()) return 0.0;
    double x = 0.0, y = 0.0, total = 0.0;
    for (std::size_t i = 0; i < angles.size() && i < weights.size(); ++i) {
        x += weights[i] * std::cos(4.0 * angles[i]);
        y += weights[i] * std::sin(4.0 * angles[i]);
        total += weights[i];
    }
    if (total <= 1e-9 || (x * x + y * y) <= 1e-12 * total * total) {
        return angles.front();
    }
    return foldedAngle90(0.25 * std::atan2(y, x));
}

// ---- 内部链结构(与 outlineRegular 的 DirectionChain 兼容) ----
struct DetChain {
    double angleRad = 0.0;
    double length = 0.0;
    double weight = 0.0;
    double rmse = 0.0;
    bool stable = false;
    int system = -1;
};

// ---- KDE 峰检测 ----
struct DetPeak {
    double angleRad = 0.0;
};

static std::vector<DetPeak> detectPeaks(const std::vector<DetChain>& chains) {
    std::vector<DetPeak> peaks;
    double totalWeight = 0.0;
    for (const auto& c : chains) totalWeight += c.weight;
    if (totalWeight < 1e-9) return peaks;

    constexpr int kBins = 180;
    std::vector<double> kde(kBins, 0.0);
    const double bw = kKdeBandwidthDeg * M_PI / 180.0;
    for (int b = 0; b < kBins; ++b) {
        const double angle = (b + 0.5) * 0.5 * M_PI / 180.0;
        double density = 0.0;
        for (const auto& c : chains) {
            double diff = std::abs(angle - c.angleRad);
            if (diff > M_PI / 4.0) diff = M_PI / 2.0 - diff;
            density += c.weight * std::exp(-0.5 * (diff / bw) * (diff / bw));
        }
        kde[b] = density;
    }

    // 循环局部极大值(0°/90° 边界正确处理)
    std::vector<int> maxima;
    for (int b = 0; b < kBins; ++b) {
        const int prev = (b + kBins - 1) % kBins;
        const int next = (b + 1) % kBins;
        if (kde[b] > kde[prev] && kde[b] >= kde[next] && kde[b] > totalWeight * 0.02) {
            maxima.push_back(b);
        }
    }
    std::sort(maxima.begin(), maxima.end(),
              [&](int a, int b) { return kde[a] > kde[b]; });
    // 分离约束
    constexpr double minSepBins = 12.0 / 0.5; // 12° 间隔
    std::vector<int> selected;
    for (int m : maxima) {
        bool ok = true;
        for (int s : selected) {
            double d = std::abs(m - s);
            if (d > kBins / 2.0) d = kBins - d;
            if (d < minSepBins) { ok = false; break; }
        }
        if (ok) selected.push_back(m);
        if ((int)selected.size() >= 4) break;
    }
    for (int bin : selected) {
        peaks.push_back({(bin + 0.5) * 0.5 * M_PI / 180.0});
    }
    return peaks;
}

// =============================================================
// 主检测函数
// =============================================================
DetectedDirectionResult DetectBuildingDirection(
    const std::vector<pcl::PointXYZ>& smoothRing,
    const std::vector<pcl::PointXYZ>& rawRing,
    double pixelSize,
    long long fid,
    int partIdx)
{
    DetectedDirectionResult result;
    const std::size_t n = smoothRing.size();
    if (n < 8) {
        result.rejectReason = "too_few_vertices";
        return result;
    }

    // ---- 1. 简化链提取(轻量版: 直接从平滑环顶点提取) ----
    // 不调用完整的 ExtractDirectionSystems, 而是直接计算每条边的方向
    // (方向检测只需要角度, 不需要完整的split-and-merge)
    std::vector<DetChain> chains;
    chains.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = smoothRing[i];
        const auto& b = smoothRing[(i + 1) % n];
        const double len = std::hypot(b.x - a.x, b.y - a.y);
        if (len < 0.1) continue;
        DetChain c;
        c.angleRad = foldedAngle90(std::atan2(b.y - a.y, b.x - a.x));
        c.length = len;
        c.weight = len; // 无 rmse 信息时用纯长度
        c.rmse = 0.0;
        c.stable = len >= kStableMinLength;
        chains.push_back(c);
    }

    // 如果有 raw 环且顶点更多, 用 raw 环的边做更精确的角度投票
    if (!rawRing.empty() && rawRing.size() > n) {
        std::vector<DetChain> rawChains;
        rawChains.reserve(rawRing.size());
        for (std::size_t i = 0; i < rawRing.size(); ++i) {
            const auto& a = rawRing[i];
            const auto& b = rawRing[(i + 1) % rawRing.size()];
            const double len = std::hypot(b.x - a.x, b.y - a.y);
            if (len < 0.05) continue;
            DetChain c;
            c.angleRad = foldedAngle90(std::atan2(b.y - a.y, b.x - a.x));
            c.length = len;
            c.weight = len * len; // raw 边很短, 平方权重让长墙主导
            c.rmse = 0.0;
            c.stable = false; // raw 边不单独作为稳定链
            rawChains.push_back(c);
        }
        // 将 raw 链的方向投票加入(权重很小, 主要由平滑环主导)
        for (auto& rc : rawChains) {
            rc.weight *= 0.1;
            chains.push_back(rc);
        }
    }

    if (chains.empty()) {
        result.rejectReason = "no_valid_edges";
        return result;
    }
    result.totalChains = (int)chains.size();
    result.stableChains = 0;
    for (const auto& c : chains) {
        if (c.stable) {
            result.stableChains++;
            result.totalStableLength += c.length;
        }
    }

    // ---- 2. KDE 峰检测 ----
    // 只用稳定链的方向投票
    std::vector<DetChain> stableOnly;
    for (const auto& c : chains) {
        if (c.stable) stableOnly.push_back(c);
    }
    if (stableOnly.empty()) {
        // 没有足够长的链 → 用所有链
        stableOnly = chains;
    }
    const auto peaks = detectPeaks(stableOnly);

    // ---- 3. 两阶段贪心聚类 ----
    // 阶段一: 种子发现(全部稳定链按权重降序)
    std::vector<DetChain>& allChains = chains;
    const std::size_t nc = allChains.size();
    const double minSep = kSystemMinSeparationDeg * M_PI / 180.0;

    // 按权重排序的稳定链索引
    std::vector<std::size_t> sortedStable;
    for (std::size_t i = 0; i < nc; ++i) {
        if (allChains[i].stable) sortedStable.push_back(i);
    }
    std::sort(sortedStable.begin(), sortedStable.end(),
              [&](std::size_t a, std::size_t b) {
                  return allChains[a].weight > allChains[b].weight;
              });

    std::vector<double> seedAngles;
    for (std::size_t si : sortedStable) {
        if ((int)seedAngles.size() >= kMaxCandidates) break;
        const double ang = allChains[si].angleRad;
        bool near = false;
        for (double s : seedAngles) {
            if (foldedDistance90(ang, s) < minSep) { near = true; break; }
        }
        if (!near) seedAngles.push_back(ang);
    }
    if (seedAngles.empty()) {
        result.rejectReason = "no_seeds";
        return result;
    }

    // 阶段二: 统一就近分配 + 验证迭代
    std::vector<int> systemOf(nc, -1);
    std::vector<double> sysAngles = seedAngles;

    auto assignAll = [&]() {
        std::fill(systemOf.begin(), systemOf.end(), -1);
        if (sysAngles.empty()) return;
        const double assignRad = kAssignDeg * M_PI / 180.0;
        const double rescueRad = kLongRescueDeg * M_PI / 180.0;
        for (std::size_t i = 0; i < nc; ++i) {
            if (!allChains[i].stable) continue;
            std::size_t nearest = 0;
            double best = foldedDistance90(allChains[i].angleRad, sysAngles[0]);
            for (std::size_t s = 1; s < sysAngles.size(); ++s) {
                const double d = foldedDistance90(allChains[i].angleRad, sysAngles[s]);
                if (d < best) { best = d; nearest = s; }
            }
            if (best <= assignRad) {
                systemOf[i] = (int)nearest;
            } else if (allChains[i].length >= 3.0 && best <= rescueRad) {
                systemOf[i] = (int)nearest;
            }
        }
    };

    auto recomputeAngles = [&]() {
        for (std::size_t s = 0; s < sysAngles.size(); ++s) {
            std::vector<double> ma, mw;
            for (std::size_t i = 0; i < nc; ++i) {
                if (systemOf[i] == (int)s && allChains[i].stable) {
                    ma.push_back(allChains[i].angleRad);
                    mw.push_back(allChains[i].weight);
                }
            }
            if (!ma.empty()) sysAngles[s] = circularMean90(ma, mw);
        }
    };

    // 删除弱系统(迭代≤3轮)
    double totalStableLen = result.totalStableLength;
    for (int iter = 0; iter < 3; ++iter) {
        assignAll();
        recomputeAngles();
        // 验证每个系统
        bool removed = false;
        for (std::size_t s = sysAngles.size(); s-- > 0;) {
            double lenSum = 0.0;
            int cnt = 0;
            for (std::size_t i = 0; i < nc; ++i) {
                if (systemOf[i] == (int)s && allChains[i].stable) {
                    lenSum += allChains[i].length;
                    ++cnt;
                }
            }
            bool credible;
            if (cnt >= 2) {
                credible = lenSum >= kSystemMinChainFrac * std::max(totalStableLen, 1e-9);
            } else if (cnt == 1) {
                credible = lenSum >= std::max(kSeedMinLength * 5.0, 0.10 * (totalStableLen + result.totalStableLength));
            } else {
                credible = false;
            }
            if (!credible) {
                // 删除系统 s, 重映射索引
                sysAngles.erase(sysAngles.begin() + s);
                for (auto& si : systemOf) {
                    if (si == (int)s) si = -1;
                    else if (si > (int)s) --si;
                }
                removed = true;
            }
        }
        if (!removed) break;
    }

    // 精化后角距过近的系统合并
    {
        const double mergeRad = kMergeDeg * M_PI / 180.0;
        bool merged = true;
        while (merged && sysAngles.size() >= 2) {
            merged = false;
            std::size_t bi = 0, bj = 1;
            double bestD = 1e9;
            for (std::size_t a = 0; a < sysAngles.size(); ++a) {
                for (std::size_t b = a + 1; b < sysAngles.size(); ++b) {
                    const double d = foldedDistance90(sysAngles[a], sysAngles[b]);
                    if (d < bestD) { bestD = d; bi = a; bj = b; }
                }
            }
            if (bestD < mergeRad) {
                const std::size_t donor = (allChains[0].weight > 0) ? bj : bi;
                const std::size_t keeper = (donor == bi) ? bj : bi;
                for (auto& si : systemOf) {
                    if (si == (int)donor) si = (int)keeper;
                    else if (si > (int)donor) --si;
                }
                sysAngles.erase(sysAngles.begin() + donor);
                assignAll();
                recomputeAngles();
                merged = true;
            }
        }
    }

    // ---- 4. 构建输出 ----
    assignAll();
    recomputeAngles();

    // 未归组长链统计
    double unassignedLen = 0.0, unassignedW = 0.0;
    double totalW = 0.0;
    for (std::size_t i = 0; i < nc; ++i) {
        if (!allChains[i].stable) continue;
        totalW += allChains[i].weight;
        if (systemOf[i] < 0 && allChains[i].length >= 3.0) {
            unassignedLen += allChains[i].length;
            unassignedW += allChains[i].weight;
        }
    }
    const double unassignedLenRatio =
        totalStableLen > 1e-9 ? unassignedLen / totalStableLen : 0.0;
    const double unassignedWRatio = totalW > 1e-9 ? unassignedW / totalW : 0.0;

    // 填充系统信息
    for (std::size_t s = 0; s < sysAngles.size() && s < (size_t)kMaxFinalSystems; ++s) {
        DetectedDirectionSystem ds;
        ds.angleRad = sysAngles[s];
        double ma[64], mw[64];
        int cnt = 0;
        for (std::size_t i = 0; i < nc && cnt < 64; ++i) {
            if (systemOf[i] == (int)s && allChains[i].stable) {
                ma[cnt] = allChains[i].angleRad;
                mw[cnt] = allChains[i].weight;
                ds.totalLength += allChains[i].length;
                ds.weight += allChains[i].weight;
                ++cnt;
            }
        }
        ds.chainCount = cnt;
        if (cnt > 0) {
            // 集中度
            double x4 = 0, y4 = 0, w4 = 0;
            for (int i = 0; i < cnt; ++i) {
                x4 += mw[i] * std::cos(4.0 * ma[i]);
                y4 += mw[i] * std::sin(4.0 * ma[i]);
                w4 += mw[i];
            }
            ds.concentration = w4 > 1e-9 ? std::hypot(x4, y4) / w4 : 0.0;
            // 置信度
            const double frac = totalW > 1e-9 ? ds.weight / totalW : 0.0;
            ds.confidence = std::min(1.0, frac * ds.concentration);
        }
        result.systems.push_back(ds);
    }

    if (result.systems.empty()) {
        result.rejectReason = "no_credible_systems";
        return result;
    }

    result.primaryAngle = result.systems[0].angleRad;
    result.concentration = result.systems[0].concentration;

    // 多方向判定
    if (result.systems.size() >= 2) {
        const auto& second = result.systems[1];
        const double secondFrac = totalW > 1e-9 ? second.weight / totalW : 0.0;
        result.multiDirection = secondFrac >= 0.20 && second.chainCount >= 2;
    }

    // 方向确定性
    if (result.systems[0].confidence < kCertaintyConf) {
        result.rejectReason = "low_confidence";
        result.valid = false;
    } else if (unassignedLenRatio > kUnassignedHighRatio ||
               unassignedWRatio > kUnassignedHighRatio) {
        result.rejectReason = "unassigned_long_chain";
        result.valid = false;
    } else {
        result.valid = true;
    }

    // 日志
    if (fid >= 0) {
        std::cerr << "[DirectionDetect] fid=" << fid << " part=" << partIdx
                  << " valid=" << (result.valid ? 1 : 0)
                  << " systems=" << result.systems.size()
                  << " multi=" << (result.multiDirection ? 1 : 0)
                  << " primary_deg=" << result.primaryAngle * 180.0 / M_PI
                  << " concentration=" << result.concentration
                  << " stable=" << result.stableChains
                  << " stable_len=" << result.totalStableLength
                  << " unassigned_ratio=" << unassignedLenRatio;
        if (!result.valid) std::cerr << " reject=" << result.rejectReason;
        std::cerr << std::endl;
        for (std::size_t s = 0; s < result.systems.size(); ++s) {
            const auto& sys = result.systems[s];
            std::cerr << "[DirectionDetect] fid=" << fid
                      << " system=" << s
                      << " angle_deg=" << sys.angleRad * 180.0 / M_PI
                      << " chains=" << sys.chainCount
                      << " length=" << sys.totalLength
                      << " concentration=" << sys.concentration
                      << " confidence=" << sys.confidence << std::endl;
        }
    }

    return result;
}

// =============================================================
// 辅助: 边到最近方向系统的分配
// =============================================================
int AssignEdgeToDirectionSystem(
    double edgeAngle,
    const std::vector<double>& systemAngles,
    double maxAssignDeg)
{
    if (systemAngles.empty()) return -1;
    const double tol = maxAssignDeg * M_PI / 180.0;
    int best = -1;
    double bestD = tol;
    for (std::size_t s = 0; s < systemAngles.size(); ++s) {
        const double d = foldedDistance90(edgeAngle, systemAngles[s]);
        if (d < bestD) { bestD = d; best = (int)s; }
    }
    return best;
}

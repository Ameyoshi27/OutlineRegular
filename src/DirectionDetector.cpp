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
constexpr double kIndependentPeakMinSeparationDeg = 7.0;
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
    double length = 0.0;    // 路径长度
    double chordLength = 0.0;
    double weight = 0.0;    // length^kLenExponent
    double rmse = 0.0;
    double maxDeviation = 0.0;
    double continuity = 1.0;  // chord/path
    double centerX = 0.0;
    double centerY = 0.0;
    int edgeCount = 0;
    // ---- 只读诊断字段(救援结构分析用, 不参与任何计算) ----
    std::size_t startEdgeIndex = 0;      // 链起点在环上的顶点索引
    std::size_t endEdgeIndex = 0;        // 链终点在环上的顶点索引
    double rawDirectedAngleRad = 0.0;    // 首尾弦有向角 [-π,π)
    int ringOrder = 0;                   // 输出序号(环序)
};

// 正交最小二乘直线拟合(TLS): 返回主轴方向与偏差统计。
// 方向折叠到 [0, π/2)。
struct TlsLineFit {
    double angleRad = 0.0;
    double rmse = 0.0;
    double maxDeviation = 0.0;
    bool valid = false;
};
TlsLineFit fitLineTLS(const pcl::PointXYZ* points, std::size_t count) {
    TlsLineFit fit;
    if (count < 2) return fit;
    double cx = 0.0, cy = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
        cx += points[i].x; cy += points[i].y;
    }
    cx /= count; cy /= count;
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
        const double dx = points[i].x - cx;
        const double dy = points[i].y - cy;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    // 2x2 协方差主轴: λmax 的特征向量
    const double halfTrace = 0.5 * (sxx + syy);
    const double disc = std::sqrt(
        0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy);
    const double lambdaMax = halfTrace + disc;
    double dirX = 0.0, dirY = 0.0;
    if (std::abs(sxy) > 1e-12) {
        dirX = sxy; dirY = lambdaMax - sxx;
    } else if (sxx >= syy) {
        dirX = 1.0; dirY = 0.0;
    } else {
        dirX = 0.0; dirY = 1.0;
    }
    const double norm = std::hypot(dirX, dirY);
    if (norm < 1e-12) return fit;
    dirX /= norm; dirY /= norm;
    double sumSq = 0.0, maxDev = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
        const double dx = points[i].x - cx;
        const double dy = points[i].y - cy;
        const double dev = std::abs(-dirY * dx + dirX * dy);
        sumSq += dev * dev;
        maxDev = std::max(maxDev, dev);
    }
    fit.angleRad = foldedAngle90(std::atan2(dirY, dirX));
    fit.rmse = std::sqrt(sumSq / count);
    fit.maxDeviation = maxDev;
    fit.valid = true;
    return fit;
}

// 全局形状 PCA: 等弧长采样 + 协方差主轴。
// 简化后顶点密度不均, 顶点等权 PCA 会被密集区带偏,
// 等弧长采样恢复边界长度的话语权。
struct ShapePcaDirection {
    bool valid = false;
    double angleRad = 0.0;
    double axisRatio = 1.0;
    double anisotropy = 0.0;
};
ShapePcaDirection computeShapePca(
    const std::vector<pcl::PointXYZ>& ring, double pixelSize) {
    ShapePcaDirection pca;
    if (ring.size() < 4) return pca;
    const double spacing = std::clamp(pixelSize, 0.3, 0.8);
    std::vector<double> sx, sy;  // 等弧长采样点
    sx.reserve(static_cast<std::size_t>(600.0 / spacing));
    sy.reserve(sx.capacity());
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto& a = ring[i];
        const auto& b = ring[(i + 1) % ring.size()];
        const double len = std::hypot(b.x - a.x, b.y - a.y);
        if (len < 1e-6) continue;  // 跳过零长边(含闭合重复点)
        const int steps = std::max(1, static_cast<int>(len / spacing));
        for (int k = 0; k < steps; ++k) {
            const double t = static_cast<double>(k) / steps;
            sx.push_back(a.x + t * (b.x - a.x));
            sy.push_back(a.y + t * (b.y - a.y));
        }
    }
    if (sx.size() < 8) return pca;
    double cx = 0.0, cy = 0.0;
    for (std::size_t i = 0; i < sx.size(); ++i) {
        cx += sx[i]; cy += sy[i];
    }
    cx /= sx.size(); cy /= sy.size();
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < sx.size(); ++i) {
        const double dx = sx[i] - cx;
        const double dy = sy[i] - cy;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    const double halfTrace = 0.5 * (sxx + syy);
    const double disc = std::sqrt(
        0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy);
    const double lambdaMax = halfTrace + disc;
    const double lambdaMin = halfTrace - disc;
    if (lambdaMax < 1e-9 || lambdaMin < 1e-9) return pca;
    double dirX = 0.0, dirY = 0.0;
    if (std::abs(sxy) > 1e-12) {
        dirX = sxy; dirY = lambdaMax - sxx;
    } else if (sxx >= syy) {
        dirX = 1.0; dirY = 0.0;
    } else {
        dirX = 0.0; dirY = 1.0;
    }
    const double norm = std::hypot(dirX, dirY);
    if (norm < 1e-12) return pca;
    pca.angleRad = foldedAngle90(std::atan2(dirY / norm, dirX / norm));
    pca.axisRatio = lambdaMax / std::max(lambdaMin, 1e-9);
    pca.anisotropy = (lambdaMax - lambdaMin) /
        std::max(lambdaMax + lambdaMin, 1e-9);
    pca.valid = true;
    return pca;
}

// ---- 边链化: 空间连续 + 转向平缓 + 整体直线性三重约束 ----
// 平滑环上斜墙是方向渐变的短边串, 逐边投票会在角度空间被摊薄,
// 而精确对齐栅格轴的楼梯边高度集中, 密度峰被伪方向占据。
// 以链为单位投票可恢复墙体的真实长度话语权;
// 整体直线性约束(TLS 拟合偏差 + 弦路比)阻止弯曲边串被合并成
// 错误长斜弦后获得 length^1.5 高权重。
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
    // 直线性容差: 超限的弯曲链不切断(保留方向证据), 但权重按
    // 弯曲度指数降权——弯曲链(渐变墙/弧段)不能获得 length^1.5
    // 的完整话语权, 否则其首尾弦方向会虚假主导 KDE。
    const double maxDevLimit = std::max(0.35, 1.2 * pixelSize);

    std::vector<DetChain> chains;
    std::size_t chainFromVertex = 0;
    std::size_t chainVertexCount = 0;  // 链覆盖的顶点数(含起点)
    double chainDir = 0.0;
    double chainLen = 0.0;
    int chainEdges = 0;

    // 环序收集链覆盖顶点 [chainFromVertex, 环序 count 个]
    auto collectVertices = [&](std::size_t count,
                               std::vector<pcl::PointXYZ>& out) {
        out.clear();
        out.reserve(count);
        for (std::size_t k = 0; k < count; ++k) {
            out.push_back(ring[(chainFromVertex + k) % ring.size()]);
        }
    };
    std::vector<pcl::PointXYZ> chainPoints;

    auto closeChain = [&](std::size_t endVertex) {
        if (chainEdges == 0) return;
        const auto& a = ring[chainFromVertex];
        const auto& b = ring[endVertex];
        const double chord = std::hypot(b.x - a.x, b.y - a.y);
        if (chord < 1e-6 || chainLen < minChainLength) return;
        DetChain chain;
        chain.length = chainLen;
        chain.chordLength = chord;
        chain.centerX = 0.5 * (a.x + b.x);
        chain.centerY = 0.5 * (a.y + b.y);
        chain.edgeCount = chainEdges;
        chain.continuity = chord / std::max(chainLen, 1e-9);
        TlsLineFit fit;
        if (chainVertexCount >= 3) {
            collectVertices(chainVertexCount, chainPoints);
            fit = fitLineTLS(chainPoints.data(), chainPoints.size());
        }
        if (fit.valid) {
            // TLS 主轴方向比首尾弦更接近真实墙方向(渐变段对称轴)
            chain.angleRad = fit.angleRad;
            chain.rmse = fit.rmse;
            chain.maxDeviation = fit.maxDeviation;
        } else {
            chain.angleRad = foldedAngle90(std::atan2(b.y - a.y, b.x - a.x));
        }
        // 弯曲度惩罚: maxDev 在容差内全额权重, 超出按指数衰减
        const double excess = std::max(0.0, chain.maxDeviation / maxDevLimit - 1.0);
        const double linearityFactor = std::exp(-excess);
        chain.weight = std::pow(chainLen, kLenExponent) * linearityFactor;
        // 只读诊断: 环序索引/有向角/序号(不参与计算)
        chain.startEdgeIndex = chainFromVertex;
        chain.endEdgeIndex = endVertex;
        chain.rawDirectedAngleRad = std::atan2(b.y - a.y, b.x - a.x);
        chain.ringOrder = static_cast<int>(chains.size());
        chains.push_back(chain);
        chainEdges = 0;
    };

    for (std::size_t step = 0; step < edgeList.size(); ++step) {
        const std::size_t idx = (startEdge + step) % edgeList.size();
        const auto& edge = edgeList[idx];
        if (chainEdges == 0) {
            chainFromVertex = edge.from;
            chainVertexCount = 2;  // 起点 + 本边终点
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
            ++chainVertexCount;
        } else {
            // 转向超限: 关闭当前链(终点=新边起点), 新链从当前边开始
            closeChain(edge.from);
            chainFromVertex = edge.from;
            chainVertexCount = 2;
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
        bool evidenceBacked = false;       // 直墙证据救回(prominence 不达标)
        const char* acceptReason = "prominence";
        double refinedAngleRad = std::numeric_limits<double>::quiet_NaN();
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
    const double independentMinSeparation =
        kIndependentPeakMinSeparationDeg * M_PI / 180.0;
    const double independentAngleTolerance = 3.0 * M_PI / 180.0;
    const double straightDeviationLimit = std::max(0.35, 1.2 * pixelSize);
    const double independentLengthLimit = std::max(6.0, 0.04 * totalPerimeter);
    std::vector<Peak> accepted;
    std::vector<bool> independentlySeparated;
    // ---- 直墙支持评估: 峰 ±5° 内高质量直链统计(救援与角度精化共用) ----
    struct StraightSupport {
        int straightCount = 0;
        double straightLength = 0.0;
        double longest = 0.0;
        const DetChain* longestChain = nullptr;
        double x4 = 0.0, y4 = 0.0, weightSum = 0.0;  // 4θ 圆均值累积
        int supportRuns = 0;              // 支持链沿环序的连续 run 数(环形)
        double longestRun = 0.0;          // 最长 run 的链长和
    };
    auto evaluateStraightSupport = [&](double peakAngle) {
        StraightSupport s;
        std::vector<const DetChain*> supportChains;
        for (const auto& chain : chains) {
            if (foldedDistance90(chain.angleRad, peakAngle) >
                5.0 * M_PI / 180.0) continue;
            // 高质量直链: 弦路比 + 垂直偏差(随弦长放宽) + 最小长度
            const double devLimit = std::max(0.60, 0.02 * chain.chordLength);
            const bool straight = chain.continuity >= 0.97 &&
                chain.maxDeviation <= devLimit &&
                chain.length >= std::max(2.0, 3.0 * pixelSize);
            if (!straight) continue;
            ++s.straightCount;
            s.straightLength += chain.length;
            if (chain.length > s.longest) {
                s.longest = chain.length;
                s.longestChain = &chain;
            }
            // 精化权重 = 投票权重(len^1.5 × 直线性折扣)
            s.x4 += chain.weight * std::cos(4.0 * chain.angleRad);
            s.y4 += chain.weight * std::sin(4.0 * chain.angleRad);
            s.weightSum += chain.weight;
            supportChains.push_back(&chain);
        }
        // 环序 run 统计: chains 输出序即 ringOrder, 连续支持链成 run,
        // 首尾 run(ringOrder 0 与 N-1 相邻)合并
        if (!supportChains.empty()) {
            const std::size_t total = chains.size();
            std::vector<char> isSup(total, 0);
            for (const DetChain* c : supportChains) {
                isSup[static_cast<std::size_t>(c->ringOrder)] = 1;
            }
            std::size_t start = 0;
            while (start < total && isSup[start]) {
                std::size_t p = (start + total - 1) % total;
                bool wrapped = false;
                for (std::size_t g = 0; g < total && isSup[p]; ++g) {
                    start = p;
                    p = (p + total - 1) % total;
                    wrapped = true;
                }
                if (!wrapped || start == 0) break;
            }
            bool inRun = false;
            double runLen = 0.0;
            for (std::size_t step = 0; step <= total; ++step) {
                const std::size_t i = (start + step) % total;
                if (step < total && isSup[i]) {
                    runLen += chains[i].length;
                    inRun = true;
                } else if (inRun) {
                    ++s.supportRuns;
                    s.longestRun = std::max(s.longestRun, runLen);
                    inRun = false;
                    runLen = 0.0;
                }
            }
        }
        return s;
    };
    // 主峰最强支持链(空间独立性参照)
    const StraightSupport mainSupport =
        evaluateStraightSupport(peaks.front().bin * binStep);
    const double rescueLengthGate = std::max(8.0, 0.04 * totalPerimeter);
    const double rescueLongestGate = std::max(10.0, 0.04 * totalPerimeter);

    int rejectedDiag = 0;
    for (auto& peak : peaks) {
        if (static_cast<int>(accepted.size()) >= kMaxFinalSystems) break;
        peak.prominence = circularProminence(density, peak.bin);
        if (!accepted.empty() &&
            peak.prominence < kPeakProminenceRatio * mainHeight) {
            // ---- 直墙证据救援 ----
            // 不降全局 prominence; 只有"高度比够 + 与已接受峰距离够 +
            // 独立直墙支撑够"的峰被救回。低而宽的锯齿隆起依旧被拒。
            const double heightRatio =
                peak.height / std::max(mainHeight, 1e-12);
            bool farEnough = true;
            for (const auto& kept : accepted) {
                if (foldedDistance90(peak.bin * binStep, kept.angleRad) <
                    12.0 * M_PI / 180.0) {
                    farEnough = false;
                    break;
                }
            }
            StraightSupport support;
            const char* rescueReason = nullptr;
            if (heightRatio < 0.12) {
                rescueReason = "insufficient_height";
            } else if (!farEnough) {
                rescueReason = "too_close";
            } else {
                support = evaluateStraightSupport(peak.bin * binStep);
                const bool enough = (support.straightCount >= 2 &&
                                     support.straightLength >= rescueLengthGate) ||
                    (support.straightCount >= 1 &&
                     support.longest >= rescueLongestGate);
                if (!enough) {
                    rescueReason = "insufficient_straight_support";
                } else if (mainSupport.longestChain && support.longestChain) {
                    // 空间独立: 与主峰最强支持链不是同一段墙
                    const double centerGap = std::hypot(
                        support.longestChain->centerX -
                            mainSupport.longestChain->centerX,
                        support.longestChain->centerY -
                            mainSupport.longestChain->centerY);
                    const double requiredGap = std::max(3.0, 0.15 *
                        (support.longestChain->chordLength +
                         mainSupport.longestChain->chordLength));
                    if (centerGap < requiredGap) {
                        rescueReason = "not_spatially_independent";
                    }
                }
            }
            const bool rescued = rescueReason == nullptr;
            if (rescued) {
                peak.evidenceBacked = true;
                peak.acceptReason = "straight_wall_rescue";
                // 角度精化: 救回峰不用 KDE bin 的粗位置, 用 ±5° 直链
                // 4θ 加权圆均值(折叠空间周期 90°, 必须四倍角)
                if (support.weightSum > 1e-9 &&
                    std::hypot(support.x4, support.y4) >
                        1e-12 * support.weightSum) {
                    peak.refinedAngleRad = foldedAngle90(
                        0.25 * std::atan2(support.y4, support.x4));
                }
            }
            if (fid >= 0 && (rescued || rejectedDiag < 3)) {
                if (!rescued) ++rejectedDiag;
                std::cerr << "[DirectionPeakFilter] fid=" << fid
                          << " part=" << partIdx
                          << " raw_deg=" << peak.bin * binStep * 180.0 / M_PI
                          << " height_ratio=" << heightRatio
                          << " prominence_ratio=" << peak.prominence /
                                 std::max(mainHeight, 1e-12)
                          << " straight_chains=" << support.straightCount
                          << " straight_length=" << support.straightLength
                          << " longest_chain=" << support.longest
                          << " spatially_independent="
                          << (rescued ? 1 : 0)
                          << " accepted=" << (rescued ? 1 : 0)
                          << " reason=" << (rescued
                              ? "straight_wall_rescue" : rescueReason)
                          << std::endl;
            }
            if (!rescued) continue;
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
        if (peak.evidenceBacked &&
            std::isfinite(peak.refinedAngleRad)) {
            peak.angleRad = peak.refinedAngleRad;
        }
        bool tooClose = false;
        bool independentClosePeak = false;
        for (const auto& kept : accepted) {
            const double peakGap = foldedDistance90(peak.angleRad, kept.angleRad);
            if (peakGap >= minSeparation) continue;

            // A 7-10 degree pair is retained only when each peak owns a
            // different long, straight wall chain in a different location.
            const DetChain* peakSupport = nullptr;
            const DetChain* keptSupport = nullptr;
            for (const auto& chain : chains) {
                const double dPeak = foldedDistance90(chain.angleRad, peak.angleRad);
                const double dKept = foldedDistance90(chain.angleRad, kept.angleRad);
                const double chainDeviationLimit = std::max(
                    straightDeviationLimit, 0.02 * chain.chordLength);
                const bool straight = chain.maxDeviation <= chainDeviationLimit &&
                    chain.continuity >= 0.97;
                if (!straight || chain.length < independentLengthLimit) continue;
                if (dPeak <= independentAngleTolerance && dPeak + 1e-6 < dKept &&
                    (!peakSupport || chain.length > peakSupport->length)) {
                    peakSupport = &chain;
                }
                if (dKept <= independentAngleTolerance && dKept + 1e-6 < dPeak &&
                    (!keptSupport || chain.length > keptSupport->length)) {
                    keptSupport = &chain;
                }
            }
            bool spatiallyIndependent = false;
            if (peakSupport && keptSupport) {
                const double centerGap = std::hypot(
                    peakSupport->centerX - keptSupport->centerX,
                    peakSupport->centerY - keptSupport->centerY);
                const double requiredGap = std::max(
                    3.0, 0.15 * (peakSupport->chordLength + keptSupport->chordLength));
                spatiallyIndependent = centerGap >= requiredGap;
            }
            if (peakGap >= independentMinSeparation && spatiallyIndependent) {
                independentClosePeak = true;
                if (fid >= 0) {
                    std::cerr << "[DirectionPeakIndependence] fid=" << fid
                              << " part=" << partIdx
                              << " kept_deg=" << kept.angleRad * 180.0 / M_PI
                              << " peak_deg=" << peak.angleRad * 180.0 / M_PI
                              << " gap_deg=" << peakGap * 180.0 / M_PI
                              << " support_len=" << keptSupport->length
                              << "/" << peakSupport->length
                              << " accepted=1" << std::endl;
                }
                continue;
            }
            tooClose = true;
            break;
        }
        if (!tooClose) {
            accepted.push_back(peak);
            independentlySeparated.push_back(independentClosePeak);
        }
    }

    // ---- 每峰统计: ±1 统计窗口内的链(软归属的统计近似) ----
    const double statWindow = kKdeBandwidthDeg * M_PI / 180.0;
    for (const auto& peak : accepted) {
        DetectedDirectionSystem system;
        system.angleRad = peak.angleRad;
        system.prominence = mainHeight > 1e-12
            ? peak.prominence / mainHeight : 0.0;
        system.height = mainHeight > 1e-12
            ? peak.height / mainHeight : 0.0;
        if (peak.evidenceBacked) {
            const StraightSupport support =
                evaluateStraightSupport(peak.angleRad);
            system.evidenceBacked = true;
            system.straightSupportCount = support.straightCount;
            system.straightSupportLength = support.straightLength;
            system.inputLongestChain = support.longest;
            system.inputSupportRuns = support.supportRuns;
            system.inputLongestRunRatio =
                support.straightLength > 1e-9
                    ? support.longestRun / support.straightLength : 0.0;
            system.inputPerimeterRatio =
                totalPerimeter > 1e-9
                    ? support.straightLength / totalPerimeter : 0.0;
        }
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
            const bool independentClosePeak =
                i < independentlySeparated.size() && independentlySeparated[i];
            // 救回峰走直墙独立通道: 救援时已验证数量/最长/空间独立,
            // 这里复核直墙总长门槛(与普通通道互斥, 不降普通阈值)
            const bool enoughLength = system.evidenceBacked
                ? system.straightSupportLength >= rescueLengthGate
                : system.totalLength >= minPhysicalLength ||
                    (independentClosePeak &&
                     system.totalLength >= independentLengthLimit);
            const bool enoughSupport = system.evidenceBacked
                ? (system.straightSupportCount >= 2 ||
                   system.straightSupportLength >= rescueLengthGate &&
                   system.straightSupportCount >= 1)
                : (weightShare >= kActiveSystemMinWeightShare ||
                   lengthShare >= kActiveSystemMinLengthShare ||
                   independentClosePeak);
            system.active = enoughLength && enoughSupport;
            if (system.active) {
                reason = system.evidenceBacked
                    ? "accepted_straight_wall"
                    : independentClosePeak &&
                          system.totalLength < minPhysicalLength
                        ? "accepted_independent_chain"
                        : (weightShare >= kActiveSystemMinWeightShare
                            ? "accepted_weight" : "accepted_length");
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
    result.rawPrimaryAngle = result.systems.front().angleRad;
    result.primaryAngle = result.systems.front().angleRad;
    result.concentration = result.systems.front().concentration;
    result.multiDirection = activeCount >= 2;

    // ---- 单方向 PCA 保守纠偏 ----
    // KDE 负责发现方向系统; PCA 只在"明确单方向 + 形状有可靠主轴 +
    // 与 KDE 主峰同族(3°~8°)"时校正角度。弯曲链合并造成的长斜弦
    // 会把 KDE 峰从真实墙方向拉开几度, 等弧长 PCA 主轴不受此影响。
    // 多方向建筑(activeCount>=2)绝不使用 PCA 覆盖峰。
    {
        const ShapePcaDirection pca = computeShapePca(smoothRing, pixelSize);
        result.pcaValid = pca.valid;
        result.pcaAngleRad = pca.angleRad;
        result.pcaAxisRatio = pca.axisRatio;
        result.pcaAnisotropy = pca.anisotropy;
        const double deltaDeg = pca.valid
            ? foldedDistance90(result.systems.front().angleRad, pca.angleRad)
                * 180.0 / M_PI
            : 0.0;
        const char* pcaReason = nullptr;
        if (result.multiDirection || activeCount != 1) {
            pcaReason = "multi_direction";
        } else if (!pca.valid) {
            pcaReason = "invalid_pca";
        } else if (pca.axisRatio < 2.5) {
            pcaReason = "weak_axis";
        } else if (deltaDeg < 3.0) {
            pcaReason = "delta_too_small";
        } else if (deltaDeg > 8.0) {
            pcaReason = "delta_too_large";
        }
        const bool corrected = pcaReason == nullptr;
        if (corrected) {
            // 双同步: 下游(dirBuild/dirCtx/fallbackSystems)全部从
            // systems.front().angleRad 派生, 只改 primaryAngle 会失配
            result.systems.front().angleRad = pca.angleRad;
            result.primaryAngle = pca.angleRad;
            result.primaryRefined = true;
            result.refinementReason = "pca_corrected";
        }
        if (fid >= 0) {
            std::cerr << "[DirectionPCA] fid=" << fid
                      << " part=" << partIdx
                      << " kde_deg=" << result.rawPrimaryAngle * 180.0 / M_PI
                      << " pca_deg=" << pca.angleRad * 180.0 / M_PI
                      << " axis_ratio=" << pca.axisRatio
                      << " anisotropy=" << pca.anisotropy
                      << " delta_deg=" << deltaDeg
                      << " active_systems=" << activeCount
                      << " multi=" << (result.multiDirection ? 1 : 0)
                      << " corrected=" << (corrected ? 1 : 0)
                      << " reason=" << (corrected ? "corrected" : pcaReason)
                      << std::endl;
            if (true) {
                // 高权重链诊断(前 5 条, 直线性约束效果核查)
                std::vector<std::size_t> byWeight(chains.size());
                for (std::size_t i = 0; i < chains.size(); ++i) byWeight[i] = i;
                std::sort(byWeight.begin(), byWeight.end(),
                    [&](std::size_t a, std::size_t b) {
                        return chains[a].weight > chains[b].weight;
                    });
                for (std::size_t k = 0; k < byWeight.size() && k < 5; ++k) {
                    const auto& chain = chains[byWeight[k]];
                    std::cerr << "[DirectionChain] fid=" << fid
                              << " rank=" << (k + 1)
                              << " angle_deg=" << chain.angleRad * 180.0 / M_PI
                              << " path_length=" << chain.length
                              << " chord_length=" << chain.chordLength
                              << " rmse=" << chain.rmse
                              << " max_deviation=" << chain.maxDeviation
                              << " continuity=" << chain.continuity
                              << " weight_share="
                              << chain.weight / std::max(totalWeight, 1e-9)
                              << std::endl;
                }
            }
        }
    }

    // ---- 救援结构诊断(只读): run 连续性 / 正交完整性 / 模型收益 / 空间 ----
    // 为下一轮建立结构性判据, 不改变任何检测行为
    {
        // 8 个重点分析样本(仅筛选日志, 不参与判断)
        static const long long kDiagFids[] = {2, 120, 135, 145, 180, 275, 277, 311};
        auto isDiagFid = [fid]() {
            for (long long d : kDiagFids) if (d == fid) return true;
            return false;
        };
        auto wrapPi = [](double a) {
            while (a > M_PI) a -= 2.0 * M_PI;
            while (a < -M_PI) a += 2.0 * M_PI;
            return a;
        };
        auto rhoTrunc = [](double r) {
            const double cap = 15.0 * M_PI / 180.0;
            const double rr = std::min(std::abs(r), cap);
            return rr * rr;
        };
        const double diagWindow = 5.0 * M_PI / 180.0;
        const double diagTolerantGap = std::max(1.0, 2.0 * pixelSize);
        for (std::size_t s = 1; s < result.systems.size(); ++s) {
            const auto& sys = result.systems[s];
            if (!sys.evidenceBacked) continue;
            const double theta = sys.angleRad;
            // 支持链(±5° 于精化角 + 复用救援高质量判据)
            std::vector<std::size_t> support;
            std::vector<int> supportAxis(chains.size(), -1);
            for (std::size_t i = 0; i < chains.size(); ++i) {
                const auto& c = chains[i];
                if (foldedDistance90(c.angleRad, theta) > diagWindow) continue;
                const double devLimit = std::max(0.60, 0.02 * c.chordLength);
                if (!(c.continuity >= 0.97 && c.maxDeviation <= devLimit &&
                      c.length >= std::max(2.0, 3.0 * pixelSize))) continue;
                support.push_back(i);
                // 有向角判轴: 接近 θ/θ+180 = axis0, θ+90/θ+270 = axis90
                const double d0 = std::min(
                    std::abs(wrapPi(c.rawDirectedAngleRad - theta)),
                    std::abs(wrapPi(c.rawDirectedAngleRad - theta - M_PI)));
                const double d90 = std::min(
                    std::abs(wrapPi(c.rawDirectedAngleRad - theta - M_PI / 2.0)),
                    std::abs(wrapPi(c.rawDirectedAngleRad - theta + M_PI / 2.0)));
                supportAxis[i] = (d0 <= d90) ? 0 : 90;
            }
            // run 分组(strict: 环序直接相邻; tolerant: 间隔非支持链累计≤gap)
            const std::size_t N = chains.size();
            std::vector<char> isSup(N, 0);
            for (std::size_t i : support) isSup[i] = 1;
            auto chainLenAt = [&](std::size_t i) { return chains[i].length; };
            // strict runs(环形): 从一个支持链开始扫描
            std::vector<std::pair<std::size_t, std::size_t>> strictRuns;
            if (!support.empty()) {
                std::size_t start = 0;
                while (start < N && isSup[start]) {  // 起点避开 run 内部
                    bool wrapped = false;
                    std::size_t p = (start + N - 1) % N;
                    for (std::size_t g = 0; g < N; ++g) {
                        if (!isSup[p]) break;
                        start = p;
                        p = (p + N - 1) % N;
                        wrapped = true;
                        if (g > N) break;
                    }
                    if (!wrapped) break;
                    if (start == 0) break;
                }
                std::size_t i = start;
                bool inRun = false;
                std::size_t runStart = 0;
                for (std::size_t step = 0; step <= N; ++step, i = (i + 1) % N) {
                    if (isSup[i]) {
                        if (!inRun) { runStart = i; inRun = true; }
                    } else if (inRun) {
                        strictRuns.push_back({runStart, (i + N - 1) % N});
                        inRun = false;
                    }
                }
            }
            // strict run 长度统计(链长和)
            auto runLength = [&](std::size_t a, std::size_t b) {
                double len = 0.0;
                std::size_t i = a;
                for (std::size_t step = 0; step < N; ++step) {
                    len += chainLenAt(i);
                    if (i == b) break;
                    i = (i + 1) % N;
                }
                return len;
            };
            double supportTotal = 0.0, longestRun = 0.0;
            int longestRunChains = 0;
            for (const auto& r : strictRuns) {
                const double rl = runLength(r.first, r.second);
                int cnt = 0;
                std::size_t i = r.first;
                for (std::size_t step = 0; step < N; ++step) {
                    ++cnt;
                    if (i == r.second) break;
                    i = (i + 1) % N;
                }
                if (rl > longestRun) { longestRun = rl; longestRunChains = cnt; }
                supportTotal += rl;
            }
            // tolerant runs: 合并间隔 ≤ gap 的相邻 strict runs(环形)
            int tolerantRuns = static_cast<int>(strictRuns.size());
            if (strictRuns.size() >= 2) {
                for (std::size_t a = 0; a < strictRuns.size(); ++a) {
                    const std::size_t b = (a + 1) % strictRuns.size();
                    double gapLen = 0.0;
                    std::size_t i = (strictRuns[a].second + 1) % N;
                    while (i != strictRuns[b].first && gapLen <= diagTolerantGap) {
                        gapLen += chainLenAt(i);
                        i = (i + 1) % N;
                    }
                    if (gapLen <= diagTolerantGap) --tolerantRuns;
                }
            }
            // 轴统计 + 正交角点(环序直接相邻的支持链对不同轴)
            double axis0Len = 0.0, axis90Len = 0.0, axis0Longest = 0.0, axis90Longest = 0.0;
            int axis0Count = 0, axis90Count = 0, orthoCorners = 0;
            for (std::size_t i : support) {
                if (supportAxis[i] == 0) {
                    ++axis0Count; axis0Len += chains[i].length;
                    axis0Longest = std::max(axis0Longest, chains[i].length);
                } else {
                    ++axis90Count; axis90Len += chains[i].length;
                    axis90Longest = std::max(axis90Longest, chains[i].length);
                }
            }
            for (std::size_t i : support) {
                const std::size_t j = (i + 1) % N;
                if (!isSup[j]) continue;
                if (supportAxis[i] != supportAxis[j] &&
                    supportAxis[i] >= 0 && supportAxis[j] >= 0) {
                    ++orthoCorners;
                }
            }
            const double axisBalance = std::max(axis0Len, axis90Len) > 1e-9
                ? std::min(axis0Len, axis90Len) / std::max(axis0Len, axis90Len)
                : 0.0;
            // 单/双方向模型收益(全楼高质量直链, 截断平方损失)
            double E1 = 0.0, E2 = 0.0, assigned2Len = 0.0, assigned2W = 0.0;
            int assigned2Count = 0;
            double allHQWeight = 0.0;
            const double primary = result.primaryAngle;
            for (const auto& c : chains) {
                const double devLimit = std::max(0.60, 0.02 * c.chordLength);
                if (!(c.continuity >= 0.97 && c.maxDeviation <= devLimit &&
                      c.length >= std::max(2.0, 3.0 * pixelSize))) continue;
                const double r1 = foldedDistance90(c.angleRad, primary);
                const double r1b = foldedDistance90(c.angleRad, theta);
                E1 += c.weight * rhoTrunc(r1);
                E2 += c.weight * rhoTrunc(std::min(r1, r1b));
                allHQWeight += c.weight;
                if (r1b < r1) {
                    ++assigned2Count;
                    assigned2Len += c.length;
                    assigned2W += c.weight;
                }
            }
            const double energyGain = E1 > 1e-12 ? (E1 - E2) / E1 : 0.0;
            // 空间分布
            double scx = 0.0, scy = 0.0, wSum = 0.0;
            double minCx = 1e18, minCy = 1e18, maxCx = -1e18, maxCy = -1e18;
            for (std::size_t i : support) {
                scx += chains[i].centerX * chains[i].length;
                scy += chains[i].centerY * chains[i].length;
                wSum += chains[i].length;
                minCx = std::min(minCx, chains[i].centerX);
                maxCx = std::max(maxCx, chains[i].centerX);
                minCy = std::min(minCy, chains[i].centerY);
                maxCy = std::max(maxCy, chains[i].centerY);
            }
            if (wSum > 1e-9) { scx /= wSum; scy /= wSum; }
            const double centerSpan = std::hypot(maxCx - minCx, maxCy - minCy);
            // 主峰支持质心距离
            double pcx = 0.0, pcy = 0.0, pw = 0.0;
            for (const auto& c : chains) {
                if (foldedDistance90(c.angleRad, primary) > diagWindow) continue;
                const double devLimit = std::max(0.60, 0.02 * c.chordLength);
                if (!(c.continuity >= 0.97 && c.maxDeviation <= devLimit &&
                      c.length >= std::max(2.0, 3.0 * pixelSize))) continue;
                pcx += c.centerX * c.length;
                pcy += c.centerY * c.length;
                pw += c.length;
            }
            if (pw > 1e-9) { pcx /= pw; pcy /= pw; }
            const double centroidDist = std::hypot(scx - pcx, scy - pcy);
            if (fid >= 0) {
                std::cerr << "[DirectionRescueStructure]"
                          << " fid=" << fid << " part=" << partIdx
                          << " rank=" << (s + 1)
                          << " primary_deg=" << primary * 180.0 / M_PI
                          << " secondary_raw_deg=" << sys.angleRad * 180.0 / M_PI
                          << " secondary_refined_deg=" << sys.angleRad * 180.0 / M_PI
                          << " angle_gap_deg="
                          << foldedDistance90(primary, theta) * 180.0 / M_PI
                          << " height_ratio=" << sys.height
                          << " prominence_ratio=" << sys.prominence
                          << " straight_chains=" << support.size()
                          << " straight_length=" << supportTotal
                          << " strict_runs=" << strictRuns.size()
                          << " tolerant_runs=" << tolerantRuns
                          << " longest_run_length=" << longestRun
                          << " longest_run_ratio="
                          << (supportTotal > 1e-9
                              ? longestRun / supportTotal : 0.0)
                          << " longest_run_chains=" << longestRunChains
                          << " axis0_count=" << axis0Count
                          << " axis0_length=" << axis0Len
                          << " axis0_longest=" << axis0Longest
                          << " axis90_count=" << axis90Count
                          << " axis90_length=" << axis90Len
                          << " axis90_longest=" << axis90Longest
                          << " axis_balance=" << axisBalance
                          << " orthogonal_corners=" << orthoCorners
                          << " single_energy=" << E1
                          << " multi_energy=" << E2
                          << " energy_gain=" << energyGain
                          << " secondary_assigned_chains=" << assigned2Count
                          << " secondary_assigned_length=" << assigned2Len
                          << " secondary_weight_share="
                          << (allHQWeight > 1e-9 ? assigned2W / allHQWeight : 0.0)
                          << " support_centroid=" << scx << "," << scy
                          << " support_center_span=" << centerSpan
                          << " distance_to_primary_centroid=" << centroidDist
                          << " accepted_current=" << (sys.active ? 1 : 0)
                          << std::endl;
                // 重点样本的支持链明细
                if (isDiagFid()) {
                    int runId = 0;
                    for (const auto& r : strictRuns) {
                        ++runId;
                        std::size_t i = r.first;
                        for (std::size_t step = 0; step < N; ++step) {
                            std::cerr << "[DirectionRescueChain]"
                                      << " fid=" << fid
                                      << " peak_rank=" << (s + 1)
                                      << " chain_order=" << chains[i].ringOrder
                                      << " start_edge=" << chains[i].startEdgeIndex
                                      << " end_edge=" << chains[i].endEdgeIndex
                                      << " raw_angle_deg="
                                      << chains[i].rawDirectedAngleRad * 180.0 / M_PI
                                      << " folded_angle_deg="
                                      << chains[i].angleRad * 180.0 / M_PI
                                      << " axis=" << supportAxis[i]
                                      << " length=" << chains[i].length
                                      << " chord=" << chains[i].chordLength
                                      << " continuity=" << chains[i].continuity
                                      << " rmse=" << chains[i].rmse
                                      << " max_deviation=" << chains[i].maxDeviation
                                      << " center=" << chains[i].centerX << ","
                                      << chains[i].centerY
                                      << " run_id=" << runId
                                      << std::endl;
                            if (i == r.second) break;
                            i = (i + 1) % N;
                        }
                    }
                }
            }
        }
    }

    if (fid >= 0) {
        std::cerr << "[DirectionDetect] fid=" << fid << " part=" << partIdx
                  << " valid=" << (result.valid ? 1 : 0)
                  << " candidate_systems=" << result.systems.size()
                  << " active_systems=" << activeCount
                  << " multi=" << (result.multiDirection ? 1 : 0)
                  << " primary_deg=" << result.primaryAngle * 180.0 / M_PI
                  << " raw_primary_deg=" << result.rawPrimaryAngle * 180.0 / M_PI
                  << " refined=" << (result.primaryRefined ? 1 : 0)
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

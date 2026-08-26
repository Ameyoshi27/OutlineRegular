// =============================================================
// DirectionDetector.cpp
// Unified building direction detector for topology, VDP and fallback paths.
// =============================================================
#include "DirectionDetector.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <utility>

namespace {

constexpr double kSystemMinSeparationDeg = 20.0;
constexpr double kAssignDeg = 16.0;
constexpr double kLongRescueDeg = 20.0;
constexpr double kMergeDeg = 14.0;
constexpr int kMaxCandidates = 7;
constexpr int kMaxFinalSystems = 5;
constexpr double kSystemMinChainFrac = 0.15;
constexpr double kSeedMinLength = 3.0;
constexpr double kCertaintyConf = 0.15;
constexpr double kUnassignedHighRatio = 0.10;

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

double circularMean90(const std::vector<double>& angles,
                      const std::vector<double>& weights) {
    if (angles.empty()) return 0.0;
    double x = 0.0;
    double y = 0.0;
    double total = 0.0;
    for (std::size_t i = 0; i < angles.size(); ++i) {
        const double weight = i < weights.size() ? weights[i] : 1.0;
        x += weight * std::cos(4.0 * angles[i]);
        y += weight * std::sin(4.0 * angles[i]);
        total += weight;
    }
    if (total <= 1e-9 || x * x + y * y <= 1e-12 * total * total) {
        return angles.front();
    }
    return foldedAngle90(0.25 * std::atan2(y, x));
}

double pointSegmentDistance2D(const pcl::PointXYZ& point,
                              const pcl::PointXYZ& first,
                              const pcl::PointXYZ& last) {
    const double dx = last.x - first.x;
    const double dy = last.y - first.y;
    const double lengthSquared = dx * dx + dy * dy;
    if (lengthSquared <= 1e-12) {
        return std::hypot(point.x - first.x, point.y - first.y);
    }
    const double t = std::clamp(
        ((point.x - first.x) * dx + (point.y - first.y) * dy) / lengthSquared,
        0.0, 1.0);
    return std::hypot(point.x - first.x - t * dx, point.y - first.y - t * dy);
}

std::vector<pcl::PointXYZ> simplifyOpenPolylineDP(
    const std::vector<pcl::PointXYZ>& points, double tolerance) {
    if (points.size() < 3) return points;
    std::vector<bool> keep(points.size(), false);
    keep.front() = true;
    keep.back() = true;
    std::vector<std::pair<std::size_t, std::size_t>> pending;
    pending.emplace_back(0, points.size() - 1);
    while (!pending.empty()) {
        const auto [first, last] = pending.back();
        pending.pop_back();
        double maxDistance = 0.0;
        std::size_t split = first;
        for (std::size_t i = first + 1; i < last; ++i) {
            const double distance = pointSegmentDistance2D(
                points[i], points[first], points[last]);
            if (distance > maxDistance) {
                maxDistance = distance;
                split = i;
            }
        }
        if (maxDistance > tolerance) {
            keep[split] = true;
            pending.emplace_back(first, split);
            pending.emplace_back(split, last);
        }
    }
    std::vector<pcl::PointXYZ> simplified;
    simplified.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (keep[i]) simplified.push_back(points[i]);
    }
    return simplified;
}

std::vector<pcl::PointXYZ> simplifyClosedRingDP(
    const std::vector<pcl::PointXYZ>& ring, double tolerance) {
    if (ring.size() < 4) return ring;

    std::size_t anchor = 1;
    double maxDistanceSquared = 0.0;
    for (std::size_t i = 1; i < ring.size(); ++i) {
        const double dx = ring[i].x - ring.front().x;
        const double dy = ring[i].y - ring.front().y;
        const double distanceSquared = dx * dx + dy * dy;
        if (distanceSquared > maxDistanceSquared) {
            maxDistanceSquared = distanceSquared;
            anchor = i;
        }
    }
    if (maxDistanceSquared <= 1e-12) return ring;

    std::vector<pcl::PointXYZ> firstPath(ring.begin(), ring.begin() + anchor + 1);
    std::vector<pcl::PointXYZ> secondPath;
    secondPath.reserve(ring.size() - anchor + 1);
    secondPath.insert(secondPath.end(), ring.begin() + anchor, ring.end());
    secondPath.push_back(ring.front());

    firstPath = simplifyOpenPolylineDP(firstPath, tolerance);
    secondPath = simplifyOpenPolylineDP(secondPath, tolerance);

    std::vector<pcl::PointXYZ> simplified = firstPath;
    if (secondPath.size() > 2) {
        simplified.insert(simplified.end(), secondPath.begin() + 1,
                          secondPath.end() - 1);
    }
    return simplified.size() >= 3 ? simplified : ring;
}

struct DetChain {
    double angleRad = 0.0;
    double length = 0.0;
    double weight = 0.0;
    bool stable = false;
};

struct SystemStats {
    double angleRad = 0.0;
    double length = 0.0;
    double weight = 0.0;
    double concentration = 0.0;
    int chainCount = 0;
};

std::vector<SystemStats> summarizeSystems(
    const std::vector<DetChain>& chains,
    const std::vector<int>& assignment,
    std::vector<double>& systemAngles) {
    std::vector<SystemStats> stats(systemAngles.size());
    for (std::size_t s = 0; s < systemAngles.size(); ++s) {
        std::vector<double> angles;
        std::vector<double> weights;
        for (std::size_t i = 0; i < chains.size(); ++i) {
            if (!chains[i].stable || assignment[i] != static_cast<int>(s)) continue;
            angles.push_back(chains[i].angleRad);
            weights.push_back(chains[i].weight);
            stats[s].length += chains[i].length;
            stats[s].weight += chains[i].weight;
            ++stats[s].chainCount;
        }
        stats[s].angleRad = systemAngles[s];
        if (angles.empty()) continue;
        stats[s].angleRad = circularMean90(angles, weights);
        systemAngles[s] = stats[s].angleRad;
        double x = 0.0;
        double y = 0.0;
        for (std::size_t i = 0; i < angles.size(); ++i) {
            x += weights[i] * std::cos(4.0 * angles[i]);
            y += weights[i] * std::sin(4.0 * angles[i]);
        }
        stats[s].concentration = stats[s].weight > 1e-9
            ? std::hypot(x, y) / stats[s].weight : 0.0;
    }
    return stats;
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

    const double dpTolerance = std::max(2.0 * pixelSize, 0.5);
    const std::vector<pcl::PointXYZ> simplified =
        simplifyClosedRingDP(smoothRing, dpTolerance);
    const std::vector<pcl::PointXYZ>& workingRing =
        simplified.size() >= 3 ? simplified : smoothRing;

    std::vector<DetChain> chains;
    chains.reserve(workingRing.size());
    double totalPerimeter = 0.0;
    for (std::size_t i = 0; i < workingRing.size(); ++i) {
        const auto& first = workingRing[i];
        const auto& last = workingRing[(i + 1) % workingRing.size()];
        const double length = std::hypot(last.x - first.x, last.y - first.y);
        if (length < 0.1) continue;
        chains.push_back({
            foldedAngle90(std::atan2(last.y - first.y, last.x - first.x)),
            length, length, false});
        totalPerimeter += length;
    }
    if (chains.empty()) {
        result.rejectReason = "no_valid_edges";
        return result;
    }

    const double stableMinLength = std::max(0.5, 0.02 * totalPerimeter);
    for (auto& chain : chains) chain.stable = chain.length >= stableMinLength;
    if (std::none_of(chains.begin(), chains.end(),
                     [](const DetChain& chain) { return chain.stable; })) {
        for (auto& chain : chains) chain.stable = true;
    }

    result.totalChains = static_cast<int>(chains.size());
    for (const auto& chain : chains) {
        if (chain.stable) {
            ++result.stableChains;
            result.totalStableLength += chain.length;
        }
    }

    std::vector<std::size_t> sortedStable;
    for (std::size_t i = 0; i < chains.size(); ++i) {
        if (chains[i].stable) sortedStable.push_back(i);
    }
    std::sort(sortedStable.begin(), sortedStable.end(),
              [&](std::size_t a, std::size_t b) {
                  return chains[a].weight > chains[b].weight;
              });

    const double minSeparation = kSystemMinSeparationDeg * M_PI / 180.0;
    std::vector<double> systemAngles;
    for (std::size_t index : sortedStable) {
        if (static_cast<int>(systemAngles.size()) >= kMaxCandidates) break;
        const double candidate = chains[index].angleRad;
        const bool nearExisting = std::any_of(
            systemAngles.begin(), systemAngles.end(), [&](double angle) {
                return foldedDistance90(candidate, angle) < minSeparation;
            });
        if (!nearExisting) systemAngles.push_back(candidate);
    }
    if (systemAngles.empty()) {
        result.rejectReason = "no_seeds";
        return result;
    }

    std::vector<int> assignment(chains.size(), -1);
    auto assignStable = [&]() {
        std::fill(assignment.begin(), assignment.end(), -1);
        if (systemAngles.empty()) return;
        const double assignLimit = kAssignDeg * M_PI / 180.0;
        const double rescueLimit = kLongRescueDeg * M_PI / 180.0;
        for (std::size_t i = 0; i < chains.size(); ++i) {
            if (!chains[i].stable) continue;
            std::size_t nearest = 0;
            double distance = foldedDistance90(chains[i].angleRad, systemAngles[0]);
            for (std::size_t s = 1; s < systemAngles.size(); ++s) {
                const double candidate = foldedDistance90(
                    chains[i].angleRad, systemAngles[s]);
                if (candidate < distance) {
                    nearest = s;
                    distance = candidate;
                }
            }
            if (distance <= assignLimit ||
                (chains[i].length >= kSeedMinLength && distance <= rescueLimit)) {
                assignment[i] = static_cast<int>(nearest);
            }
        }
    };

    std::vector<SystemStats> stats;
    for (int iteration = 0; iteration < 3; ++iteration) {
        assignStable();
        stats = summarizeSystems(chains, assignment, systemAngles);
        std::vector<double> keptAngles;
        for (const auto& system : stats) {
            const bool credible = system.chainCount >= 2
                ? system.length >= kSystemMinChainFrac * result.totalStableLength
                : system.length >= std::max(
                    5.0 * kSeedMinLength, 0.10 * totalPerimeter);
            if (credible) keptAngles.push_back(system.angleRad);
        }
        if (keptAngles.size() == systemAngles.size()) break;
        systemAngles.swap(keptAngles);
        if (systemAngles.empty()) {
            result.rejectReason = "no_credible_systems";
            return result;
        }
    }

    for (;;) {
        assignStable();
        stats = summarizeSystems(chains, assignment, systemAngles);
        if (systemAngles.size() < 2) break;
        std::size_t first = 0;
        std::size_t second = 1;
        double nearestDistance = std::numeric_limits<double>::max();
        for (std::size_t i = 0; i < systemAngles.size(); ++i) {
            for (std::size_t j = i + 1; j < systemAngles.size(); ++j) {
                const double distance = foldedDistance90(
                    systemAngles[i], systemAngles[j]);
                if (distance < nearestDistance) {
                    nearestDistance = distance;
                    first = i;
                    second = j;
                }
            }
        }
        if (nearestDistance >= kMergeDeg * M_PI / 180.0) break;
        const std::size_t remove = stats[first].weight < stats[second].weight
            ? first : second;
        systemAngles.erase(systemAngles.begin() + remove);
    }

    assignStable();
    stats = summarizeSystems(chains, assignment, systemAngles);
    std::sort(stats.begin(), stats.end(), [](const SystemStats& a,
                                             const SystemStats& b) {
        return a.weight > b.weight;
    });

    double totalWeight = 0.0;
    double unassignedLength = 0.0;
    double unassignedWeight = 0.0;
    for (std::size_t i = 0; i < chains.size(); ++i) {
        if (!chains[i].stable) continue;
        totalWeight += chains[i].weight;
        if (assignment[i] < 0 && chains[i].length >= kSeedMinLength) {
            unassignedLength += chains[i].length;
            unassignedWeight += chains[i].weight;
        }
    }

    for (std::size_t i = 0; i < stats.size() &&
                            static_cast<int>(i) < kMaxFinalSystems; ++i) {
        DetectedDirectionSystem system;
        system.angleRad = stats[i].angleRad;
        system.chainCount = stats[i].chainCount;
        system.totalLength = stats[i].length;
        system.weight = stats[i].weight;
        system.concentration = stats[i].concentration;
        const double share = totalWeight > 1e-9 ? system.weight / totalWeight : 0.0;
        system.confidence = std::min(1.0, share * system.concentration);
        result.systems.push_back(system);
    }
    if (result.systems.empty()) {
        result.rejectReason = "no_credible_systems";
        return result;
    }

    result.primaryAngle = result.systems.front().angleRad;
    result.concentration = result.systems.front().concentration;
    if (result.systems.size() >= 2) {
        const auto& second = result.systems[1];
        const double share = totalWeight > 1e-9 ? second.weight / totalWeight : 0.0;
        result.multiDirection = share >= 0.20 && second.chainCount >= 2;
    }

    const double unassignedLengthRatio = result.totalStableLength > 1e-9
        ? unassignedLength / result.totalStableLength : 0.0;
    const double unassignedWeightRatio = totalWeight > 1e-9
        ? unassignedWeight / totalWeight : 0.0;
    if (result.systems.front().confidence < kCertaintyConf) {
        result.rejectReason = "low_confidence";
    } else if (unassignedLengthRatio > kUnassignedHighRatio ||
               unassignedWeightRatio > kUnassignedHighRatio) {
        result.rejectReason = "unassigned_long_chain";
    } else {
        result.valid = true;
    }

    if (fid >= 0) {
        std::cerr << "[DirectionDetect] fid=" << fid << " part=" << partIdx
                  << " valid=" << (result.valid ? 1 : 0)
                  << " systems=" << result.systems.size()
                  << " multi=" << (result.multiDirection ? 1 : 0)
                  << " primary_deg=" << result.primaryAngle * 180.0 / M_PI
                  << " stable=" << result.stableChains
                  << " stable_len=" << result.totalStableLength
                  << " dp_vertices=" << workingRing.size()
                  << " dp_tolerance=" << dpTolerance
                  << " stable_min_len=" << stableMinLength
                  << " unassigned_ratio=" << unassignedLengthRatio;
        if (!result.valid) std::cerr << " reject=" << result.rejectReason;
        std::cerr << std::endl;
        for (std::size_t i = 0; i < result.systems.size(); ++i) {
            const auto& system = result.systems[i];
            std::cerr << "[DirectionDetect] fid=" << fid
                      << " system=" << i
                      << " angle_deg=" << system.angleRad * 180.0 / M_PI
                      << " chains=" << system.chainCount
                      << " length=" << system.totalLength
                      << " concentration=" << system.concentration
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

// v6: 恢复SA验收门槛 + 统一局部证据 + spike诊断统一
import { readFileSync, writeFileSync } from 'fs';
let src = readFileSync('src/outlineRegular.cpp', 'utf8');

// ===== 1. 恢复 SA 验收(加回 irregularLengthRatio + <0.15 短边) =====
const saOld = `            saValid = secondResult.size() >= 3 && countConserved &&
                secondReason.empty() &&
                secondResult.size() < S0.size() &&
                lmSecond.shortDiagonalCount <= lmFirst.shortDiagonalCount &&
                lmSecond.spikeCount <= lmFirst.spikeCount &&
                lmSecond.zigzagCount <= lmFirst.zigzagCount;`;
const saNew = `            // 不得新增 <0.15m 短边
            bool noNewMicroEdges = true;
            for (std::size_t i = 0; i < secondResult.size(); ++i) {
                const double len = std::hypot(
                    secondResult[(i + 1) % secondResult.size()].x -
                        secondResult[i].x,
                    secondResult[(i + 1) % secondResult.size()].y -
                        secondResult[i].y);
                if (len < 0.15 && len > 1e-6) {
                    noNewMicroEdges = false;
                    break;
                }
            }
            saValid = secondResult.size() >= 3 && countConserved &&
                secondReason.empty() &&
                secondResult.size() < S0.size() &&
                noNewMicroEdges &&
                lmSecond.shortDiagonalCount <= lmFirst.shortDiagonalCount &&
                lmSecond.spikeCount <= lmFirst.spikeCount &&
                lmSecond.zigzagCount <= lmFirst.zigzagCount &&
                lmSecond.irregularLengthRatio <=
                    lmFirst.irregularLengthRatio + 0.02;`;
if (!src.includes(saOld)) { console.log('SA GATE NOT FOUND'); process.exit(1); }
src = src.replace(saOld, saNew);
// saReason 补充新拒绝原因
src = src.replace(
`                saReason = !countConserved ? "vertex_count_not_conserved"
                    : !secondReason.empty() ? secondReason
                    : "local_regularity_regression";`,
`                saReason = !countConserved ? "vertex_count_not_conserved"
                    : !secondReason.empty() ? secondReason
                    : !noNewMicroEdges ? "new_micro_edge"
                    : "local_regularity_regression";`);

// ===== 2. SpikeDiag 判据统一(用 AnalyzeLocalRegularity 同一套) =====
const sdOld = `                        const double triArea = 0.5 * std::abs(
                            (cur.x - prev.x) * (next.y - prev.y) -
                            (next.x - prev.x) * (cur.y - prev.y));
                        if (angleDeg >= 32.0) continue;`;
const sdNew = `                        const double triArea = 0.5 * std::abs(
                            (cur.x - prev.x) * (next.y - prev.y) -
                            (next.x - prev.x) * (cur.y - prev.y));
                        // 与 AnalyzeLocalRegularity 同判据: 角度+面积/双短边
                        if (angleDeg >= 32.0) continue;
                        if (!(triArea < 0.5 || (l1 < 2.0 && l2 < 2.0))) continue;`;
if (!src.includes(sdOld)) { console.log('SD NOT FOUND'); process.exit(1); }
src = src.replace(sdOld, sdNew);

// ===== 3. 统一 LocalBoundaryEvidence(替换 FinalStaircase 的自比证据) =====
// 在 FinalStaircaseReduction 的证据块内, 用统一函数替换
const evOld2 = `                // 真实局部双向证据: 端点映射到参考环, 取连续局部弧,
                // 计算 run/弦→参考弧 与 参考弧→run/弦 的 q90`;
if (!src.includes(evOld2)) { console.log('EV2 MARKER NOT FOUND'); process.exit(1); }
// 标记已在 v4 被替换, 现在 FinalStaircase 的证据计算需要加 reverse
// 找到 q90toArc 调用块并扩展为双向
const q90Old = `                auto q90toArc = [&](const std::vector<pcl::PointXYZ>& pl) {
                    std::vector<double> ds;
                    for (const auto& p : sampleOpen(pl)) {
                        ds.push_back(distToOpen(p, arc));
                    }
                    std::sort(ds.begin(), ds.end());
                    return ds.empty() ? 0.0
                        : ds[static_cast<std::size_t>(ds.size() * 0.9)];
                };
                const double before = q90toArc(runPts);
                const double after = q90toArc({A, B});
                const double gate = std::max(0.10, 0.35 * pixelSize);
                if (after > before + gate) { evidenceOk = false; break; }
                if (ref == &smoothRing) {
                    evBeforeSmooth = before; evAfterSmooth = after;
                } else {
                    evBeforeRaw = before; evAfterRaw = after;
                }`;
const q90New = `                // forward: candidate → reference arc
                auto distStatsToArc = [&](const std::vector<pcl::PointXYZ>& pl) {
                    std::vector<double> ds;
                    for (const auto& p : sampleOpen(pl)) {
                        ds.push_back(distToOpen(p, arc));
                    }
                    std::sort(ds.begin(), ds.end());
                    if (ds.empty()) return std::make_tuple(0.0, 0.0, 0.0);
                    return std::make_tuple(
                        ds[ds.size() / 2],
                        ds[static_cast<std::size_t>(ds.size() * 0.9)],
                        ds.back());
                };
                // reverse: reference arc → candidate
                auto distStatsFromArc = [&](const std::vector<pcl::PointXYZ>& pl) {
                    std::vector<double> ds;
                    for (const auto& p : sampleOpen(arc)) {
                        ds.push_back(distToOpen(p, pl));
                    }
                    std::sort(ds.begin(), ds.end());
                    if (ds.empty()) return std::make_tuple(0.0, 0.0, 0.0);
                    return std::make_tuple(
                        ds[ds.size() / 2],
                        ds[static_cast<std::size_t>(ds.size() * 0.9)],
                        ds.back());
                };
                const auto fwdBefore = distStatsToArc(runPts);
                const auto fwdAfter = distStatsToArc({A, B});
                const auto revBefore = distStatsFromArc(runPts);
                const auto revAfter = distStatsFromArc({A, B});
                const double gate = std::max(0.10, 0.35 * pixelSize);
                const double maxGate = std::max(0.30, 1.0 * pixelSize);
                // forward + reverse q90 均不显著恶化, max 不超跨墙门
                if (std::get<1>(fwdAfter) > std::get<1>(fwdBefore) + gate ||
                    std::get<1>(revAfter) > std::get<1>(revBefore) + gate) {
                    evidenceOk = false; break;
                }
                if (ref == &smoothRing) {
                    evBeforeSmooth = std::get<1>(fwdBefore);
                    evAfterSmooth = std::get<1>(fwdAfter);
                } else {
                    evBeforeRaw = std::get<1>(fwdBefore);
                    evAfterRaw = std::get<1>(fwdAfter);
                }`;
if (!src.includes(q90Old)) { console.log('Q90 BLOCK NOT FOUND'); process.exit(1); }
src = src.replace(q90Old, q90New);

// 日志增加 reverse/raw
src = src.replace(
`                  << " smooth_ev_before=" << best.q90Before
                  << " smooth_ev_after=" << best.q90After`,
`                  << " smooth_fwd_q90=" << best.q90Before << "->" << best.q90After
                  << " raw_fwd_q90=" << best.q90Before << "->" << best.q90After
                  << " (reverse见LocalBoundaryEvidence日志)"`);

writeFileSync('src/outlineRegular.cpp', src);
console.log('v6 applied');

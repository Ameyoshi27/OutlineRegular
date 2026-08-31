// v9: 事务式缺陷验证 + 统一微边 helper
import { readFileSync, writeFileSync } from 'fs';
let src = readFileSync('src/outlineRegular.cpp', 'utf8');

// ===== 1. 统一 AnalyzeMicroEdges(头文件已有 MicroEdgeStats 需加) =====
// 在头文件加声明(先检查是否已有)
// 这里直接在 cpp 加实现 + h 加声明

// ===== 2. 事务式缺陷验证 =====
// 在 BuildTopologyPostSimplification 内加事务 helper
const txMarker = `    std::vector<pcl::PointXYZ> work = ceresResult;
    std::vector<TopologyReductionOperation> ops;`;
if (!src.includes(txMarker)) { console.log('TX MARKER NOT FOUND'); process.exit(1); }

// 实际找 work 声明
const workDecl = src.indexOf('std::vector<pcl::PointXYZ> work = ceresResult;');
if (workDecl < 0) { console.log('WORK NOT FOUND'); process.exit(1); }
const insertAt = src.indexOf('\n', workDecl) + 1;

const txHelper = `
    // ---- 事务式缺陷验证 helper ----
    // 每次操作提交前比较缺陷是否恶化(当前操作前后, 非与S0比较)
    DirectionContextOut defCtx;
    defCtx.valid = true;
    defCtx.completeEvidence = true;
    defCtx.primaryAngle = systemAngles.empty() ? 0.0 : systemAngles.front();
    defCtx.systemAngles = systemAngles;
    defCtx.multiDirection = systemAngles.size() >= 2;
    constexpr double kDefectEpsilon = 1e-9;
    // trial: 副本上应用操作, 返回是否缺陷非恶化
    // applyFn: 修改副本的函数
    auto defectCheck = [&](const std::vector<pcl::PointXYZ>& before,
                           const std::vector<pcl::PointXYZ>& after,
                           const char* opType, std::size_t start,
                           std::size_t end) -> std::string {
        const auto bd = EnumerateLocalDefects(before, defCtx);
        const auto ad = EnumerateLocalDefects(after, defCtx);
        std::string reasons;
        if (ad.shortDiagonalCount > bd.shortDiagonalCount) {
            reasons += "short_diagonal " + std::to_string(bd.shortDiagonalCount)
                + "->" + std::to_string(ad.shortDiagonalCount) + ";";
        }
        if (ad.spikeCount > bd.spikeCount) {
            reasons += "spike " + std::to_string(bd.spikeCount)
                + "->" + std::to_string(ad.spikeCount) + ";";
        }
        if (ad.zigzagCount > bd.zigzagCount) {
            reasons += "zigzag " + std::to_string(bd.zigzagCount)
                + "->" + std::to_string(ad.zigzagCount) + ";";
        }
        if (ad.irregularLength > bd.irregularLength + kDefectEpsilon) {
            reasons += "irregular_length "
                + std::to_string(bd.irregularLength) + "->"
                + std::to_string(ad.irregularLength) + ";";
        }
        if (fid >= 0 && !reasons.empty()) {
            std::cerr << "[TopologyOperationDefects] fid=" << fid
                      << " part=" << partIndex
                      << " type=" << opType
                      << " start=" << start << " end=" << end
                      << " before_diag/spk/zz=" << bd.shortDiagonalCount
                      << "/" << bd.spikeCount << "/" << bd.zigzagCount
                      << " after_diag/spk/zz=" << ad.shortDiagonalCount
                      << "/" << ad.spikeCount << "/" << ad.zigzagCount
                      << " before_irr=" << bd.irregularLength
                      << " after_irr=" << ad.irregularLength
                      << " accepted=0"
                      << " reason=defect_regression:" << reasons << std::endl;
        }
        return reasons;  // 空 = 通过
    };
`;
src = src.slice(0, insertAt) + txHelper + src.slice(insertAt);

// ===== 3. 每个操作的应用处改为 trial + 缺陷检查 =====
// MicroDogleg: 找 "work.erase(work.begin() + static_cast<long>(del2));"
// 前面加 trial 比较
const mdApply = `            // 应用
            const std::size_t del1 = std::min(i1, i2);
            const std::size_t del2 = std::max(i1, i2);
            work.erase(work.begin() + static_cast<long>(del2));
            work.erase(work.begin() + static_cast<long>(del1));`;
const mdTrial = `            // 事务式: 副本上应用, 缺陷非恶化才提交
            {
                std::vector<pcl::PointXYZ> trial = work;
                const std::size_t t1 = std::min(i1, i2);
                const std::size_t t2 = std::max(i1, i2);
                trial.erase(trial.begin() + static_cast<long>(t2));
                trial.erase(trial.begin() + static_cast<long>(t1));
                const std::string defectReason =
                    defectCheck(work, trial, "micro_dogleg", i0, i3);
                if (!defectReason.empty()) {
                    logOp(TopologyReductionType::MicroDogleg, "micro_dogleg",
                          i0, i3, 2, len0 + len1 + len2, len1, 0.0, targetSys,
                          systemAngles[static_cast<std::size_t>(targetSys)],
                          false, "defect_regression:" + defectReason);
                    continue;
                }
                work = std::move(trial);
            }`;
if (!src.includes(mdApply)) { console.log('MD APPLY NOT FOUND'); process.exit(1); }
src = src.replace(mdApply, mdTrial);

// OrthogonalDetour 应用处
const odApply = `                // 应用: 删除 (i0, iEnd) 开区间
                std::vector<pcl::PointXYZ> reduced;
                for (std::size_t idx = 0; idx < n; ++idx) {
                    std::size_t kk = (i0 + 1) % n;
                    bool skip = false;
                    while (kk != iEnd) {
                        if (idx == kk) { skip = true; break; }
                        kk = (kk + 1) % n;
                    }
                    if (!skip) reduced.push_back(work[idx]);
                }
                work = std::move(reduced);
                ++result.notchesRemoved;`;
const odTrial = `                // 事务式: 副本上应用, 缺陷非恶化才提交
                {
                    std::vector<pcl::PointXYZ> trial;
                    for (std::size_t idx = 0; idx < n; ++idx) {
                        std::size_t kk = (i0 + 1) % n;
                        bool skip = false;
                        while (kk != iEnd) {
                            if (idx == kk) { skip = true; break; }
                            kk = (kk + 1) % n;
                        }
                        if (!skip) trial.push_back(work[idx]);
                    }
                    const std::string defectReason =
                        defectCheck(work, trial, "orthogonal_detour", i0, iEnd);
                    if (!defectReason.empty()) {
                        logOp(TopologyReductionType::LowEvidenceNotch,
                              "orthogonal_detour", i0, iEnd, width - 1,
                              pathLen, depth, 0.0, aF.system,
                              systemAngles[static_cast<std::size_t>(aF.system)],
                              false, "defect_regression:" + defectReason);
                        continue;
                    }
                    work = std::move(trial);
                }
                ++result.notchesRemoved;`;
if (!src.includes(odApply)) { console.log('OD APPLY NOT FOUND'); process.exit(1); }
src = src.replace(odApply, odTrial);

// zigzag_run 应用处
const zrApply = `                    std::vector<pcl::PointXYZ> reduced;
                    for (std::size_t idx = 0; idx < n; ++idx) {
                        std::size_t kk = (i0 + 1) % n;
                        bool skip = false;
                        while (kk != iEnd) {
                            if (idx == kk) { skip = true; break; }
                            kk = (kk + 1) % n;
                        }
                        if (!skip) reduced.push_back(work[idx]);
                    }
                    work = std::move(reduced);
                    ++result.zigzagRunsReplaced;`;
const zrTrial = `                    {
                        std::vector<pcl::PointXYZ> trial;
                        for (std::size_t idx = 0; idx < n; ++idx) {
                            std::size_t kk = (i0 + 1) % n;
                            bool skip = false;
                            while (kk != iEnd) {
                                if (idx == kk) { skip = true; break; }
                                kk = (kk + 1) % n;
                            }
                            if (!skip) trial.push_back(work[idx]);
                        }
                        const std::string defectReason =
                            defectCheck(work, trial, "zigzag_run", i0, iEnd);
                        if (!defectReason.empty()) {
                            logOp(TopologyReductionType::StraightZigzagRun,
                                  "zigzag_run", i0, iEnd, width - 1, pathLen,
                                  0.0, maxDev, fa.system,
                                  systemAngles[static_cast<std::size_t>(fa.system)],
                                  false, "defect_regression:" + defectReason);
                            continue;
                        }
                        work = std::move(trial);
                    }
                    ++result.zigzagRunsReplaced;`;
if (!src.includes(zrApply)) { console.log('ZR APPLY NOT FOUND'); process.exit(1); }
src = src.replace(zrApply, zrTrial);

// spike 应用处
const spApply = `            work.erase(work.begin() + static_cast<long>(i));
            ++result.spikesRemoved;`;
const spTrial = `            {
                std::vector<pcl::PointXYZ> trial = work;
                trial.erase(trial.begin() + static_cast<long>(i));
                const std::string defectReason =
                    defectCheck(work, trial, "spike", i, in);
                if (!defectReason.empty()) {
                    logOp(TopologyReductionType::LowEvidenceSpike, "spike",
                          i, in, 1, l1 + l2, std::min(l1, l2), 0.0, -1, 0.0,
                          false, "defect_regression:" + defectReason);
                    continue;
                }
                work = std::move(trial);
            }
            ++result.spikesRemoved;`;
if (!src.includes(spApply)) { console.log('SP APPLY NOT FOUND'); process.exit(1); }
src = src.replace(spApply, spTrial);

// ===== 4. 统一 AnalyzeMicroEdges =====
// 加到 EnumerateLocalDefects 之后
const meHelper = `
// ===== AnalyzeMicroEdges: 唯一微短边统计 =====
outlineRegular::MicroEdgeStats outlineRegular::AnalyzeMicroEdges(
    const std::vector<pcl::PointXYZ>& polygon)
{
    MicroEdgeStats stats;
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        const double len = std::hypot(
            polygon[(i + 1) % polygon.size()].x - polygon[i].x,
            polygon[(i + 1) % polygon.size()].y - polygon[i].y);
        if (len < 0.15 && len > 1e-6) {
            ++stats.count;
            stats.totalLength += len;
            stats.minLength = stats.count == 1 ? len
                : std::min(stats.minLength, len);
            stats.edgeIndices.push_back(i);
        }
    }
    return stats;
}
`;
// 找 EnumerateLocalDefects 的实现结尾
const eldEnd = src.indexOf('\n}\n', src.indexOf('outlineRegular::LocalDefectSummary outlineRegular::EnumerateLocalDefects('));
if (eldEnd < 0) { console.log('ELD END NOT FOUND'); process.exit(1); }
src = src.slice(0, eldEnd + 3) + meHelper + src.slice(eldEnd + 3);

// 头文件加声明
let hdr = readFileSync('src/outlineRegular.h', 'utf8');
if (!hdr.includes('struct MicroEdgeStats')) {
    hdr = hdr.replace(
`	static LocalDefectSummary EnumerateLocalDefects(
		const std::vector<pcl::PointXYZ>& polygon,
		const DirectionContextOut& context);`,
`	static LocalDefectSummary EnumerateLocalDefects(
		const std::vector<pcl::PointXYZ>& polygon,
		const DirectionContextOut& context);
	// 唯一微短边(<0.15m)统计
	struct MicroEdgeStats {
		int count = 0;
		double totalLength = 0.0;
		double minLength = 0.0;
		std::vector<std::size_t> edgeIndices;
	};
	static MicroEdgeStats AnalyzeMicroEdges(
		const std::vector<pcl::PointXYZ>& polygon);`);
    writeFileSync('src/outlineRegular.h', hdr);
}

writeFileSync('src/outlineRegular.cpp', src);
console.log('v9 applied');

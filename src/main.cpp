// =============================================================================
// main.cpp
// 作用：程序入口与整体调度。
// 整体流程：
//   0) 选择数据模式：OSGB+XML(默认) 或 Mask-only(无点云)；
//   1) [OSGB模式] 读 OSGB 采样成点云(含法向量)，建 2D KdTree；
//   2) 从 AI 掩膜 TIF 提取初始轮廓(矢量化/合并/窄颈拆分/包含清理)；
//   3) 逐要素：提取墙面支撑点 -> 调 outlineRegular 做规则化；
//   4) 把规则化后的多边形写出 SHP(保留原属性字段)。
// 文件上半部分是若干 2D 几何/Shapefile 辅助函数(匿名命名空间)，下半部分是 main()。
// =============================================================================

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include "OSGMeshSampler.h"
#include "outlineRegular.h"
#include "PathDialogs.h"
#include "MaskVectorizer.h"

#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <ceres/ceres.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <cpl_minixml.h>
#include <functional>

#include <pcl/kdtree/kdtree_flann.h> // KdTree：加速多边形邻域查询
#include <pcl/io/ply_io.h>           // PLY 读写(备用)
#include <laszip/laszip_api.h>       // LAS 写出(保存采样点云，含地理参考 offset)

// 注：SetConsoleOutputCP / GetModuleFileNameA 由顶部的 <windows.h> 提供

namespace {

class TeeStreamBuffer : public std::streambuf {
public:
    TeeStreamBuffer(std::streambuf* console_buffer, std::streambuf* file_buffer)
        : console_buffer_(console_buffer), file_buffer_(file_buffer)
    {
    }

protected:
    int overflow(int ch) override
    {
        if (ch == EOF) return !EOF;
        const int console_result = console_buffer_ ? console_buffer_->sputc(static_cast<char>(ch)) : ch;
        const int file_result = file_buffer_ ? file_buffer_->sputc(static_cast<char>(ch)) : ch;
        return (console_result == EOF || file_result == EOF) ? EOF : ch;
    }

    int sync() override
    {
        const int console_result = console_buffer_ ? console_buffer_->pubsync() : 0;
        const int file_result = file_buffer_ ? file_buffer_->pubsync() : 0;
        return (console_result == 0 && file_result == 0) ? 0 : -1;
    }

private:
    std::streambuf* console_buffer_ = nullptr;
    std::streambuf* file_buffer_ = nullptr;
};

class ConsoleLogTee {
public:
    explicit ConsoleLogTee(const std::filesystem::path& log_path)
    {
        log_file_.open(log_path, std::ios::out | std::ios::trunc | std::ios::binary);
        if (!log_file_.is_open()) {
            std::cerr << "[Log] Cannot open console log file: " << log_path.string() << std::endl;
            return;
        }

        cout_buffer_ = std::cout.rdbuf();
        cerr_buffer_ = std::cerr.rdbuf();
        tee_cout_ = std::make_unique<TeeStreamBuffer>(cout_buffer_, log_file_.rdbuf());
        tee_cerr_ = std::make_unique<TeeStreamBuffer>(cerr_buffer_, log_file_.rdbuf());
        std::cout.rdbuf(tee_cout_.get());
        std::cerr.rdbuf(tee_cerr_.get());
        enabled_ = true;

        std::cout << "[Log] Console output is also saved to: "
            << log_path.string() << std::endl;
    }

    ~ConsoleLogTee()
    {
        if (!enabled_) return;
        std::cout.flush();
        std::cerr.flush();
        std::cout.rdbuf(cout_buffer_);
        std::cerr.rdbuf(cerr_buffer_);
        log_file_.flush();
    }

private:
    bool enabled_ = false;
    std::ofstream log_file_;
    std::streambuf* cout_buffer_ = nullptr;
    std::streambuf* cerr_buffer_ = nullptr;
    std::unique_ptr<TeeStreamBuffer> tee_cout_;
    std::unique_ptr<TeeStreamBuffer> tee_cerr_;
};

// 作用：获取 exe 所在目录(失败时退回当前工作目录)。
std::filesystem::path GetExecutableDirectory()
{
    char module_path[MAX_PATH] = {};
    const DWORD length = GetModuleFileNameA(nullptr, module_path, MAX_PATH);
    if (length == 0 || length >= MAX_PATH) {
        return std::filesystem::current_path();
    }
    return std::filesystem::path(module_path).parent_path();
}

// ---- 全局算法参数(可按数据情况调整) ----
const float kSampleDensity = 30.0f;              // OSGB 缃戞牸閲囨牱瀵嗗害(鐐?骞虫柟绫?
const int kMaxPagedLODDepth = 10;                // PagedLOD 最大递归深度
const double kBoundarySupportBuffer = 3.0;      // 边界支撑点搜索的缓冲宽度(米)，仅作 KD-tree 粗搜索范围
const double kDensifyStep = 0.5;                 // 兜底加密边界时的步长(米)
const double kVerticalNormalMaxAbsZ = 0.45;      // 判定为墙面点的法向 z 分量阈值；超过则视为水平面点
                                                 // 超过则视为水平面点，在只取墙面点时被排除
const float kSupportVoxelLeaf = 0.15f;           // Support 2D voxel leaf size, meters.
const double kSupportDensityRadius = 1.0;        // XY radius for local support-density reliability weight.
const double kSupportDensityMinWeight = 0.40;    // Sparse support is downweighted, not discarded.
const double kSupportDensityMaxWeight = 1.15;    // Dense support gets only a mild boost.
const double kMinPolygonBBoxArea = 20.0;
// 小建筑阈值已移至 outlineRegular.cpp
// (kSmallBuildingSingleDirectionArea): small footprints no longer take a
// rectangle fast path, they run the normal pipeline with a forced single
// main direction.
const double kMinOutputPolygonArea = 15.0;      // Final write-out safety floor in square meters.
const double kMinOutputPolygonBBoxArea = 20.0;  // Final write-out bbox safety floor in square meters.
const double kModelCoverageBuffer = 3.0;
const std::size_t kMinModelEvidencePoints = 10;
const double kMinSharedBoundaryLength = 0.6;
const double kMinSharedPerimeterRatio = 0.025;
const double kShapeMergeMinSharedLength = 2.0;
const double kShapeMergeMinSharedPerimeterRatio = 0.015;
const double kStrongSeamEvidenceRatio = 0.85;  // 只有非常强的证据才阻止合并；调高以放宽被误分的大建筑
const double kMergeSeamLength = 3.0;           // 宽缝阈值：共享边 >= 此值才合并(窄缝连接不同栋，不合)
const double kLargeBuildingMergeArea = 1500.0; // 大建筑兜底合并：min(面积)>此值且共享边>=kMergeSeamLength 时
                                                // 跳过证据间直接合并(两栋各自独立的大楼几乎不可能共享几十米长边)
const double kMaxSplitCompactnessGainForMerge = 0.4;
const double kContainedFootprintMaxArea = 1000.0;
const double kOutputOverlapMinArea = 0.5;
const double kInitialMergeSimplifyTolerance = 0.35; // meters. Smooths raster/union seams before regularization.
const double kRobustSharedBoundaryTolerance = 0.65;
const double kRobustSharedAngleToleranceDeg = 8.0;
const double kOutputOverlapMinRatio = 0.005;
const double kSeamSampleStep = 0.5;
const double kSeamWallSearchRadius = 0.6;
const double kSeamRoofSideOffset = 1.2;
const double kSeamRoofSearchRadius = 0.8;
const double kSeamRoofHeightDifference = 1.2;
const int kMinSeamEvidenceSamples = 6;
const int kMinSeamHeightPairs = 3;
const double kSupportOwnershipTieTolerance = 0.20; // meters. Prevents dense-neighbor wall points from supporting both footprints.
const double kSupportOwnershipGridSize = 20.0;     // meters. Spatial index cell size for initial footprint ownership tests.
// ---- Effective wall-support association (dense-neighbour protection) ----
// kBoundarySupportBuffer stays the KD-tree COARSE search radius (tolerating
// local mask offset). Points entering direction estimation, support evidence
// and Ceres must additionally pass narrow-band / extended-band association,
// unambiguous ownership and wall-normal/edge-normal consistency.
const double kSupportNarrowBandDistance = 0.85;      // m, high-confidence wall band
const double kSupportExtendedBandDistance = 1.50;    // m, only via continuous chains
const double kSupportAmbiguousDistanceGap = 0.40;    // m, nearest-vs-second gap below this = ambiguous
const double kSupportNormalEdgeMaxAngleDeg = 25.0;   // deg, wall XY normal vs edge normal
const double kSupportExtendedMinRunLength = 2.0;     // m, chain length required for extended band
const double kSupportExtendedChainGap = 0.9;         // m, arc gap linking extended-band points
const double kSupportFallbackWeight = 0.30;          // fallback (non-wall) support weight
const double kSupportWideWallWeight = 0.50;          // wide-band wall tier weight (below effective)
const double kSupportDensifyWeight = 0.25;           // last-resort densified boundary weight
const double kNarrowNeckMaxWidth = 3.0;            // meters. Post-vectorization split threshold for missed building separations.
const double kNarrowNeckMinBoundarySeparation = 4.0;
const double kNarrowNeckBoundarySeparationRatio = 0.02;
const double kNarrowNeckMinPartArea = 20.0;
const double kNarrowNeckCutBuffer = 0.12;
const int kNarrowNeckMaxCutsPerFeature = 12;
// Experimental Mask-only branch. It is isolated from the OSGB pipeline and
// falls back to the existing VDP path whenever topology checks fail.
constexpr bool kUseTopologyPreservingResidualRegularization = true;
constexpr bool kMaskCurveDetectionDebugOnly = false;

// ---- 计时索引(秒，用于定位规则化各阶段耗时) ----
double g_supportTime = 0.0;   // 支撑点提取(含 KdTree 查询)累计
double g_optimizeTime = 0.0;
std::size_t g_removedSmallPolygons = 0;
std::size_t g_removedOutsideModel = 0;
std::size_t g_supportOwnershipRejected = 0;
std::size_t g_supportOwnershipAmbiguousRejected = 0;
std::size_t g_supportFilterCoarseWallTotal = 0;
std::size_t g_supportFilterEffectiveWallTotal = 0;
std::size_t g_supportFilterFallbackUsedTotal = 0;
std::size_t g_supportFilterWideWallTotal = 0;

// 作用：删除一个 Shapefile 的全部伴随文件(.shp/.shx/.dbf/.prj/.cpg 等)。返回是否删掉了主文件。
bool RemoveShapefileFamily(const std::filesystem::path& shpPath, bool verbose = true)
{
    static const char* const extensions[] = {
        ".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix",
        ".sbn", ".sbx", ".fbn", ".fbx", ".ain", ".aih", ".atx"
    };
    const std::filesystem::path stem = shpPath.parent_path() / shpPath.stem();
    for (const char* extension : extensions) {
        std::filesystem::path file = stem;
        file += extension;
        std::error_code ec;
        if (!std::filesystem::remove(file, ec) && ec) {
            if (!verbose) return false;
            std::cerr << "[SHP] Cannot remove " << file.string()
                      << ": " << ec.message() << std::endl;
            return false;
        }
    }
    return true;
}

// 作用：判断某个 Shapefile 家族(任一伴随文件)是否已存在。
bool ShapefileFamilyExists(const std::filesystem::path& shpPath)
{
    static const char* const extensions[] = {
        ".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix",
        ".sbn", ".sbx", ".fbn", ".fbx", ".ain", ".aih", ".atx"
    };
    const std::filesystem::path stem = shpPath.parent_path() / shpPath.stem();
    for (const char* extension : extensions) {
        std::filesystem::path file = stem;
        file += extension;
        std::error_code ec;
        if (std::filesystem::exists(file, ec)) return true;
    }
    return false;
}

// 作用：给路径加 _1/_2/... 后缀，返回第一个不存在的同族路径。
std::filesystem::path MakeUniqueShapefilePath(const std::filesystem::path& desired)
{
    const std::filesystem::path parent = desired.parent_path();
    const std::string stem = desired.stem().string();
    const std::string extension = desired.extension().empty()
        ? std::string(".shp")
        : desired.extension().string();
    for (int i = 1; i < 10000; ++i) {
        std::filesystem::path candidate = parent / (stem + "_" + std::to_string(i) + extension);
        if (!ShapefileFamilyExists(candidate)) return candidate;
    }
    return parent / (stem + "_new" + extension);
}

// 作用：优先删除旧族；文件被占用时改写到唯一后缀路径并提示。
std::filesystem::path PrepareWritableShapefilePath(
    const std::filesystem::path& desired,
    const std::string& label)
{
    if (RemoveShapefileFamily(desired, false)) return desired;

    const std::filesystem::path fallback = MakeUniqueShapefilePath(desired);
    std::cout << "[SHP] " << label << " is open or locked, writing to: "
              << fallback.string() << std::endl;
    return fallback;
}

// 作用：从 OSGB 的 metadata.xml 解析 SRSOrigin 地理参考偏移。
bool ReadMetadataOffset(const std::string& xmlPath, Eigen::Vector3d& offset)
{
    CPLXMLNode* tree = CPLParseXMLFile(xmlPath.c_str());
    if (!tree) {
        std::cerr << "Cannot parse metadata XML: " << xmlPath << std::endl;
        return false;
    }

    std::function<const char*(CPLXMLNode*)> findSRSOrigin = [&](CPLXMLNode* node) -> const char* {
        for (CPLXMLNode* current = node; current; current = current->psNext) {
            if (current->eType == CXT_Element &&
                current->pszValue &&
                std::string(current->pszValue) == "SRSOrigin") {
                return CPLGetXMLValue(current, nullptr, nullptr);
            }
            if (current->psChild) {
                if (const char* value = findSRSOrigin(current->psChild)) {
                    return value;
                }
            }
        }
        return nullptr;
    };

    const char* value = findSRSOrigin(tree);
    if (!value) {
        std::cerr << "metadata XML has no SRSOrigin: " << xmlPath << std::endl;
        CPLDestroyXMLNode(tree);
        return false;
    }

    double x = 0.0, y = 0.0, z = 0.0;
    const int parsed = std::sscanf(value, "%lf,%lf,%lf", &x, &y, &z);
    CPLDestroyXMLNode(tree);
    if (parsed < 2) {
        std::cerr << "Invalid SRSOrigin value in metadata XML: " << value << std::endl;
        return false;
    }

    offset = Eigen::Vector3d(x, y, parsed >= 3 ? z : 0.0);
    return true;
}

// 将采样点云保存为 LAS(点数据格式 3 = GPS+RGB，二进制 .las)。
// 内存点云为相对坐标，这里按 mc.offset 叠加成世界坐标后写入；laszip 按 header 的
// offset/scale 自动把 double 转为整型存储，读回时再加回。
bool SaveSampledCloudAsLas(const MyCloud& mc, const std::string& path)
{
    if (!mc.cloud || mc.cloud->empty()) {
        std::cerr << "[LAS] point cloud is empty, skip saving." << std::endl;
        return false;
    }

    laszip_POINTER writer = nullptr;
    if (laszip_create(&writer) != 0) {
        std::cerr << "[LAS] laszip_create failed." << std::endl;
        return false;
    }

    // 点数据格式 3(GPS time + RGB)，记录长度 34 字节
    if (laszip_set_point_type_and_size(writer, 3, 34) != 0) {
        laszip_CHAR* err = nullptr; laszip_get_error(writer, &err);
        std::cerr << "[LAS] laszip_set_point_type_and_size 失败: " << (err ? err : "?") << std::endl;
        laszip_destroy(writer);
        return false;
    }

    laszip_header_struct* header = nullptr;
    laszip_get_header_pointer(writer, &header);
    header->x_scale_factor = 0.0001;
    header->y_scale_factor = 0.0001;
    header->z_scale_factor = 0.0001;
    header->x_offset = mc.offset.x();
    header->y_offset = mc.offset.y();
    header->z_offset = mc.offset.z();
    header->number_of_point_records = static_cast<laszip_U32>(mc.cloud->size());

    // FALSE = 不压缩(.las)；TRUE 则为压缩(.laz)
    if (laszip_open_writer(writer, path.c_str(), FALSE) != 0) {
        laszip_CHAR* err = nullptr; laszip_get_error(writer, &err);
        std::cerr << "[LAS] laszip_open_writer 失败: " << (err ? err : "?") << std::endl;
        laszip_destroy(writer);
        return false;
    }

    laszip_point_struct* point = nullptr;
    laszip_get_point_pointer(writer, &point);

    const laszip_F64 ox = mc.offset.x();
    const laszip_F64 oy = mc.offset.y();
    const laszip_F64 oz = mc.offset.z();
    const std::size_t n = mc.cloud->size();
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = mc.cloud->points[i];
        // 世界坐标 = 相对坐标 + offset
        const laszip_F64 coords[3] = {
            static_cast<laszip_F64>(p.x) + ox,
            static_cast<laszip_F64>(p.y) + oy,
            static_cast<laszip_F64>(p.z) + oz
        };
        laszip_set_coordinates(writer, coords);
        // 8bit RGB -> 16bit LAS RGB(乘 257，与读取时的 /257 对称)
        point->rgb[0] = static_cast<laszip_U16>(p.r) * 257;
        point->rgb[1] = static_cast<laszip_U16>(p.g) * 257;
        point->rgb[2] = static_cast<laszip_U16>(p.b) * 257;
        point->classification = 1; // Unclassified
        laszip_set_point(writer, point);
        laszip_write_point(writer);
    }

    laszip_close_writer(writer);
    laszip_destroy(writer);
    return true;
}

// 作用：对支撑点云做 XY 体素降采样，返回聚合后的点与权重。
pcl::PointCloud<pcl::PointXYZ>::Ptr DownsampleSupport2D(
    const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
    double leaf)
{
    auto result = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);
    if (!cloud || cloud->empty() || leaf <= 0.0) return result;

    struct Accumulator {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        std::size_t count = 0;
    };
    std::map<std::pair<long long, long long>, Accumulator> cells;
    for (const auto& p : cloud->points) {
        const long long ix = static_cast<long long>(std::floor(static_cast<double>(p.x) / leaf));
        const long long iy = static_cast<long long>(std::floor(static_cast<double>(p.y) / leaf));
        auto& cell = cells[{ix, iy}];
        cell.x += p.x;
        cell.y += p.y;
        cell.z += p.z;
        ++cell.count;
    }

    result->reserve(cells.size());
    for (const auto& item : cells) {
        const auto& cell = item.second;
        pcl::PointXYZ p;
        p.x = static_cast<float>(cell.x / cell.count);
        p.y = static_cast<float>(cell.y / cell.count);
        p.z = static_cast<float>(cell.z / cell.count);
        result->push_back(p);
    }
    return result;
}

struct WeightedSupportCloud {
    pcl::PointCloud<pcl::PointXYZ>::Ptr points{new pcl::PointCloud<pcl::PointXYZ>};
    std::vector<double> weights;
};

// 作用：把数值钳制到 [lo, hi] 区间。
double ClampDoubleLocal(double value, double min_value, double max_value)
{
    return std::max(min_value, std::min(max_value, value));
}

// 作用：由法向 z 分量计算墙面可信度权重(越竖直权重越高)。
double NormalSupportWeightFromZ(double normal_z)
{
    const double verticality = ClampDoubleLocal(
        1.0 - std::abs(normal_z) / std::max(kVerticalNormalMaxAbsZ, 1e-6),
        0.0, 1.0);
    return ClampDoubleLocal(0.15 + 0.85 * verticality * verticality, 0.15, 1.0);
}

// 作用：带权重的 XY 体素降采样，每个体素输出代表点与平均权重。
WeightedSupportCloud DownsampleSupport2DWithWeights(
    const WeightedSupportCloud& support,
    double leaf)
{
    WeightedSupportCloud result;
    if (!support.points || support.points->empty() || leaf <= 0.0) return result;

    struct Representative {
        pcl::PointXYZ point;
        double weight = -1.0;
    };
    std::map<std::pair<long long, long long>, Representative> cells;
    for (std::size_t i = 0; i < support.points->size(); ++i) {
        const auto& p = support.points->points[i];
        const double weight = i < support.weights.size()
            ? ClampDoubleLocal(support.weights[i], 0.10, 1.0)
            : 1.0;
        const long long ix = static_cast<long long>(std::floor(static_cast<double>(p.x) / leaf));
        const long long iy = static_cast<long long>(std::floor(static_cast<double>(p.y) / leaf));
        auto& cell = cells[{ix, iy}];
        if (weight > cell.weight) {
            cell.point = p;
            cell.weight = weight;
        }
    }

    result.points->reserve(cells.size());
    result.weights.reserve(cells.size());
    for (const auto& item : cells) {
        result.points->push_back(item.second.point);
        result.weights.push_back(ClampDoubleLocal(item.second.weight, 0.10, 1.0));
    }
    return result;
}

// 作用：按局部点密度调整支撑点权重(稀疏降权、密集轻微加权)。
void ApplySupportDensityWeights(WeightedSupportCloud& support)
{
    if (!support.points || support.points->size() < 8) return;
    if (support.weights.size() != support.points->size()) {
        support.weights.assign(support.points->size(), 1.0);
    }

    pcl::PointCloud<pcl::PointXYZ>::Ptr support2d(new pcl::PointCloud<pcl::PointXYZ>);
    support2d->resize(support.points->size());
    for (std::size_t i = 0; i < support.points->size(); ++i) {
        const auto& p = support.points->points[i];
        support2d->points[i].x = p.x;
        support2d->points[i].y = p.y;
        support2d->points[i].z = 0.0f;
    }

    pcl::KdTreeFLANN<pcl::PointXYZ> densityTree;
    densityTree.setInputCloud(support2d);

    std::vector<int> neighborCounts(support.points->size(), 0);
    std::vector<int> indices;
    std::vector<float> distances;
    for (std::size_t i = 0; i < support2d->size(); ++i) {
        indices.clear();
        distances.clear();
        densityTree.radiusSearch(
            support2d->points[i],
            static_cast<float>(kSupportDensityRadius),
            indices,
            distances);
        neighborCounts[i] = std::max(0, static_cast<int>(indices.size()) - 1);
    }

    std::vector<int> sortedCounts = neighborCounts;
    std::nth_element(
        sortedCounts.begin(),
        sortedCounts.begin() + sortedCounts.size() / 2,
        sortedCounts.end());
    const double expectedCount = std::max(
        1.0,
        static_cast<double>(sortedCounts[sortedCounts.size() / 2]));

    double densityWeightSum = 0.0;
    double densityWeightMin = std::numeric_limits<double>::max();
    double densityWeightMax = 0.0;
    for (std::size_t i = 0; i < support.weights.size(); ++i) {
        const double raw = static_cast<double>(neighborCounts[i]) / expectedCount;
        const double densityWeight = ClampDoubleLocal(
            raw, kSupportDensityMinWeight, kSupportDensityMaxWeight);
        support.weights[i] = ClampDoubleLocal(
            support.weights[i] * densityWeight, 0.05, 1.15);
        densityWeightSum += densityWeight;
        densityWeightMin = std::min(densityWeightMin, densityWeight);
        densityWeightMax = std::max(densityWeightMax, densityWeight);
    }

    std::cerr << "[SupportDensity] radius=" << kSupportDensityRadius
              << " expected_neighbors=" << expectedCount
              << " avg_weight=" << densityWeightSum / support.weights.size()
              << " min=" << densityWeightMin
              << " max=" << densityWeightMax << std::endl;
}

// ===== Distance2D =====
// 作用：两个三维点的水平(XY)距离(忽略 Z)。
double Distance2D(const pcl::PointXYZ& a, const pcl::PointXYZ& b)
{
    const double dx = static_cast<double>(a.x) - b.x;
    const double dy = static_cast<double>(a.y) - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

// ===== PointToSegmentDistance2D =====
// 作用：点 p 到线段 a-b 的水平最短距离(垂足被 clamp 到线段范围内)。
double PointToSegmentDistance2D(const pcl::PointXYZ& p, const pcl::PointXYZ& a, const pcl::PointXYZ& b)
{
    const double vx = static_cast<double>(b.x) - a.x;
    const double vy = static_cast<double>(b.y) - a.y;
    const double wx = static_cast<double>(p.x) - a.x;
    const double wy = static_cast<double>(p.y) - a.y;
    const double len2 = vx * vx + vy * vy;
    double t = len2 > 0.0 ? (wx * vx + wy * vy) / len2 : 0.0;
    t = std::max(0.0, std::min(1.0, t));
    const double px = static_cast<double>(a.x) + t * vx;
    const double py = static_cast<double>(a.y) + t * vy;
    const double dx = static_cast<double>(p.x) - px;
    const double dy = static_cast<double>(p.y) - py;
    return std::sqrt(dx * dx + dy * dy);
}

// ===== PointInPolygon2D =====
// 作用：判断点 p 是否在多边形 poly 内部(水平面，射线法)。
bool PointInPolygon2D(const pcl::PointXYZ& p, const std::vector<pcl::PointXYZ>& poly)
{
    bool inside = false;
    const std::size_t n = poly.size();
    if (n < 3) return false;

    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const auto& pi = poly[i];
        const auto& pj = poly[j];
        const bool intersect = ((pi.y > p.y) != (pj.y > p.y)) &&
            (p.x < (pj.x - pi.x) * (p.y - pi.y) / ((pj.y - pi.y) + 1e-12f) + pi.x);
        if (intersect) inside = !inside;
    }
    return inside;
}

// 作用：点到多边形边界的水平最短距离。
double DistanceToRingBoundary2D(const pcl::PointXYZ& p, const std::vector<pcl::PointXYZ>& ring)
{
    if (ring.empty()) return std::numeric_limits<double>::max();
    double best = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < ring.size(); ++i) {
        best = std::min(best, PointToSegmentDistance2D(p, ring[i], ring[(i + 1) % ring.size()]));
    }
    return best;
}

struct RingEnvelope2D {
    double minX = 0.0;
    double maxX = 0.0;
    double minY = 0.0;
    double maxY = 0.0;
};

RingEnvelope2D ComputeRingEnvelope2D(const std::vector<pcl::PointXYZ>& ring)
{
    RingEnvelope2D env;
    if (ring.empty()) return env;
    env.minX = env.maxX = ring.front().x;
    env.minY = env.maxY = ring.front().y;
    for (const auto& p : ring) {
        env.minX = std::min(env.minX, static_cast<double>(p.x));
        env.maxX = std::max(env.maxX, static_cast<double>(p.x));
        env.minY = std::min(env.minY, static_cast<double>(p.y));
        env.maxY = std::max(env.maxY, static_cast<double>(p.y));
    }
    return env;
}

// 作用：判断点是否在包围盒(可加边距)内。
bool EnvelopeContainsPoint(const RingEnvelope2D& env, const pcl::PointXYZ& p, double buffer)
{
    return p.x >= env.minX - buffer && p.x <= env.maxX + buffer &&
           p.y >= env.minY - buffer && p.y <= env.maxY + buffer;
}

// ===== RemoveClosingDuplicate =====
// 作用：如果多边形首尾点几乎重合，去掉末尾的重复闭合点。
void RemoveClosingDuplicate(std::vector<pcl::PointXYZ>& points)
{
    if (points.size() >= 2 && Distance2D(points.front(), points.back()) < 1e-6) {
        points.pop_back();
    }
}

// ===== ExtractExteriorRing =====
// 作用：从 OGR 多边形中取出外环顶点，转为 pcl::PointXYZ 点列，并去闭合重复点。
std::vector<pcl::PointXYZ> ExtractExteriorRing(
    OGRPolygon* polygon, const Eigen::Vector3d& metadataOffset)
{
    std::vector<pcl::PointXYZ> ring;
    if (!polygon) return ring;

    OGRLinearRing* ogr_ring = polygon->getExteriorRing();
    if (!ogr_ring) return ring;

    for (int i = 0; i < ogr_ring->getNumPoints(); ++i) {
        pcl::PointXYZ p;
        p.x = static_cast<float>(ogr_ring->getX(i) - metadataOffset.x());
        p.y = static_cast<float>(ogr_ring->getY(i) - metadataOffset.y());
        p.z = static_cast<float>(ogr_ring->getZ(i) - metadataOffset.z());
        ring.push_back(p);
    }
    RemoveClosingDuplicate(ring);
    return ring;
}

struct InitialRingRecord {
    std::size_t id = 0;
    GIntBig sourceFid = OGRNullFID;
    int ringIndex = 0;
    std::vector<pcl::PointXYZ> ring;
    RingEnvelope2D envelope;
};

struct CellKey {
    int x = 0;
    int y = 0;

    bool operator==(const CellKey& other) const {
        return x == other.x && y == other.y;
    }
};

struct CellKeyHash {
    std::size_t operator()(const CellKey& key) const {
        const std::uint64_t ux = static_cast<std::uint32_t>(key.x);
        const std::uint64_t uy = static_cast<std::uint32_t>(key.y);
        return static_cast<std::size_t>((ux << 32) ^ uy);
    }
};

struct SupportOwnershipContext {
    std::vector<InitialRingRecord> records;
    std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> grid;
    double cellSize = kSupportOwnershipGridSize;

    CellKey cellFor(double x, double y) const {
        return {
            static_cast<int>(std::floor(x / cellSize)),
            static_cast<int>(std::floor(y / cellSize))
        };
    }

    void addRecord(InitialRingRecord record) {
        if (record.ring.size() < 3) return;
        record.id = records.size();
        record.envelope = ComputeRingEnvelope2D(record.ring);
        const std::size_t id = record.id;
        records.push_back(std::move(record));
        const auto minCell = cellFor(records.back().envelope.minX - kBoundarySupportBuffer,
                                     records.back().envelope.minY - kBoundarySupportBuffer);
        const auto maxCell = cellFor(records.back().envelope.maxX + kBoundarySupportBuffer,
                                     records.back().envelope.maxY + kBoundarySupportBuffer);
        for (int cy = minCell.y; cy <= maxCell.y; ++cy) {
            for (int cx = minCell.x; cx <= maxCell.x; ++cx) {
                grid[{cx, cy}].push_back(id);
            }
        }
    }

    // 三态归属：一个点位于两个轮廓之间，且到两者的距离
    // 差值小于 kSupportAmbiguousDistanceGap 时为"归属不明确"，必须
    // 禁止进入任何一栋的独立规则化(旧的二值
    // ownsSupportPoint 按 ring id 强制分配这类点，导致一片混合
    // 点云被交错地分给两个轮廓)。
    enum class SupportOwnership { Owned, Ambiguous, OwnedByOther };

    SupportOwnership classifySupportPoint(const pcl::PointXYZ& p,
                                          std::size_t currentId,
                                          double currentDistance) const {
        if (currentId >= records.size()) return SupportOwnership::Owned;
        const auto cell = cellFor(p.x, p.y);
        const auto found = grid.find(cell);
        if (found == grid.end()) return SupportOwnership::Owned;

        double nearestOther = std::numeric_limits<double>::max();
        for (std::size_t otherId : found->second) {
            if (otherId == currentId || otherId >= records.size()) continue;
            const auto& other = records[otherId];
            if (!EnvelopeContainsPoint(other.envelope, p, kBoundarySupportBuffer)) continue;
            const double otherDistance = DistanceToRingBoundary2D(p, other.ring);
            if (otherDistance > kBoundarySupportBuffer) continue;
            nearestOther = std::min(nearestOther, otherDistance);
        }
        if (nearestOther < currentDistance) {
            return SupportOwnership::OwnedByOther;
        }
        if (nearestOther - currentDistance < kSupportAmbiguousDistanceGap) {
            return SupportOwnership::Ambiguous;
        }
        return SupportOwnership::Owned;
    }

    bool ownsSupportPoint(const pcl::PointXYZ& p,
                          std::size_t currentId,
                          double currentDistance) const {
        return classifySupportPoint(p, currentId, currentDistance) ==
               SupportOwnership::Owned;
    }
};

// 作用：按 sourceFid+ringIndex 在归属索引里查环 id。
std::size_t FindInitialRingRecordId(
    const SupportOwnershipContext* ownership,
    GIntBig sourceFid,
    int ringIndex)
{
    if (!ownership) return static_cast<std::size_t>(-1);
    for (const auto& record : ownership->records) {
        if (record.sourceFid == sourceFid && record.ringIndex == ringIndex) {
            return record.id;
        }
    }
    return static_cast<std::size_t>(-1);
}

// ===== DensifyBoundary =====
// 作用：沿多边形边界按步长 step 插值加密，生成密集点云作为支撑点的兜底来源。
pcl::PointCloud<pcl::PointXYZ>::Ptr DensifyBoundary(const std::vector<pcl::PointXYZ>& ring, double step)
{
    auto cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);
    if (ring.size() < 2) return cloud;

    for (std::size_t i = 0; i < ring.size(); ++i) {
        const pcl::PointXYZ& a = ring[i];
        const pcl::PointXYZ& b = ring[(i + 1) % ring.size()];
        const double len = Distance2D(a, b);
        const int segments = std::max(1, static_cast<int>(std::ceil(len / step)));
        for (int s = 0; s < segments; ++s) {
            const double t = static_cast<double>(s) / segments;
            pcl::PointXYZ p;
            p.x = static_cast<float>(a.x * (1.0 - t) + b.x * t);
            p.y = static_cast<float>(a.y * (1.0 - t) + b.y * t);
            p.z = static_cast<float>(a.z * (1.0 - t) + b.z * t);
            cloud->push_back(p);
        }
    }
    return cloud;
}

// ===== MakePolygon =====
// 作用：把点列重新组装成一个 OGRPolygon(自动闭合外环)。点数<3 返回 nullptr。
OGRPolygon* MakePolygon(
    const std::vector<pcl::PointXYZ>& ring, const Eigen::Vector3d& metadataOffset)
{
    if (ring.size() < 3) return nullptr;

    auto* polygon = new OGRPolygon();
    OGRLinearRing ogr_ring;
    for (const auto& p : ring) {
        ogr_ring.addPoint(p.x + metadataOffset.x(),
            p.y + metadataOffset.y(),
            p.z + metadataOffset.z());
    }
    ogr_ring.addPoint(ring.front().x + metadataOffset.x(),
        ring.front().y + metadataOffset.y(),
        ring.front().z + metadataOffset.z());
    polygon->addRing(&ogr_ring);
    return polygon;
}

// ===== ExtractBoundarySupportFromOSGB =====
// 作用：从采样的点云里，挑出支撑边界规则化用的点。
//       [KdTree 版]用预先建好的 kdtree 做一次 2D 半径查询，只取多边形扩展包围盒范围内的
//       候选点，再按模式筛选(有效墙面窄带/宽带墙面/非墙面兜底)。
// 参数：sampled - OSGB 采样点云；ring - 多边形边界；mode - 提取模式；
//       kdtree - 预建的 2D KdTree(Z 已置 0)。返回：支撑点云(点+权重)。
// 有效墙面支撑管线的逐栋过滤统计。
struct SupportFilterStats {
    std::size_t coarseWall = 0;        // 3 m coarse candidates passing the wall-normal gate
    std::size_t narrowBand = 0;        // effective points inside the narrow band
    std::size_t extendedBand = 0;      // effective points from surviving chains
    std::size_t extendedRaw = 0;       // extended-band candidates before chain filter
    std::size_t ambiguousRejected = 0; // ownership == Ambiguous
    std::size_t otherOwnerRejected = 0;
    std::size_t normalRejected = 0;    // wall/edge normal inconsistency or degenerate XY normal
    std::size_t tRejected = 0;         // projection parameter outside [0.05, 0.95]
};

// 点到线段距离，同时返回原始投影参数(未钳制)。
double PointSegmentDistanceWithT(
    const pcl::PointXYZ& p, const pcl::PointXYZ& a, const pcl::PointXYZ& b, double& t);

// 提取模式：
//  EffectiveWall — narrow band (<=0.85 m) + chain-verified extended band
//                  (<=1.5 m), unambiguous ownership, normal/edge consistency.
//  WideWall      — the legacy 3 m wall band (wall normals + unambiguous
//                  ownership only): fallback tier for footprints whose mask
//                  offset exceeds the narrow/extended bands, so their edges
//                  still get REAL wall anchoring instead of roof points.
//  NonWall       — 3 m / interior points without the wall-normal gate
//                      (Ceres 最后兜底；绝不参与方向估计)。
enum class SupportExtractMode { EffectiveWall, WideWall, NonWall };

// 作用：按模式从采样点云提取边界支撑点(有效墙面/宽带墙面/非墙面)，含归属三态判定、法向-边一致性校验与扩展带链过滤。
WeightedSupportCloud ExtractWeightedBoundarySupportFromOSGB(
    const MyCloudPtr& sampled,
    const std::vector<pcl::PointXYZ>& ring,
    SupportExtractMode mode,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree,
    const SupportOwnershipContext* ownership = nullptr,
    std::size_t currentRingId = static_cast<std::size_t>(-1),
    SupportFilterStats* filterStats = nullptr)
{
    WeightedSupportCloud support;
    if (!sampled || !sampled->cloud || sampled->cloud->empty() || ring.size() < 3 || !kdtree) return support;

    double min_x = ring[0].x, max_x = ring[0].x;
    double min_y = ring[0].y, max_y = ring[0].y;
    for (const auto& p : ring) {
        min_x = std::min(min_x, static_cast<double>(p.x));
        max_x = std::max(max_x, static_cast<double>(p.x));
        min_y = std::min(min_y, static_cast<double>(p.y));
        max_y = std::max(max_y, static_cast<double>(p.y));
    }

    // 以扩展包围盒中心为圆心、半对角线为半径做一次 2D 查询，得到候选点索引
    const double cx = (min_x + max_x) * 0.5;
    const double cy = (min_y + max_y) * 0.5;
    const double half_w = (max_x - min_x) * 0.5 + kBoundarySupportBuffer;
    const double half_h = (max_y - min_y) * 0.5 + kBoundarySupportBuffer;
    const float radius = static_cast<float>(std::sqrt(half_w * half_w + half_h * half_h) + 1e-3);

    pcl::PointXYZ center;
    center.x = static_cast<float>(cx);
    center.y = static_cast<float>(cy);
    center.z = 0.0f;

    std::vector<int> candidates;
    std::vector<float> dists2;
    kdtree->radiusSearch(center, radius, candidates, dists2);

    // 最近边关联所需的辅助量：前缀弧长 + XY 边法向。
    const std::size_t n = ring.size();
    std::vector<double> prefix(n + 1, 0.0);
    std::vector<double> edgeNx(n, 0.0);
    std::vector<double> edgeNy(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = ring[i];
        const auto& b = ring[(i + 1) % n];
        prefix[i + 1] = prefix[i] + Distance2D(a, b);
        const double len = Distance2D(a, b);
        if (len > 1e-9) {
            edgeNx[i] = -(b.y - a.y) / len;
            edgeNy[i] = (b.x - a.x) / len;
        }
    }
    const double perimeter = prefix[n];
    const double normalEdgeCos =
        std::cos(kSupportNormalEdgeMaxAngleDeg * M_PI / 180.0);
    struct ExtendedCandidate {
        pcl::PointXYZ p;
        double weight = 1.0;
        double arc = 0.0;
    };
    std::vector<ExtendedCandidate> extendedCandidates;
    auto bump = [&filterStats](std::size_t SupportFilterStats::*field) {
        if (filterStats) ++(filterStats->*field);
    };

    for (int idx : candidates) {
        if (idx < 0 || static_cast<std::size_t>(idx) >= sampled->cloud->size()) continue;
        const auto& src = sampled->cloud->points[idx];

        // 包围盒粗筛保留，保证结果与全量扫描一致
        if (src.x < min_x - kBoundarySupportBuffer || src.x > max_x + kBoundarySupportBuffer ||
            src.y < min_y - kBoundarySupportBuffer || src.y > max_y + kBoundarySupportBuffer) {
            continue;
        }

        // 只取墙面点时：法向 z 分量过大(接近平面)的点跳过
        const bool wallMode =
            mode == SupportExtractMode::EffectiveWall || mode == SupportExtractMode::WideWall;
        double normalWeight = 1.0;
        if (sampled->normal && static_cast<std::size_t>(idx) < sampled->normal->size()) {
            const auto& n = sampled->normal->points[idx];
            const double absNormalZ = std::abs(n.normal_z);
            if (wallMode && absNormalZ > kVerticalNormalMaxAbsZ) {
                continue;
            }
            if (wallMode) {
                normalWeight = NormalSupportWeightFromZ(n.normal_z);
            }
        }

        pcl::PointXYZ p;
        p.x = src.x;
        p.y = src.y;
        p.z = src.z;

        // 最近边、距离和投影参数。
        double boundaryDistance = std::numeric_limits<double>::max();
        std::size_t bestEdge = n;
        double bestT = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double t = 0.0;
            const double d = PointSegmentDistanceWithT(p, ring[i], ring[(i + 1) % n], t);
            if (d < boundaryDistance) {
                boundaryDistance = d;
                bestEdge = i;
                bestT = t;
            }
        }
        const bool nearBoundary =
            bestEdge < n && boundaryDistance <= kBoundarySupportBuffer;

        if (mode == SupportExtractMode::EffectiveWall) {
            // ---- Effective wall-support pipeline ----
            if (!nearBoundary) continue;
            bump(&SupportFilterStats::coarseWall);
            if (ownership) {
                const auto own = ownership->classifySupportPoint(
                    p, currentRingId, boundaryDistance);
                if (own == SupportOwnershipContext::SupportOwnership::OwnedByOther) {
                    ++g_supportOwnershipRejected;
                    bump(&SupportFilterStats::otherOwnerRejected);
                    continue;
                }
                if (own == SupportOwnershipContext::SupportOwnership::Ambiguous) {
                    ++g_supportOwnershipAmbiguousRejected;
                    bump(&SupportFilterStats::ambiguousRejected);
                    continue;
                }
            }
            if (bestT < 0.05 || bestT > 0.95) {
                bump(&SupportFilterStats::tRejected);
                continue;
            }
            // 墙面 XY 法向与最近边法向的一致性校验。
            // XY 法向退化的点同样不算有效墙面支撑。
            bool normalOk = true;
            if (sampled->normal &&
                static_cast<std::size_t>(idx) < sampled->normal->size()) {
                const auto& nv = sampled->normal->points[idx];
                const double nLen = std::hypot(nv.normal_x, nv.normal_y);
                if (nLen < 1e-6) {
                    normalOk = false;
                } else {
                    const double dot = std::abs(
                        (nv.normal_x * edgeNx[bestEdge] +
                         nv.normal_y * edgeNy[bestEdge]) / nLen);
                    if (dot < normalEdgeCos) normalOk = false;
                }
            }
            if (!normalOk) {
                bump(&SupportFilterStats::normalRejected);
                continue;
            }
            if (boundaryDistance <= kSupportNarrowBandDistance) {
                support.points->push_back(p);
                support.weights.push_back(normalWeight);
                bump(&SupportFilterStats::narrowBand);
            } else if (boundaryDistance <= kSupportExtendedBandDistance) {
                ExtendedCandidate c;
                c.p = p;
                c.weight = normalWeight;
                c.arc = prefix[bestEdge] +
                    bestT * (prefix[bestEdge + 1] - prefix[bestEdge]);
                extendedCandidates.push_back(c);
                bump(&SupportFilterStats::extendedRaw);
            }
            continue;
        }

        if (mode == SupportExtractMode::WideWall) {
            // ---- Legacy 3 m wall band: wall normals + unambiguous ownership
            // only (no band/t/normal-edge gates). Restores real-wall
            // anchoring for footprints whose mask offset exceeds the
            // narrow/extended bands. Weights are lowered by the caller. ----
            if (!nearBoundary) continue;
            if (ownership) {
                const auto own = ownership->classifySupportPoint(
                    p, currentRingId, boundaryDistance);
                if (own != SupportOwnershipContext::SupportOwnership::Owned) {
                    if (own == SupportOwnershipContext::SupportOwnership::Ambiguous) {
                        ++g_supportOwnershipAmbiguousRejected;
                    } else {
                        ++g_supportOwnershipRejected;
                    }
                    continue;
                }
            }
            support.points->push_back(p);
            support.weights.push_back(normalWeight);
            continue;
        }

        // ---- Non-wall fallback extraction (wide band allowed, ambiguous
        // ownership still excluded; weights are lowered by the caller) ----
        const bool inside = PointInPolygon2D(p, ring);
        if (nearBoundary || inside) {
            if (nearBoundary && ownership) {
                const auto own = ownership->classifySupportPoint(
                    p, currentRingId, boundaryDistance);
                if (own == SupportOwnershipContext::SupportOwnership::OwnedByOther) {
                    ++g_supportOwnershipRejected;
                    bump(&SupportFilterStats::otherOwnerRejected);
                    continue;
                }
                if (own == SupportOwnershipContext::SupportOwnership::Ambiguous) {
                    ++g_supportOwnershipAmbiguousRejected;
                    bump(&SupportFilterStats::ambiguousRejected);
                    continue;
                }
            }
            support.points->push_back(p);
            support.weights.push_back(normalWeight);
        }
    }

    // 扩展带链过滤：只保留投影跨度能证明存在
    // 连续墙面支撑的链(跨环缝的链会被
    // 切开——最小实现的取舍)。
    if (mode == SupportExtractMode::EffectiveWall && !extendedCandidates.empty()) {
        std::sort(extendedCandidates.begin(), extendedCandidates.end(),
                  [](const ExtendedCandidate& a, const ExtendedCandidate& b) {
                      return a.arc < b.arc;
                  });
        std::size_t start = 0;
        for (std::size_t k = 1; k <= extendedCandidates.size(); ++k) {
            const bool chainEnd =
                k == extendedCandidates.size() ||
                extendedCandidates[k].arc - extendedCandidates[k - 1].arc >
                    kSupportExtendedChainGap;
            if (!chainEnd) continue;
            const double span =
                extendedCandidates[k - 1].arc - extendedCandidates[start].arc;
            if (span >= kSupportExtendedMinRunLength) {
                for (std::size_t m = start; m < k; ++m) {
                    support.points->push_back(extendedCandidates[m].p);
                    support.weights.push_back(extendedCandidates[m].weight);
                    if (filterStats) ++filterStats->extendedBand;
                }
            }
            start = k;
        }
    }

    return support;
}

struct DebugHypothesisRecord {
    GIntBig sourceFid = OGRNullFID;
    int ringIndex = 0;
    std::vector<pcl::PointXYZ> points;
};

struct RegularizationDebugCollector {
    std::vector<DebugHypothesisRecord> hypotheses;
    PointCloudT::Ptr support{new PointCloudT};
};

double BoundingBoxArea2D(const std::vector<pcl::PointXYZ>& ring);
double PolygonArea2D(const std::vector<pcl::PointXYZ>& ring);
std::vector<pcl::PointXYZ> OrientedBoundingRectangle(
    const std::vector<pcl::PointXYZ>& ring, double angle);
bool OutputRingPassesSizeFloor(const std::vector<pcl::PointXYZ>& ring);

// ===== RegularizeRing =====
// 作用：对单条多边形边界做规则化。
//       优先用有效墙面支撑点，不足退化为宽带墙面/全部点，再不足用边界自身加密点兜底；
//       然后构造 outlineRegular 并执行 regular_Contour() 得到规则化结果。
// ===== Support evidence (pseudo-footprint filter) =====
// 掩膜尾巴经窄颈切分后可能形成大于 20m² 的独立多边形，
// 背后却没有真实建筑。证据必须来自"原始"墙面专用
// 支撑集(提取后立即计算，先于任何回退、
// 加密或体素降采样)——最终的
// "support" 数混入了邻楼点和边界加密点，
// 不能作为建筑是否存在的依据。
struct SupportEvidence {
    std::size_t wallPointCount = 0;
    std::size_t associatedWallPointCount = 0;
    std::size_t coveredBins = 0;
    std::size_t totalBins = 0;
    double boundaryCoverage = 0.0;
    double maxEmptyGap = 0.0;        // metres, circular longest unoccupied arc
    double perimeter = 0.0;
    double associatedDensity = 0.0;  // associated points per metre
    bool hasEvidence = false;
};

constexpr bool kSupportEvidenceFilterDryRun = true;   // calibration phase: log only
constexpr double kSupportEvidenceSmallArea = 100.0;         // m2
constexpr double kSupportEvidenceMinCoverage = 0.25;        // occupied-bin fraction
constexpr double kSupportEvidenceMaxGapRatio = 0.50;        // x perimeter
constexpr std::size_t kSupportEvidenceMinAssociatedPoints = 8;
constexpr double kSupportEvidenceBinLength = 1.8;           // m per boundary bin
constexpr double kSupportEvidenceAssociationDistance = 0.8; // m, point-to-edge

// 点到线段距离，同时输出原始投影参数
// (未钳制)，供 [0.05, 0.95] 关联测试使用。
double PointSegmentDistanceWithT(
    const pcl::PointXYZ& p, const pcl::PointXYZ& a, const pcl::PointXYZ& b, double& t)
{
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double lenSq = dx * dx + dy * dy;
    if (lenSq < 1e-12) {
        t = 0.0;
        return std::hypot(p.x - a.x, p.y - a.y);
    }
    t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / lenSq;
    const double tc = std::max(0.0, std::min(1.0, t));
    return std::hypot(p.x - a.x - tc * dx, p.y - a.y - tc * dy);
}

// 作用：在(预过滤的)有效墙面支撑上计算边界覆盖率/最大空缺/关联密度等存在性证据。
SupportEvidence ComputeSupportEvidence(
    const std::vector<pcl::PointXYZ>& ring,
    const WeightedSupportCloud& wallOnlySupport)
{
    SupportEvidence evidence;
    const std::size_t n = ring.size();
    if (n < 3 || !wallOnlySupport.points || wallOnlySupport.points->empty()) {
        return evidence;
    }
    std::vector<double> prefix(n + 1, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        prefix[i + 1] = prefix[i] + Distance2D(ring[i], ring[(i + 1) % n]);
    }
    const double perimeter = prefix[n];
    evidence.perimeter = perimeter;
    evidence.wallPointCount = wallOnlySupport.points->size();
    if (perimeter < 2.0) return evidence;   // degenerate guard

    const std::size_t totalBins = std::max<std::size_t>(
        4, static_cast<std::size_t>(std::ceil(perimeter / kSupportEvidenceBinLength)));
    std::vector<bool> occupied(totalBins, false);

    // 证据输入是已预过滤的有效墙面支撑(窄带或
    // 已验证的扩展链)，关联只是再次确认分桶。
    const double associationDistance = kSupportExtendedBandDistance;
    std::size_t associated = 0;
    for (const auto& p : *wallOnlySupport.points) {
        double bestDistance = std::numeric_limits<double>::max();
        std::size_t bestEdge = n;
        double bestT = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double t = 0.0;
            const double d = PointSegmentDistanceWithT(
                p, ring[i], ring[(i + 1) % n], t);
            if (d < bestDistance) {
                bestDistance = d;
                bestEdge = i;
                bestT = t;
            }
        }
        if (bestEdge >= n) continue;
        if (bestDistance > associationDistance) continue;
        // 忽略紧贴顶点的投影：不携带任何边级证据。
        if (bestT < 0.05 || bestT > 0.95) continue;
        ++associated;
        const double arc = prefix[bestEdge] +
            bestT * (prefix[bestEdge + 1] - prefix[bestEdge]);
        const std::size_t bin = std::min(
            totalBins - 1,
            static_cast<std::size_t>(std::floor(arc / perimeter * totalBins)));
        occupied[bin] = true;
    }
    evidence.associatedWallPointCount = associated;
    evidence.totalBins = totalBins;
    evidence.coveredBins = static_cast<std::size_t>(
        std::count(occupied.begin(), occupied.end(), true));
    evidence.boundaryCoverage =
        static_cast<double>(evidence.coveredBins) / static_cast<double>(totalBins);

    // 未占用桶的最长环形连续段，换算回弧长。
    if (evidence.coveredBins == totalBins) {
        evidence.maxEmptyGap = 0.0;
    } else {
        std::size_t bestRun = 0;
        std::size_t run = 0;
        bool wrapped = false;
        for (std::size_t k = 0; k < 2 * totalBins; ++k) {
            if (!occupied[k % totalBins]) {
                ++run;
                if (run > bestRun) bestRun = run;
            } else {
                if (k >= totalBins) break;
                if (k == 0) continue;
                run = 0;
            }
        }
        bestRun = std::min(bestRun, totalBins - evidence.coveredBins);
        evidence.maxEmptyGap =
            static_cast<double>(bestRun) * perimeter / static_cast<double>(totalBins);
    }
    evidence.associatedDensity =
        static_cast<double>(associated) / std::max(perimeter, 1.0);
    evidence.hasEvidence =
        associated >= kSupportEvidenceMinAssociatedPoints ||
        evidence.boundaryCoverage >= kSupportEvidenceMinCoverage;
    return evidence;
}

// 保守的三级判定。第三级保留其余所有：大
// 建筑绝不在此时删除(无模型覆盖的要素
// 已被上游 RegularizeGeometry 的 RingHasModelCoverage 过滤)。
bool ShouldDropForLackOfSupport(
    const SupportEvidence& evidence,
    double ringArea,
    std::string& reason)
{
    if (evidence.associatedWallPointCount == 0 &&
        evidence.boundaryCoverage < 0.10 &&
        evidence.maxEmptyGap > 0.60 * evidence.perimeter &&
        ringArea < kSupportEvidenceSmallArea) {
        reason = "no_support";
        return true;
    }
    if (ringArea < kSupportEvidenceSmallArea &&
        evidence.associatedWallPointCount < kSupportEvidenceMinAssociatedPoints &&
        evidence.boundaryCoverage < kSupportEvidenceMinCoverage &&
        evidence.maxEmptyGap > kSupportEvidenceMaxGapRatio * evidence.perimeter) {
        reason = "small_low_coverage";
        return true;
    }
    if (ringArea < kSupportEvidenceSmallArea && !evidence.hasEvidence) {
        reason = "low_confidence_keep";
    } else {
        reason = "ok";
    }
    return false;
}

// 作用：OSGB 模式的单环规则化入口：四级支撑阶梯(effective->宽带墙面->非墙面->加密) + 方向提示 + outlineRegular 全管线。
std::vector<pcl::PointXYZ> RegularizeRing(
    const std::vector<pcl::PointXYZ>& ring,
    const MyCloudPtr& sampled,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree,
    const SupportOwnershipContext* ownership = nullptr,
    std::size_t currentRingId = static_cast<std::size_t>(-1),
    std::vector<pcl::PointXYZ>* debugBestHypothesis = nullptr,
    std::vector<pcl::PointXYZ>* debugSupport = nullptr,
    long long sourceFid = -1,
    SupportEvidence* supportEvidence = nullptr)
{
    if (ring.size() < 3) return ring;

    // ---- 有效墙面支撑提取(3m 粗搜 + 窄带/
    // 扩展带 + 归属 + 法向-边过滤在内部完成) ----
    SupportFilterStats filterStats;
    auto wallOnlySupport = ExtractWeightedBoundarySupportFromOSGB(
        sampled, ring, SupportExtractMode::EffectiveWall,
        kdtree, ownership, currentRingId, &filterStats);
    // 原始有效集快照：下面的降采样会替换
    // weightedSupport 的点云，而方向提示必须使用原始的
    // 有效墙面点。
    const WeightedSupportCloud effectiveWallSupportRaw = wallOnlySupport;
    g_supportFilterCoarseWallTotal += filterStats.coarseWall;
    g_supportFilterEffectiveWallTotal +=
        filterStats.narrowBand + filterStats.extendedBand;

    // ---- 支撑点提取(计时) ----
    auto t0 = std::chrono::steady_clock::now();
    // ---- 伪轮廓支撑证据：在"原始"墙面专用支撑集上计算
    // (先于任何回退/加密/降采样)。判定删除时
    // 返回空环，要素不会进入 VDP/Ceres。
    {
        SupportEvidence evidence = ComputeSupportEvidence(ring, wallOnlySupport);
        if (supportEvidence) *supportEvidence = evidence;
        const double ringArea = PolygonArea2D(ring);
        std::string reason;
        const bool drop = ShouldDropForLackOfSupport(evidence, ringArea, reason);
        // 此处 source_fid 是稳定属性 id(栅格模式下等于
        // 初始轮廓的 "id" 字段，否则为 OGR FID)，见 main() 循环。
        std::cerr << "[SupportEvidence] source_fid=" << sourceFid
                  << " evidence_source=effective_wall_support"
                  << " area=" << ringArea
                  << " perimeter=" << evidence.perimeter
                  << " wall_points=" << evidence.wallPointCount
                  << " associated_wall_points=" << evidence.associatedWallPointCount
                  << " coverage=" << evidence.boundaryCoverage
                  << " max_gap=" << evidence.maxEmptyGap
                  << " density=" << evidence.associatedDensity
                  << " decision=" << (drop ? (kSupportEvidenceFilterDryRun
                                                 ? "would_drop" : "drop")
                                           : "keep")
                  << " reason=" << reason << std::endl;
        if (ringArea < kSupportEvidenceSmallArea) {
            std::cerr << "[SupportEvidence] small_area_detail source_fid=" << sourceFid
                      << " covered_bins=" << evidence.coveredBins
                      << "/" << evidence.totalBins << std::endl;
        }
        if (drop && !kSupportEvidenceFilterDryRun) {
            return {};
        }
    }

    auto weightedSupport = wallOnlySupport;
    std::size_t wideWallCount = 0;
    std::size_t fallbackCount = 0;
    std::size_t densifyCount = 0;
    if (weightedSupport.points->size() < 20) {
        // 宽带墙面层：旧版 3m 墙面带(无歧义归属、
        // 墙面法向)。掩膜偏移超出窄带/
        // 扩展带的轮廓在这里仍能获得真实墙面锚定，
        // 而不是屋面点，从而把边拉回墙面。
        auto wideWallSupport = ExtractWeightedBoundarySupportFromOSGB(
            sampled, ring, SupportExtractMode::WideWall,
            kdtree, ownership, currentRingId);
        wideWallCount = wideWallSupport.points->size();
        g_supportFilterWideWallTotal += wideWallCount;
        for (std::size_t i = 0; i < wideWallSupport.points->size(); ++i) {
            weightedSupport.points->push_back(wideWallSupport.points->points[i]);
            weightedSupport.weights.push_back(kSupportWideWallWeight);
        }
    }
    if (weightedSupport.points->size() < 20) {
        // 非墙面回退：歧义归属已在提取器内部排除；
        // 权重明显低于墙面点。
        auto fallbackSupport = ExtractWeightedBoundarySupportFromOSGB(
            sampled, ring, SupportExtractMode::NonWall,
            kdtree, ownership, currentRingId);
        fallbackCount = fallbackSupport.points->size();
        g_supportFilterFallbackUsedTotal += fallbackCount;
        for (std::size_t i = 0; i < fallbackSupport.points->size(); ++i) {
            weightedSupport.points->push_back(fallbackSupport.points->points[i]);
            weightedSupport.weights.push_back(kSupportFallbackWeight);
        }
    }
    if (weightedSupport.points->size() < 20) {
        // 仅作为 Ceres 的最后兜底：轮廓自身的
        // 加密边界，权重最低。
        auto densified = DensifyBoundary(ring, kDensifyStep);
        densifyCount = densified->size();
        for (const auto& p : *densified) {
            weightedSupport.points->push_back(p);
            weightedSupport.weights.push_back(kSupportDensifyWeight);
        }
    }
    // 支撑点过多时降采样：computeModelResolution 是 O(N²)，Ceres/拓扑修复也随点数变贵
    // Ceres 只用 XY。三维体素会保留同一墙面不同高度
    // 的重复点，无意中放大高墙的影响。
    auto downsampled = DownsampleSupport2DWithWeights(weightedSupport, kSupportVoxelLeaf);
    if (downsampled.points && !downsampled.points->empty()) weightedSupport = downsampled;
    ApplySupportDensityWeights(weightedSupport);
    auto support = weightedSupport.points;
    if (debugSupport) {
        debugSupport->assign(support->points.begin(), support->points.end());
    }
    auto t1 = std::chrono::steady_clock::now();
    g_supportTime += std::chrono::duration<double>(t1 - t0).count();

    std::cerr << "[SupportFilter] fid=" << sourceFid
              << " coarse_wall=" << filterStats.coarseWall
              << " effective_wall=" << filterStats.narrowBand + filterStats.extendedBand
              << " narrow_band=" << filterStats.narrowBand
              << " extended_band=" << filterStats.extendedBand
              << " wall_wide=" << wideWallCount
              << " fallback=" << fallbackCount
              << " densify=" << densifyCount
              << " ambiguous_rejected=" << filterStats.ambiguousRejected
              << " other_owner_rejected=" << filterStats.otherOwnerRejected
              << " normal_rejected=" << filterStats.normalRejected
              << " t_rejected=" << filterStats.tRejected
              << " direction_input=" << effectiveWallSupportRaw.points->size()
              << " ceres_input=" << support->size() << std::endl;

    // 小建筑不再在此走方向矩形快通道；
    // 它们走正常管线，由 outlineRegular 强制单一主方向
    // (kSmallBuildingSingleDirectionArea)，避免多方向
    // 规则化在小轮廓上产生斜边。
    // footprints.

    // ---- 规则化优化(计时) ----
    // AI 掩膜即初始轮廓，不是独立的 DLG 观测。
    outlineRegular regularizer(ring, support, weightedSupport.weights);
    regularizer.setSourceFeatureId(sourceFid);
    double wallDirection = 0.0;
    double wallPeakRatio = 0.0;
    std::size_t wallPairCount = 0;
    {
        // 方向直方图的强度门槛(kSupportDirectionMinPairs)
        // 按体素降采样后的墙面点校准。直接喂入原始
        // 有效集(30点/m²)会让配对数虚高约20倍，使一段
        // 短异常墙就能冒充"强"证据，把建筑
        // 旋转(实测：400-8000配对导致15-45°的"纠偏")。
        // 提示前先降采样以恢复尺度一致性。
        auto downsampledEffective =
            DownsampleSupport2DWithWeights(effectiveWallSupportRaw, kSupportVoxelLeaf);
        const auto directionInput =
            (downsampledEffective.points && !downsampledEffective.points->empty())
                ? downsampledEffective.points
                : effectiveWallSupportRaw.points;
        outlineRegular::estimateSupportDirection2D(
            directionInput, wallDirection, wallPeakRatio, wallPairCount);
    }
    regularizer.setSupportDirectionHint(wallDirection, wallPeakRatio, wallPairCount);
    auto t2 = std::chrono::steady_clock::now();
    regularizer.regular_Contour();
    auto t3 = std::chrono::steady_clock::now();
    g_optimizeTime += std::chrono::duration<double>(t3 - t2).count();

    if (debugBestHypothesis) {
        *debugBestHypothesis = regularizer.getBestEnergyHypothesis();
        if (debugBestHypothesis->size() < 3) *debugBestHypothesis = ring;
        RemoveClosingDuplicate(*debugBestHypothesis);
    }

    std::vector<pcl::PointXYZ> result;
    if (regularizer.final_points && regularizer.final_points->size() >= 3) {
        result.assign(regularizer.final_points->points.begin(), regularizer.final_points->points.end());
    } else {
        result = ring;
    }
    RemoveClosingDuplicate(result);
    return result.size() >= 3 ? result : ring;
}

// 作用：计算点列的轴对齐包围盒面积(m²)。
double BoundingBoxArea2D(const std::vector<pcl::PointXYZ>& ring)
{
    if (ring.empty()) return 0.0;
    double min_x = ring.front().x;
    double max_x = ring.front().x;
    double min_y = ring.front().y;
    double max_y = ring.front().y;
    for (const auto& p : ring) {
        min_x = std::min(min_x, static_cast<double>(p.x));
        max_x = std::max(max_x, static_cast<double>(p.x));
        min_y = std::min(min_y, static_cast<double>(p.y));
        max_y = std::max(max_y, static_cast<double>(p.y));
    }
    return std::max(0.0, max_x - min_x) * std::max(0.0, max_y - min_y);
}

// 作用：计算 2D 多边形面积(鞋带公式，绝对值)。
double PolygonArea2D(const std::vector<pcl::PointXYZ>& ring)
{
    if (ring.size() < 3) return 0.0;
    double area = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto& a = ring[i];
        const auto& b = ring[(i + 1) % ring.size()];
        area += static_cast<double>(a.x) * b.y - static_cast<double>(b.x) * a.y;
    }
    return std::abs(area) * 0.5;
}

// 作用：按指定方向生成点列的最小旋转外接矩形(4 顶点)。
std::vector<pcl::PointXYZ> OrientedBoundingRectangle(
    const std::vector<pcl::PointXYZ>& ring,
    double angle)
{
    std::vector<pcl::PointXYZ> result;
    if (ring.size() < 3) return result;
    const double ux = std::cos(angle);
    const double uy = std::sin(angle);
    const double vx = -uy;
    const double vy = ux;
    double minU = std::numeric_limits<double>::max();
    double maxU = -std::numeric_limits<double>::max();
    double minV = std::numeric_limits<double>::max();
    double maxV = -std::numeric_limits<double>::max();
    for (const auto& p : ring) {
        const double u = p.x * ux + p.y * uy;
        const double v = p.x * vx + p.y * vy;
        minU = std::min(minU, u);
        maxU = std::max(maxU, u);
        minV = std::min(minV, v);
        maxV = std::max(maxV, v);
    }
    const float z = ring.front().z;
    auto point = [&](double u, double v) {
        return pcl::PointXYZ(static_cast<float>(u * ux + v * vx),
            static_cast<float>(u * uy + v * vy), z);
    };
    result = { point(minU, minV), point(maxU, minV),
               point(maxU, maxV), point(minU, maxV) };
    return result;
}

// 作用：检查输出环是否通过面积/bbox 下限安全地板。
bool OutputRingPassesSizeFloor(const std::vector<pcl::PointXYZ>& ring)
{
    return ring.size() >= 3 &&
        PolygonArea2D(ring) >= kMinOutputPolygonArea &&
        BoundingBoxArea2D(ring) >= kMinOutputPolygonBBoxArea;
}

// 作用：判断环内及环边 3m 范围内是否有足够的采样点(模型覆盖证据)。
bool RingHasModelCoverage(
    const std::vector<pcl::PointXYZ>& ring,
    const MyCloudPtr& sampled,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree)
{
    if (ring.size() < 3 || !sampled || !sampled->cloud || sampled->cloud->empty() || !kdtree) {
        return false;
    }

    double minX = ring.front().x;
    double maxX = ring.front().x;
    double minY = ring.front().y;
    double maxY = ring.front().y;
    for (const auto& p : ring) {
        minX = std::min(minX, static_cast<double>(p.x));
        maxX = std::max(maxX, static_cast<double>(p.x));
        minY = std::min(minY, static_cast<double>(p.y));
        maxY = std::max(maxY, static_cast<double>(p.y));
    }

    const Eigen::Vector3f cloudMin = sampled->boundingBox.min();
    const Eigen::Vector3f cloudMax = sampled->boundingBox.max();
    if (maxX + kModelCoverageBuffer < cloudMin.x() ||
        minX - kModelCoverageBuffer > cloudMax.x() ||
        maxY + kModelCoverageBuffer < cloudMin.y() ||
        minY - kModelCoverageBuffer > cloudMax.y()) {
        return false;
    }

    const double centerX = (minX + maxX) * 0.5;
    const double centerY = (minY + maxY) * 0.5;
    const double halfWidth = (maxX - minX) * 0.5 + kModelCoverageBuffer;
    const double halfHeight = (maxY - minY) * 0.5 + kModelCoverageBuffer;
    const float radius = static_cast<float>(std::hypot(halfWidth, halfHeight) + 1e-3);

    pcl::PointXYZ center;
    center.x = static_cast<float>(centerX);
    center.y = static_cast<float>(centerY);
    center.z = 0.0f;
    std::vector<int> candidates;
    std::vector<float> distances;
    kdtree->radiusSearch(center, radius, candidates, distances);

    std::size_t evidenceCount = 0;
    for (int index : candidates) {
        if (index < 0 || static_cast<std::size_t>(index) >= sampled->cloud->size()) continue;
        const auto& source = sampled->cloud->points[static_cast<std::size_t>(index)];
        if (source.x < minX - kModelCoverageBuffer || source.x > maxX + kModelCoverageBuffer ||
            source.y < minY - kModelCoverageBuffer || source.y > maxY + kModelCoverageBuffer) {
            continue;
        }

        pcl::PointXYZ point;
        point.x = source.x;
        point.y = source.y;
        point.z = source.z;
        bool supported = PointInPolygon2D(point, ring);
        if (!supported) {
            for (std::size_t i = 0; i < ring.size(); ++i) {
                if (PointToSegmentDistance2D(point, ring[i], ring[(i + 1) % ring.size()]) <=
                    kModelCoverageBuffer) {
                    supported = true;
                    break;
                }
            }
        }
        if (supported && ++evidenceCount >= kMinModelEvidencePoints) return true;
    }
    return false;
}

SupportOwnershipContext BuildSupportOwnershipContext(
    OGRLayer* layer,
    const Eigen::Vector3d& metadataOffset,
    const MyCloudPtr& sampled,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree)
{
    SupportOwnershipContext context;
    if (!layer) return context;

    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (!geometry) {
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
        if (type == wkbPolygon) {
            auto ring = ExtractExteriorRing(geometry->toPolygon(), metadataOffset);
            if (BoundingBoxArea2D(ring) >= kMinPolygonBBoxArea &&
                RingHasModelCoverage(ring, sampled, kdtree)) {
                InitialRingRecord record;
                record.sourceFid = feature->GetFID();
                record.ringIndex = 0;
                record.ring = std::move(ring);
                context.addRecord(std::move(record));
            }
        } else if (type == wkbMultiPolygon) {
            auto* multi = geometry->toMultiPolygon();
            int partIndex = 0;
            for (auto&& part : *multi) {
                auto ring = ExtractExteriorRing(part->toPolygon(), metadataOffset);
                if (BoundingBoxArea2D(ring) >= kMinPolygonBBoxArea &&
                    RingHasModelCoverage(ring, sampled, kdtree)) {
                    InitialRingRecord record;
                    record.sourceFid = feature->GetFID();
                    record.ringIndex = partIndex;
                    record.ring = std::move(ring);
                    context.addRecord(std::move(record));
                }
                ++partIndex;
            }
        }

        OGRFeature::DestroyFeature(feature);
    }
    layer->ResetReading();
    return context;
}

struct InitialOutlineMergeStats {
    long long adjacentCandidates = 0;
    long long keptByModelEvidence = 0;
    long long keptByShape = 0;
    long long keptNarrowSeam = 0;     // 共享边 < kMergeSeamLength(窄缝连接，不合并)
    long long mergedPairs = 0;
    long long mergedByBigPair = 0;    // 大建筑兜底合并(跳过证据间的计数，便于调 kLargeBuildingMergeArea)
    long long mergedByShape = 0;      // shared-boundary geometry strongly suggests an over-segmentation seam
    long long removedFeatures = 0;
    // 形状判定(keptByShape)的子原因拆分(诊断用)
    long long shapeUnionFailed = 0;       // A∪B 返回空
    long long shapeNotSinglePolygon = 0;  // 并集不是单个 wkbPolygon(常见于栅格小缺口→MultiPolygon)
    long long shapePerimeterRatio = 0;    // 共享边占周长比 < kMinSharedPerimeterRatio
    long long shapeCompactness = 0;       // 绱у噾搴︽敹鐩?> kMaxSplitCompactnessGainForMerge
    double maxCompactnessGain = 0.0;      // 观察到的最大紧凑度增益(判断阈值要放宽多少)
    long long mergedByParent = 0;
};

struct OutlineFeatureRecord {
    GIntBig fid = OGRNullFID;
    std::unique_ptr<OGRGeometry> geometry;
    OGREnvelope envelope = {};
    double area = 0.0;
    double perimeter = 0.0;
    int parent = -1;   // 所属原始连通分量(<=0 表示未知，不参与同源合并)
};

struct OutputOverlapRepairStats {
    long long candidatePairs = 0;
    long long overlapPairs = 0;
    long long resolvedPairs = 0;
    long long unresolvedPairs = 0;
    long long shiftedFeatures = 0;
    long long optimizedGroups = 0;
    double maxOverlapArea = 0.0;
    double maxShiftDistance = 0.0;
};

struct SeamSegment {
    pcl::PointXYZ a;
    pcl::PointXYZ b;
};

struct SeamEvidence {
    double wallRatio = 0.0;
    double heightRatio = 0.0;
    int modelSampleCount = 0;
    int wallCount = 0;
    int heightPairCount = 0;
    int heightDifferenceCount = 0;
};

struct NarrowNeckSplitStats {
    long long inspectedFeatures = 0;
    long long splitFeatures = 0;
    long long cuts = 0;
    long long createdParts = 0;
    long long candidateOnly = 0;
    long long rejectedInvalidRing = 0;
    long long rejectedSmallPartArea = 0;
    long long rejectedInvalidPolygon = 0;
    long long debugRejectLogs = 0;
};

constexpr int kNarrowNeckDebugRejectLogLimit = 24;

void CopyFieldValues(OGRFeature* src, OGRFeature* dst);

// 作用：计算 OGR 几何的长度(线)。
double GeometryLength(const OGRGeometry* geometry)
{
    if (!geometry || geometry->IsEmpty()) return 0.0;
    const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    // 线性几何直接求长; 集合递归累加; 面走 Boundary 后再求;
    // 点等非曲线类型不得调用 OGR_G_Length(会触发 non-curve 警告)
    if (type == wkbLineString || type == wkbLinearRing) {
        return OGR_G_Length(OGRGeometry::ToHandle(
            const_cast<OGRGeometry*>(geometry)));
    }
    if (type == wkbMultiLineString || type == wkbGeometryCollection) {
        const auto* collection = geometry->toGeometryCollection();
        double length = 0.0;
        for (int i = 0; collection && i < collection->getNumGeometries(); ++i) {
            length += GeometryLength(collection->getGeometryRef(i));
        }
        return length;
    }
    if (type == wkbPolygon || type == wkbTriangle || type == wkbMultiPolygon) {
        // 周长: 先取 Boundary(线几何) 再递归
        std::unique_ptr<OGRGeometry> boundary(geometry->Boundary());
        return boundary ? GeometryLength(boundary.get()) : 0.0;
    }
    return 0.0;
}

// 作用：计算 OGR 几何的面积(面)。
double GeometryArea(const OGRGeometry* geometry)
{
    if (!geometry || geometry->IsEmpty()) return 0.0;
    const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    if (type == wkbPolygon || type == wkbTriangle) {
        return OGR_G_Area(OGRGeometry::ToHandle(
            const_cast<OGRGeometry*>(geometry)));
    }
    if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
        const auto* collection = geometry->toGeometryCollection();
        double area = 0.0;
        for (int i = 0; collection && i < collection->getNumGeometries(); ++i) {
            area += GeometryArea(collection->getGeometryRef(i));
        }
        return area;
    }
    // Intersection/union may also contain line or point components. They do
    // not contribute to area and must not call OGR_G_Area (which emits a
    // misleading non-surface warning).
    return 0.0;
}

// 作用：计算 OGR 几何的周长(边界长度)。
double GeometryPerimeter(OGRGeometry* geometry)
{
    if (!geometry) return 0.0;
    std::unique_ptr<OGRGeometry> boundary(geometry->Boundary());
    return GeometryLength(boundary.get());
}

struct NarrowNeckCandidate {
    double ax = 0.0;
    double ay = 0.0;
    double bx = 0.0;
    double by = 0.0;
    double width = std::numeric_limits<double>::max();
    double boundarySeparation = 0.0;
    std::size_t vertexA = 0;
    std::size_t vertexB = 0;
};

struct Point2D64 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// 作用：两个 Point2D64 的水平距离。
double Distance2D64(const Point2D64& a, const Point2D64& b)
{
    return std::hypot(a.x - b.x, a.y - b.y);
}

// 作用：取出 OGR 多边形外环并转为 Point2D64 列(去闭合点)。
std::vector<Point2D64> ExtractExteriorRing64(OGRPolygon* polygon)
{
    std::vector<Point2D64> ring;
    if (!polygon || !polygon->getExteriorRing()) return ring;
    OGRLinearRing* ogrRing = polygon->getExteriorRing();
    for (int i = 0; i < ogrRing->getNumPoints(); ++i) {
        ring.push_back({ogrRing->getX(i), ogrRing->getY(i), ogrRing->getZ(i)});
    }
    if (ring.size() >= 2 && Distance2D64(ring.front(), ring.back()) < 1e-9) {
        ring.pop_back();
    }
    return ring;
}

// 作用：Point2D64 环的有向面积(带符号)。
double PolygonArea2D64(const std::vector<Point2D64>& ring)
{
    if (ring.size() < 3) return 0.0;
    double area = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto& a = ring[i];
        const auto& b = ring[(i + 1) % ring.size()];
        area += a.x * b.y - b.x * a.y;
    }
    return std::abs(area) * 0.5;
}

// 作用：射线法判断点是否在 Point2D64 多边形内。
bool PointInPolygon2D64(const Point2D64& p, const std::vector<Point2D64>& poly)
{
    bool inside = false;
    const std::size_t n = poly.size();
    if (n < 3) return false;
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const auto& pi = poly[i];
        const auto& pj = poly[j];
        const bool intersect = ((pi.y > p.y) != (pj.y > p.y)) &&
            (p.x < (pj.x - pi.x) * (p.y - pi.y) / ((pj.y - pi.y) + 1e-12) + pi.x);
        if (intersect) inside = !inside;
    }
    return inside;
}

// 作用：三点方向叉积。
double Orient2D64(const Point2D64& a, const Point2D64& b, const Point2D64& c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// 作用：判断点是否在线段上(含共线与范围检查)。
bool OnSegment2D64(const Point2D64& a, const Point2D64& b, const Point2D64& p)
{
    constexpr double eps = 1e-8;
    if (std::abs(Orient2D64(a, b, p)) > eps) return false;
    return p.x >= std::min(a.x, b.x) - eps && p.x <= std::max(a.x, b.x) + eps &&
           p.y >= std::min(a.y, b.y) - eps && p.y <= std::max(a.y, b.y) + eps;
}

// 作用：判断两条线段是否真正相交(不含端点触碰)。
bool SegmentsIntersect2D64(
    const Point2D64& a,
    const Point2D64& b,
    const Point2D64& c,
    const Point2D64& d)
{
    constexpr double eps = 1e-8;
    const double o1 = Orient2D64(a, b, c);
    const double o2 = Orient2D64(a, b, d);
    const double o3 = Orient2D64(c, d, a);
    const double o4 = Orient2D64(c, d, b);

    if (((o1 > eps && o2 < -eps) || (o1 < -eps && o2 > eps)) &&
        ((o3 > eps && o4 < -eps) || (o3 < -eps && o4 > eps))) {
        return true;
    }
    if (std::abs(o1) <= eps && OnSegment2D64(a, b, c)) return true;
    if (std::abs(o2) <= eps && OnSegment2D64(a, b, d)) return true;
    if (std::abs(o3) <= eps && OnSegment2D64(c, d, a)) return true;
    if (std::abs(o4) <= eps && OnSegment2D64(c, d, b)) return true;
    return false;
}

// 作用：判断边是否与顶点重合接触。
bool EdgeTouchesVertex(std::size_t edgeIndex, std::size_t vertexIndex, std::size_t vertexCount)
{
    return edgeIndex == vertexIndex || ((edgeIndex + 1) % vertexCount) == vertexIndex;
}

// 作用：判断窄颈切线(除端点外)是否在多边形内部。
bool NeckCutLineInsidePolygon2D64(
    const std::vector<Point2D64>& ring,
    std::size_t vertexA,
    std::size_t vertexB)
{
    const std::size_t n = ring.size();
    if (n < 4 || vertexA >= n || vertexB >= n || vertexA == vertexB) return false;
    const Point2D64& a = ring[vertexA];
    const Point2D64& b = ring[vertexB];

    // 内部采样点拒绝穿过外部空间的弦。
    constexpr int sampleCount = 7;
    for (int s = 1; s < sampleCount; ++s) {
        const double t = static_cast<double>(s) / sampleCount;
        Point2D64 p;
        p.x = a.x * (1.0 - t) + b.x * t;
        p.y = a.y * (1.0 - t) + b.y * t;
        p.z = a.z * (1.0 - t) + b.z * t;
        if (!PointInPolygon2D64(p, ring)) return false;
    }

    // A valid internal chord may only touch the original boundary at its two
    // endpoints. Any other boundary intersection means the line leaves the
    // polygon or crosses another lobe.
    for (std::size_t edge = 0; edge < n; ++edge) {
        if (EdgeTouchesVertex(edge, vertexA, n) ||
            EdgeTouchesVertex(edge, vertexB, n)) {
            continue;
        }
        if (SegmentsIntersect2D64(a, b, ring[edge], ring[(edge + 1) % n])) {
            return false;
        }
    }
    return true;
}

// 作用：把 Point2D64 环组装成 OGRPolygon。
OGRPolygon* MakePolygon64(const std::vector<Point2D64>& ring)
{
    if (ring.size() < 3) return nullptr;
    auto* polygon = new OGRPolygon();
    OGRLinearRing ogrRing;
    for (const auto& p : ring) {
        ogrRing.addPoint(p.x, p.y, p.z);
    }
    ogrRing.addPoint(ring.front().x, ring.front().y, ring.front().z);
    polygon->addRing(&ogrRing);
    return polygon;
}

struct ClosestSegmentPoints2D {
    double ax = 0.0;
    double ay = 0.0;
    double bx = 0.0;
    double by = 0.0;
    double ta = 0.0;
    double tb = 0.0;
    double distance = std::numeric_limits<double>::max();
};

// 作用：钳制到 [0,1]。
double Clamp01(double value)
{
    return std::max(0.0, std::min(1.0, value));
}

ClosestSegmentPoints2D ClosestPointsOnSegments2D(
    const pcl::PointXYZ& a0,
    const pcl::PointXYZ& a1,
    const pcl::PointXYZ& b0,
    const pcl::PointXYZ& b1)
{
    const double ux = a1.x - a0.x;
    const double uy = a1.y - a0.y;
    const double vx = b1.x - b0.x;
    const double vy = b1.y - b0.y;
    const double wx = a0.x - b0.x;
    const double aa = ux * ux + uy * uy;
    const double bb = ux * vx + uy * vy;
    const double cc = vx * vx + vy * vy;
    const double dd = ux * wx;
    const double ee = vx * wx;
    const double denom = aa * cc - bb * bb;

    double s = 0.0;
    double t = 0.0;
    if (aa <= 1e-12 && cc <= 1e-12) {
        s = t = 0.0;
    } else if (aa <= 1e-12) {
        s = 0.0;
        t = Clamp01(ee / cc);
    } else if (cc <= 1e-12) {
        t = 0.0;
        s = Clamp01(-dd / aa);
    } else if (std::abs(denom) > 1e-12) {
        s = Clamp01((bb * ee - cc * dd) / denom);
        t = Clamp01((aa * ee - bb * dd) / denom);
    } else {
        s = 0.0;
        t = Clamp01(ee / cc);
    }

    // 固定一个参数后重新钳制一次。对短
    // raster-outline segments and avoids pulling in a larger geometry library.
    const double px = a0.x + s * ux;
    const double py = a0.y + s * uy;
    t = cc > 1e-12 ? Clamp01(((px - b0.x) * vx + (py - b0.y) * vy) / cc) : 0.0;
    const double qx = b0.x + t * vx;
    const double qy = b0.y + t * vy;
    s = aa > 1e-12 ? Clamp01(((qx - a0.x) * ux + (qy - a0.y) * uy) / aa) : 0.0;

    ClosestSegmentPoints2D result;
    result.ta = s;
    result.tb = t;
    result.ax = a0.x + s * ux;
    result.ay = a0.y + s * uy;
    result.bx = b0.x + t * vx;
    result.by = b0.y + t * vy;
    result.distance = std::hypot(result.ax - result.bx, result.ay - result.by);
    return result;
}

// 作用：判断切线是否大部分位于多边形内。
bool CutLineMostlyInsidePolygon(
    const NarrowNeckCandidate& candidate,
    const std::vector<pcl::PointXYZ>& ring)
{
    for (int i = 1; i <= 5; ++i) {
        const double t = static_cast<double>(i) / 6.0;
        pcl::PointXYZ p;
        p.x = static_cast<float>(candidate.ax * (1.0 - t) + candidate.bx * t);
        p.y = static_cast<float>(candidate.ay * (1.0 - t) + candidate.by * t);
        p.z = 0.0f;
        if (!PointInPolygon2D(p, ring)) return false;
    }
    return true;
}

bool FindNarrowNeckCandidates(
    const std::vector<Point2D64>& ring,
    std::vector<NarrowNeckCandidate>& candidates);

// 作用：从候选中选最优窄颈(宽度最小、边界弧最长)。
bool FindBestNarrowNeckCandidate(
    const std::vector<Point2D64>& ring,
    NarrowNeckCandidate& best)
{
    std::vector<NarrowNeckCandidate> candidates;
    if (!FindNarrowNeckCandidates(ring, candidates)) return false;
    best = candidates.front();
    return true;
}

// 作用：枚举所有满足宽度/弧长/面积门槛的窄颈候选顶点对。
bool FindNarrowNeckCandidates(
    const std::vector<Point2D64>& ring,
    std::vector<NarrowNeckCandidate>& candidates)
{
    candidates.clear();
    if (ring.size() < 4) return false;

    const std::size_t n = ring.size();
    std::vector<double> edgeLengths(n, 0.0);
    std::vector<double> prefix(n + 1, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        edgeLengths[i] = Distance2D64(ring[i], ring[(i + 1) % n]);
        prefix[i + 1] = prefix[i] + edgeLengths[i];
    }
    const double perimeter = prefix.back();
    if (perimeter <= 2.0 * kNarrowNeckMinBoundarySeparation) return false;

    // A narrow neck is defined by two non-adjacent boundary vertices.
    // 两点之间的两段弧成为两个新的外环。
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const std::size_t forwardEdges = j - i;
            const std::size_t backwardEdges = n - forwardEdges;
            if (forwardEdges < 2 || backwardEdges < 2) continue;

            const double arcA = prefix[j] - prefix[i];
            const double arcB = perimeter - arcA;
            const double boundarySeparation = std::min(arcA, arcB);
            if (arcA < kNarrowNeckMinBoundarySeparation ||
                arcB < kNarrowNeckMinBoundarySeparation) {
                continue;
            }

            const double width = Distance2D64(ring[i], ring[j]);
            if (width <= 1e-6 || width > kNarrowNeckMaxWidth) continue;
            if (!NeckCutLineInsidePolygon2D64(ring, i, j)) continue;

            std::vector<Point2D64> ringA;
            std::vector<Point2D64> ringB;
            for (std::size_t index = i;; index = (index + 1) % n) {
                ringA.push_back(ring[index]);
                if (index == j) break;
            }
            for (std::size_t index = j;; index = (index + 1) % n) {
                ringB.push_back(ring[index]);
                if (index == i) break;
            }
            std::unique_ptr<OGRPolygon> polygonA(MakePolygon64(ringA));
            std::unique_ptr<OGRPolygon> polygonB(MakePolygon64(ringB));
            if (!polygonA || !polygonB) continue;
            const double areaA = GeometryArea(polygonA.get());
            const double areaB = GeometryArea(polygonB.get());
            if (areaA < kNarrowNeckMinPartArea || areaB < kNarrowNeckMinPartArea) {
                continue;
            }

            NarrowNeckCandidate candidate;
            candidate.ax = ring[i].x;
            candidate.ay = ring[i].y;
            candidate.bx = ring[j].x;
            candidate.by = ring[j].y;
            candidate.width = width;
            candidate.boundarySeparation = boundarySeparation;
            candidate.vertexA = i;
            candidate.vertexB = j;
            candidates.push_back(candidate);
        }
    }

    std::sort(candidates.begin(), candidates.end(),
        [](const NarrowNeckCandidate& lhs, const NarrowNeckCandidate& rhs) {
            if (std::abs(lhs.width - rhs.width) > 1e-6) {
                return lhs.width < rhs.width;
            }
            return lhs.boundarySeparation > rhs.boundarySeparation;
        });
    return !candidates.empty();
}

// 作用：把(多)多边形几何拆成面积达标的单多边形部件列表。
void CollectPolygonParts(
    OGRGeometry* geometry,
    std::vector<std::unique_ptr<OGRGeometry>>& parts)
{
    if (!geometry || geometry->IsEmpty()) return;
    const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    if (type == wkbPolygon) {
        if (GeometryArea(geometry) >= kNarrowNeckMinPartArea) {
            parts.emplace_back(geometry->clone());
        }
    } else if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
        auto* collection = geometry->toGeometryCollection();
        for (int i = 0; collection && i < collection->getNumGeometries(); ++i) {
            CollectPolygonParts(collection->getGeometryRef(i), parts);
        }
    }
}

// 作用：用一个窄颈候选切分单个多边形，输出各部件。
bool SplitOnePolygonAtNarrowNeck(
OGRPolygon* polygon,
    std::vector<std::unique_ptr<OGRGeometry>>& splitParts,
    bool* foundCandidate = nullptr,
    NarrowNeckSplitStats* stats = nullptr,
    GIntBig fid = -1)
{
    splitParts.clear();
    if (!polygon || GeometryArea(polygon) < kNarrowNeckMinPartArea * 2.0) return false;

    auto ring = ExtractExteriorRing64(polygon);
    std::vector<NarrowNeckCandidate> candidates;
    if (!FindNarrowNeckCandidates(ring, candidates)) return false;
    if (foundCandidate) *foundCandidate = true;

    const std::size_t n = ring.size();
    for (const auto& candidate : candidates) {
        splitParts.clear();
        const std::size_t a = candidate.vertexA;
        const std::size_t b = candidate.vertexB;
        if (a >= n || b >= n || a == b) continue;

        std::vector<Point2D64> ringA;
        std::vector<Point2D64> ringB;
        for (std::size_t index = a;; index = (index + 1) % n) {
            ringA.push_back(ring[index]);
            if (index == b) break;
        }
        for (std::size_t index = b;; index = (index + 1) % n) {
            ringB.push_back(ring[index]);
            if (index == a) break;
        }

        if (ringA.size() < 3 || ringB.size() < 3) {
            if (stats) {
                ++stats->rejectedInvalidRing;
                if (stats->debugRejectLogs < kNarrowNeckDebugRejectLogLimit) {
                    ++stats->debugRejectLogs;
                    std::cout << "[Mask neck split reject] fid=" << fid
                              << " reason=invalid_ring"
                              << " width=" << candidate.width
                              << " sep=" << candidate.boundarySeparation
                              << " vA=" << candidate.vertexA
                              << " vB=" << candidate.vertexB << std::endl;
                }
            }
            continue;
        }

        std::unique_ptr<OGRPolygon> polygonA(
            MakePolygon64(ringA));
        std::unique_ptr<OGRPolygon> polygonB(
            MakePolygon64(ringB));
        if (!polygonA || !polygonB) {
            if (stats) {
                ++stats->rejectedInvalidPolygon;
                if (stats->debugRejectLogs < kNarrowNeckDebugRejectLogLimit) {
                    ++stats->debugRejectLogs;
                    std::cout << "[Mask neck split reject] fid=" << fid
                              << " reason=invalid_polygon"
                              << " width=" << candidate.width
                              << " sep=" << candidate.boundarySeparation
                              << " vA=" << candidate.vertexA
                              << " vB=" << candidate.vertexB << std::endl;
                }
            }
            continue;
        }
        const double areaA = GeometryArea(polygonA.get());
        const double areaB = GeometryArea(polygonB.get());
        if (areaA < kNarrowNeckMinPartArea ||
            areaB < kNarrowNeckMinPartArea) {
            if (stats) {
                ++stats->rejectedSmallPartArea;
                if (stats->debugRejectLogs < kNarrowNeckDebugRejectLogLimit) {
                    ++stats->debugRejectLogs;
                    std::cout << "[Mask neck split reject] fid=" << fid
                              << " reason=small_part_area"
                              << " width=" << candidate.width
                              << " sep=" << candidate.boundarySeparation
                              << " areaA=" << areaA
                              << " areaB=" << areaB
                              << " vA=" << candidate.vertexA
                              << " vB=" << candidate.vertexB << std::endl;
                }
            }
            continue;
        }

        splitParts.emplace_back(polygonA.release());
        splitParts.emplace_back(polygonB.release());
        return true;
    }

    if (stats) {
        ++stats->rejectedInvalidPolygon;
        if (stats->debugRejectLogs < kNarrowNeckDebugRejectLogLimit) {
            ++stats->debugRejectLogs;
            std::cout << "[Mask neck split reject] fid=" << fid
                      << " reason=all_candidates_failed"
                      << " candidates=" << candidates.size() << std::endl;
        }
    }
    splitParts.clear();
    return false;
}

void SplitGeometryAtNarrowNecks(
    OGRGeometry* geometry,
    std::vector<std::unique_ptr<OGRGeometry>>& outputParts,
    long long& cutCount,
    long long& candidateOnlyCount,
    NarrowNeckSplitStats* stats = nullptr,
    GIntBig fid = -1)
{
    if (!geometry || geometry->IsEmpty()) return;
    const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
        auto* collection = geometry->toGeometryCollection();
        for (int i = 0; collection && i < collection->getNumGeometries(); ++i) {
            SplitGeometryAtNarrowNecks(
                collection->getGeometryRef(i), outputParts, cutCount, candidateOnlyCount, stats, fid);
        }
        return;
    }
    if (type != wkbPolygon) {
        outputParts.emplace_back(geometry->clone());
        return;
    }

    std::vector<std::unique_ptr<OGRGeometry>> pending;
    pending.emplace_back(geometry->clone());
    int localCuts = 0;
    for (int pass = 0;
         pass < kNarrowNeckMaxCutsPerFeature &&
         localCuts < kNarrowNeckMaxCutsPerFeature &&
         !pending.empty();
         ++pass) {
        bool splitThisPass = false;
        std::vector<std::unique_ptr<OGRGeometry>> next;
        for (auto& item : pending) {
            auto* polygon = item->toPolygon();
            std::vector<std::unique_ptr<OGRGeometry>> splitParts;
            bool foundCandidate = false;
            if (localCuts < kNarrowNeckMaxCutsPerFeature &&
                SplitOnePolygonAtNarrowNeck(polygon, splitParts, &foundCandidate, stats, fid)) {
                ++cutCount;
                ++localCuts;
                splitThisPass = true;
                for (auto& part : splitParts) next.push_back(std::move(part));
            } else {
                if (foundCandidate) ++candidateOnlyCount;
                next.push_back(std::move(item));
            }
        }
        pending = std::move(next);
        if (!splitThisPass) break;
    }
    for (auto& item : pending) outputParts.push_back(std::move(item));
}

// 作用：初始轮廓窄颈拆分阶段入口：就地修改 Shapefile。
bool SplitInitialOutlinesAtNarrowNecks(
    const std::string& shpPath,
    NarrowNeckSplitStats& stats)
{
    stats = {};
    GDALDataset* dataset = static_cast<GDALDataset*>(
        GDALOpenEx(shpPath.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
                   nullptr, nullptr, nullptr));
    if (!dataset) return false;
    OGRLayer* layer = dataset->GetLayer(0);
    if (!layer) {
        GDALClose(dataset);
        return false;
    }

    std::vector<GIntBig> fids;
    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        fids.push_back(feature->GetFID());
        OGRFeature::DestroyFeature(feature);
    }

    OGRFeatureDefn* defn = layer->GetLayerDefn();
    for (GIntBig fid : fids) {
        OGRFeature* feature = layer->GetFeature(fid);
        if (!feature) continue;
        ++stats.inspectedFeatures;
        OGRGeometry* geometry = feature->GetGeometryRef();
        std::vector<std::unique_ptr<OGRGeometry>> parts;
        long long cutsBefore = stats.cuts;
        SplitGeometryAtNarrowNecks(geometry, parts, stats.cuts, stats.candidateOnly, &stats, fid);
        if (parts.size() <= 1 || stats.cuts == cutsBefore) {
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        if (layer->DeleteFeature(fid) != OGRERR_NONE) {
            OGRFeature::DestroyFeature(feature);
            continue;
        }
        ++stats.splitFeatures;
        stats.createdParts += static_cast<long long>(parts.size());
        for (auto& part : parts) {
            OGRFeature* out = OGRFeature::CreateFeature(defn);
            CopyFieldValues(feature, out);
            out->SetGeometry(part.get());
            layer->CreateFeature(out);
            OGRFeature::DestroyFeature(out);
        }
        OGRFeature::DestroyFeature(feature);
    }

    layer->SyncToDisk();
    GDALClose(dataset);
    return true;
}

struct BoundarySegment2D {
    double ax = 0.0;
    double ay = 0.0;
    double bx = 0.0;
    double by = 0.0;
    double length = 0.0;
};

// 作用：把外环每条边加入边界段列表。
void AddExteriorRingSegments(const OGRLinearRing* ring, std::vector<BoundarySegment2D>& segments)
{
    if (!ring || ring->getNumPoints() < 2) return;
    for (int i = 1; i < ring->getNumPoints(); ++i) {
        BoundarySegment2D segment;
        segment.ax = ring->getX(i - 1);
        segment.ay = ring->getY(i - 1);
        segment.bx = ring->getX(i);
        segment.by = ring->getY(i);
        segment.length = std::hypot(segment.bx - segment.ax, segment.by - segment.ay);
        if (segment.length > 1e-6) segments.push_back(segment);
    }
}

// 作用：收集几何的全部外环边界段(支持多部件)。
void CollectExteriorBoundarySegments(const OGRGeometry* geometry, std::vector<BoundarySegment2D>& segments)
{
    if (!geometry) return;
    const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    if (type == wkbPolygon) {
        const auto* polygon = geometry->toPolygon();
        if (polygon) AddExteriorRingSegments(polygon->getExteriorRing(), segments);
    } else if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
        const auto* collection = geometry->toGeometryCollection();
        for (int i = 0; collection && i < collection->getNumGeometries(); ++i) {
            CollectExteriorBoundarySegments(collection->getGeometryRef(i), segments);
        }
    }
}

// 作用：点到直线(用边界段表示)的距离。
double PointLineDistance2D(double x, double y, const BoundarySegment2D& segment)
{
    const double dx = segment.bx - segment.ax;
    const double dy = segment.by - segment.ay;
    return std::abs((x - segment.ax) * dy - (y - segment.ay) * dx) /
        std::max(segment.length, 1e-9);
}

// 作用：计算两边界段在平行条件下的投影重叠长度。
double SegmentProjectionOverlap2D(
    const BoundarySegment2D& a,
    const BoundarySegment2D& b,
    double tolerance)
{
    const double aux = (a.bx - a.ax) / a.length;
    const double auy = (a.by - a.ay) / a.length;
    const double bux = (b.bx - b.ax) / b.length;
    const double buy = (b.by - b.ay) / b.length;
    const double angleCos = std::abs(aux * bux + auy * buy);
    const double minCos = std::cos(kRobustSharedAngleToleranceDeg * M_PI / 180.0);
    if (angleCos < minCos) return 0.0;

    const double distB0 = PointLineDistance2D(b.ax, b.ay, a);
    const double distB1 = PointLineDistance2D(b.bx, b.by, a);
    const double midAx = 0.5 * (a.ax + a.bx);
    const double midAy = 0.5 * (a.ay + a.by);
    const double distAmid = PointLineDistance2D(midAx, midAy, b);
    if (std::min(distB0, distB1) > tolerance && distAmid > tolerance) {
        return 0.0;
    }

    const double b0 = (b.ax - a.ax) * aux + (b.ay - a.ay) * auy;
    const double b1 = (b.bx - a.ax) * aux + (b.by - a.ay) * auy;
    const double lo = std::max(0.0, std::min(b0, b1));
    const double hi = std::min(a.length, std::max(b0, b1));
    return std::max(0.0, hi - lo);
}

// 作用：鲁棒估计两轮廓的共享边界长度(精确值+容差投影)。
double RobustSharedBoundaryLength(
    const OGRGeometry* a,
    const OGRGeometry* b,
    double exactSharedLength)
{
    std::vector<BoundarySegment2D> segmentsA;
    std::vector<BoundarySegment2D> segmentsB;
    CollectExteriorBoundarySegments(a, segmentsA);
    CollectExteriorBoundarySegments(b, segmentsB);
    if (segmentsA.empty() || segmentsB.empty()) return exactSharedLength;

    double robustLength = 0.0;
    for (const auto& sa : segmentsA) {
        for (const auto& sb : segmentsB) {
            robustLength += SegmentProjectionOverlap2D(
                sa, sb, kRobustSharedBoundaryTolerance);
        }
    }
    return std::max(exactSharedLength, robustLength);
}

// 作用：清理合并后的几何(去自交、平滑细缝、保持拓扑)。
std::unique_ptr<OGRGeometry> CleanMergedGeometry(const OGRGeometry* geometry)
{
    if (!geometry) return nullptr;
    std::unique_ptr<OGRGeometry> cleaned(geometry->clone());
    if (!cleaned) return nullptr;

    // 先去除 Union/Buffer(0) 产生的自交伪影。
    if (!cleaned->IsValid()) {
        std::unique_ptr<OGRGeometry> fixed(cleaned->Buffer(0.0));
        if (fixed) cleaned = std::move(fixed);
    }

    // 再在保持拓扑的前提下平滑细小栅格缝。
    std::unique_ptr<OGRGeometry> simplified(
        cleaned->SimplifyPreserveTopology(kInitialMergeSimplifyTolerance));
    if (simplified && !simplified->IsEmpty()) {
        cleaned = std::move(simplified);
    }

    if (!cleaned->IsValid()) {
        std::unique_ptr<OGRGeometry> fixed(cleaned->Buffer(0.0));
        if (fixed) cleaned = std::move(fixed);
    }
    return cleaned;
}

// 作用：按形状判断两个过分割轮廓是否该合并(缺口互补等)。
bool ShouldMergeByFootprintShape(
    const OutlineFeatureRecord& a,
    const OutlineFeatureRecord& b,
    double sharedLength)
{
    if (!a.geometry || !b.geometry) return false;
    const double smallerPerimeter = std::max(1e-9, std::min(a.perimeter, b.perimeter));
    const double sharedRatio = sharedLength / smallerPerimeter;
    if (sharedLength < kShapeMergeMinSharedLength &&
        sharedRatio < kShapeMergeMinSharedPerimeterRatio) {
        return false;
    }

    std::unique_ptr<OGRGeometry> merged(a.geometry->Union(b.geometry.get()));
    merged = CleanMergedGeometry(merged.get());
    if (!merged || merged->IsEmpty() ||
        wkbFlatten(merged->getGeometryType()) != wkbPolygon) {
        return false;
    }

    const double mergedArea = GeometryArea(merged.get());
    const double areaSum = a.area + b.area;
    if (mergedArea <= 0.0 || mergedArea > areaSum * 1.02) return false;

    const double mergedPerimeter = GeometryPerimeter(merged.get());
    const double expectedPerimeter = std::max(0.0, a.perimeter + b.perimeter - 2.0 * sharedLength);
    if (expectedPerimeter > 1e-6 &&
        std::abs(mergedPerimeter - expectedPerimeter) / expectedPerimeter > 0.25) {
        return false;
    }

    return true;
}

// 作用：紧凑度 = 4πA/P²(圆为 1)。
double Compactness(double area, double perimeter)
{
    return area > 0.0 && perimeter > 0.0
        ? 4.0 * M_PI * area / (perimeter * perimeter) : 0.0;
}

// 作用：判断两几何是否存在有意义的重叠(超面积/比例阈值)。
bool HasMeaningfulGeometryOverlap(
    const OGRGeometry* a,
    const OGRGeometry* b,
    double& overlapArea)
{
    overlapArea = 0.0;
    if (!a || !b) return false;
    if (!a->Intersects(b)) return false;

    std::unique_ptr<OGRGeometry> overlap(a->Intersection(b));
    overlapArea = GeometryArea(overlap.get());
    if (overlapArea < kOutputOverlapMinArea) return false;

    const double areaA = GeometryArea(a);
    const double areaB = GeometryArea(b);
    const double smaller = std::max(1e-9, std::min(areaA, areaB));
    return overlapArea / smaller >= kOutputOverlapMinRatio;
}

// 作用：计算几何质心(面优先，失败退回包围盒中心)。
bool GeometryCentroid2D(const OGRGeometry* geometry, double& x, double& y)
{
    x = y = 0.0;
    if (!geometry || geometry->IsEmpty()) return false;
    OGRPoint centroid;
    if (const auto* surface = dynamic_cast<const OGRSurface*>(geometry)) {
        if (const_cast<OGRSurface*>(surface)->Centroid(&centroid) == OGRERR_NONE) {
            x = centroid.getX();
            y = centroid.getY();
            return true;
        }
    }
    OGREnvelope env;
    geometry->getEnvelope(&env);
    x = 0.5 * (env.MinX + env.MaxX);
    y = 0.5 * (env.MinY + env.MaxY);
    return true;
}

// 作用：把几何整体平移 (dx, dy)。
std::unique_ptr<OGRGeometry> TranslateGeometry2D(
    const OGRGeometry* geometry,
    double dx,
    double dy)
{
    if (!geometry) return nullptr;
    std::unique_ptr<OGRGeometry> moved(geometry->clone());
    if (!moved) return nullptr;

    std::function<void(OGRGeometry*)> translate = [&](OGRGeometry* g) {
        if (!g) return;
        const OGRwkbGeometryType type = wkbFlatten(g->getGeometryType());
        if (type == wkbPoint) {
            auto* p = g->toPoint();
            p->setX(p->getX() + dx);
            p->setY(p->getY() + dy);
        } else if (type == wkbLineString || type == wkbLinearRing) {
            auto* line = g->toLineString();
            for (int i = 0; i < line->getNumPoints(); ++i) {
                line->setPoint(i, line->getX(i) + dx, line->getY(i) + dy, line->getZ(i));
            }
        } else if (type == wkbPolygon) {
            auto* polygon = g->toPolygon();
            if (polygon->getExteriorRing()) translate(polygon->getExteriorRing());
            for (int i = 0; i < polygon->getNumInteriorRings(); ++i) {
                translate(polygon->getInteriorRing(i));
            }
        } else if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
            auto* collection = g->toGeometryCollection();
            for (int i = 0; i < collection->getNumGeometries(); ++i) {
                translate(collection->getGeometryRef(i));
            }
        }
    };
    translate(moved.get());
    return moved;
}

// 作用：尝试沿质心连线方向平移"动方"使其脱离"定方"(带新冲突检查)。
bool TryResolvePairByTranslation(
    const OGRGeometry* fixed,
    const OGRGeometry* moving,
    double overlapArea,
    std::unique_ptr<OGRGeometry>& resolved,
    double& shiftDistance)
{
    resolved.reset();
    shiftDistance = 0.0;
    if (!fixed || !moving) return false;

    double fixedCx = 0.0, fixedCy = 0.0;
    double movingCx = 0.0, movingCy = 0.0;
    if (!GeometryCentroid2D(fixed, fixedCx, fixedCy) ||
        !GeometryCentroid2D(moving, movingCx, movingCy)) {
        return false;
    }

    double vx = movingCx - fixedCx;
    double vy = movingCy - fixedCy;
    double len = std::hypot(vx, vy);
    if (len < 1e-6) {
        OGREnvelope fixedEnv, movingEnv;
        fixed->getEnvelope(&fixedEnv);
        moving->getEnvelope(&movingEnv);
        const double gapX = std::min(std::abs(movingEnv.MaxX - fixedEnv.MinX),
                                     std::abs(fixedEnv.MaxX - movingEnv.MinX));
        const double gapY = std::min(std::abs(movingEnv.MaxY - fixedEnv.MinY),
                                     std::abs(fixedEnv.MaxY - movingEnv.MinY));
        if (gapX <= gapY) {
            vx = movingCx >= fixedCx ? 1.0 : -1.0;
            vy = 0.0;
        } else {
            vx = 0.0;
            vy = movingCy >= fixedCy ? 1.0 : -1.0;
        }
        len = 1.0;
    }
    vx /= len;
    vy /= len;

    const double areaScale = std::sqrt(std::max(std::min(GeometryArea(fixed), GeometryArea(moving)), 1.0));
    const double baseStep = std::clamp(std::sqrt(std::max(overlapArea, 0.0)) * 0.35, 0.05, 0.6);
    const double maxShift = std::clamp(0.08 * areaScale, 0.5, 3.0);
    const int maxSteps = std::max(4, static_cast<int>(std::ceil(maxShift / baseStep)));
    for (int step = 1; step <= maxSteps; ++step) {
        const double distance = std::min(maxShift, baseStep * step);
        std::unique_ptr<OGRGeometry> candidate = TranslateGeometry2D(moving, vx * distance, vy * distance);
        if (!candidate || candidate->IsEmpty()) continue;
        double remainingOverlap = 0.0;
        if (!HasMeaningfulGeometryOverlap(fixed, candidate.get(), remainingOverlap)) {
            resolved = std::move(candidate);
            shiftDistance = distance;
            return true;
        }
    }
    return false;
}

struct GroupEdgeLine {
    double theta = 0.0;
    double nx = 0.0;
    double ny = 0.0;
    double initialD = 0.0;
    double d = 0.0;
    double insideSign = 1.0;
};

struct GroupBuildingModel {
    GIntBig fid = OGRNullFID;
    std::vector<pcl::PointXYZ> initialRing;
    std::vector<pcl::PointXYZ> ring;
    std::vector<GroupEdgeLine> edges;
    double area = 0.0;
};

struct DRegularizerResidual {
    explicit DRegularizerResidual(double initialD, double weight)
        : initialD_(initialD), weight_(weight) {}

    template <typename T>
    bool operator()(const T* const d, T* residual) const {
        residual[0] = T(weight_) * (d[0] - T(initialD_));
        return true;
    }

    static ceres::CostFunction* Create(double initialD, double weight) {
        return new ceres::AutoDiffCostFunction<DRegularizerResidual, 1, 1>(
            new DRegularizerResidual(initialD, weight));
    }

    double initialD_;
    double weight_;
};

struct IntrusionHalfspaceResidual {
    IntrusionHalfspaceResidual(double x, double y, double nx, double ny,
                               double insideSign, double margin, double weight)
        : x_(x), y_(y), nx_(nx), ny_(ny),
          insideSign_(insideSign), margin_(margin), weight_(weight) {}

    template <typename T>
    bool operator()(const T* const d, T* residual) const {
        const T signedInside = T(insideSign_) * (T(nx_) * T(x_) + T(ny_) * T(y_) - d[0]);
        residual[0] = T(weight_) * (signedInside + T(margin_));
        return true;
    }

    static ceres::CostFunction* Create(double x, double y, double nx, double ny,
                                       double insideSign, double margin, double weight) {
        return new ceres::AutoDiffCostFunction<IntrusionHalfspaceResidual, 1, 1>(
            new IntrusionHalfspaceResidual(x, y, nx, ny, insideSign, margin, weight));
    }

    double x_;
    double y_;
    double nx_;
    double ny_;
    double insideSign_;
    double margin_;
    double weight_;
};

// 作用：两条 Hesse 法线式直线的交点。
bool LineIntersection2D(
    const GroupEdgeLine& a,
    const GroupEdgeLine& b,
    double& x,
    double& y)
{
    const double det = a.nx * b.ny - b.nx * a.ny;
    if (std::abs(det) < 1e-8) return false;
    x = (a.d * b.ny - b.d * a.ny) / det;
    y = (a.nx * b.d - b.nx * a.d) / det;
    return std::isfinite(x) && std::isfinite(y);
}

// 作用：把输出多边形转为组平差模型(边参数化 θ/d)。
bool BuildGroupModelFromGeometry(const OutlineFeatureRecord& feature,
    const Eigen::Vector3d& metadataOffset,
    GroupBuildingModel& model)
{
    model = {};
    if (!feature.geometry || wkbFlatten(feature.geometry->getGeometryType()) != wkbPolygon) {
        return false;
    }
    OGRPolygon* polygon = feature.geometry->toPolygon();
    if (!polygon) return false;

    model.fid = feature.fid;
    model.initialRing = ExtractExteriorRing(polygon, metadataOffset);
    RemoveClosingDuplicate(model.initialRing);
    if (model.initialRing.size() < 3) return false;
    model.ring = model.initialRing;
    model.area = std::max(GeometryArea(feature.geometry.get()), PolygonArea2D(model.ring));

    double cx = 0.0;
    double cy = 0.0;
    for (const auto& p : model.ring) {
        cx += p.x;
        cy += p.y;
    }
    cx /= static_cast<double>(model.ring.size());
    cy /= static_cast<double>(model.ring.size());

    model.edges.reserve(model.ring.size());
    for (std::size_t i = 0; i < model.ring.size(); ++i) {
        const auto& a = model.ring[i];
        const auto& b = model.ring[(i + 1) % model.ring.size()];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double len = std::hypot(dx, dy);
        if (len < 1e-6) return false;

        GroupEdgeLine edge;
        edge.theta = std::atan2(dy, dx);
        edge.nx = -std::sin(edge.theta);
        edge.ny = std::cos(edge.theta);
        edge.initialD = edge.nx * a.x + edge.ny * a.y;
        edge.d = edge.initialD;
        const double centroidSigned = edge.nx * cx + edge.ny * cy - edge.d;
        edge.insideSign = centroidSigned >= 0.0 ? 1.0 : -1.0;
        model.edges.push_back(edge);
    }
    return model.edges.size() == model.ring.size();
}

// 前向声明: ReconstructGroupRing 的简单性护栏用到(定义在下方)
bool MaskOnlyRingIsSimple(const std::vector<pcl::PointXYZ>& ring);

// 作用：由优化后的边参数重建环(面积比护栏)。
bool ReconstructGroupRing(GroupBuildingModel& model)
{
    if (model.edges.size() < 3) return false;
    std::vector<pcl::PointXYZ> rebuilt;
    rebuilt.reserve(model.edges.size());
    const float z = model.initialRing.empty() ? 0.0f : model.initialRing.front().z;
    for (std::size_t i = 0; i < model.edges.size(); ++i) {
        const std::size_t previous = (i + model.edges.size() - 1) % model.edges.size();
        double x = 0.0;
        double y = 0.0;
        if (!LineIntersection2D(model.edges[previous], model.edges[i], x, y)) {
            return false;
        }
        // 交点距离护栏: 近平行边优化后可能产生远距离交点(飞点),
        // 交点必须落在初始环对应顶点附近
        if (i < model.initialRing.size() &&
            std::hypot(x - model.initialRing[i].x, y - model.initialRing[i].y) > 0.8) {
            return false;
        }
        pcl::PointXYZ p;
        p.x = static_cast<float>(x);
        p.y = static_cast<float>(y);
        p.z = z;
        rebuilt.push_back(p);
    }
    RemoveClosingDuplicate(rebuilt);
    if (rebuilt.size() < 3) return false;
    // 简单多边形(组平差不允许引入自交)
    if (!MaskOnlyRingIsSimple(rebuilt)) return false;
    const double area = PolygonArea2D(rebuilt);
    const double initialArea = std::max(PolygonArea2D(model.initialRing), 1e-6);
    if (area < 1e-6 || area / initialArea < 0.5 || area / initialArea > 1.8) {
        return false;
    }
    // 质心位移护栏: 组平差不允许整体搬动建筑
    {
        double cx = 0.0, cy = 0.0, ix = 0.0, iy = 0.0;
        for (const auto& p : rebuilt) { cx += p.x; cy += p.y; }
        for (const auto& p : model.initialRing) { ix += p.x; iy += p.y; }
        cx /= static_cast<double>(rebuilt.size());
        cy /= static_cast<double>(rebuilt.size());
        ix /= static_cast<double>(model.initialRing.size());
        iy /= static_cast<double>(model.initialRing.size());
        if (std::hypot(cx - ix, cy - iy) > 0.5) return false;
    }
    model.ring.swap(rebuilt);
    return true;
}

// 作用：沿环采样顶点/中点用于重叠检测。
std::vector<pcl::PointXYZ> SampleRingForOverlap(const std::vector<pcl::PointXYZ>& ring)
{
    std::vector<pcl::PointXYZ> samples;
    if (ring.size() < 3) return samples;
    samples.reserve(ring.size() * 3);
    double cx = 0.0;
    double cy = 0.0;
    for (const auto& p : ring) {
        samples.push_back(p);
        cx += p.x;
        cy += p.y;
    }
    cx /= static_cast<double>(ring.size());
    cy /= static_cast<double>(ring.size());
    pcl::PointXYZ center;
    center.x = static_cast<float>(cx);
    center.y = static_cast<float>(cy);
    center.z = ring.front().z;
    samples.push_back(center);

    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto& a = ring[i];
        const auto& b = ring[(i + 1) % ring.size()];
        pcl::PointXYZ mid;
        mid.x = 0.5f * (a.x + b.x);
        mid.y = 0.5f * (a.y + b.y);
        mid.z = 0.5f * (a.z + b.z);
        samples.push_back(mid);
    }
    return samples;
}

// 作用：查点在环上最近的边下标。
int NearestEdgeIndex(const pcl::PointXYZ& p, const std::vector<pcl::PointXYZ>& ring)
{
    int best = -1;
    double bestDistance = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const double distance = PointToSegmentDistance2D(p, ring[i], ring[(i + 1) % ring.size()]);
        if (distance < bestDistance) {
            bestDistance = distance;
            best = static_cast<int>(i);
        }
    }
    return best;
}

// 作用：Ceres 组平差：相邻建筑互斥的软约束迭代(角度固定、只调边偏移)。
bool OptimizeCeresConflictGroup(
    std::vector<GroupBuildingModel>& group,
    OutputOverlapRepairStats& stats)
{
    if (group.size() < 2) return false;
    bool anyImproved = false;
    constexpr int kOuterIterations = 4;
    const double margin = 0.08;
    const double overlapWeight = 12.0;
    const double priorWeight = 1.0;

    for (int outer = 0; outer < kOuterIterations; ++outer) {
        ceres::Problem problem;
        int intrusionResiduals = 0;

        for (auto& building : group) {
            for (auto& edge : building.edges) {
                problem.AddResidualBlock(
                    DRegularizerResidual::Create(edge.initialD, priorWeight),
                    new ceres::HuberLoss(0.3),
                    &edge.d);
            }
        }

        for (std::size_t a = 0; a < group.size(); ++a) {
            for (std::size_t b = a + 1; b < group.size(); ++b) {
                const auto samplesA = SampleRingForOverlap(group[a].ring);
                for (const auto& p : samplesA) {
                    if (!PointInPolygon2D(p, group[b].ring)) continue;
                    const int edgeIndex = NearestEdgeIndex(p, group[b].ring);
                    if (edgeIndex < 0) continue;
                    auto& edge = group[b].edges[static_cast<std::size_t>(edgeIndex)];
                    problem.AddResidualBlock(
                        IntrusionHalfspaceResidual::Create(
                            p.x, p.y, edge.nx, edge.ny, edge.insideSign, margin, overlapWeight),
                        new ceres::HuberLoss(0.5),
                        &edge.d);
                    ++intrusionResiduals;
                }

                const auto samplesB = SampleRingForOverlap(group[b].ring);
                for (const auto& p : samplesB) {
                    if (!PointInPolygon2D(p, group[a].ring)) continue;
                    const int edgeIndex = NearestEdgeIndex(p, group[a].ring);
                    if (edgeIndex < 0) continue;
                    auto& edge = group[a].edges[static_cast<std::size_t>(edgeIndex)];
                    problem.AddResidualBlock(
                        IntrusionHalfspaceResidual::Create(
                            p.x, p.y, edge.nx, edge.ny, edge.insideSign, margin, overlapWeight),
                        new ceres::HuberLoss(0.5),
                        &edge.d);
                    ++intrusionResiduals;
                }
            }
        }

        if (intrusionResiduals == 0) break;

        ceres::Solver::Options options;
        options.max_num_iterations = 40;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        // 逐建筑重建: 单体失败(护栏拦截/求交失败)只回退该建筑到
        // 初始环与初始边参数, 不回退整个组; 全部失败才终止
        int rebuiltCount = 0;
        for (auto& building : group) {
            if (ReconstructGroupRing(building)) {
                ++rebuiltCount;
            } else {
                building.ring = building.initialRing;
                for (auto& edge : building.edges) edge.d = edge.initialD;
            }
        }
        if (rebuiltCount == 0) return anyImproved;
        anyImproved = true;
    }
    if (anyImproved) {
        ++stats.optimizedGroups;
    }
    return anyImproved;
}

// 作用：把(多)线几何拆成接缝段列表。
void CollectSeamSegments(const OGRGeometry* geometry,
                         const Eigen::Vector3d& metadataOffset,
                         std::vector<SeamSegment>& segments)
{
    if (!geometry) return;
    const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    if (type == wkbLineString) {
        const auto* line = geometry->toLineString();
        for (int i = 1; i < line->getNumPoints(); ++i) {
            SeamSegment segment;
            segment.a.x = static_cast<float>(line->getX(i - 1) - metadataOffset.x());
            segment.a.y = static_cast<float>(line->getY(i - 1) - metadataOffset.y());
            segment.a.z = 0.0f;
            segment.b.x = static_cast<float>(line->getX(i) - metadataOffset.x());
            segment.b.y = static_cast<float>(line->getY(i) - metadataOffset.y());
            segment.b.z = 0.0f;
            if (Distance2D(segment.a, segment.b) > 1e-6) segments.push_back(segment);
        }
        return;
    }
    if (type == wkbMultiLineString || type == wkbGeometryCollection) {
        const auto* collection = geometry->toGeometryCollection();
        for (int i = 0; i < collection->getNumGeometries(); ++i) {
            CollectSeamSegments(collection->getGeometryRef(i), metadataOffset, segments);
        }
    }
}

// 作用：查询点附近是否存在墙面点。
bool HasWallPointNear(
    const pcl::PointXYZ& query,
    const MyCloudPtr& sampled,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree,
    bool* hasAnyPoint = nullptr)
{
    std::vector<int> indices;
    std::vector<float> distances;
    const int found = kdtree->radiusSearch(query, static_cast<float>(kSeamWallSearchRadius),
                                           indices, distances);
    if (hasAnyPoint) *hasAnyPoint = found > 0;
    if (found <= 0) {
        return false;
    }
    for (int index : indices) {
        if (index < 0 || static_cast<std::size_t>(index) >= sampled->normal->size()) continue;
        if (std::abs(sampled->normal->points[index].normal_z) <= kVerticalNormalMaxAbsZ) {
            return true;
        }
    }
    return false;
}

// 作用：取点附近屋顶高度中位数。
bool MedianRoofHeightNear(
    const pcl::PointXYZ& query,
    const MyCloudPtr& sampled,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree,
    double& height)
{
    std::vector<int> indices;
    std::vector<float> distances;
    kdtree->radiusSearch(query, static_cast<float>(kSeamRoofSearchRadius), indices, distances);
    std::vector<float> heights;
    heights.reserve(indices.size());
    for (int index : indices) {
        if (index < 0 || static_cast<std::size_t>(index) >= sampled->cloud->size() ||
            static_cast<std::size_t>(index) >= sampled->normal->size()) {
            continue;
        }
        if (std::abs(sampled->normal->points[index].normal_z) >= 0.65f) {
            heights.push_back(sampled->cloud->points[index].z);
        }
    }
    if (heights.size() < 3) return false;
    const std::size_t middle = heights.size() / 2;
    std::nth_element(heights.begin(), heights.begin() + middle, heights.end());
    height = heights[middle];
    return true;
}

// 作用：沿共享边采样统计墙面比例与两侧屋顶高差(合并/切分的证据判据)。
SeamEvidence EvaluateSeamEvidence(
    const OGRGeometry* sharedBoundary,
    const MyCloudPtr& sampled,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree,
    const Eigen::Vector3d& metadataOffset)
{
    SeamEvidence evidence;
    if (!sharedBoundary || !sampled || !sampled->cloud || !sampled->normal || !kdtree) {
        return evidence;
    }
    std::vector<SeamSegment> segments;
    CollectSeamSegments(sharedBoundary, metadataOffset, segments);
    int sampleCount = 0;
    int wallCount = 0;
    int heightDifferenceCount = 0;
    int heightPairCount = 0;
    for (const SeamSegment& segment : segments) {
        const double dx = segment.b.x - segment.a.x;
        const double dy = segment.b.y - segment.a.y;
        const double length = std::hypot(dx, dy);
        const int steps = std::max(1, static_cast<int>(std::ceil(length / kSeamSampleStep)));
        const double nx = -dy / length;
        const double ny = dx / length;
        for (int step = 0; step <= steps; ++step) {
            const double t = static_cast<double>(step) / steps;
            pcl::PointXYZ center;
            center.x = static_cast<float>(segment.a.x + t * dx);
            center.y = static_cast<float>(segment.a.y + t * dy);
            center.z = 0.0f;
            bool hasAnyPoint = false;
            const bool hasWallPoint = HasWallPointNear(center, sampled, kdtree, &hasAnyPoint);
            if (!hasAnyPoint) continue;
            ++sampleCount;
            if (hasWallPoint) ++wallCount;

            pcl::PointXYZ sideA = center;
            pcl::PointXYZ sideB = center;
            sideA.x += static_cast<float>(nx * kSeamRoofSideOffset);
            sideA.y += static_cast<float>(ny * kSeamRoofSideOffset);
            sideB.x -= static_cast<float>(nx * kSeamRoofSideOffset);
            sideB.y -= static_cast<float>(ny * kSeamRoofSideOffset);
            double heightA = 0.0;
            double heightB = 0.0;
            if (MedianRoofHeightNear(sideA, sampled, kdtree, heightA) &&
                MedianRoofHeightNear(sideB, sampled, kdtree, heightB)) {
                ++heightPairCount;
                if (std::abs(heightA - heightB) >= kSeamRoofHeightDifference) {
                    ++heightDifferenceCount;
                }
            }
        }
    }
    evidence.modelSampleCount = sampleCount;
    evidence.wallCount = wallCount;
    evidence.heightPairCount = heightPairCount;
    evidence.heightDifferenceCount = heightDifferenceCount;
    if (sampleCount > 0) evidence.wallRatio = static_cast<double>(wallCount) / sampleCount;
    if (heightPairCount > 0) evidence.heightRatio = static_cast<double>(heightDifferenceCount) / heightPairCount;
    return evidence;
}

class DisjointSet {
public:
    explicit DisjointSet(std::size_t size) : parent_(size), rank_(size, 0)
    {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    std::size_t find(std::size_t value)
    {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }

    void unite(std::size_t a, std::size_t b)
    {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (rank_[a] < rank_[b]) std::swap(a, b);
        parent_[b] = a;
        if (rank_[a] == rank_[b]) ++rank_[a];
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<int> rank_;
};

#if 0
// 作用：输出重叠修复：组平差 + 平移回退 + 硬裁剪(当前默认禁用)。
bool ResolveOutputOverlaps(OGRLayer* layer,
    const Eigen::Vector3d& metadataOffset,
    OutputOverlapRepairStats& stats)
{
    stats = {};
    if (!layer) return false;

    std::vector<OutlineFeatureRecord> features;
    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry && !geometry->IsEmpty()) {
            OutlineFeatureRecord record;
            record.fid = feature->GetFID();
            record.geometry.reset(geometry->clone());
            if (record.geometry && !record.geometry->IsValid()) {
                std::unique_ptr<OGRGeometry> fixed(record.geometry->Buffer(0.0));
                if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                    record.geometry = std::move(fixed);
                } else {
                    record.validForCeres = false;
                    std::cerr << "[Output overlap] fid=" << record.fid
                              << " remains invalid after Buffer(0); skip Ceres"
                              << std::endl;
                }
            }
            record.geometry->getEnvelope(&record.envelope);
            record.area = GeometryArea(record.geometry.get());
            record.perimeter = GeometryPerimeter(record.geometry.get());
            if (record.area > 0.0 && record.perimeter > 0.0) {
                features.push_back(std::move(record));
            }
        }
        OGRFeature::DestroyFeature(feature);
    }
    if (features.size() < 2) return true;

    DisjointSet groups(features.size());
    long long initialOverlapPairs = 0;
    for (std::size_t i = 0; i < features.size(); ++i) {
        for (std::size_t j = i + 1; j < features.size(); ++j) {
            const auto& a = features[i];
            const auto& b = features[j];
            if (a.envelope.MaxX < b.envelope.MinX || b.envelope.MaxX < a.envelope.MinX ||
                a.envelope.MaxY < b.envelope.MinY || b.envelope.MaxY < a.envelope.MinY) {
                continue;
            }
            double overlapArea = 0.0;
            if (!HasMeaningfulGeometryOverlap(a.geometry.get(), b.geometry.get(), overlapArea)) {
                continue;
            }
            ++initialOverlapPairs;
            stats.maxOverlapArea = std::max(stats.maxOverlapArea, overlapArea);
            groups.unite(i, j);
        }
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> members;
    for (std::size_t i = 0; i < features.size(); ++i) {
        members[groups.find(i)].push_back(i);
    }

    for (const auto& item : members) {
        const auto& indices = item.second;
        if (indices.size() < 2) continue;

        std::vector<GroupBuildingModel> group;
        group.reserve(indices.size());
        bool groupOk = true;
        for (std::size_t index : indices) {
            const OGRPolygon* sourcePolygon = features[index].geometry &&
                wkbFlatten(features[index].geometry->getGeometryType()) == wkbPolygon
                ? features[index].geometry->toPolygon() : nullptr;
            const OGRLinearRing* sourceRing = sourcePolygon
                ? sourcePolygon->getExteriorRing() : nullptr;
            if (!features[index].validForCeres || !sourceRing ||
                sourceRing->getNumPoints() - 1 > 40) {
                groupOk = false;
                std::cerr << "[Output overlap] skip Ceres group fid="
                          << features[index].fid << " valid="
                          << (features[index].validForCeres ? 1 : 0)
                          << " vertices=" << (sourceRing ? sourceRing->getNumPoints() - 1 : 0)
                          << std::endl;
                break;
            }
            GroupBuildingModel model;
            if (!BuildGroupModelFromGeometry(features[index], metadataOffset, model)) {
                groupOk = false;
                break;
            }
            group.push_back(std::move(model));
        }
        if (!groupOk || group.size() < 2) continue;

        if (!OptimizeCeresConflictGroup(group, stats)) continue;

        for (std::size_t local = 0; local < group.size(); ++local) {
            const std::size_t featureIndex = indices[local];
            std::unique_ptr<OGRPolygon> polygon(MakePolygon(group[local].ring, metadataOffset));
            if (!polygon) continue;
            features[featureIndex].geometry.reset(polygon.release());
            features[featureIndex].geometry->getEnvelope(&features[featureIndex].envelope);
            features[featureIndex].area = GeometryArea(features[featureIndex].geometry.get());
            features[featureIndex].perimeter = GeometryPerimeter(features[featureIndex].geometry.get());
            ++stats.shiftedFeatures;
            for (std::size_t p = 0; p < group[local].ring.size() &&
                    p < group[local].initialRing.size(); ++p) {
                stats.maxShiftDistance = std::max(stats.maxShiftDistance,
                    Distance2D(group[local].ring[p], group[local].initialRing[p]));
            }
        }
    }

    // 在不改动较大轮廓的前提下解决残余冲突。平移
    // is accepted only when it does not create a new conflict with any third building.
    for (std::size_t pass = 0; pass < features.size(); ++pass) {
        bool changed = false;
        for (std::size_t i = 0; i < features.size(); ++i) {
            for (std::size_t j = i + 1; j < features.size(); ++j) {
                double overlapArea = 0.0;
                if (!HasMeaningfulGeometryOverlap(
                        features[i].geometry.get(), features[j].geometry.get(), overlapArea)) {
                    continue;
                }

                const std::size_t fixedIndex = features[i].area >= features[j].area ? i : j;
                const std::size_t movingIndex = fixedIndex == i ? j : i;
                std::unique_ptr<OGRGeometry> candidate;
                double shiftDistance = 0.0;
                bool translated = TryResolvePairByTranslation(
                    features[fixedIndex].geometry.get(),
                    features[movingIndex].geometry.get(),
                    overlapArea, candidate, shiftDistance);
                if (translated) {
                    for (std::size_t k = 0; k < features.size(); ++k) {
                        if (k == movingIndex || k == fixedIndex) continue;
                        double otherOverlap = 0.0;
                        if (HasMeaningfulGeometryOverlap(
                                candidate.get(), features[k].geometry.get(), otherOverlap)) {
                            translated = false;
                            break;
                        }
                    }
                }

                if (translated) {
                    features[movingIndex].geometry = std::move(candidate);
                    stats.maxShiftDistance = std::max(stats.maxShiftDistance, shiftDistance);
                    ++stats.translatedFeatures;
                } else {
                    if (!features[movingIndex].geometry->IsValid()) {
                        std::unique_ptr<OGRGeometry> fixed(
                            features[movingIndex].geometry->Buffer(0.0));
                        if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                            features[movingIndex].geometry = std::move(fixed);
                        }
                    }
                    if (!features[fixedIndex].geometry->IsValid()) {
                        std::unique_ptr<OGRGeometry> fixed(
                            features[fixedIndex].geometry->Buffer(0.0));
                        if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                            features[fixedIndex].geometry = std::move(fixed);
                        }
                    }
                    candidate.reset(features[movingIndex].geometry->Difference(
                        features[fixedIndex].geometry.get()));
                    if (!candidate || candidate->IsEmpty()) {
                        std::cerr << "[Output overlap] Difference removed fid="
                                  << features[movingIndex].fid << std::endl;
                        features[movingIndex].geometry.reset();
                    } else {
                        if (!candidate->IsValid()) {
                            std::unique_ptr<OGRGeometry> fixed(candidate->Buffer(0.0));
                            if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                                candidate = std::move(fixed);
                            } else {
                                std::cerr << "[Output overlap] Difference invalid fid="
                                          << features[movingIndex].fid << std::endl;
                                candidate.reset();
                            }
                        }
                        features[movingIndex].geometry = std::move(candidate);
                    }
                    ++stats.clippedFeatures;
                }

                auto& moving = features[movingIndex];
                if (moving.geometry) {
                    moving.geometry->getEnvelope(&moving.envelope);
                    moving.area = GeometryArea(moving.geometry.get());
                    moving.perimeter = GeometryPerimeter(moving.geometry.get());
                    const double bboxArea =
                        (moving.envelope.MaxX - moving.envelope.MinX) *
                        (moving.envelope.MaxY - moving.envelope.MinY);
                    if (moving.area < kMinOutputPolygonArea ||
                        bboxArea < kMinOutputPolygonBBoxArea) {
                        std::cerr << "[Output filter] drop clipped fid=" << moving.fid
                                  << " area=" << moving.area
                                  << " bbox=" << bboxArea << std::endl;
                        moving.geometry.reset();
                        moving.area = moving.perimeter = 0.0;
                    }
                } else {
                    moving.area = moving.perimeter = 0.0;
                }
                changed = true;
            }
        }
        if (!changed) break;
    }

    if (!layer->TestCapability(OLCRandomWrite)) {
        std::cerr << "[Output overlap] layer does not support random writes" << std::endl;
        return false;
    }
    bool writeOk = true;
    for (const auto& feature : features) {
        OGRFeature* target = layer->GetFeature(feature.fid);
        if (!target) {
            std::cerr << "[Output overlap] GetFeature failed fid=" << feature.fid << std::endl;
            writeOk = false;
            continue;
        }
        const OGRErr setGeometryErr = feature.geometry
            ? target->SetGeometry(feature.geometry.get())
            : target->SetGeometry(nullptr);
        const OGRErr setFeatureErr = setGeometryErr == OGRERR_NONE
            ? layer->SetFeature(target) : setGeometryErr;
        OGRFeature::DestroyFeature(target);
        if (setFeatureErr != OGRERR_NONE) {
            std::cerr << "[Output overlap] SetFeature failed fid=" << feature.fid
                      << " err=" << static_cast<int>(setFeatureErr) << std::endl;
            writeOk = false;
        }
    }

    stats.maxOverlapArea = 0.0;
    for (std::size_t i = 0; i < features.size(); ++i) {
        if (!features[i].geometry) continue;
        for (std::size_t j = i + 1; j < features.size(); ++j) {
            if (!features[j].geometry) continue;
            const auto& a = features[i];
            const auto& b = features[j];
            if (a.envelope.MaxX < b.envelope.MinX || b.envelope.MaxX < a.envelope.MinX ||
                a.envelope.MaxY < b.envelope.MinY || b.envelope.MaxY < a.envelope.MinY) {
                continue;
            }
            ++stats.candidatePairs;
            double overlapArea = 0.0;
            if (HasMeaningfulGeometryOverlap(a.geometry.get(), b.geometry.get(), overlapArea)) {
                ++stats.overlapPairs;
                stats.maxOverlapArea = std::max(stats.maxOverlapArea, overlapArea);
            }
        }
    }
    stats.unresolvedPairs = stats.overlapPairs;
    stats.resolvedPairs = std::max<long long>(0, initialOverlapPairs - stats.overlapPairs);

    if (layer->SyncToDisk() != OGRERR_NONE) {
        std::cerr << "[Output overlap] SyncToDisk failed" << std::endl;
        writeOk = false;
    }
    struct WrittenGeometry {
        GIntBig fid = OGRNullFID;
        std::unique_ptr<OGRGeometry> geometry;
        OGREnvelope envelope = {};
    };
    std::vector<WrittenGeometry> writtenGeometries;
    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        const OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry && !geometry->IsEmpty()) {
            WrittenGeometry written;
            written.fid = feature->GetFID();
            written.geometry.reset(geometry->clone());
            if (written.geometry && !written.geometry->IsValid()) {
                std::unique_ptr<OGRGeometry> fixed(written.geometry->Buffer(0.0));
                if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                    written.geometry = std::move(fixed);
                } else {
                    std::cerr << "[Output overlap audit] invalid fid="
                              << written.fid << std::endl;
                    writeOk = false;
                }
            }
            if (written.geometry) written.geometry->getEnvelope(&written.envelope);
            writtenGeometries.push_back(std::move(written));
        }
        OGRFeature::DestroyFeature(feature);
    }
    stats.overlapPairs = 0;
    stats.maxOverlapArea = 0.0;
    for (std::size_t i = 0; i < writtenGeometries.size(); ++i) {
        for (std::size_t j = i + 1; j < writtenGeometries.size(); ++j) {
            const auto& a = writtenGeometries[i];
            const auto& b = writtenGeometries[j];
            if (!a.geometry || !b.geometry ||
                a.envelope.MaxX < b.envelope.MinX || b.envelope.MaxX < a.envelope.MinX ||
                a.envelope.MaxY < b.envelope.MinY || b.envelope.MaxY < a.envelope.MinY) {
                continue;
            }
            double overlapArea = 0.0;
            if (HasMeaningfulGeometryOverlap(
                    a.geometry.get(), b.geometry.get(), overlapArea)) {
                ++stats.overlapPairs;
                stats.maxOverlapArea = std::max(stats.maxOverlapArea, overlapArea);
                std::cerr << "[Output overlap audit] unresolved fid=" << a.fid
                          << " x " << b.fid << " area=" << overlapArea << std::endl;
            }
        }
    }
    stats.unresolvedPairs = stats.overlapPairs;
    stats.resolvedPairs = std::max<long long>(0, initialOverlapPairs - stats.overlapPairs);
    return writeOk && stats.unresolvedPairs == 0;
}
#endif

// ===== Mask-only 输出重叠解决器 =====
// 在全部单体规则化完成后统一处理建筑间相交(不在单体过程中掩盖错误)。
// 冲突处理顺序: 更好的单体候选(初始轮廓) → 带护栏的 Ceres 组平差 →
// 有位移上限的平移回退 → Difference 裁剪低优先级建筑。
// 小建筑矩形获得优先级加成: 邻居让位, 矩形本体只被裁掉超出部分,
// 不会因首次相交整体回退。
struct MaskOnlyOverlapStats {
    long long candidates = 0;
    long long conflictPairs = 0;
    long long groups = 0;
    long long repaired = 0;
    long long clipped = 0;
    long long unresolved = 0;
    long long bufferRepaired = 0;   // Buffer(0) 修复的几何
    long long candidateSwaps = 0;   // 阶段一: 换回初始轮廓解决的对
    long long groupAdjusted = 0;    // 阶段二: 组平差改写的建筑
    long long translatedCount = 0;  // 阶段三: 平移解决的对
    std::vector<long long> clippedFids;
};

// 保留几何中面积最大的合法面部件(MultiPolygon/集合拆解)
std::unique_ptr<OGRGeometry> KeepLargestValidPolygonPart(
    std::unique_ptr<OGRGeometry> geometry)
{
    if (!geometry) return nullptr;
    const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    if (type == wkbPolygon || type == wkbTriangle) {
        return geometry->IsValid() ? std::move(geometry) : nullptr;
    }
    if (type != wkbMultiPolygon && type != wkbGeometryCollection) {
        return nullptr;
    }
    const auto* collection = geometry->toGeometryCollection();
    double bestArea = 0.0;
    std::unique_ptr<OGRGeometry> best;
    for (int i = 0; collection && i < collection->getNumGeometries(); ++i) {
        const OGRGeometry* part = collection->getGeometryRef(i);
        if (!part || wkbFlatten(part->getGeometryType()) != wkbPolygon) continue;
        if (!part->IsValid() || part->IsEmpty()) continue;
        const double a = GeometryArea(part);
        if (a > bestArea) {
            bestArea = a;
            best.reset(part->clone());
        }
    }
    return best;
}

// 把环上偏离方向系统的 >=0.8m 边绕中点旋回最近合法方向(方案A简化版)
std::vector<pcl::PointXYZ> SnapEdgesToDirectionSystems(
    const std::vector<pcl::PointXYZ>& ring,
    const outlineRegular::DirectionContextOut& dirCtx)
{
    if (ring.size() < 3) return {};
    if (!dirCtx.valid || dirCtx.systemAngles.empty()) return ring;
    std::vector<pcl::PointXYZ> out = ring;
    for (std::size_t i = 0; i < out.size(); ++i) {
        const std::size_t j = (i + 1) % out.size();
        const double len = std::hypot(out[j].x - out[i].x, out[j].y - out[i].y);
        if (len < 0.8) continue;
        const double ang = std::atan2(out[j].y - out[i].y, out[j].x - out[i].x);
        // 折叠角空间最近合法角
        double bestAng = dirCtx.systemAngles.front();
        double bestDist = 1e9;
        for (double a : dirCtx.systemAngles) {
            // 折叠空间中 a 代表 a/a+90 两个真实方向; 取展开后最近者
            for (int k = 0; k < 2; ++k) {
                const double cand = a + k * M_PI / 2.0;
                double d = std::abs(ang - cand);
                while (d > M_PI) d = std::abs(d - 2.0 * M_PI);
                d = std::min(d, M_PI - d);
                if (d < bestDist) { bestDist = d; bestAng = cand; }
            }
        }
        if (bestDist < 0.5 * M_PI / 180.0) continue;  // 已合法
        if (bestDist > 45.0 * M_PI / 180.0) continue;  // 过远旋转会破坏几何
        // 绕边中点旋转两端点
        const double mx = 0.5 * (out[i].x + out[j].x);
        const double my = 0.5 * (out[i].y + out[j].y);
        const double d = std::cos(bestAng), e = std::sin(bestAng);
        out[i].x = static_cast<float>(mx - d * len / 2.0);
        out[i].y = static_cast<float>(my - e * len / 2.0);
        out[j].x = static_cast<float>(mx + d * len / 2.0);
        out[j].y = static_cast<float>(my + e * len / 2.0);
    }
    RemoveClosingDuplicate(out);
    return out;
}

bool ResolveMaskOnlyOutputOverlaps(
    OGRLayer* layer,
    const Eigen::Vector3d& originOffset,
    const std::unordered_map<long long, double>& priorityByFid,
    const std::unordered_map<long long, std::vector<pcl::PointXYZ>>& alternateByFid,
    const std::unordered_map<long long, outlineRegular::DirectionContextOut>& directionByFid,
    MaskOnlyOverlapStats& stats)
{
    stats = {};
    if (!layer) return false;
    // ---- 1. 读取 + 几何有效性检查(必要时 Buffer(0) 修复) ----
    std::vector<OutlineFeatureRecord> features;
    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry && !geometry->IsEmpty()) {
            OutlineFeatureRecord record;
            record.fid = feature->GetFID();
            record.geometry.reset(geometry->clone());
            if (record.geometry && !record.geometry->IsValid()) {
                std::unique_ptr<OGRGeometry> fixed(record.geometry->Buffer(0.0));
                if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                    record.geometry = std::move(fixed);
                    ++stats.bufferRepaired;
                }
            }
            if (!record.geometry || record.geometry->IsEmpty() ||
                !record.geometry->IsValid()) {
                OGRFeature::DestroyFeature(feature);
                continue;
            }
            record.geometry->getEnvelope(&record.envelope);
            record.area = GeometryArea(record.geometry.get());
            record.perimeter = GeometryPerimeter(record.geometry.get());
            if (record.area > 0.0 && record.perimeter > 0.0) {
                features.push_back(std::move(record));
            }
        }
        OGRFeature::DestroyFeature(feature);
    }
    stats.candidates = static_cast<long long>(features.size());
    if (features.size() < 2) {
        std::cout << "[MaskOnlyOverlap] candidates=" << stats.candidates
                  << " conflict_pairs=0 groups=0 repaired=0 unresolved=0" << std::endl;
        return true;
    }

    auto priorityOf = [&](const OutlineFeatureRecord& f) {
        const auto it = priorityByFid.find(static_cast<long long>(f.fid));
        return it == priorityByFid.end() ? 0.5 : it->second;
    };
    // 严格重叠面积(绝对面积, 不带比例豁免; 验收阈值 0.05m²)
    auto strictOverlapArea = [](const OGRGeometry* a, const OGRGeometry* b) {
        if (!a || !b || !a->Intersects(b)) return 0.0;
        std::unique_ptr<OGRGeometry> inter(a->Intersection(b));
        return inter ? GeometryArea(inter.get()) : 0.0;
    };
    constexpr double kStrictOverlapArea = 0.05;
    auto bboxDisjoint = [](const OutlineFeatureRecord& a, const OutlineFeatureRecord& b) {
        return a.envelope.MaxX < b.envelope.MinX || b.envelope.MaxX < a.envelope.MinX ||
               a.envelope.MaxY < b.envelope.MinY || b.envelope.MaxY < a.envelope.MinY;
    };
    auto refreshEnvelope = [](OutlineFeatureRecord& f) {
        if (f.geometry) {
            f.geometry->getEnvelope(&f.envelope);
            f.area = GeometryArea(f.geometry.get());
            f.perimeter = GeometryPerimeter(f.geometry.get());
        } else {
            f.area = f.perimeter = 0.0;
        }
    };
    auto collectConflicts = [&](std::vector<std::pair<std::size_t, std::size_t>>& pairs) {
        pairs.clear();
        for (std::size_t i = 0; i < features.size(); ++i) {
            if (!features[i].geometry) continue;
            for (std::size_t j = i + 1; j < features.size(); ++j) {
                if (!features[j].geometry) continue;
                if (bboxDisjoint(features[i], features[j])) continue;
                if (strictOverlapArea(features[i].geometry.get(),
                                      features[j].geometry.get()) > kStrictOverlapArea) {
                    pairs.emplace_back(i, j);
                }
            }
        }
    };

    // ---- 2. bbox 预筛 + 冲突图 ----
    std::vector<std::pair<std::size_t, std::size_t>> conflicts;
    collectConflicts(conflicts);
    stats.conflictPairs = static_cast<long long>(conflicts.size());
    DisjointSet dset(features.size());
    for (const auto& pr : conflicts) dset.unite(pr.first, pr.second);
    {
        std::unordered_map<std::size_t, int> groupSizes;
        for (std::size_t i = 0; i < features.size(); ++i) ++groupSizes[dset.find(i)];
        for (const auto& gs : groupSizes) {
            if (gs.second >= 2) ++stats.groups;
        }
    }
    std::cout << "[MaskOnlyOverlap] candidates=" << stats.candidates << std::endl;
    std::cout << "[MaskOnlyOverlap] conflict_pairs=" << stats.conflictPairs << std::endl;
    std::cout << "[MaskOnlyOverlap] groups=" << stats.groups << std::endl;
    if (conflicts.empty()) return true;
    // ---- 阶段一(已删除): 换回初始轮廓的备用候选 ----
    // 该通道会用未规则化轮廓替换规则化结果, 违反"初始轮廓不得
    // 直出"约束; 备用候选必须是已过完整方向/拓扑/局部规则性检查
    // 的规则化候选(当前无此类候选, candidate_swaps 恒为 0)。


    // ---- 4. 阶段二: 带护栏的 Ceres 组平差 ----
    {
        std::vector<std::pair<std::size_t, std::size_t>> current;
        collectConflicts(current);
        DisjointSet ceresDset(features.size());
        for (const auto& pr : current) ceresDset.unite(pr.first, pr.second);
        std::unordered_map<std::size_t, std::vector<std::size_t>> members;
        for (std::size_t i = 0; i < features.size(); ++i) {
            members[ceresDset.find(i)].push_back(i);
        }
        for (const auto& item : members) {
            const auto& indices = item.second;
            if (indices.size() < 2) continue;
            std::vector<GroupBuildingModel> group;
            group.reserve(indices.size());
            bool groupOk = true;
            for (std::size_t index : indices) {
                const OGRPolygon* sourcePolygon =
                    features[index].geometry &&
                    wkbFlatten(features[index].geometry->getGeometryType()) == wkbPolygon
                        ? features[index].geometry->toPolygon() : nullptr;
                const OGRLinearRing* sourceRing =
                    sourcePolygon ? sourcePolygon->getExteriorRing() : nullptr;
                if (!sourceRing || sourceRing->getNumPoints() - 1 > 40) {
                    groupOk = false;
                    break;
                }
                GroupBuildingModel model;
                if (!BuildGroupModelFromGeometry(features[index], originOffset, model)) {
                    groupOk = false;
                    break;
                }
                group.push_back(std::move(model));
            }
            if (!groupOk || group.size() < 2) continue;
            // 局部坐标确认: 进入 pcl/Ceres 的坐标必须是测试区域尺度,
            // 出现 3.84e7 量级说明 offset 丢失(float 分辨率仅 ~4m)
            {
                double maxAbsXY = 0.0;
                for (const auto& g : group) {
                    for (const auto& p : g.initialRing) {
                        maxAbsXY = std::max(maxAbsXY,
                            std::max(std::abs((double)p.x), std::abs((double)p.y)));
                    }
                }
                std::cerr << "[MaskOnlyOverlapCeres] fid=" << group.front().fid
                          << " local_coord=1 group_size=" << group.size()
                          << " max_abs_xy=" << maxAbsXY << std::endl;
            }
            OutputOverlapRepairStats ceresStats;
            if (!OptimizeCeresConflictGroup(group, ceresStats)) continue;
            for (std::size_t local = 0; local < group.size(); ++local) {
                const std::size_t featureIndex = indices[local];
                std::unique_ptr<OGRPolygon> polygon(
                    MakePolygon(group[local].ring, originOffset));
                if (!polygon || polygon->IsEmpty() || !polygon->IsValid()) continue;
                features[featureIndex].geometry.reset(polygon.release());
                refreshEnvelope(features[featureIndex]);
                ++stats.groupAdjusted;
            }
        }
    }

    // ---- 5/6. 阶段三+四: 平移回退 → Difference 裁剪(按优先级) ----
    for (int pass = 0; pass < 4; ++pass) {
        std::vector<std::pair<std::size_t, std::size_t>> current;
        collectConflicts(current);
        if (current.empty()) break;
        bool changed = false;
        for (const auto& pr : current) {
            auto& fa = features[pr.first];
            auto& fb = features[pr.second];
            if (!fa.geometry || !fb.geometry) continue;
            const bool aYields = priorityOf(fa) < priorityOf(fb);
            OutlineFeatureRecord& loser = aYields ? fa : fb;
            OutlineFeatureRecord& winner = aYields ? fb : fa;

            // 阶段三: 有位移上限的平移(不得引入新冲突)
            const double overlapArea =
                strictOverlapArea(loser.geometry.get(), winner.geometry.get());
            std::unique_ptr<OGRGeometry> candidate;
            double shiftDistance = 0.0;
            bool translated = TryResolvePairByTranslation(
                winner.geometry.get(), loser.geometry.get(),
                overlapArea, candidate, shiftDistance);
            if (translated) {
                OGREnvelope candEnv;
                candidate->getEnvelope(&candEnv);
                for (std::size_t k = 0; k < features.size() && translated; ++k) {
                    if (k == pr.first || k == pr.second || !features[k].geometry) continue;
                    const auto& env = features[k].envelope;
                    if (candEnv.MaxX < env.MinX || env.MaxX < candEnv.MinX ||
                        candEnv.MaxY < env.MinY || env.MaxY < candEnv.MinY) continue;
                    if (strictOverlapArea(candidate.get(), features[k].geometry.get()) >
                        kStrictOverlapArea) {
                        translated = false;
                    }
                }
                // 平移器内部用较宽的重叠阈值验收, 残余重叠可能仍在
                // 严格阈值之上——用严格阈值复核本对是否真正解决
                if (translated &&
                    strictOverlapArea(candidate.get(), winner.geometry.get()) >
                        kStrictOverlapArea) {
                    translated = false;
                }
            }
            if (translated) {
                loser.geometry = std::move(candidate);
                refreshEnvelope(loser);
                ++stats.translatedCount;
                changed = true;
                continue;
            }

            // 阶段四: Difference 裁剪(低优先级让位; 前后有效性检查 +
            // 方向完整性: 裁剪引入的固定建筑边可能不属于被裁建筑的
            // 方向系统; 违规边旋回本建筑合法方向(方案A简化版),
            // 2轮仍失败则删除低优先级残片
            if (!loser.geometry->IsValid()) {
                std::unique_ptr<OGRGeometry> fixed(loser.geometry->Buffer(0.0));
                if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                    loser.geometry = std::move(fixed);
                }
            }
            if (!winner.geometry->IsValid()) {
                std::unique_ptr<OGRGeometry> fixed(winner.geometry->Buffer(0.0));
                if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                    winner.geometry = std::move(fixed);
                }
            }
            std::unique_ptr<OGRGeometry> clipped(
                loser.geometry->Difference(winner.geometry.get()));
            if (!clipped || clipped->IsEmpty()) {
                loser.geometry.reset();
            } else {
                if (!clipped->IsValid()) {
                    std::unique_ptr<OGRGeometry> fixed(clipped->Buffer(0.0));
                    if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                        clipped = std::move(fixed);
                    }
                }
                // 只保留面积最大的合法面部件(MultiPolygon 拆分)
                clipped = KeepLargestValidPolygonPart(std::move(clipped));
                // 方向完整性: >=0.8m 边须属于本建筑方向系统;
                // 违规边绕中点旋回最近合法方向, 最多 2 轮
                if (clipped) {
                    const auto dirIt = directionByFid.find(
                        static_cast<long long>(loser.fid));
                    if (dirIt != directionByFid.end()) {
                        for (int round = 0; round < 2 && clipped; ++round) {
                            auto localRing = ExtractExteriorRing(
                                clipped->toPolygon(), originOffset);
                            if (localRing.size() < 3) break;
                            const std::string viol =
                                outlineRegular::CheckFallbackLocalRegularity(
                                    localRing, dirIt->second, -1, 0, "clip");
                            if (viol.empty()) break;
                            std::cerr << "[MaskOnlyOverlap] clip_fix fid="
                                      << loser.fid << " round=" << round
                                      << " reason=" << viol << std::endl;
                            auto fixedRing = SnapEdgesToDirectionSystems(
                                localRing, dirIt->second);
                            if (fixedRing.size() < 3) { clipped.reset(); break; }
                            std::unique_ptr<OGRPolygon> fixedPoly(
                                MakePolygon(fixedRing, originOffset));
                            if (!fixedPoly || fixedPoly->IsEmpty() ||
                                !fixedPoly->IsValid() ||
                                strictOverlapArea(fixedPoly.get(),
                                    winner.geometry.get()) > 0.05) {
                                if (round == 1) clipped.reset();
                                continue;
                            }
                            clipped = std::move(fixedPoly);
                        }
                    }
                }
                if (clipped) {
                    loser.geometry = std::move(clipped);
                } else {
                    loser.geometry.reset();
                }
            }
            refreshEnvelope(loser);
            ++stats.clipped;
            stats.clippedFids.push_back(static_cast<long long>(loser.fid));
            // 裁剪后面积地板: <15m² 或 bbox<20m² 的残片删除
            if (loser.geometry) {
                const double bboxArea =
                    (loser.envelope.MaxX - loser.envelope.MinX) *
                    (loser.envelope.MaxY - loser.envelope.MinY);
                if (loser.area < 15.0 || bboxArea < 20.0) {
                    loser.geometry.reset();
                    refreshEnvelope(loser);
                }
            }
            changed = true;
        }
        if (!changed) break;
    }
    stats.repaired = std::max<long long>(0, stats.conflictPairs -
        [&] {
            std::vector<std::pair<std::size_t, std::size_t>> current;
            collectConflicts(current);
            return static_cast<long long>(current.size());
        }());

    // ---- 7. 写回 ----
    if (!layer->TestCapability(OLCRandomWrite)) {
        std::cerr << "[MaskOnlyOverlap] layer does not support random writes" << std::endl;
        return false;
    }
    for (const auto& feature : features) {
        OGRFeature* target = layer->GetFeature(feature.fid);
        if (!target) continue;
        const OGRErr setGeometryErr = feature.geometry
            ? target->SetGeometry(feature.geometry.get())
            : target->SetGeometry(nullptr);
        if (setGeometryErr == OGRERR_NONE) {
            const OGRErr writeErr = layer->SetFeature(target);
            if (writeErr != OGRERR_NONE) {
                std::cerr << "[MaskOnlyOverlap] SetFeature failed fid="
                          << feature.fid << " err=" << writeErr << std::endl;
            }
        } else {
            std::cerr << "[MaskOnlyOverlap] SetGeometry failed fid="
                      << feature.fid << " err=" << setGeometryErr << std::endl;
        }
        OGRFeature::DestroyFeature(target);
    }

    // Shapefile 的 SetFeature 可能延迟写入 .shp/.shx；必须在最终审计前
    // 强制落盘，否则同一 layer 的读缓存可能看到“已修复”、关闭后文件
    // 却仍保留旧几何。
    if (layer->SyncToDisk() != OGRERR_NONE) {
        std::cerr << "[MaskOnlyOverlap] SyncToDisk failed after write-back"
                  << std::endl;
    }
    layer->ResetReading();

    // ---- 8. 最终严格检查: 先有效性修复, 再绝对面积相交检测 ----
    {
        std::vector<OutlineFeatureRecord> finalFeats;
        layer->ResetReading();
        while (OGRFeature* feature = layer->GetNextFeature()) {
            OGRGeometry* geometry = feature->GetGeometryRef();
            if (geometry && !geometry->IsEmpty()) {
                OutlineFeatureRecord record;
                record.fid = feature->GetFID();
                record.geometry.reset(geometry->clone());
                if (record.geometry && !record.geometry->IsValid()) {
                    std::unique_ptr<OGRGeometry> fixed(record.geometry->Buffer(0.0));
                    if (fixed && !fixed->IsEmpty() && fixed->IsValid()) {
                        record.geometry = std::move(fixed);
                        OGRFeature* target = layer->GetFeature(feature->GetFID());
                        if (target) {
                            target->SetGeometry(record.geometry.get());
                            layer->SetFeature(target);
                            OGRFeature::DestroyFeature(target);
                        }
                    }
                }
                if (record.geometry && !record.geometry->IsEmpty()) {
                    record.geometry->getEnvelope(&record.envelope);
                    finalFeats.push_back(std::move(record));
                }
            }
            OGRFeature::DestroyFeature(feature);
        }
        for (std::size_t i = 0; i < finalFeats.size(); ++i) {
            for (std::size_t j = i + 1; j < finalFeats.size(); ++j) {
                if (bboxDisjoint(finalFeats[i], finalFeats[j])) continue;
                const double area = strictOverlapArea(
                    finalFeats[i].geometry.get(), finalFeats[j].geometry.get());
                if (area > kStrictOverlapArea) {
                    ++stats.unresolved;
                    std::cerr << "[MaskOnlyOverlap] unresolved fid="
                              << finalFeats[i].fid << " x " << finalFeats[j].fid
                              << " area=" << area << std::endl;
                }
            }
        }
    }

    std::cout << "[MaskOnlyOverlap] repaired=" << stats.repaired << std::endl;
    std::cout << "[MaskOnlyOverlap] unresolved=" << stats.unresolved << std::endl;
    std::string clippedList;
    for (std::size_t i = 0; i < stats.clippedFids.size(); ++i) {
        clippedList += (i ? "," : "") + std::to_string(stats.clippedFids[i]);
    }
    std::cout << "[MaskOnlyOverlap] clipped_fids=" << clippedList << std::endl;
    return stats.unresolved == 0;
}

// 作用：初始轮廓过分割合并阶段：按接缝证据/形状/大建筑兜底就地合并 Shapefile。
// GeometryOnly 模式仅执行同 parent + 共享边下限的确定性合并，不调 OSGB 证据。
enum class InitialMergeMode {
    GeometryOnly,
    GeometryAndOSGB
};

bool MergeOversegmentedInitialOutlines(
    const std::string& shpPath,
    const MyCloudPtr& sampled,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree,
    const Eigen::Vector3d& metadataOffset,
    InitialOutlineMergeStats& stats,
    InitialMergeMode mode = InitialMergeMode::GeometryAndOSGB)
{
    stats = {};
    GDALDataset* dataset = static_cast<GDALDataset*>(
        GDALOpenEx(shpPath.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
                   nullptr, nullptr, nullptr));
    if (!dataset) return false;
    OGRLayer* layer = dataset->GetLayer(0);
    if (!layer) {
        GDALClose(dataset);
        return false;
    }

    std::vector<OutlineFeatureRecord> features;
    const int parentFieldIdx = layer->GetLayerDefn()->GetFieldIndex("parent");
    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry && !geometry->IsEmpty()) {
            OutlineFeatureRecord record;
            record.fid = feature->GetFID();
            record.geometry.reset(geometry->clone());
            record.geometry->getEnvelope(&record.envelope);
            record.area = GeometryArea(record.geometry.get());
            record.perimeter = GeometryPerimeter(record.geometry.get());
            if (parentFieldIdx >= 0) record.parent = feature->GetFieldAsInteger(parentFieldIdx);
            if (record.area > 0.0 && record.perimeter > 0.0) {
                features.push_back(std::move(record));
            }
        }
        OGRFeature::DestroyFeature(feature);
    }

    DisjointSet groups(features.size());
    for (std::size_t i = 0; i < features.size(); ++i) {
        for (std::size_t j = i + 1; j < features.size(); ++j) {
            const auto& a = features[i];
            const auto& b = features[j];
            if (a.envelope.MaxX < b.envelope.MinX || b.envelope.MaxX < a.envelope.MinX ||
                a.envelope.MaxY < b.envelope.MinY || b.envelope.MaxY < a.envelope.MinY) {
                continue;
            }
            std::unique_ptr<OGRGeometry> boundaryA(a.geometry->Boundary());
            std::unique_ptr<OGRGeometry> boundaryB(b.geometry->Boundary());
            std::unique_ptr<OGRGeometry> shared(
                boundaryA && boundaryB ? boundaryA->Intersection(boundaryB.get()) : nullptr);
            const double exactSharedLength = GeometryLength(shared.get());
            const double sharedLength = RobustSharedBoundaryLength(
                a.geometry.get(), b.geometry.get(), exactSharedLength);
            if (sharedLength < kMinSharedBoundaryLength) continue;   // 邻接下限(点/角接触不算)
            ++stats.adjacentCandidates;
            if (sharedLength > exactSharedLength + 0.5) {
                std::cerr << "[Merge diag] ROBUST exactShared=" << exactSharedLength
                          << " robustShared=" << sharedLength
                          << " areaA=" << a.area << " areaB=" << b.area << std::endl;
            }

            if (mode == InitialMergeMode::GeometryOnly) {
                // 保守几何合并：同 parent(同实例色)直接合并。
                // 兜底两条纯几何判据(不依赖 OSGB 点云)，覆盖 parent 因
                // 颜色漂移/连通域断裂而失真的过分割：
                //   1) 形状互补(并集面积/周长吻合) —— 中小建筑的主判据；
                //   2) 大建筑长共享边 —— min面积>1500m² 且共享边>=3m。
                // wrapper 伪影会被形状互补校验(并集须为单多边形)拦住。
                if (a.parent > 0 && a.parent == b.parent) {
                    groups.unite(i, j);
                    ++stats.mergedPairs;
                    ++stats.mergedByParent;
                    continue;
                }
                if (ShouldMergeByFootprintShape(a, b, sharedLength)) {
                    groups.unite(i, j);
                    ++stats.mergedPairs;
                    ++stats.mergedByShape;
                    continue;
                }
                if (std::min(a.area, b.area) > kLargeBuildingMergeArea &&
                    sharedLength >= kMergeSeamLength) {
                    groups.unite(i, j);
                    ++stats.mergedPairs;
                    ++stats.mergedByBigPair;
                    continue;
                }
                continue;
            }

            const bool largePair = std::min(a.area, b.area) > 500.0;  // 澶у潡瀵?鐤戜技澶у缓绛戯紝鎵撹瘖鏂?
bool evidenceComputed = false;
            SeamEvidence evidence;
            auto getEvidence = [&]() -> const SeamEvidence& {
                if (!evidenceComputed) {
                    evidence = EvaluateSeamEvidence(
                        shared.get(), sampled, kdtree, metadataOffset);
                    evidenceComputed = true;
                }
                return evidence;
            };

            if (ShouldMergeByFootprintShape(a, b, sharedLength)) {
                if (largePair) std::cerr << "[Merge diag] SHAPE sharedLen=" << sharedLength
                                         << " areaA=" << a.area << " areaB=" << b.area << std::endl;
                groups.unite(i, j);
                ++stats.mergedPairs;
                ++stats.mergedByShape;
                continue;
            }

            // 宽缝：共享边 >= kMergeSeamLength 才算"该合"。
            // 窄缝(如分水岭切开的粘连)共享边短 → 属不同栋 → 不合并。
if (sharedLength < kMergeSeamLength) {
                if (largePair) std::cerr << "[Merge diag] NARROW sharedLen=" << sharedLength
                                         << " areaA=" << a.area << " areaB=" << b.area << std::endl;
                ++stats.keptNarrowSeam;
                continue;
            }
            const bool sameParent = a.parent > 0 && a.parent == b.parent;
            if (sameParent) {
                groups.unite(i, j);
                ++stats.mergedPairs;
                ++stats.mergedByParent;
                continue;
            }
            // 大建筑兜底：两侧都是大块且共享长边 → 视为同一综合体，跳过证据间直接合并。
            // 注(kLargeBuildingMergeArea 处)：长边本身是同一栋的最强几何证据，而证据间在
            // 大建筑屋面会因女儿墙/伸缩缝等触发 wallRatio 饱和而错误阻止合并。
if (std::min(a.area, b.area) > kLargeBuildingMergeArea) {
                if (largePair) std::cerr << "[Merge diag] BIGMERGE sharedLen=" << sharedLength
                                         << " areaA=" << a.area << " areaB=" << b.area << std::endl;
                groups.unite(i, j);
                ++stats.mergedPairs;
                ++stats.mergedByBigPair;
                continue;
            }
            const SeamEvidence& seamEvidence = getEvidence();
            const double maxEvidence = std::max(seamEvidence.wallRatio, seamEvidence.heightRatio);
            if (maxEvidence >= kStrongSeamEvidenceRatio) {
                if (largePair) std::cerr << "[Merge diag] EVIDENCE sharedLen=" << sharedLength
                                         << " wall=" << seamEvidence.wallRatio << " height=" << seamEvidence.heightRatio
                                         << " areaA=" << a.area << " areaB=" << b.area << std::endl;
                ++stats.keptByModelEvidence;
                continue;
            }
            if (largePair) std::cerr << "[Merge diag] MERGE sharedLen=" << sharedLength
                                     << " wall=" << seamEvidence.wallRatio << " height=" << seamEvidence.heightRatio
                                     << " areaA=" << a.area << " areaB=" << b.area << std::endl;
            groups.unite(i, j);
            ++stats.mergedPairs;
        }
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> members;
    for (std::size_t i = 0; i < features.size(); ++i) members[groups.find(i)].push_back(i);
    for (const auto& item : members) {
        const auto& indices = item.second;
        if (indices.size() < 2) continue;
        std::unique_ptr<OGRGeometry> merged(features[indices.front()].geometry->clone());
        for (std::size_t k = 1; k < indices.size(); ++k) {
            std::unique_ptr<OGRGeometry> next(merged->Union(features[indices[k]].geometry.get()));
            if (next) merged = std::move(next);
        }
        std::unique_ptr<OGRGeometry> mergeCleaned = CleanMergedGeometry(merged.get());
        if (mergeCleaned) merged = std::move(mergeCleaned);
        // 点接触/小缝隙会让并集变成 MultiPolygon；用 Buffer(0) 融合成单多边形
if (merged && wkbFlatten(merged->getGeometryType()) != wkbPolygon) {
            std::unique_ptr<OGRGeometry> cleaned(merged->Buffer(0.0));
            if (cleaned && wkbFlatten(cleaned->getGeometryType()) == wkbPolygon) {
                merged = std::move(cleaned);
            } else {
                continue;  // 实在合不成单多边形，放弃这组
            }
        }
        if (!merged || wkbFlatten(merged->getGeometryType()) != wkbPolygon) continue;

        OGRFeature* target = layer->GetFeature(features[indices.front()].fid);
        if (!target) continue;
        target->SetGeometry(merged.get());
        const bool updated = layer->SetFeature(target) == OGRERR_NONE;
        OGRFeature::DestroyFeature(target);
        if (!updated) continue;
        for (std::size_t k = 1; k < indices.size(); ++k) {
            if (layer->DeleteFeature(features[indices[k]].fid) == OGRERR_NONE) {
                ++stats.removedFeatures;
            }
        }
    }
    // Convergence guard: if pairs were identified but no features were actually
    GDALClose(dataset);
    return true;
}

// 把多边形(可能带洞)的所有 interior ring 去掉，只保留外环(实心化)。changed 置为是否真的改动过。
std::unique_ptr<OGRGeometry> SolidGeometry(const OGRGeometry* g, bool& changed)
{
    changed = false;
    if (!g) return nullptr;
    const OGRwkbGeometryType t = wkbFlatten(g->getGeometryType());
    if (t == wkbPolygon) {
        const OGRPolygon* p = g->toPolygon();
        if (p->getNumInteriorRings() == 0) return std::unique_ptr<OGRGeometry>(g->clone());
        auto* np = new OGRPolygon();
        if (const OGRLinearRing* ext = p->getExteriorRing())
            np->addRingDirectly(new OGRLinearRing(*ext));
        np->assignSpatialReference(g->getSpatialReference());
        changed = true;
        return std::unique_ptr<OGRGeometry>(np);
    }
    if (t == wkbMultiPolygon) {
        const OGRMultiPolygon* mp = g->toMultiPolygon();
        bool anyHole = false;
        for (int i = 0; i < mp->getNumGeometries(); ++i) {
            const OGRPolygon* p = mp->getGeometryRef(i)->toPolygon();
            if (p && p->getNumInteriorRings() > 0) { anyHole = true; break;
        }
        }
        if (!anyHole) return std::unique_ptr<OGRGeometry>(g->clone());
        auto* nmp = new OGRMultiPolygon();
        for (int i = 0; i < mp->getNumGeometries(); ++i) {
            const OGRPolygon* p = mp->getGeometryRef(i)->toPolygon();
            if (!p || !p->getExteriorRing()) continue;
            auto* np = new OGRPolygon();
            np->addRingDirectly(new OGRLinearRing(*(p->getExteriorRing())));
            nmp->addGeometryDirectly(np);
        }
        nmp->assignSpatialReference(g->getSpatialReference());
        changed = true;
        return std::unique_ptr<OGRGeometry>(nmp);
    }
    return std::unique_ptr<OGRGeometry>(g->clone());
}

// 合并收敛后：先把所有多边形去 interior ring(实心化)，再删除被较大轮廓完全包含的小轮廓。
// 大建筑合并后，内部残留的分色/栅格化小碎片常以"大轮廓带洞+洞里套小轮廓"的形式存在；
// 带洞的 A 对洞里的 B 调 Contains 返回 false 会漏检，故先全局去洞，再按实际几何做包含清理。
// maxArea：只删面积 <= maxArea 的被包含者。
// solidifyHoles：true(默认/OSGB)=先实心化再判包含(旧行为)；
//                false(Mask-only)=保留原始孔洞，用外环壳体做候选判定。
long long RemoveContainedSmallFootprints(
    const std::string& shpPath, double maxArea, bool solidifyHoles = true)
{
    GDALDataset* dataset = static_cast<GDALDataset*>(
        GDALOpenEx(shpPath.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
                   nullptr, nullptr, nullptr));
    if (!dataset) return 0;
    OGRLayer* layer = dataset->GetLayer(0);
    if (!layer) { GDALClose(dataset); return 0; }

    struct ContainedRec {
        GIntBig fid = OGRNullFID;
        std::unique_ptr<OGRGeometry> geom;
        OGREnvelope env = {};
        double area = 0.0;
        bool solidified = false;              // 原带孔、已实心化，需写回
        bool remove = false;
    };
    std::vector<ContainedRec> recs;
    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry && !geometry->IsEmpty()) {
            ContainedRec rec;
            rec.fid = feature->GetFID();
            bool changed = false;
            if (solidifyHoles) {
                rec.geom = SolidGeometry(geometry, changed);   // 先实心化，让洞里的小轮廓可被 Contains
                rec.solidified = changed;
            } else {
                // Mask-only：保留原始几何不实心化；用外环壳体副本做包含候选判定。
                if (wkbFlatten(geometry->getGeometryType()) == wkbPolygon) {
                    OGRPolygon* polygon = geometry->toPolygon();
                    if (polygon->getExteriorRing()) {
                        auto* shell = new OGRPolygon();
                        shell->addRing(new OGRLinearRing(*polygon->getExteriorRing()));
                        rec.geom.reset(shell);
                    }
                }
                if (!rec.geom) rec.geom.reset(geometry->clone());
            }
            rec.geom->getEnvelope(&rec.env);
            rec.area = GeometryArea(rec.geom.get());
            if (rec.area > 0.0) recs.push_back(std::move(rec));
        }
        OGRFeature::DestroyFeature(feature);
    }

    std::sort(recs.begin(), recs.end(),
              [](const ContainedRec& a, const ContainedRec& b) { return a.area > b.area; });
    for (std::size_t i = 1; i < recs.size(); ++i) {
        ContainedRec& b = recs[i];
        if (b.area > maxArea) continue;          // 只删"小"轮廓
        for (std::size_t j = 0; j < i; ++j) {
            const ContainedRec& a = recs[j];
            if (a.area <= b.area) continue;
            if (a.env.MinX > b.env.MinX || a.env.MinY > b.env.MinY ||
                a.env.MaxX < b.env.MaxX || a.env.MaxY < b.env.MaxY) {
                continue;                         // envelope 预筛：A 包不住 B 的外接矩形则跳过
            }
            if (a.geom->Contains(b.geom.get())) {
                b.remove = true;
                break;
            }
        }
    }

    long long removed = 0;
    long long holesPreserved = 0;
    long long solidified = 0;
    for (const ContainedRec& rec : recs) {
        if (rec.remove) {
            if (rec.fid != OGRNullFID &&
                layer->DeleteFeature(rec.fid) == OGRERR_NONE) ++removed;
        } else if (rec.solidified && rec.fid != OGRNullFID) {
            if (OGRFeature* feature = layer->GetFeature(rec.fid)) {
                feature->SetGeometry(rec.geom.get());
                if (layer->SetFeature(feature) == OGRERR_NONE) ++solidified;
                OGRFeature::DestroyFeature(feature);
            }
        }
    }
    if (removed > 0 || solidified > 0) {
        layer->SyncToDisk();
        std::cout << "[Mask] solidified (removed holes) in " << solidified
                  << " footprints." << std::endl;
    }
        if (!solidifyHoles) {
        std::cout << "[MaskOnly] holes_preserved=" << holesPreserved
                  << " contained_removed=" << removed << std::endl;
    }
GDALClose(dataset);
    return removed;
}

// 作用：把支撑点按要素 fid 着色后追加进调试收集器。
void AppendDebugSupportPoints(
    const std::vector<pcl::PointXYZ>& support,
    GIntBig sourceFid,
    int ringIndex,
    RegularizationDebugCollector* debug)
{
    if (!debug || support.empty()) return;
    const unsigned int hash = static_cast<unsigned int>(
        (sourceFid == OGRNullFID ? 0 : sourceFid) * 2654435761u +
        static_cast<unsigned int>(ringIndex) * 97531u);
    const std::uint8_t r = static_cast<std::uint8_t>(80 + (hash & 0x7F));
    const std::uint8_t g = static_cast<std::uint8_t>(80 + ((hash >> 8) & 0x7F));
    const std::uint8_t b = static_cast<std::uint8_t>(80 + ((hash >> 16) & 0x7F));
    debug->support->reserve(debug->support->size() + support.size());
    for (const auto& p : support) {
        PointT q;
        q.x = p.x;
        q.y = p.y;
        q.z = p.z;
        q.r = r;
        q.g = g;
        q.b = b;
        q.a = 255;
        debug->support->push_back(q);
    }
}

// ===== RegularizeGeometry =====
// 作用：对一个 OGR 几何(单多边形/多多边形)做规则化，返回新的几何。
//       其它类型直接克隆返回。
std::unique_ptr<OGRGeometry> RegularizeGeometry(
    OGRGeometry* geometry,
    const MyCloudPtr& sampled,
    const pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr& kdtree,
    const Eigen::Vector3d& metadataOffset,
    const SupportOwnershipContext* ownership = nullptr,
    GIntBig sourceFid = OGRNullFID,
    RegularizationDebugCollector* debug = nullptr)
{
    if (!geometry) return nullptr;

    OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
    if (type == wkbPolygon) {
        auto* polygon = geometry->toPolygon();
        auto ring = ExtractExteriorRing(polygon, metadataOffset);
        if (BoundingBoxArea2D(ring) < kMinPolygonBBoxArea) {
            ++g_removedSmallPolygons;
            return nullptr;
        }
        if (!RingHasModelCoverage(ring, sampled, kdtree)) {
            ++g_removedOutsideModel;
            return nullptr;
        }
        std::vector<pcl::PointXYZ> bestHypothesis;
        std::vector<pcl::PointXYZ> support;
        const std::size_t ringId = FindInitialRingRecordId(ownership, sourceFid, 0);
        auto regularized = RegularizeRing(ring, sampled, kdtree,
            ownership, ringId,
            debug ? &bestHypothesis : nullptr,
            debug ? &support : nullptr, sourceFid);
        if (debug && bestHypothesis.size() >= 3) {
            debug->hypotheses.push_back({sourceFid, 0, bestHypothesis});
            AppendDebugSupportPoints(support, sourceFid, 0, debug);
        }
        if (!OutputRingPassesSizeFloor(regularized)) {
            ++g_removedSmallPolygons;
            std::cerr << "[Output filter] drop fid=" << sourceFid
                      << " area=" << PolygonArea2D(regularized)
                      << " bbox=" << BoundingBoxArea2D(regularized) << std::endl;
            return nullptr;
        }
        return std::unique_ptr<OGRGeometry>(MakePolygon(regularized, metadataOffset));
    }

    if (type == wkbMultiPolygon) {
        auto* out = new OGRMultiPolygon();
        auto* multi = geometry->toMultiPolygon();
        int partIndex = 0;
        for (auto&& part : *multi) {
            auto ring = ExtractExteriorRing(part->toPolygon(), metadataOffset);
            if (BoundingBoxArea2D(ring) < kMinPolygonBBoxArea) {
                ++g_removedSmallPolygons;
                ++partIndex;
                continue;
            }
            if (!RingHasModelCoverage(ring, sampled, kdtree)) {
                ++g_removedOutsideModel;
                ++partIndex;
                continue;
            }
            std::vector<pcl::PointXYZ> bestHypothesis;
            std::vector<pcl::PointXYZ> support;
            const std::size_t ringId = FindInitialRingRecordId(ownership, sourceFid, partIndex);
            auto regularized = RegularizeRing(ring, sampled, kdtree,
                ownership, ringId,
                debug ? &bestHypothesis : nullptr,
                debug ? &support : nullptr, sourceFid);
            if (debug && bestHypothesis.size() >= 3) {
                debug->hypotheses.push_back({sourceFid, partIndex, bestHypothesis});
                AppendDebugSupportPoints(support, sourceFid, partIndex, debug);
            }
            if (!OutputRingPassesSizeFloor(regularized)) {
                ++g_removedSmallPolygons;
                std::cerr << "[Output filter] drop fid=" << sourceFid
                          << " part=" << partIndex
                          << " area=" << PolygonArea2D(regularized)
                          << " bbox=" << BoundingBoxArea2D(regularized) << std::endl;
                ++partIndex;
                continue;
            }
            std::unique_ptr<OGRPolygon> polygon(MakePolygon(regularized, metadataOffset));
            if (polygon) {
                out->addGeometryDirectly(polygon.release());
            }
            ++partIndex;
        }
        if (out->getNumGeometries() == 0) {
            delete out;
            return nullptr;
        }
        return std::unique_ptr<OGRGeometry>(out);
    }

    return std::unique_ptr<OGRGeometry>(geometry->clone());
}

// ===== CopyFields =====
// 作用：把输入图层的所有属性字段结构复制到输出图层(只复制字段定义，不复制值)。
void CopyFields(OGRLayer* inputLayer, OGRLayer* outputLayer)
{
    OGRFeatureDefn* defn = inputLayer->GetLayerDefn();
    for (int i = 0; i < defn->GetFieldCount(); ++i) {
        OGRFieldDefn* field = defn->GetFieldDefn(i);
        outputLayer->CreateField(field);
    }
}

// ===== CopyFieldValues =====
// 作用：把一个要素的所有属性值复制到另一个要素(按字段下标对应)。
void CopyFieldValues(OGRFeature* src, OGRFeature* dst)
{
    for (int i = 0; i < src->GetFieldCount(); ++i) {
        dst->SetField(i, src->GetRawFieldRef(i));
    }
}

// ===== Mask-only 初始轮廓保拓扑平滑 =====
// 用 OGR 的 SimplifyPreserveTopology 消除栅格楼梯，配合面积/有效性/
// 邻接回退检查，保证建筑数量、孔洞和相邻关系不被破坏。
// 窄颈拆分必须在平滑之前执行(平滑会删除窄颈候选顶点)。

struct InitialOutlineSmoothStats {
    long long inspected = 0;
    long long changed = 0;
    long long beforeVertices = 0;
    long long afterVertices = 0;
    long long revertedInvalid = 0;
    long long revertedType = 0;
    long long revertedArea = 0;
    long long revertedDeviation = 0;
    long long revertedAdjacency = 0;
    long long revertedNoReduction = 0;
    long long holesPreserved = 0;
};

const double kInitialOutlineSmoothMinTolerance = 0.20;      // m
const double kInitialOutlineSmoothMaxTolerance = 0.45;      // m
const double kInitialOutlineSmoothAreaChangeRatio = 0.05;   // 5%
const double kInitialOutlineSmoothMaxDeviation = 0.60;      // m, 双向边界最大偏移
const double kInitialOutlineSmoothAdjacencyTolerance = 0.60; // m, 共享边邻接距离
const double kInitialOutlineSmoothOverlapArea = 0.05;       // m², 新增重叠面积下限

// 采样环边界并计算到另一环边界的最大距离(近似 Hausdorff)
double MaxBoundaryDeviation(
    const OGRGeometry* geomA,
    const OGRGeometry* geomB,
    double sampleStep)
{
    if (!geomA || !geomB) return 0.0;
    std::unique_ptr<OGRGeometry> boundaryB(geomB->Boundary());
    if (!boundaryB) return 0.0;
    std::unique_ptr<OGRGeometry> boundaryA(geomA->Boundary());
    if (!boundaryA) return 0.0;

    // 沿 A 的边界采样点，测到 B 边界的距离
    OGRPoint sample;
    double maxDist = 0.0;
    const OGRGeometryCollection* colA =
        dynamic_cast<const OGRGeometryCollection*>(boundaryA.get());
    const int nGeomA = colA ? colA->getNumGeometries() : 0;
    for (int g = 0; g < (nGeomA > 0 ? nGeomA : 1); ++g) {
        const OGRGeometry* line =
            (colA && nGeomA > 0) ? colA->getGeometryRef(g) : boundaryA.get();
        if (!line || wkbFlatten(line->getGeometryType()) != wkbLineString) continue;
        const auto* ls = line->toLineString();
        const int nPts = ls->getNumPoints();
        for (int i = 0; i < nPts; ++i) {
            ls->getPoint(i, &sample);
            const double d = boundaryB->Distance(&sample);
            if (d > maxDist) maxDist = d;
        }
    }
    return maxDist;
}

int CountTotalVertices(const OGRGeometry* geom)
{
    if (!geom) return 0;
    const OGRwkbGeometryType type = wkbFlatten(geom->getGeometryType());
    if (type == wkbPolygon) {
        const auto* poly = geom->toPolygon();
        int count = poly->getExteriorRing() ? poly->getExteriorRing()->getNumPoints() : 0;
        for (int i = 0; i < poly->getNumInteriorRings(); ++i) {
            count += poly->getInteriorRing(i)->getNumPoints();
        }
        return count;
    }
    if (type == wkbMultiPolygon) {
        const auto* mp = geom->toMultiPolygon();
        int count = 0;
        for (int i = 0; i < mp->getNumGeometries(); ++i) {
            count += CountTotalVertices(mp->getGeometryRef(i));
        }
        return count;
    }
    return 0;
}

int CountInteriorRings(const OGRGeometry* geom)
{
    if (!geom) return 0;
    const OGRwkbGeometryType type = wkbFlatten(geom->getGeometryType());
    if (type == wkbPolygon) {
        return geom->toPolygon()->getNumInteriorRings();
    }
    if (type == wkbMultiPolygon) {
        const auto* mp = geom->toMultiPolygon();
        int count = 0;
        for (int i = 0; i < mp->getNumGeometries(); ++i) {
            count += CountInteriorRings(mp->getGeometryRef(i));
        }
        return count;
    }
    return 0;
}

// 对整个初始轮廓 Shapefile 执行保拓扑平滑(就地修改)。
// 返回 true 表示至少处理了一个要素。
bool SmoothInitialOutlinesTopologyPreserving(
    const std::string& shpPath,
    double pixelSizeX,
    double pixelSizeY,
    InitialOutlineSmoothStats& stats)
{
    stats = {};
    GDALDataset* dataset = static_cast<GDALDataset*>(
        GDALOpenEx(shpPath.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
                   nullptr, nullptr, nullptr));
    if (!dataset) return false;
    OGRLayer* layer = dataset->GetLayer(0);
    if (!layer) {
        GDALClose(dataset);
        return false;
    }

    const double tolerance = std::clamp(
        1.0 * std::max(pixelSizeX, pixelSizeY),
        kInitialOutlineSmoothMinTolerance,
        kInitialOutlineSmoothMaxTolerance);

    struct SmoothRec {
        GIntBig fid = OGRNullFID;
        std::unique_ptr<OGRGeometry> original;
        std::unique_ptr<OGRGeometry> smoothed;
        OGREnvelope env = {};
        bool changed = false;
        bool reverted = false;
    };
    std::vector<SmoothRec> recs;
    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        SmoothRec rec;
        rec.fid = feature->GetFID();
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry && !geometry->IsEmpty()) {
            rec.original.reset(geometry->clone());
            rec.original->getEnvelope(&rec.env);
        }
        if (rec.original) recs.push_back(std::move(rec));
        OGRFeature::DestroyFeature(feature);
    }
    if (recs.empty()) {
        GDALClose(dataset);
        return false;
    }

    // 第一遍：逐要素 SimplifyPreserveTopology + 安全检查
    for (auto& rec : recs) {
        ++stats.inspected;
        stats.beforeVertices += CountTotalVertices(rec.original.get());
        if (!rec.original) continue;

        const OGRwkbGeometryType origType = wkbFlatten(rec.original->getGeometryType());
        const int origHoles = CountInteriorRings(rec.original.get());
        const double origArea = GeometryArea(rec.original.get());
        const int origVerts = CountTotalVertices(rec.original.get());

        std::unique_ptr<OGRGeometry> simplified(rec.original->SimplifyPreserveTopology(tolerance));
        if (!simplified || simplified->IsEmpty()) {
            ++stats.revertedInvalid;
            continue;
        }

        const OGRwkbGeometryType simpType = wkbFlatten(simplified->getGeometryType());
        if (simpType != origType) {
            ++stats.revertedType;
            continue;
        }

        const int simpHoles = CountInteriorRings(simplified.get());
        if (simpHoles < origHoles) {
            ++stats.revertedInvalid;
            continue;
        }

        const double simpArea = GeometryArea(simplified.get());
        if (std::abs(simpArea - origArea) >
            kInitialOutlineSmoothAreaChangeRatio * std::max(origArea, 1.0)) {
            ++stats.revertedArea;
            continue;
        }

        const double deviation = MaxBoundaryDeviation(
            rec.original.get(), simplified.get(), tolerance);
        if (deviation > kInitialOutlineSmoothMaxDeviation) {
            ++stats.revertedDeviation;
            continue;
        }

        const int simpVerts = CountTotalVertices(simplified.get());
        if (simpVerts >= origVerts) {
            ++stats.revertedNoReduction;
            continue;
        }

        if (!simplified->IsValid()) {
            // 不用 Buffer(0) 修复，直接回退
            ++stats.revertedInvalid;
            continue;
        }

        rec.smoothed = std::move(simplified);
        rec.changed = true;
        ++stats.changed;
        if (simpHoles > 0) ++stats.holesPreserved;
    }

    // 第二遍：跨建筑邻接检查——平滑后不产生新重叠或异常缝隙
    for (std::size_t i = 0; i < recs.size(); ++i) {
        auto& a = recs[i];
        if (!a.changed) continue;
        for (std::size_t j = 0; j < recs.size(); ++j) {
            if (i == j) continue;
            auto& b = recs[j];
            if (!b.original) continue;
            // 平滑后几何的 envelope
            OGREnvelope envSmoothedA;
            if (a.smoothed) a.smoothed->getEnvelope(&envSmoothedA);
            else envSmoothedA = a.env;
            // 原始几何的 envelope(用于判断原始是否相邻)
            if (envSmoothedA.MaxX < b.env.MinX - kInitialOutlineSmoothAdjacencyTolerance ||
                b.env.MaxX < envSmoothedA.MinX - kInitialOutlineSmoothAdjacencyTolerance ||
                envSmoothedA.MaxY < b.env.MinY - kInitialOutlineSmoothAdjacencyTolerance ||
                b.env.MaxY < envSmoothedA.MinY - kInitialOutlineSmoothAdjacencyTolerance) {
                continue;
            }

            // 计算原始对是否有重叠
            const bool hadOverlap = a.original->Intersects(b.original.get());
            // 平滑后是否有重叠
            const bool nowOverlap = a.smoothed->Intersects(b.original.get());
            if (!hadOverlap && nowOverlap) {
                // 新增重叠：检查面积
                std::unique_ptr<OGRGeometry> overlap(a.smoothed->Intersection(b.original.get()));
                const double overlapArea = GeometryArea(overlap.get());
                if (overlapArea > kInitialOutlineSmoothOverlapArea) {
                    a.changed = false;
                    a.smoothed.reset();
                    ++stats.revertedAdjacency;
                    break;
                }
            }
        }
    }

    // 写回：只更新 changed 且未回退的要素
    for (const auto& rec : recs) {
        if (!rec.changed || !rec.smoothed || rec.reverted) continue;
        OGRFeature* feature = layer->GetFeature(rec.fid);
        if (!feature) continue;
        feature->SetGeometry(rec.smoothed.get());
        layer->SetFeature(feature);
        OGRFeature::DestroyFeature(feature);
        stats.afterVertices += CountTotalVertices(rec.smoothed.get());
    }
    // 未平滑的要素也计入 after 统计
    for (const auto& rec : recs) {
        if (rec.changed && rec.smoothed && !rec.reverted) continue;
        stats.afterVertices += CountTotalVertices(rec.original.get());
    }

    layer->SyncToDisk();
    GDALClose(dataset);

    std::cout << "[Mask smooth] inspected=" << stats.inspected
              << " changed=" << stats.changed
              << " before_vertices=" << stats.beforeVertices
              << " after_vertices=" << stats.afterVertices
              << " reverted_invalid=" << stats.revertedInvalid
              << " reverted_type=" << stats.revertedType
              << " reverted_area=" << stats.revertedArea
              << " reverted_deviation=" << stats.revertedDeviation
              << " reverted_adjacency=" << stats.revertedAdjacency
              << " reverted_no_reduction=" << stats.revertedNoReduction
              << " holes_preserved=" << stats.holesPreserved << std::endl;
    return true;
}

// 给 Shapefile 的每个要素新增或刷新 "area" 属性(平方米)。
// 掩膜各阶段会重建要素，中途写的面积会过期，
// 故在全部阶段结束后统一盖章。
bool StampAreaField(const std::string& shpPath)
{
    GDALDataset* dataset = static_cast<GDALDataset*>(
        GDALOpenEx(shpPath.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
                   nullptr, nullptr, nullptr));
    if (!dataset) return false;
    OGRLayer* layer = dataset->GetLayer(0);
    if (!layer) {
        GDALClose(dataset);
        return false;
    }
    OGRFeatureDefn* defn = layer->GetLayerDefn();
    int areaIdx = defn->GetFieldIndex("area");
    if (areaIdx < 0) {
        OGRFieldDefn areaField("area", OFTReal);
        if (layer->CreateField(&areaField) != OGRERR_NONE) {
            GDALClose(dataset);
            return false;
        }
        areaIdx = layer->GetLayerDefn()->GetFieldIndex("area");
    }
    long long updated = 0;
    layer->ResetReading();
    while (OGRFeature* feature = layer->GetNextFeature()) {
        OGRGeometry* geometry = feature->GetGeometryRef();
        if (geometry && !geometry->IsEmpty()) {
            feature->SetField(areaIdx, GeometryArea(geometry));
            if (layer->SetFeature(feature) == OGRERR_NONE) ++updated;
        }
        OGRFeature::DestroyFeature(feature);
    }
    layer->SyncToDisk();
    GDALClose(dataset);
    return updated > 0;
}

// 作用：把每栋的最优假设写出为 debug_best_hypothesis.shp。
bool SaveDebugBestHypotheses(
    const std::filesystem::path& outputPath,
    const RegularizationDebugCollector& debug,
    OGRSpatialReference* spatialRef,
    const Eigen::Vector3d& metadataOffset)
{
    // 删除旧文件，避免本次没有有效记录时残留上一次调试结果。
if (!RemoveShapefileFamily(outputPath)) return false;
    if (debug.hypotheses.empty()) return true;
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    if (!driver) return false;

    GDALDataset* dataset = driver->Create(outputPath.string().c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!dataset) return false;
    OGRLayer* layer = dataset->CreateLayer("debug_best_hypothesis", spatialRef, wkbPolygon25D, nullptr);
    if (!layer) {
        GDALClose(dataset);
        return false;
    }
    OGRFieldDefn srcField("src_fid", OFTInteger64);
    OGRFieldDefn idField("id", OFTInteger64);
    OGRFieldDefn ringField("ring_id", OFTInteger);
    OGRFieldDefn ptsField("npoints", OFTInteger);
    OGRFieldDefn areaField("area", OFTReal);
    layer->CreateField(&srcField);
    layer->CreateField(&idField);
    layer->CreateField(&ringField);
    layer->CreateField(&ptsField);
    layer->CreateField(&areaField);
    OGRFeatureDefn* defn = layer->GetLayerDefn();

    for (const auto& record : debug.hypotheses) {
        if (record.points.size() < 3) continue;
        std::unique_ptr<OGRPolygon> polygon(MakePolygon(record.points, metadataOffset));
        if (!polygon) continue;
        OGRFeature* feature = OGRFeature::CreateFeature(defn);
        feature->SetField("src_fid", static_cast<GIntBig>(record.sourceFid));
        // 统一的稳定建筑 id(可与 initial_building_outline 和
        // 规则化输出关联；栅格模式下与 src_fid 相同)。
        feature->SetField("id", static_cast<GIntBig>(record.sourceFid));
        feature->SetField("ring_id", record.ringIndex);
        feature->SetField("npoints", static_cast<int>(record.points.size()));
        feature->SetField("area", polygon->get_Area());
        feature->SetGeometry(polygon.get());
        layer->CreateFeature(feature);
        OGRFeature::DestroyFeature(feature);
    }
    GDALClose(dataset);
    return true;
}

// 作用：把调试支撑点云保存为 LAS。
bool SaveDebugSupportLas(
    const std::filesystem::path& outputPath,
    const RegularizationDebugCollector& debug,
    const Eigen::Vector3d& metadataOffset)
{
    std::error_code ec;
    std::filesystem::remove(outputPath, ec);
    if (ec) {
        std::cerr << "[LAS] Cannot remove previous debug file "
                  << outputPath.string() << ": " << ec.message() << std::endl;
        return false;
    }
    if (!debug.support || debug.support->empty()) return true;
    MyCloud cloud;
    cloud.cloud = debug.support;
    cloud.offset = metadataOffset;
    cloud.hasoffset = true;
    return SaveSampledCloudAsLas(cloud, outputPath.string());
}

// ============================================================================
// Mask-only 模式：完全不依赖 OSGB 的规则化。
// TIF -> 矢量化初始轮廓 -> 等弧长轮廓残差
// 点 -> 复用同一个 outlineRegular::regular_Contour 管线(VDP、
// 以轮廓弦方向直方图为证据的方向判定、Ceres) -> 输出。
// ============================================================================


// 沿初始轮廓的均匀残差点间距(米)。等弧长
// 均匀布点消除栅格轮廓的楼梯顶点密度偏置；
// 点间距 1-2.5m 的弦方向同时充当
// 去楼梯化的边方向直方图。
const double kMaskResidualSpacing = 0.5;

// 作用：判断 Mask-only 结果环是否自交(线段两两相交测试)。
bool MaskOnlyRingIsSimple(const std::vector<pcl::PointXYZ>& ring)
{
    const std::size_t n = ring.size();
    if (n < 4) return true;
    std::vector<Point2D64> pts;
    pts.reserve(n);
    for (const auto& p : ring) pts.push_back({p.x, p.y, p.z});
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 2; j < n; ++j) {
            if (i == 0 && j == n - 1) continue;
            if (SegmentsIntersect2D64(pts[i], pts[(i + 1) % n],
                                      pts[j], pts[(j + 1) % n])) {
                return false;
            }
        }
    }
    return true;
}

enum class MaskOnlyFallback { Final = 0, Hypothesis, Initial };
// 规则化最终路径(统计与日志用)
enum class MaskOnlyPath {
    Topology = 0,
    Vdp,
    BestHypothesis, // diagnostic-only legacy value; never emitted
    InitialRing,    // diagnostic-only legacy value; never emitted
    StrictFallback
};

// Mask-only 规则化入口：轮廓点既是残差点云，
// 也是方向证据，供原封不动的 regular_Contour 管线使用。
std::vector<pcl::PointXYZ> RegularizeRingFromMaskOnly(
    const std::vector<pcl::PointXYZ>& ring,
    long long sourceFid,
    std::vector<pcl::PointXYZ>* bestHypothesisOut,
    MaskOnlyFallback* fallbackLevel,
    MaskOnlyPath* pathOut = nullptr,
    int partIndex = 0,
    const std::vector<pcl::PointXYZ>* rawRing = nullptr,
    outlineRegular::DirectionContextOut* dirCtxOut = nullptr,
    double maskPixelSize = 0.5)
{
    if (fallbackLevel) *fallbackLevel = MaskOnlyFallback::Final;
    if (pathOut) *pathOut = MaskOnlyPath::InitialRing;
    if (ring.size() < 3) {
        std::cerr << "[RegularizationPath] fid=" << sourceFid
                  << " failed_too_few_input_vertices" << std::endl;
        return {};
    }

    auto residual = DensifyBoundary(ring, kMaskResidualSpacing);
    pcl::PointCloud<pcl::PointXYZ>::Ptr contourCloud(new pcl::PointCloud<pcl::PointXYZ>);
    contourCloud->points.assign(residual->begin(), residual->end());
    if (contourCloud->size() < 8) {
        const auto strict = OrientedBoundingRectangle(ring, 0.0);
        if (strict.size() >= 3) {
            if (pathOut) *pathOut = MaskOnlyPath::StrictFallback;
            std::cerr << "[RegularizationPath] fid=" << sourceFid
                      << " strict_direction_fallback angle_deg=0 vertices="
                      << strict.size() << std::endl;
            return strict;
        }
        std::cerr << "[RegularizationPath] fid=" << sourceFid
                  << " failed_insufficient_residual_points" << std::endl;
        return {};
    }

    double contourDir = 0.0;
    double contourRatio = 0.0;
    std::size_t contourPairs = 0;
    outlineRegular::estimateSupportDirection2D(
        contourCloud, contourDir, contourRatio, contourPairs);

    std::cerr << "[MaskOnlySupport] fid=" << sourceFid
              << " original_vertices=" << ring.size()
              << " resampled_points=" << contourCloud->size()
              << " spacing=" << kMaskResidualSpacing
              << " weights=contour_uniform" << std::endl;
    std::cerr << "[MaskOnlyDirection] fid=" << sourceFid
              << " contour_hist_deg=" << contourDir * 180.0 / M_PI
              << " peak_ratio=" << contourRatio
              << " pairs=" << contourPairs << std::endl;

    // 方向诊断已并入 TopologyPreservingRegularize(单次调用统一用于
    // 决策/吸附/上下文/日志), 不再单独执行不同输入的第二遍诊断

    std::vector<double> weights(contourCloud->size(), 1.0);
    outlineRegular regularizer(ring, contourCloud, weights);
    regularizer.setSourceFeatureId(sourceFid);
    regularizer.setSupportDirectionHint(contourDir, contourRatio, contourPairs);
    // 拓扑通道的方向结论(传给 VDP 备用结果做方向一致性检查)
    outlineRegular::DirectionContextOut dirCtx;

    // ---- 拓扑保持通道(主流程) ----
    // 从边链出发而非VDP假设，保留真实凹凸/窄颈拓扑。
    // 合格即直接采用; 失败(链不足/方向不确定/候选不合格/Ceres 大位移/
    // 斜边/圆形轮廓)转 VDP 备用, 方向上下文随之传出。
    // ---- Mask-only 局部圆弧检测(所有路径共用) ----
    // 检测在平滑环上确定区间; rawRing 可用时提供拟合支撑;
    // 恢复在各路径的直线候选通过质量检查后执行
    const auto maskCurves = DetectMaskConicArcs(
        ring, rawRing ? *rawRing : std::vector<pcl::PointXYZ>{},
        maskPixelSize > 0.0 ? maskPixelSize : 0.3, sourceFid, partIndex);

        if (kUseTopologyPreservingResidualRegularization) {
        bool topoFallback = false;
        auto topoResult = regularizer.TopologyPreservingRegularize(
            ring, kMaskResidualSpacing, topoFallback, &dirCtx, partIndex,
            rawRing);
        if (!topoResult.empty()) {
            if (bestHypothesisOut) *bestHypothesisOut = topoResult;
            if (pathOut) *pathOut = MaskOnlyPath::Topology;
            std::cerr << "[RegularizationPath] fid=" << sourceFid
                      << " topology" << std::endl;
            std::cerr << "[MaskOnlyTopology] fid=" << sourceFid
                      << " vertices=" << topoResult.size() << std::endl;
            // 曲线恢复: 直线候选已通过质量检查, 回贴圆弧采样
            if (!kMaskCurveDetectionDebugOnly) topoResult = RestoreMaskConicArcs(
                topoResult, maskCurves, kMaskResidualSpacing, sourceFid, partIndex);
            return topoResult;
        }
        if (topoFallback) {
            std::cerr << "[MaskOnlyTopology] fid=" << sourceFid
                      << " fallback to VDP pipeline" << std::endl;
        }
    }
    // 方向上下文传出(通道内已填好, VDP 路径不再改动)
    if (dirCtxOut) *dirCtxOut = dirCtx;

    // 无正射证据仲裁栅格楼梯上的曲线检测；
    // Mask-only 不得恢复伪曲线。
    regularizer.setCurveRestorationEnabled(false);
    regularizer.regular_Contour();

    std::vector<pcl::PointXYZ> bestHypothesis;
    if (bestHypothesisOut) {
        *bestHypothesisOut = regularizer.getBestEnergyHypothesis();
        if (bestHypothesisOut->size() < 3) *bestHypothesisOut = ring;
        RemoveClosingDuplicate(*bestHypothesisOut);
        bestHypothesis = *bestHypothesisOut;
    } else {
        bestHypothesis = regularizer.getBestEnergyHypothesis();
        if (bestHypothesis.size() < 3) bestHypothesis = ring;
        RemoveClosingDuplicate(bestHypothesis);
    }

    std::vector<pcl::PointXYZ> result;
    if (regularizer.final_points && regularizer.final_points->size() >= 3) {
        result.assign(regularizer.final_points->points.begin(),
                      regularizer.final_points->points.end());
        RemoveClosingDuplicate(result);
    }

    // ---- VDP 备用结果质量闸门: 与拓扑结果同一标准 ----
    // 不合格的 VDP 结果不能覆盖合格的拓扑结果(拓扑合格已在上方直接
    // 返回), 也不能直接进入输出; VDP 结果额外过兜底局部规则性检查
    // (短斜边/尖刺/锯齿); 逐级退到最优假设(保守清理后重验)、初始轮廓。
    if (result.size() >= 3) {
        const std::string vdpReason =
            outlineRegular::CheckRingQuality(result, ring, dirCtx, 2.5, sourceFid, partIndex);
        const std::string vdpLocal = vdpReason.empty()
            ? outlineRegular::CheckFallbackLocalRegularity(
                  result, dirCtx, sourceFid, partIndex, "vdp")
            : std::string();
        if (!vdpReason.empty() || !vdpLocal.empty()) {
            std::cerr << "[VDPReject] fid=" << sourceFid
                      << " stage=final reason="
                      << (!vdpReason.empty() ? vdpReason : vdpLocal) << std::endl;
            result.clear();
        }
    }
    if (result.size() >= 3) {
        if (pathOut) *pathOut = MaskOnlyPath::Vdp;
        std::cerr << "[RegularizationPath] fid=" << sourceFid << " vdp" << std::endl;
        std::cerr << "[MaskOnlyHypothesis] fid=" << sourceFid
                  << " initial_vertices=" << ring.size()
                  << " final_vertices=" << result.size()
                  << " initial_area=" << PolygonArea2D(ring)
                  << " final_area=" << PolygonArea2D(result) << std::endl;
        // 曲线恢复: VDP直线候选通过质量检查后回贴
        if (!kMaskCurveDetectionDebugOnly) result = RestoreMaskConicArcs(
            result, maskCurves, kMaskResidualSpacing, sourceFid, partIndex);
        return result;
    }
    // The VDP hypothesis is diagnostic data and a source for the strict
    // fallback only. It must never be emitted as a final mask-only building.
    std::vector<pcl::PointXYZ> strictInput = bestHypothesis;
    if (bestHypothesis.size() >= 3) {
        const std::string hypReason =
            outlineRegular::CheckRingQuality(bestHypothesis, ring, dirCtx, 2.5, sourceFid, partIndex);
        const std::string hypLocal = hypReason.empty()
            ? outlineRegular::CheckFallbackLocalRegularity(
                  bestHypothesis, dirCtx, sourceFid, partIndex, "best")
            : std::string();
        if (!hypReason.empty() || !hypLocal.empty()) {
            // 最优假设不合格: 保守清理低面积尖刺/近共线锯齿后重验
            // (只删确定性噪声, 真实 90° 凹凸不受影响)
            auto cleaned = outlineRegular::CleanLowEvidenceIrregularities(bestHypothesis);
            const std::string cleanReason =
                outlineRegular::CheckRingQuality(cleaned, ring, dirCtx, 2.5, sourceFid, partIndex);
            const std::string cleanLocal = cleanReason.empty()
                ? outlineRegular::CheckFallbackLocalRegularity(
                      cleaned, dirCtx, sourceFid, partIndex, "best_cleaned")
                : std::string();
            if (!cleanReason.empty() || !cleanLocal.empty() || cleaned.size() < 3) {
                std::cerr << "[VDPReject] fid=" << sourceFid
                          << " stage=best_hypothesis reason="
                          << (!hypReason.empty() ? hypReason : hypLocal)
                          << " (cleaned 也未通过: "
                          << (!cleanReason.empty() ? cleanReason : cleanLocal)
                          << ")" << std::endl;
            } else {
                strictInput = cleaned;
                std::cerr << "[VDPUsableForStrictFallback] fid=" << sourceFid
                          << " cleaned "
                          << bestHypothesis.size() << "->" << cleaned.size()
                          << " verts)" << std::endl;
            }
        } else {
            std::cerr << "[VDPUsableForStrictFallback] fid=" << sourceFid
                      << " best_vertices=" << bestHypothesis.size()
                      << " (direct emission disabled)" << std::endl;
        }
    }

    // Never write the noisy initial ring as a final regularized building.
    // Use the strongest known direction to build one last strict candidate;
    // when no direction evidence exists, an axis-aligned OBR is still a
    // regularized result and is preferable to exposing pixel stair steps.
    const double fallbackAngle = dirCtx.valid ? dirCtx.primaryAngle : 0.0;
    std::vector<pcl::PointXYZ> strictFallback;
    if (strictInput.size() >= 3) {
        regularizer.BuildStrictDirectionalFallback(
            strictInput, fallbackAngle, strictFallback);
    }
    if (strictFallback.size() < 3) {
        strictFallback = OrientedBoundingRectangle(ring, fallbackAngle);
    }
    RemoveClosingDuplicate(strictFallback);
    if (strictFallback.size() >= 3) {
        if (pathOut) *pathOut = MaskOnlyPath::StrictFallback;
        std::cerr << "[RegularizationPath] fid=" << sourceFid
                  << " strict_direction_fallback angle_deg="
                  << fallbackAngle * 180.0 / M_PI
                  << " vertices=" << strictFallback.size() << std::endl;
        // 曲线恢复: 严格方向候选通过质量检查后回贴
        if (!kMaskCurveDetectionDebugOnly) strictFallback = RestoreMaskConicArcs(
            strictFallback, maskCurves, kMaskResidualSpacing, sourceFid, partIndex);
        return strictFallback;
    }

    // This is an actual processing failure. Do not silently claim that the
    // initial ring was regularized; the caller will omit the invalid result.
    if (fallbackLevel) *fallbackLevel = MaskOnlyFallback::Final;
    std::cerr << "[RegularizationPath] fid=" << sourceFid
              << " failed_no_regularized_candidate" << std::endl;
    return {};
}

// 作用：Mask-only 模式主入口：TIF→初始轮廓→等弧长残差点→复用规则化→输出。
int RunMaskOnlyMode(const std::string& inputRaster, const std::string& outputVectorIn)
{
    std::cout << "[Mode] Mask-only" << std::endl;
    std::cout << "[MaskOnly] input_tif=" << inputRaster << std::endl;
    std::cout << "[MaskOnly] output_shp=" << outputVectorIn << std::endl;
    GDALAllRegister();

    // 仿射变换 + 空间参考(必需：输出保留 TIF 的坐标系)。
    GDALDataset* tif = static_cast<GDALDataset*>(
        GDALOpenEx(inputRaster.c_str(), GDAL_OF_RASTER, nullptr, nullptr, nullptr));
    if (!tif) {
        std::cerr << "[MaskOnly] cannot open input raster: " << inputRaster << std::endl;
        return 1;
    }
    double gt[6] = {0};
    if (tif->GetGeoTransform(gt) != CE_None) {
        std::cerr << "[MaskOnly] input raster has no GeoTransform." << std::endl;
        GDALClose(tif);
        return 1;
    }
    const char* projection = tif->GetProjectionRef();
    if (!projection || projection[0] == '\0') {
        std::cerr << "[MaskOnly] input raster has no spatial reference." << std::endl;
        GDALClose(tif);
        return 1;
    }
    const double originX = std::round(gt[0] + gt[1] * tif->GetRasterXSize() * 0.5);
    const double originY = std::round(gt[3] + gt[5] * tif->GetRasterYSize() * 0.5);
    const Eigen::Vector3d originOffset(originX, originY, 0.0);
    std::cout << "[MaskOnly] origin_x=" << originX << " origin_y=" << originY << std::endl;
    std::cout << "[MaskOnly] geotransform=" << gt[0] << "," << gt[1] << "," << gt[2]
              << "," << gt[3] << "," << gt[4] << "," << gt[5] << std::endl;
    GDALClose(tif);

    // 输出准备：调试文件放在输出 Shapefile 旁边。
    std::filesystem::path outPath(outputVectorIn);
    std::error_code dirEc;
    std::filesystem::create_directories(outPath.parent_path(), dirEc);
    if (dirEc && !outPath.parent_path().empty()) {
        std::cerr << "[MaskOnly] cannot create output directory: "
                  << outPath.parent_path().string() << " (" << dirEc.message() << ")" << std::endl;
        return 1;
    }
    const std::string outputVector =
        PrepareWritableShapefilePath(outPath, "output Shapefile").string();
    const std::filesystem::path debugDir =
        std::filesystem::path(outputVector).parent_path();

    std::filesystem::path initialPath = debugDir / "initial_building_outline.shp";
    if (std::filesystem::exists(initialPath) && !RemoveShapefileFamily(initialPath, false)) {
        initialPath = MakeUniqueShapefilePath(initialPath);
        std::cout << "[SHP] initial outline path locked, writing to: "
                  << initialPath.string() << std::endl;
    }

    // ---- Initial outlines: full mask geometry chain ----
    MaskVectorizationStats maskStats;
    if (!VectorizeBuildingMask(inputRaster, initialPath.string(), maskStats)) {
        std::cerr << "[MaskOnly] mask vectorization failed (no valid building pixels?)."
                  << std::endl;
        return 1;
    }

    // ---- Stage 1: GeometryOnly 合并(同连通域 + 共享边>=3m) ----
    long long totalMergePasses = 0;
    long long totalMergedPairs = 0;
    for (int pass = 1; pass <= 20; ++pass) {
        InitialOutlineMergeStats passStats;
        const bool mergeOk = MergeOversegmentedInitialOutlines(
            initialPath.string(), nullptr, nullptr, Eigen::Vector3d::Zero(),
            passStats, InitialMergeMode::GeometryOnly);
        if (!mergeOk && pass == 1) {
            std::cerr << "[MaskOnly merge] failed to open initial outlines." << std::endl;
            break;
        }
        ++totalMergePasses;
        totalMergedPairs += passStats.mergedPairs;
        std::cout << "[MaskOnly merge] pass " << pass
                  << " candidates=" << passStats.adjacentCandidates
                  << " same_parent=" << passStats.mergedByParent
                  << " merged=" << passStats.mergedPairs
                  << " removed=" << passStats.removedFeatures << std::endl;
        if (passStats.mergedPairs == 0) break;
        if (passStats.removedFeatures == 0 && pass > 1) break;
    }

    // ---- Stage 2: 包含清理(Mask-only 用 20m² 阈值，保留孔洞) ----
    const long long containedRemoved1 =
        RemoveContainedSmallFootprints(initialPath.string(),
            kNarrowNeckMinPartArea, /*solidifyHoles=*/false);
    std::cout << "[MaskOnly contained] pass1 removed=" << containedRemoved1 << std::endl;

    // ---- Stage 3: 窄颈递归切分 ----
    NarrowNeckSplitStats neckStats;
    if (!SplitInitialOutlinesAtNarrowNecks(initialPath.string(), neckStats)) {
        std::cerr << "[MaskOnly] neck split failed." << std::endl;
    } else {
        std::cout << "[MaskOnly neck split] split=" << neckStats.splitFeatures
                  << " cuts=" << neckStats.cuts
                  << " parts=" << neckStats.createdParts << std::endl;
    }

    // ---- Stage 4: 合并后二次包含清理 ----
    const long long containedRemoved2 =
        RemoveContainedSmallFootprints(initialPath.string(),
            kNarrowNeckMinPartArea, /*solidifyHoles=*/false);
    std::cout << "[MaskOnly contained] pass2 removed=" << containedRemoved2 << std::endl;

    // 保拓扑平滑：在窄颈拆分和包含清理之后、StampAreaField 之前执行。
    // 先保存原始副本用于对比诊断。
    {
        const std::filesystem::path rawPath = debugDir / "initial_building_outline_raw.shp";
        if (std::filesystem::exists(initialPath)) {
            // 复制整个 Shapefile 族
            for (const char* ext : {".shp", ".shx", ".dbf", ".prj", ".cpg"}) {
                std::error_code copyEc;
                std::filesystem::path src = initialPath;
                src.replace_extension(ext);
                std::filesystem::path dst = rawPath;
                dst.replace_extension(ext);
                if (std::filesystem::exists(src)) {
                    std::filesystem::copy_file(src, dst,
                        std::filesystem::copy_options::overwrite_existing, copyEc);
                }
            }
            std::cout << "[MaskOnly] raw initial saved: " << rawPath.string() << std::endl;
        }
        InitialOutlineSmoothStats smoothStats;
        SmoothInitialOutlinesTopologyPreserving(
            initialPath.string(), maskStats.pixelSizeX, maskStats.pixelSizeY, smoothStats);
    }

    StampAreaField(initialPath.string());
    std::cout << "[MaskOnly] initial outlines saved: " << initialPath.string() << std::endl;

    // ---- Feature loop (no OSGB: no coverage filter, no ownership) ----
    GDALDataset* inDataset = static_cast<GDALDataset*>(
        GDALOpenEx(initialPath.string().c_str(), GDAL_OF_VECTOR,
                   nullptr, nullptr, nullptr));
    if (!inDataset) {
        std::cerr << "[MaskOnly] cannot reopen initial outlines." << std::endl;
        return 1;
    }
    OGRLayer* inLayer = inDataset->GetLayer(0);
    if (!inLayer) {
        std::cerr << "[MaskOnly] initial outline layer missing." << std::endl;
        GDALClose(inDataset);
        return 1;
    }
    // 关闭数据集前克隆 SRS：调试假设写出器
    // 在 inDataset 销毁后运行，仍需携带 TIF 参考系。
    OGRSpatialReference* srsClone = inLayer->GetSpatialRef()
        ? inLayer->GetSpatialRef()->Clone()
        : nullptr;
    const int initialCount = inLayer->GetFeatureCount(TRUE);
    std::cout << "[MaskOnly] initial feature count=" << initialCount << std::endl;

    // ---- 原始(平滑前)轮廓加载: 按 FID 精确关联, 供双残差规则化 ----
    // raw 是平滑前的同一图层副本，SmoothInitialOutlinesTopologyPreserving
    // 只用 SetFeature 写回，因此 FID 在平滑前后保持稳定；id 可能因窄颈
    // 切分而重复，不能作为几何残差的唯一关联键。
    std::unordered_map<long long, std::vector<std::vector<pcl::PointXYZ>>> rawRingsByFid;
    {
        const std::filesystem::path rawPath = debugDir / "initial_building_outline_raw.shp";
        if (std::filesystem::exists(rawPath)) {
            GDALDataset* rawDataset = static_cast<GDALDataset*>(
                GDALOpenEx(rawPath.string().c_str(), GDAL_OF_VECTOR,
                           nullptr, nullptr, nullptr));
            if (rawDataset) {
                OGRLayer* rawLayer = rawDataset->GetLayer(0);
                if (rawLayer) {
                    rawLayer->ResetReading();
                    long long rawRingCount = 0;
                    while (OGRFeature* feature = rawLayer->GetNextFeature()) {
                        OGRGeometry* geometry = feature->GetGeometryRef();
                        const long long key = feature->GetFID();
                        if (geometry && !geometry->IsEmpty()) {
                            const OGRwkbGeometryType rawType =
                                wkbFlatten(geometry->getGeometryType());
                            std::vector<OGRPolygon*> rawParts;
                            if (rawType == wkbPolygon) {
                                rawParts.push_back(geometry->toPolygon());
                            } else if (rawType == wkbMultiPolygon ||
                                       rawType == wkbGeometryCollection) {
                                auto* collection = geometry->toGeometryCollection();
                                for (int i = 0; collection && i < collection->getNumGeometries(); ++i) {
                                    if (wkbFlatten(collection->getGeometryRef(i)->getGeometryType()) == wkbPolygon) {
                                        rawParts.push_back(collection->getGeometryRef(i)->toPolygon());
                                    }
                                }
                            }
                            for (OGRPolygon* part : rawParts) {
                                auto rawRing = ExtractExteriorRing(part, originOffset);
                                RemoveClosingDuplicate(rawRing);
                                if (rawRing.size() >= 3) {
                                    rawRingsByFid[key].push_back(std::move(rawRing));
                                    ++rawRingCount;
                                }
                            }
                        }
                        OGRFeature::DestroyFeature(feature);
                    }
                    std::cout << "[MaskOnly] raw rings loaded=" << rawRingCount
                              << " features=" << rawRingsByFid.size() << std::endl;
                }
                GDALClose(rawDataset);
            }
        } else {
            std::cout << "[MaskOnly] raw outline missing; dual residual disabled"
                      << std::endl;
        }
    }

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    if (!driver) {
        std::cerr << "[MaskOnly] ESRI Shapefile driver unavailable." << std::endl;
        GDALClose(inDataset);
        return 1;
    }
    GDALDataset* outDataset = driver->Create(
        outputVector.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!outDataset) {
        std::cerr << "[MaskOnly] cannot create output vector: " << outputVector << std::endl;
        GDALClose(inDataset);
        return 1;
    }
    OGRLayer* outLayer = outDataset->CreateLayer(
        "regularized_building", inLayer->GetSpatialRef(), wkbPolygon25D, nullptr);
    if (!outLayer) {
        std::cerr << "[MaskOnly] cannot create output layer." << std::endl;
        GDALClose(outDataset);
        GDALClose(inDataset);
        return 1;
    }
    CopyFields(inLayer, outLayer);
    const int inIdFieldIdx = inLayer->GetLayerDefn()->GetFieldIndex("id");
    int outIdFieldIdx = outLayer->GetLayerDefn()->GetFieldIndex("id");
    if (outIdFieldIdx < 0) {
        OGRFieldDefn idField("id", OFTInteger64);
        if (outLayer->CreateField(&idField) == OGRERR_NONE) {
            outIdFieldIdx = outLayer->GetLayerDefn()->GetFieldIndex("id");
        }
    }
    int outAreaFieldIdx = outLayer->GetLayerDefn()->GetFieldIndex("area");
    if (outAreaFieldIdx < 0) {
        OGRFieldDefn areaField("area", OFTReal);
        if (outLayer->CreateField(&areaField) == OGRERR_NONE) {
            outAreaFieldIdx = outLayer->GetLayerDefn()->GetFieldIndex("area");
        }
    }
    OGRFeatureDefn* outDefn = outLayer->GetLayerDefn();

    RegularizationDebugCollector debugCollector;
    // Legacy best/initial slots remain for backward-compatible enum values;
    // current mask-only output uses topology, vdp, or strict_fallback only.
    struct PathStatAccum {
        long long count = 0;
        double sumVerts = 0.0;
        long long sumShortEdges = 0;
        double sumAreaRatio = 0.0;
    };
    std::vector<PathStatAccum> pathStats(5);
    // 重叠解决器输入: 建筑优先级与"更好的单体候选"(初始轮廓, 局部坐标)
    const double maskPixelSize = std::max(maskStats.pixelSizeX, maskStats.pixelSizeY);
    std::unordered_map<long long, double> overlapPriorityByFid;
    std::unordered_map<long long, std::vector<pcl::PointXYZ>> overlapAlternateByFid;
    // 输出 FID → 该建筑的方向上下文(Difference 后方向检查用)
    std::unordered_map<long long, outlineRegular::DirectionContextOut> overlapDirectionByFid;
    long long total = 0;
    long long okCount = 0;
    long long emptyGeom = 0;
    long long smallSkipped = 0;
    long long holesDropped = 0;
    long long selfIntersect = 0;
    long long fallbackHypothesis = 0;
    long long fallbackInitial = 0;
    double sumInitVerts = 0.0;
    double sumBestVerts = 0.0;
    double sumFinalVerts = 0.0;
    double sumAreaRatio = 0.0;
    long long areaSamples = 0;
    auto loopStart = std::chrono::steady_clock::now();

    inLayer->ResetReading();
    while (OGRFeature* inFeature = inLayer->GetNextFeature()) {
        ++total;
        OGRGeometry* geometry = inFeature->GetGeometryRef();
        if (!geometry || geometry->IsEmpty()) {
            std::cerr << "[MaskOnly] skip empty geometry fid=" << inFeature->GetFID() << std::endl;
            ++emptyGeom;
            OGRFeature::DestroyFeature(inFeature);
            continue;
        }
        GIntBig buildingId = inFeature->GetFID();
        if (inIdFieldIdx >= 0 && inFeature->IsFieldSetAndNotNull(inIdFieldIdx)) {
            buildingId = inFeature->GetFieldAsInteger64(inIdFieldIdx);
        }

        // Collect exterior rings of every polygon part.
        std::vector<OGRPolygon*> parts;
        const OGRwkbGeometryType type = wkbFlatten(geometry->getGeometryType());
        if (type == wkbPolygon) {
            parts.push_back(geometry->toPolygon());
        } else if (type == wkbMultiPolygon || type == wkbGeometryCollection) {
            auto* collection = geometry->toGeometryCollection();
            for (int i = 0; collection && i < collection->getNumGeometries(); ++i) {
                const OGRwkbGeometryType partType =
                    wkbFlatten(collection->getGeometryRef(i)->getGeometryType());
                if (partType == wkbPolygon) {
                    parts.push_back(collection->getGeometryRef(i)->toPolygon());
                }
            }
        }

        // 要素级重叠处理输入累积: 按 part 最大优先级 + 对应方向上下文,
        // CreateFeature 成功后以真实输出 FID 写入映射(属性 id 与输出 FID
        // 因 smallSkipped 等不再一致, 禁止再用 buildingId 当键)
        double featurePriority = -1e9;
        MaskOnlyPath featurePath = MaskOnlyPath::InitialRing;
        outlineRegular::DirectionContextOut featureDirCtx;
        std::vector<std::unique_ptr<OGRPolygon>> outParts;
        int ringIdx = 0;
        for (OGRPolygon* part : parts) {
            if (part->getNumInteriorRings() > 0) {
                ++holesDropped;
            }
            auto ring = ExtractExteriorRing(part, originOffset);
            ++ringIdx;
            if (BoundingBoxArea2D(ring) < kMinPolygonBBoxArea) {
                ++smallSkipped;
                continue;
            }
            const std::vector<pcl::PointXYZ>* rawRingForPart = nullptr;
            const auto rawIt = rawRingsByFid.find(
                static_cast<long long>(inFeature->GetFID()));
            if (rawIt != rawRingsByFid.end() && !rawIt->second.empty()) {
                // SimplifyPreserveTopology normally preserves part order. For
                // defensive handling of reordered multipart geometries, choose
                // the raw part with the largest bbox overlap with this smooth part.
                double bestScore = -1.0;
                for (const auto& rawCandidate : rawIt->second) {
                    double sMinX = 1e18, sMinY = 1e18, sMaxX = -1e18, sMaxY = -1e18;
                    double rMinX = 1e18, rMinY = 1e18, rMaxX = -1e18, rMaxY = -1e18;
                    for (const auto& p : ring) {
                        sMinX = std::min(sMinX, static_cast<double>(p.x));
                        sMinY = std::min(sMinY, static_cast<double>(p.y));
                        sMaxX = std::max(sMaxX, static_cast<double>(p.x));
                        sMaxY = std::max(sMaxY, static_cast<double>(p.y));
                    }
                    for (const auto& p : rawCandidate) {
                        rMinX = std::min(rMinX, static_cast<double>(p.x));
                        rMinY = std::min(rMinY, static_cast<double>(p.y));
                        rMaxX = std::max(rMaxX, static_cast<double>(p.x));
                        rMaxY = std::max(rMaxY, static_cast<double>(p.y));
                    }
                    const double overlapW = std::max(
                        0.0, std::min(sMaxX, rMaxX) - std::max(sMinX, rMinX));
                    const double overlapH = std::max(
                        0.0, std::min(sMaxY, rMaxY) - std::max(sMinY, rMinY));
                    const double minBoxArea = std::max(1e-9, std::min(
                        (sMaxX - sMinX) * (sMaxY - sMinY),
                        (rMaxX - rMinX) * (rMaxY - rMinY)));
                    const double score = overlapW * overlapH / minBoxArea;
                    if (score > bestScore) {
                        bestScore = score;
                        rawRingForPart = &rawCandidate;
                    }
                }
            }
            outlineRegular::DirectionContextOut ringDirCtx;
            std::vector<pcl::PointXYZ> bestHypothesis;
            MaskOnlyFallback fallback = MaskOnlyFallback::Final;
            MaskOnlyPath path = MaskOnlyPath::InitialRing;
            auto result = RegularizeRingFromMaskOnly(
                ring, static_cast<long long>(buildingId), &bestHypothesis,
                &fallback, &path, ringIdx - 1,
                rawRingForPart, &ringDirCtx, maskPixelSize);
            if (fallback == MaskOnlyFallback::Hypothesis) ++fallbackHypothesis;
            if (fallback == MaskOnlyFallback::Initial) ++fallbackInitial;
            // 路径统计: 顶点数/短边数(<0.5m)/面积变化
            {
                auto& st = pathStats[static_cast<std::size_t>(path)];
                ++st.count;
                st.sumVerts += static_cast<double>(result.size());
                for (std::size_t i = 0; i < result.size(); ++i) {
                    const double len = std::hypot(
                        result[(i + 1) % result.size()].x - result[i].x,
                        result[(i + 1) % result.size()].y - result[i].y);
                    if (len < 0.5) ++st.sumShortEdges;
                }
                const double initA = std::abs(PolygonArea2D(ring));
                const double outA = std::abs(PolygonArea2D(result));
                if (initA > 1e-6) st.sumAreaRatio += outA / initA;
            }
            // 重叠解决器输入采集:
            // 优先级 = 路径置信度 + 面积规模 + 小建筑矩形语义加成
            {
                double prio = 0.5;
                switch (path) {
                    case MaskOnlyPath::Topology: prio += 0.35; break;
                    case MaskOnlyPath::Vdp: prio += 0.15; break;
                    case MaskOnlyPath::BestHypothesis: break;
                    case MaskOnlyPath::InitialRing: prio -= 0.25; break;
                    case MaskOnlyPath::StrictFallback: prio += 0.05; break;
                }
                prio += 0.20 * std::min(1.0, std::abs(PolygonArea2D(result)) / 200.0);
                // 小建筑矩形(≤6顶点的 VDP 矩形拟合): 提高优先级,
                // 相交时邻居让位, 矩形只被裁掉超出部分
                if (path == MaskOnlyPath::Vdp && result.size() <= 6) prio += 0.15;
                // 要素级累积(取 part 最大优先级; 方向上下文随最高优先级 part)
                if (prio > featurePriority) {
                    featurePriority = prio;
                    featurePath = path;
                    featureDirCtx = ringDirCtx;
                }
            }

            if (bestHypothesis.size() >= 3) {
                debugCollector.hypotheses.push_back(
                    {buildingId, ringIdx - 1, bestHypothesis});
            }
            sumInitVerts += static_cast<double>(ring.size());
            sumBestVerts += static_cast<double>(bestHypothesis.size());
            sumFinalVerts += static_cast<double>(result.size());
            const double initArea = PolygonArea2D(ring);
            if (initArea > 1e-6) {
                sumAreaRatio += PolygonArea2D(result) / initArea;
                ++areaSamples;
            }
            if (!MaskOnlyRingIsSimple(result)) ++selfIntersect;

            // 面积地板: <15m² 或 bbox<20m² 的残片不写出
            {
                const double partArea = std::abs(PolygonArea2D(result));
                double bbMinX=1e18,bbMinY=1e18,bbMaxX=-1e18,bbMaxY=-1e18;
                for (const auto& p : result) {
                    bbMinX=std::min(bbMinX,(double)p.x); bbMaxX=std::max(bbMaxX,(double)p.x);
                    bbMinY=std::min(bbMinY,(double)p.y); bbMaxY=std::max(bbMaxY,(double)p.y);
                }
                const double partBBox = (bbMaxX-bbMinX)*(bbMaxY-bbMinY);
                if (partArea < 15.0 || partBBox < 20.0) {
                    std::cerr << "[MaskOnlySizeFloor] fid=" << buildingId
                              << " area=" << partArea << " bbox=" << partBBox
                              << " dropped=1" << std::endl;
                    ++smallSkipped;
                    continue;
                }
            }
            std::unique_ptr<OGRPolygon> outPolygon(MakePolygon(result, originOffset));
            if (outPolygon) outParts.push_back(std::move(outPolygon));
        }

        if (outParts.empty()) {
            std::cerr << "[MaskOnly] no valid rings for fid=" << buildingId << std::endl;
            OGRFeature::DestroyFeature(inFeature);
            continue;
        }
        OGRFeature* outFeature = OGRFeature::CreateFeature(outDefn);
        CopyFieldValues(inFeature, outFeature);
        if (outIdFieldIdx >= 0) outFeature->SetField(outIdFieldIdx, buildingId);
        if (outParts.size() == 1) {
            outFeature->SetGeometry(outParts.front().get());
        } else {
            OGRMultiPolygon multi;
            for (const auto& part : outParts) multi.addGeometry(part.get());
            outFeature->SetGeometry(&multi);
        }
        if (outAreaFieldIdx >= 0) {
            outFeature->SetField(outAreaFieldIdx,
                GeometryArea(outFeature->GetGeometryRef()));
        }
        // CreateFeature 成功后以真实输出 FID 写入重叠处理映射
        // (属性 id 与输出 FID 因 smallSkipped 等不再一致)
        if (outLayer->CreateFeature(outFeature) == OGRERR_NONE) {
            ++okCount;
            if (featurePriority > -1e8) {
                const GIntBig outFid = outFeature->GetFID();
                overlapPriorityByFid[static_cast<long long>(outFid)] =
                    featurePriority;
                overlapDirectionByFid[static_cast<long long>(outFid)] =
                    featureDirCtx;
                const char* pathNames[] = { "topology", "vdp",
                    "best_hypothesis", "initial_ring", "strict_fallback" };
                std::cerr << "[MaskOnlyOverlapKey] output_fid=" << outFid
                          << " building_id=" << buildingId
                          << " priority=" << featurePriority
                          << " path=" << pathNames[static_cast<int>(featurePath)]
                          << std::endl;
            }
        } else {
            std::cerr << "[MaskOnly] CreateFeature failed fid=" << buildingId << std::endl;
        }
        OGRFeature::DestroyFeature(outFeature);
        OGRFeature::DestroyFeature(inFeature);
        if (total % 200 == 0) {
            std::cout << "  ...mask-only processing feature " << total << std::endl;
        }
    }
    auto loopEnd = std::chrono::steady_clock::now();

    // ---- 全局重叠解决: 所有单体规则化完成后统一处理建筑间相交 ----
    {
        MaskOnlyOverlapStats overlapStats;
        ResolveMaskOnlyOutputOverlaps(
            outLayer, originOffset,
            overlapPriorityByFid, overlapAlternateByFid,
            overlapDirectionByFid, overlapStats);
        std::cout << "[MaskOnlyOverlap] summary buffer_repaired="
                  << overlapStats.bufferRepaired
                  << " candidate_swaps=" << overlapStats.candidateSwaps
                  << " group_adjusted=" << overlapStats.groupAdjusted
                  << " translated=" << overlapStats.translatedCount
                  << " clipped=" << overlapStats.clipped << std::endl;
    }
    // Shapefile 的同一 GDALDataset 可能继续返回写入前的要素缓存。
    // 关闭并重新打开后再做一次最终修复，确保日志与磁盘文件一致。
    outDataset->FlushCache();
    GDALClose(outDataset);
    outDataset = nullptr;
    GDALDataset* reopenedOutput = static_cast<GDALDataset*>(GDALOpenEx(
        outputVectorIn.c_str(), GDAL_OF_VECTOR | GDAL_OF_UPDATE,
        nullptr, nullptr, nullptr));
    if (reopenedOutput) {
        OGRLayer* reopenedLayer = reopenedOutput->GetLayer(0);
        if (reopenedLayer) {
            MaskOnlyOverlapStats finalOverlapStats;
            ResolveMaskOnlyOutputOverlaps(
                reopenedLayer, originOffset,
                overlapPriorityByFid, overlapAlternateByFid,
                overlapDirectionByFid,
                finalOverlapStats);
            std::cout << "[MaskOnlyOverlap] reopened_final unresolved="
                      << finalOverlapStats.unresolved << std::endl;
            // 最终审计: 面积地板复查 + area 字段与最终几何一致
            {
                reopenedLayer->ResetReading();
                std::vector<GIntBig> toDelete;
                while (OGRFeature* f = reopenedLayer->GetNextFeature()) {
                    OGRGeometry* g = f->GetGeometryRef();
                    if (!g || g->IsEmpty()) {
                        toDelete.push_back(f->GetFID());
                        OGRFeature::DestroyFeature(f);
                        continue;
                    }
                    OGREnvelope env;
                    g->getEnvelope(&env);
                    const double a = GeometryArea(g);
                    const double bboxA =
                        (env.MaxX - env.MinX) * (env.MaxY - env.MinY);
                    if (a < 15.0 || bboxA < 20.0) {
                        std::cerr << "[MaskOnlySizeFloor] fid="
                                  << f->GetFID() << " area=" << a
                                  << " bbox=" << bboxA
                                  << " dropped=1 stage=final_audit"
                                  << std::endl;
                        toDelete.push_back(f->GetFID());
                    } else {
                        const int areaIdx = f->GetFieldIndex("area");
                        if (areaIdx >= 0) f->SetField(areaIdx, a);
                        reopenedLayer->SetFeature(f);
                    }
                    OGRFeature::DestroyFeature(f);
                }
                for (GIntBig fid : toDelete) {
                    reopenedLayer->DeleteFeature(fid);
                }
                if (!toDelete.empty()) {
                    std::cout << "[MaskOnlySizeFloor] final_audit dropped="
                              << toDelete.size() << std::endl;
                }
            }
        }
        reopenedOutput->FlushCache();
        GDALClose(reopenedOutput);
    } else {
        std::cerr << "[MaskOnlyOverlap] failed to reopen output for final audit"
                  << std::endl;
    }
    GDALClose(inDataset);

    // 路径统计汇总
    {
        const char* pathNames[5] = {
            "topology", "vdp", "best_hypothesis", "initial_ring",
            "strict_fallback"
        };
        for (int p = 0; p < 5; ++p) {
            const auto& st = pathStats[static_cast<std::size_t>(p)];
            if (st.count == 0) continue;
            std::cout << "[PathStats] path=" << pathNames[p]
                      << " count=" << st.count
                      << " avg_verts=" << st.sumVerts / st.count
                      << " avg_short_edges=" << static_cast<double>(st.sumShortEdges) / st.count
                      << " avg_area_ratio=" << st.sumAreaRatio / st.count
                      << std::endl;
        }
    }

    const std::filesystem::path debugBestPath = debugDir / "debug_best_hypothesis.shp";
    SaveDebugBestHypotheses(debugBestPath, debugCollector, srsClone, originOffset);
    std::cout << "[MaskOnly] best hypotheses saved: " << debugBestPath.string()
              << " (" << debugCollector.hypotheses.size() << " rings)" << std::endl;
    std::cout << "[MaskOnly] support LAS skipped: no OSGB source" << std::endl;

    // 双残差调试点输出(参与 Ceres 的 raw/smooth 残差采样)
    {
        const std::filesystem::path rawDbgPath =
            debugDir / "debug_mask_raw_residual_points.shp";
        if (outlineRegular::SaveRawResidualDebugDump(
                rawDbgPath.string(), originOffset, srsClone)) {
            std::cout << "[MaskOnly] raw residual debug points saved: "
                      << rawDbgPath.string() << std::endl;
        }
    }
    if (srsClone) srsClone->Release();

    const double loopSec =
        std::chrono::duration<double>(loopEnd - loopStart).count();
    std::cout << "[MaskOnly] summary initial_features=" << initialCount
              << " processed=" << total
              << " output_features=" << okCount
              << " empty_geom=" << emptyGeom
              << " small_area_skipped=" << smallSkipped
              << " holes_dropped=" << holesDropped
              << " self_intersections=" << selfIntersect
              << " fallback_hypothesis=" << fallbackHypothesis
              << " fallback_initial=" << fallbackInitial << std::endl;
    if (areaSamples > 0) {
        std::cout << "[MaskOnly] avg_vertices initial="
                  << (total > 0 ? sumInitVerts / total : 0.0)
                  << " best=" << (total > 0 ? sumBestVerts / total : 0.0)
                  << " final=" << (total > 0 ? sumFinalVerts / total : 0.0)
                  << " area_ratio_mean=" << sumAreaRatio / areaSamples << std::endl;
    }
    std::cout << "[MaskOnly] feature loop: " << loopSec << " s" << std::endl;
    std::cout << "[MaskOnly] Output: " << outputVector << std::endl;
    return okCount > 0 ? 0 : 2;
}


} // namespace

// ===== main =====
// 作用：程序入口。读 OSGB -> 从 TIF 提取初始轮廓 -> 按模型过滤 -> 规则化 -> 写出 SHP。
int main(int argc, char* argv[])
{
    SetConsoleOutputCP(65001u);
    ConsoleLogTee console_log(GetExecutableDirectory() / "console_output.txt");

    if (argc == 4 && std::string(argv[1]) == "--vectorize-mask") {
        GDALAllRegister();
        MaskVectorizationStats stats;
        if (!VectorizeBuildingMask(argv[2], argv[3], stats)) return 1;
        std::cout << "[Mask] raster=" << stats.width << "x" << stats.height
                  << ", source bands=" << stats.sourceBands
                  << ", used bands=" << stats.usedBands
                  << ", building pixels=" << stats.buildingPixels
                  << ", polygons=" << stats.polygonCount
                  << ", contained removed=" << stats.containedPolygonsRemoved << std::endl;
        std::cout << "[Mask separation] colors=" << stats.sourceColorCount
                  << (stats.colorLabelsPreserved ? " (preserved)" : " (binary fallback)")
                  << ", pixel size=" << stats.pixelSizeX << " x " << stats.pixelSizeY
                  << " m, seed erosion radius=" << stats.erosionRadiusPixels
                  << " px, seeds=" << stats.seedCount
                  << ", split components=" << stats.splitComponentCount
                  << ", narrow waist splits=" << stats.narrowWaistSplitCount << std::endl;
        std::cout << "[Mask extent] X=" << stats.minX << ".." << stats.maxX
                  << ", Y=" << stats.minY << ".." << stats.maxY << std::endl;
        std::cout << "[Mask] Output: " << argv[3] << std::endl;
        return 0;
    }

    if (argc == 4 && std::string(argv[1]) == "--mask-only") {
        return RunMaskOnlyMode(argv[2], argv[3]);
    }

    std::cout << "Select data mode: (1) OSGB + XML (default)  (2) Mask-only: ";
    std::string dataModeChoice;
    std::cin >> dataModeChoice;
    if (dataModeChoice == "2" || dataModeChoice == "m" || dataModeChoice == "M") {
        std::string maskOnlyTif;
        std::string maskOnlyOut;
        std::cout << "Select AI mask GeoTIFF..." << std::endl;
        if (!PickOpenTifFile(maskOnlyTif)) {
            std::cerr << "No mask GeoTIFF selected. Exiting." << std::endl;
            return 1;
        }
        std::cout << "Select output Shapefile (.shp) to save..." << std::endl;
        if (!PickSaveFile(maskOnlyOut)) {
            std::cerr << "No output Shapefile selected. Exiting." << std::endl;
            return 1;
        }
        return RunMaskOnlyMode(maskOnlyTif, maskOnlyOut);
    }

    // ---- 1) Select OSGB input ----
    std::string osgbDir;
    std::string metadataXml;
    std::string inputRaster;
    std::string inputVector;
    std::string outputVector;

    std::cout << "Select input OSGB folder..." << std::endl;
    if (!PickFolder(osgbDir)) {
        std::cerr << "No input OSGB folder selected. Exiting." << std::endl;
        return 1;
    }

    // ---- 2) 初始化 GDAL，确保输出目录存在，打印路径 ----
    GDALAllRegister();
    std::cout << "Input OSGB dir: " << osgbDir << std::endl;

    // ---- 2.5) 读 metadata XML，解析地理参考(SRSOrigin) ----
    //   OSGB 是相对坐标，点云也按相对坐标采样；这里先拿到 offset，供稍后保存 LAS 时
    //   叠加成世界坐标，以及后续把世界坐标的多边形换算回相对坐标。
std::cout << "Select metadata XML file..." << std::endl;
    Eigen::Vector3d metadataOffset = Eigen::Vector3d::Zero();
    if (PickOpenXmlFile(metadataXml)) {
        if (!ReadMetadataOffset(metadataXml, metadataOffset)) {
            return 1;
        }
    } else {
        std::cout << "No metadata XML selected. Using metadata offset 0,0,0." << std::endl;
    }

    // ---- 3) 读 OSGB，采样成点云(含法向量) ----
    auto s0 = std::chrono::steady_clock::now();
    MyCloudPtr sampled = std::make_shared<MyCloud>();
    if (!OSGMeshSampler::convertOSGBToMyCloud(osgbDir, sampled, kSampleDensity, kMaxPagedLODDepth) ||
        !sampled || !sampled->cloud || sampled->cloud->empty()) {
        std::cerr << "Failed to sample OSGB mesh." << std::endl;
        return 1;
    }
    auto s1 = std::chrono::steady_clock::now();
    std::cout << "[Timing] OSGB sampling: " << std::chrono::duration<double>(s1 - s0).count() << " s" << std::endl;

    if (!sampled->normal || sampled->normal->size() != sampled->cloud->size()) {
        std::cerr << "Sampled point and normal counts do not match." << std::endl;
        return 1;
    }

    // ---- 3.5) 建一次 2D KdTree(Z 置 0)，供后续每个多边形做邻域查询 ----
    auto kt0 = std::chrono::steady_clock::now();
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2d(new pcl::PointCloud<pcl::PointXYZ>);
    cloud2d->resize(sampled->cloud->size());
    for (std::size_t i = 0; i < sampled->cloud->size(); ++i) {
        const auto& p = sampled->cloud->points[i];
        cloud2d->points[i].x = p.x;
        cloud2d->points[i].y = p.y;
        cloud2d->points[i].z = 0.0f;
    }
    pcl::KdTreeFLANN<pcl::PointXYZ>::Ptr kdtree(new pcl::KdTreeFLANN<pcl::PointXYZ>);
    kdtree->setInputCloud(cloud2d);
    auto kt1 = std::chrono::steady_clock::now();
    std::cout << "[Timing] KdTree built in "
              << std::chrono::duration<double, std::milli>(kt1 - kt0).count()
              << " ms (" << cloud2d->size() << " points)" << std::endl;

    // ---- 3.6) 可选：保存采样点云为 LAS(点格式 3=GPS+RGB，二进制，含地理参考 offset) ----
    //          默认放在 build_deps_release\src\sampled_cloud.las；
    //          世界坐标 = 相对坐标 + SRSOrigin(offset 写入 LAS 文件头)。
    //          内存里的 sampled 仍保持相对坐标，后续处理不变。
std::cout << "閲囨牱鐐逛簯鍏?" << sampled->cloud->size()
              << " 个点。是否保存为 LAS? (y/n): ";
    char saveCloud = 'n';
    std::cin >> saveCloud;
    if (saveCloud == 'y' || saveCloud == 'Y') {
        char exePath[1024] = {0};
        GetModuleFileNameA(nullptr, exePath, 1024);
        // exe 在 build_deps_release\Release\，上一级再进 src = build_deps_release\src
        std::filesystem::path lasDir = std::filesystem::path(exePath).parent_path() / ".." / "src";
        std::filesystem::create_directories(lasDir);
        std::string lasPath = (lasDir / "sampled_cloud.las").string();

        // LAS 文件头 offset = 地理参考(SRSOrigin)；SaveSampledCloudAsLas 内部按
        // 世界坐标 = 相对坐标 + offset 写入，不改动 sampled->cloud 的相对坐标。
        sampled->offset = metadataOffset;
        sampled->hasoffset = true;

        std::cout << "姝ｅ湪淇濆瓨鐐逛簯鍒? " << lasPath << " ..." << std::endl;
        auto p0 = std::chrono::steady_clock::now();
        const bool ok = SaveSampledCloudAsLas(*sampled, lasPath);
        auto p1 = std::chrono::steady_clock::now();
        if (ok) {
            std::error_code ec;
            auto sz = std::filesystem::file_size(lasPath, ec);
            std::cout << "[LAS] 已保存 " << sampled->cloud->size() << " points, 耗时 "
                      << std::chrono::duration<double>(p1 - p0).count() << " s"
                      << (ec ? "" : (", 文件 " + std::to_string(sz) + " bytes")) << std::endl;
        } else {
            std::cerr << "[LAS] save failed." << std::endl;
        }
    }

    // ---- 4) Select footprint source: AI mask raster or existing vector ----
    std::cout << "Choose footprint source: image mask or existing vector? "
              << "(i=image, v=vector): ";
    char sourceMode = 'i';
    std::cin >> sourceMode;
    const bool useRasterMode = (sourceMode == 'i' || sourceMode == 'I');

    if (useRasterMode) {
    // ---- 4a) Extract initial outlines from the AI building mask ----
    std::cout << "Select AI building mask GeoTIFF..." << std::endl;
    if (!PickOpenTifFile(inputRaster)) {
        std::cerr << "No input GeoTIFF selected. Exiting." << std::endl;
        return 1;
    }

    char executablePath[MAX_PATH] = {0};
    GetModuleFileNameA(nullptr, executablePath, MAX_PATH);
    inputVector = PrepareWritableShapefilePath(
        std::filesystem::path(executablePath).parent_path() /
            "initial_building_outline.shp",
        "initial outline Shapefile").string();

    std::cout << "[Mask] Extracting all non-black building regions..." << std::endl;
    auto maskStart = std::chrono::steady_clock::now();
    MaskVectorizationStats maskStats;
    if (!VectorizeBuildingMask(inputRaster, inputVector, maskStats)) {
        std::cerr << "Failed to extract building outlines from GeoTIFF." << std::endl;
        return 1;
    }
    InitialOutlineMergeStats initialMergeStats;
    int mergePasses = 0;
    while (true) {
        InitialOutlineMergeStats passStats;
        const bool ok = MergeOversegmentedInitialOutlines(
            inputVector, sampled, kdtree, metadataOffset, passStats);
        if (!ok) {
            if (mergePasses == 0) {
                std::cerr << "Failed to evaluate adjacent initial outlines." << std::endl;
                return 1;
            }
            break;
        }
        ++mergePasses;
        initialMergeStats.adjacentCandidates += passStats.adjacentCandidates;
        initialMergeStats.keptByModelEvidence += passStats.keptByModelEvidence;
        initialMergeStats.keptNarrowSeam += passStats.keptNarrowSeam;
        initialMergeStats.mergedPairs += passStats.mergedPairs;
        initialMergeStats.mergedByBigPair += passStats.mergedByBigPair;
        initialMergeStats.mergedByShape += passStats.mergedByShape;
        initialMergeStats.mergedByParent += passStats.mergedByParent;
        initialMergeStats.removedFeatures += passStats.removedFeatures;
        std::cout << "[Mask merge pass " << mergePasses << "] accepted="
                  << passStats.mergedPairs << ", removed=" << passStats.removedFeatures
                  << ", narrow=" << passStats.keptNarrowSeam
                  << ", evidence=" << passStats.keptByModelEvidence
                  << ", big=" << passStats.mergedByBigPair
                  << ", shape=" << passStats.mergedByShape
                  << ", parent=" << passStats.mergedByParent << std::endl;
        if (passStats.mergedPairs == 0) break;
        if (mergePasses >= 20) {
            std::cerr << "[Mask merge] reached max 20 passes." << std::endl;
            break;
        }
    }

    NarrowNeckSplitStats neckSplitStats;
    if (!SplitInitialOutlinesAtNarrowNecks(inputVector, neckSplitStats)) {
        std::cerr << "[Mask neck split] failed to evaluate initial outlines." << std::endl;
    } else {
        std::cout << "[Mask neck split] inspected=" << neckSplitStats.inspectedFeatures
                  << ", split features=" << neckSplitStats.splitFeatures
                  << ", cuts=" << neckSplitStats.cuts
                  << ", created parts=" << neckSplitStats.createdParts
                  << ", candidate-only=" << neckSplitStats.candidateOnly
                  << ", reject-small-area=" << neckSplitStats.rejectedSmallPartArea
                  << ", reject-invalid-ring=" << neckSplitStats.rejectedInvalidRing
                  << ", reject-invalid-polygon=" << neckSplitStats.rejectedInvalidPolygon
                  << std::endl;
    }

    // ---- 合并收敛后：吸收被大轮廓完全包含的内部小轮廓(分色/栅格化残留) ----
    const long long containedFootprintsRemoved =
        RemoveContainedSmallFootprints(inputVector, kContainedFootprintMaxArea);
    if (containedFootprintsRemoved > 0) {
        std::cout << "[Mask] removed " << containedFootprintsRemoved
                  << " contained small footprints (<= " << kContainedFootprintMaxArea
                  << " m^2) inside larger outlines." << std::endl;
    }

    maskStats.polygonCount = std::max<long long>(
        0, maskStats.polygonCount - initialMergeStats.removedFeatures
               + (neckSplitStats.createdParts > 0
                    ? neckSplitStats.createdParts - neckSplitStats.splitFeatures
                    : 0)
               - containedFootprintsRemoved);
    auto maskEnd = std::chrono::steady_clock::now();
    std::cout << "[Mask] raster=" << maskStats.width << "x" << maskStats.height
              << ", source bands=" << maskStats.sourceBands
              << ", used bands=" << maskStats.usedBands
              << ", building pixels=" << maskStats.buildingPixels
              << ", polygons=" << maskStats.polygonCount
              << ", contained removed=" << maskStats.containedPolygonsRemoved << std::endl;
    std::cout << "[Mask separation] colors=" << maskStats.sourceColorCount
              << (maskStats.colorLabelsPreserved ? " (preserved)" : " (binary fallback)")
              << ", pixel size=" << maskStats.pixelSizeX << " x " << maskStats.pixelSizeY
              << " m, seed erosion radius=" << maskStats.erosionRadiusPixels
              << " px, seeds=" << maskStats.seedCount
              << ", split components=" << maskStats.splitComponentCount
              << ", narrow waist splits=" << maskStats.narrowWaistSplitCount << std::endl;
    std::cout << "[Mask merge] passes=" << mergePasses
              << ", total adjacent candidates=" << initialMergeStats.adjacentCandidates
              << ", narrow seams (not merged)=" << initialMergeStats.keptNarrowSeam
              << ", kept by evidence=" << initialMergeStats.keptByModelEvidence
              << ", accepted pairs=" << initialMergeStats.mergedPairs
              << ", shape merges=" << initialMergeStats.mergedByShape
              << ", parent merges=" << initialMergeStats.mergedByParent
              << ", removed features=" << initialMergeStats.removedFeatures << std::endl;
    std::cout << "[Mask neck split] split features=" << neckSplitStats.splitFeatures
              << ", cuts=" << neckSplitStats.cuts
              << ", created parts=" << neckSplitStats.createdParts
              << ", candidate-only=" << neckSplitStats.candidateOnly
              << ", reject-small-area=" << neckSplitStats.rejectedSmallPartArea
              << ", reject-invalid-ring=" << neckSplitStats.rejectedInvalidRing
              << ", reject-invalid-polygon=" << neckSplitStats.rejectedInvalidPolygon
              << std::endl;
    std::cout << "[Mask extent] X=" << maskStats.minX << ".." << maskStats.maxX
              << ", Y=" << maskStats.minY << ".." << maskStats.maxY << std::endl;
    // Stamp the area attribute after all mask stages (merge/split recreate
    // features, so any earlier area values would be stale).
    StampAreaField(inputVector);
    std::cout << "[Mask] Initial outlines saved: " << inputVector << std::endl;
    std::cout << "[Timing] mask vectorization: "
              << std::chrono::duration<double>(maskEnd - maskStart).count() << " s" << std::endl;
    } else {
        std::cout << "Select input footprint Shapefile..." << std::endl;
        if (!PickOpenFile(inputVector)) {
            std::cerr << "No input Shapefile selected. Exiting." << std::endl;
            return 1;
        }
    }

    std::cout << "Select output Shapefile (.shp) to save..." << std::endl;
    if (!PickSaveFile(outputVector)) {
        std::cerr << "No output Shapefile selected. Exiting." << std::endl;
        return 1;
    }

    if (std::filesystem::path(outputVector).lexically_normal() ==
        std::filesystem::path(inputVector).lexically_normal()) {
        std::cerr << "Output Shapefile must differ from the input footprint Shapefile." << std::endl;
        return 1;
    }

    std::filesystem::create_directories(std::filesystem::path(outputVector).parent_path());
    outputVector = PrepareWritableShapefilePath(
        std::filesystem::path(outputVector),
        "output Shapefile").string();
    std::cout << "Metadata XML:   " << (metadataXml.empty() ? "(none)" : metadataXml) << std::endl;
    std::cout << "SRSOrigin:      " << metadataOffset.x() << ", "
              << metadataOffset.y() << ", " << metadataOffset.z() << std::endl;
    std::cout << "Footprint mode: " << (useRasterMode ? "image mask" : "existing vector") << std::endl;
    if (useRasterMode) {
        std::cout << "Input raster:   " << inputRaster << std::endl;
        std::cout << "Initial vector: " << inputVector << std::endl;
    } else {
        std::cout << "Input vector:   " << inputVector << std::endl;
    }
    std::cout << "Output vector:  " << outputVector << std::endl;

    const Eigen::Vector3f cloudMin = sampled->boundingBox.min();
    const Eigen::Vector3f cloudMax = sampled->boundingBox.max();
    std::cout << "[Model extent local] X=" << cloudMin.x() << ".." << cloudMax.x()
              << ", Y=" << cloudMin.y() << ".." << cloudMax.y() << std::endl;
    std::cout << "[Model extent world] X=" << cloudMin.x() + metadataOffset.x()
              << ".." << cloudMax.x() + metadataOffset.x()
              << ", Y=" << cloudMin.y() + metadataOffset.y()
              << ".." << cloudMax.y() + metadataOffset.y() << std::endl;

    GDALDataset* inDataset = static_cast<GDALDataset*>(
        GDALOpenEx(inputVector.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!inDataset) {
        std::cerr << "Cannot open input vector: " << inputVector << std::endl;
        return 1;
    }

    OGRLayer* inLayer = inDataset->GetLayer(0);
    if (!inLayer) {
        std::cerr << "Input vector has no layer." << std::endl;
        GDALClose(inDataset);
        return 1;
    }

    // ---- 5) 创建输出 Shapefile(ESRI Shapefile 驱动，三维面 wkbPolygon25D) ----
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
    if (!driver) {
        std::cerr << "Cannot find ESRI Shapefile driver." << std::endl;
        GDALClose(inDataset);
        return 1;
    }

    GDALDataset* outDataset = driver->Create(outputVector.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (!outDataset) {
        std::cerr << "Cannot create output vector: " << outputVector << std::endl;
        GDALClose(inDataset);
        return 1;
    }

    OGRLayer* outLayer = outDataset->CreateLayer(
        "regularized_building", inLayer->GetSpatialRef(), wkbPolygon25D, nullptr);
    if (!outLayer) {
        std::cerr << "Cannot create output layer." << std::endl;
        GDALClose(outDataset);
        GDALClose(inDataset);
        return 1;
    }

    CopyFields(inLayer, outLayer);   // 复制字段结构
    // Ensure the output carries the stable building id (assigned once at mask
    // vectorization, propagated through merge/split stages); fall back to the
    // input FID when the field is absent (existing-vector input mode).
    const int inIdFieldIdx = inLayer->GetLayerDefn()->GetFieldIndex("id");
    int outIdFieldIdx = outLayer->GetLayerDefn()->GetFieldIndex("id");
    if (outIdFieldIdx < 0) {
        OGRFieldDefn idField("id", OFTInteger64);
        if (outLayer->CreateField(&idField) == OGRERR_NONE) {
            outIdFieldIdx = outLayer->GetLayerDefn()->GetFieldIndex("id");
        }
    }
    int outAreaFieldIdx = outLayer->GetLayerDefn()->GetFieldIndex("area");
    if (outAreaFieldIdx < 0) {
        OGRFieldDefn areaField("area", OFTReal);
        if (outLayer->CreateField(&areaField) == OGRERR_NONE) {
            outAreaFieldIdx = outLayer->GetLayerDefn()->GetFieldIndex("area");
        }
    }
    OGRFeatureDefn* outDefn = outLayer->GetLayerDefn();

    auto ownershipStart = std::chrono::steady_clock::now();
    SupportOwnershipContext supportOwnership =
        BuildSupportOwnershipContext(inLayer, metadataOffset, sampled, kdtree);
    auto ownershipEnd = std::chrono::steady_clock::now();
    std::cout << "[Support ownership] indexed rings=" << supportOwnership.records.size()
              << ", grid cells=" << supportOwnership.grid.size()
              << ", build time="
              << std::chrono::duration<double>(ownershipEnd - ownershipStart).count()
              << " s" << std::endl;

    // ---- 6) 遍历每个建筑物要素：规则化几何 + 复制属性 + 写出 ----
    int total = 0;
    int ok = 0;
    RegularizationDebugCollector debugCollector;
    auto loopStart = std::chrono::steady_clock::now();
    inLayer->ResetReading();
    while (OGRFeature* inFeature = inLayer->GetNextFeature()) {
        ++total;
        if (total % 50 == 0) {
            std::cout << "  ...processing feature " << total << std::endl;
        }
        GIntBig buildingId = inFeature->GetFID();
        if (inIdFieldIdx >= 0 && inFeature->IsFieldSetAndNotNull(inIdFieldIdx)) {
            buildingId = inFeature->GetFieldAsInteger64(inIdFieldIdx);
        }
        std::unique_ptr<OGRGeometry> outGeometry = RegularizeGeometry(
            inFeature->GetGeometryRef(), sampled, kdtree, metadataOffset,
            &supportOwnership,
            buildingId, &debugCollector);
        if (!outGeometry || outGeometry->IsEmpty()) {
            std::cerr << "[Output] skip empty geometry for fid=" << inFeature->GetFID() << std::endl;
            OGRFeature::DestroyFeature(inFeature);
            continue;
        }

        OGRFeature* outFeature = OGRFeature::CreateFeature(outDefn);
        if (!outFeature) {
            std::cerr << "[Output] failed to create feature for fid=" << inFeature->GetFID() << std::endl;
            OGRFeature::DestroyFeature(inFeature);
            continue;
        }

        CopyFieldValues(inFeature, outFeature);
        if (outIdFieldIdx >= 0) {
            outFeature->SetField(outIdFieldIdx, buildingId);
        }
        const OGRErr geomErr = outFeature->SetGeometry(outGeometry.get());
        if (outAreaFieldIdx >= 0) {
            outFeature->SetField(outAreaFieldIdx, GeometryArea(outGeometry.get()));
        }
        if (geomErr != OGRERR_NONE) {
            std::cerr << "[Output] SetGeometry failed for fid=" << inFeature->GetFID()
                      << ", geomType=" << outGeometry->getGeometryName()
                      << ", err=" << static_cast<int>(geomErr) << std::endl;
            OGRFeature::DestroyFeature(outFeature);
            OGRFeature::DestroyFeature(inFeature);
            continue;
        }

        const OGRErr createErr = outLayer->CreateFeature(outFeature);
        if (createErr != OGRERR_NONE) {
            std::cerr << "[Output] CreateFeature failed for fid=" << inFeature->GetFID()
                      << ", geomType=" << outGeometry->getGeometryName()
                      << ", err=" << static_cast<int>(createErr) << std::endl;
        } else {
            ++ok;
        }
        OGRFeature::DestroyFeature(outFeature);
        OGRFeature::DestroyFeature(inFeature);
    }
    auto loopEnd = std::chrono::steady_clock::now();
    outlineRegular::PrintHypothesisRepairSummary();

    const std::filesystem::path debugDir = std::filesystem::path(inputVector).parent_path();
    const std::filesystem::path debugBestPath = debugDir / "debug_best_hypothesis.shp";
    const std::filesystem::path debugSupportPath = debugDir / "debug_support_points.las";
    if (SaveDebugBestHypotheses(debugBestPath, debugCollector,
            inLayer->GetSpatialRef(), metadataOffset)) {
        std::cout << "[Debug] best hypotheses saved: " << debugBestPath.string()
                  << " (" << debugCollector.hypotheses.size() << " rings)" << std::endl;
    } else {
        std::cerr << "[Debug] failed to save best hypotheses: "
                  << debugBestPath.string() << std::endl;
    }
    if (SaveDebugSupportLas(debugSupportPath, debugCollector, metadataOffset)) {
        std::cout << "[Debug] support points saved: " << debugSupportPath.string()
                  << " (" << (debugCollector.support ? debugCollector.support->size() : 0)
                  << " points)" << std::endl;
    } else {
        std::cerr << "[Debug] failed to save support points: "
                  << debugSupportPath.string() << std::endl;
    }

    // Temporarily skip output overlap repair until write-out stability is verified.
    // OutputOverlapRepairStats overlapStats;
    // if (!ResolveOutputOverlaps(outLayer, overlapStats)) {
    //     std::cerr << "[Output overlap] failed to repair overlapping polygons." << std::endl;
    // } else if (overlapStats.resolvedPairs > 0 || overlapStats.overlapPairs > 0) {
    //     std::cout << "[Output overlap] candidate pairs=" << overlapStats.candidatePairs
    //               << ", remaining overlap pairs=" << overlapStats.overlapPairs
    //               << ", resolved pairs=" << overlapStats.resolvedPairs
    //               << ", shifted features=" << overlapStats.shiftedFeatures
    //               << ", optimized groups=" << overlapStats.optimizedGroups
    //               << ", unresolved pairs=" << overlapStats.unresolvedPairs
    //               << ", max remaining overlap area=" << overlapStats.maxOverlapArea
    //               << ", max shift=" << overlapStats.maxShiftDistance << std::endl;
    // }

    // ---- 计时汇总：定位"支撑点提取"还是"规则化优化"是瓶颈 ----
    double loopSec = std::chrono::duration<double>(loopEnd - loopStart).count();
    std::cout << "[Timing] feature loop: " << loopSec << " s (" << total << " features)" << std::endl;
    std::cout << "[Timing]   support extraction: " << g_supportTime << " s" << std::endl;
    std::cout << "[Timing]   regularization opt: " << g_optimizeTime << " s" << std::endl;
    std::cout << "[Filter] removed polygons with bbox area < " << kMinPolygonBBoxArea
              << " m^2: " << g_removedSmallPolygons << std::endl;
    std::cout << "[Filter] removed polygons outside sampled model: "
              << g_removedOutsideModel << std::endl;
    std::cout << "[Support ownership] rejected_other=" << g_supportOwnershipRejected
              << ", ambiguous_rejected=" << g_supportOwnershipAmbiguousRejected << std::endl;
    std::cout << "[Support filter] coarse_wall=" << g_supportFilterCoarseWallTotal
              << ", effective_wall=" << g_supportFilterEffectiveWallTotal
              << ", wall_wide=" << g_supportFilterWideWallTotal
              << ", fallback_used=" << g_supportFilterFallbackUsedTotal << std::endl;
    if (total > 0) {
        std::cout << "[Timing]   avg per feature: " << (loopSec / total) << " s" << std::endl;
    }

    GDALClose(outDataset);
    GDALClose(inDataset);

    // ---- 7) 打印结果统计 ----
    std::cout << "Regularization finished. success=" << ok << " / total=" << total << std::endl;
    std::cout << "Output: " << outputVector << std::endl;
    return ok > 0 ? 0 : 2;
}

// =============================================================================
// OSGMeshSampler.cpp
// 作用：实现 OSGB -> 点云 的转换。
// 整体流程(convertOSGBToMyCloud)：
//   1) 检查 OSG 是否有读 .osgb 的插件；
//   2) 递归枚举目录下所有 .osgb 文件；
//   3) 用 osgDB::readNodeFile 逐个读入，并用 TriangleExtractor(一个 NodeVisitor)
//      遍历场景图，抽出所有三角面(同时递归加载 PagedLOD 的外部子瓦片)；
//   4) 对收集到的三角面做"面积加权随机采样"，得到带颜色和法向量的点云。
// =============================================================================

#include "OSGMeshSampler.h"

#include <osgDB/ReadFile>
#include <osgDB/Registry>
#include <osg/Node>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/PagedLOD>
#include <osg/NodeVisitor>
#include <osg/Vec3>
#include <osg/Vec4>
#include <osg/Vec2>
#include <osg/Texture2D>
#include <osg/Image>
#include <osg/StateSet>
#include <osg/Matrix>

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <filesystem>
#include <unordered_set>
#include <limits>
#include <random>

namespace OSGMeshSampler {

// ---------------------------------------------------------------------------
// Triangle：单个三角面的数据(三个顶点、三个颜色、面积、编号)
// ---------------------------------------------------------------------------
struct Triangle {
    osg::Vec3d v0, v1, v2;
    osg::Vec4 c0, c1, c2;
    float    area;
    uint32_t triId;
};

// ---------------------------------------------------------------------------
// TriangleExtractor：一个 osg::NodeVisitor，遍历场景图抽取所有三角面。
// 关键能力：遇到 PagedLOD 时，会把其引用、尚未加载的外部子瓦片(.osgb)也读入并继续遍历，
//           从而把分页 LOD 的多级瓦片"摊平"成全部三角面。
// ---------------------------------------------------------------------------
class TriangleExtractor : public osg::NodeVisitor
{
public:
    std::vector<Triangle> tris;                 // 收集到的所有三角面
    std::string baseDir;                        // 数据根目录(用于解析 PagedLOD 的相对文件名)
    int maxDepth;                               // PagedLOD 最大递归深度
    int currentDepth;                           // 当前递归深度
    uint32_t triCounter;                        // 三角面编号计数器
    std::unordered_set<std::string> visitedFiles; // 已访问过的瓦片文件(避免重复加载)

    TriangleExtractor(const std::string& dir, int depth)
        : osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN)
        , baseDir(dir)
        , maxDepth(depth)
        , currentDepth(0)
        , triCounter(0)
    {}

    // ===== apply(Geode) =====
    // 作用：访问 Geode(叶子几何节点)时，取出其中的 Geometry，抽取三角面。
    //       顶点先用 localToWorld 矩阵变换到世界坐标。
    void apply(osg::Geode& geode) override
    {
		const osg::Matrixd localToWorld = osg::computeLocalToWorld(getNodePath());
        for (unsigned int i = 0; i < geode.getNumDrawables(); ++i) {
            osg::Geometry* geom = geode.getDrawable(i)->asGeometry();
            if (!geom) continue;
            processGeometry(geom, localToWorld);
        }
        traverse(geode);
    }

    // ===== apply(PagedLOD) =====
    // 作用：访问 PagedLOD 时，除了遍历已加载的子节点，还要把它"按需引用"的
    //       外部子瓦片文件读进来并继续遍历，实现分页 LOD 的递归展开。
    void apply(osg::PagedLOD& plod) override
    {
        if (currentDepth >= maxDepth) { traverse(plod); return; }

        // 先遍历已经加载进来的子节点
        for (unsigned int i = 0; i < plod.getNumChildren(); ++i) {
            osg::Node* child = plod.getChild(i);
            if (child) child->accept(*this);
        }

        // 再加载 PagedLOD 引用、但尚未加载的外部文件
        std::string dbPath = plod.getDatabasePath();
        for (unsigned int i = 0; i < plod.getNumFileNames(); ++i) {
            std::string fname = plod.getFileName(i);
            if (fname.empty()) continue;

            // 解析子瓦片的完整路径：优先用 PagedLOD 自带的 databasePath，否则用 baseDir
            std::string fullPath;
            if (!dbPath.empty()) {
                fullPath = (std::filesystem::path(dbPath) / fname).string();
            } else {
                fullPath = (std::filesystem::path(baseDir) / fname).string();
            }

            // 已访问过的瓦片跳过，避免重复加载/死循环
            std::string canonical = std::filesystem::weakly_canonical(fullPath).string();
            if (visitedFiles.count(canonical)) continue;
            visitedFiles.insert(canonical);

            osg::ref_ptr<osg::Node> childNode = osgDB::readNodeFile(fullPath);
            if (childNode) {
                std::cout << "  [PagedLOD] Loading: " << fullPath << std::endl;
                currentDepth++;
                childNode->accept(*this);   // 递归展开子瓦片
                currentDepth--;
            }
        }
    }

private:
    // ===== processGeometry =====
    // 作用：从一个 osg::Geometry 中抽取所有三角面。
    //       颜色来源优先级：逐顶点颜色 > 整体颜色 > 纹理采样 > 默认白色。
    //       顶点用 localToWorld 变换到世界坐标，并计算每个三角面的面积。
    void processGeometry(osg::Geometry* geom, const osg::Matrixd& localToWorld)
    {
        osg::Vec3Array* verts = dynamic_cast<osg::Vec3Array*>(geom->getVertexArray());
        if (!verts || verts->empty()) return;

        // 颜色：先判断是逐顶点色、整体色，还是没有颜色(需从纹理采样)
        osg::Vec4Array* colors = dynamic_cast<osg::Vec4Array*>(geom->getColorArray());
        bool hasPerVertexColor = (colors && colors->size() >= verts->size());
        bool hasOverallColor = false;
        osg::Vec4 overallColor(1,1,1,1);
        if (!hasPerVertexColor && colors && !colors->empty()) {
            hasOverallColor = true;
            overallColor = (*colors)[0];
        }

        // 没有顶点/整体颜色时，尝试从纹理 + UV 采样得到颜色
        bool sampleTexColor = (!hasPerVertexColor && !hasOverallColor);
        osg::Image* texImage = nullptr;
        osg::Vec2Array* texCoords = nullptr;
        if (sampleTexColor) {
            texCoords = dynamic_cast<osg::Vec2Array*>(geom->getTexCoordArray(0));
            if (texCoords) {
                osg::StateSet* ss = geom->getStateSet();
                if (ss) {
                    osg::Texture2D* tex = dynamic_cast<osg::Texture2D*>(
                        ss->getTextureAttribute(0, osg::StateAttribute::TEXTURE));
                    if (tex) texImage = tex->getImage();
                }
            }
            if (!texImage) sampleTexColor = false;
        }

        // 遍历所有图元集，只处理 GL_TRIANGLES，每 3 个索引一个三角面
        for (unsigned int ps = 0; ps < geom->getNumPrimitiveSets(); ++ps) {
            osg::PrimitiveSet* prim = geom->getPrimitiveSet(ps);
            if (!prim) continue;
            GLenum mode = prim->getMode();
            if (mode != GL_TRIANGLES) continue;

            unsigned int cnt = prim->getNumIndices();
            for (unsigned int j = 0; j + 2 < cnt; j += 3) {
                unsigned int i0 = prim->index(j);
                unsigned int i1 = prim->index(j+1);
                unsigned int i2 = prim->index(j+2);
                if (i0 >= verts->size() || i1 >= verts->size() || i2 >= verts->size())
                    continue;

                Triangle tri;
				// 顶点变换到世界坐标
				const osg::Vec3d worldV0 = osg::Vec3d((*verts)[i0]) * localToWorld;
				const osg::Vec3d worldV1 = osg::Vec3d((*verts)[i1]) * localToWorld;
				const osg::Vec3d worldV2 = osg::Vec3d((*verts)[i2]) * localToWorld;
				tri.v0 = worldV0;
				tri.v1 = worldV1;
				tri.v2 = worldV2;
                tri.triId = triCounter++;

                // 按颜色来源填充三个顶点的颜色
                if (hasPerVertexColor) {
                    tri.c0 = (*colors)[i0]; tri.c1 = (*colors)[i1]; tri.c2 = (*colors)[i2];
                } else if (hasOverallColor) {
                    tri.c0 = tri.c1 = tri.c2 = overallColor;
                } else if (sampleTexColor && texCoords) {
                    tri.c0 = sampleTex(texImage, (*texCoords)[i0]);
                    tri.c1 = sampleTex(texImage, (*texCoords)[i1]);
                    tri.c2 = sampleTex(texImage, (*texCoords)[i2]);
                } else {
                    tri.c0 = tri.c1 = tri.c2 = osg::Vec4(1,1,1,1);
                }

                // 三角面面积 = 0.5 * |叉积|；面积过小(退化三角形)则丢弃
                osg::Vec3d e1 = tri.v1 - tri.v0;
                osg::Vec3d e2 = tri.v2 - tri.v0;
                tri.area = 0.5f * (e1 ^ e2).length();
                if (tri.area > 1e-10f) tris.push_back(tri);
            }
        }
    }

    // ===== sampleTex =====
    // 作用：从纹理图像 img 中，按 uv 坐标采样一个像素颜色，返回归一化的 RGBA。
    static osg::Vec4 sampleTex(osg::Image* img, const osg::Vec2& uv)
    {
        if (!img || !img->data()) return osg::Vec4(1,1,1,1);
        int s = img->s(), t = img->t();
        if (s <= 0 || t <= 0) return osg::Vec4(1,1,1,1);
        float u = std::max(0.0f, std::min(0.9999f, uv.x()));
        float v = std::max(0.0f, std::min(0.9999f, uv.y()));
        int x = static_cast<int>(u * (s-1));
        int y = static_cast<int>(v * (t-1));
        unsigned char* px = (unsigned char*)img->data(x, y);
        float a = (img->getPixelSizeInBits() >= 32) ? px[3]/255.0f : 1.0f;
        return osg::Vec4(px[0]/255.0f, px[1]/255.0f, px[2]/255.0f, a);
    }
};

// ---------------------------------------------------------------------------
// sampleTriangles：对收集到的三角面做"面积加权随机采样"，生成点云。
// 每个三角面采样点数 ≈ 面积 * density；点在三角形内按重心坐标均匀分布，
// 颜色按重心坐标插值，法向量取三角面法向。
// ---------------------------------------------------------------------------
static void sampleTriangles(
    const std::vector<Triangle>& triangles,
    float density,
    PointCloudT::Ptr& outCloud,
    pcl::PointCloud<pcl::Normal>::Ptr& outNormal)
{
    std::mt19937 rng(std::random_device{}());          // 随机数引擎
    std::uniform_real_distribution<float> dist01(0.0f, 1.0f);

    // 预估总点数并预分配
    size_t est = 0;
    for (const auto& tri : triangles)
        est += static_cast<size_t>(std::max(1, static_cast<int>(tri.area * density)));
    outCloud->reserve(est);
    outNormal->reserve(est);

    for (const auto& tri : triangles) {
        int n = std::max(1, static_cast<int>(tri.area * density)); // 该三角面采样点数

        // 计算三角面法向量(叉积归一化)
        osg::Vec3d e1 = tri.v1 - tri.v0;
        osg::Vec3d e2 = tri.v2 - tri.v0;
        osg::Vec3d cross = e1 ^ e2;
        double len = cross.length();
        if (len < 1e-10f) continue;
        osg::Vec3d nOSG = cross / len;

        pcl::Normal pn;
        pn.normal_x = nOSG.x(); pn.normal_y = nOSG.y(); pn.normal_z = nOSG.z();
        pn.curvature = 0.0f;

        // 在三角形内随机撒 n 个点(用重心坐标保证均匀分布)
        for (int s = 0; s < n; ++s) {
            float r1 = dist01(rng), r2 = dist01(rng);
            float sq = std::sqrt(r1);
            float b0 = 1.0f - sq;
            float b1 = sq * (1.0f - r2);
            float b2 = sq * r2;

            osg::Vec3d p = tri.v0*b0 + tri.v1*b1 + tri.v2*b2 ;
            osg::Vec4 c = tri.c0*b0 + tri.c1*b1 + tri.c2*b2;   // 颜色(重心插值)

            PointT pt;
            pt.x = p.x(); pt.y = p.y(); pt.z = p.z();
            pt.r = (uint8_t)std::min(255.0f, std::max(0.0f, c.r()*255.0f));
            pt.g = (uint8_t)std::min(255.0f, std::max(0.0f, c.g()*255.0f));
            pt.b = (uint8_t)std::min(255.0f, std::max(0.0f, c.b()*255.0f));
            pt.a = (uint8_t)std::min(255.0f, std::max(0.0f, c.a()*255.0f));

            outCloud->push_back(pt);
            outNormal->push_back(pn);
        }
    }
    outCloud->points.resize(outCloud->size());
    outNormal->points.resize(outNormal->size());
}

// ---------------------------------------------------------------------------
// 对外接口：convertOSGBToMyCloud(见头文件注释)
// ---------------------------------------------------------------------------
bool convertOSGBToMyCloud(
    const std::string& osgbDir,
    MyCloudPtr& outputCloud,
    float sampleDensity,
    int maxDepth)
{
    if (!outputCloud) outputCloud = std::make_shared<MyCloud>();

    // 1) 确认 OSG 能读 .osgb(需要 osgdb_osg 插件已加载)
    if (!osgDB::Registry::instance()->getReaderWriterForExtension("osgb")) {
        std::cerr << "[OSGMeshSampler] No OSGB plugin!" << std::endl;
        return false;
    }

    // 2) 检查目录存在
    std::filesystem::path dir(osgbDir);
    if (!std::filesystem::exists(dir) || !std::filesystem::is_directory(dir)) {
        std::cerr << "[OSGMeshSampler] Dir not found: " << osgbDir << std::endl;
        return false;
    }

    // 3) 递归收集目录下所有 .osgb 文件(兼容 Tile_xxx/ 子目录结构)
    std::vector<std::filesystem::path> osgbFiles;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() == ".osgb") {
            osgbFiles.push_back(entry.path());
        }
    }
    if (osgbFiles.empty()) {
        std::cerr << "[OSGMeshSampler] No .osgb in " << osgbDir << std::endl;
        return false;
    }

    std::cout << "[OSGMeshSampler] " << osgbFiles.size() << " .osgb files, density="
              << sampleDensity << std::endl;

    // 4) 逐个读入 .osgb，用 TriangleExtractor 抽取三角面
    std::vector<Triangle> allTris;
    std::unordered_set<std::string> globalVisited;

    for (const auto& fp : osgbFiles) {
        std::string canon = std::filesystem::weakly_canonical(fp).string();
        if (globalVisited.count(canon)) continue;   // 跨文件去重
        globalVisited.insert(canon);

        std::cout << "  Loading: " << fp.string() << std::endl;
        osg::ref_ptr<osg::Node> root = osgDB::readNodeFile(fp.string());
        if (!root) { std::cerr << "  Failed: " << fp.string() << std::endl; continue; }

        TriangleExtractor ext(osgbDir, maxDepth);
        ext.visitedFiles = globalVisited;           // 共享已访问集合
        root->accept(ext);
        for (const auto& v : ext.visitedFiles) globalVisited.insert(v);

        std::cout << "    " << ext.tris.size() << " triangles" << std::endl;
        allTris.insert(allTris.end(), ext.tris.begin(), ext.tris.end());
    }

    std::cout << "[OSGMeshSampler] Total triangles: " << allTris.size() << std::endl;
    if (allTris.empty()) { std::cerr << "[OSGMeshSampler] No triangles!" << std::endl; return false; }


    // 5) 面积加权采样，生成点云 + 法向量
    std::cout << "[OSGMeshSampler] Sampling..." << std::endl;
    sampleTriangles(allTris, sampleDensity, outputCloud->cloud, outputCloud->normal);

    // 6) 填充 MyCloud 的状态字段并计算包围盒
    outputCloud->points_num = static_cast<int>(outputCloud->cloud->size());
    outputCloud->hasCloud = true;
    outputCloud->hasMesh = false;
    outputCloud->isValid = true;

    if (!outputCloud->cloud->empty()) {
        for (const auto& pt : outputCloud->cloud->points)
            outputCloud->boundingBox.extend(Eigen::Vector3f(pt.x, pt.y, pt.z));
    }

    std::cout << "[OSGMeshSampler] Done. " << outputCloud->cloud->size() << " points" << std::endl;
    return true;
}

} // namespace OSGMeshSampler

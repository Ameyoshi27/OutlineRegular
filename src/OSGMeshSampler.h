// =============================================================================
// OSGMeshSampler.h
// 作用：声明 OSGMeshSampler 命名空间下的核心接口 convertOSGBToMyCloud。
//       功能是把一个目录下的倾斜摄影 OSGB 瓦片读入 OpenSceneGraph，
//       把网格三角面"面积加权随机采样"成带颜色和法向量的点云，存入 MyCloud。
// =============================================================================

#ifndef OSG_MESH_SAMPLER_H
#define OSG_MESH_SAMPLER_H

#include "MyCloud.h"
#include <string>
#include <vector>
#include <random>
#include <iostream>

namespace osg { class Node; }

namespace OSGMeshSampler {

// convertOSGBToMyCloud：读取 osgbDir 下所有 .osgb(含 PagedLOD 子瓦片)，
// 把网格采样为点云写入 outputCloud。
// 参数：osgbDir      - OSGB 数据所在目录(会递归查找 .osgb)；
//       outputCloud  - 输出的点云容器(MyCloud)；
//       sampleDensity- 采样密度(每平方米采样点数)，越大点越密；
//       maxDepth     - PagedLOD 递归加载的最大深度，防止无限下钻。
// 返回：成功(采样到点)返回 true，否则 false。
bool convertOSGBToMyCloud(
    const std::string& osgbDir,
    MyCloudPtr& outputCloud,
    float sampleDensity = 50.0f,
    int maxDepth = 10);

} // namespace OSGMeshSampler

#endif // OSG_MESH_SAMPLER_H

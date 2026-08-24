// =============================================================================
// MyCloud.h
// 作用：定义本项目通用的"点云容器"数据结构 MyCloud。
//       把一次采样的结果打包在一起：带颜色的点云、法向量、网格(可选)、
//       轴对齐包围盒以及若干状态标志，方便在 OSGMeshSampler / outlineRegular /
//       main 等模块之间统一传递。
// =============================================================================

#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/PolygonMesh.h>

#include <memory>
#include <vector>

using PointT = pcl::PointXYZRGBA;              // 带颜色(RGBA)的三维点类型
using PointCloudT = pcl::PointCloud<PointT>;   // 对应的点云类型

// MyCloud：项目内部统一的点云数据封装
class MyCloud {
public:
    // 构造时为点云 / 法向量 / 网格分配空对象，并将包围盒置空
    MyCloud()
        : cloud(new PointCloudT),
          normal(new pcl::PointCloud<pcl::Normal>),
          mesh(new pcl::PolygonMesh)
    {
        boundingBox.setEmpty();
    }

    PointCloudT::Ptr cloud;                     // 采样得到的带颜色点云
    pcl::PointCloud<pcl::Normal>::Ptr normal;   // 与 cloud 一一对应的法向量
    pcl::PolygonMesh::Ptr mesh;                 // 网格(本项目主要用点云，网格保留备用)

    int points_num = 0;       // 点数
    bool hasCloud = false;    // 是否已有点云数据
    bool hasMesh = false;     // 是否已有网格数据
    bool isValid = false;     // 数据是否有效(采样成功)
    bool haslabel = false;    // 是否带标签(预留)

    Eigen::AlignedBox3f boundingBox; // 轴对齐包围盒(AABB)，用于空间范围判断

    // 点云的地理参考偏移(LAS header offset)：世界坐标 = 相对坐标 + offset。
    // OSGB 按相对坐标采样，保存 LAS 时把 offset 写入文件头并叠加到坐标上；
    // 内存里的点云始终保持相对坐标，后续处理不变。移植自 E:\jt 的 MyCloud。
    Eigen::Vector3d offset{Eigen::Vector3d::Zero()};
    bool hasoffset = false;
};

using MyCloudPtr = std::shared_ptr<MyCloud>; // MyCloud 的共享指针，模块间传递用

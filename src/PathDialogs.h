// =============================================================================
// PathDialogs.h
// 作用：声明三个原生 Windows 文件/文件夹选择对话框函数。
//       程序运行时弹窗让用户选择：输入 OSGB 文件夹、输入 Shapefile、输出 Shapefile 路径。
//       仅使用 Windows SDK(IFileDialog)，不引入任何第三方库。
// =============================================================================

#pragma once

#include <string>

bool PickOpenXmlFile(std::string& out);
bool PickOpenTifFile(std::string& out);

// 弹出文件夹选择对话框，选中后把路径(UTF-8)写入 out。
// 返回 true=成功，false=用户取消。
bool PickFolder(std::string& out);    // 选择文件夹(输入 OSGB 目录)

// 弹出文件打开对话框(过滤 *.shp)，选中后把路径(UTF-8)写入 out。
// 返回 true=成功，false=用户取消。
bool PickOpenFile(std::string& out);  // 选择已存在的 .shp(输入)

// 弹出文件保存对话框(过滤 *.shp，默认文件名 regularized_building.shp)，选中后把路径(UTF-8)写入 out。
// 返回 true=成功，false=用户取消。
bool PickSaveFile(std::string& out);  // 选择要保存的 .shp 路径(输出)

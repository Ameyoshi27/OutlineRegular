# outlineRegularTool

Windows C++ 控制台工具，当前处理链为：

```text
倾斜 OSGB 模型 + AI 彩色建筑掩膜 GeoTIFF
  -> 初始建筑轮廓 Shapefile
  -> 按 OSGB 采样点云过滤模型范围外轮廓
  -> 规则化建筑轮廓 Shapefile
```

第三方库已打包在 `deps/x64-windows`，不需要另行安装 vcpkg。

## 运行流程

1. 选择输入 OSGB 文件夹，程序采样网格并建立二维 KdTree。
2. 可选保存采样点云 PLY。
3. 可选选择 OSGB 的 `metadata.xml`；取消时使用零 offset。
4. 选择多波段建筑掩膜 GeoTIFF。任一有效颜色波段非零的像素视为建筑，黑色视为背景。
5. 程序使用 GDAL 将二值掩膜矢量化，并在 exe 目录保存完整的
   `initial_building_outline.shp`。该文件保留 GeoTIFF 的投影和绝对坐标。
6. 选择规则化结果的输出 Shapefile 路径。
7. 对每个初始轮廓先减去 metadata offset，再检查二维包围盒面积和采样点云证据。
   小于 20 平方米或不在当前 OSGB 模型覆盖范围内的轮廓不参与规则化。
8. 对保留轮廓提取墙面/边界支撑点并执行规则化，输出时加回 metadata offset。

彩色掩膜会先合并为“建筑/非建筑”二值结果，因此相接且中间没有黑色间隔的建筑会形成一个连通面。

## 构建

双击 `build_vs2019.bat`，或执行：

```bat
cmake -S . -B build_deps_release -G "Visual Studio 16 2019" -A x64 -T host=x64
cmake --build build_deps_release --config Release --target outlineRegularTool
```

产物为 `build_deps_release/Release/outlineRegularTool.exe`。

## 独立验证掩膜

不运行 OSGB 采样时，可以直接验证 TIF 矢量化：

```bat
outlineRegularTool.exe --vectorize-mask input.tif output.shp
```

控制台会打印影像尺寸、使用波段数、建筑像素数、面数量和地理范围。

## 运行依赖

`run_release.bat` 会把 `deps/x64-windows/bin` 和 OSG 插件目录加入运行环境。
如果提示 `[OSGMeshSampler] No OSGB plugin!`，请检查
`deps/x64-windows/plugins/osgPlugins-3.6.5/osgdb_osgb.dll`。

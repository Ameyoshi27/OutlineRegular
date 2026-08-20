cmake_minimum_required(VERSION 3.16)

project(outlineRegularTool LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "" FORCE)

set(LOCAL_DEPS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/deps/x64-windows")
if(EXISTS "${LOCAL_DEPS_ROOT}")
    list(APPEND CMAKE_PREFIX_PATH "${LOCAL_DEPS_ROOT}")
    list(APPEND CMAKE_PREFIX_PATH "${LOCAL_DEPS_ROOT}/share")
    list(APPEND CMAKE_PREFIX_PATH "${LOCAL_DEPS_ROOT}/lib/cmake")
    list(APPEND CMAKE_INCLUDE_PATH "${LOCAL_DEPS_ROOT}/include")
    list(APPEND CMAKE_LIBRARY_PATH "${LOCAL_DEPS_ROOT}/lib")
elseif(EXISTS "F:/vcpkg/vcpkg/installed/x64-windows")
    list(APPEND CMAKE_PREFIX_PATH "F:/vcpkg/vcpkg/installed/x64-windows")
    list(APPEND CMAKE_PREFIX_PATH "F:/vcpkg/vcpkg/installed/x64-windows/share")
elseif(EXISTS "C:/vcpkg-master/vcpkg/installed/x64-windows")
    list(APPEND CMAKE_PREFIX_PATH "C:/vcpkg-master/vcpkg/installed/x64-windows")
    list(APPEND CMAKE_PREFIX_PATH "C:/vcpkg-master/vcpkg/installed/x64-windows/share")
endif()

if(MSVC)
    add_compile_options(/utf-8 /bigobj /MP /wd4996 /wd4819 /wd4251 /wd4267 /wd4244 /wd4275 /wd4005 /wd5033 /wd4018)
    add_compile_options("$<$<CONFIG:Release>:/Zi>")
    add_link_options("$<$<CONFIG:Release>:/DEBUG>")
endif()

find_package(PCL REQUIRED COMPONENTS common io kdtree surface search filters)
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
set(GFLAGS_USE_TARGET_NAMESPACE TRUE)
find_package(gflags REQUIRED)
find_package(glog REQUIRED)
find_package(Ceres REQUIRED)
find_package(GDAL REQUIRED)
if(EXISTS "${LOCAL_DEPS_ROOT}/lib/opencv_core4.lib" AND
   EXISTS "${LOCAL_DEPS_ROOT}/lib/opencv_imgproc4.lib")
    # The bundled OpenCV config imports every module and therefore requires a
    # missing protoc executable. Morphology needs only these two existing libs.
    set(OPENCV_MORPH_INCLUDE_DIRS "${LOCAL_DEPS_ROOT}/include")
    set(OPENCV_MORPH_LIBS
        "${LOCAL_DEPS_ROOT}/lib/opencv_core4.lib"
        "${LOCAL_DEPS_ROOT}/lib/opencv_imgproc4.lib")
else()
    find_package(OpenCV REQUIRED COMPONENTS core imgproc)
    set(OPENCV_MORPH_INCLUDE_DIRS ${OpenCV_INCLUDE_DIRS})
    set(OPENCV_MORPH_LIBS ${OpenCV_LIBS})
endif()
find_package(OpenSceneGraph REQUIRED COMPONENTS osgDB osgUtil osg)

add_executable(outlineRegularTool
    src/main.cpp
    src/PathDialogs.cpp
    src/PathDialogs.h
    src/MaskVectorizer.cpp
    src/MaskVectorizer.h
    src/outlineRegular.cpp
    src/outlineRegular.h
    src/modularFunction.cpp
    src/modularFunction.h
    src/OSGMeshSampler.cpp
    src/OSGMeshSampler.h
    src/MyCloud.h
)

target_include_directories(outlineRegularTool PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${PCL_INCLUDE_DIRS}
    ${OPENCV_MORPH_INCLUDE_DIRS}
    ${LOCAL_DEPS_ROOT}/include      # laszip/laszip_api.h
)

target_compile_definitions(outlineRegularTool PRIVATE
    _USE_MATH_DEFINES
    NOMINMAX
    GLOG_NO_ABBREVIATED_SEVERITIES
    LASZIP_DYN_LINK                 # laszip 动态链接：使 laszip_api.h 走 __declspec(dllimport)
)

target_link_directories(outlineRegularTool PRIVATE ${PCL_LIBRARY_DIRS})
target_compile_definitions(outlineRegularTool PRIVATE ${PCL_DEFINITIONS})

target_link_libraries(outlineRegularTool PRIVATE
    ${PCL_LIBRARIES}
    Eigen3::Eigen
    CGAL::CGAL
    Ceres::ceres
    ${OPENCV_MORPH_LIBS}
    ${OPENSCENEGRAPH_LIBRARIES}
    ${LOCAL_DEPS_ROOT}/lib/laszip3.lib   # LAS 写出(laszip 动态库的导入库)
)

target_include_directories(outlineRegularTool PRIVATE ${OPENSCENEGRAPH_INCLUDE_DIRS})

if(TARGET GDAL::GDAL)
    target_link_libraries(outlineRegularTool PRIVATE GDAL::GDAL)
else()
    target_include_directories(outlineRegularTool PRIVATE ${GDAL_INCLUDE_DIR})
    target_link_libraries(outlineRegularTool PRIVATE ${GDAL_LIBRARY})
endif()

# 运行时依赖：把 laszip3.dll 复制到可执行文件旁，避免运行时报"找不到 laszip3.dll"。
if(EXISTS "${LOCAL_DEPS_ROOT}/bin/laszip3.dll")
    add_custom_command(TARGET outlineRegularTool POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${LOCAL_DEPS_ROOT}/bin/laszip3.dll"
                "$<TARGET_FILE_DIR:outlineRegularTool>")
endif()

# Standalone mask nesting inspector (read-only diagnostic).
add_executable(inspect_mask_nesting tools/inspect_mask_nesting.cpp)
target_compile_definitions(inspect_mask_nesting PRIVATE _USE_MATH_DEFINES NOMINMAX)
if(TARGET GDAL::GDAL)
    target_link_libraries(inspect_mask_nesting PRIVATE GDAL::GDAL)
else()
    target_include_directories(inspect_mask_nesting PRIVATE ${GDAL_INCLUDE_DIR})
    target_link_libraries(inspect_mask_nesting PRIVATE ${GDAL_LIBRARY})
endif()

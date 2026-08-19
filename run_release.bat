@echo off
setlocal
cd /d "%~dp0"

set PATH=%~dp0deps\x64-windows\bin;%~dp0build_deps_release\Release;%PATH%
set OSG_LIBRARY_PATH=%~dp0deps\x64-windows\plugins\osgPlugins-3.6.5

build_deps_release\Release\outlineRegularTool.exe
pause

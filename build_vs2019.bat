@echo off
setlocal

cd /d "%~dp0"

REM Add VS 2019 bundled CMake to PATH (system PATH usually lacks cmake when double-clicking)
set "VS_CMAKE=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin"
if exist "%VS_CMAKE%\cmake.exe" set "PATH=%VS_CMAKE%;%PATH%"

where cmake >nul 2>nul || ( echo [ERROR] cmake not found. Install "C++ CMake tools for Windows" via VS Installer. & pause & exit /b 1 )

cmake -S . -B build_deps_release -G "Visual Studio 16 2019" -A x64 -T host=x64
if errorlevel 1 ( echo [ERROR] CMake configure failed. & pause & exit /b %errorlevel% )

cmake --build build_deps_release --config Release --target outlineRegularTool
if errorlevel 1 ( echo [ERROR] Build failed. & pause & exit /b %errorlevel% )

echo.
echo Build finished:
echo %cd%\build_deps_release\Release\outlineRegularTool.exe
pause

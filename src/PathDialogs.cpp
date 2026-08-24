// =============================================================================
// PathDialogs.cpp
// 作用：用 Win32 原生 IFileOpenDialog / IFileSaveDialog 实现三个路径选择对话框。
//       不新增任何第三方库；COM 库通过 #pragma comment(lib,"ole32") 链接。
// 流程：每个函数 -> 初始化 COM -> 创建对话框 -> 配置选项/标题/过滤器 ->
//       显示并取回用户所选路径(宽字符)-> 转 UTF-8 std::string 返回。
// =============================================================================

#include "PathDialogs.h"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <shobjidl.h>

#include <vector>

// 用 #pragma 链接 COM 库，无需改 CMakeLists。
#pragma comment(lib, "ole32.lib")

namespace {

// RAII：在本线程初始化/反初始化 COM。
class ScopedCOM {
public:
    ScopedCOM() : hr_(CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED)) {}
    ~ScopedCOM()
    {
        // 只有当本次调用确实增加了引用计数时才反初始化。
        // RPC_E_CHANGED_MODE 表示 COM 已被别处以 MTA 方式初始化，那种情况不动它、直接继续。
        if (hr_ == S_OK || hr_ == S_FALSE)
            CoUninitialize();
    }
private:
    HRESULT hr_;
};

// 宽字符串(以 \0 结尾) -> UTF-8 std::string
std::string WideToUtf8(const wchar_t* w)
{
    if (!w) return std::string();
    int len = WideCharToMultiByte(CP_UTF8, 0, w, -1, nullptr, 0, nullptr, nullptr);
    if (len <= 0) return std::string();
    std::vector<char> buf(static_cast<size_t>(len));
    WideCharToMultiByte(CP_UTF8, 0, w, -1, buf.data(), len, nullptr, nullptr);
    return std::string(buf.data());
}

// ===== RunDialog =====
// 作用：把已经配置好的 IFileDialog 显示出来，取回用户所选的文件系统路径。
// 参数：dlg - 已创建并配置好的对话框；out - 输出路径(UTF-8)。
// 返回：true=用户确认并选到了路径，false=取消或失败。
bool RunDialog(IFileDialog* dlg, std::string& out)
{
    // 以控制台窗口作为属主；没有控制台时返回 NULL，也允许。
    HWND owner = GetConsoleWindow();
    if (FAILED(dlg->Show(owner))) return false;

    IShellItem* item = nullptr;
    bool ok = false;
    if (SUCCEEDED(dlg->GetResult(&item)) && item) {
        PWSTR path = nullptr;
        if (SUCCEEDED(item->GetDisplayName(SIGDN_FILESYSPATH, &path)) && path) {
            out = WideToUtf8(path);
            CoTaskMemFree(path);
            ok = !out.empty();
        }
        item->Release();
    }
    return ok;
}

} // namespace

// ===== PickFolder =====
// 作用：弹出"选择文件夹"对话框，用于让用户选输入 OSGB 目录。
bool PickFolder(std::string& out)
{
    ScopedCOM com;

    IFileOpenDialog* dlg = nullptr;
    HRESULT hr = CoCreateInstance(CLSID_FileOpenDialog, nullptr, CLSCTX_INPROC_SERVER,
                                  IID_PPV_ARGS(&dlg));
    if (FAILED(hr) || !dlg) return false;

    DWORD opts = 0;
    dlg->GetOptions(&opts);
    dlg->SetOptions(opts | FOS_PICKFOLDERS | FOS_FORCEFILESYSTEM); // 只选文件夹、限定文件系统
    dlg->SetTitle(L"选择输入 OSGB 文件夹");

    bool ok = RunDialog(dlg, out);
    dlg->Release();
    return ok;
}

// ===== PickOpenFile =====
// 作用：弹出"打开文件"对话框，让用户选输入 Shapefile(*.shp)。
bool PickOpenFile(std::string& out)
{
    ScopedCOM com;

    IFileOpenDialog* dlg = nullptr;
    HRESULT hr = CoCreateInstance(CLSID_FileOpenDialog, nullptr, CLSCTX_INPROC_SERVER,
                                  IID_PPV_ARGS(&dlg));
    if (FAILED(hr) || !dlg) return false;

    const COMDLG_FILTERSPEC filters[] = {
        { L"Shapefile (*.shp)", L"*.shp" },
        { L"All files (*.*)",   L"*.*" }
    };
    dlg->SetFileTypes(2, filters);
    dlg->SetDefaultExtension(L"shp");
    dlg->SetTitle(L"选择输入 Shapefile (.shp)");
    DWORD opts = 0;
    dlg->GetOptions(&opts);
    dlg->SetOptions(opts | FOS_FILEMUSTEXIST); // 文件必须已存在

    bool ok = RunDialog(dlg, out);
    dlg->Release();
    return ok;
}

bool PickOpenXmlFile(std::string& out)
{
    ScopedCOM com;

    IFileOpenDialog* dlg = nullptr;
    HRESULT hr = CoCreateInstance(CLSID_FileOpenDialog, nullptr, CLSCTX_INPROC_SERVER,
                                  IID_PPV_ARGS(&dlg));
    if (FAILED(hr) || !dlg) return false;

    const COMDLG_FILTERSPEC filters[] = {
        { L"Metadata XML (*.xml)", L"*.xml" },
        { L"All files (*.*)",     L"*.*" }
    };
    dlg->SetFileTypes(2, filters);
    dlg->SetDefaultExtension(L"xml");
    dlg->SetTitle(L"Select metadata XML file");
    DWORD opts = 0;
    dlg->GetOptions(&opts);
    dlg->SetOptions(opts | FOS_FILEMUSTEXIST);

    bool ok = RunDialog(dlg, out);
    dlg->Release();
    return ok;
}

bool PickOpenTifFile(std::string& out)
{
    ScopedCOM com;

    IFileOpenDialog* dlg = nullptr;
    HRESULT hr = CoCreateInstance(CLSID_FileOpenDialog, nullptr, CLSCTX_INPROC_SERVER,
                                  IID_PPV_ARGS(&dlg));
    if (FAILED(hr) || !dlg) return false;

    const COMDLG_FILTERSPEC filters[] = {
        { L"GeoTIFF mask (*.tif;*.tiff)", L"*.tif;*.tiff" },
        { L"All files (*.*)",             L"*.*" }
    };
    dlg->SetFileTypes(2, filters);
    dlg->SetDefaultExtension(L"tif");
    dlg->SetTitle(L"Select AI building mask GeoTIFF");
    DWORD opts = 0;
    dlg->GetOptions(&opts);
    dlg->SetOptions(opts | FOS_FILEMUSTEXIST);

    bool ok = RunDialog(dlg, out);
    dlg->Release();
    return ok;
}

// ===== PickSaveFile =====
// 作用：弹出"保存文件"对话框，让用户选输出 Shapefile 的保存路径(默认名 regularized_building.shp)。
bool PickSaveFile(std::string& out)
{
    ScopedCOM com;

    IFileSaveDialog* dlg = nullptr;
    HRESULT hr = CoCreateInstance(CLSID_FileSaveDialog, nullptr, CLSCTX_INPROC_SERVER,
                                  IID_PPV_ARGS(&dlg));
    if (FAILED(hr) || !dlg) return false;

    const COMDLG_FILTERSPEC filters[] = {
        { L"Shapefile (*.shp)", L"*.shp" }
    };
    dlg->SetFileTypes(1, filters);
    dlg->SetDefaultExtension(L"shp");
    dlg->SetFileName(L"regularized_building.shp"); // 默认文件名
    dlg->SetTitle(L"选择输出 Shapefile 保存路径");
    DWORD opts = 0;
    dlg->GetOptions(&opts);
    dlg->SetOptions(opts | FOS_OVERWRITEPROMPT | FOS_PATHMUSTEXIST); // 覆盖时提示、路径须存在

    bool ok = RunDialog(dlg, out);
    dlg->Release();
    return ok;
}

// Dear ImGui: standalone example application for Win32 + OpenGL 3

// Learn about Dear ImGui:
// - FAQ                  https://dearimgui.com/faq
// - Getting Started      https://dearimgui.com/getting-started
// - Documentation        https://dearimgui.com/docs (same as your local docs/ folder).
// - Introduction, links and more at the top of imgui.cpp

// This is provided for completeness, however it is strongly recommended you use OpenGL with SDL or GLFW.

#define _CRT_SECURE_NO_WARNINGS

#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_win32.h"
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <GL/GL.h>
#include <tchar.h>

#include  <commdlg.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <string>
#include <map>
using std::string;
using std::map;

// Data stored per platform window
struct WGL_WindowData { HDC hDC; };

// Data
static HGLRC            g_hRC;
static WGL_WindowData   g_MainWindow;
static int              g_Width;
static int              g_Height;

#define GL_CLAMP_TO_EDGE 0x812F


// Forward declarations of helper functions
bool CreateDeviceWGL(HWND hWnd, WGL_WindowData* data);
void CleanupDeviceWGL(HWND hWnd, WGL_WindowData* data);
void ResetDeviceWGL();
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

struct texture_struct {
    GLuint* data;
    int width;
    int height;
    int channels;
};

texture_struct* LoadTextureFromFile(const char* filename)
{

    texture_struct* TS = new texture_struct();
    // Load from file
    int image_width = 0;
    int image_height = 0;
    unsigned char* image_data = stbi_load(filename, &image_width, &image_height, NULL, 4);
    if (image_data == NULL)
        return TS;

    // Create a OpenGL texture identifier
    GLuint image_texture;
    glGenTextures(1, &image_texture);
    glBindTexture(GL_TEXTURE_2D, image_texture);

    // Setup filtering parameters for display
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); // This is required on WebGL for non power-of-two textures
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE); // Same

    // Upload pixels into texture
#if defined(GL_UNPACK_ROW_LENGTH) && !defined(__EMSCRIPTEN__)
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
#endif
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image_width, image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);
    stbi_image_free(image_data);

    TS->data = new GLuint(image_texture);
    TS->width= image_width;
    TS->height= image_height;

    return TS;
}

map<string, texture_struct*> load_assets() {
    map<string, texture_struct*> lookup;

    lookup["file_splash"] = LoadTextureFromFile(".\\Assets\\splash_file.png");
    lookup["file_add"] = LoadTextureFromFile(".\\Assets\\file_add.png");
    lookup["file_load"] = LoadTextureFromFile(".\\Assets\\file_load.png");
    return lookup;
}

void resetStyling() {
    ImGui::StyleColorsClassic();
    ImGuiStyle& style = ImGui::GetStyle();

    style.WindowPadding = ImVec2(10, 10);
    style.WindowRounding = 10;

    style.Colors[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 1);
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.1, 0.1, 0.1, 1);

}

// Main code
int _main(int, char**)
{
    // Create application window
    //ImGui_ImplWin32_EnableDpiAwareness();
    WNDCLASSEXW wc = { sizeof(wc), CS_OWNDC, WndProc, 0L, 0L, GetModuleHandle(nullptr), nullptr, nullptr, nullptr, nullptr, L"ImGui Example", nullptr };
    ::RegisterClassExW(&wc);
    HWND hwnd = ::CreateWindowW(wc.lpszClassName, L"Dear ImGui Win32+OpenGL3 Example", WS_OVERLAPPEDWINDOW, 100, 100, 1280, 800, nullptr, nullptr, wc.hInstance, nullptr);

    // Initialize OpenGL
    if (!CreateDeviceWGL(hwnd, &g_MainWindow))
    {
        CleanupDeviceWGL(hwnd, &g_MainWindow);
        ::DestroyWindow(hwnd);
        ::UnregisterClassW(wc.lpszClassName, wc.hInstance);
        return 1;
    }
    wglMakeCurrent(g_MainWindow.hDC, g_hRC);

    // Show the window
    ::ShowWindow(hwnd, SW_SHOWDEFAULT);
    ::UpdateWindow(hwnd);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;   // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;    // Enable Gamepad Controls

    // Setup Dear ImGui style
    //ImGui::StyleColorsClassic();

    // Setup Platform/Renderer backends
    ImGui_ImplWin32_InitForOpenGL(hwnd);
    ImGui_ImplOpenGL3_Init();

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return a nullptr. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
    // - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use Freetype for higher quality font rendering.
    // - Read 'docs/FONTS.md' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf", 18.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, nullptr, io.Fonts->GetGlyphRangesJapanese());
    //IM_ASSERT(font != nullptr);

    resetStyling();
    ImGuiStyle& style = ImGui::GetStyle();
    ImGuiStyle StyleTextOnlyButton = ImGui::GetStyle();


    StyleTextOnlyButton.Colors[ImGuiCol_Button] = style.Colors[ImGuiCol_WindowBg];
    StyleTextOnlyButton.Colors[ImGuiCol_ButtonHovered] = ImVec4(0.3, 0.3, 0.3, 1);
    // Our state
    io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf", 18.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    io.Fonts->Build();
    ImFont* arial_font_18 = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\Arial.ttf", 18.0f, nullptr, io.Fonts->GetGlyphRangesJapanese());
    ImFont* arial_font_13 = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\Arial.ttf", 13.0f, nullptr, io.Fonts->GetGlyphRangesJapanese());
    IM_ASSERT(arial_font_18 != nullptr);
    io.Fonts->Build();

    // Our state
    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    map<string, texture_struct*> assets = load_assets();

    // Main loop
    bool done = false;

    char fileInputBuf[240];
    int resX = 1080;
    int resY = 1080;
    for (int i = 0; i < 240; i++) fileInputBuf[i] = 0;


    while (!done)
    {
        // Poll and handle messages (inputs, window resize, etc.)
        // See the WndProc() function below for our to dispatch events to the Win32 backend.
        MSG msg;
        while (::PeekMessage(&msg, nullptr, 0U, 0U, PM_REMOVE))
        {
            ::TranslateMessage(&msg);
            ::DispatchMessage(&msg);
            if (msg.message == WM_QUIT)
                done = true;
        }
        if (done)
            break;

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplWin32_NewFrame();
        ImGui::NewFrame();


        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        auto viewport = ImGui::GetMainViewport();
        auto& viewport_size = viewport->Size;

        ImGuiWindowFlags static_window =
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize;

        ImGuiWindowFlags no_scroll = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;

        bool show_welcome_window = false;

        resetStyling();
        ImGui::SetNextWindowPos(ImVec2(0, 0), 0, ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(1920, 1080), 0);
        style.Colors[ImGuiCol_WindowBg] = ImVec4(0.3, 0.3, 0.3, 1);
        ImGui::Begin("Backpane", NULL, static_window | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoBringToFrontOnFocus);
        ImGui::Text("asdfasdfasd");
        ImGui::End();
        resetStyling();
        ImGui::PushFont(arial_font_18);
        if (show_welcome_window) {
            int ypadding = 25;
            int xpadding = 15;
            style.WindowPadding = ImVec2(xpadding, ypadding);
            ImVec2 welcome_size = { 700, 500 };
            ImVec2 welcome_pos = ImVec2(((viewport_size.x - welcome_size.x) / 2.0), (viewport_size.y - welcome_size.y) / 2.0);

            
            ImGui::SetNextWindowPos(welcome_pos, 0);
            ImGui::SetNextWindowSize(welcome_size);
            ImGui::Begin("welcome window", NULL, static_window | no_scroll | ImGuiWindowFlags_NoTitleBar); {
                style.WindowPadding = ImVec2(10, 10);

                ImGui::BeginChild("##", ImVec2(200 - xpadding * 2, welcome_size.y - ypadding * 2), 0, 0); {
                    style = StyleTextOnlyButton;
                    style.ItemSpacing = ImVec2(0 ,4);
                    ImGui::Text("File Options");
                    ImGui::Separator();
                    ImGui::Image((ImTextureID)(*assets["file_add"]->data), ImVec2(15, 18),ImVec2(),ImVec2(1,1), ImVec4(0.9, 0.9, 0.9, 0.9));
                    ImGui::SameLine();
                    ImGui::SmallButton("New File");
                    ImGui::Image((ImTextureID)(*assets["file_load"]->data), ImVec2(16, 18), ImVec2(), ImVec2(1, 1), ImVec4(0.9, 0.9, 0.9, 0.9));
                    ImGui::SameLine();
                    ImGui::SmallButton("Load File");
                    ImGui::Text("Recent Files");
                    ImGui::Separator();
                    resetStyling();
                }ImGui::EndChild();
                ImGui::SameLine();
                ImDrawList* dl = ImGui::GetForegroundDrawList();
                ImVec2 p_min = ImGui::GetCursorScreenPos();
                p_min = { p_min.x + xpadding,p_min.y - ypadding };
                ImVec2 p_max = ImVec2(p_min.x + welcome_size.x - 200, p_min.y + welcome_size.y);
                dl->AddImageRounded((ImTextureID)(*assets["file_splash"]->data), p_min, p_max,
                    ImVec2(0, 0), ImVec2(1, 1), ImGui::GetColorU32(ImVec4(1, 1, 1, 1)),
                    10);
            }ImGui::End();
        }
        resetStyling();
        ImGui::SetNextWindowPos({ 0,0 });
        ImGui::SetNextWindowSize({ 250,viewport_size.y });
        ImGui::Begin("Render Settings", NULL, static_window | no_scroll | ImGuiWindowFlags_NoTitleBar);
        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File"))
            {
                if (ImGui::MenuItem("Open..", "Ctrl+O")) { /* Do stuff */ }
                if (ImGui::MenuItem("Save", "Ctrl+S")) { /* Do stuff */ }
                if (ImGui::MenuItem("Close", "Ctrl+W")) {}
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }
        ImGui::PushFont(arial_font_13);
        ImGui::Text("Model File");
        ImGui::PopFont();
        ImGui::PushItemWidth(230);
        ImGui::InputTextWithHint("##","path to file", fileInputBuf, IM_ARRAYSIZE(fileInputBuf));
        ImGui::Separator();
        style.ItemSpacing = (ImVec2(4, 2));
        ImGui::BeginChild("##Resolution", ImVec2(80, 200), 0, 0); {
            ImGui::Text("Resolution");
            ImGui::PopFont();
            ImGui::PushFont(arial_font_13);
            ImGui::Text("X: ");
            ImGui::SameLine();
            ImGui::InputInt("##resx", &resX, 0);
            ImGui::Text("Y: ");
            ImGui::SameLine();
            ImGui::InputInt("##resy", &resY, 0);
        } ImGui::EndChild();
        ImGui::End();
        // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
        
        ImGui::PopFont();
        // Rendering
        ImGui::Render();
        glViewport(0, 0, g_Width, g_Height);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // Present
        ::SwapBuffers(g_MainWindow.hDC);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplWin32_Shutdown();
    ImGui::DestroyContext();

    CleanupDeviceWGL(hwnd, &g_MainWindow);
    wglDeleteContext(g_hRC);
    ::DestroyWindow(hwnd);
    ::UnregisterClassW(wc.lpszClassName, wc.hInstance);

    return 0;
}

// Helper functions
bool CreateDeviceWGL(HWND hWnd, WGL_WindowData* data)
{
    HDC hDc = ::GetDC(hWnd);
    PIXELFORMATDESCRIPTOR pfd = { 0 };
    pfd.nSize = sizeof(pfd);
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 32;

    const int pf = ::ChoosePixelFormat(hDc, &pfd);
    if (pf == 0)
        return false;
    if (::SetPixelFormat(hDc, pf, &pfd) == FALSE)
        return false;
    ::ReleaseDC(hWnd, hDc);

    data->hDC = ::GetDC(hWnd);
    if (!g_hRC)
        g_hRC = wglCreateContext(data->hDC);
    return true;
}

void CleanupDeviceWGL(HWND hWnd, WGL_WindowData* data)
{
    wglMakeCurrent(nullptr, nullptr);
    ::ReleaseDC(hWnd, data->hDC);
}

// Forward declare message handler from imgui_impl_win32.cpp
extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

// Win32 message handler
// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
LRESULT WINAPI WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    if (ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam))
        return true;

    switch (msg)
    {
    case WM_SIZE:
        if (wParam != SIZE_MINIMIZED)
        {
            g_Width = LOWORD(lParam);
            g_Height = HIWORD(lParam);
        }
        return 0;
    case WM_SYSCOMMAND:
        if ((wParam & 0xfff0) == SC_KEYMENU) // Disable ALT application menu
            return 0;
        break;
    case WM_DESTROY:
        ::PostQuitMessage(0);
        return 0;
    }
    return ::DefWindowProcW(hWnd, msg, wParam, lParam);
}

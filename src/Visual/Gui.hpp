#ifndef MEDYAN_Visual_Gui_hpp
#define MEDYAN_Visual_Gui_hpp

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include "Visual/DisplaySettings.hpp"

namespace medyan::visual {


// Note:
//   - This must be used while OpenGL and GLFW environments are live
class ImguiGuard {

private:
    ImGuiContext* context_ = nullptr;

public:
    ImguiGuard(GLFWwindow* window) {
        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        context_ = ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO(); (void)io;
        //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
        //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

        // Setup Dear ImGui style
        ImGui::StyleColorsDark();
        //ImGui::StyleColorsClassic();

        // Setup Platform/Renderer bindings
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330 core");

        // Load Fonts
        // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
        // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
        // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
        // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
        // - Read 'docs/FONTS.md' for more instructions and details.
        // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
        //io.Fonts->AddFontDefault();
        //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
        //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
        //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
        //io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
        //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
        //IM_ASSERT(font != NULL);
    }

    ~ImguiGuard() {
        // Cleanup
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext(context_);

    }
};



inline void guiColorPicker4Popup(const char* strId, float* pColor) {
    const bool bgColorPopup = ImGui::ColorButton(
        strId,
        [&] { ImVec4 to; std::memcpy(&to, pColor, sizeof(to)); return to; }(),
        0
    );
    ImGui::SameLine();
    ImGui::Text(strId);

    if(bgColorPopup) {
        ImGui::OpenPopup(strId);
    }

    if(ImGui::BeginPopup(strId)) {
        ImGui::ColorPicker4(strId, pColor, ImGuiColorEditFlags_PickerHueWheel);

        ImGui::EndPopup();
    }
}

inline void guiMainWindow(DisplaySettings& displaySettings) {
    // Exceptionally add an extra assert here for people confused about initial Dear ImGui setup
    // Most ImGui functions would normally just crash if the context is missing.
    IM_ASSERT(ImGui::GetCurrentContext() != NULL && "Missing dear imgui context.");

    ImGuiWindowFlags windowFlags =
        ImGuiWindowFlags_MenuBar |
        ImGuiWindowFlags_NoCollapse;

    // We specify a default position/size in case there's no data in the .ini file.
    // We only do it to make the demo applications a little more welcoming, but typically this isn't required.
    // ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(500, 700), ImGuiCond_FirstUseEver);

    // Main body of the Demo window starts here.
    if (!ImGui::Begin("medyan control", nullptr, windowFlags))
    {
        // Early out if the window is collapsed, as an optimization.
        ImGui::End();
        return;
    }

    // Most "big" widgets share a common width settings by default. See 'Demo->Layout->Widgets Width' for details.

    // e.g. Use 2/3 of the space for widgets and 1/3 for labels (default)
    //ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.65f);

    // e.g. Leave a fixed amount of width for labels (by passing a negative value), the rest goes to widgets.
    ImGui::PushItemWidth(ImGui::GetFontSize() * -12);

    // Menu Bar
    if (ImGui::BeginMenuBar())
    {
        if (ImGui::BeginMenu("file"))
        {
            ImGui::MenuItem("(empty menu)", NULL, false, false);

            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    // Main display settings
    if(ImGui::BeginCombo("display mode", text(displaySettings.displayMode), 0)) {
        for (int i = 0; i < underlying(DisplayMode::last_); ++i) {

            const DisplayMode mode = static_cast<DisplayMode>(i);

            const bool isSelected = (displaySettings.displayMode == mode);
            if (ImGui::Selectable(text(mode), isSelected)) {
                if(!isSelected) {
                    // Do stuff to switch to new display mode
                }
                displaySettings.displayMode = mode;
            }

            // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
            if (isSelected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }

    ImGui::Spacing();


    if (ImGui::CollapsingHeader("info")) {

        // Function to print view object information
        const auto printView = [](const ObjectViewSettings& viewSettings) {
            ImGui::Text("view");
            ImGui::BulletText(
                "camera position: (%.1f, %.1f, %.1f)",
                viewSettings.camera.position[0],
                viewSettings.camera.position[1],
                viewSettings.camera.position[2]
            );
            ImGui::BulletText(
                "camera target: (%.1f, %.1f, %.1f)",
                viewSettings.camera.target[0],
                viewSettings.camera.target[1],
                viewSettings.camera.target[2]
            );
            ImGui::BulletText("projection: %s", text(viewSettings.projection.type));
            if(viewSettings.projection.type == ObjectViewSettings::Projection::Type::perspective) {
                ImGui::BulletText("fov: %.1f", viewSettings.projection.fov);
            }
            ImGui::BulletText(
                "z range: (%.1f, %.1f)",
                viewSettings.projection.zNear,
                viewSettings.projection.zFar
            );
        };

        printView(displaySettings.mainView);
        // ImGui::Separator();
    }

    if(ImGui::CollapsingHeader("settings")) {

        guiColorPicker4Popup("background", displaySettings.mainView.canvas.bgColor.data());
        ImGui::Separator();

        ImGui::SliderFloat("camera key speed", &displaySettings.mainView.control.cameraKeyPositionPerFrame, 50.0, 500.0, "%.1f");
    }


    // End of main window
    ImGui::End();
}

// Note:
//   - This should only be used in GLFW main loop
inline void imguiLoopRender(DisplaySettings& displaySettings) {
    ImGui_ImplGlfw_NewFrame();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui::NewFrame();

    if(displaySettings.gui.enabled) {
        guiMainWindow(displaySettings);
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

} // namespace medyan::visual

#endif

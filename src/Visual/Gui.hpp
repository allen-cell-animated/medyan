#ifndef MEDYAN_Visual_Gui_hpp
#define MEDYAN_Visual_Gui_hpp

#include <cstdio>
#include <type_traits>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <imgui_internal.h>
#include <imgui_stdlib.h>
#include <nfd.h>

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"
#include "Visual/FrameData.hpp"

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


class ImguiDisableHelperGuard {
    bool yesPleaseReallyDisableItForReal_;

public:
    ImguiDisableHelperGuard(bool forReal) : yesPleaseReallyDisableItForReal_(forReal) {
        if(yesPleaseReallyDisableItForReal_) {
            ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
            ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
        }
    }

    ~ImguiDisableHelperGuard() {
        if(yesPleaseReallyDisableItForReal_) {
            ImGui::PopItemFlag();
            ImGui::PopStyleVar();
        }
    }
};


//-----------------------------------------------------------------------------
// Auxiliary functions
//-----------------------------------------------------------------------------

// Returns whether the color is changed.
template<
    int colorDof,
    std::enable_if_t< colorDof == 3 || colorDof == 4 >* = nullptr
>
inline bool guiAuxColorPickerPopup(
    const char*        strId, 
    float*             pColor
) {
    bool changed = false;

    const bool bgColorPopup = ImGui::ColorButton(
        strId,
        ImVec4(
            pColor[0],
            pColor[1],
            pColor[2],
            colorDof == 3 ? 1.0f : pColor[3]
        ),
        0
    );
    ImGui::SameLine();
    ImGui::Text(strId);

    if(bgColorPopup) {
        ImGui::OpenPopup(strId);
    }

    if(ImGui::BeginPopup(strId)) {
        if(colorDof == 3) {
            changed = ImGui::ColorPicker3(strId, pColor, ImGuiColorEditFlags_PickerHueWheel);
        } else {
            changed = ImGui::ColorPicker4(strId, pColor, ImGuiColorEditFlags_PickerHueWheel);
        }

        ImGui::EndPopup();
    }

    return changed;
}

inline bool guiAuxColorPicker3Popup(const char* strId, float* pColor) {
    return guiAuxColorPickerPopup<3>(strId, pColor);
}
inline bool guiAuxColorPicker4Popup(const char* strId, float* pColor) {
    return guiAuxColorPickerPopup<4>(strId, pColor);
}

// Function to build combo box automatically for an enumeration type.
//
// Returns whether a new value is selected.
//
// Notes:
//   - The elements in the enum type must be automatically valued (ie the
//     values are automatically 0, 1, 2, ...).
//   - The type must have "last_" as the last element.
//   - A "text(Enum)" function must be implemented to display the elements.
template<
    typename Enum,
    typename Reselect,              // function void(Enum old, Enum new) to execute when selected
    std::enable_if_t<
        std::is_enum_v< Enum > &&
        std::is_invocable_r_v< void, Reselect, Enum, Enum > // Reselect: (Enum, Enum) -> void
    >* = nullptr                    // type requirements
>
inline bool guiAuxEnumComboBox(
    const char*          name,
    Enum&                value,
    Reselect&&           reselect,
    ImGuiSelectableFlags flags = 0
) {
    bool changed = false;

    if(ImGui::BeginCombo(name, text(value), 0)) {
        for (int i = 0; i < underlying(Enum::last_); ++i) {

            const Enum valueI = static_cast<Enum>(i);

            const bool isSelected = (value == valueI);
            if (ImGui::Selectable(text(valueI), isSelected, flags)) {
                const auto oldValue = value;
                value = valueI;
                reselect(oldValue, valueI);

                if(!isSelected) {
                    // selected one item that was previously not selected
                    changed = true;
                }
            }

            // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
            if (isSelected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }

    return changed;
}
// Where Reselect function is a no-op.
template<
    typename Enum,
    std::enable_if_t< std::is_enum_v< Enum > >* = nullptr     // type requirements
>
inline bool guiAuxEnumComboBox(
    const char*          name,
    Enum&                value,
    ImGuiSelectableFlags flags = 0
) {
    return guiAuxEnumComboBox(name, value, [](Enum, Enum) {}, flags);
}

// Tristate checkbox
//
// State: 0 -> unchecked, 1 -> checked, -1 -> intermediate
// Clicking the checkbox will never set the state to -1.
//
// adapted from https://github.com/ocornut/imgui/issues/2644#issuecomment-507023896
inline bool guiAuxCheckboxTristate(const char* label, int* pState) {
    bool clicked = false;
    if (*pState == -1)
    {
        ImGui::PushItemFlag(ImGuiItemFlags_MixedValue, true);
        bool b = false;
        clicked = ImGui::Checkbox(label, &b);
        if (clicked)
            *pState = 1;
        ImGui::PopItemFlag();
    }
    else
    {
        bool b = (*pState != 0);
        clicked = ImGui::Checkbox(label, &b);
        if (clicked)
            *pState = (int)b;
    }
    return clicked;
}

// check box to toggle flags
//
// Returns whether the flag is changed.
template< typename UInt >
inline bool guiAuxCheckboxFlags(const char* label, UInt& flag, UInt setTo) {
    const UInt oldFlag     = flag;
    const bool atLeast1Set = flag & setTo;
    const bool allSet      = (flag & setTo) == setTo;
    int        state       = allSet ? 1 : (atLeast1Set ? -1 : 0);

    if(guiAuxCheckboxTristate(label, &state)) {
        if(state == 1) {
            flag |= setTo;
        } else if(state == 0) {
            flag &= ~setTo;
        } // cannot be -1
    }

    return flag != oldFlag;
}

// select file
inline void guiAuxSelectFile(std::filesystem::path& file) {
    nfdchar_t* filepath = nullptr;
    auto result = NFD_OpenDialog(nullptr, std::filesystem::current_path().string().c_str(), &filepath);

    if ( result == NFD_OKAY ) {
        file = filepath;
        std::free(filepath);
    }
    else if ( result == NFD_CANCEL ) {
        // do nothing
    }
    else {
        LOG(ERROR) << "Error: " << NFD_GetError();
    }
}

//-----------------------------------------------------------------------------
// Main GUI functions
//-----------------------------------------------------------------------------

inline void guiHelpWindow(DisplaySettings& displaySettings) {
    if(!displaySettings.gui.helpWindow) return;

    ImGui::SetNextWindowSize(ImVec2(360, 400), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("help", &displaySettings.gui.helpWindow)) {

        if(ImGui::CollapsingHeader("controls", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Text("key controls:");
            ImGui::BulletText("f: take snapshot of current screen");
            ImGui::BulletText("g: toggle gui on/off");
            ImGui::BulletText("w/a/s/d: move camera horizontally");
        }
    }


    // End of help window, no matter Begin returns true or false.
    ImGui::End();
}

// Functions for profile configuration
// Returns whether the profile is changed.
inline bool guiGeometryDisplaySettings(SurfaceDisplaySettings& settings) {
    bool changed = false;

    changed |= ImGui::Checkbox("enabled", &settings.enabled);
    changed |= guiAuxEnumComboBox("polygon", settings.polygonMode);

    return changed;
}
inline bool guiGeometryDisplaySettings(LineDisplaySettings& settings) {
    bool changed = false;

    changed |= ImGui::Checkbox("enabled", &settings.enabled);

    return changed;
}

inline bool guiActiveProfileConfig(MembraneProfile& profile) {
    bool changed = false;

    changed |= guiGeometryDisplaySettings(profile.displaySettings.surface);

    changed |= guiAuxColorPicker3Popup("fixed color", profile.displaySettings.colorFixed.value.data());

    return changed;
}
inline bool guiActiveProfileConfig(FilamentProfile& profile) {
    bool changed = false;

    // path mode specific
    changed |= guiAuxEnumComboBox("path mode", profile.displaySettings.pathMode);

    switch(profile.displaySettings.pathMode) {
        case FilamentDisplaySettings::PathMode::line:
            break;
        case FilamentDisplaySettings::PathMode::extrude:
            changed |= ImGui::SliderFloat("extrude radius", &profile.displaySettings.pathExtrudeRadius, 1.0f, 24.0f, "%.1f");
            changed |= ImGui::SliderInt("extrude sides", &profile.displaySettings.pathExtrudeSides, 4, 20);
            break;
        case FilamentDisplaySettings::PathMode::bead:
            changed |= ImGui::SliderFloat("bead radius", &profile.displaySettings.beadRadius, 1.0f, 50.0f, "%.1f");
            break;
    }

    if(displayGeometryType(profile.displaySettings) == DisplayGeometryType::line) {
        changed |= guiGeometryDisplaySettings(profile.displaySettings.line);
    } else {
        changed |= guiGeometryDisplaySettings(profile.displaySettings.surface);
    }

    changed |= guiAuxColorPicker3Popup("fixed color", profile.displaySettings.colorFixed.value.data());

    return changed;
}
inline bool guiActiveProfileConfig(LinkerProfile& profile) {
    bool changed = false;

    // selector
    {
        // on/off checkbox
        bool filterOn = profile.selector.name.has_value();
        const bool filterOnChanged = ImGui::Checkbox(filterOn ? "##filter check" : "filter##filter check", &filterOn);
        changed |= filterOnChanged;

        if(filterOnChanged) {
            if(filterOn) profile.selector.name.emplace();
            else         profile.selector.name.reset();
        }

        // show text box
        if(filterOn) {
            ImGui::SameLine();
            changed |= ImGui::InputText("filter##filter text", &*profile.selector.name);
        }
    }

    // path mode specific
    changed |= guiAuxEnumComboBox("path mode", profile.displaySettings.pathMode);

    switch(profile.displaySettings.pathMode) {
        case LinkerDisplaySettings::PathMode::line:
            break;
        case LinkerDisplaySettings::PathMode::extrude:
            changed |= ImGui::SliderFloat("extrude radius", &profile.displaySettings.pathExtrudeRadius, 1.0f, 24.0f, "%.1f");
            changed |= ImGui::SliderInt("extrude sides", &profile.displaySettings.pathExtrudeSides, 4, 20);
            break;
    }

    if(displayGeometryType(profile.displaySettings) == DisplayGeometryType::line) {
        changed |= guiGeometryDisplaySettings(profile.displaySettings.line);
    } else {
        changed |= guiGeometryDisplaySettings(profile.displaySettings.surface);
    }

    // color
    changed |= guiAuxColorPicker3Popup("fixed color", profile.displaySettings.colorFixed.value.data());

    return changed;
}
inline bool guiActiveProfileConfig(AuxLineProfile& profile) {
    bool changed = false;

    changed |= guiAuxCheckboxFlags("compartment border", profile.displaySettings.flag, AuxLineDisplaySettings::targetCompartmentBorder);
    changed |= guiAuxCheckboxFlags("compartment grid", profile.displaySettings.flag, AuxLineDisplaySettings::targetCompartmentAll);

    changed |= guiGeometryDisplaySettings(profile.displaySettings.line);

    return changed;
}

// Returns whether the profile is changed
inline bool guiActiveProfileConfig(ElementProfile& profile) {
    bool changed = false;

    // Combo box to select target
    if(ImGui::BeginCombo("target", elementProfileDisplayName(profile.index()), 0)) {
        for (int i = 0; i < std::variant_size_v< ElementProfile >; ++i) {

            const bool isSelected = (profile.index() == i);
            if (ImGui::Selectable(elementProfileDisplayName(i), isSelected, 0)) {
                if(!isSelected) {
                    // A new profile is chosen
                    profile = makeProfileWithIndex(i);

                    // Because we selected a profile that was not previously selected
                    changed = true;
                }
            }

            // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
            if (isSelected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }

    // Display actual settings
    changed |= std::visit([](auto& actualProfile) { return guiActiveProfileConfig(actualProfile); }, profile);

    return changed;
}

inline void guiProfileWindow(
    DisplaySettings& displaySettings,
    DisplayStates&   displayStates
) {
    if(!displaySettings.gui.profileWindow) return;

    ImGui::SetNextWindowSize(ImVec2(750, 300), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("profile", &displaySettings.gui.profileWindow)) {

        static int selectedTrajectoryIndex = -1; // used only in trajectory mode. -1 means not selected
        static int selectedProfileIndex = -1; // used only when profiles are listed. -1 means not selected

        // left pane, for trajectory list
        {
            ImGui::BeginChild("trajectory pane", ImVec2(200, 0), true);

            if(displaySettings.displayMode == DisplayMode::trajectory) {
                if(ImGui::CollapsingHeader("new file", ImGuiTreeNodeFlags_Framed)) {

                    static DisplayTrajectoryFileSettings newFile;

                    const auto fileSelection = [&](std::filesystem::path& file, const char* buttonName) {
                        if(ImGui::Button(buttonName)) {
                            guiAuxSelectFile(file);
                        }
                        if(!file.empty()) {
                            ImGui::SameLine();
                            char clearButtonLabelId[128];
                            std::snprintf(clearButtonLabelId, 128, "x##clear %s", buttonName);
                            if(ImGui::Button(clearButtonLabelId)) {
                                file.clear();
                            }

                            ImGui::TextWrapped(file.string().c_str());
                        }
                    };

                    fileSelection(newFile.trajSnapshot, "snapshot.traj");

                    // The "add" button
                    {
                        ImguiDisableHelperGuard guard(displayStates.sync.busy.load());

                        if(ImGui::Button("add")) {
                            pushAnAsyncTask(
                                displayStates.sync,
                                // file information is copied but not referenced
                                [&, file = newFile] {
                                    backgroundTaskReadTrajectory(displayStates.sync, file);
                                }
                            );
                        }
                    }
                }
            }

            if(displaySettings.displayMode == DisplayMode::realtime) {
                // always select realtime
                ImGui::Selectable("realtime", true);
            }
            else {
                const int numTrajectories = displayStates.trajectoryDataStates.trajectories.size();
                for (int i = 0; i < numTrajectories; ++i)
                {
                    auto& thisTraj = displayStates.trajectoryDataStates.trajectories[i];

                    char label[64];
                    std::snprintf(label, 64, "%d %s", i, thisTraj.displayMasterSwitch ? "" : "disabled");
                    if (ImGui::Selectable(label, selectedTrajectoryIndex == i)) {
                        selectedTrajectoryIndex = i;
                    }

                    if(ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                        thisTraj.displayMasterSwitch ^= true;
                        if(thisTraj.displayMasterSwitch) {
                            // mark mesh data to be updated
                            for(auto& profileData : thisTraj.profileData) {
                                profileData.shouldUpdateMeshData = true;
                            }
                        }
                    }
                }
            }

            ImGui::EndChild();
        }
        ImGui::SameLine();

        // middle pane, for trajectory configuration and profiles
        {
            ImGui::BeginChild("config pane", ImVec2(250, 0), true);

            // The function to display all profiles
            const auto displayProfiles = [&](std::vector< ProfileWithMeshData >& vecProfileData) {

                // add profiles
                if(ImGui::Button("add profile")) {
                    vecProfileData.emplace_back();
                }

                const int numProfiles = vecProfileData.size();

                for(int i = 0; i < numProfiles; ++i) {
                    char label[128];
                    std::snprintf(label, 128, "%d %s", i, elementProfileDisplayName(vecProfileData[i].profile.index()));
                    if(ImGui::Selectable(label, selectedProfileIndex == i)) {
                        selectedProfileIndex = i;
                    }
                }
            };

            if(displaySettings.displayMode == DisplayMode::trajectory) {
                auto& trajectories = displayStates.trajectoryDataStates.trajectories;
                selectedTrajectoryIndex = std::min(selectedTrajectoryIndex, (int)trajectories.size() - 1);

                if(selectedTrajectoryIndex >= 0) {
                    // files
                    if(ImGui::CollapsingHeader("file")) {
                        ImGui::TextWrapped("snapshot file: %s", trajectories[selectedTrajectoryIndex].inputs.trajSnapshot.string().c_str());
                        if(ImGui::Button("close")) {
                            // close this file
                            trajectories.erase(trajectories.begin() + selectedTrajectoryIndex);
                        }
                    }

                    // profiles
                    if(ImGui::CollapsingHeader("profiles", ImGuiTreeNodeFlags_DefaultOpen)) {
                        if(selectedTrajectoryIndex < displayStates.trajectoryDataStates.trajectories.size()) {
                            displayProfiles(displayStates.trajectoryDataStates.trajectories[selectedTrajectoryIndex].profileData);
                        }
                    }
                }
            }
            else {
                if(ImGui::CollapsingHeader("profiles", ImGuiTreeNodeFlags_DefaultOpen)) {
                    displayProfiles(displayStates.realtimeDataStates.profileData);
                }
            }

            ImGui::EndChild();
        }
        ImGui::SameLine();

        // right pane, for profile configuration
        {
            ImGui::BeginChild("profile pane", ImVec2(0, 0), true);

            if(displaySettings.displayMode == DisplayMode::trajectory) {
                const bool valid =
                    selectedTrajectoryIndex >= 0 && selectedTrajectoryIndex < displayStates.trajectoryDataStates.trajectories.size() &&
                    selectedProfileIndex    >= 0 && selectedProfileIndex    < displayStates.trajectoryDataStates.trajectories[selectedTrajectoryIndex].profileData.size();

                if(valid) {
                    auto& traj = displayStates.trajectoryDataStates.trajectories[selectedTrajectoryIndex];
                    if(ImGui::Button("remove profile")) {
                        traj.profileData.erase(traj.profileData.begin() + selectedProfileIndex);
                    }
                    else {
                        auto& profileData = traj.profileData[selectedProfileIndex];
                        if(guiActiveProfileConfig(profileData.profile)) {
                            profileData.shouldUpdateMeshData = true;
                        }
                    }
                }
            }
            else {
                const bool valid =
                    selectedProfileIndex >= 0 && selectedProfileIndex < displayStates.realtimeDataStates.profileData.size();

                if(valid) {
                    if(ImGui::Button("remove profile")) {
                        displayStates.realtimeDataStates.profileData.erase(displayStates.realtimeDataStates.profileData.begin() + selectedProfileIndex);
                    }
                    else {
                        auto& profileData = displayStates.realtimeDataStates.profileData[selectedProfileIndex];
                        if(guiActiveProfileConfig(profileData.profile)) {
                            profileData.shouldUpdateMeshData = true;
                        }
                    }
                }
            }

            ImGui::EndChild();
        }
    }

    // End of help window, no matter Begin returns true or false.
    ImGui::End();
}

inline void guiTrajectoryControlPanel(
    DisplaySettings& displaySettings,
    DisplayStates&   displayStates
) {
    // Playback
    ImGui::Checkbox("play", &displayStates.playback.isPlaying);
    ImGui::SliderFloat("speed", &displaySettings.playback.fps, 1, 24, "%.1f");
    ImGui::SliderInt("playback", &displayStates.playback.currentFrame, 0, displayStates.playback.maxFrame);
}

inline void guiViewSettings(ObjectViewSettings& viewSettings) {
    using PT = ObjectViewSettings::Projection::Type;

    guiAuxColorPicker4Popup("background", viewSettings.canvas.bgColor.data());
    guiAuxEnumComboBox("projection", viewSettings.projection.type);

    ImGui::SliderFloat("z near", &viewSettings.projection.zNear, 1.0, 200.0, "%.1f");
    ImGui::SliderFloat("z far", &viewSettings.projection.zFar, 1000.0, 20000.0, "%.1f");
    ImGui::SliderFloat("camera key speed", &viewSettings.control.cameraKeyPositionPerFrame, 50.0, 500.0, "%.1f");

    guiAuxEnumComboBox("mouse mode", viewSettings.control.cameraMouseMode);

}

inline void guiMainWindow(
    DisplaySettings& displaySettings,
    DisplayStates  & displayStates
) {
    const bool busy = displayStates.sync.busy.load();

    // Exceptionally add an extra assert here for people confused about initial Dear ImGui setup
    // Most ImGui functions would normally just crash if the context is missing.
    IM_ASSERT(ImGui::GetCurrentContext() != NULL && "Missing dear imgui context.");

    ImGuiWindowFlags windowFlags =
        ImGuiWindowFlags_MenuBar |
        ImGuiWindowFlags_NoCollapse;

    // We specify a default position/size in case there's no data in the .ini file.
    // We only do it to make the demo applications a little more welcoming, but typically this isn't required.
    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(400, 600), ImGuiCond_FirstUseEver);

    // Main body of the Demo window starts here.
    if (!ImGui::Begin("medyan control", nullptr, windowFlags))
    {
        // Early out if the window is collapsed, as an optimization.
        ImGui::End();
        return;
    }

    // Most "big" widgets share a common width settings by default. See 'Demo->Layout->Widgets Width' for details.

    // e.g. Use 2/3 of the space for widgets and 1/3 for labels (default)
    // ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.65f);

    // e.g. Leave a fixed amount of width for labels (by passing a negative value), the rest goes to widgets.
    // ImGui::PushItemWidth(ImGui::GetFontSize() * -12);

    // Menu Bar
    if (ImGui::BeginMenuBar())
    {
        if (ImGui::BeginMenu("file")) {
            ImGui::MenuItem("(empty menu)", NULL, false, false);
            ImGui::EndMenu();
        }

        if(ImGui::BeginMenu("window")) {
            ImGui::MenuItem("profile", nullptr, &displaySettings.gui.profileWindow, true);
            ImGui::EndMenu();
        }

        if(ImGui::BeginMenu("help"))
        {
            ImGui::MenuItem("help window", nullptr, &displaySettings.gui.helpWindow, true);
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    // Main display settings
    guiAuxEnumComboBox(
        "display mode",
        displaySettings.displayMode,
        [&](DisplayMode modeOld, DisplayMode modeNew) {
            switch(modeNew) {
                case DisplayMode::trajectory:
                    if(modeOld != modeNew && displayStates.trajectoryDataStates.trajectories.empty()) {
                        pushAnAsyncTask(
                            displayStates.sync,
                            [&] {
                                backgroundTaskReadTrajectory(displayStates.sync, DisplayTrajectoryFileSettings {});
                            }
                        );
                    }
                    break;
            }
        },
        busy ? ImGuiSelectableFlags_Disabled : 0
    );

    if(displaySettings.displayMode == DisplayMode::trajectory) {
        guiTrajectoryControlPanel(displaySettings, displayStates);
    }

    ImGui::Spacing();


    if (ImGui::CollapsingHeader("info", ImGuiTreeNodeFlags_DefaultOpen)) {

        // Function to print view object information
        const auto printView = [](const ObjectViewSettings& viewSettings) {
            ImGui::Text("view");
            ImGui::BulletText(
                "window size: (%d, %d)",
                viewSettings.canvas.width,
                viewSettings.canvas.height
            );
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
                ImGui::BulletText("fov: %.2f", viewSettings.projection.fov);
            } else {
                ImGui::BulletText("scale: %.1f", viewSettings.projection.scale);
            }
            ImGui::BulletText(
                "z range: (%.1f, %.1f)",
                viewSettings.projection.zNear,
                viewSettings.projection.zFar
            );
        };

        // busy information
        if(busy) {
            ImGui::Text("busy...");
            ImGui::Separator();
        }
        // fps information
        ImGui::Text("fps: %.1f", displayStates.timing.fps);
        ImGui::Separator();
        // main view information
        printView(displaySettings.mainView);

    }

    if(ImGui::CollapsingHeader("settings")) {

        // main view
        guiViewSettings(displaySettings.mainView);
    }


    // End of main window
    ImGui::End();
}

// Note:
//   - This should only be used in GLFW main loop
inline void imguiLoopRender(
    DisplaySettings& displaySettings,
    DisplayStates  & displayStates
) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if(displaySettings.gui.enabled) {
        guiMainWindow(displaySettings, displayStates);
        guiProfileWindow(displaySettings, displayStates);
        guiHelpWindow(displaySettings);
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

} // namespace medyan::visual

#endif

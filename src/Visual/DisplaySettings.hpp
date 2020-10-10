#ifndef MEDYAN_Visual_DisplaySettings_hpp
#define MEDYAN_Visual_DisplaySettings_hpp

#include <array>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "Util/Io/Log.hpp"

namespace medyan::visual {

//-------------------------------------
// Global display mode
//-------------------------------------

enum class DisplayMode {
    realtime,        // Display system data from memory during simulation
    trajectory,      // Read and display trajectory from output files

    // Reflection hack for counting elements
    last_
};

constexpr auto text(DisplayMode mode) {
    switch(mode) {
        case DisplayMode::realtime:   return "realtime";
        case DisplayMode::trajectory: return "trajectory";
        default:                      return "";
    }
}

//-------------------------------------
// GUI settings
//-------------------------------------

struct GuiSettings {
    // the master switch
    bool enabled = false;

    // windows switches
    bool helpWindow = false;
    bool profileWindow = true;
};

//-------------------------------------
// Graphics settings
//-------------------------------------

struct ObjectViewSettings {

    struct Canvas {
        std::array< float, 4 > bgColor { 0.0f, 0.0f, 0.0f, 0.0f };

        int width = 1200;
        int height = 800;
    };

    // The Camera struct stores the states and parameters of a camera, used in view
    // transformation.
    struct Camera {
        glm::mat4 view() const {
            return glm::lookAt(position, target, up);
        }

        // Positions
        glm::vec3 position = glm::vec3(0.0f, 0.0f, 2000.0f);
        glm::vec3 target   = glm::vec3(0.0f, 0.0f, 0.0f);
        glm::vec3 right    = glm::vec3(1.0f, 0.0f, 0.0f);
        glm::vec3 up       = glm::vec3(0.0f, 1.0f, 0.0f);
    };

    struct Projection {
        enum class Type { orthographic, perspective, last_ };

        Type  type       = Type::perspective;
        float fov        = glm::radians(45.0f); // perspective
        float scale      = 1.0f;                // orthographic
        float zNear      = 10.0f;
        float zFar       = 5000.0f;

        // Helper functions

        auto projOrtho(float width, float height) const {
            return glm::ortho(-width * 0.5f * scale, width * 0.5f * scale, -height * 0.5f * scale, height * 0.5f * scale, zNear, zFar);
        }
        auto projPerspective(float width, float height) const {
            return glm::perspective(fov, width / height, zNear, zFar);
        }
        auto proj(float width, float height) const {
            return type == Type::orthographic ?
                projOrtho      (width, height):
                projPerspective(width, height);
        }
    };

    struct Control {
        enum class CameraMouseMode {
            none,
            rotate,
            pan,
            last_
        };
        enum class SnapshotResolution {
            scaleWithScreen,
            manual,
            last_
        };

        // camera
        //-----------------------------
        float cameraKeyPositionPerFrame = 200.0f;

        CameraMouseMode cameraMouseMode = CameraMouseMode::rotate;
        float cameraRotatePositionPerCursorPixel = 2.0f;   // camera move distance per pixel (rotate)
        float cameraPanPositionPerCursorPixel = 2.0f;      // camera move distance per pixel (pan)

        // snapshot saving section
        //-----------------------------
        std::filesystem::path snapshotSavePath = ".";
        std::string           snapshotBaseName = "snapshot";
        SnapshotResolution    snapshotResolution = SnapshotResolution::scaleWithScreen;
        float                 snapshotScale = 1.0f;   // Used with "scale with screen"
        // Used with "scale with screen".
        // If set to true, it will undo the scaling when calculating orthographic projection.
        bool                  snapshotUndoScaleOrtho = true;
        int                   snapshotWidth = 1920;   // Used with "manual"
        int                   snapshotHeight = 1080;   // Used with "manual"

        // This function gives the snapshot size in integers (width, height).
        auto snapshotSize(int screenViewWidth, int screenViewHeight) const {
            int width, height;
            if(snapshotResolution == SnapshotResolution::scaleWithScreen) {
                width = screenViewWidth * snapshotScale;
                height = screenViewHeight * snapshotScale;
            } else {
                width = snapshotWidth;
                height = snapshotHeight;
            }

            return std::pair { width, height };
        }

        // Finds the snapshot file name.
        //
        // The file name will be displayed as "baseName-index.ext", where the
        // index will be displayed with padding zeros.
        //
        // The index should be -1 or non-negative, where -1 will not append
        // "-index" after the base name.
        auto snapshotPngFile(int index = -1) const {
            if(index == -1) {
                return snapshotSavePath / (snapshotBaseName + ".png");
            }
            else {
                if(index < 0 || index > 999999) {
                    LOG(ERROR) << "Snapshot index " << index << " out of range.";
                    throw std::runtime_error("Snapshot index out of range");
                }
                // 11 chars after base name: "-012345.png"
                char appendStr[12] {};
                std::snprintf(appendStr, 12, "-%06d.png", index);
                return snapshotSavePath / (snapshotBaseName + appendStr);
            }
        }
    };

    // Canvas
    Canvas canvas;

    // View
    Camera camera;

    // Projection
    Projection projection;

    // Control
    Control control;
};

constexpr auto text(ObjectViewSettings::Projection::Type proj) {
    switch(proj) {
        case ObjectViewSettings::Projection::Type::orthographic: return "orthographic";
        case ObjectViewSettings::Projection::Type::perspective:  return "perspective";
        default:                                                 return "";
    }
}
constexpr auto text(ObjectViewSettings::Control::CameraMouseMode mode) {
    switch(mode) {
        case ObjectViewSettings::Control::CameraMouseMode::none:   return "none";
        case ObjectViewSettings::Control::CameraMouseMode::rotate: return "rotate";
        case ObjectViewSettings::Control::CameraMouseMode::pan:    return "pan";
        default:                                                   return "";
    }
}
constexpr auto text(ObjectViewSettings::Control::SnapshotResolution value) {
    switch(value) {
        case ObjectViewSettings::Control::SnapshotResolution::scaleWithScreen: return "scale with screen";
        case ObjectViewSettings::Control::SnapshotResolution::manual:          return "manual";
        default:                                                               return "";
    }
}

//-------------------------------------
// Playback settings
//-------------------------------------

struct TrajectoryPlaybackSettings {
    float fps = 6;

};

//-------------------------------------
// All display settings
//-------------------------------------

struct DisplaySettings {
    DisplayMode displayMode = DisplayMode::realtime;

    // GUI
    GuiSettings gui;

    // Main display
    ObjectViewSettings mainView;

    // Trajectory playback
    TrajectoryPlaybackSettings playback;
};

} // namespace medyan::visual

#endif

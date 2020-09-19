#ifndef MEDYAN_Visual_DisplaySettings_hpp
#define MEDYAN_Visual_DisplaySettings_hpp

#include <array>
#include <type_traits>

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

constexpr auto underlying(DisplayMode mode) {
    return static_cast<std::underlying_type_t<DisplayMode>>(mode);
}
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
    bool enabled = true;
};

//-------------------------------------
// Graphics settings
//-------------------------------------

struct ObjectViewSettings {

    struct Canvas {
        std::array< float, 4 > bgColor { 0.0f, 0.0f, 0.0f, 0.0f };
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
        enum class Type { orthographic, perspective };

        Type  type       = Type::perspective;
        float fov        = glm::radians(45.0f); // perspective
        float zNear      = 10.0f;
        float zFar       = 5000.0f;

        // Helper functions

        auto projOrtho(float width, float height) const {
            return glm::ortho(-width / 2, width / 2, -height / 2, height / 2, zNear, zFar);
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
        // camera
        float cameraKeyPositionPerFrame = 150.0f;
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

//-------------------------------------
// All display settings
//-------------------------------------

struct DisplaySettings {
    DisplayMode displayMode = DisplayMode::realtime;

    // GUI
    GuiSettings gui;

    // Main display
    ObjectViewSettings mainView;
};

} // namespace medyan::visual

#endif

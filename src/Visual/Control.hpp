#ifndef MEDYAN_Visual_Control_hpp
#define MEDYAN_Visual_Control_hpp

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"

namespace medyan::visual {

// mouse camera drag

inline void mouseDragRotateCameraAroundTarget(
    ObjectViewSettings::Camera& camera,
    double mouseLastX,
    double mouseLastY,
    double mouseX,
    double mouseY,
    float  cameraRotatePositionPerCursorPixel
) {

    const auto dist = glm::distance(camera.target, camera.position);

    camera.position -= (
        camera.right * float(mouseX - mouseLastX) +
        camera.up    * float(mouseLastY - mouseY)
    ) * cameraRotatePositionPerCursorPixel;
    camera.position = camera.target + glm::normalize(camera.position - camera.target) * dist;
    
    // Update direction
    camera.right = glm::normalize(glm::cross(camera.target - camera.position, camera.up));
    camera.up = glm::normalize(glm::cross(camera.right, camera.target - camera.position));
}

inline void mouseDragPanCamera(
    ObjectViewSettings::Camera& camera,
    double mouseLastX,
    double mouseLastY,
    double mouseX,
    double mouseY,
    float  cameraPanPositionPerCursorPixel
) {
    const auto objectDelta = (
        camera.right * float(mouseX - mouseLastX) +
        camera.up    * float(mouseLastY - mouseY)
    ) * cameraPanPositionPerCursorPixel;

    camera.position -= objectDelta;
    camera.target   -= objectDelta;
}

inline void mouseDragCamera(
    ObjectViewSettings::Camera& camera,
    double mouseLastX,
    double mouseLastY,
    double mouseX,
    double mouseY,
    const ObjectViewSettings::Control& control
) {
    using MM = ObjectViewSettings::Control::CameraMouseMode;

    switch(control.cameraMouseMode) {
        case MM::rotate:
            mouseDragRotateCameraAroundTarget(
                camera,
                mouseLastX, mouseLastY, mouseX, mouseY,
                control.cameraRotatePositionPerCursorPixel
            );
            break;

        case MM::pan:
            mouseDragPanCamera(
                camera,
                mouseLastX, mouseLastY, mouseX, mouseY,
                control.cameraPanPositionPerCursorPixel
            );
            break;

        default:
            break;
    }
}

} // namespace medyan::visual

#endif

#ifndef MEDYAN_Visual_Camera_Hpp
#define MEDYAN_Visual_Camera_Hpp

#include "Util/Math/Vec.hpp"
#include "Visual/Common.hpp"

namespace visual {

// The Camera struct stores the states and parameters of a camera, used in view
// transformation.
// The Camera struct also manages parameters for window control.
struct Camera {
    glm::mat4 view() const {
        return glm::lookAt(position, target, up);
    }

    // Positions
    glm::vec3 position = glm::vec3(0.0f, 0.0f, 2000.0f);
    glm::vec3 target   = glm::vec3(0.0f, 0.0f, 0.0f);
    glm::vec3 right    = glm::vec3(1.0f, 0.0f, 0.0f);
    glm::vec3 up       = glm::vec3(0.0f, 1.0f, 0.0f);

    // Control settings
    float keyControlSpeed     = 150.0f;
    double mouseControlSpeed  = 2.0;

};

} // namespace visual

#endif

#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include "util/environment.h"
#include "Visual/Common.hpp"

#ifdef VISUAL

namespace visual {

class Window {
    void run() const {
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    }
};

} // namespace visual

#endif // ifdef VISUAL

#endif

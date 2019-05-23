#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include <iostream> // cout, endl

#include "util/environment.h"
#include "Visual/Common.hpp"

#ifdef VISUAL

namespace visual {

class Window {
public:
    void run() const {
        // Initialization
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef PLATFORM_MACOS
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

        // Create a new window
        unsigned int windowWidth = 1200;
        unsigned int windowHeight = 800;
        GLFWwindow* window = glfwCreateWindow(windowWidth, windowHeight, "MEDYAN", NULL, NULL);

        if(window == NULL) {
            std::cout << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            return;
        }
        glfwMakeContextCurrent(window);
    }
};

} // namespace visual

#endif // ifdef VISUAL

#endif

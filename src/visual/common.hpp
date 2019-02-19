#ifndef MEDYAN_VISUAL_COMMON_HPP
#define MEDYAN_VISUAL_COMMON_HPP

#include "util/environment.h"

#ifdef VISUAL
    #ifdef PLATFORM_WINDOWS
        #define APIENTRY __stdcall
    #endif
    #include <glad/glad.h>
    #include <GLFW/glfw3.h>
#endif

namespace visual {

#ifdef VISUAL
constexpr bool enabled = true;
#else
constexpr bool enabled = false;
#endif

} // namespace visual

#endif

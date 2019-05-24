#ifndef MEDYAN_Visual_Common_Hpp
#define MEDYAN_Visual_Common_Hpp

#include "Util/Environment.h"

#ifdef VISUAL
    #ifdef PLATFORM_WINDOWS
        #define APIENTRY __stdcall
    #endif
    #include <glad/glad.h>
    #include <GLFW/glfw3.h>

    #include <glm/glm.hpp>
    #include <glm/gtc/matrix_transform.hpp>
    #include <glm/gtc/type_ptr.hpp>
#endif

namespace visual {

#ifdef VISUAL
constexpr bool enabled = true;
#else
constexpr bool enabled = false;
#endif

} // namespace visual

#endif

#ifndef MEDYAN_Visual_Common_Hpp
#define MEDYAN_Visual_Common_Hpp

#ifdef VISUAL
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

#ifndef MEDYAN_VISUAL_COMMON_HPP
#define MEDYAN_VISUAL_COMMON_HPP

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

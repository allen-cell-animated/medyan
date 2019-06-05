#ifndef MEDYAN_Visual_VisualElement_Hpp
#define MEDYAN_Visual_VisualElement_Hpp

#include <mutex>
#include <vector>

#include "Visual/Common.hpp"

namespace visual {

class Profile {
public:
    enum class Target { Filament, Membrane, Force, Compartment };

    // User settings
    bool enabled = false;

    Target target;

    GLenum polygonMode = GL_LINE;
};

struct GlState {
    // vao, vbo, ebo
    GLuint vao, vbo, ebo;

    // Element info
    std::vector< float >    vertexCoords;
    std::vector< unsigned > vertexIndices;
    bool coordChanged = true;
    bool indexChanged = true;
};

// Each VisualElement consists of the following parts:
//   - A mutex that locks all the members
//   - OpenGL compatible data structure as output
//   - All the related settings
struct VisualElement {
    std::mutex me;

    // Set when init
    //-------------------------------------------------------------------------
    GLenum eleMode = GL_TRIANGLES;

    // settings
    //-------------------------------------------------------------------------
    Profile profile;

    // OpenGL data
    //-------------------------------------------------------------------------
    GlState state;

};

} // namespace visual

#endif

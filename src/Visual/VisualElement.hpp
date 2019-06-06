#ifndef MEDYAN_Visual_VisualElement_Hpp
#define MEDYAN_Visual_VisualElement_Hpp

#include <cstdint>
#include <mutex>
#include <vector>

#include "Visual/Common.hpp"

namespace visual {

struct Profile {

    using FlagType = std::uint_fast8_t;

    static constexpr FlagType targetFilament      = 1 << 0;
    static constexpr FlagType targetMembrane      = 1 << 1;
    static constexpr FlagType displayForce        = 1 << 2;
    static constexpr FlagType targetConcentration = 1 << 3;

    // User settings
    bool enabled = false;

    FlagType flag = 0;

    GLenum polygonMode = GL_LINE;
};

struct GlState {
    // vao, vbo, ebo
    GLuint vao, vbo, ebo;

    // Element info
    std::vector< float >    vertexAttribs;
    std::vector< unsigned > vertexIndices;
    bool attribChanged = true;
    bool indexChanged  = true;

    // Set by init and helper
    GLenum eleMode = GL_TRIANGLES;
};

// Each VisualElement consists of the following parts:
//   - A mutex that locks all the members
//   - All the related settings
//   - OpenGL compatible data structure as output
// Note:
//   - Due to constructor of the visual element, the object must be created
//     within an OpenGL context.
struct VisualElement {
    std::mutex me;

    // settings
    //-------------------------------------------------------------------------
    Profile profile;

    // OpenGL data
    //-------------------------------------------------------------------------
    GlState state;

    // ctor and dtor
    VisualElement() {
        glGenBuffers(1, &state.vbo);
        glGenBuffers(1, &state.ebo);
        glGenVertexArrays(1, &state.vao);

        glBindVertexArray(state.vao);
        glBindBuffer(GL_ARRAY_BUFFER, state.vbo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state.ebo);

        // Vertex attribute
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        // temporarily retarget
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    ~VisualElement() {
        glDeleteVertexArrays(1, &state.vao);
        glDeleteBuffers(1, &state.vbo);
        glDeleteBuffers(1, &state.ebo);
    }

};

} // namespace visual

#endif

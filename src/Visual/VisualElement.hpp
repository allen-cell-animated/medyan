#ifndef MEDYAN_Visual_VisualElement_Hpp
#define MEDYAN_Visual_VisualElement_Hpp

#include <cstdint>
#include <mutex>
#include <vector>

#include "Visual/Common.hpp"

namespace visual {

struct Profile {

    using FlagType = std::uint_fast8_t;

    enum class PathMode { Line, Extrude };

    static constexpr FlagType targetFilament      = 1 << 0;
    static constexpr FlagType targetMembrane      = 1 << 1;
    static constexpr FlagType displayForce        = 1 << 2;
    static constexpr FlagType targetConcentration = 1 << 3;

    // User settings
    //-------------------------------------------------------------------------
    bool enabled = false;

    // flag for choosing data
    FlagType flag = 0;

    // shape settings
    PathMode          pathMode = PathMode::Line; // TODO implement this in visual helper
    float             pathExtrudeRadius = 10.0f;
    std::uint_fast8_t pathExtrudeSides = 12;
    float             forceScale = 5.0f;

    GLenum polygonMode = GL_LINE;

    // color settings
};

struct GlState {
    static constexpr unsigned int vaStride = 6;
    static constexpr unsigned int vaPosStart = 0;
    static constexpr unsigned int vaPosSize = 3;
    static constexpr unsigned int vaNormalStart = 3;
    static constexpr unsigned int vaNormalSize = 3;

    // vao, vbo, ebo
    GLuint vao, vbo, ebo;

    // Element info
    std::vector< float >    vertexAttribs;
    std::vector< unsigned > vertexIndices;
    bool attribChanged = true;
    bool indexChanged  = true;

    // Set by init and helper
    GLenum eleMode = GL_TRIANGLES;

    // ctor and dtor
    GlState() {
        glGenBuffers(1, &vbo);
        glGenBuffers(1, &ebo);
        glGenVertexArrays(1, &vao);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

        // Vertex attribute
        //---------------------------------------------------------------------
        // Position
        glVertexAttribPointer(0, vaPosSize,    GL_FLOAT, GL_FALSE, vaStride * sizeof(float), static_cast<const char*>(0) + vaPosStart   );
        glEnableVertexAttribArray(0);
        // Normal
        glVertexAttribPointer(1, vaNormalSize, GL_FLOAT, GL_FALSE, vaStride * sizeof(float), static_cast<const char*>(0) + vaNormalStart);
        glEnableVertexAttribArray(1);

        // temporarily retarget
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    ~GlState() {
        glDeleteVertexArrays(1, &vao);
        glDeleteBuffers(1, &vbo);
        glDeleteBuffers(1, &ebo);
    }

};

// Each VisualElement consists of the following parts:
//   - A mutex that locks all the members
//   - All the related settings
//   - OpenGL compatible data structure as output
// Note:
//   - Due to constructor of the GlState element, the object must be created
//     within an OpenGL context.
struct VisualElement {
    std::mutex me;

    // settings
    //-------------------------------------------------------------------------
    Profile profile;

    // OpenGL data
    //-------------------------------------------------------------------------
    GlState state;

};

} // namespace visual

#endif

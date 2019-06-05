#ifndef MEDYAN_Visual_Profile_Hpp
#define MEDYAN_Visual_Profile_Hpp

#include <functional>

#include "Visual/Common.hpp"
#include "Visual/Shader.hpp"

#ifdef VISUAL

namespace visual {

struct VisualElementSettings {
    bool enabled;

    // Shared data usage
    //-------------------------------------------------------------------------
    // Vertex attribute
    GLuint        vertexAttribIndex;
    GLint         vertexAttribSize;
    GLenum        vertexAttribType = GL_FLOAT;
    GLsizei       vertexAttribStride;
    const GLvoid* vertexAttribPointer = (void*)0;

    // Draw mode
    GLenum polygonMode = GL_LINE;

    GLenum eleMode;
    GLsizei eleCount;
    GLenum indexType = GL_UNSIGNED_INT;

    Shader shader;

    // Data processing and update

    void init() {
        // Set up vertex
        glGenBuffers(1, &vbo);
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &ebo);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

        // Vertex attribute
        glVertexAttribPointer(vertexAttribIndex, vertexAttribSize, vertexAttribType, GL_FALSE, vertexAttribStride, vertexAttribPointer);
        // ^^^ also register vbo as bound
        glEnableVertexAttribArray(0);

        // temporarily retarget
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    } // init()

};

// Each profile is a preset
// TODO

} // namespace visual

#endif // VISUAL

#endif

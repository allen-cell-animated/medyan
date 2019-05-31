#ifndef MEDYAN_Visual_Profile_Hpp
#define MEDYAN_Visual_Profile_Hpp

#include <functional>

#include "Visual/Common.hpp"
#include "Visual/Shader.hpp"
#include "Visual/SharedData.hpp"

#ifdef VISUAL

namespace visual {

// Each profile defines how to render a set of objects
struct Profile {
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

    // Obtained after init
    // TODO: this part should go to the renderer
    unsigned vao;
    unsigned vbo;
    unsigned ebo;

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

    void render() const {
        glPolygonMode(GL_FRONT_AND_BACK, polygonMode);

        glUseProgram(shader.id);
        glBindVertexArray(vao);
        // Update data
        // TODO: This is the job of the renderer
        {
            std::lock_guard<std::mutex> guard(shared::dataMutex);
            if(shared::coordChanged) {
                glBindBuffer(GL_ARRAY_BUFFER, state::vbo[0]);
                glBufferData(GL_ARRAY_BUFFER, sizeof(float) * shared::vertexCoords.size(), &shared::vertexCoords[0], GL_DYNAMIC_DRAW);
                shared::coordChanged = false;
            }
            if(shared::indexChanged) {
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state::ebo[0]);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * shared::triangleVertexIndices.size(), &shared::triangleVertexIndices[0], GL_DYNAMIC_DRAW);
                elementCount = shared::triangleVertexIndices.size();
                shared::indexChanged = false;
            }
        }

        glDrawElements(eleMode, elementCount, indexType, 0);
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);


    } // render()
};

} // namespace visual

#endif // VISUAL

#endif

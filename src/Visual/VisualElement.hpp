#ifndef MEDYAN_Visual_VisualElement_Hpp
#define MEDYAN_Visual_VisualElement_Hpp

#include <cstdint>
#include <memory> // shared_ptr
#include <mutex>
#include <vector>

#include "Visual/Common.hpp"
#include "Visual/Shader.hpp"

namespace visual {

struct Profile {

    using FlagType = std::uint_fast16_t;

    enum class PathMode { Line, Extrude, Bead };
    enum class GridMode { Boundary, Mesh };
    enum class ColorMode { Fixed };

    static constexpr FlagType targetFilament      = 1 << 0;
    static constexpr FlagType targetMembrane      = 1 << 1;
    static constexpr FlagType targetLinker        = 1 << 2;
    static constexpr FlagType targetMotor         = 1 << 3;
    static constexpr FlagType targetBrancher      = 1 << 4;
    static constexpr FlagType displayForce        = 1 << 5;
    static constexpr FlagType forceUseSearchDir   = 1 << 6;
    static constexpr FlagType targetCompartment   = 1 << 7;
    static constexpr FlagType targetConcentration = 1 << 8;

    // User settings
    //-------------------------------------------------------------------------
    bool enabled = false;

    // flag for choosing data
    FlagType flag = 0;

    // shape settings
    // -------- Pathable objects --------
    PathMode          pathMode = PathMode::Line; // TODO implement this in visual helper
    float             pathExtrudeRadius = 7.5f;
    std::uint_fast8_t pathExtrudeSides = 10;
    // -------- Bead-like objects --------
    float             beadRadius = 12.0f;
    std::uint_fast8_t beadLongitudeSegs = 10;
    std::uint_fast8_t beadLatitudeSegs = 5;
    // -------- Forces --------
    float             forceScale = 5.0f;
    // -------- Grid-like objects --------
    GridMode          gridMode = GridMode::Boundary;

    GLenum polygonMode = GL_LINE;

    // color settings
    glm::vec3 colorAmbient; // Currently not used
    glm::vec3 colorDiffuse;
    glm::vec3 colorSpecular  = { 0.7, 0.7, 0.7 };
    float     colorShininess = 32.0f;
};

struct GlSize {
    unsigned int vaStride      = 9;
    unsigned int vaPosStart    = 0;
    unsigned int vaPosSize     = 3;
    unsigned int vaNormalStart = 3;
    unsigned int vaNormalSize  = 3;
    unsigned int vaColorStart  = 6;
    unsigned int vaColorSize   = 3;
};

// Note:
//   - The object of this type must be created within an OpenGL context.
struct GlState {
    const GlSize size; // This is a redundant GlSize for thread safe access

    // vao, vbo, ebo
    GLuint vao;
    GLuint vbo;
    // GLuint ebo;

    // Element info
    std::vector< float >    vertexAttribs;
    // std::vector< unsigned > vertexIndices;
    bool attribChanged = true;
    // bool indexChanged  = true;

    // Set by init and helper
    GLenum eleMode = GL_TRIANGLES;

    // ctor and dtor
    GlState(GlSize newSize) : size(newSize) {
        glGenBuffers(1, &vbo);
        // glGenBuffers(1, &ebo);
        glGenVertexArrays(1, &vao);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

        // Vertex attribute
        //---------------------------------------------------------------------
        GLuint idx = 0;
        // Position
        if(size.vaPosSize) {
            glVertexAttribPointer(idx, size.vaPosSize,    GL_FLOAT, GL_FALSE, size.vaStride * sizeof(float), static_cast<const char*>(0) + sizeof(float) * size.vaPosStart   );
            glEnableVertexAttribArray(idx);
            ++idx;
        }
        // Normal
        if(size.vaNormalSize) {
            glVertexAttribPointer(idx, size.vaNormalSize, GL_FLOAT, GL_FALSE, size.vaStride * sizeof(float), static_cast<const char*>(0) + sizeof(float) * size.vaNormalStart);
            glEnableVertexAttribArray(idx);
            ++idx;
        }
        // Color
        if(size.vaColorSize) {
            glVertexAttribPointer(idx, size.vaColorSize,  GL_FLOAT, GL_FALSE, size.vaStride * sizeof(float), static_cast<const char*>(0) + sizeof(float) * size.vaColorStart );
            glEnableVertexAttribArray(idx);
            ++idx;
        }

        // temporarily retarget
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    ~GlState() {
        glDeleteVertexArrays(1, &vao);
        glDeleteBuffers(1, &vbo);
        // glDeleteBuffers(1, &ebo);
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

    VisualElement(GlSize size) : state(size) {}

};

// Preset configurations for shaders and vector structures
//
// Thread safety: Not thread safe. Creation or modification of objects of this
//     type should only happen in the opengl main thread.
// Note:
//   - The creation/modification/destruction of objects of this type must be
//     within an OpenGL context
struct VisualPreset {
    const GlSize size;
    const Shader shader;

    // Read or modification of visualElements vector requires lock.
    std::mutex veMutex;
    std::vector< std::shared_ptr< VisualElement > > visualElements;

    VisualPreset(GlSize size, const char* vertexShaderSrc, const char* fragmentShaderSrc)
        : size(size), shader(vertexShaderSrc, fragmentShaderSrc) { }

    ~VisualPreset() {
        {
            // Delete visual elements
            std::lock_guard< std::mutex > guard(veMutex);
            visualElements.clear();
        }
    } // ~VisualPreset()
}; // VisualPreset

} // namespace visual

#endif

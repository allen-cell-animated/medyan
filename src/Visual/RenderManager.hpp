#ifndef MEDYAN_Visual_RenderManager_Hpp
#define MEDYAN_Visual_RenderManager_Hpp

#include <thread>

#include "Visual/Common.hpp"

namespace visual {

#ifdef VISUAL

// The render manager is used to manage initiating/destroying the rendering
// context and threads.
// It should be able to store/update the settings for rendering and
// update it.
// It should be able to access system data (as well as settings) and update
// the rendering result.
class RenderManager {

public:
    RenderManager() {}
    ~RenderManager() {
        // Deallocate resources
        glDeleteVertexArrays(1, &_vao);
        glDeleteBuffers(1, &_vbo);

        glfwTerminate();

    }

private:
    std::thread _windowThread;

    unsigned int _vao;
    unsigned int _vbo;

};

#endif

} // namespace visual

#endif

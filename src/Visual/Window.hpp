#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include <array>
#include <iostream> // cout, endl
#include <stdexcept> // runtime_error
#include <vector>

#include "Util/Environment.h"
#include "Util/Io/Log.hpp"
#include "Visual/Camera.hpp"
#include "Visual/Common.hpp"
#include "Visual/Shader.hpp"
#include "Visual/ShaderSrc.hpp"
#include "Visual/SharedData.hpp"
#include "Visual/VisualElement.hpp"

#ifdef VISUAL

// For best portability, the window signal handling could only be done from the
// main thread (due to MacOS Cocoa framework).

namespace visual {

struct Window {
    GLFWwindow* window;
}; // Currently not used

inline void glfwError(int id, const char* description) {
    LOG(ERROR) << description;
    throw std::runtime_error("Error in GLFW environment");
}


template< typename T >
inline void replaceBuffer(GLenum target, const std::vector<T>& source) {
    GLint prevSize;
    glGetBufferParameteriv(target, GL_BUFFER_SIZE, &prevSize);

    const std::size_t newSize = sizeof(T) * source.size();

    if(newSize > prevSize) {
        glBufferData(target, newSize, source.data(), GL_DYNAMIC_DRAW);
    } else {
        glBufferSubData(target, 0, newSize, source.data());
    }
} 


// The RAII object for managing visualization context and window
//
// Note:
//   - Only one object of this type can be created at a time.
//   - The window signal handling should only be done from the main thread (due
//     to MacOS Cocoa framework).
class VisualContext {
public:
    struct WindowStates {
        struct Transformation {
            enum class ProjectionType { Orthographic, Perspective };

            ProjectionType projType = ProjectionType::Perspective;
            float fov               = glm::radians(45.0f); // perspective
            float nearDistance      = 0.1f;
            float farDistance       = 20000.0f;
            glm::mat4 projection;
            glm::mat4 model;

            Camera camera;
        };

        // Size
        int width = 1200;
        int height = 800;

        // Model, view, perspective...
        Transformation trans;

        // Other states
        float deltaTime = 0.01f;
        float lastTime  = 0.0f;
        bool mouseLeftAlreadyPressed = false;
        double mouseLastX;
        double mouseLastY;
    };

    VisualContext() {
        // GLFW initializing
        LOG(DEBUG) << "Initializing GLFW";
        glfwSetErrorCallback(&glfwError);
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

        // Window initializing
        window_ = glfwCreateWindow(windowStates_.width, windowStates_.height, "MEDYAN", NULL, NULL);
        if(window_ == NULL) {
            LOG(ERROR) << "Failed to create GLFW window";
            glfwTerminate();
            return;
        }
        glfwMakeContextCurrent(window_);

        // GLAD initializing
        LOG(DEBUG) << "initializing GLAD";
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            LOG(ERROR) << "Failed to initialize GLAD";
            return;
        }

        glViewport(0, 0, windowStates_.width, windowStates_.height);
        glfwSetWindowUserPointer(window_, &windowStates_);

        // Set window callbacks
        const auto framebufferSizeCallback = [](GLFWwindow* window, int width, int height) {
            auto ws = static_cast< WindowStates* >(glfwGetWindowUserPointer(window));
            ws->width = width;
            ws->height = height;
            glViewport(0, 0, width, height);
        };
        const auto cursorPositionCallback = [](GLFWwindow* window, double xpos, double ypos) {
            auto ws = static_cast< WindowStates* >(glfwGetWindowUserPointer(window));
            auto& camera = ws->trans.camera;

            const int mouseState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
            if(mouseState == GLFW_PRESS) {
                if(ws->mouseLeftAlreadyPressed) {
                    // Transform
                    const double dist = glm::distance(camera.target, camera.position);

                    camera.position -= (camera.right * float(xpos - ws->mouseLastX) + camera.up * float(ws->mouseLastY - ypos)) * (float)camera.mouseControlSpeed;
                    camera.position = camera.target + glm::normalize(camera.position - camera.target) * (float)dist;
                    
                    // Update direction
                    camera.right = glm::normalize(glm::cross(camera.target - camera.position, camera.up));
                    camera.up = glm::normalize(glm::cross(camera.right, camera.target - camera.position));

                } else {
                    ws->mouseLeftAlreadyPressed = true;
                }
                ws->mouseLastX = xpos;
                ws->mouseLastY = ypos;
            } else {
                ws->mouseLeftAlreadyPressed = false;
            }
        };
        const auto scrollCallback = [](GLFWwindow* window, double xoffset, double yoffset) {
            auto ws = static_cast< WindowStates* >(glfwGetWindowUserPointer(window));
            auto& fov = ws->trans.fov;

            fov -= 0.02 * yoffset;
            if(fov < 0.01f) fov = 0.01f;
            if(fov > 3.13f) fov = 3.13f;
        };

        glfwSetFramebufferSizeCallback(window_, framebufferSizeCallback);
        glfwSetCursorPosCallback(window_, cursorPositionCallback);
        glfwSetScrollCallback(window_, scrollCallback);

    } // VisualContext()

    ~VisualContext() {
        glfwTerminate();
    }

    auto window() const { return window_; }
    auto&       windowStates()       { return windowStates_; }
    const auto& windowStates() const { return windowStates_; }

    // Helper function to process window inputs
    void processInput() {

        auto& camera = windowStates_.trans.camera;

        const float currentTime = glfwGetTime();
        windowStates_.deltaTime = currentTime - windowStates_.lastTime;
        windowStates_.lastTime = currentTime;
        const float cameraMove = camera.keyControlSpeed * windowStates_.deltaTime;

        if (glfwGetKey(window_, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            LOG(INFO) << "Escape key hit!";
            glfwSetWindowShouldClose(window_, true);
        }

        if(glfwGetKey(window_, GLFW_KEY_W) == GLFW_PRESS) {
            const auto change = cameraMove * glm::normalize(camera.target - camera.position);
            camera.position += change;
            camera.target += change;
        }
        if(glfwGetKey(window_, GLFW_KEY_S) == GLFW_PRESS) {
            const auto change = cameraMove * glm::normalize(camera.target - camera.position);
            camera.position -= change;
            camera.target -= change;
        }
        if(glfwGetKey(window_, GLFW_KEY_A) == GLFW_PRESS) {
            const auto change = glm::normalize(glm::cross(glm::normalize(camera.target - camera.position), camera.up)) * cameraMove;
            camera.position -= change;
            camera.target -= change;
        }
        if(glfwGetKey(window_, GLFW_KEY_D) == GLFW_PRESS) {
            const auto change = glm::normalize(glm::cross(glm::normalize(camera.target - camera.position), camera.up)) * cameraMove;
            camera.position += change;
            camera.target += change;
        }
    }

private:
    GLFWwindow* window_;
    WindowStates windowStates_;
}; // VisualContext


// The RAII object for all the rendering process
struct VisualDisplay {
    // The overall opengl context. Must be at top
    VisualContext vc;

    // Visual presets
    std::array< VisualPreset, 2 > vps {{
        { GlSize {9, 0, 3, 3, 3, 6, 3}, shader::VertexElementLight, shader::FragElementLight },
        { GlSize {6, 0, 3, 3, 0, 3, 3}, shader::VertexElementLine,  shader::FragElementLine  }
    }};
    VisualPreset& vpLight = vps[0];
    VisualPreset& vpLine  = vps[1];

    VisualDisplay() {
        // Configure global opengl state
        glEnable(GL_DEPTH_TEST);

        // Setup visual with light
        {
            auto& vp = vpLight;
            std::lock_guard< std::mutex > guard(vp.veMutex);

            const auto newVe = [&]() -> auto& {
                vp.visualElements.emplace_back(new VisualElement(vp.size));
                return vp.visualElements.back();
            };

            {
                auto& ve = newVe();
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetMembrane;
                ve->profile.colorAmbient = glm::vec3(0.4f, 0.6f, 0.95f);
                ve->profile.colorDiffuse = glm::vec3(0.4f, 0.6f, 0.95f);
            }
            {
                auto& ve = newVe();
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetMembrane | Profile::displayForce;
                ve->profile.colorAmbient = glm::vec3(0.3f, 0.6f, 0.95f);
                ve->profile.colorDiffuse = glm::vec3(0.3f, 0.6f, 0.95f);
            }
            {
                auto& ve = newVe();
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetFilament;
                ve->profile.pathMode = Profile::PathMode::Extrude;
                ve->profile.colorAmbient = glm::vec3(0.95f, 0.1f, 0.15f);
                ve->profile.colorDiffuse = glm::vec3(0.95f, 0.1f, 0.15f);
            }
            {
                auto& ve = newVe();
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetFilament | Profile::displayForce;
                ve->profile.colorAmbient = glm::vec3(0.95f, 0.1f, 0.15f);
                ve->profile.colorDiffuse = glm::vec3(0.95f, 0.1f, 0.15f);
            }
            {
                auto& ve = newVe();
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetLinker;
                ve->profile.colorAmbient = glm::vec3(0.1f, 0.9f, 0.0f);
                ve->profile.colorDiffuse = glm::vec3(0.1f, 0.9f, 0.0f);
            }
            {
                auto& ve = newVe();
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetMotor;
                ve->profile.colorAmbient = glm::vec3(0.1f, 0.1f, 0.99f);
                ve->profile.colorDiffuse = glm::vec3(0.1f, 0.1f, 0.99f);
            }
            {
                auto& ve = newVe();
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetBrancher;
                ve->profile.colorAmbient = glm::vec3(0.95f, 0.9f, 0.05f);
                ve->profile.colorDiffuse = glm::vec3(0.95f, 0.9f, 0.05f);
            }
        } // ~lock_guard (vpLight)

        // Setup visual for non-lights
        {
            auto& vp = vpLine;
            std::lock_guard< std::mutex > guard(vp.veMutex);

            const auto newVe = [&]() -> auto& {
                vp.visualElements.emplace_back(new VisualElement(vp.size));
                return vp.visualElements.back();
            };

            {
                auto& ve = newVe();
                ve->state.eleMode = GL_LINES;
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetCompartment;
                ve->profile.colorAmbient = glm::vec3(0.8f, 0.8f, 0.8f);
                ve->profile.colorDiffuse = glm::vec3(0.8f, 0.8f, 0.8f);
            }
        } // ~lock_guard (vpLine)

    } // VisualDisplay()

    void run() {

        while (!glfwWindowShouldClose(vc.window())) {
            // input
            vc.processInput();

            auto& ws = vc.windowStates();

            // rendering
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // transform
            ws.trans.projection = glm::perspective(ws.trans.fov, (float)ws.width / (float)ws.height, ws.trans.nearDistance, ws.trans.farDistance);
            ws.trans.model      = glm::mat4(1.0f);
            glm::mat3 modelInvTrans3(glm::transpose(glm::inverse(ws.trans.model)));

            for(auto& vp : vps) {
                std::lock_guard< std::mutex > guard(vp.veMutex);

                if(!vp.visualElements.empty()) {
                    glUseProgram(vp.shader.id());

                    vp.shader.setMat4("projection", ws.trans.projection);
                    vp.shader.setMat4("model",      ws.trans.model);
                    vp.shader.setMat3("modelInvTrans3", modelInvTrans3);
                    vp.shader.setMat4("view",       ws.trans.camera.view());
                    vp.shader.setVec3("CameraPos",  ws.trans.camera.position);

                    for(const auto& ve : vp.visualElements) {
                        std::lock_guard< std::mutex > guard(ve->me);

                        glBindVertexArray(ve->state.vao);

                        glPolygonMode(GL_FRONT_AND_BACK, ve->profile.polygonMode);

                        // Update data
                        if(ve->state.attribChanged) {
                            glBindBuffer(GL_ARRAY_BUFFER, ve->state.vbo);
                            replaceBuffer(GL_ARRAY_BUFFER, ve->state.vertexAttribs);
                            ve->state.attribChanged = false;
                        }

                        // Draw
                        glDrawArrays(ve->state.eleMode, 0, ve->state.vertexAttribs.size() / ve->state.size.vaStride);
                        // glDrawElements(ve->state.eleMode, ve->state.vertexIndices.size(), GL_UNSIGNED_INT, (void*)0);
                    }
                }
            } // End loop visual presets
            glBindVertexArray(0);

            // check
            glfwSwapBuffers(vc.window());
            glfwPollEvents();

        } // End main loop

    } // void run() const
}; // VisualDisplay

} // namespace visual

#endif // ifdef VISUAL

#endif

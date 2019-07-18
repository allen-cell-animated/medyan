#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include <array>
#include <iostream> // cout, endl
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

namespace state {
// Defines variables used in the main thread

enum class ProjectionType { Orthographic, Perspective };

// settings
int windowWidth = 1200;
int windowHeight = 800;
ProjectionType projType = ProjectionType::Perspective;
float fov = glm::radians(45.0f); // perspective
float nearDistance = 0.1f;
float farDistance = 20000.0f;
glm::mat4 projection;
glm::mat4 model;

Camera camera;

float deltaTime = 0.01f;
float lastTime = 0.0f;
bool mouseLeftAlreadyPressed = false;
double mouseLastX;
double mouseLastY;

// gl data
Shader sd(shader::VertexElementLight, shader::FragElementLight); // FIXME: should not be init'ed here

} // namespace state

inline void glfwError(int id, const char* description) {
    LOG(ERROR) << description;
    system("pause");
}

inline void processInput(GLFWwindow* window) {
    using namespace state;

    float currentTime = glfwGetTime();
    deltaTime = currentTime - lastTime;
    lastTime = currentTime;
    float cameraMove = camera.keyControlSpeed * deltaTime;

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        LOG(INFO) << "Escape key hit!";
        glfwSetWindowShouldClose(window, true);
    }

    if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        const auto change = cameraMove * glm::normalize(camera.target - camera.position);
        camera.position += change;
        camera.target += change;
    }
    if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        const auto change = cameraMove * glm::normalize(camera.target - camera.position);
        camera.position -= change;
        camera.target -= change;
    }
    if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        const auto change = glm::normalize(glm::cross(glm::normalize(camera.target - camera.position), camera.up)) * cameraMove;
        camera.position -= change;
        camera.target -= change;
    }
    if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        const auto change = glm::normalize(glm::cross(glm::normalize(camera.target - camera.position), camera.up)) * cameraMove;
        camera.position += change;
        camera.target += change;
    }
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

// The main loop for all windows.
// Note:
//   - This function must be called from the main thread.
inline void mainLoop(GLFWwindow* window) {
    using namespace state;

    // Loop
    LOG(INFO) << "Entering main loop";

    while (!glfwWindowShouldClose(window)) {
        // input
        processInput(window);

        // rendering
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // transform
        state::projection = glm::perspective(state::fov, (float)state::windowWidth / (float)state::windowHeight, state::nearDistance, state::farDistance);
        state::model = glm::mat4(1.0f);
        glm::mat3 modelInvTrans3(glm::transpose(glm::inverse(state::model)));

        state::sd.setMat4("projection", state::projection);
        state::sd.setMat4("model", state::model);
        state::sd.setMat3("modelInvTrans3", modelInvTrans3);
        state::sd.setMat4("view", camera.view());
        state::sd.setVec3("CameraPos", camera.position);

        glUseProgram(state::sd.id());

        {
            std::lock_guard< std::mutex > guard(shared::veMutex);

            for(const auto& ve : shared::visualElements) {
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
        glBindVertexArray(0);

        // check
        glfwSwapBuffers(window);
        glfwPollEvents();
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
    struct WindowSettings {
        int width = 1200;
        int height = 800;
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

        // GLAD initializing
        LOG(DEBUG) << "initializing GLAD";
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            LOG(ERROR) << "Failed to initialize GLAD";
            return;
        }

        // Window initializing
        window_ = glfwCreateWindow(windowSettings_.width, windowSettings_.height, "MEDYAN", NULL, NULL);
        if(window_ == NULL) {
            LOG(ERROR) << "Failed to create GLFW window";
            glfwTerminate();
            return;
        }
        glfwMakeContextCurrent(window_);

        glViewport(0, 0, windowSettings_.width, windowSettings_.height);
        glfwSetWindowUserPointer(window_, &windowSettings_);

        // Set window callbacks
        const auto framebufferSizeCallback = [](GLFWwindow* window, int width, int height) {
            auto ws = static_cast< WindowSettings* >(glfwGetWindowUserPointer(window));
            ws->width = width;
            ws->height = height;
            glViewport(0, 0, width, height);
        };
        const auto cursorPositionCallback = [](GLFWwindow* window, double xpos, double ypos) {
            using namespace state;

            int mouseState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
            if(mouseState == GLFW_PRESS) {
                if(mouseLeftAlreadyPressed) {
                    // Transform
                    double dist = glm::distance(camera.target, camera.position);

                    camera.position -= (camera.right * float(xpos - mouseLastX) + camera.up * float(mouseLastY - ypos)) * (float)camera.mouseControlSpeed;
                    camera.position = camera.target + glm::normalize(camera.position - camera.target) * (float)dist;
                    
                    // Update direction
                    camera.right = glm::normalize(glm::cross(camera.target - camera.position, camera.up));
                    camera.up = glm::normalize(glm::cross(camera.right, camera.target - camera.position));

                } else {
                    mouseLeftAlreadyPressed = true;
                }
                mouseLastX = xpos;
                mouseLastY = ypos;
            } else {
                mouseLeftAlreadyPressed = false;
            }
        };
        const auto scrollCallback = [](GLFWwindow* window, double xoffset, double yoffset) {
            using namespace state;
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

private:
    GLFWwindow* window_;
    WindowSettings windowSettings_;
};

struct VisualDisplay {
    // The overall opengl context. Must be at top
    VisualContext vc;

    // Visual presets
    std::array< VisualPreset, 1 > vps {{
        { GlSize {}, shader::VertexElementLight, shader::FragElementLight }
    }};
    VisualPreset& vpLight = vps[0];

    VisualDisplay() {
        // Configure global opengl state
        glEnable(GL_DEPTH_TEST);

        // Setup visual with light
        {
            std::lock_guard< std::mutex > guard(vpLight.veMutex);

            const auto newVe = [this]() -> auto& {
                vpLight.visualElements.emplace_back(new VisualElement(vpLight.size));
                return vpLight.visualElements.back();
            };

            {
                auto& ve = newVe();
                ve->profile.enabled = true;
                ve->profile.flag = Profile::targetCompartment;
                ve->profile.colorAmbient = glm::vec3(0.8f, 0.8f, 0.8f);
                ve->profile.colorDiffuse = glm::vec3(0.8f, 0.8f, 0.8f);
            }
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
        } // ~lock_guard

    } // VisualDisplay()

    void run() const { mainLoop(vc.window()); }
};

} // namespace visual

#endif // ifdef VISUAL

#endif

#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include <iostream> // cout, endl
#include <vector>

#include "Util/Environment.h"
#include "util/io/log.h"
#include "Visual/Camera.hpp"
#include "Visual/Common.hpp"
#include "Visual/Shader.hpp"
#include "Visual/SharedData.hpp"

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
GLFWwindow* window;
unsigned int vao[2]; // 0 for coord, 1 for force
unsigned int vbo[2];
unsigned int ebo[2];
Shader sd;
unsigned int elementCount;
unsigned int forceElementCount;

} // namespace state

inline void glfwError(int id, const char* description) {
    LOG(ERROR) << description;
    system("pause");
}

inline void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    state::windowWidth = width;
    state::windowHeight = height;
    glViewport(0, 0, width, height);
}
inline void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    using namespace state;
    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    if(state == GLFW_PRESS) {
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
}
inline void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    using namespace state;
    fov -= 0.02 * yoffset;
    if(fov < 0.01f) fov = 0.01f;
    if(fov > 3.13f) fov = 3.13f;
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

inline void createWindow() {
    LOG(INFO) << "Initializing GLFW";
    glfwSetErrorCallback(&glfwError);
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    state::window = glfwCreateWindow(state::windowWidth, state::windowHeight, "MEDYAN", NULL, NULL);
    if(state::window == NULL) {
        LOG(ERROR) << "Failed to create GLFW window";
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(state::window);

    LOG(INFO) << "initializing GLAD";
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        LOG(ERROR) << "Failed to initialize GLAD";
        return;
    }
    glViewport(0, 0, state::windowWidth, state::windowHeight);
    glfwSetFramebufferSizeCallback(state::window, framebuffer_size_callback);
    glfwSetCursorPosCallback(state::window, cursor_position_callback);
    glfwSetScrollCallback(state::window, scroll_callback);

    // Configure global opengl state
    glEnable(GL_DEPTH_TEST);

    // Shader
    const char* const vertexshader = R"(
#version 330 core
layout (location = 0) in vec3 aPos;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main() {
    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)";

    const char* const fragmentshader = R"(
#version 330 core
out vec4 FragColor;

void main() {
    FragColor = vec4(0.5f, 0.25f, 0.1f, 1.0f);
}
)";
    state::sd.init(vertexshader, fragmentshader);

    // Set up vertex
    glGenBuffers(2, state::vbo);
    glGenVertexArrays(2, state::vao);
    glGenBuffers(2, state::ebo);
    glBindVertexArray(state::vao[0]); // Bind this first!

    // Coord vertex array
    glBindBuffer(GL_ARRAY_BUFFER, state::vbo[0]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state::ebo[0]);
    // Vertex attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    // ^^^ also register vbo as bound
    glEnableVertexAttribArray(0);

    // temporarily retarget
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // Force vertex array
    glBindVertexArray(state::vao[1]);
    glBindBuffer(GL_ARRAY_BUFFER, state::vbo[1]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state::ebo[1]);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // Draw wireframe
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

}

// The main loop for all windows.
// Note:
//   - This function must be called from the main thread.
inline void mainLoop() {
    using namespace state;

    // Loop
    LOG(INFO) << "Entering main loop";

    while (!glfwWindowShouldClose(state::window)) {
        // input
        processInput(state::window);

        // rendering
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // transform
        state::projection = glm::perspective(state::fov, (float)state::windowWidth / (float)state::windowHeight, state::nearDistance, state::farDistance);
        state::sd.setMat4("projection", state::projection);
        state::model = glm::mat4(1.0f);
        //state::model = glm::rotate(state::model, 10.0f * (float)glfwGetTime(), glm::vec3(0.6f, 0.8f, 0.0f));
        state::sd.setMat4("model", state::model);

        state::sd.setMat4("view", camera.view());

        glUseProgram(state::sd.id);
        glBindVertexArray(state::vao[0]);
        // Update data
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

        glDrawElements(GL_TRIANGLES, elementCount, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        glBindVertexArray(state::vao[1]);
        {
            std::lock_guard<std::mutex> guard(shared::dataMutex);
            if(shared::forceChanged) {
                glBindBuffer(GL_ARRAY_BUFFER, state::vbo[1]);
                glBufferData(GL_ARRAY_BUFFER, sizeof(float) * shared::arrowVertexCoords.size(), &shared::arrowVertexCoords[0], GL_DYNAMIC_DRAW);
                shared::forceChanged = false;
            }
            if(shared::forceIndexChanged) {
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state::ebo[1]);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * shared::lineVertexIndices.size(), &shared::lineVertexIndices[0], GL_DYNAMIC_DRAW);
                forceElementCount = shared::lineVertexIndices.size();
                shared::forceIndexChanged = false;
            }
        }
        glDrawElements(GL_LINES, forceElementCount, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0); // no need to unbind every time

        // check
        glfwSwapBuffers(state::window);
        glfwPollEvents();
    }
}

inline void deallocate() {

    // Deallocate resources
    glDeleteVertexArrays(2, state::vao);
    glDeleteBuffers(2, state::vbo);
    glDeleteBuffers(2, state::ebo);

    glfwTerminate();

}

} // namespace visual

#endif // ifdef VISUAL

#endif

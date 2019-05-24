#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include <iostream> // cout, endl
#include <vector>

#include "Util/Environment.h"
#include "util/io/log.h"
#include "Visual/Common.hpp"
#include "Visual/Shader.hpp"

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
float farDistance = 100.0f;
glm::mat4 projection;
glm::mat4 view;
glm::mat4 model;

glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 50.0f);
glm::vec3 cameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 cameraRight = glm::vec3(1.0f, 0.0f, 0.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
float cameraSpeed = 6.0f;
float deltaTime = 0.01f;
float lastTime = 0.0f;
double mouseSpeed = 0.5;
bool mouseLeftAlreadyPressed = false;
double mouseLastX;
double mouseLastY;

GLFWwindow* window;
unsigned int vao;
unsigned int vbo;
unsigned int ebo;
Shader sd;

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
            double dist = glm::distance(cameraTarget, cameraPos);

            cameraPos -= cameraRight * float(xpos - mouseLastX) + cameraUp * float(mouseLastY - ypos);
            cameraPos = cameraTarget + glm::normalize(cameraPos - cameraTarget) * (float)dist;
            
            // Update direction
            cameraRight = glm::normalize(glm::cross(cameraTarget - cameraPos, cameraUp));
            cameraUp = glm::normalize(glm::cross(cameraRight, cameraTarget - cameraPos));

        } else {
            mouseLeftAlreadyPressed = true;
        }
        mouseLastX = xpos;
        mouseLastY = ypos;
    } else {
        mouseLeftAlreadyPressed = false;
    }
}
inline void processInput(GLFWwindow* window) {
    using namespace state;

    float currentTime = glfwGetTime();
    deltaTime = currentTime - lastTime;
    lastTime = currentTime;
    float cameraMove = cameraSpeed * deltaTime;

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        LOG(INFO) << "Escape key hit!";
        glfwSetWindowShouldClose(window, true);
    }

    if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        const auto change = cameraMove * glm::normalize(cameraTarget - cameraPos);
        cameraPos += change;
        cameraTarget += change;
    }
    if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        const auto change = cameraMove * glm::normalize(cameraTarget - cameraPos);
        cameraPos -= change;
        cameraTarget -= change;
    }
    if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        const auto change = glm::normalize(glm::cross(glm::normalize(cameraTarget - cameraPos), cameraUp)) * cameraMove;
        cameraPos -= change;
        cameraTarget -= change;
    }
    if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        const auto change = glm::normalize(glm::cross(glm::normalize(cameraTarget - cameraPos), cameraUp)) * cameraMove;
        cameraPos += change;
        cameraTarget += change;
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
    // glfwSetInputMode(state::window, GLFW_STICKY_MOUSE_BUTTONS, GLFW_TRUE);

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
    glGenBuffers(1, &state::vbo);
    glGenVertexArrays(1, &state::vao);
    glGenBuffers(1, &state::ebo);
    glBindVertexArray(state::vao); // Bind this first!

    glBindBuffer(GL_ARRAY_BUFFER, state::vbo);
    float vertices[]{
        -15.5f, -25.9f, 8.0f,
        30.0f,-0.5f,0.0f,
        0.0f, 30.0f,-7.0f,
        -9.9f,0.4f,0.1f
    };
    unsigned int indices[]{
        0, 1, 2,
        1, 2, 3
    };
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state::ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    // Vertex attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    // ^^^ also register vbo as bound
    glEnableVertexAttribArray(0);

    // temporarily retarget
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Draw wireframe
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

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

        state::view = glm::lookAt(cameraPos, cameraTarget, cameraUp);
        state::sd.setMat4("view", state::view);

        glUseProgram(state::sd.id);
        glBindVertexArray(state::vao);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        // glBindVertexArray(0); // no need to unbind every time

        // check
        glfwSwapBuffers(state::window);
        glfwPollEvents();
    }
}

inline void deallocate() {

    // Deallocate resources
    glDeleteVertexArrays(1, &state::vao);
    glDeleteBuffers(1, &state::vbo);

    glfwTerminate();

}

} // namespace visual

#endif // ifdef VISUAL

#endif

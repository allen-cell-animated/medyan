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
    //LOG(INFO) << "Changing window size: " << width << " x " << height;
    glViewport(0, 0, width, height);
}
inline void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        LOG(INFO) << "Escape key hit!";
        glfwSetWindowShouldClose(window, true);
    }
}

inline void createWindow(unsigned int width, unsigned int height) {
    LOG(INFO) << "Initializing GLFW";
    glfwSetErrorCallback(&glfwError);
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    state::window = glfwCreateWindow(width, height, "MEDYAN", NULL, NULL);
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
    glViewport(0, 0, width, height);
    glfwSetFramebufferSizeCallback(state::window, framebuffer_size_callback);

    // Shader
    const char* const vertexshader = R"(
#version 330 core
layout (location = 0) in vec3 aPos;

void main() {
    gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);
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
        -0.5f, -0.9f, 0.0f,
        1.2f,-0.5f,0.0f,
        0.0f, 0.5f,0.0f,
        -0.9f,0.4f,0.1f
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
    // Loop
    LOG(INFO) << "Entering main loop";

    while (!glfwWindowShouldClose(state::window)) {
        // input
        processInput(state::window);

        // rendering
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

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

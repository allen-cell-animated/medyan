#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include <iostream> // cout, endl

#include "Util/Environment.h"
#include "util/io/log.h"
#include "Visual/Common.hpp"

#ifdef VISUAL

namespace visual {

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

inline void init() {
    LOG(INFO) << "Initializing GLFW";
    glfwSetErrorCallback(&glfwError);
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    unsigned int windowWidth = 1200;
    unsigned int windowHeight = 800;

    GLFWwindow* window = glfwCreateWindow(windowWidth, windowHeight, "MEDYAN", NULL, NULL);
    if(window == NULL) {
        LOG(ERROR) << "Failed to create GLFW window";
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(window);

    LOG(INFO) << "initializing GLAD";
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        LOG(ERROR) << "Failed to initialize GLAD";
        return;
    }
    glViewport(0, 0, windowWidth, windowHeight);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Shader
    const char* const vertexshader = R"(
#version 330 core
layout (location = 0) in vec3 aPos;

void main() {
    gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);
}
)";
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexshader, NULL);
    glCompileShader(vertexShader);
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        LOG(ERROR) << "Shader compile failed: " << infoLog;
    }

    const char* const fragmentshader = R"(
#version 330 core
out vec4 FragColor;

void main() {
    FragColor = vec4(0.5f, 0.25f, 0.1f, 1.0f);
}
)";
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentshader, NULL);
    glCompileShader(fragmentShader);

    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        LOG(ERROR) << "Shader program link failed: " << infoLog;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // Set up vertex
    unsigned int vbo;
    glGenBuffers(1, &vbo);
    unsigned int vao, ebo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &ebo);
    glBindVertexArray(vao); // Bind this first!

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
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
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
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


    // Loop
    LOG(INFO) << "Entering main loop";
    while (!glfwWindowShouldClose(window)) {
        // input
        processInput(window);

        // rendering
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        glBindVertexArray(vao);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        // glBindVertexArray(0); // no need to unbind every time

        // check
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Deallocate resources
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);

    glfwTerminate();

}


} // namespace visual

#endif // ifdef VISUAL

#endif

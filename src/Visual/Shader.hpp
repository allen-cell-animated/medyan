#ifndef MEDYAN_Visual_Shader_Hpp
#define MEDYAN_Visual_Shader_Hpp

#ifdef VISUAL

#include "util/io/log.h"
#include "Visual/Common.hpp"

namespace visual {

struct Shader {
    unsigned int id;

    void init(const char* vertexShaderSrc, const char* fragmentShaderSrc) {
        int success;
        char infoLog[512];

        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSrc, NULL);
        glCompileShader(vertexShader);
        glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
            LOG(ERROR) << "Vertex shader compile failed: " << infoLog;
        }

        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSrc, NULL);
        glCompileShader(fragmentShader);
        glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
        if (!success) {
            glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
            LOG(ERROR) << "Fragment shader compile failed: " << infoLog;
        }

        id = glCreateProgram();
        glAttachShader(id, vertexShader);
        glAttachShader(id, fragmentShader);
        glLinkProgram(id);
        glGetProgramiv(id, GL_LINK_STATUS, &success);
        if (!success) {
            glGetProgramInfoLog(id, 512, NULL, infoLog);
            LOG(ERROR) << "Shader program link failed: " << infoLog;
        }

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

    }
};

} // namespace visual

#endif // VISUAL

#endif

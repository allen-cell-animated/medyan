#ifndef MEDYAN_Visual_ShaderSrc_Hpp
#define MEDYAN_Visual_ShaderSrc_Hpp

namespace visual {
namespace shader {

constexpr const char* VertexElement = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

out vec3 ModelPos;
out vec3 Normal;

uniform mat4 model;
uniform mat3 modelInvTrans3;
uniform mat4 view;
uniform mat4 projection;

void main() {
    ModelPos = vec3(model * vec4(aPos, 1.0));
    Normal = modelInvTrans3 * aNormal;

    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)";

constexpr const char* FragElement = R"(
#version 330 core
out vec4 FragColor;

in vec3 ModelPos;
in vec3 Normal;

uniform vec3 CameraPos;

void main() {
    FragColor = vec4(0.2f, 0.5f, 0.95f, 1.0f);
}
)";

} // namespace shader
} // namespace visual

#endif

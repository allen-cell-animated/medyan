#ifndef MEDYAN_Visual_ShaderSrc_Hpp
#define MEDYAN_Visual_ShaderSrc_Hpp

namespace visual {
namespace shader {

constexpr const char* VertexElement = R"(
#version 330 core
layout (location = 0) in vec3 aPos;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main() {
    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)";

constexpr const char* FragElement = R"(
#version 330 core
out vec4 FragColor;

void main() {
    FragColor = vec4(0.5f, 0.25f, 0.1f, 1.0f);
}
)";

} // namespace shader
} // namespace visual

#endif

#ifndef MEDYAN_Visual_ShaderSrc_Hpp
#define MEDYAN_Visual_ShaderSrc_Hpp

namespace visual {
namespace shader {

constexpr const char* VertexElementLight = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec3 aColor;

out vec3 ModelPos;
out vec3 Normal;
out vec3 Color;

uniform mat4 model;
uniform mat3 modelInvTrans3;
uniform mat4 view;
uniform mat4 projection;

void main() {
    ModelPos = vec3(model * vec4(aPos, 1.0));
    Normal = modelInvTrans3 * aNormal;

    // Simplified color (used in ambient and diffuse color)
    Color = aColor;

    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)";

constexpr const char* FragElementLight = R"(
#version 330 core
out vec4 FragColor;

in vec3 ModelPos;
in vec3 Normal;
in vec3 Color;

struct Material {
    vec3 diffuse; // also for ambient
    vec3 specular;
    float shininess;
};

struct DirLight {
    vec3 direction;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

struct PointLight {
    vec3 position;

    float constant;
    float linear;
    float quadratic;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

#define NUM_DIR_LIGHTS 1
#define NUM_POINT_LIGHTS 4

uniform vec3 CameraPos;
uniform DirLight dirLights[NUM_DIR_LIGHTS];
uniform PointLight pointLights[NUM_POINT_LIGHTS];
uniform Material material;

vec3 calcDirLight(DirLight light, vec3 normal, vec3 viewDir);
vec3 calcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir);

void main() {
    FragColor = vec4(Color, 1.0f);
}

vec3 calcDirLight(DirLight light, vec3 normal, vec3 viewDir) {
    vec3 lightDir = normalize(-light.direction);
    // Diffuse
    float diffuseFac = max(dot(normal, lightDir), 0.0);
    // Specular
    vec3 reflectDir = reflect(-lightDir, normal);
    float specularFac = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);
    // Combine
    vec3 ambient = light.ambient * material.diffuse;
    vec3 diffuse = light.diffuse * diffuseFac * material.diffuse;
    vec3 specular = light.specular * specularFac * material.specular;
    return ambient + diffuse + specular;
}

vec3 calcPointLight(PointLight light, vec3 normal, vec3 fragPos, vec3 viewDir) {
    vec3 lightDir = normalize(light.position - fragPos);
    // Diffuse
    float diffuseFac = max(dot(normal, lightDir), 0.0);
    // Specular
    vec3 reflectDir = reflect(-lightDir, normal);
    float specularFac = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);
    // Attenuation
    float distance = length(light.position - fragPos);
    float attenuation = 1.0 / (light.constant + light.linear * distance + light.quadratic * distance * distance);
    // Combine
    vec3 ambient = light.ambient * material.diffuse;
    vec3 diffuse = light.diffuse * diffuseFac * material.diffuse;
    vec3 specular = light.specular * specularFac * material.specular;
    return ambient + diffuse + specular;
}
)";

constexpr const char* VertexElementLine = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;

out vec3 Color;

uniform mat4 model;
uniform mat3 modelInvTrans3;
uniform mat4 view;
uniform mat4 projection;

void main() {
    // Simplified color (used in ambient and diffuse color)
    Color = aColor;

    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)";

constexpr const char* FragElementLine = R"(
#version 330 core
out vec4 FragColor;

in vec3 Color;

void main() {
    FragColor = vec4(Color, 1.0f);
}
)";

} // namespace shader
} // namespace visual

#endif

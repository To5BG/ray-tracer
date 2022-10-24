#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    float totalArea = glm::length(glm::cross((v0 - v2), (v1 - v2))) / 2;
    float alpha = (glm::length(glm::cross(p - v2, v1 - v2))) / (2 * totalArea);
    float beta = (glm::length(glm::cross(p - v2, v0 - v2))) / (2 * totalArea);
    float gamma = (glm::length(glm::cross(p - v1, v0 - v1))) / (2 * totalArea);
    return glm::vec3(alpha, beta, gamma);
}

glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    // TODO: implement this function.
    glm::vec3 interpolatedNormal = barycentricCoord.x * n0 + barycentricCoord.y * n1 + barycentricCoord.z * n2;
    return glm::normalize(interpolatedNormal);
}

glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    // TODO: implement this function.
    glm::vec2 interpolatedCoord;
    float interpolatedCoordX = t0.x * barycentricCoord.x + t1.x * barycentricCoord.y + t2.x * barycentricCoord.z;
    float interpolatedCoordY = t0.y * barycentricCoord.x + t1.y * barycentricCoord.y + t2.y * barycentricCoord.z;
    interpolatedCoord = glm::vec2(interpolatedCoordX, interpolatedCoordY);
    return interpolatedCoord;
}

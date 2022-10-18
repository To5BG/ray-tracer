#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // TODO: implement the Phong shading model.
    return hitInfo.material.kd;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    Ray reflectionRay {};
    float epsilon = 0.00001f;
    reflectionRay.direction = ray.direction - (2 * glm::dot(ray.direction, hitInfo.normal)) * hitInfo.normal;
    reflectionRay.origin = ray.t * ray.direction*(1-epsilon) + ray.origin;
    
    return reflectionRay;
}
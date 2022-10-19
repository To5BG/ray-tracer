#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <iostream>
const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // basic vectors
    glm::vec3 VertexPos = (ray.origin + (ray.direction * ray.t));
    glm::vec3 N = normalize(hitInfo.normal);
    glm::vec3 L = glm::normalize(lightPosition - VertexPos);

    // diffuse
    float lambertian = glm::dot(N, L) > 0 ? glm::dot(N, L) : 0;
    glm::vec3 diffuse;

    // if light behind object
    if (lambertian <= 0) 
        diffuse = glm::vec3{ 0 };
    
    diffuse = lightColor * hitInfo.material.kd * lambertian;

    // specular
    float specular = 0.0;

    // if visible
    if (lambertian > 0) {
        glm::vec3 R = glm::reflect(L, N);
        glm::vec3 V = glm::normalize(-ray.origin + VertexPos);
        float specAngle = glm::dot(R, V) > 0 ? glm::dot(R, V) : 0;
        specular = pow(specAngle, hitInfo.material.shininess);
    }

    glm::vec3 specularVector = lightColor * hitInfo.material.ks * specular;

  //  std::cout << diffuse.x << " " << diffuse.y << " " <<diffuse.z << std::endl;

    return specularVector + diffuse;
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
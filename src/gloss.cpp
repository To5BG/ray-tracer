#include "gloss.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <random>
#include "light.h"
#include <framework/opengl_includes.h>

int extr_glossy_filterSize = 64;
float extr_glossy_sigma = 1;
std::normal_distribution<> normal(0.0, 1.0);
std::mt19937 mtGen(std::random_device {}());

std::vector<Ray> getRaySamples(HitInfo& hitInfo, Ray& ref)
{
    // Calculate square basis
    glm::vec3 w = glm::normalize(ref.direction);
    // Follow textbook's suggestion - replace smallest magnitude 
    // dimension with 1 to get a vector sufficiently different from w
    int minTerm = fabs(w.x) <= fabs(w.y) && fabs(w.x) <= fabs(w.z) ? 0 : (fabs(w.y) <= fabs(w.x) && fabs(w.y) <= fabs(w.z) ? 1 : 2);
    glm::vec3 t = w;
    t[minTerm] = 1;
    // Better performance, but edge case - if incident angle is exactly 90, it will result in a zero vector
    // glm::vec3 u = glm::normalize(glm:cross(hitInfo.normal, w));
    glm::vec3 u = glm::normalize(glm::cross(t, w));
    glm::vec3 v = glm::cross(w, u);

    // Visual debug -> square that covers all possible ray perturbations
    glColor3f( 0.5f, 0.0f, 0.5f );
    // Side of sample square -> proportional to sigma, inversely proportional to shininess
    float side = extr_glossy_sigma / (hitInfo.material.shininess * 3.0f);
    drawTriangle(
        Vertex { ref.origin + 0.3f * w - side * (u + v) }, 
        Vertex { ref.origin + 0.3f * w + side * (u - v) },
        Vertex { ref.origin + 0.3f * w + side * (v - u) });
    drawTriangle(
        Vertex { ref.origin + 0.3f * w + side * (u + v) }, 
        Vertex { ref.origin + 0.3f * w + side * (v - u) }, 
        Vertex { ref.origin + 0.3f * w + side * (u - v) });

    std::vector<Ray> rays;
    // 0,0 on distribution space is 0.5, 0.5 on square's space -> offset by center
    float offset = -side / 2.0f;
    for (int i = 0; i < extr_glossy_filterSize; i++) {
        glm::vec3 rPrime = ref.direction + float(offset + normal(mtGen) * side) * u + float(offset + normal(mtGen) * side) * v;
        // Make sure perturbed ray is not below the surface its reflected from (at large sigmas and/or small incident angle)
        if (glm::dot(rPrime, hitInfo.normal) > 0)
            rays.push_back({ ref.origin, glm::normalize(rPrime), std::numeric_limits<float>::max() });
    }
    return rays;
}

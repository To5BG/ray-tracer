#include "gloss.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <random>
#include "light.h"

int extr_glossy_filterSize = 20;
float extr_glossy_sigma = 1;
std::normal_distribution<> normal(0.0, 1.0);
std::mt19937 mtGen(std::random_device {}());

std::vector<Ray> getRaySamples(HitInfo& hitInfo, Ray& ref)
{
    // Calculate square basis
    glm::vec3 w = glm::normalize(ref.direction);
    glm::vec3 u = glm::normalize(glm::cross(w, hitInfo.normal));
    glm::vec3 v = glm::cross(w, u);

    float side = extr_glossy_sigma / (hitInfo.material.shininess * 4.0f);
    drawTriangle(Vertex { ref.origin + 0.25f * w - side * (u + v) }, 
        Vertex { ref.origin + 0.25f * w + side * (u - v) },
        Vertex { ref.origin + 0.25f * w + side * (v - u) });
    drawTriangle(Vertex { ref.origin + 0.2f * w + side * (u + v) }, 
        Vertex { ref.origin + 0.25f * w + side * (v - u) }, 
        Vertex { ref.origin + 0.25f * w + side * (u - v) });

    std::vector<Ray> rays;
    float offset = -side / 2.0f;
    for (int i = 0; i < extr_glossy_filterSize; i++) {
        glm::vec3 rPrime = ref.direction + float(offset + normal(mtGen) * side) * u + float(offset + normal(mtGen) * side) * v;
        rays.push_back({ ref.origin, glm::normalize(rPrime), std::numeric_limits<float>::max() });
    }
    return rays;
}

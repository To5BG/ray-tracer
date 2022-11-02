#include "gloss.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <random>

int extr_glossy_filterSize = 20;
float extr_glossy_sigma = 0.03;

std::vector<Ray> getRaySamples(HitInfo& hitInfo, Ray& ref)
{
    // Calculate square basis
    glm::vec3 w = glm::normalize(ref.direction);
    glm::vec3 u = glm::normalize(glm::cross(w, hitInfo.normal));
    glm::vec3 v = glm::cross(w, u);

    std::vector<Ray> rays;
    std::random_device random;
    std::mt19937 mtGen(random());
    std::normal_distribution<> normal(0.0, 1.0);

    float a = extr_glossy_sigma / 100.0f;
    float alpha = -a / 2.0f + normal(mtGen) * a;
    float beta = -a / 2.0f + normal(mtGen) * a;

    for (int i = 0; i < extr_glossy_filterSize; i++) {
        alpha = -a / 2.0f + normal(mtGen) * a;
        beta = -a / 2.0f + normal(mtGen) * a;
        rays.push_back({ ref.origin, glm::normalize(ref.direction + alpha * u + beta * v), std::numeric_limits<float>::max() });
    }
    return rays;
}

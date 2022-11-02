#include "gloss.h"
#include <iostream>
#include <algorithm>
#include <map>

//std::map<int, ContinuousGaussianFilter> glossy_gaussian_filters;

AxisAlignedBox getBlurSquare(Ray& ray) 
{
    glm::vec3 w = glm::normalize(ray.direction);
    // p.334 of textbook - ensure we choose a t that is different enough from w to form a basis
    float minDimension = std::fmin(std::fmin(w.x, w.y), w.z);
    glm::vec3 t = (minDimension == w.x) ? glm::vec3 { 1.0f, w.y, w.z } : (minDimension == w.y ? glm::vec3 { w.x, 1.0f, w.z } : glm::vec3 { w.x, w.y, 1.0f });
    // Calculate basis
    glm::vec3 u = glm::normalize(glm::cross(w, t));
    glm::vec3 v = glm::normalize(glm::cross(w, u));

    return { ray.origin + w - (u + v), ray.origin + w + (u + v) };
}

// Helper for returning pi
const float pi = std::atan(1) * 4;

ContinuousGaussianFilter::ContinuousGaussianFilter()
{
    ContinuousGaussianFilter(3, 0.8f, 1.0f);
}
    
// Create Gaussian filter
ContinuousGaussianFilter::ContinuousGaussianFilter(int f, float s, float d)
{
    filterSize = f;
    //dist = d;
    sigma = pow(s, 2);

    kernel.reserve(pow(2 * filterSize + 1, 2));
    float divTerm = 1.0f / (2 * pi * sigma);
    float sum = 0;

    for (float i = -filterSize; i <= filterSize; i++) {
        float xTerm = std::exp(-(i * i) / (2.0f * sigma)) * divTerm;
        for (int j = -filterSize; j <= filterSize; j++) {
            float term = xTerm * std::exp(-(j * j) / (2.0f * sigma));
            sum += term;
            kernel.push_back(term);
        }
    }
    // Normalize
    std::for_each(kernel.begin(), kernel.end(), [&](float& e) { e /= sum; std::cout << e << std::endl; });
}

ContinuousGaussianFilter glossy_filter = ContinuousGaussianFilter {};
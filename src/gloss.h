#pragma once
#include "common.h"
#include <framework/ray.h>
#include <glm/glm.hpp>
#include <map>

AxisAlignedBox getBlurSquare(Ray& ray);

struct ContinuousGaussianFilter {
    int filterSize;
    float sigma;
    //float dist;
    std::vector<float> kernel;

    ContinuousGaussianFilter();

    ContinuousGaussianFilter(int f, float s, float d);
};

extern ContinuousGaussianFilter glossy_filter;
// extern std::map <int, ContinuousGaussianFilter> glossy_gaussian_filters;
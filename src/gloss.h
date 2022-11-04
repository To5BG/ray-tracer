#pragma once
#include "common.h"
#include <framework/ray.h>
#include <glm/glm.hpp>
#include <map>

extern int extr_glossy_filterSize;
extern float extr_glossy_sigma;

// Returns an array with normally distributed ray samples from another ray
std::vector<Ray> getRaySamples(HitInfo& hitInfo, Ray& ray);

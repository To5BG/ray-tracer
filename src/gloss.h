#pragma once
#include "common.h"
#include <framework/ray.h>
#include <glm/glm.hpp>
#include <map>

extern int extr_glossy_filterSize;
extern float extr_glossy_sigma;

std::vector<Ray> getRaySamples(HitInfo& hitInfo, Ray& ray);

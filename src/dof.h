#pragma once
#include "common.h"
#include <framework/ray.h>
#include <glm/glm.hpp>
#include <map>
#include <framework/trackball.h>

extern int extr_dof_samples;
extern float extr_dof_aperture;
extern float extr_dof_f;
extern bool draw_random_rays;

extern std::vector<Ray> dofRays;
extern std::vector<Ray> camFrame;

std::vector<Ray> getEyeFrame(Ray& ray, const Trackball& camera);


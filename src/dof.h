#pragma once
#include "common.h"
#include <framework/ray.h>
#include <glm/glm.hpp>
#include <map>
#include <framework/trackball.h>

extern int extr_dof_samples;
extern float extr_dof_aperture;
extern float extr_dof_f;

std::vector<Ray> getEyeFrame(Trackball& camera);

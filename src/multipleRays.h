#pragma once
#include "common.h"

extern int rayMultiplier;
extern std::vector<Ray> rays; 

glm::vec3 calculateColor(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features, float x, float y, std::optional<glm::vec2> win = std::nullopt, std::optional<int> debug = 0);

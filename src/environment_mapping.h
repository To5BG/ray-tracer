#pragma once
#include "common.h"
#include <glm/glm.hpp>

extern bool extr_enabledSkyBox;
extern bool extr_enabledReflMap;

glm::vec3 environment_lookup(glm::vec3 v);
glm::vec3 environment_map(const Image& img, glm::vec3 texCoord, Features& feats);
glm::vec3 acquireTexelClamp(const Image& image, const glm::vec2& texCoord, const Features& feats);
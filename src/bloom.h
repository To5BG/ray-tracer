#pragma once
#include "common.h"

extern float threshold;
extern int bloomsize;
extern bool debugBloom;

void addBloom(std::vector<glm::vec3>& pixels, int width, int height);
glm::vec3 filterPixel(int idx, int width, int height);

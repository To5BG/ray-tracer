#pragma once
#include "common.h"

extern float threshold;
extern int bloomsize;
void addBloom(std::vector<glm::vec3>& pixels, int width, int height);

glm::vec3 filterPixel(std::vector<glm::vec3> source, glm::vec2 filter, int i, int j, glm::vec3 col);
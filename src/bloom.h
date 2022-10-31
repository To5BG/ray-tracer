#pragma once
#include "common.h"

extern float threshold; // color threshold
extern int bloomsize; // the size of the filter
extern bool debugBloom; // for visual debug
extern float sigma; // to compute standard deviation
extern int gaussian; // either box filter or gaussian filter
extern float scale; // bloom scale to make it brighter

std::vector<std::vector<float>> gaussianKernel(float sigma);
glm::vec3 filterPixel(int idx, int width, int height);
void addBloom(std::vector<glm::vec3>& pixels, int width, int height);

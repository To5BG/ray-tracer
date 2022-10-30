#include "bloom.h"
#include <algorithm> 
#include <iostream>
float threshold = 0.7f;
int bloomsize = 2;
std::vector<glm::vec3> filterPixels = {};


glm::vec3 filterPixel(int j)
{
    glm::vec3 avg = glm::vec3 { 0 };
    for (int i = 0; i < 2*bloomsize+1; i++) {
        avg += filterPixels[j + (i - bloomsize)] / (float)(2*bloomsize + 1);
    }
    return avg;
}

void addBloom(std::vector<glm::vec3>& pixels, int width, int height) {

    filterPixels.clear();
	for (int i = 0; i < pixels.size(); i++) {

		glm::vec3 color = pixels[i];

        if (color.x > threshold || color.y > threshold || color.z > threshold) {
            filterPixels.push_back(color);
        } else {
            filterPixels.push_back(glm::vec3{0});
        }

  	}

    for (int i = 0; i < pixels.size(); i++) {
     if (i > bloomsize && i < pixels.size() - bloomsize)
            pixels[i] += filterPixel(i);
    
    }
}

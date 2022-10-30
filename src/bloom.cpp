#include "bloom.h"
#include <algorithm> 
#include <iostream>

float threshold = 0.7f;
int bloomsize = 2;
std::vector<glm::vec3> filterPixels = {};
bool debugBloom = false;

//boxfilter size (bloomsize*2+1) for both width and height
glm::vec3 filterPixel(int index, int width, int height)
{
    // calculate the average for each pixel, when applying the bloom filter
    glm::vec3 avg = glm::vec3 { 0 };
    int filtersize = 2 * bloomsize + 1;
    // loop over all corresponding pixels
    for (int i = 0; i < filtersize; i++) {
        for (int j = 0; j < filtersize; j++) {
            int heightpos = (j - bloomsize)* height; 
            int widthpos = (i - bloomsize);
            // TODO check if the filter is not outside the range of the image, else consider the pixel black. This doesn't work yet
            if (heightpos >= 0 && heightpos <= height * width &&  widthpos >= 0 && widthpos <= width)
                avg += filterPixels[index + widthpos + heightpos] / (float)(pow(filtersize, 2));
            else
                avg += glm::vec3 { 0 };
        }
    }
    return avg;
}

void addBloom(std::vector<glm::vec3>& pixels, int width, int height) {

    // clear the filter
    filterPixels.clear();
    // loop over each pixel and check if the color has a value above threshold
    for (int i = 0; i < pixels.size(); i++) {

		glm::vec3 color = pixels[i];

        if (color.x > threshold || color.y > threshold || color.z > threshold) {
            filterPixels.push_back(color);
        } else {
            filterPixels.push_back(glm::vec3{0});
        }

  	}

    // loop over each pixel again, and apply the boxfilter
    for (int i = 0; i < pixels.size(); i++) {
        if (i > bloomsize && i < pixels.size() - bloomsize) {
            if (debugBloom) {
                pixels[i] = filterPixel(i, width, height);
            } else {
                pixels[i] += filterPixel(i, width, height);
            }
         }
    }
}

#include "bloom.h"
#include <algorithm> 
#include <iostream>

float threshold = 0.7f;
int bloomsize = 2;
float sigma = 1.0f;
std::vector<glm::vec3> filterPixels = {};
bool debugBloom = false;
int gaussian = 1;
float scale = 1.0f;

// cr
std::vector<std::vector<float>> gaussianKernel(float sigma)
{
    // standard deviation and size
    float sd = pow(sigma, 2);
    int size = bloomsize * 2 + 1;
    // create the 2D vector holding all values
    std::vector<std::vector<float>> gaus;
    // will be the total sum to devide every value by
    float denom = 0;
    const auto spread = 1.0 / (2 * sd);

    //create the gaussian kernel
    for (int i = 0; i < size; i++) {
        const auto x = i - bloomsize;
        float valx = std::exp(-x * x * spread);

        std::vector<float> gausx;
        for (int j = 0; j < size; j++) {
            const auto y = j - bloomsize;
            float valy = std::exp(-y * y * spread);
            float value = valy * valx;
            gausx.push_back(value);
            denom += valy * valx;
        }
        gaus.push_back(gausx);
    }

    // devide everything by the denominator
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            gaus[i][j] /= denom;
        }
    }
    return gaus;
}


//boxfilter size (bloomsize*2+1) for both width and height
glm::vec3 filterPixel(int index, int width, int height)
{
    // calculate the average for each pixel, when applying the bloom filter
    glm::vec3 avg = glm::vec3 { 0 };
    int filtersize = 2 * bloomsize + 1;

    std::vector<std::vector<float>> gaussianFilter;

    if (gaussian) {
        gaussianFilter = gaussianKernel(sigma);
    }
    // loop over all corresponding pixels
    for (int i = 0; i < filtersize; i++) {
        for (int j = 0; j < filtersize; j++) {
            int heightpos = (j - bloomsize)* height; 
            int widthpos = (i - bloomsize);

            // check if the range is correct
            int x = widthpos + (index % width);
            int y = index + heightpos + widthpos;
            if (x >= 0 && x < width && y >= 0 && y <height*width) {
             
                // do either gaussian or box filter
                if (gaussian) {
                    // multiply each pixel with corresponding gaussian index
                    avg += filterPixels[index + widthpos + heightpos] * gaussianFilter[j][i];
                } else {
                    // average each pixel (box filter)
                    avg += filterPixels[index + widthpos + heightpos] / (float)(pow(filtersize, 2));
                }
            }
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

        //check threshold
        if (color.x > threshold || color.y > threshold || color.z > threshold) {
            filterPixels.push_back(color);
        } else {
            filterPixels.push_back(glm::vec3{0.0f});
        }

  	}

    // loop over each pixel again, and apply the boxfilter
    for (int i = 0; i < pixels.size(); i++) {
        if (debugBloom) {
            // calculate each pixel and multiply with scale for extra brightness
            pixels[i] = filterPixel(i, width, height) * scale;
        } else {
            pixels[i] += filterPixel(i, width, height) * scale;
        }
    }
}

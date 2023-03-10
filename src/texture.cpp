#include "texture.h"
#include <framework/image.h>
#include <iostream>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    float xTex = texCoord.x;
    float yTex = texCoord.y;
    float i = std::fmod(std::round((xTex * image.width) - 0.5), image.width);
    float j = std::fmod(std::round(((1 - yTex) * image.height) - 0.5), image.height);
    int toIndex = (j * image.width + i);
    return image.pixels[toIndex];
}
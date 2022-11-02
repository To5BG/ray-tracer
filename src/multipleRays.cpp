#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "multipleRays.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <iostream>
int rayMultiplier = 1;
std::vector<Ray> rays; 

// very long function signature
glm::vec3 calculateColor(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features, float x, float y, std::optional<glm::vec2> win, std::optional<int> debug)
{
    // calculate window resolution, if in debug mode, it needs window resulution somehow, so that's why i have an optional check for debug
    glm::vec2 windowResolution;
    if (std::move(*debug)) {
        rays.clear();
        windowResolution = std::move(*win);
    } else {
        windowResolution = screen.resolution();
    }

    // calculate the color
    glm::vec3 color = glm::vec3{0};
    for (int i = 0; i < rayMultiplier; i++) {
        for (int j = 0; j < rayMultiplier; j++) {

            // create random offset
            float randomx = (float)((rand() % 100) / 100.0f);
            float randomy = (float)((rand() % 100) / 100.0f);

            // calc the pixelpos
            glm::vec2 normalizedPixelPos {
                float(x + ((i) / float(rayMultiplier) + randomx)) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y + ((j) / float(rayMultiplier) + randomy)) / float(windowResolution.y) * 2.0f - 1.0f
            };

            // calculate the ray
            Ray cameraRay = camera.generateRay(normalizedPixelPos);
            // if debug mode, push the ray to the vector which holds the rays, so we can draw them every frame
            if (std::move(*debug))  rays.push_back(cameraRay);
            color += getFinalColor(scene, bvh, cameraRay, features);
        }
    }

    // calculate the average color
    color /= pow(rayMultiplier, 2.0f);

    return color;

}
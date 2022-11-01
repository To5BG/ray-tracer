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


glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    
   
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            
            // TODO: put your own implementation of recursive ray tracing here.
            if (rayDepth > 0 && hitInfo.material.ks != glm::vec3({0,0,0})) {
                drawRay(ray, glm::vec3(1.0f));
                Ray reflection = computeReflectionRay(ray, hitInfo);
                
                return hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth - 1);
            } else {
                drawRay(ray, Lo);

        // Set the color of the pixel to white if the ray hits.
        return Lo;
            }
        }

        // Draw a white debug ray if the ray hits.
        if (features.enableShading) {
            drawRay(ray, Lo);
        }else {
            drawRay(ray, glm::vec3(1.0f));
        }

        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);

            glm::vec3 color = getFinalColor(scene, bvh, cameraRay, features);

            if (features.extra.enableMultipleRaysPerPixel) {

                 const glm::vec2 normalizedPixelPos1 {
                    float(x - 0.5f) / float(windowResolution.x) * 2.0f - 1.0f,
                    float(y -0.5f) / float(windowResolution.y) * 2.0f - 1.0f
                };

                  const glm::vec2 normalizedPixelPos2 {
                     float(x + 0.5f) / float(windowResolution.x) * 2.0f - 1.0f,
                     float(y + 0.5f) / float(windowResolution.y) * 2.0f - 1.0f
                 };

                const Ray cameraRay1 = camera.generateRay(normalizedPixelPos1);
                const Ray cameraRay2 = camera.generateRay(normalizedPixelPos2);

                glm::vec3 color1 = getFinalColor(scene, bvh, cameraRay1, features);
                glm::vec3 color2 = getFinalColor(scene, bvh, cameraRay2, features);
                color = (color1 + color2)/2.0f;
            }

            screen.setPixel(x, y, color);
        }
    }
}
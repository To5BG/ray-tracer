#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "bloom.h"
#include "multipleRays.h"
#include <framework/trackball.h>
#include "gloss.h"
#include "dof.h"
#ifdef NDEBUG
#include <omp.h>
#endif
#include <iostream>

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    float t = ray.t;
   
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.extra.enableTransparency) {
            drawRay(ray, Lo);
            if (hitInfo.material.transparency != 1.0) {
                Ray newRay = {};
                newRay.direction = ray.direction;
                newRay.origin = (ray.t+0.000001f) * ray.direction + ray.origin;
                newRay.t = std::numeric_limits<float>::max();
                return hitInfo.material.transparency * Lo + (1.0f - hitInfo.material.transparency) * getFinalColor(scene, bvh, newRay, features, rayDepth);
               
            } 
        }

        if (features.enableRecursive || features.extra.enableGlossyReflection) {
            
            if (rayDepth > 0 && hitInfo.material.ks != glm::vec3({ 0, 0, 0 })) {
                drawRay(ray, glm::vec3 { 1.0f });
                Ray reflection = computeReflectionRay(ray, hitInfo);
                if (features.extra.enableGlossyReflection) {
                    glm::vec3 avg = glm::vec3 { 0.0f };
                    std::vector<Ray> rays = getRaySamples(hitInfo, reflection);
                    std::for_each(rays.begin(), rays.end(), [&](Ray r) {
                        avg += getFinalColor(scene, bvh, r, features, rayDepth - 1);
                    });
                    return hitInfo.material.ks * avg / float(extr_glossy_filterSize);
                }
                return hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth - 1);
            }
        }
        // Draw a white debug ray if the ray hits, but if shading/multiple rays make the ray the same color as the hit
        drawRay(ray, features.enableShading || features.extra.enableMultipleRaysPerPixel ? Lo : glm::vec3 { 1.0f });
        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } 
    // Draw a red debug ray if the ray missed.
    drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
    // Set the color of the pixel to black if the ray misses.
    return glm::vec3(0.0f);
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
           
            glm::vec3 color;

            // calculate the ray, (multiple if boolean set to true)
            if (features.extra.enableMultipleRaysPerPixel) {
                 color = calculateColor(scene, camera, bvh, screen, features, float(x), float(y) , std::nullopt, 0);
            } else {
                const glm::vec2 normalizedPixelPos {
                    float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                    float(y) / float(windowResolution.y) * 2.0f - 1.0f
                };
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                color = getFinalColor(scene, bvh, cameraRay, features);
            }

            screen.setPixel(x, y, color);
        }
    }
    if (features.extra.enableBloomEffect) {
        addBloom(screen.pixels(), windowResolution.x, windowResolution.y);
    }
}
#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "bloom.h"
#include <framework/trackball.h>
#include "gloss.h"
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
                if (features.enableRecursive) {
                    return hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth - 1);
                }
                int i = 0;
                glm::vec3 avg = glm::vec3 { 0.0f };
                //AxisAlignedBox sq = getBlurSquare(reflection);
                //drawAABB(sq, DrawMode::Wireframe, { 1.0f, 1.0f, 1.0f }, 0.5f);
                //glm::vec3 U = sq.lower;
                glm::vec3 U = ray.origin + ray.direction - 1.0f / hitInfo.material.shininess;
                glm::vec3 V = ray.origin + ray.direction + 1.0f / hitInfo.material.shininess;
                for (float u = -glossy_filter.filterSize; u <= glossy_filter.filterSize; u++) {
                    for (float v = -glossy_filter.filterSize; v <= glossy_filter.filterSize; v++) {
                // 
                //        glm::vec3 V = sq.upper;
                // 
                        Ray pertr_refl = { reflection.origin, reflection.direction + u * U + v * V, reflection.t };
                        i++;
                        avg += (hitInfo.material.ks * glossy_filter.kernel[i] * getFinalColor(scene, bvh, pertr_refl, features, rayDepth - 1));
                    }
                }
                return avg;
            }
            drawRay(ray, Lo);
            // Set the color of the pixel to white if the ray hits.
            return Lo;
        }
        // Draw a white debug ray if the ray hits.
        drawRay(ray, features.enableShading ? Lo : glm::vec3 { 1.0f });
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
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
            
        }
    }
    if (features.extra.enableBloomEffect) {
        addBloom(screen.pixels(), windowResolution.x, windowResolution.y);
    }
}
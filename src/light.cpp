#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>

int samplesPerUnit = 10;

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    float distance = glm::distance(segmentLight.endpoint0, segmentLight.endpoint1);
    float patition0 = 0.5f;
    float patition1 = 0.5f;
   
    if (distance != 0) {
        patition1 = glm::distance(segmentLight.endpoint0, position) / distance;
        patition0 = 1 - patition1;
    }

    color = patition0 * segmentLight.color0 + patition1 * segmentLight.color1;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    position = glm::vec3(0.0);
    color = glm::vec3(0.0);
    // TODO: implement this function.
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    // to prevent the ray from interacting with the object at the origin
    float epsilon = 0.00001f;

    if (features.enableHardShadow || features.enableSoftShadow) {

        // the intersection point to calculate the visibility for, minus a really small offset for floating point issues
        glm::vec3 intersectPos = ray.origin + (1-epsilon) * ray.t * ray.direction;

        // calculate the light ray 
        Ray light = Ray {};
        light.origin = intersectPos;
        light.direction = samplePos - intersectPos;
        light.t = 1.0f;
        
        // if there is an object between the light and the hitpos
        if (bvh.intersect(light, hitInfo, features)) {
            drawRay(light, { 1.0f, 0.0f, 0.0f });
            return 0.0f;
        }

        // draw a ray in the lightcolor if no intersection is found
        drawRay(light, debugColor);

    }
    return 1.0;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{

    // for hard shadow this will indiciate if we should light the hitpoint or not
    float lighted = features.enableHardShadow ? 0.0f : 1.0f;

    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.
        glm::vec3 shading = glm::vec3 { 0.0f };

        // iterate over all lights
        for (auto& light : scene.lights) {
            
            if (std::holds_alternative<PointLight>(light)) {

                // pointLight
                PointLight pointLight = std::get<PointLight>(light);

                // add shadow
                lighted = features.enableHardShadow ? testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo) : 1.0f;
               
                // shade
                shading += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo) * lighted;

            } else if (std::holds_alternative<SegmentLight>(light)) {

                if (features.enableSoftShadow) {

                    // segmentLight
                    const SegmentLight segmentLight = std::get<SegmentLight>(light);

                    // the amount of samples per distance, multiplied by the distance to get the total amount of needed samples
                    float samples = (float)std::max((float)samplesPerUnit * glm::distance(segmentLight.endpoint0, segmentLight.endpoint1),1.0f);
                
                    // sets the seed so that every time we run the for loop we get the same result for rand() (otherwise there is noise)
                    srand(1);

                    // create the samples
                    for (int i = 0; i < samples; i++) {

                        glm::vec3 color;

                        // random offset [-0.5,0.5)
                        float randomUniform = (float)(((rand() % 100) - 50) / 100.0f);

                        // randomize the sample position
                        glm::vec3 segmentDistance = ((segmentLight.endpoint1 - segmentLight.endpoint0) / samples);
                        glm::vec3 randomAddition = segmentDistance * randomUniform;

                        // get the position based on sampling
                        glm::vec3 interpolatedPos = segmentLight.endpoint0 * float(1 - i / samples) + segmentLight.endpoint1 * float(i / samples);
                        glm::vec3 position = interpolatedPos - randomAddition;
                    
                        sampleSegmentLight(segmentLight, position, color);
                        lighted = testVisibilityLightSample(position, color, bvh, features, ray, hitInfo);

                        // generate a ray for the sample
                        Ray sampleRay = Ray {};
                        sampleRay.origin = position;
                        sampleRay.direction = ray.direction * ray.t + ray.origin - position;
                        sampleRay.t = 1;

                        // check if there is no shadow
                        if (lighted)
                            drawRay(sampleRay, color);

                        // shade each sample, and divide by the total amount of samples 
                        shading += computeShading(position, color, features, ray, hitInfo) / samples * lighted;
                    }
                }


            } else if(std::holds_alternative<ParallelogramLight>(light)) {

                // parallelogramLight

                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                shading += computeShading(parallelogramLight.v0, parallelogramLight.color0, features, ray, hitInfo);
                shading += computeShading(parallelogramLight.v0 + parallelogramLight.edge01, parallelogramLight.color1, features, ray, hitInfo);
                shading += computeShading(parallelogramLight.v0 + parallelogramLight.edge02, parallelogramLight.color2, features, ray, hitInfo);
                shading += computeShading(parallelogramLight.v0 + parallelogramLight.edge01 + parallelogramLight.edge02, parallelogramLight.color3, features, ray, hitInfo);
            }
        }

        return shading;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd * lighted;
    }
}

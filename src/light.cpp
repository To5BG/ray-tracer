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
float sampleCountRec;
float segmentLengthRec;
glm::vec3 sampleDistanceSegment;

int samplesPerUnitParallel = 10;
float dirISamplesRec;
float dirJSamplesRec;

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    // random offset [0,1)
    float randomUniform = (float)((rand() % 50) / 50.0f);
    // randomize the sample position
    glm::vec3 randomAddition = sampleDistanceSegment * randomUniform;
    // get the position based on sampling
    float alpha = position.x * sampleCountRec;
    glm::vec3 interpolatedPos = segmentLight.endpoint0 * (1 - alpha) + segmentLight.endpoint1 * alpha;
    position = interpolatedPos + randomAddition;
    
    // interpolation
    float patition0 = 0.5f;
    float patition1 = 0.5f;
    // if the endpoints are at the same position make the lights both 50%
    if (segmentLengthRec != 0) {
        patition1 = glm::distance(segmentLight.endpoint0, position) * segmentLengthRec;
        patition0 = 1 - patition1;
    }
    // calculate the color
    color = patition0 * segmentLight.color0 + patition1 * segmentLight.color1;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    // random position in sample 
    float randomUniformI = (float)(((rand() % 50)) / 50.0f) + position.x;
    float randomUniformJ = (float)(((rand() % 50)) / 50.0f) + position.y;
    // calculate the partitions for each color
    float x = randomUniformI * dirISamplesRec;
    float y = randomUniformJ * dirJSamplesRec;
    // calculate offset from both edges
    glm::vec3 Ipos = parallelogramLight.edge01 * x;
    glm::vec3 Jpos = parallelogramLight.edge02 * y;

    // calc position
    position = parallelogramLight.v0 + Ipos + Jpos;
    // bilinearly interpolate the color
    color = (1 - y) * (x * parallelogramLight.color1 + (1 - x) * parallelogramLight.color0) + 
        y * (x * parallelogramLight.color3 + (1 - x) * parallelogramLight.color2);
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
    // If shading is enabled, compute the contribution from all lights.
    glm::vec3 shading = glm::vec3 { 0.0f };
    // iterate over all lights
    for (auto& light : scene.lights) {
            
        if (std::holds_alternative<PointLight>(light)) {
            // pointLight
            PointLight pointLight = std::get<PointLight>(light);
            // add shadow
            lighted = features.enableHardShadow 
                ? testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo) : 1.0f;
            // shade
            glm::vec3 shade = features.enableShading 
                ? computeShading(pointLight.position, pointLight.color, features, ray, hitInfo) : hitInfo.material.kd;
            shading += shade * lighted;

        } else if (std::holds_alternative<SegmentLight>(light)) {

            if (features.enableSoftShadow || features.enableShading) {
                // segmentLight
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                // length of segment
                segmentLengthRec = 1.0f / glm::distance(segmentLight.endpoint0, segmentLight.endpoint1);
                // the amount of samples per distance, multiplied by the distance to get the total amount of needed samples
                sampleCountRec = 1.0f / std::floor(std::fmax(samplesPerUnit / segmentLengthRec, 1.0f));
                // distance between two samples
                sampleDistanceSegment = ((segmentLight.endpoint1 - segmentLight.endpoint0) * sampleCountRec);
                // sets the seed so that every time we run the for loop we get the same result for rand() (otherwise there is noise)
                if (enableDebugDraw) srand(1);
                // create the samples
                float sampleCount = 1.0f / sampleCountRec;
                for (int i = 0; i < sampleCount; i++) {
                    // initialize color and position
                    glm::vec3 color;
                    glm::vec3 position;
                    // pass the index of the sample
                    position.x = i;
                    // create samples for position and color
                    sampleSegmentLight(segmentLight, position, color);
                    // check if in shadow
                    lighted = features.enableSoftShadow ? testVisibilityLightSample(position, color, bvh, features, ray, hitInfo) : 1.0f;
                    // check if there is no shadow
                    if (lighted)
                        // dray ray for sample
                        drawRay(Ray { position, ray.direction * ray.t + ray.origin - position, 1 }, color);
                    // shade
                    glm::vec3 shade = features.enableShading 
                        ? computeShading(position, color, features, ray, hitInfo) * sampleCountRec 
                        : hitInfo.material.kd;
                    shading += shade * lighted;
                }
            }

        } else if(std::holds_alternative<ParallelogramLight>(light)) {

            if (features.enableSoftShadow || features.enableShading) {
                // parallelogramLight
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                // amount of samples over both edges
                dirISamplesRec = 1.0f / std::floor(std::fmax(samplesPerUnitParallel * glm::length(parallelogramLight.edge01), 1.0f));
                dirJSamplesRec = 1.0f / std::floor(std::fmax(samplesPerUnitParallel * glm::length(parallelogramLight.edge02), 1.0f));
                // sets the seed so that every time we run the for loop we get the same result for rand() (otherwise there is noise)
                if (enableDebugDraw) srand(1);
                // iterate over both axis of the light
                float dirISampleCount = 1.0f / dirISamplesRec;
                float dirJSampleCount = 1.0f / dirJSamplesRec;
                for (int i = 0; i < dirISampleCount; i++) {
                    for (int j = 0; j < dirJSampleCount; j++) {
                        // create pointers for position and color
                        glm::vec3 position = { i, j, 0.0f };
                        glm::vec3 color;
                        // sample
                        sampleParallelogramLight(parallelogramLight, position, color);
                        // calculate the light
                        lighted = features.enableSoftShadow ? testVisibilityLightSample(position, color, bvh, features, ray, hitInfo) : 1.0f;
                        // check if there is no shadow
                        if (lighted)
                            // generate a ray for the sample
                            drawRay(Ray { position, ray.direction * ray.t + ray.origin - position, 1.0f }, color);
                        // shade
                        glm::vec3 shade = features.enableShading 
                            ? computeShading(position, color, features, ray, hitInfo) * dirISamplesRec * dirJSamplesRec 
                            : hitInfo.material.kd;
                        shading += shade * lighted;
                    }
                }
            }
        }
    }
    // return the shading amount
    return shading;
}

#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    float totalArea = glm::length(glm::cross((v0 - v2), (v1 - v2))) / 2;
    float alpha = (glm::length(glm::cross(p - v2, v1 - v2))) / (2 * totalArea);
    float beta = (glm::length(glm::cross(p - v2, v0 - v2))) / (2 * totalArea);
    float gamma = (glm::length(glm::cross(p - v1, v0 - v1))) / (2 * totalArea);
    if (alpha < 0 || beta < 0 || alpha + beta > 1.0 || alpha + gamma > 1.0 || beta + gamma > 1.0) {
        return false;
    } else {
        return true;
    }
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / glm::dot(ray.direction, plane.normal);
    glm::vec3 point = ray.origin + t * ray.direction;
    if (t < 0 || glm::dot(plane.normal, ray.direction) == 0) {
        return false;
    } else if (t < ray.t && t >= 0) {
        ray.t = t;
        return true;
    }
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    float bottom = glm::length(glm::cross((v0 - v2), (v1 - v2)));
    if (bottom == 0) {
        bottom = 1, 0;
    }
    glm::vec3 normal = glm::cross((v0 - v2), (v1 - v2)) / bottom;
    plane.normal = normal;
    float D = glm::dot(normal, v0);
    plane.D = D;
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    float oldT = ray.t;
    Plane plane = trianglePlane(v0, v1, v2);
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / glm::dot(ray.direction, plane.normal);
    glm::vec3 point = ray.origin + t * ray.direction;
    if (!intersectRayWithPlane(plane, ray)) {
        ray.t = oldT;
        return false;
    } else {
        if (pointInTriangle(v0, v1, v2, plane.normal, point)) {
            if (ray.t < t) {
                ray.t = oldT;
                return false;
            } else {
                ray.t = t;
                return true;
            }
        } else {
            ray.t = oldT;
            return false;
        }
    }
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    float oldT = ray.t;
    ray.origin = ray.origin - sphere.center;
    float A = glm::pow(ray.direction.x, 2) + glm::pow(ray.direction.y, 2) + glm::pow(ray.direction.z, 2);
    float B = 2 * (ray.direction.x * ray.origin.x + ray.direction.y * ray.origin.y + ray.direction.z * ray.origin.z);
    float C = glm::pow(ray.origin.x, 2) + glm::pow(ray.origin.y, 2) + glm::pow(ray.origin.z, 2) - glm::pow(sphere.radius, 2);
    float discriminant = glm::pow(B, 2) - 4 * A * C;

    if (discriminant < 0) {
        ray.origin = ray.origin + sphere.center;
        return false;
    } else {
        if (discriminant == 0) {
            float uniqueT = -B / (2 * A);
            if (ray.t < uniqueT) {
                ray.t = oldT;
                ray.origin = ray.origin + sphere.center;
                return false;
            } else {
                ray.t = uniqueT;
                ray.origin = ray.origin + sphere.center;
                return true;
            }
        } else {
            if (discriminant > 0) {
                float t0 = (-B + glm::sqrt(discriminant)) / 2 * A;
                float t1 = (-B - glm::sqrt(discriminant)) / 2 * A;
                if (t1 < 0) {
                    if (t0 < oldT && t0 > 0) {
                        ray.t = t0;
                        ray.origin = ray.origin + sphere.center;
                        return true;
                    } else {
                        ray.t = oldT;
                        ray.origin = ray.origin + sphere.center;
                        return false;
                    }
                } else {
                    float t = std::min(t0, t1);
                    if (t < ray.t && t > 0) {
                        ray.t = t;
                        ray.origin = ray.origin + sphere.center;
                        return true;
                    } else {
                        ray.t = oldT;
                        ray.origin = ray.origin + sphere.center;
                        return false;
                    }
                }
            }
        }
    }
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float oldT = ray.t;
    float txmin = (box.lower.x - ray.origin.x) / ray.direction.x;
    float txmax = (box.upper.x - ray.origin.x) / ray.direction.x;
    float tymin = (box.lower.y - ray.origin.y) / ray.direction.y;
    float tymax = (box.upper.y - ray.origin.y) / ray.direction.y;
    float tzmin = (box.lower.z - ray.origin.z) / ray.direction.z;
    float tzmax = (box.upper.z - ray.origin.z) / ray.direction.z;
    float tinx = std::min(txmin, txmax);
    float toutx = std::max(txmin, txmax);
    float tiny = std::min(tymin, tymax);
    float touty = std::max(tymin, tymax);
    float tinz = std::min(tzmin, tzmax);
    float toutz = std::max(tzmin, tzmax);
    float tin = std::max(tinx, std::max(tiny, tinz));
    float tout = std::min(toutx, std::min(touty, toutz));
    if (tin < 0) {
        if (tout < 0) {
            ray.t = oldT;
            return true;
        } else {
            ray.t = tout;
            return true;
        }
    }
    if (tin > tout || tout < 0) {
        return false;
    } else {
        if (oldT < tin) {
            ray.t = oldT;
            return false;
        } else {
            ray.t = tin;
            return true;
        }
    }
}
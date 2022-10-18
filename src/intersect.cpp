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
#include <array>

bool isZero(float a, float epsilon = 0.000001f)
{
    return std::fabs(a) <= epsilon;
}

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    // Check if point lies on plane
    if (!isZero(glm::dot(n, p - v0)))
        return false;
    glm::vec3 a = glm::cross(p - v2, p - v0);
    glm::vec3 b = glm::cross(p - v1, p - v2);
    glm::vec3 c = glm::cross(p - v0, p - v1);
    // Check if all three cross-products point in the same direction
    return glm::dot(a, b) >= 0 && 
        glm::dot(b, c) >= 0 && 
        glm::dot(a, c) >= 0;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    // Check if ray is parallel with/lies on the plane
    if (isZero(glm::dot(plane.normal, ray.direction)))
        return false;
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / glm::dot(ray.direction, plane.normal);
    // Check if t is positive
    bool valid = t > 0;
    ray.t = valid ? t : ray.t;
    return valid;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 n = glm::cross(v2 - v0, v1 - v0);
    // Choose arbitrary normal if triangle is degenerate
    plane.normal = isZero(glm::length(n)) ? glm::vec3 { 1.0f, 0.0f, 0.0f } : glm::normalize(n);
    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    Plane p = trianglePlane(v0, v1, v2);
    float oldT = ray.t;
    if (!intersectRayWithPlane(p, ray))
        return false;
    // Update t if it lies on triangle and is smaller than curr ray.t
    bool pit = pointInTriangle(v0, v1, v2, p.normal, ray.origin + ray.direction * ray.t) && ray.t < oldT;
    ray.t = pit ? ray.t : oldT;
    if(pit) hitInfo.normal = -p.normal;
    return pit;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    float a = glm::dot(ray.direction, ray.direction);
    float b = 2.0f * glm::dot(ray.direction, ray.origin - sphere.center);
    float c = glm::dot(ray.origin - sphere.center, ray.origin - sphere.center) - pow(sphere.radius, 2.0f);
    float D = b * b - 4 * a * c;
    // Check if ray misses sphere or origin lies on the sphere's surface
    if (D < 0 || isZero(c))
        return false;
    float t1 = (-b - sqrt(D)) / (2 * a);
    float t2 = (-b + sqrt(D)) / (2 * a);
    // Check if ray points opposite to sphere
    if (t1 < 0 && t2 < 0)
        return false;
    // newT is either max value if inside of sphere, or min value if outside
    float newT = t1 * t2 < 0 ? std::fmax(t1, t2) : std::fmin(t1, t2);
    // Check if t is smaller than oldT
    bool valid = newT < ray.t;
    ray.t = valid ? newT : ray.t;
    return valid;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    // slight optimization by applying divisions once
    glm::vec3 invDir = { 1 / ray.direction.x, 1 / ray.direction.y, 1 / ray.direction.z };
    // Store all Ts, even idx -> min, odd idx -> max
    std::array<float, 6> ts = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    // For each dimension (i) -> set min and max to either float.min and float.max 
    // if ray does not intersect dim plane, or to closer vertex based on ray dir
    for (int i = 0; i < 3; i++) {
        bool dir = ray.direction[i] > 0;
        //Min
        ts[2 * i] = isZero(ray.direction[i]) ? std::numeric_limits<float>::min() 
            : ((dir ? box.lower[i] : box.upper[i]) - ray.origin[i]) * invDir[i];
        //Max
        ts[2 * i + 1] = isZero(ray.direction[i]) ? std::numeric_limits<float>::max()
            : ((dir ? box.upper[i] : box.lower[i]) - ray.origin[i]) * invDir[i];
    }
    // Calculate t_in and t_out
    float t_in = glm::max(glm::max(ts[0], ts[2]), ts[4]);
    float t_out = glm::min(glm::min(ts[1], ts[3]), ts[5]);
    // Check if ray misses AABB, if ray is pointing opposite to AABB, or if origin lies on AABB
    if (t_in > t_out || t_out <= 0 || isZero(t_in))
        return false;
    // Check if origin inside of AABB
    ray.t = t_in < 0 ? t_out : t_in;
    return true;
}

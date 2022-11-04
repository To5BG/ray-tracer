#include "dof.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <random>
#include "light.h"
#include <framework/opengl_includes.h>

int extr_dof_samples = 3;
float extr_dof_aperture = 0.005f;
float extr_dof_f = 5;
bool draw_random_rays = true;
std::normal_distribution<> nnormal(0.0, 1.0f);
std::mt19937 mtGenn(std::random_device {}());

std::vector<Ray> DOFrays = {};
std::vector<Ray> CamRays = {};

std::vector<Ray> getEyeFrame(Ray& ray)
{
    // Calculate square basis
    glm::vec3 w = glm::normalize(ray.direction);
    // Follow textbook's suggestion - replace smallest magnitude 
    // dimension with 1 to get a vector sufficiently different from w
    int minTerm = fabs(w.x) <= fabs(w.y) && fabs(w.x) <= fabs(w.z) ? 0 : (fabs(w.y) <= fabs(w.x) && fabs(w.y) <= fabs(w.z) ? 1 : 2);
    glm::vec3 t = w;
    t[minTerm] = 1;
    // Better performance, but edge case - if incident angle is exactly 90, it will result in a zero vector
    // glm::vec3 u = glm::normalize(glm:cross(hitInfo.normal, w));
    glm::vec3 u = glm::normalize(glm::cross(t, w));
    glm::vec3 v = glm::cross(w, u);

    std::vector<Ray> rays;
    rays.push_back(ray);
    // 0,0 on distribution space is 0.5, 0.5 on square's space -> offset by center
    float side = extr_dof_aperture;
    float offset = -side / 2.0f;
    glm::vec3 focusPoint = ray.origin + w * extr_dof_f;

    CamRays.clear();
    CamRays.push_back(Ray { ray.origin + glm::vec3 { extr_dof_aperture, extr_dof_aperture, 0 }, focusPoint - ray.origin, std::numeric_limits<float>::max() });
    CamRays.push_back(Ray { ray.origin + glm::vec3 { -extr_dof_aperture, extr_dof_aperture, 0 }, focusPoint - ray.origin, std::numeric_limits<float>::max() });
    CamRays.push_back(Ray { ray.origin + glm::vec3 { -extr_dof_aperture, -extr_dof_aperture, 0 }, focusPoint - ray.origin, std::numeric_limits<float>::max() });
    CamRays.push_back(Ray { ray.origin + glm::vec3 { extr_dof_aperture, -extr_dof_aperture, 0 }, focusPoint - ray.origin, std::numeric_limits<float>::max() });

    for (int i = 0; i < extr_dof_samples; i++) {
        glm::vec3 origin = ray.origin + float(offset + nnormal(mtGenn) * side) * u + float(offset + nnormal(mtGenn) * side) * v;
        // Make sure perturbed ray is not below the surface its reflected from (at large sigmas and/or small incident angle)
        rays.push_back({ origin, glm::normalize(focusPoint - origin), std::numeric_limits<float>::max() });
    }
    return rays;
}

#include "dof.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <random>
#include "light.h"
#include <framework/opengl_includes.h>

int extr_dof_samples = 3;
float extr_dof_aperture = 1.4f;
float extr_dof_f = 2.5f;
float extr_dof = 0.0f;

bool draw_random_rays = true;
std::uniform_real_distribution<> uniform(0.0, 1.0f);
std::mt19937 mtGenn(std::random_device {}());

std::vector<Ray> dofRays;
std::vector<Ray> camFrame = std::vector<Ray>(4);

std::vector<Ray> getEyeFrame(Ray& ray, const Trackball& camera)
{
    // Calculate square basis
    glm::vec3 w = glm::normalize(ray.direction);
    // Follow textbook's suggestion - replace smallest magnitude 
    // dimension with 1 to get a vector sufficiently different from w
    //int minTerm = fabs(w.x) <= fabs(w.y) && fabs(w.x) <= fabs(w.z) ? 0 : (fabs(w.y) <= fabs(w.x) && fabs(w.y) <= fabs(w.z) ? 1 : 2);
    //glm::vec3 t = w;
    //t[minTerm] = 1;
    // Better performance, but edge case - if incident angle is exactly 90, it will result in a zero vector
    // glm::vec3 u = glm::normalize(glm:cross(hitInfo.normal, w));
    glm::vec3 u = camera.up();
    glm::vec3 v = camera.left();

    std::vector<Ray> rays;
    float side = extr_dof_f / (extr_dof_aperture * 6.0f);
    // 0,0 on distribution space is 0.5, 0.5 on square's space -> offset by center
    float offset = -side / 2.0f;
    glm::vec3 focusPoint = ray.origin + w * extr_dof_f;
    
    camFrame[0] = { ray.origin + side * (u + v), glm::normalize(focusPoint - ray.origin - side * (u + v)) };
    camFrame[1] = { ray.origin + side * (v - u), glm::normalize(focusPoint - ray.origin - side * (v - u)) };
    camFrame[2] = { ray.origin + side * (-u - v), glm::normalize(focusPoint - ray.origin - side * (-u - v)) };
    camFrame[3] = { ray.origin + side * (u - v), glm::normalize(focusPoint - ray.origin - side * (u - v)) };

    for (int i = 0; i < extr_dof_samples; i++) {
        glm::vec3 origin = ray.origin + float(offset + uniform(mtGenn) * side) * u + float(offset + uniform(mtGenn) * side) * v;
        rays.push_back({ origin, glm::normalize(focusPoint - origin) });
    }
    return rays;
}

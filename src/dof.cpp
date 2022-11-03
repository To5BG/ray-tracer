#include "dof.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <random>
#include "light.h"
#include <framework/opengl_includes.h>

int extr_dof_samples = 1;
float extr_dof_aperture = 0.03f;
float extr_dof_f = 5;
std::normal_distribution<> nnormal(0.0, 1.0f);
std::mt19937 mtGenn(std::random_device {}());

std::vector<Ray> getEyeFrame(Trackball& camera)
{
    glm::vec3 w = glm::normalize(camera.lookAt());
    glm::vec3 u = glm::normalize(camera.left());
    glm::vec3 v = glm::normalize(camera.up());

    std::vector<Ray> rays;
    // 0,0 on distribution space is 0.5, 0.5 on square's space -> offset by center
    float side = extr_dof_aperture;
    float offset = -side / 2.0f;
    glm::vec3 focusPoint = camera.position() + w * extr_dof_f;
    for (int i = 0; i < extr_dof_samples; i++) {
        glm::vec3 origin = camera.position() + float(offset + nnormal(mtGenn) * side) * u + float(offset + nnormal(mtGenn) * side) * v;
        // Make sure perturbed ray is not below the surface its reflected from (at large sigmas and/or small incident angle)
        rays.push_back({ origin, glm::normalize(focusPoint - origin), std::numeric_limits<float>::max() });
    }
    return rays;
}

#include "light.h"
#include "environment_mapping.h"
#include <glm/geometric.hpp>

bool extr_enabledSkyBox = false;
bool extr_enabledReflMap = false;

    // Find the texCoords and index of face for cube mapping 
glm::vec3 environment_lookup(glm::vec3 v)
{
    glm::vec3 absV = glm::abs(v);
    float maxA, ut, vt, idx;

    // Go through each box to figure out for which box to take coords
    if (absV.x >= absV.y && absV.x >= absV.z) {
        // multiplied by two for simpler calculations later
        maxA = 2 * absV.x;
        ut = (v.x >= 0) ? -v.z : v.z;
        vt = -v.y;
        idx = (v.x < 0);
    } else if (absV.y >= absV.z) {
        maxA = 2 * absV.y;
        ut = v.x;
        vt = (v.y >= 0) ? v.z : -v.z;
        idx = (v.y < 0) + 2;
    } else {
        maxA = 2 * absV.z;
        ut = (v.z >= 0) ? v.x : -v.x;
        vt = -v.y;
        idx = (v.z < 0) + 4;
    }

    // Invert x planes due to project's config
    if (idx <= 1)
        idx = 1 - idx;

    return glm::vec3 { - (ut / maxA) + 0.5, - (vt / maxA) + 0.5, idx };
}

glm::vec3 acquireTexelClamp(const Image& image, const glm::vec2& texCoord, const Features& feats)
{
    int i = texCoord.x * image.width;
    int j = std::clamp((1.0f - texCoord.y) * image.height, 0.0f, image.height - 1.0f);
    return image.pixels[j * image.width + i];
}

// Reverse operation to environment_lookUp
glm::vec3 environment_map(glm::vec3 texCoord) 
{
    float ut = 2.0f * texCoord[0] - 1.0f;
    float vt = 2.0f * texCoord[1] - 1.0f;
    glm::vec3 res = {};
    int idx = texCoord[2];
    switch (idx) {
    case 0:
        res = { 1.0f, vt, -ut };
        break;
    case 1:
        res = { -1.0f, vt, ut };
        break;
    case 2:
        res = { ut, 1.0f, -vt };
        break;
    case 3:
        res = { ut, -1.0f, vt };
        break;
    case 4:
        res = { ut, vt, 1.0f };
        break;
    case 5:
        res = { -ut, vt, -1.0f };
        break;
    }
    return res;
}
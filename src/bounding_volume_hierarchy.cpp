#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <numeric>
#include <deque>

AxisAlignedBox calculateAABB(std::vector<Prim>& prims, std::vector<int>& prim_ids) 
{
    glm::vec3 min = glm::vec3 { std::numeric_limits<float>::max() };
    glm::vec3 max = glm::vec3 { std::numeric_limits<float>::min() };
    std::for_each(prim_ids.begin(), prim_ids.end(), [&](int i) {
        Prim p = prims[i];
        min = { std::fmin(min.x, p.min.x), std::fmin(min.y, p.min.y), std::fmin(min.z, p.min.z) };
        max = { std::fmax(max.x, p.max.x), std::fmax(max.y, p.max.y), std::fmax(max.z, p.max.z) };
    });
    return { min, max };
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    //Initial values
    this->m_numLeaves = 0;
    this->m_numLevels = 0;
    std::vector<BVHNode> nodes;

    // Create prim vector
    std::vector<Prim> prims;
    for (int i = 0; i < pScene->meshes.size(); i++) 
    {
        Mesh& mesh = pScene->meshes[i];
        for (int j = 0; j < mesh.triangles.size(); j++) {
            glm::uvec3 t = mesh.triangles[j];
            glm::vec3 v0 = mesh.vertices[t.x].position,
                      v1 = mesh.vertices[t.y].position,
                      v2 = mesh.vertices[t.z].position;
            prims.push_back(Prim { /*std::vector { v0, v1, v2 },*/ {
                    std::fmin(std::fmin(v0.x, v1.x), v2.x), std::fmin(std::fmin(v0.y, v1.y), v2.y), std::fmin(std::fmin(v0.z, v1.z), v2.z)
                }, { 
                    std::fmax(std::fmax(v0.x, v1.x), v2.x), std::fmax(std::fmax(v0.y, v1.y), v2.y), std::fmax(std::fmax(v0.z, v1.z), v2.z) }, 
                    (v0 + v1 + v2) / glm::vec3 { 3.0f }, /*t,*/ j, i });
        }
    }
    std::vector<int> i(prims.size());
    std::iota(i.begin(), i.end(), 0);

    //Recursive
    ConstructorHelper(prims, i, nodes, 0, -1, 0);

    //Iterative
    //int idx = 0;
    //std::deque<std::vector<int>> queue;
    //queue.push_back(i);

    //while (queue.size() != 0) 
    //{
    //    this->m_numLevels = std::max(this->m_numLevels, int(std::log2(idx - 1)));
    //    std::vector<int> prim_ids = queue.front();
    //    int parentIdx = (idx - 1) / 2;
    //    queue.pop_front();

    //    BVHNode current;
    //    current.n_id = idx;
    //    current.box = calculateAABB(prims, prim_ids);

    //    if (this->m_numLevels == max_level || prim_ids.size() == 1) {
    //        current.isLeafNode = true;
    //        std::for_each(prim_ids.begin(), prim_ids.end(), [&](int i) {
    //            current.ids.push_back(prims[i].t_id);
    //            current.ids.push_back(prims[i].m_id);
    //        });
    //        this->m_numLeaves++;
    //        nodes.push_back(current);
    //        if (parentIdx != -0)
    //            nodes[parentIdx].ids.push_back(idx);
    //    } else {
    //        current.isLeafNode = false;
    //        // Sort by centroid
    //        std::sort(prim_ids.begin(), prim_ids.end(), [&](int i, int j) {
    //            return prims[i].centr[this->m_numLevels % 3] < prims[j].centr[this->m_numLevels % 3];
    //        });
    //        nodes.push_back(current);
    //        if (parentIdx != -0)
    //            nodes[parentIdx].ids.push_back(idx);

    //        queue.push_back({ prim_ids.begin(), prim_ids.begin() + prim_ids.size() / 2 });
    //        queue.push_back({ prim_ids.begin() + prim_ids.size() / 2, prim_ids.end() });
    //    }
    //    idx++;
    //}
    this->nodes = nodes;
}

void BoundingVolumeHierarchy::ConstructorHelper(std::vector<Prim>& prims, std::vector<int> prim_ids,
    std::vector<BVHNode>& nodes, int currLevel, int parentIdx, int idx)
{
    this->m_numLevels = std::max(this->m_numLevels, currLevel);
    BVHNode current;
    current.n_id = idx;
    current.box = calculateAABB(prims, prim_ids);

    if (currLevel == max_level || prim_ids.size() == 1) {
        current.isLeafNode = true;
        std::for_each(prim_ids.begin(), prim_ids.end(), [&](int i) {
            current.ids.push_back(prims[i].t_id);
            current.ids.push_back(prims[i].m_id);
        });
        this->m_numLeaves++;
        nodes.push_back(current);
        if (parentIdx != -1) 
            nodes[parentIdx].ids.push_back(nodes.size() - 1);
    } else {
        current.isLeafNode = false;

        // Sort by centroid
        std::sort(prim_ids.begin(), prim_ids.end(), [&](int i, int j) {
            return prims[i].centr[currLevel % 3] < prims[j].centr[currLevel % 3];
        });

        nodes.push_back(current);
        if(parentIdx != -1) 
            nodes[parentIdx].ids.push_back(nodes.size() - 1);

        // Recursive call
        ConstructorHelper(prims, { prim_ids.begin(), prim_ids.begin() + prim_ids.size() / 2 }, nodes, currLevel + 1, nodes.size() - 1, idx + 1);
        ConstructorHelper(prims, { prim_ids.begin() + prim_ids.size() / 2, prim_ids.end() }, nodes, currLevel + 1, nodes.size() - 1, idx + 2);
    }
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return this->m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return this->m_numLeaves;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}
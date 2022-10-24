#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include "config.h"
#include <glm/glm.hpp>
#include <numeric>
#include <deque>

// Helper method for calcualating the new bounding volume based on prims and the ids of prims to calculate for
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
    this->max_level = 10; // Hardcoded for now, add slider later
    this->m_numLeaves = 0;
    this->m_numLevels = 0;
    std::vector<BVHNode> nodes;

    // Create prim vector -> store for each primitive ->
    // triangle: *opt{vertices}, min_vector, max_vector, centroid, *opt{indexes of vertices}, index of triangle, index of mesh
    // sphere: *opt{center x3}, min_vector, max_vector, center, *opt{zero * 3}, index of sphere, -1
    // *opt{...} -> commented out, uncomment if needed for other features
    std::vector<Prim> prims;
    // Add meshes
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
    // Add spheres
    for (int i = 0; i < pScene->spheres.size(); i++) 
    {
        Sphere s = pScene->spheres[i];
        prims.push_back(Prim { /*std::vector { s.center },*/
            { s.center - s.radius }, { s.center + s.radius }, s.center, /*glm::vec3 { 0 },*/ i, -1 });
    }
    // Vector with ints from 0 to prim.size()
    std::vector<int> i(prims.size());
    std::iota(i.begin(), i.end(), 0);

    //Recursive
    ConstructorHelper(prims, i, nodes, 0, -1, 0);

    //Iterative
    //int idx = 0, level = 0;
    //std::deque<std::vector<int>> queue;
    //queue.push_back(i);

    //while (queue.size() != 0) 
    //{
    //    level = int(std::log2(idx + 1));
    //    std::vector<int> prim_ids = queue.front();
    //    int parentIdx = (idx - 1) / 2;
    //    queue.pop_front();

    //    BVHNode current;
    //    current.n_id = idx;
    //    current.level = level;
    //    current.box = calculateAABB(prims, prim_ids);

    //    if (level == max_level || prim_ids.size() == 1) {
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
    //            return prims[i].centr[level % 3] < prims[j].centr[level % 3];
    //        });
    //        nodes.push_back(current);
    //        if (parentIdx != -0)
    //            nodes[parentIdx].ids.push_back(idx);

    //        queue.push_back({ prim_ids.begin(), prim_ids.begin() + prim_ids.size() / 2 });
    //        queue.push_back({ prim_ids.begin() + prim_ids.size() / 2, prim_ids.end() });
    //    }
    //    idx++;
    //}
    // this->m_numLevels = nodes[nodes.size() - 1].level + 1;
    this->nodes = nodes;
}

// Constructor helper for recursion
void BoundingVolumeHierarchy::ConstructorHelper(std::vector<Prim>& prims, std::vector<int> prim_ids, 
    std::vector<BVHNode>& nodes, int currLevel, int parentIdx, int idx)
{
    // Update number of levels
    this->m_numLevels = std::max(this->m_numLevels, currLevel + 1);
    BVHNode current;
    //current.n_id = idx;
    current.level = currLevel;
    current.box = calculateAABB(prims, prim_ids);

    // Handle leaf node
    if (currLevel == this->max_level || prim_ids.size() == 1) {
        current.isLeafNode = true;
        std::for_each(prim_ids.begin(), prim_ids.end(), [&](int i) {
            current.ids.push_back(prims[i].t_id);
            current.ids.push_back(prims[i].m_id);
        });
        this->m_numLeaves++;
        nodes.push_back(current);
        if (parentIdx != -1) 
            nodes[parentIdx].ids.push_back(idx);
    // Handle internal node
    } else {
        current.isLeafNode = false;

        // Sort by centroid/sphere center
        std::sort(prim_ids.begin(), prim_ids.end(), [&](int i, int j) {
            return prims[i].centr[currLevel % 3] < prims[j].centr[currLevel % 3];
        });

        nodes.push_back(current);
        if(parentIdx != -1) 
            nodes[parentIdx].ids.push_back(idx);

        // Recursive call
        ConstructorHelper(prims, { prim_ids.begin(), prim_ids.begin() + prim_ids.size() / 2 }, nodes, currLevel + 1, idx, nodes.size());
        ConstructorHelper(prims, { prim_ids.begin() + prim_ids.size() / 2, prim_ids.end() }, nodes, currLevel + 1, idx, nodes.size());
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

int BoundingVolumeHierarchy::maxLevel() const
{
    return this->max_level;
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
    std::for_each(this->nodes.begin(), this->nodes.end(), [&](BVHNode n) {
        if (n.level == level)
            drawAABB(n.box, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 1.0f);
    });
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
    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
    int acc = 0;
    for (BVHNode curr : this->nodes) 
    {
        if (curr.isLeafNode) acc++;
        if (acc == leafIdx && leafIdx > 0 && leafIdx < this->nodes.size()) {
            drawAABB(curr.box, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 1.0f);
            for (int j = 0; j < curr.ids.size(); j += 2) 
            {
                int mesh_id = curr.ids[j + 1];
                // sphere
                if (mesh_id == -1) {
                    Sphere& sphere = this->m_pScene->spheres[curr.ids[j]];
                    drawSphere(sphere);
                } else {
                    Mesh& mesh = this->m_pScene->meshes[mesh_id];
                    glm::uvec3 t = mesh.triangles[curr.ids[j]];
                    drawTriangle(mesh.vertices[t.x], mesh.vertices[t.y], mesh.vertices[t.z]);
                }
            }
            return;
        }
    }
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
        bool foundIntersection = false;
        Vertex v0Debug;
        Vertex v1Debug;
        Vertex v2Debug;
        float smallestT = ray.t;
        
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    if (features.enableTextureMapping) {
                        if (mesh.material.kdTexture != nullptr) {
                            glm::vec2 texCoords = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                            glm::vec3 tex = acquireTexel(*mesh.material.kdTexture, texCoords, features);
                            hitInfo.material = mesh.material;
                            hitInfo.material.kd = tex;
                            hit = true;
                        }
                    } else {
                        hitInfo.material = mesh.material;
                        hit = true;
                    }

                    if (features.enableNormalInterp){
                        if (smallestT > ray.t){
                            // update all debug rays and toggle the foundIntersection boolean
                            foundIntersection = true;
                            smallestT = ray.t;
                            v0Debug = v0;
                            v1Debug = v1;
                            v2Debug = v2;
                        }
                    }
                }
            }
        }

        if (features.enableNormalInterp && foundIntersection) {
            glm::vec3 point = ray.origin + ray.direction * ray.t;
            float length = 0.5f;

            // draw the rays of each vertex of the triangle
            drawRay(Ray { v0Debug.position, v0Debug.normal, length });
            drawRay(Ray { v1Debug.position, v1Debug.normal, length });
            drawRay(Ray { v2Debug.position, v2Debug.normal, length });

            // get the interpolated normal
            glm::vec3 color = glm::vec3 { 0.0f, 1.0f, 0.0f };
            glm::vec3 barycentric = computeBarycentricCoord(v0Debug.position, v1Debug.position, v2Debug.position, point);
            glm::vec3 interpolatedNormal = interpolateNormal(v0Debug.normal, v1Debug.normal, v2Debug.normal, barycentric);

            // draw the interpolated ray
            drawRay(Ray {point, interpolatedNormal, length}, color);
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
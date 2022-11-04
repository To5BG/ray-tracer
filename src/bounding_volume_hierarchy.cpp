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
#include <stack>
#include <iostream>

extern bool intersectedButNotTraversed;

int extr_max_level = 32;
int extr_sah_bins = 64;
bool extr_debugSAH = false;

// Helper method for calculating the new bounding volume based on prims and the ids of prims to calculate for
AxisAlignedBox calculateAABB(std::vector<Prim>& prims, std::vector<int>& prim_ids, int start = 0, int end = -1) 
{
    glm::vec3 min = glm::vec3 { std::numeric_limits<float>::max() };
    glm::vec3 max = glm::vec3 { -std::numeric_limits<float>::max() };
    for (int i = start; i < (end == -1 ? prim_ids.size() : end); i++) {
        Prim p = prims[prim_ids[i]];
        min = { std::fmin(min.x, p.min.x), std::fmin(min.y, p.min.y), std::fmin(min.z, p.min.z) };
        max = { std::fmax(max.x, p.max.x), std::fmax(max.y, p.max.y), std::fmax(max.z, p.max.z) };
    }
    return { min, max };
}

// Helper method for calculating surface area of an AABB, used for the surface-area heuristic
float surfaceArea(AxisAlignedBox& a)
{
    glm::vec3 b = a.upper - a.lower;
    b = { b.y, b.z, b.x };
    return glm::dot(a.upper - a.lower, b);
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene)
{
    // Start clock for benchmarking
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
    // Initial values
    this->max_level = extr_max_level;
    this->m_numLeaves = 0;
    this->m_numLevels = 0;
    std::vector<BVHNode> nodes;
    // Save bool for visual debug
    this->feat = features;

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
    ConstructorHelper(prims, i, nodes, 0, -1, 0, this->feat.extra.enableBvhSahBinning);
    this->nodes = nodes;
    const auto end = clock::now();
    // check if SAH or regular BVH was generated
    if (this->feat.extra.enableBvhSahBinning) {
        std::cout << "Time to create BVH with SAH and ";
        // if enough primitives to bin, consider that
        if (extr_sah_bins < prims.size())
            std::cout << std::setprecision(4) << extr_sah_bins << "-way binning: ";
        // otherwise sah with every centroid checked
        else
            std::cout << "no binning: ";
        std::cout << std::setprecision(-1) << std::chrono::duration<float, std::milli>(end - start).count() << " milliseconds " << std::endl;
    } else
        std::cout << "Time to create basic BVH: " << std::chrono::duration<float, std::milli>(end - start).count() << " milliseconds " << std::endl;
}

// Constructor helper for recursion
void BoundingVolumeHierarchy::ConstructorHelper(std::vector<Prim>& prims, std::vector<int> prim_ids, 
    std::vector<BVHNode>& nodes, int currLevel, int parentIdx, int idx, bool enabledSAHBinning)
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
        // Store child idx to parent
        if (parentIdx != -1) 
            nodes[parentIdx].ids.push_back(idx);
    // Handle internal node
    } else {
        current.isLeafNode = false;

        int splitPoint = 0;
        if (enabledSAHBinning) 
        {
            // Sort by all axes and find minimal SAH cost
            // SAH cost defined as (SA(left) * prim_num(left) + SA(right) * prim_num(right)) / SA(full)
            // Explanation given in the report for the above formula
            float minCost = std::numeric_limits<float>::max();
            // Micro-optimization by calculating current box volume only once
            float cur_surfaceAreaRec = 1.0f / surfaceArea(current.box);
            int axis = -1;
            for (int a = 0; a < 3; a++) 
            {
                // Sort by current axis
                std::sort(prim_ids.begin(), prim_ids.end(), [&](int i, int j) {
                    return prims[i].centr[a] < prims[j].centr[a];
                });
                // if more or bins than primitives or equal > check every centroid split
                if (extr_sah_bins >= prim_ids.size()) {
                    // Find best split by the aforementioned SAH criteria
                    for (int i = 1; i < prim_ids.size(); i++) {
                        AxisAlignedBox left = calculateAABB(prims, prim_ids, 0, i);
                        AxisAlignedBox right = calculateAABB(prims, prim_ids, i);
                        float currCost = (surfaceArea(left) * i + surfaceArea(right) * (prim_ids.size() - i)) * cur_surfaceAreaRec;
                        if (currCost < minCost) {
                            // Store min cost, index of split point, and best axis
                            minCost = currCost;
                            axis = a;
                            splitPoint = i;
                        }
                    }
                } else {
                    // Get centroid range, then split on even bins
                    float centroidrange = prims[prim_ids[prim_ids.size() - 1]].centr[a] - prims[prim_ids[0]].centr[a];
                    // Precompute distance between two bin boundaries
                    float dist = centroidrange / extr_sah_bins;
                    // Start from first primitive, for each bin boundary increment to encapsulate all prim ids left of the boundary
                    int countLeft = 1;
                    for (int i = 0; i < extr_sah_bins; i++) 
                    {
                        float currSplit = prims[prim_ids[0]].centr[a] + i * dist; 
                        // increment countLeft until prim with centroid proj larger than current bin boundary is found
                        while (prims[prim_ids[countLeft]].centr[a] < currSplit)
                            countLeft++;
                        AxisAlignedBox left = calculateAABB(prims, prim_ids, 0, countLeft);
                        AxisAlignedBox right = calculateAABB(prims, prim_ids, countLeft);
                        float currCost = (surfaceArea(left) * countLeft + surfaceArea(right) * (prim_ids.size() - countLeft)) * cur_surfaceAreaRec;
                        if (currCost < minCost) {
                            // Store min cost, index of split point, and best axis
                            minCost = currCost;
                            axis = a;
                            splitPoint = countLeft;
                        }
                    }
                }
            }
            // Sort by best axis for further recursive calls
            std::sort(prim_ids.begin(), prim_ids.end(), [&](int i, int j) {
                return prims[i].centr[axis] < prims[j].centr[axis];
            });
        } else {
            // Sort by centroid/sphere center
            std::sort(prim_ids.begin(), prim_ids.end(), [&](int i, int j) {
                return prims[i].centr[currLevel % 3] < prims[j].centr[currLevel % 3];
            });
            // Take median as split point
            splitPoint = prim_ids.size() / 2;
        }
        // Store child idx to parent
        nodes.push_back(current);
        if(parentIdx != -1) 
            nodes[parentIdx].ids.push_back(idx);

        // Recursive call
        ConstructorHelper(prims, { prim_ids.begin(), prim_ids.begin() + splitPoint }, nodes, currLevel + 1, idx, nodes.size(), enabledSAHBinning);
        ConstructorHelper(prims, { prim_ids.begin() + splitPoint, prim_ids.end() }, nodes, currLevel + 1, idx, nodes.size(), enabledSAHBinning);
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
    if (extr_debugSAH) 
    {
        Features newFeat = feat;
        newFeat.extra.enableBvhSahBinning = false;
        BoundingVolumeHierarchy bvh_two = { m_pScene, newFeat };
        std::for_each(bvh_two.nodes.begin(), bvh_two.nodes.end(), [&](BVHNode n) {
            if (n.level == level)
                drawAABB(n.box, DrawMode::Wireframe, glm::vec3(1.0f, 0.0f, 0.0f), 0.5f);
        });
    }
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


bool BoundingVolumeHierarchy::traversal(HitInfo& hitInfo, Ray& ray, const Features& features, BVHNode& node, bool& hit, float& absoluteT, int& finalMesh, int& finalTriangle, bool& sphereInt) const
{
    float oldT = ray.t; // ray distance
    if (node.level == 0 && !intersectRayWithShape(node.box, ray)) {
        ray.t = oldT;
        return false;
    } else {
        ray.t = oldT;
    }

    if (!features.enableRecursive && !features.extra.enableTransparency) {
        ray.t = std::numeric_limits<float>::max();
        if (intersectRayWithShape(node.box, ray)) {
            if (ray.t > absoluteT) {
                ray.t = oldT;
                if (intersectedButNotTraversed) {
                    drawAABB(node.box, DrawMode::Wireframe, glm::vec3(0.5f, 0.0f, 0.7f), 1.0f); // purple
                }
                return hit;
            } else {
                drawAABB(node.box, DrawMode::Wireframe, glm::vec3(0.0f, 0.7f, 0.0f), 1.0f); // green
                ray.t = oldT;
            }
        }
    } else {
        drawAABB(node.box, DrawMode::Wireframe, glm::vec3(0.0f, 0.7f, 0.0f), 1.0f); // green
    }

    if (node.isLeafNode) { // If leaf
        bool foundIntersection = false;
        Vertex v0Debug;
        Vertex v1Debug;
        Vertex v2Debug;
        float smallestT = ray.t;
        int i = 0;
        while (i < node.ids.size()) { // For each triangle mesh pair in ids
            int triangleID = node.ids[i]; // Get triangle ID
            int meshID = node.ids[i + 1]; // Get mesh ID
            if (meshID == -1) { 
                 Sphere& sphere = m_pScene->spheres[triangleID];
                 if (intersectRayWithShape(sphere, ray, hitInfo)) {
                     hitInfo.material = sphere.material;
                     hit = true;
                     sphereInt = true;
                }
            } else {
                Mesh mesh = m_pScene->meshes[meshID]; // Get mesh
                glm::uvec3 triangle = mesh.triangles[triangleID]; // Get triangle
                const auto v0 = mesh.vertices[triangle[0]];
                const auto v1 = mesh.vertices[triangle[1]];
                const auto v2 = mesh.vertices[triangle[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    if (ray.t < absoluteT) {
                        absoluteT = ray.t;
                    }
                    finalMesh = meshID;
                    finalTriangle = triangleID;

                    if (features.enableTextureMapping) {
                        if (mesh.material.kdTexture != nullptr) {
                            glm::vec2 texCoords = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                            glm::vec3 tex = acquireTexel(*mesh.material.kdTexture, texCoords, features);
                            hitInfo.material = mesh.material;
                            hitInfo.material.kd = tex;
                            hit = true;
                        } else {
                            hitInfo.material = mesh.material;
                            hit = true;
                        }
                    } else {

                        hitInfo.material = mesh.material;
                        hit = true;
                    }

                    if (features.enableNormalInterp) {

                        if (smallestT > ray.t) {
                            // update all debug rays and toggle the foundIntersection boolean
                            foundIntersection = true;
                            smallestT = ray.t;
                            v0Debug = v0;
                            v1Debug = v1;
                            v2Debug = v2;
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
                    hitInfo.normal = interpolatedNormal;
                    // draw the interpolated ray
                    drawRay(Ray { point, interpolatedNormal, length }, color);
                }
            }
            
            i += 2; // Go to next pair
        }

        return hit;
    } else // If internal
    {
        int left = node.ids[0];
        int right = node.ids[1];
        BVHNode leftNode = nodes[left];
        BVHNode rightNode = nodes[right];
        bool intersectsLeft = false;
        bool intersectsRight = false;
        float leftT; // Used to decide which node to push on stack first; node with closest t gets pushed on stack first
        float rightT;
        ray.t = std::numeric_limits<float>::max();
        if (intersectRayWithShape(leftNode.box, ray)) { // If left box is intersected, add to stack
            leftT = ray.t;
            intersectsLeft = true;
            ray.t = oldT;
        }
        ray.t = std::numeric_limits<float>::max();
        if (intersectRayWithShape(rightNode.box, ray)) { // If right box is intersected, add to stack
            rightT = ray.t;
            intersectsRight = true;
            ray.t = oldT;
        }
        ray.t = oldT;

        if (intersectsLeft && intersectsRight) { // If both left and right node intersect, check which is closest. If needed, swap the two
            if (leftT < rightT) {
                traversal(hitInfo, ray, features, leftNode, hit, absoluteT, finalMesh, finalTriangle, sphereInt);
                traversal(hitInfo, ray, features, rightNode, hit, absoluteT, finalMesh, finalTriangle, sphereInt);
                return hit;
            } else {
                traversal(hitInfo, ray, features, rightNode, hit, absoluteT, finalMesh, finalTriangle, sphereInt);
                traversal(hitInfo, ray, features, leftNode, hit, absoluteT, finalMesh, finalTriangle, sphereInt);
                return hit;
            }
        } else {
            if (!intersectsLeft && intersectsRight) {
                traversal(hitInfo, ray, features, rightNode, hit, absoluteT, finalMesh, finalTriangle, sphereInt);
                return hit;
            } else if (intersectsLeft && !intersectsRight) {
                traversal(hitInfo, ray, features, leftNode, hit, absoluteT, finalMesh, finalTriangle, sphereInt);
                return hit;
            }
        }

        return hit;
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
                        } else {
                            hitInfo.material = mesh.material;
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
            hitInfo.normal = interpolatedNormal;
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

        BVHNode root = nodes[0];
        float absoluteT = std::numeric_limits<float>::max();
        int finalMesh;
        int finalTriangle;
        bool hit = false;
        bool sphereInt = false;
        traversal(hitInfo, ray, features, root, hit, absoluteT, finalMesh, finalTriangle, sphereInt);
        if (hit && !sphereInt) {
            Mesh& mesh = m_pScene->meshes[finalMesh];
            glm::uvec3 t = mesh.triangles[finalTriangle];
            drawTriangle(mesh.vertices[t.x], mesh.vertices[t.y], mesh.vertices[t.z]);
        }
        return hit;
    }
}
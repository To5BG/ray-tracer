#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;

struct BVHNode {
    bool isLeafNode;
    int n_id;
    // isLeafNode ? (triangles in AABB) : (2 child nodes);
    // triangles are given as (triangle_id, mesh_id) pairs
    std::vector<int> ids;
    AxisAlignedBox box;
};

struct Prim {
    // positions of vertices
    std::vector<glm::vec3> vs;
    // centroid pos
    glm::vec3 centr;
    // ids of vertices
    glm::uvec3 v_id;
    // id of triangle
    int t_id;
    // id of linked mesh
    int m_id;
};

AxisAlignedBox calculateAABB(std::vector<Prim>& prims, std::vector<int>& prim_ids);

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    // constructor helper for recursion
    void ConstructorHelper(std::vector<Prim>& prims, std::vector<int> prim_ids, std::vector<BVHNode>& nodes, int currLevel, int parentIdx, int idx);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;


private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<BVHNode> nodes;
    int max_level = 24;
};
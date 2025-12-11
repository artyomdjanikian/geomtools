#pragma once

#include <eigen3/Eigen/Core>
#include <vector>
#include "tools3d.h"

typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

// axis aligned bounding box tree, working on OpenMesh polygon meshes
// supporting ray tracing, nearest point search, computing signed distance to a point

class AABBTree
{
    public:

    void Build(MyMesh &mesh, uint32_t maxDepth);
    void Print() const;

    struct sRayX {
        double dist;
        int facetIndex;
        Eigen::Vector3d pnt;
    };

    bool TraceRay(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayDir, std::vector<sRayX> &rayXVec) const;
    sRayX FindNearestPoint(Eigen::Vector3d queryPnt) const;
    double ComputeSignedDistance(MyMesh &mesh, Eigen::Vector3d pos, AABBTree::sRayX &result) const;

    Bounds3d GetBounds() const;

    // TODO : std::vector<int> FindNearbyFacetCandidates(double dist);
    // TODO : std::vector<std::pair<int, int>> FindCollisionCandidates();

private:
    class AABBTreeNode
    {
        AABBTreeNode(size_t _iLow, size_t _iHigh) : iLow(_iLow), iHigh(_iHigh) {pChildren[0] = pChildren[1] = nullptr;}
 
        friend class AABBTree;

        private:
        AABBTreeNode * pChildren[2];
        Bounds3d bounds;

        size_t iLow;
        size_t iHigh;
    };

private:

    struct sFacetRecord {
        std::array<Eigen::Vector3d, 3> pnts = { };
        Eigen::Vector3d centroid = {0.0, 0.0, 0.0};
        int facetIndex = -1;
    };
    
    Bounds3d computeBounds(size_t iLow, size_t iHigh);
    void subdivide(AABBTree::AABBTreeNode * pNode, uint32_t maxDepth);
    void print(AABBTreeNode * pNode, int depth) const;

    std::vector<sFacetRecord> m_facets;
    AABBTreeNode * m_pRootNode;
};

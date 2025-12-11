 #include <iostream>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <eigen3/Eigen/Core>
#include "AABBTree.h"

typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

Bounds3d AABBTree::computeBounds(size_t iLow, size_t iHigh)
{
    Bounds3d bounds; 

    for(size_t ii = iLow; ii < iHigh; ii++) {
        bounds += m_facets[ii].pnts[0];
        bounds += m_facets[ii].pnts[1];
        bounds += m_facets[ii].pnts[2];
    }

    return bounds;
}

void AABBTree::Build(MyMesh &mesh, uint32_t maxDepth) 
{
    for(auto facetHandle : mesh.faces()) {
        sFacetRecord rec;

        for(size_t iP = 0; iP < 3; iP++)
            rec.pnts[iP] = {0.0, 0.0, 0.0};

        size_t iV = 0;
        for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
            // only triangles allowed in AABBTree
            if(iV == 3)
                throw std::logic_error("only triangles in AABBTree");

            auto pnt = mesh.point(*fv_it);

            rec.pnts[iV++] = {pnt[0], pnt[1], pnt[2]};
        }

        rec.centroid = (rec.pnts[0] + rec.pnts[1] + rec.pnts[2])/3.0;
        rec.facetIndex = facetHandle.idx();

        m_facets.emplace_back(std::move(rec));
    }

    m_pRootNode = new AABBTreeNode(0, m_facets.size());

    subdivide(m_pRootNode, maxDepth);
}

void AABBTree::Print() const
{
    print(m_pRootNode, 0);
}

void AABBTree::print(AABBTree::AABBTreeNode * pNode, int depth) const
{
    for(int iP = 0; iP < depth; iP++)
        printf(" ");

    int nChildren = 0;
    if(pNode->pChildren[0] != nullptr)
        nChildren++;
    if(pNode->pChildren[1] != nullptr)
        nChildren++;

    double diag = (pNode->bounds.maxBound-pNode->bounds.minBound).norm();

    if(nChildren == 0)
        printf("diag %.6f, %lu facets\n", diag, pNode->iHigh-pNode->iLow);
    else
        printf("diag %.6f, %d children\n", diag, nChildren);

    if(pNode->pChildren[0] != nullptr)
        print(pNode->pChildren[0], depth+1);

    if(pNode->pChildren[1] != nullptr)
        print(pNode->pChildren[1], depth+1);
}

bool AABBTree::TraceRay(Eigen::Vector3d rayOrigin, Eigen::Vector3d rayDir, std::vector<AABBTree::sRayX> &rayXVec) const
{
    AABBTree::sRayX result;
    result.dist = std::numeric_limits<double>::max();
    result.facetIndex = -1;
 
    std::vector<std::pair<AABBTreeNode *, double>> stack;

    stack.push_back({m_pRootNode, -1.0});

    Line3d line(rayOrigin, rayOrigin + rayDir);

    int nNodesTested = 0;
    int nFacetsTested = 0;

    while(stack.size()) {
        auto nodeRec = stack.back();
        auto pNode = nodeRec.first;
        auto nodeDist = nodeRec.second;
        stack.pop_back();

        nNodesTested++;

        // check if ray intersects bounds
        if(pNode == nullptr || result.dist < nodeDist)
            continue;
    
        if(pNode->pChildren[0] != nullptr || pNode->pChildren[1] != nullptr) {

            bool child1X = false;
            double child1Dist = -1.0;
            if(pNode->pChildren[0] != nullptr) {
                std::pair<double, double> params;

                if(pNode->pChildren[0]->bounds.IntersectLine(line, &params)) {
                    child1Dist = params.first;
                    if(params.first >= 0.0 || params.second >= 0.0)
                        child1X = true;
                }
            }

            bool child2X = false;
            double child2Dist = -1.0;
            if(pNode->pChildren[1] != nullptr) {
                std::pair<double, double> params;

                if(pNode->pChildren[1]->bounds.IntersectLine(line, &params)) {
                    child2Dist = params.first;
                    if(params.first >= 0.0 || params.second >= 0.0)
                        child2X = true;
                }
            }
    
            int iFirstChild = 0;
            int iSecondChild = 1;

            // since it's lifo, push the largest distance first
            if(child2Dist < child1Dist) {
                std::swap(iFirstChild, iSecondChild);
                std::swap(child1Dist, child2Dist);
                std::swap(child1X, child2X);
            }

            if(pNode->pChildren[iSecondChild] != nullptr && child2X) 
                stack.push_back({pNode->pChildren[iSecondChild], child2Dist});
            
            if(pNode->pChildren[iFirstChild] != nullptr && child1X)
                stack.push_back({pNode->pChildren[iFirstChild], child1Dist});
        }
        else {
            for(size_t iF = pNode->iLow; iF < pNode->iHigh; iF++) {
                Triangle3d triangle(m_facets[iF].pnts[0], m_facets[iF].pnts[1],m_facets[iF].pnts[2]);

                nFacetsTested++;

                std::pair<double, double> params;
                int nInt = triangle.IntersectLine(line, &params);

//                    if(nInt > 0)
//                        printf("   test facet %d : %d intersections (%.6f, %.6f)\n", m_facets[iF].facetIndex, nInt, params.first, params.second);
//                    else
//                        printf("   test facet %d\n", m_facets[iF].facetIndex);

                for(int iInt = 0; iInt < nInt; iInt++) {
                    double param = iInt == 0 ? params.first : params.second;
                    if(param > eps) {
                        if(param < result.dist) {
                            result.dist = param;
                            result.facetIndex = m_facets[iF].facetIndex;
                            result.pnt = line.Eval(param);
                        }
                    }
                }
            }
        }
    }

//    printf("%d nodes tested, %d facets tested\n", nNodesTested, nFacetsTested);

    if(result.facetIndex != -1) {
//        printf(" *** found %d\n", result.facetIndex);
        rayXVec.push_back(result);

        return true;
    }

    return false;
}

#if 0
class cStack {

public:
    void push_back(const std::pair<AABBTreeNode *, double> &pair) {
        m_stack[count++] = pair;
    }

    const std::pair<AABBTreeNode *, double> & back() const {
        return m_stack[count-1];
    }

    void pop_back() {
        count--;
    }

private:
    std::array<<std::pair<AABBTreeNode *, double>, 20> m_stack;
    size_t count = 0;
};
#endif

AABBTree::sRayX AABBTree::FindNearestPoint(Eigen::Vector3d queryPnt) const
{
    int nNodesTested = 0;
    int nFacetsTested = 0;

    AABBTree::sRayX result;
    result.dist = std::numeric_limits<double>::max();
    result.facetIndex = -1;

    // TODO : store a pair of tree node and bbox distance already computed
    std::vector<std::pair<AABBTreeNode *, double>> stack;

    stack.push_back({m_pRootNode, -1.0});

    while(stack.size()) {
        nNodesTested++;

        auto nodeRec = stack.back();
        auto pNode = nodeRec.first;
        auto nodeDist = nodeRec.second;
        stack.pop_back();

        if(pNode != nullptr && (result.dist == std::numeric_limits<double>::max() || nodeDist <= fabs(result.dist))) {
//            pNode->bounds.PointDistanceSq(queryPnt) <= fabs(result.dist)) {
               auto pChild1 = pNode->pChildren[0];
               auto pChild2 = pNode->pChildren[1];

               if(pChild1 != nullptr || pChild2 != nullptr) {
                    double d1Sq = pChild1 != nullptr ? pChild1->bounds.PointDistanceSq(queryPnt) : -1.0;
                    double d2Sq = pChild2 != nullptr ? pChild2->bounds.PointDistanceSq(queryPnt) : -1.0;

                    // since it's lifo, push bigger distance first
                    if(d2Sq > d1Sq) {
                        if(d2Sq <= fabs(result.dist))
                            stack.push_back({pChild2, d2Sq});
                        if(d1Sq <= fabs(result.dist))
                            stack.push_back({pChild1, d1Sq});
                    }
                    else {
                        if(d1Sq <= fabs(result.dist))
                            stack.push_back({pChild1, d1Sq});
                        if(d2Sq <= fabs(result.dist))
                            stack.push_back({pChild2, d2Sq});
                    }
                }
                else {
                    for(size_t iF = pNode->iLow; iF < pNode->iHigh; iF++) {

                        nFacetsTested++;
                        Triangle3d triangle(m_facets[iF].pnts[0], m_facets[iF].pnts[1],m_facets[iF].pnts[2]);

                        std::pair<Eigen::Vector3d, double> pair = triangle.PointDistanceSq(queryPnt);

//                        printf(" facet %d, dist %.6f\n", m_facets[iF].facetIndex, sqrt(pair.second));

                        if(pair.second < fabs(result.dist)) {

//                            printf(" * take %d (%.6f, %.6f, %.6f)\n", m_facets[iF].facetIndex, pair.first[0], pair.first[1], pair.first[2]);
                            result.dist = pair.second;
                            result.pnt = pair.first;
                            result.facetIndex = m_facets[iF].facetIndex;

                            double projOrient = triangle.GetNormal().dot(queryPnt-result.pnt);
                            if(projOrient <= 0.0)
                                result.dist *= -1.0;
                        }
                    }
                }
            }
    }

    // compute proper sign of distance

    if(result.dist > 0.0)
        result.dist = sqrt(result.dist);
    else
        result.dist = -sqrt(-result.dist);

    //    auto iF = result.facetIndex;
    //    Plane3d plane(m_facets[iF].pnts[0], m_facets[iF].pnts[1],m_facets[iF].pnts[2]);

    //    auto normal = plane.GetNormal();

    //    printf("query point (%.6f, %.6f, %.6f)\n", queryPnt[0], queryPnt[1], queryPnt[2]);
    //    printf("facet %d : (%.6f, %.6f, %.6f)\n", iF, normal[0], normal[1], normal[2]);
    //    printf("nearest point (%.6f, %.6f, %.6f)\n\n", result.pnt[0], result.pnt[1], result.pnt[2]);

    //    if(plane.SignedDistance(queryPnt) < 0.0)
    //        result.dist *= -1.0;

    //    printf("%d nodes tested, %d facets tested\n", nNodesTested, nFacetsTested);

    return result;
}

// TODO : make it a method of AABBTree
double AABBTree::ComputeSignedDistance(MyMesh &mesh, Eigen::Vector3d pos, AABBTree::sRayX &result) const
{
        // check if projPnt is on the edge/vertex
        auto face = MyMesh::FaceHandle(result.facetIndex);

        auto getTriangle = [&mesh](MyMesh::FaceHandle fh)
        {
            std::vector<Eigen::Vector3d> triPnts;

            for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it)
            {
                auto pnt = mesh.point(*fv_it);
                triPnts.push_back(toVec(pnt));
            }

            assert(triPnts.size() >= 3);
            return Triangle3d(triPnts[0], triPnts[1], triPnts[2]);
        };

        // iterate over vertices of face
        // check if point lies on vertex
        for (MyMesh::FaceHalfedgeIter fhe_it = mesh.fh_iter(face); fhe_it.is_valid(); ++fhe_it)
        {
            auto fromVertex = mesh.from_vertex_handle(*fhe_it);
            auto fromVertexPnt = mesh.point(fromVertex);

            if ((toVec(fromVertexPnt) - result.pnt).norm() < eps)
            {

                // circulate halfedges around fromVertex
                for (MyMesh::VertexIHalfedgeIter vih_it = mesh.vih_iter(fromVertex); vih_it.is_valid(); ++vih_it)
                {
                    auto adjacentFace = mesh.face_handle(*vih_it);
                    auto oppositeFace = mesh.opposite_face_handle(*vih_it);

                    // write a function to build Triangle3d from adjacentFace vertex points
                    // circulate around adjacentFace to find projection onto triangle plane

                    auto facetTriangle = getTriangle(adjacentFace);
                    auto oppTriangle = getTriangle(oppositeFace);

                    auto oppCentroid = oppTriangle.EvalBarycentric({0.33, 0.33, 0.34});

                    double signedDistance = facetTriangle.GetPlane().SignedDistance(pos);
                    double oppCentroidDist = facetTriangle.GetPlane().SignedDistance(oppCentroid);

                    ///            printf("dist %.6f, facet dist %.6f, opp centroid dist %.6f\n", dist, signedDistance, oppCentroidDist);

                    if ((signedDistance > 0.0 && oppCentroidDist < 0.0) ||
                        (signedDistance < 0.0 && oppCentroidDist > 0.0))
                    {
                        //                printf("take %.6f\n", signedDistance);
                        return signedDistance;
                        // characteristic edge projection
                    }
                }
            }
        }

        // iterate over edges
        // check if resuot.pnt lies on the edge

        for (MyMesh::FaceHalfedgeIter fhe_it = mesh.fh_iter(face); fhe_it.is_valid(); ++fhe_it)
        {
            auto fromVertex = mesh.from_vertex_handle(*fhe_it);
            auto toVertex = mesh.to_vertex_handle(*fhe_it);
            Line3d edgeLine(toVec(mesh.point(fromVertex)), toVec(mesh.point(toVertex)));

            // interior edge projection
            double projParam = edgeLine.ProjectParam(result.pnt);

            //          if(projParam >= 0 && projParam <= 1.0) {
            //            printf("projParam %.6f, dist %.6f\n", projParam, (edgeLine.Eval(projParam) - result.pnt).norm());
            //          }

            if (projParam >= 0.0 && projParam <= 1.0 && (edgeLine.Eval(projParam) - result.pnt).norm() < 1e-10)
            {

                auto hedgeHandle = *fhe_it;

                for (int iHedge = 0; iHedge < 2; iHedge++)
                {

                    auto adjacentFace = mesh.face_handle(hedgeHandle);
                    auto oppositeFace = mesh.opposite_face_handle(hedgeHandle);

                    auto facetTriangle = getTriangle(adjacentFace);
                    auto oppTriangle = getTriangle(oppositeFace);

                    auto oppCentroid = oppTriangle.EvalBarycentric({0.33, 0.33, 0.34});

                    double signedDistance = facetTriangle.GetPlane().SignedDistance(pos);
                    double oppCentroidDist = facetTriangle.GetPlane().SignedDistance(oppCentroid);

                    //              printf(" dist %.6f, facet %d dist %.6f, opp centroid %d dist %.6f\n", dist, adjacentFace.idx(), signedDistance, oppositeFace.idx(), oppCentroidDist);

                    if ((signedDistance > 0.0 && oppCentroidDist < 0.0) ||
                        (signedDistance < 0.0 && oppCentroidDist > 0.0))
                    {
                        //                printf("take %.6f\n", signedDistance);
                        return signedDistance;
                    }

                    hedgeHandle = mesh.opposite_halfedge_handle(hedgeHandle);
                }
            }
        }

//        auto facetTriangle = getTriangle(face);
//        double signedDistance = facetTriangle.GetPlane().SignedDistance(pos);

        return result.dist;
}

Bounds3d AABBTree::GetBounds() const
{
    Bounds3d bounds;

    if (m_pRootNode != nullptr)
        bounds = m_pRootNode->bounds;

    return bounds;
}

void AABBTree::subdivide(AABBTree::AABBTreeNode * pNode, uint32_t maxDepth)
{
    static int count = 0;

    auto iLow = pNode->iLow;
    auto iHigh = pNode->iHigh;

    pNode->bounds = computeBounds(iLow, iHigh);

    //    printf("%d, depth %d, %lu facets, %.6f span", count++, 10-maxDepth, iHigh-iLow, (pNode->bounds.maxBound-pNode->bounds.minBound).norm());

    // depth limit or triangle count limit
    if (maxDepth == 0 || iLow + 2 >= iHigh)
    {
        //        printf(", leaf\n");
        return;
    }

    //    printf("\n");

    auto getSplitCoord = [this, bounds = pNode->bounds](size_t iLow, size_t iHigh)
    {
        double edgeXSq = 0.0;
        double edgeYSq = 0.0;
        double edgeZSq = 0.0;

        for (size_t ii = iLow; ii < iHigh; ii++)
        {
            const auto &facet = m_facets[ii];

            auto e1 = facet.pnts[0] - facet.pnts[1];
            auto e2 = facet.pnts[1] - facet.pnts[2];
            auto e3 = facet.pnts[2] - facet.pnts[0];

            edgeXSq += e1[0] * e1[0] + e2[0] * e2[0] + e3[0] * e3[0];
            edgeYSq += e1[1] * e1[1] + e2[1] * e2[1] + e3[1] * e3[1];
            edgeZSq += e1[2] * e1[2] + e2[2] * e2[2] + e3[2] * e3[2];
        }

        edgeXSq /= bounds.maxBound[0] - bounds.minBound[0];
        edgeYSq /= bounds.maxBound[1] - bounds.minBound[1];
        edgeZSq /= bounds.maxBound[2] - bounds.minBound[2];

        // choose axis for ordering
        int iCoord = 0;

        if (edgeXSq <= edgeYSq && edgeXSq <= edgeZSq)
            iCoord = 0;
        else if (edgeYSq <= edgeXSq && edgeYSq <= edgeZSq)
            iCoord = 1;
        else if (edgeZSq <= edgeXSq && edgeZSq <= edgeYSq)
            iCoord = 2;

        return iCoord;
    };

    auto iCoord = getSplitCoord(iLow, iHigh);

    double midCoord = (pNode->bounds.maxBound[iCoord] + pNode->bounds.minBound[iCoord]) * 0.5;

    // each node reorders the portion of facets owned by subtree rooted in the node

    // TODO : pick pivot element that's closest to midCoord
    // TODO : quicksort partition followed by quicksort of only one half
    std::sort(m_facets.begin() + iLow, m_facets.begin() + iHigh,
                [this, iCoord](const auto &f1, const auto &f2)
                { return f1.centroid[iCoord] < f2.centroid[iCoord]; });

    // allocate 2 nodes, split facets in two arrays, put in 2 nodes
    size_t iMid = (iHigh + iLow) / 2;

    pNode->pChildren[0] = new AABBTreeNode(iLow, iMid);
    pNode->pChildren[1] = new AABBTreeNode(iMid, iHigh);

    subdivide(pNode->pChildren[0], maxDepth - 1);
    subdivide(pNode->pChildren[1], maxDepth - 1);
}

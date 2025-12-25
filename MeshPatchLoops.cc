#include "common.h"
#include "tools3d.h"
#include "AABBTree.h"
#include "MeshSampling.h"
#include "MeshTools.h"
#include "PointSet.h"
#include "cGrid3d.h"
#include "happly.h"

void PatchLoops(MyMesh &mesh, double loopStep)
{
    AABBTree tree;
    tree.Build(mesh, 20);
  
    int nVertices = 0;
    Bounds3d bounds;

    happly::PLYExport bPly;

    // compute loop border vertices, facets, and points

    std::vector<int> borderFacets;
    std::set<MyMesh::VertexHandle> borderVertices;

    std::vector<Eigen::Vector3d> borderPoints;

    for(auto hedge : mesh.halfedges()) {
      if(mesh.is_boundary(hedge)) {

        auto toVertex = mesh.to_vertex_handle(hedge);
        auto fromVertex = mesh.from_vertex_handle(hedge);

        borderVertices.insert(toVertex);
        borderVertices.insert(fromVertex);

        borderPoints.push_back(toVec(mesh.point(fromVertex)));

        bounds += toVec(mesh.point(toVertex));
        bounds += toVec(mesh.point(fromVertex));
        nVertices += 2;

        auto oppHedge = mesh.opposite_halfedge_handle(hedge);
        auto nextHedge = mesh.next_halfedge_handle(oppHedge);
        auto oppVertex = mesh.to_vertex_handle(nextHedge);

        borderFacets.push_back(mesh.opposite_face_handle(hedge).idx());

        bPly.AddVertex(toVec(mesh.point(toVertex)), 1.0, 0.0, 0.0);
        bPly.AddVertex(toVec(mesh.point(fromVertex)), 0.0, 0.0, 1.0);
        bPly.AddVertex(toVec(mesh.point(oppVertex)), 0.0, 1.0, 0.0);

        printf("border face %d\n", mesh.opposite_face_handle(hedge).idx());
      }
    }

    if(nVertices == 0)
      return;

    printf("%d boundary vertices\n", nVertices);

    // 3d grid of points in the extended bounds of border loop vertices
    Eigen::Vector3d minPnt = bounds.minBound - Eigen::Vector3d(5.0*loopStep, 5.0*loopStep, 5.0*loopStep);
    Eigen::Vector3d maxPnt = bounds.maxBound + Eigen::Vector3d(5.0*loopStep, 5.0*loopStep, 5.0*loopStep);

    Eigen::Vector3d diag = maxPnt-minPnt;

    int nPnts[3];
    for(int ii = 0; ii < 3; ii++ ) {
      nPnts[ii] = static_cast<int>(ceil(diag[ii]/loopStep));
    }

    printf("%d x %d x %d grid\n", nPnts[0], nPnts[1], nPnts[2]);

    bounds += minPnt;
    bounds += maxPnt;

    auto isOnBorderFacet = [](MyMesh &mesh, MyMesh::FaceHandle facetHandle, Eigen::Vector3d projPnt ) {
      bool isOnBorderFacet = false;

      for (MyMesh::FaceHalfedgeIter fhe_it = mesh.fh_begin(facetHandle); fhe_it != mesh.fh_end(facetHandle); ++fhe_it)
      {
        MyMesh::HalfedgeHandle heh = mesh.opposite_halfedge_handle(*fhe_it);

        if( mesh.is_boundary(heh) ) {
          auto from = mesh.from_vertex_handle(heh);
          auto to = mesh.to_vertex_handle(heh);

          Line3d lineSeg(toVec(mesh.point(from)), toVec(mesh.point(to)));

          double distSq = lineSeg.SegmentDistanceSq(projPnt);
          if(distSq <= eps)
            isOnBorderFacet = true;
        }
      }

      return isOnBorderFacet;
    };

    // use grid to find points "filling" the gap
    auto fillLoopPoints_interpolateNearestPointVectorField = [&isOnBorderFacet](MyMesh &mesh, const AABBTree &tree, std::set<MyMesh::VertexHandle> &borderVertices, Bounds3d & bounds, double loopStep)
    {
      std::vector<Eigen::Vector3d> allpoints;

      happly::PLYExport cPly;

      cGrid3d grid;
      grid.init(bounds, loopStep);

      for (size_t ii = 0; ii < grid.GetNi(); ii++)
        for (size_t jj = 0; jj < grid.GetNj(); jj++)
          for (size_t kk = 0; kk < grid.GetNk(); kk++)
          {
            auto pos = grid.GetPosition(ii, jj, kk);

            auto np = tree.FindNearestPoint(pos);

            if (np.facetIndex != -1)
            {
              // TODO : check normals
              Eigen::Vector3d normal = toVec(mesh.normal(MyMesh::FaceHandle(np.facetIndex)));
              Eigen::Vector3d vec = np.pnt - pos;

              grid.GetVector(ii, jj, kk) = vec;

              vec.stableNormalize();
              double dotProd = vec.dot(normal);

              if (fabs(dotProd) > 0.95) // TODO : if projects not on boundary edges
                grid.GetPinned(ii, jj, kk) = true;
            }
          }

      grid.Smooth(10000);

      // TODO : for(const auto &gridPoint : grid.Points();
      happly::PLYExport sPly;

      for (size_t ii = 0; ii < grid.GetNi(); ii++)
        for (size_t jj = 0; jj < grid.GetNj(); jj++)
          for (size_t kk = 0; kk < grid.GetNk(); kk++)
          {

            Eigen::Vector3d pnt1 = grid.GetPosition(ii, jj, kk);
            Eigen::Vector3d vec1 = grid.GetVector(ii, jj, kk);

            if (!grid.GetPinned(ii, jj, kk))
            {
              Eigen::Vector3d interSurfacePoint = pnt1+vec1;
              // TODO : filter out allpoints that don't project onto border facets

              auto np = tree.FindNearestPoint(interSurfacePoint);

              MyMesh::FaceHandle facetHandle(np.facetIndex);

              if(isOnBorderFacet(mesh, facetHandle, np.pnt)) {
                sPly.AddVertex(interSurfacePoint, 1.0, 0.1, 0.2);
                allpoints.push_back(interSurfacePoint);
              }
            }

            // TODO : for(const auto &neighborPoint : grid.Neighbors(gridPoint))
            for (int iNeighbor = 0; iNeighbor < 6; iNeighbor++)
            {

              size_t iN = ii;
              size_t jN = jj;
              size_t kN = kk;

              switch (iNeighbor)
              {
              case 0:
                iN--;
                break;
              case 1:
                iN++;
                break;
              case 2:
                jN--;
                break;
              case 3:
                jN++;
                break;
              case 4:
                kN--;
                break;
              case 5:
                kN++;
                break;
              default:
                printf("!!! wrong neighbor index !!!\n");
              };

              if (grid.IsValid(iN, jN, kN) && !grid.GetPinned(ii, jj, kk) && !grid.GetPinned(iN, jN, kN))
              {
                Eigen::Vector3d vec2 = grid.GetVector(iN, jN, kN);
                Eigen::Vector3d pnt2 = grid.GetPosition(iN, jN, kN);

                // find the parameter value that minimizes the squared length of linear interpolation of vec1 and vec2

                double AA = vec1.dot(vec1);
                double BB = vec2.dot(vec2);
                double AB = vec1.dot(vec2);

                double denom = AA + BB - 2.0 * AB;

                if (denom > eps)
                {
                  double numer = BB - AB;

                  double param = numer / denom;

                  if (param > 0.0 && param < 1.0)
                  {
                    Eigen::Vector3d midPnt = pnt1 * param + pnt2 * (1.0 - param);
                    Eigen::Vector3d midVec = vec1 * param + vec2 * (1.0 - param);

                    auto np = tree.FindNearestPoint(midPnt);

                    MyMesh::FaceHandle facetHandle(np.facetIndex);

                    if(isOnBorderFacet(mesh, facetHandle, np.pnt))
                    {
                      //                    cPly.AddVertex(pnt1, 0.0, 1.0, 0.0);
                      //                    cPly.AddVertex(pnt2, 0.0, 1.0, 0.0);
                      cPly.AddVertex(midPnt + midVec, 0.6, 0.2, 0.7);

                      //                    happly::PLYExport dPly;

                      //                    dPly.AddVertex(pnt1, 0.0, 0.0, 1.0);
                      //                    dPly.AddVertex(pnt2, 0.0, 0.0, 1.0);

                      //                    dPly.AddVertex(pnt1+vec1, 0.0, 1.0, 0.0);
                      //                    dPly.AddVertex(pnt2+vec2, 0.0, 1.0, 0.0);
                      //                    dPly.AddVertex(midPnt+midVec, 1.0, 0.0, 1.0);

                      //                    dPly.WritePLY("debug.ply");

                      //                    printf("param = %.6f, lengths %.6f, %.6f, mid %.6f\n", param, vec1.norm(), vec2.norm(), midVec.norm());

                      //                      points.push_back(midPnt+midVec);
                      //                    getchar();
                      // TODO : debug this
                    }
                  }
                }
              }
            }
          }

      sPly.WritePLY("interpoints.ply");

      cPly.WritePLY("fillerpoints.ply");

      return allpoints;
    };

    auto fillLoopPoints_interpolateSignedDistance = [](MyMesh &mesh, const AABBTree &tree, std::set<MyMesh::VertexHandle> &borderVertices, Bounds3d &bounds, double loopStep)
    {
      std::vector<Eigen::Vector3d> allpoints;

      cGrid3d grid;
      grid.init(bounds, loopStep);

      for (size_t ii = 0; ii < grid.GetNi(); ii++)
        for (size_t jj = 0; jj < grid.GetNj(); jj++)
          for (size_t kk = 0; kk < grid.GetNk(); kk++)
          {

            auto pos = grid.GetPosition(ii, jj, kk);

            auto np = tree.FindNearestPoint(pos);

            if (np.facetIndex != -1)
            {
              // TODO : check normals
              Eigen::Vector3d normal = toVec(mesh.normal(MyMesh::FaceHandle(np.facetIndex)));
              Eigen::Vector3d vec = np.pnt - pos;

              vec.stableNormalize();
              double dotProd = vec.dot(normal);

              if (fabs(dotProd) > 0.95)
                grid.GetPinned(ii, jj, kk) = true;

              grid.GetVector(ii, jj, kk) = {np.dist, 0.0, 0.0};
            }
          }

      grid.Smooth(10000);

      double maxDist = 0.0;
      double minDist = 0.0;
      bool isInit = false;

      for (size_t ii = 0; ii < grid.GetNi(); ii++)
        for (size_t jj = 0; jj < grid.GetNj(); jj++)
          for (size_t kk = 0; kk < grid.GetNk(); kk++)
          {
            if (!grid.GetPinned(ii, jj, kk))
            {

              Eigen::Vector3d vector = grid.GetVector(ii, jj, kk);

              if (!isInit)
              {
                isInit = true;
                maxDist = vector[0];
                minDist = vector[0];
              }
              else
              {
                maxDist = std::max(maxDist, vector[0]);
                minDist = std::min(minDist, vector[0]);
              }
            }
          }

      printf("min dist %.6f, max dist %.6f\n", minDist, maxDist);

      happly::PLYExport sPly;

      for (size_t ii = 0; ii < grid.GetNi(); ii++)
        for (size_t jj = 0; jj < grid.GetNj(); jj++)
          for (size_t kk = 0; kk < grid.GetNk(); kk++)
          {
            if( true || !grid.GetPinned(ii, jj, kk)) {
              Eigen::Vector3d pnt = grid.GetPosition(ii, jj, kk);
              Eigen::Vector3d vector = grid.GetVector(ii, jj, kk);

              double dist = vector[0];
              double weight = (dist - minDist) / (maxDist - minDist);
              sPly.AddVertex(pnt, weight, 4.0 * weight * (1.0 - weight), (1.0 - weight));
            }
          }

      sPly.WritePLY("distpoints.ply");

      happly::PLYExport cPly;

      // TODO : marching cubes
      for (size_t ii = 0; ii < grid.GetNi(); ii++)
        for (size_t jj = 0; jj < grid.GetNj(); jj++)
          for (size_t kk = 0; kk < grid.GetNk(); kk++)
          {

            Eigen::Vector3d vec1 = grid.GetVector(ii, jj, kk);
            Eigen::Vector3d pnt1 = grid.GetPosition(ii, jj, kk);

            for (int iNeighbor = 0; iNeighbor < 6; iNeighbor++)
            {

              size_t iN = ii;
              size_t jN = jj;
              size_t kN = kk;

              switch (iNeighbor)
              {
              case 0:
                iN--;
                break;
              case 1:
                iN++;
                break;
              case 2:
                jN--;
                break;
              case 3:
                jN++;
                break;
              case 4:
                kN--;
                break;
              case 5:
                kN++;
                break;
              default:
                printf("!!! wrong neighbor index !!!\n");
              };

              if (grid.IsValid(iN, jN, kN))
              {
                Eigen::Vector3d vec2 = grid.GetVector(iN, jN, kN);
                Eigen::Vector3d pnt2 = grid.GetPosition(iN, jN, kN);

                double distProd = vec1[0]*vec2[0];

                if(distProd < eps && (!grid.GetPinned(ii, jj, kk) || !grid.GetPinned(iN, jN, kN)))
                {

                  double paramRange = vec1[0] - vec2[0];
                  assert(fabs(paramRange) > 0.0);

//                  Eigen::Vector3d marchingCubePoint = (vec2[0]*pnt1 - vec1[0]*pnt2) * 1.0/paramRange;

                  Eigen::Vector3d marchingCubePoint = (pnt1 + pnt2) * 0.5;

                  cPly.AddVertex(marchingCubePoint, 1.0, 0.1, 0.2);

                  allpoints.push_back(marchingCubePoint);
                }
              }
            }
          }
      

      // TODO : filter out strangepoints that don't project onto border facets
      cPly.WritePLY("interpoints.ply");

      getchar();

      return allpoints;
    };

    // TODO : save all points with signed distances
    //auto allpoints = fillLoopPoints_interpolateSignedDistance(mesh, tree, borderVertices, bounds, loopStep);
    auto allpoints = fillLoopPoints_interpolateNearestPointVectorField(mesh, tree, borderVertices, bounds, loopStep);

    happly::PLYExport exportPly;

    // cluster filler points
    printf("%lu raw points\n", allpoints.size());
    auto points = PointSetCluster(allpoints, 200, 500);

    printf("%lu border points, %lu new vertex points\n", borderPoints.size(), points.size());

// add border points to the point set
//    for(auto bpoint : borderPoints)
//      points.push_back(bpoint);


    // remove duplicate points
    PointSetRemoveDuplicates(points, 1e-5);

//    PointSetSmooth(points, 20, 100);

    // filter out interior points
    auto weights = PointSetGetExteriorWeights(points, 20);

    {
      std::vector<Eigen::Vector3d> newPoints;

      happly::PLYExport exportPly;

      for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

        if(weights[iPnt] < 0.2) {
          newPoints.push_back(points[iPnt]);

          const auto &point = points[iPnt];
          exportPly.AddVertex(point, weights[iPnt], 4.0*weights[iPnt]*(1.0-weights[iPnt]), (1.0-weights[iPnt]));
        }
      }

      exportPly.WritePLY("clusterpoints.ply");
      points = std::move(newPoints);
    }

    // for each pair of border points, compute the shortest path along the point set
    // mark points along the shortest path to use those in triangle creation

    // TODO : extract into a function
    if(borderPoints.size() > 3) {
      auto patchPoints = borderPoints;
      size_t nBorderPoints = patchPoints.size();

      printf("%lu border points, %lu exterior points\n", nBorderPoints, points.size());

      patchPoints.insert(patchPoints.end(), points.begin(), points.end());

      std::vector<bool> pathPoints(patchPoints.size(), false);

      // initialize distance matrix, use square of point distances
      Eigen::MatrixXd distanceMatrix(patchPoints.size(), patchPoints.size());

      for(size_t ii = 0; ii < patchPoints.size(); ii++) {

        double shortestEdge = -1.0;

        for(size_t jj = 0; jj < patchPoints.size(); jj++) {
          double thisDist = (patchPoints[ii]-patchPoints[jj]).norm();

          if(jj != ii) {
            if(shortestEdge < 0.0 || thisDist < shortestEdge)
              shortestEdge = thisDist;
          }

          // TODO : square distance
          distanceMatrix(ii, jj) = thisDist*thisDist;
        }

        printf("%.6f ", shortestEdge);
      }

      printf("closest neighbors\n");

      // shortest path from ii'th to jj'th border point

      std::vector<std::set<size_t>> dijedges;
      dijedges.resize(patchPoints.size());

//      nBorderPoints = patchPoints.size();
      for(size_t ii = 0; ii < nBorderPoints-1; ii++)
        for(size_t jj = ii+1; jj < nBorderPoints; jj++) {
      
          std::vector<std::size_t> shortestPath = findShortestPath(patchPoints, distanceMatrix, ii, jj);


          //  mark shortest path points 

          if(shortestPath.size() > 2) {
          pathPoints[ii] = true;
          pathPoints[jj] = true;

          double averLen = 0.0;

          for(size_t ii = 1; ii < shortestPath.size(); ii++) {
            averLen += (patchPoints[shortestPath[ii]]-patchPoints[shortestPath[ii-1]]).norm();
          }

          if(averLen > 0.0)
            averLen /= (shortestPath.size()-1);

          printf("%lu path points, aver. length %.6f\n", shortestPath.size(), averLen);

          for(size_t iStart = 1; iStart < shortestPath.size(); iStart++) {
            dijedges[shortestPath[iStart]].insert(shortestPath[iStart-1]);
            dijedges[shortestPath[iStart-1]].insert(shortestPath[iStart]);
          }

          for(auto iPath : shortestPath) {
          //  printf("   %lu", iPath); 
            pathPoints[iPath] = true;
          }
        }
          // printf("\n")
        }

//      MeshSamplingSave(patchPoints, dijedges, "dijxtraedges.ply");

          //      PointSetSmooth(patchPoints, pathPoints, 6, 3);

          // save points that have pathPoints flags set;

      points.clear();
      happly::PLYExport eply;

      bool takeAll = true;
      for(size_t ii = 0; ii < patchPoints.size(); ii++) {
        uint8_t blue = ii < nBorderPoints ? 255 : 0;
        uint8_t green = ii < nBorderPoints ? 0 : pathPoints[ii] ? 255 : 0;
        uint8_t red = ii < nBorderPoints ? 0 : pathPoints[ii] ? 0 : 255;

        eply.AddVertex(patchPoints[ii], red, green, blue);

        if(ii >= nBorderPoints && (pathPoints[ii]||takeAll))
          points.push_back(patchPoints[ii]);
      }

      eply.WritePLY("patchpathpoints.ply");

      printf("%lu points left\n", points.size());
    }

    // TODO : Taubin smooth point cloud while pinning path points 

    std::vector<MyMesh::FaceHandle> addedTriangles;

    auto lookupHedge = [](const MyMesh& mesh, MyMesh::VertexHandle vh_from, MyMesh::VertexHandle vh_to)
    {
      MyMesh::HalfedgeHandle hedgeHandle;

      // Check if both vertex handles are valid
      if (vh_from.is_valid() && vh_to.is_valid()) {
        // Iterate over the outgoing halfedges of vh_from
        for (auto he_it = mesh.cvoh_iter(vh_from); he_it.is_valid(); ++he_it) {
          if (mesh.to_vertex_handle(*he_it) == vh_to) {
            hedgeHandle = *he_it; // Found the halfedge from vh_from to vh_to
            break;
          }
        }
      }

      return hedgeHandle;
    };

    // add triangles to the mesh using a border edge and one of new vertices

    // TODO : mesh patcher class,
    // keep input points, keep output triangles as data member
    // two methods : GrowFromBorderEdges, GrowByRemovingEar, checking for self intersections, each does one loop over boundary edges
    auto GrowByRemovingEar = [&lookupHedge, &addedTriangles](MyMesh &mesh)
    {
      int nAdded = 0;
      bool oneAdded = false;

      for (auto hedgeHandle : mesh.halfedges())
      {
        if (mesh.is_boundary(hedgeHandle))
        {
          auto nextHedge = mesh.next_halfedge_handle(hedgeHandle);

          auto pntA = toVec(mesh.point(mesh.from_vertex_handle(hedgeHandle)));
          auto pntB = toVec(mesh.point(mesh.to_vertex_handle(hedgeHandle)));
          auto pntC = toVec(mesh.point(mesh.to_vertex_handle(nextHedge)));

          MyMesh::FaceHandle oppFacetHandle = mesh.opposite_face_handle(hedgeHandle);
          MyMesh::FaceHandle nextOppHandle = mesh.opposite_face_handle(nextHedge);

          if (!oppFacetHandle.is_valid() || !nextOppHandle.is_valid())
            continue;

          auto facetNormal = toVec(mesh.normal(oppFacetHandle));
          auto nextFacetNormal = toVec(mesh.normal(nextOppHandle));

          Triangle3d earTriangle(pntA, pntB, pntC);

          auto triNormal = earTriangle.GetNormal();

//          printf("%.6f, %.6f, %.6f\n", triNormal.dot(facetNormal), triNormal.dot(nextFacetNormal), earTriangle.AspectRatio() );

          if (triNormal.dot(facetNormal) > 0.9 && triNormal.dot(nextFacetNormal) > 0.9 && earTriangle.AspectRatio() > 0.5) {
////            printf("ear candidate found\n");
              // iterate over halfedges of vertex
              // check if edge is boundary
            std::vector<MyMesh::VertexHandle> face_vhandles;
            face_vhandles.push_back(mesh.from_vertex_handle(hedgeHandle));
            face_vhandles.push_back(mesh.to_vertex_handle(hedgeHandle));
            face_vhandles.push_back(mesh.to_vertex_handle(nextHedge));

            auto faceHandle = mesh.add_face(face_vhandles);
            if (faceHandle.is_valid())
            {

              // debug begin
 
              addedTriangles.push_back(faceHandle);
              mesh.update_face_normals();
              mesh.update_vertex_normals();
              mesh.update_normals();
              nAdded++;
              oneAdded = true;
            }
          }


          // Try triangle hedgeHandle.from_vertex, hedgeHandle.to_vertex, nextHedge.to_vertex
          // compute normal
          
          // get opposite facet of hedgeHandle, get normal
          // get opposite facet of nextHedge, get normal
          // get triangle normal

          // check if dot products are good.

          // TODO : check if intersects with other facets
          // sort by aspect ratio
        }
      }

      printf("%d ear facets added\n", nAdded);
      return oneAdded;
    };

    auto GrowFromBorderEdges = [&lookupHedge, &addedTriangles](MyMesh &mesh, const std::vector<Eigen::Vector3d> &points)
    {
      // fill holes in current mesh by adding triangles with new vertices

      // find best candidate
      std::vector<MyMesh::VertexHandle> pointVertices(points.size(), MyMesh::VertexHandle(-1));

      // TODO : extract into a function
      bool oneAdded = true;
      bool avoidSplitVertices = true;
      int iTry = 0;

      while (oneAdded)
      {

        printf("\n*************** %d ***************\n", iTry);

        for (size_t iPnt = 0; iPnt < points.size(); iPnt++)
        {
          auto pnt = points[iPnt];
          printf(" %lu. (%.6f, %.6f, %.6f)\n", iPnt, pnt[0], pnt[1], pnt[2]);
        }

        printf("\n\n");

        iTry++;
        int nBoundary = 0;
        int nAdded = 0;
        oneAdded = false;
        happly::PLYExport ccPly;

        happly::PLYExport bPly;
        for (auto hedgeHandle : mesh.halfedges())
        {
          if (mesh.is_boundary(hedgeHandle))
          {
            auto toVertex = mesh.to_vertex_handle(hedgeHandle);

            auto pnt = mesh.point(toVertex);

            bPly.AddVertex(pnt, 255, 0, 0);
          }
        }
        bPly.WritePLY("boundaryvertices.ply");


        for (auto hedgeHandle : mesh.halfedges())
        {

          double green = 0.0;

          if (mesh.is_boundary(hedgeHandle))
          {
            nBoundary++;

            auto toVertex = mesh.to_vertex_handle(hedgeHandle);
            auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

            Eigen::Vector3d toPnt = toVec(mesh.point(toVertex));
            Eigen::Vector3d fromPnt = toVec(mesh.point(fromVertex));

            MyMesh::FaceHandle oppFacetHandle = mesh.opposite_face_handle(hedgeHandle);

            if (!oppFacetHandle.is_valid())
              continue;

            auto facetNormal = toVec(mesh.normal(oppFacetHandle));

            printf("\n edge %d(%.6f, %.6f, %.6f) - %d(%.6f, %.6f, %.6f), normal (%.6f, %.6f, %.6f)\n",
                   toVertex.idx(), toPnt[0], toPnt[1], toPnt[2],
                   fromVertex.idx(), fromPnt[0], fromPnt[1], fromPnt[2],
                   facetNormal[0], facetNormal[1], facetNormal[2]);

            Line3d edgeSeg(toPnt, fromPnt);

            std::vector<std::pair<double, size_t>> candidates;

            // test all points
            // build triangle hedgeHandle.to_Vertex, hedgehandle.from_vertex, currPnt
            // check if point projects in the middle of segment  (to_vertex, from_Vertex)
            // check normals between this triangle and neighboring triangles
            // check aspect ratio

            // find best triangle, add to the surface

            // helper function lookupHedge
            // helper function for edge vector
            // helper function for facet triangle

            for (size_t iPnt = 0; iPnt < points.size(); iPnt++)
            {

              auto pnt = points[iPnt];

              auto param = edgeSeg.ProjectParam(pnt);

              //          if(iPnt == 99)
              //            printf("     %d, param %.6f\n", iPnt, param);

              if (param >= 0.0 && param <= 1.0)
              {
                Eigen::Vector3d projPnt = edgeSeg.Project(pnt);

                double dist = (projPnt - pnt).norm();

                Triangle3d triangle(fromPnt, toPnt, pnt);

                // TODO : sort by shortest distance

                double aspectRatio = triangle.AspectRatio();
                auto triNormal = triangle.GetNormal();

                std::vector<MyMesh::HalfedgeHandle> oppEdges;
                oppEdges.push_back(hedgeHandle);
                oppEdges.push_back(lookupHedge(mesh, toVertex, pointVertices[iPnt]));
                oppEdges.push_back(lookupHedge(mesh, pointVertices[iPnt], fromVertex));

                double dihedCos = 0.0;
                int nDihed = 0;

                for (auto oppEdge : oppEdges)
                {
                  if (oppEdge.is_valid() && mesh.is_boundary(oppEdge))
                  {

                    MyMesh::FaceHandle oppFacetHandle = mesh.opposite_face_handle(oppEdge);

                    if (oppFacetHandle.is_valid())
                    {
                      auto facetNormal = toVec(mesh.normal(oppFacetHandle));

                      //                 printf("%d, %.6f\n", oppFacetHandle.idx(), facetNormal.dot(triNormal));

                      dihedCos += facetNormal.dot(triNormal);
                      nDihed++;
                    }
                  }
                }

                if (nDihed)
                  dihedCos /= nDihed;

                double distThresh = 1.0;

                if (dist > eps && dist < distThresh && dihedCos > eps)
                {

                  printf(" %lu. (%.6f, %.6f, %.6f), h %.6f, ar %.6f, %d dihed cos %.6f, value %.6f\n", iPnt, pnt[0], pnt[1], pnt[2], dist, aspectRatio, nDihed, dihedCos, dist / dihedCos/*, 0.01 + weights[iPnt]*/);

                  candidates.push_back({dist / dihedCos, iPnt});

                  //              candidates.push_back({(0.01+weights[iPnt])*dist/dihedCos, iPnt});
                  //             candidates.push_back({dihedCos*aspectRatio, iPnt});
                }
                // loop for facets incident on (fromVertex, toVertex), (toVertex,
                //
                //  thisVertex), (thisVertex, fromVertex)
                // if exist then compute normals

                // now, each candidate has those criteria : average dihedral angle cos (best is 1.0), aspect ratio (best is 1.0)
              }
            }

            if (candidates.size())
            {
              std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b)
                        { return a.first < b.first; });

//              double minVal = -1.0;
//              size_t iMinCand = 0;

//              for (const auto &ccd : candidates)
//                printf("(%.6f, %lu) ", ccd.first, ccd.second);
//              printf("\n");

              for (size_t iMinCand = 0; iMinCand < candidates.size(); iMinCand++)
              {
//                if (minVal < 0.0 || candidates[iCand].first < minVal)
//                {
//                  minVal = candidates[iCand].first;
//                  iMinCand = iCand;
//                }
 //             }

                if (candidates[iMinCand].first > 0.0)
                {
                  size_t iVpoint = candidates[iMinCand].second;

                  green += 0.025;

                  auto testPnt = Eigen::Vector3d(37.2, 12.2, 10.9);
                  if ((testPnt - points[iVpoint]).norm() < 0.1)
                  {
                    printf("opanki %lu, %.6f\n", candidates[iMinCand].second, candidates[iMinCand].first);
                    ////              getchar();
                  }

                  if (!pointVertices[iVpoint].is_valid())
                    pointVertices[iVpoint] = mesh.add_vertex(MyMesh::Point(points[iVpoint][0], points[iVpoint][1], points[iVpoint][2]));

                  printf("take %lu (%.6f, %.6f, %.6f)\n", iVpoint, points[iVpoint][0], points[iVpoint][1], points[iVpoint][2]);

                  // TODO : check if any of vertices in {fromVertex, toVertex, pointVertices[iVpoint]} become split
                  std::vector<std::pair<MyMesh::VertexHandle, MyMesh::VertexHandle>> newEdges = 
                    {{fromVertex, toVertex}, {toVertex, pointVertices[iVpoint]}, {pointVertices[iVpoint], fromVertex}};
                  
                  printf("  new facet (%d, %d, %d)\n", fromVertex.idx(), toVertex.idx(), pointVertices[iVpoint].idx());

                  auto countBoundaryEdges = [](MyMesh &mesh, MyMesh::VertexHandle vertexHandle,
                                              const std::vector<std::pair<MyMesh::VertexHandle, MyMesh::VertexHandle>> &newEdges)
                  {

                    int nBoundary = 0;
                    int nEdges = 0;

                    // iterate over half edges of a vertex
                    printf("  vertex %d\n", vertexHandle.idx());

                    for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(vertexHandle); voh_it.is_valid(); ++voh_it)
                    {
                      nEdges++;
                      MyMesh::HalfedgeHandle hedge = *voh_it;

                      if (mesh.is_boundary(hedge) || mesh.is_boundary(mesh.opposite_halfedge_handle(hedge)))
                      {
                        nBoundary++;
                        auto fromVertex = mesh.from_vertex_handle(hedge);
                        auto toVertex = mesh.to_vertex_handle(hedge);

                        printf("   boundary (%d, %d)\n", fromVertex.idx(), toVertex.idx());

                        // get tovertex
                        // get fromvertex

                        // boundary edge becoming two sided edge after adding facet consisting of newEdges
                        for(const auto &newEdge : newEdges)
                        {
                          printf("    new edge (%d, %d)\n", newEdge.first.idx(), newEdge.second.idx());

                          if( (newEdge.second.idx() == fromVertex.idx() && newEdge.first.idx() == toVertex.idx()) ||             
                              (newEdge.first.idx() == fromVertex.idx() && newEdge.second.idx() == toVertex.idx()) ) {
                            printf("    cancels out\n");
                            nBoundary -= 2;
                          }
                        }
                      }
                    }

                    if(nEdges)
                      nBoundary += 2;

                    printf("  %d edges, %d boundary edges\n", nEdges, nBoundary);
                    return nBoundary;
                  };

                  int nFrom = countBoundaryEdges(mesh, fromVertex, newEdges);
                  int nTo   = countBoundaryEdges(mesh, toVertex, newEdges);
                  int nMid  = countBoundaryEdges(mesh, pointVertices[iVpoint], newEdges);

                  int n4Plus = 0;
                  int n6Plus = 0;

                  if(nFrom > 2)
                    n4Plus++;
                  if(nTo > 2)
                    n4Plus++;
                  if(nMid > 2)
                    n4Plus++;

                  if(nFrom > 4)
                    n6Plus++;
                  if(nTo > 4)
                    n6Plus++;
                  if(nMid > 4)
                    n6Plus++;

                  bool splitOk = n4Plus == 0;

                  if(!avoidSplitVertices && n6Plus == 0)
                    splitOk = true;

                  if (splitOk) {

                      // iterate over halfedges of vertex
                      // check if edge is boundary
                    std::vector<MyMesh::VertexHandle> face_vhandles;
                    face_vhandles.push_back(fromVertex);
                    face_vhandles.push_back(toVertex);
                    face_vhandles.push_back(pointVertices[iVpoint]);

                    if( !DoTrianglesIntersect(mesh, face_vhandles, addedTriangles) ) {
                      auto faceHandle = mesh.add_face(face_vhandles);

                      if (faceHandle.is_valid())
                      {
             if(faceHandle.idx() == 13263 || faceHandle.idx() == 13322) {

              printf("facet %d\n", faceHandle.idx());
                happly::PLYExport ccPly;

                auto v0 = ccPly.AddVertex(mesh.point(pointVertices[iVpoint])/*points[iVpoint]*/, 1.0, green, 1.0);
                auto v1 = ccPly.AddVertex(mesh.point(fromVertex), 1.0, green, 0.5);
                auto v2 = ccPly.AddVertex(mesh.point(toVertex), 1.0, green, 0.5);
                ccPly.AddFacet(v0, v1, v2);

                ccPly.WritePLY("debug.ply");

                // why triangles not intersecting?
                getchar();
              }



                        auto v0 = ccPly.AddVertex(points[iVpoint], 1.0, green, 1.0);
                        auto v1 = ccPly.AddVertex(fromPnt, 1.0, green, 0.5);
                        auto v2 = ccPly.AddVertex(toPnt, 1.0, green, 0.5);
                        ccPly.AddFacet(v0, v1, v2);

                        addedTriangles.push_back(faceHandle);
                        mesh.update_face_normals();
                        mesh.update_vertex_normals();
                        mesh.update_normals();
                        nAdded++;
                        oneAdded = true;
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        printf("%d boundary edges, %d triangles\n", nBoundary, nAdded);

        bool saveCandidates = true;

        if(saveCandidates) {
          static int nCandidates = 0;
          ccPly.WritePLY(std::string("candidates").append(std::to_string(nCandidates++)).append(std::string(".ply")));
          printf("write out candidates\n");
        }

        if(!oneAdded && avoidSplitVertices) {
          oneAdded = true;
          avoidSplitVertices = false;
          printf("next cycle with split vertices\n");
        }
        else if (oneAdded)
          avoidSplitVertices = true;

//        getchar();
      }
    };

    GrowFromBorderEdges(mesh, points);
//    while(GrowByRemovingEar(mesh)) {}
//    GrowFromBorderEdges(mesh, points);

    getchar();

    OpenMesh::IO::write_mesh(mesh, "allcandidates.ply");

    auto adjMatrix = MeshSamplingEdges(points);

    auto normals = PointSetComputeNormals(points, 8);
    // TODO : exclude edges "too parallel" to point set normals

    auto copyAdjMatrix = adjMatrix;

    int nEdges = 0;
    int nSkip = 0;
    for(size_t iPnt = 0; iPnt < adjMatrix.size(); iPnt++)
      for(size_t jPnt : adjMatrix[iPnt]) {
        nEdges++;
        Eigen::Vector3d edgeVec = points[iPnt]-points[jPnt];
        edgeVec.stableNormalize();

        double dot1 = fabs(edgeVec.dot(normals[iPnt]));
        double dot2 = fabs(edgeVec.dot(normals[jPnt]));

        if(dot1 > 0.8 || dot2 > 0.8) {
          auto iErase = copyAdjMatrix[iPnt].erase(jPnt);
          auto jErase = copyAdjMatrix[jPnt].erase(iPnt);
          if(iErase || jErase)
            nSkip++;
        }
      }

    printf("%d edges, skip %d\n", nEdges, nSkip);

    MeshSamplingSave(points, adjMatrix, "origedges.ply");
    MeshSamplingSave(points, copyAdjMatrix, "normaledges.ply");

    MyMesh patchMesh = MeshSamplingCycles(points, copyAdjMatrix);

    OpenMesh::IO::write_mesh(patchMesh, "patchmesh.ply");
}


  // PatchLoops(mesh, 0.125);
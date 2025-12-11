#include "common.h"
#include "tools3d.h"
#include "AABBTree.h"
#include "MeshSampling.h"
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

    std::vector<Eigen::Vector3d> allpoints;

    // use grid to find points "filling" the gap

// TODO :     auto fillLoopPoints = []()
    {
    happly::PLYExport cPly;

    cGrid3d grid;
    grid.init(bounds, loopStep);

    for(size_t ii = 0; ii < grid.GetNi(); ii++)
      for(size_t jj = 0; jj < grid.GetNj(); jj++)
        for(size_t kk = 0; kk < grid.GetNk(); kk++) {
      
          auto pos = grid.GetPosition(ii, jj, kk);

          auto np = tree.FindNearestPoint(pos);

          if(np.facetIndex != -1) {
          // TODO : check normals
            Eigen::Vector3d normal = toVec(mesh.normal(MyMesh::FaceHandle(np.facetIndex)));
            Eigen::Vector3d vec = np.pnt-pos;

            grid.GetVector(ii, jj, kk) = vec;

            vec.stableNormalize();
            double dotProd = vec.dot(normal);

            if(fabs(dotProd) > 0.95)
              grid.GetPinned(ii, jj, kk) = true;
          }
        }
  
    grid.Smooth(10000);

    happly::PLYExport sPly;

    for(size_t ii = 0; ii < grid.GetNi(); ii++)
      for(size_t jj = 0; jj < grid.GetNj(); jj++)
        for(size_t kk = 0; kk < grid.GetNk(); kk++) {

          Eigen::Vector3d vec1 = grid.GetVector(ii, jj, kk);
          Eigen::Vector3d pnt1 = grid.GetPosition(ii, jj, kk);

          if(!grid.GetPinned(ii, jj, kk)) {
            sPly.AddVertex(pnt1+vec1, 1.0, 0.1, 0.2);
            allpoints.push_back(pnt1+vec1);
          }

          for(int iNeighbor = 0; iNeighbor < 6; iNeighbor++) {

            size_t iN = ii;
            size_t jN = jj;
            size_t kN = kk;

            switch(iNeighbor) {
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

            if(grid.IsValid(iN, jN, kN) && !grid.GetPinned(ii, jj, kk) && !grid.GetPinned(iN, jN, kN)) {
                Eigen::Vector3d vec2 = grid.GetVector(iN, jN, kN);
                Eigen::Vector3d pnt2 = grid.GetPosition(iN, jN, kN);

                // find the parameter value that minimizes the squared length of linear interpolation of vec1 and vec2

                double AA = vec1.dot(vec1);
                double BB = vec2.dot(vec2);
                double AB = vec1.dot(vec2);

                double denom = AA+BB-2.0*AB;

                if(denom > eps) {
                  double numer = BB-AB;

                  double param = numer/denom;

                  if(param > 0.0 && param < 1.0) {
                    Eigen::Vector3d midPnt = pnt1*param + pnt2*(1.0-param);
                    Eigen::Vector3d midVec = vec1*param + vec2*(1.0-param);

                    auto np = tree.FindNearestPoint(midPnt);

                    MyMesh::FaceHandle facetHandle(np.facetIndex);

                    bool isOnBorderFacet = false;
                    for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
                      if(borderVertices.find(*fv_it) != borderVertices.end())
                        isOnBorderFacet = true;
                    }
                              
                    if(isOnBorderFacet) {
//                    cPly.AddVertex(pnt1, 0.0, 1.0, 0.0);
//                    cPly.AddVertex(pnt2, 0.0, 1.0, 0.0);
                      cPly.AddVertex(midPnt+midVec, 0.6, 0.2, 0.7);

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
  
      // TODO : filter out strangepoints that don't project onto border facets
      cPly.WritePLY("fillerpoints.ply");
    }

    happly::PLYExport exportPly;

    // cluster filler points
    auto points = PointSetCluster(allpoints, 100, 2000);

    printf("%lu border points, %lu new vertex points\n", borderPoints.size(), points.size());

// add border points to the point set
//    for(auto bpoint : borderPoints)
//      points.push_back(bpoint);


    // remove duplicate points
    PointSetRemoveDuplicates(points, 1e-5);

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
      for(size_t ii = 0; ii < nBorderPoints-1; ii++)
        for(size_t jj = ii+1; jj < nBorderPoints; jj++) {
      
          std::vector<std::size_t> shortestPath = findShortestPath(patchPoints, distanceMatrix, ii, jj);


          //  mark shortest path points 
          pathPoints[ii] = true;
          pathPoints[jj] = true;

          double averLen = 0.0;

          for(size_t ii = 1; ii < shortestPath.size(); ii++) {
            averLen += (patchPoints[shortestPath[ii]]-patchPoints[shortestPath[ii-1]]).norm();
          }

          if(averLen > 0.0)
            averLen /= (shortestPath.size()-1);

          printf("%lu path points, aver. length %.6f\n", shortestPath.size(), averLen);

          for(auto iPath : shortestPath) {
          //  printf("   %lu", iPath); 
            pathPoints[iPath] = true;
          }
          //printf("\n");
        }

      // save points that have pathPoints flags set;

      points.clear();
      happly::PLYExport eply;

      for(size_t ii = 0; ii < patchPoints.size(); ii++) {
        uint8_t blue = ii < nBorderPoints ? 255 : 0;
        uint8_t green = ii < nBorderPoints ? 0 : pathPoints[ii] ? 255 : 0;
        uint8_t red = ii < nBorderPoints ? 0 : pathPoints[ii] ? 0 : 255;

        eply.AddVertex(patchPoints[ii], red, green, blue);

        if(ii >= nBorderPoints && pathPoints[ii])
          points.push_back(patchPoints[ii]);
      }

      eply.WritePLY("patchpathpoints.ply");

      printf("%lu points left\n", points.size());
    }

    // TODO : Taubin smooth point cloud while pinning path points 

    // fill holes in current mesh by adding triangles with new vertices

    // find best candidate
    std::vector<MyMesh::VertexHandle> pointVertices(points.size(), MyMesh::VertexHandle(-1));

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

    // TODO : extract into a function
    bool oneAdded = true;
    int iTry = 0;

    while(oneAdded) {

     printf("\n*************** %d ***************\n", iTry);


     for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {
      auto pnt = points[iPnt];
      printf(" %lu. (%.6f, %.6f, %.6f)\n", iPnt, pnt[0], pnt[1], pnt[2]);
     }

     printf("\n\n");


     iTry++;
     oneAdded = false;
     happly::PLYExport ccPly;

     for(auto hedgeHandle : mesh.halfedges()) {

      double green = 0.0;

      if(mesh.is_boundary(hedgeHandle)) {

        auto toVertex = mesh.to_vertex_handle(hedgeHandle);
        auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

        Eigen::Vector3d toPnt = toVec(mesh.point(toVertex));
        Eigen::Vector3d fromPnt = toVec(mesh.point(fromVertex));


        MyMesh::FaceHandle oppFacetHandle = mesh.opposite_face_handle(hedgeHandle);

        if(!oppFacetHandle.is_valid())
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


        for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

          auto pnt = points[iPnt];

          auto param = edgeSeg.ProjectParam(pnt);

//          if(iPnt == 99)
//            printf("     %d, param %.6f\n", iPnt, param);

          if(param >= 0.0 && param <= 1.0) {
            Eigen::Vector3d projPnt = edgeSeg.Project(pnt);

            double dist = (projPnt-pnt).norm();

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

            for(auto oppEdge : oppEdges) {
              if(oppEdge.is_valid() && mesh.is_boundary(oppEdge)) {

                MyMesh::FaceHandle oppFacetHandle = mesh.opposite_face_handle(oppEdge);

                if(oppFacetHandle.is_valid()) {
                  auto facetNormal = toVec(mesh.normal(oppFacetHandle));




 //                 printf("%d, %.6f\n", oppFacetHandle.idx(), facetNormal.dot(triNormal));

                  dihedCos += facetNormal.dot(triNormal);
                  nDihed++;
                }
              }
            }

            if(nDihed)
              dihedCos /= nDihed;

            double distThresh = 1.0;

            if(dist > eps && dist < distThresh && dihedCos > eps) {

              printf(" %lu. (%.6f, %.6f, %.6f), h %.6f, ar %.6f, %d dihed cos %.6f, value %.6f, weight %.6f\n", iPnt, pnt[0], pnt[1], pnt[2], dist, aspectRatio, nDihed, dihedCos, dist/dihedCos, 0.01+weights[iPnt]);

                candidates.push_back({dist/dihedCos, iPnt});

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

        if(candidates.size()) {
          std::sort(candidates.begin(), candidates.end(), [](const auto &a, const auto &b) { return a.first < b.first;});

          double minVal = -1.0;
          size_t iMinCand = 0;

          // TODO : why unsorted???
          for(const auto &ccd : candidates)
            printf("(%.6f, %lu) ", ccd.first, ccd.second);
          printf("\n");

          for(size_t iCand = 0; iCand < candidates.size(); iCand++) {
            if(minVal < 0.0 || candidates[iCand].first < minVal) {
              minVal = candidates[iCand].first;
              iMinCand = iCand;
            }
          }

          if(candidates[iMinCand].first > 0.0) {
            size_t iVpoint = candidates[iMinCand].second;

            auto v0 = ccPly.AddVertex(points[iVpoint], 1.0, green, 1.0);
            auto v1 = ccPly.AddVertex(fromPnt, 1.0, green, 0.5);
            auto v2 = ccPly.AddVertex(toPnt, 1.0, green, 0.5);

            green += 0.025;

            auto testPnt = Eigen::Vector3d(37.2, 12.2, 10.9);
            if((testPnt-points[iVpoint]).norm() < 0.1) {
              printf("opanki %lu, %.6f\n", candidates[iMinCand].second, candidates[iMinCand].first);
////              getchar();
            }

            ccPly.AddFacet(v0, v1, v2);

            if(!pointVertices[iVpoint].is_valid())
              pointVertices[iVpoint] = mesh.add_vertex(MyMesh::Point(points[iVpoint][0], points[iVpoint][1], points[iVpoint][2]));

            printf("take %lu (%.6f, %.6f, %.6f)\n", iVpoint, points[iVpoint][0], points[iVpoint][1], points[iVpoint][2]);

            std::vector<MyMesh::VertexHandle> face_vhandles;
            face_vhandles.push_back(fromVertex);
            face_vhandles.push_back(toVertex);
            face_vhandles.push_back(pointVertices[iVpoint]);

            auto faceHandle = mesh.add_face(face_vhandles);
            if(faceHandle.is_valid()) {
              mesh.update_face_normals();
              mesh.update_vertex_normals();
              mesh.update_normals();
              oneAdded = true;
            }
          }
        }
      }
    }

//    {
//      static int nCandidates = 0;
//      ccPly.WritePLY(std::string("candidates").append(std::to_string(nCandidates++)).append(std::string(".ply")));
//      printf("write out candidates\n");
//    }
  }

  OpenMesh::IO::write_mesh(mesh, "allcandidates.ply");

  auto adjMatrix = MeshSamplingEdges(points);
  MyMesh patchMesh = MeshSamplingCycles(points, adjMatrix);

  OpenMesh::IO::write_mesh(patchMesh, "patchmesh.ply");
}


  // PatchLoops(mesh, 0.125);
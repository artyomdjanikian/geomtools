#include "common.h"
#include "tools3d.h"
#include "cOctree.h"
#include "AABBTree.h"
#include "happly.h"

// tests octree by storing mesh triangles in the octree of limited depth,
// then extracting boundary faces of the octree, them creating a triangle mesh
// out of the boundary faces, then projecting them onto the original mesh

void OctreeTest(MyMesh &mesh, const AABBTree &tree, int nIter)
{
    printf("Octree test\n");

    cOctree octree;

    std::vector<Eigen::Vector3d> pnts;

    for(auto facetHandle : mesh.faces()) {
        pnts.clear();

        for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {

  //        auto vcolor = mesh.color(*fv_it);

  //        printf("%d, %d, %d\n", vcolor[0], vcolor[1], vcolor[2]);

          auto pnt = mesh.point(*fv_it);

          pnts.push_back(toVec(pnt));
        }

        if(pnts.size() >= 3)
          for(size_t iTria = 0; iTria < pnts.size()-2; iTria++)
            octree.AddTriangle(Triangle3d(pnts[0], pnts[iTria+1], pnts[iTria+2]));
    }

    octree.Start();

    size_t limitDepth = nIter;

    int ii = 0;
    for(auto iter = octree.begin(); iter.getDepth(); ++iter) {
    //    printf(">>> %d\n", ii);

    //    printf(">>>   %d/%d\n", iter.getedgeCountDepth(), iter.isValid());

    //    for(int iC = 0; iC < iter.getDepth(); iC++)
    //      printf(" ");

    //    printf("%d. \n", ii++);

      if(iter.getDepth() < limitDepth) {
        if( iter.isLeaf() && !octree.Subdivide(iter) )
          printf("!!! subdivide failed !!!\n");
      }
    }

    // TODO : extract into printTree
  #if 0 
    ii = 0;
    for(auto iter = octree.begin(); iter.getDepth(); ++iter) {
      printf(">>>");
      for(size_t iC = 0; iC < iter.getDepth(); iC++)
        printf(" ");

      printf("%d/%d. ", ii++, iter.isValid());

      if(iter.isLeaf())
        printf("%lu facets\n", iter.getFacets().size());
      else {
        int nChildren = 0;

        for(int iC = 0; iC < 8; iC++)
          if(iter.hasChild(iC))
            nChildren++;

        printf("%d children\n", nChildren);

        if(iter.getFacets().size() > 0)
          printf("    !!! facets !!!\n");

      }
    }
  #endif

    printf("detect boundary faces\n");

    int nBoundaryFaces = 0;
    int nLeaves = 0;

    MyMesh voxelMesh;
    std::vector<MyMesh::VertexHandle> voxelMeshVertices;

    std::map<std::pair<uint64_t, uint64_t>, int> edgeCount;
    std::map<uint64_t, Eigen::Vector3d> vertexCoords;
    std::map<uint64_t, uint64_t> vertexIndices;
    std::vector<std::array<uint64_t, 4>> faces;

    for(auto iter = octree.begin(); iter.getDepth(); ++iter) {
      if(iter.isLeaf()) {

        // iterate over 6 faces
        nLeaves++;

        const auto &bounds = iter.getBounds();
        auto center = (bounds.minBound+bounds.maxBound)*0.5;

        for(uint8_t iFace = 0; iFace < 6; iFace++) {

          std::array<Eigen::Vector3d, 4> facePoints = iter.getFaceCoords(iFace);

          Eigen::Vector3d faceCenter(0.0, 0.0, 0.0);

          for(int iV = 0; iV < 4; iV++)
            faceCenter += facePoints[iV]*0.25;

          faceCenter += (faceCenter - center)*0.1;

          // check if face center is inside a tree leaf
          bool inLeaf = false;
          for(auto iiter = octree.begin(); !inLeaf && iiter.getDepth(); ++iiter) {
            if(iiter.getBounds().IsPointInside(faceCenter)) {
              if(iiter.isLeaf())
                inLeaf = true;
            }
            else
              iiter.skip();
          }

          if(!inLeaf) {
            auto faceVertices = iter.getFaceIndices(iFace);
            faces.push_back(faceVertices);

            for(int iEdge = 0; iEdge < 4; iEdge++) {

              vertexCoords[faceVertices[iEdge]] = facePoints[iEdge];

              auto vertexId = vertexIndices.size();

              auto findIter = vertexIndices.find(faceVertices[iEdge]);

              if(findIter == vertexIndices.end()) {
                voxelMeshVertices.push_back(voxelMesh.add_vertex(MyMesh::Point(facePoints[iEdge][0], facePoints[iEdge][1], facePoints[iEdge][2])));

                auto insertIter = vertexIndices.insert({faceVertices[iEdge], vertexId});

                // was inserted
  //             meshVertexPositions.push_back({facePoints[iEdge][0], facePoints[iEdge][1], facePoints[iEdge][2]});
              }

              int iNextEdge = iEdge+1;
              if(iNextEdge == 4)
                iNextEdge = 0;

              edgeCount[{faceVertices[iEdge], faceVertices[iNextEdge]}]++;
            }

            nBoundaryFaces++;

          }
        }
      }
    }

    printf("%d leaves, %d boundary faces\n", nLeaves, nBoundaryFaces);
  //  printf("%lu/%lu vertex positions\n", vertexCoords.size(), meshVertexPositions.size());

    std::map<uint64_t, bool> isNonManifold;

    int nNmHalfedges = 0;
    for(const auto &edgeC : edgeCount) {
  //    printf("(%lu, %lu) => %d\n", edgeC.first.first, edgeC.first.second, edgeC.second);
      if(edgeC.second > 1) {
        nNmHalfedges++;
        isNonManifold[edgeC.first.first] = true;
        isNonManifold[edgeC.first.second] = true;
      }
    }

    // add vertices
    int nNonManifold = 0;

    // add faces
    std::vector<MyMesh::VertexHandle> face_vhandles;

    for(const auto &face : faces) {
      std::vector<size_t> faceIndices;

      // get 4 vertices 
      for(int ii = 0; ii < 4; ii++) {
        // if manifold then iV = vertexIndices[face[ii]]
        // if non-manifold then add with duplicate coordinate
        if(isNonManifold[face[ii]]) {
            nNonManifold++;
            auto vertexCoord = vertexCoords[face[ii]];

  //          printf("%.6f, %.6f, %.6f\n", vertexCoord[0], vertexCoord[1], vertexCoord[2]);
            faceIndices.push_back(voxelMeshVertices.size());
  //          meshVertexPositions.push_back({vertexCoord[0], vertexCoord[1], vertexCoord[2]});

            voxelMeshVertices.push_back(voxelMesh.add_vertex(MyMesh::Point(vertexCoord[0], vertexCoord[1], vertexCoord[2])));
        }
        else {
          auto iV = vertexIndices[face[ii]];
          faceIndices.push_back(iV);
        }
      }

      // add face with 4 indices
  //    meshFaceIndices.push_back(faceIndices);

      face_vhandles.clear();
      for(auto iV : faceIndices)
        face_vhandles.push_back(voxelMeshVertices[iV]);

      voxelMesh.add_face(face_vhandles);
    }

    OpenMesh::HPropHandleT<int> opaqueEdges;
    voxelMesh.add_property(opaqueEdges);

    OpenMesh::FPropHandleT<int> pids;
    voxelMesh.add_property(pids);

    int nPids = FloodFill(voxelMesh, opaqueEdges, -1, pids);

    printf("%d pids\n", nPids);

    std::vector<int> pidCount;
    pidCount.resize(nPids);

    for(auto facetHandle : voxelMesh.faces()) {
      auto facetPid = voxelMesh.property(pids, facetHandle);

      pidCount[facetPid]++;
    }

    int maxCount = 0;
    for(auto pCount : pidCount)
      maxCount = std::max(maxCount, pCount);

    std::vector<bool> isActive(voxelMesh.n_vertices(), false);

    int nFacetsToDelete = 0;
    for(auto facetHandle : voxelMesh.faces()) {
      auto facetPid = voxelMesh.property(pids, facetHandle);

      if(pidCount[facetPid] < maxCount)
        nFacetsToDelete++;
      else {
        // iterate over vertices
        for(MyMesh::FaceVertexIter fv_it = voxelMesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
          int vId = fv_it->idx();
          isActive[vId] = true;;
        }
      }
    }

    printf("%d facets to delete\n", nFacetsToDelete);
    // TODO : discard faces of small pids

    std::vector<Eigen::Vector3d> meshPoints(voxelMesh.n_vertices());
    std::vector<Eigen::Vector3d> closestPoints(voxelMesh.n_vertices());
    std::vector<bool> isTouched(voxelMesh.n_vertices(), true);

    for(int iTry = 0; iTry < 10; iTry++) {

      printf("%d. project %lu vertices\n", iTry, voxelMesh.n_vertices());
  
      // collect vertices, project onto original mesh
      for(auto vertexHandle : voxelMesh.vertices()) {
        auto pnt = voxelMesh.point(vertexHandle);
        meshPoints[vertexHandle.idx()] = Eigen::Vector3d(pnt[0], pnt[1], pnt[2]);
      }

      closestPoints.resize(meshPoints.size());

      double totalDist = 0.0;
  #pragma omp parallel for reduction(+ : totalDist)
      for(int iPnt = 0; iPnt < static_cast<int>(meshPoints.size()); iPnt++) {
        if(isTouched[iPnt]) {
          auto closestPoint = tree.FindNearestPoint(meshPoints[iPnt]);
          totalDist += fabs(closestPoint.dist);
          closestPoints[iPnt] = closestPoint.pnt;
        }
        else
          closestPoints[iPnt] = meshPoints[iPnt];
      }

      printf("   %.6f aver. distance\n", totalDist/meshPoints.size());

      bool shouldOffset = iTry > 0 && iTry%3 == 0;
      bool shouldSmooth = !shouldOffset;

      if(shouldSmooth) {

        printf("%d. smooth\n", iTry);

        // smooth if brings closer to original mesh
        int nCloser = 0;

        int nMoved = 0;
        double moveDist = 0.0;

        for(auto vertexHandle : voxelMesh.vertices()) {
          int thisIdx = vertexHandle.idx();
          if(isActive[thisIdx]) {

            int valence = 0;
            Eigen::Vector3d averPnt(0.0, 0.0, 0.0);

            for (auto vv_it=voxelMesh.vv_iter( vertexHandle ); vv_it.is_valid(); ++vv_it) {
              auto pnt = voxelMesh.point(*vv_it);
              averPnt += Eigen::Vector3d(pnt[0], pnt[1], pnt[2]);
              ++valence;
            }

            if(valence)
              averPnt *= 1.0/valence;

            double initDist = (meshPoints[thisIdx]-closestPoints[thisIdx]).norm();
            double smoothDist = (averPnt-closestPoints[thisIdx]).norm();

            if(initDist <= smoothDist) {
                auto closestPoint = tree.FindNearestPoint(averPnt);
                smoothDist = fabs(closestPoint.dist);
            }

            if(initDist > smoothDist) {
              moveDist += (meshPoints[thisIdx] - averPnt).norm();
              meshPoints[thisIdx] = averPnt;
      //        meshVertexPositions[thisIdx] = {averPnt[0], averPnt[1], averPnt[2]};
              nCloser++;
            }
          }
        }

        printf("   smoothing brings %d vertices closer, aver. move %.6f\n", nCloser, moveDist/nCloser);

        for(auto vertexHandle : voxelMesh.vertices()) {
          const auto &smoothPnt = meshPoints[vertexHandle.idx()];
          MyMesh::Point pnt(smoothPnt[0], smoothPnt[1], smoothPnt[2]);
          voxelMesh.point(vertexHandle) = pnt;
        }
      }

      if(shouldOffset) {

        printf("%d. offset\n", iTry);

        voxelMesh.request_vertex_normals();
  //     voxelMesh.update_vertex_normals();

        double totalOffset = 0.0;

        for(auto vertexHandle : voxelMesh.vertices()) {
          auto meshPnt = meshPoints[vertexHandle.idx()];
          auto closestPnt = closestPoints[vertexHandle.idx()];

  //        auto vpoint = voxelMesh.point(vertexHandle);
  //       auto vnormal = voxelMesh.normal(vertexHandle);

  ///        printf("%.6f, %.6f, %.6f\n", vnormal[0], vnormal[1], vnormal[2]);

  //        Eigen::Vector3d normal(vnormal[0], vnormal[1], vnormal[2]);
          Eigen::Vector3d surfVec = closestPnt-meshPnt;
          Eigen::Vector3d normal = surfVec;
          normal.stableNormalize();

          double dotProd = normal.dot(surfVec);
          totalOffset += 0.5*dotProd;

          normal *= 0.5*dotProd;

          meshPnt += normal;

  //        printf("%.6f, %.6f, %.6f\n", meshPnt[0], meshPnt[1], meshPnt[2]);

          MyMesh::Point pnt(meshPnt[0], meshPnt[1], meshPnt[2]);
          voxelMesh.point(vertexHandle) = pnt;
        }

        printf("   %lu vertices offset, %.6f aver. dist\n", voxelMesh.n_vertices(), totalOffset/voxelMesh.n_vertices());
      }
    }

    printf("%d facets to delete\n", nFacetsToDelete);

    // collect project onto mesh with tree
    OpenMesh::IO::write_mesh(voxelMesh, "voxelmesh.ply");

    printf("%lu vertices, %lu half edges, %d non manifold half edges\n", vertexCoords.size(), edgeCount.size(), nNmHalfedges);
    printf("%d non-manifold vertices\n", nNonManifold);

    // Write the object to file
    happly::PLYData plyOut;

    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::array<double, 3>> meshVertexColors;
    std::vector<std::vector<size_t>> meshFaceIndices;

    for(auto vertexHandle : voxelMesh.vertices()) {
      if(meshVertexPositions.size() <= static_cast<size_t>(vertexHandle.idx()))
        meshVertexPositions.resize(vertexHandle.idx()+1);

      auto pnt = voxelMesh.point(vertexHandle);
      meshVertexPositions[vertexHandle.idx()] = {pnt[0], pnt[1], pnt[2]};
    }

    meshVertexColors.resize(meshVertexPositions.size());

    // average edge length

    printf("compute average edge length\n");
    double averEdgeLen = 0.0;
    int nEdges = 0;

    for(auto hedgeHandle : voxelMesh.halfedges()) {
      nEdges++;

      auto toVertex = voxelMesh.to_vertex_handle(hedgeHandle);
      auto fromVertex = voxelMesh.from_vertex_handle(hedgeHandle);

      auto edgeLen = (voxelMesh.point(toVertex) - voxelMesh.point(fromVertex)).norm();

      averEdgeLen += edgeLen;
    }

    if(nEdges)
      averEdgeLen /= nEdges;

  #pragma omp parallel for
    for(int iPnt = 0; iPnt < static_cast<int>(meshVertexPositions.size()); iPnt++) {
      Eigen::Vector3d vertexPoint(meshVertexPositions[iPnt][0], meshVertexPositions[iPnt][1], meshVertexPositions[iPnt][2]);
      auto closestPoint = tree.FindNearestPoint(vertexPoint);
      double dist = fabs(closestPoint.dist)/averEdgeLen;

      meshVertexColors[iPnt][0] = dist;
      meshVertexColors[iPnt][1] = 4.0*dist*(1.0-dist);
      meshVertexColors[iPnt][2] = 1.0-dist;
    }

    printf("%lu vertices\n", meshVertexPositions.size());

    for(auto facetHandle : voxelMesh.faces()) {

      auto pid = voxelMesh.property(pids, facetHandle);

      if(pidCount[pid] == maxCount) {
        std::vector<size_t> faceVertices;

        for(MyMesh::FaceVertexIter fv_it = voxelMesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
          int vId = fv_it->idx();
//          auto vcolor = mesh.color(*fv_it);

//          printf("%d, %d, %d\n", vcolor[0], vcolor[1], vcolor[2]);

          auto pnt = mesh.point(*fv_it);

          pnts.push_back(toVec(pnt));
          faceVertices.push_back(vId);
        }

        meshFaceIndices.push_back(faceVertices);
      }
    }

    // Add mesh data (elements are created automatically)
    plyOut.addVertexPositions(meshVertexPositions);
    plyOut.addVertexColors(meshVertexColors);
    plyOut.addFaceIndices(meshFaceIndices);

    // TODO : color code by distance


    printf("%lu ply vertices, %lu ply faces\n", meshVertexPositions.size(), meshFaceIndices.size());

    plyOut.write("voxelboundary.ply", happly::DataFormat::ASCII);
}

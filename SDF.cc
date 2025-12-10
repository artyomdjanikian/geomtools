// computing shape diameter function by tracing and reflecting rays inside the model (billiard style)

#include "common.h"
#include "tools3d.h"
#include "AABBTree.h"
#include "happly.h"

void ShapeDiameterFuncTest(MyMesh &mesh, AABBTree &tree, int nIter)
{
    happly::PLYData plyOut;

    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::array<double, 3>> meshVertexColors;

    struct sRay {
      Eigen::Vector3d origin;
      Eigen::Vector3d dir;
      double dist;
      int facetIndex;
    };

    std::vector<sRay> rays;
    std::vector<std::pair<double, int>> facetSDF(rays.size());
    std::vector<Eigen::Vector3d> pnts;

    for(auto facetHandle : mesh.faces()) {
      Eigen::Vector3d centroid(0.0, 0.0, 0.0);

      size_t iV = 0;
      for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
          auto pnt = mesh.point(*fv_it);

          centroid += Eigen::Vector3d(pnt[0], pnt[1], pnt[2]);
          iV++;
      }

      if(iV)
        centroid *= 1.0/iV;

      pnts.push_back(centroid);

      auto norm = mesh.normal(facetHandle);
      Eigen::Vector3d normal(norm[0], norm[1], norm[2]); 

      centroid -= normal * 1e-05;

      rays.push_back({centroid, normal, 0.0, -1});

      facetSDF.push_back({0.0, 0});
    }

    // TODO : debug TraceRay performance`
    for(int iTry = 0; iTry < nIter; iTry++) {

    //    printf("%.6f, %.6f, %.6f\n", origin[0], origin[1], origin[2]);
    int nMissed = 0;

    #pragma omp parallel for reduction( + : nMissed)
    for(int iPnt = 0; iPnt < static_cast<int>(rays.size()); iPnt++) {

        auto &origin = rays[iPnt].origin;
        auto &dir = rays[iPnt].dir;

        std::vector<AABBTree::sRayX> rayXVec;
        auto found = tree.TraceRay(origin, dir, rayXVec);

        if(!found)
          nMissed++;
  //      if(!found) {
  //        printf("\n%d. !!! missed !!!\n", iTry);
  //      }
      
        if(found) {
    //      printf("\n%d. facet %d, dist %.6f, pnt (%.6f, %.6f, %.6f)\n",
    //        iTry, rayXVec[0].facetIndex, rayXVec[0].dist, rayXVec[0].pnt[0], rayXVec[0].pnt[1], rayXVec[0].pnt[2]);

          auto norm = mesh.normal(MyMesh::FaceHandle(rayXVec[0].facetIndex));

    //      printf("normal (%.6f, %.6f, %.6f)\n", norm[0], norm[1], norm[2]);
        
          Eigen::Vector3d normal(norm[0], norm[1], norm[2]); 

          Eigen::Vector3d shiftVec = rayXVec[0].pnt-origin;

          rays[iPnt].dist = shiftVec.norm();
          rays[iPnt].facetIndex = rayXVec[0].facetIndex;

          double dotProd = shiftVec.dot(normal);

          shiftVec -= normal*(2.0*dotProd);
          shiftVec.stableNormalize();

          dir = shiftVec;
          origin = rayXVec[0].pnt;

          Eigen::Vector3d centroid(0.0, 0.0, 0.0);
          int nPnts = 0;

          for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(MyMesh::FaceHandle(rayXVec[0].facetIndex)); fv_it.is_valid(); ++fv_it) {
            auto point = mesh.point(*fv_it);
            centroid += Eigen::Vector3d(point[0], point[1], point[2]);
            nPnts++;
          }

    //      printf("%d points in facet\n", nPnts);
          if(nPnts)
            origin = origin * 0.5 + centroid * (0.5/nPnts);

          origin += dir*0.00001;
    //      printf("new origin (%.6f, %.6f, %.6f)\n", origin[0], origin[1], origin[2]);
    //      printf("new direction (%.6f, %.6f, %.6f)\n", dir[0], dir[1], dir[2]);
        }
      }

      printf("%lu rays, %d missed\n", rays.size(), nMissed);

      for(int iPnt = 0; iPnt < static_cast<int>(rays.size()); iPnt++) {
        const auto &origin = rays[iPnt].origin;
        const auto &dir = rays[iPnt].dir;

        meshVertexPositions.push_back({origin[0], origin[1], origin[2]});
        meshVertexColors.push_back({0, 1.0, 0});

        facetSDF[rays[iPnt].facetIndex].first += rays[iPnt].dist;
        facetSDF[rays[iPnt].facetIndex].second++; 
      }
    }

    plyOut.addVertexPositions(meshVertexPositions);
    plyOut.addVertexColors(meshVertexColors);


    // Write the object to file
    plyOut.write("raypnts.ply", happly::DataFormat::ASCII);

    if(false){
      happly::PLYData plyOut;

      std::vector<std::array<double, 3>> meshVertexPositions;
      std::vector<std::array<double, 3>> meshVertexColors;

      double maxDist = -1.0;

      for(auto &facet : facetSDF) {

        if(facet.second > 0)
          facet.first /= facet.second;

        maxDist = std::max(maxDist, facet.first);
      }

      printf("max dist %.6f\n", maxDist);

      for(int iPnt = 0; iPnt < static_cast<int>(facetSDF.size()); iPnt++) {
        const auto &origin = pnts[iPnt];

        double weight = facetSDF[iPnt].first/maxDist;

        meshVertexPositions.push_back({origin[0], origin[1], origin[2]});
        meshVertexColors.push_back({weight, 4.0*weight*(1.0-weight), 1.0-weight});

        facetSDF[rays[iPnt].facetIndex].first += rays[iPnt].dist;
        facetSDF[rays[iPnt].facetIndex].second++; 
      }

      plyOut.addVertexPositions(meshVertexPositions);
      plyOut.addVertexColors(meshVertexColors);


      // Write the object to file
      plyOut.write("centroids.ply", happly::DataFormat::ASCII);
    }
}

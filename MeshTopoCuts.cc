#include <vector>
#include "common.h"
#include "MeshTools.h"
#include "tools3d.h"
#include "happly.h"

std::vector<std::vector<MyMesh::HalfedgeHandle>> ComputeTopoCuts(MyMesh &mesh, int nIter)
{
    // computes the lowest value eigenvector of mesh laplacian matrix
    // uses Fiedler vector method, orthogonalizes it with respect to orthValsVec
    // which are eigenvectors previously found by the same method

    std::vector<double> vals;
    std::vector<std::vector<double>> orthVals;

    for(int iTry = 0, nTry = 1; iTry < nTry; iTry++) {
      vals = ComputeEigenvector(mesh, orthVals, nIter*(1+iTry));

      if(iTry < nTry-1)
        orthVals.emplace_back(std::move(vals));
    }

//  auto gradv = ComputeGradient(mesh, vals);

//    vals.swap(gradv);

    double maxVal = vals[0];
    double minVal = vals[1];
    double averVal = 0.0;

    for(const auto &val : vals) {
      maxVal = std::max(val, maxVal);
      minVal = std::min(val, minVal);
      averVal += val;
    }

    averVal /= vals.size();
    printf("aver %.6f, from %.6f to %.6f\n", averVal, minVal, maxVal);

    {
      happly::PLYExport eply;

      for(auto edgeHandle : mesh.edges()) {

        auto vh1 = mesh.to_vertex_handle(mesh.halfedge_handle(edgeHandle, 0));
        auto vh2 = mesh.from_vertex_handle(mesh.halfedge_handle(edgeHandle, 0));

        auto toPnt = toVec(mesh.point(vh1));
        auto fromPnt = toVec(mesh.point(vh2));

        if((toPnt-fromPnt).norm() > 1.5) {

          auto pnts = Sample(Line3d(toPnt, fromPnt), 1.0);

          for(auto &pnt : pnts) {
            eply.AddVertex(pnt, 0.5, 0.1, 0.5);
          }
        }
      }

      eply.WritePLY("longedgesamplepnts.ply");

    }

    // save eigenvector as colored mesh
    {
    // Suppose these hold your data
    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::array<double, 3>> meshVertexColors;
    std::vector<std::vector<size_t>> meshFaceIndices;

    // Create an empty object

    // TODO : request_face_normals
    // TODO : request_vertex_colors, release_vertex_colors, has_vertex_colors

    happly::PLYData plyOut;

    meshVertexPositions.resize(mesh.n_vertices());
    meshVertexColors.resize(mesh.n_vertices());
    meshFaceIndices.resize(mesh.n_faces());

    for (auto v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) {  
        int thisIdx = v_it->idx();

        auto pnt = mesh.point( *v_it );
        meshVertexPositions[thisIdx] = {pnt[0], pnt[1], pnt[2]};

        double weight = (vals[thisIdx]-minVal)/(maxVal-minVal+1e-10);

        meshVertexColors[thisIdx] = {weight, 4.0*weight*(1.0-weight), 1.0-weight};
    }

    for (auto f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it) {  
      int thisIdx = f_it->idx();

      MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it);
      for(; fv_it.is_valid(); ++fv_it) {
        int vId = fv_it->idx();
        meshFaceIndices[thisIdx].push_back(vId);      }
    }

    // Add mesh data (elements are created automatically)
    plyOut.addVertexPositions(meshVertexPositions);
    plyOut.addVertexColors(meshVertexColors);
    plyOut.addFaceIndices(meshFaceIndices);


    // Write the object to file
    plyOut.write("my_output_mesh_file.ply", happly::DataFormat::ASCII);
    }
//    getchar();

//    gradv.swap(vals);

    OpenMesh::HPropHandleT<int> opaqueEdges;
    mesh.add_property(opaqueEdges);

//    meshVertexPositions.clear();
//    meshVertexColors.clear();

    std::vector<std::vector<MyMesh::HalfedgeHandle>> simpleCuts;
    int nSimpleCuts = 0;

    // 1. compute fiedler vector isolines, look for simple cuts formed by isolines
    for(int iVec = 1, nVec = 100; iVec < nVec-1; iVec++) {

      double val = (maxVal*iVec + minVal*(nVec-iVec))/nVec;

      int nLoops = ComputeOpaqueEdges(mesh, vals, val, opaqueEdges);

      printf("%d. isovector %.6f, %d loops\n", iVec, val, nLoops);

      OpenMesh::FPropHandleT<int> pids;
      mesh.add_property(pids);

      for(int iLoop = 0; iLoop < nLoops; iLoop++) {

        int nPids = FloodFill(mesh, opaqueEdges, iLoop, pids);
        if(nPids == 1) {

          printf("simple cut found\n");

          // add loop to simple cuts
          // TODO : check if equivalent to already found cuts

          // TODO : check number of edges, if < 3 then skip
          nSimpleCuts++;

          std::set<MyMesh::VertexHandle> opaqueVertices;
          std::vector<MyMesh::HalfedgeHandle> opaqueEdgeVec;

          double r = 1.0*rand()/RAND_MAX;
          double g = 1.0*rand()/RAND_MAX;
          double b = 1.0*rand()/RAND_MAX;

          int nHedgesInLoop = 0;

          for(auto hedgeHandle : mesh.halfedges()) {

            if(mesh.property(opaqueEdges, hedgeHandle) == iLoop) {
              opaqueEdgeVec.push_back(hedgeHandle);

              nHedgesInLoop++;

              auto toVertex = mesh.to_vertex_handle(hedgeHandle);
              auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

              opaqueVertices.insert(toVertex);
              opaqueVertices.insert(fromVertex);

              auto midPnt = (mesh.point(toVertex) + mesh.point(fromVertex))*0.5;

              //  meshVertexPositions.push_back({midPnt[0], midPnt[1], midPnt[2]});
              //  meshVertexColors.push_back({r, g, b});
            }
          }

          if(opaqueEdgeVec.size() > 2)
            simpleCuts.push_back(opaqueEdgeVec);
        }
      }

      mesh.remove_property(pids);
    }

    // 2. filter out equivalent cuts
    ComputeUniqueCuts(mesh, simpleCuts);

    std::vector<std::vector<MyMesh::HalfedgeHandle>> topoCuts(std::move(simpleCuts));

    printf("%lu topo cuts found\n", topoCuts.size());

    int eulerChar = mesh.n_vertices() + mesh.n_faces() - mesh.n_edges();

    printf("%.6f, %.6f, euler char %d, %d simple cuts, %lu topo cuts found\n", minVal, maxVal, eulerChar, nSimpleCuts, topoCuts.size());

    // 3. compute cross cuts
    std::vector<std::vector<MyMesh::HalfedgeHandle>> crossCuts;

//    SaveLoopsPLY(mesh, topoCuts, "topocuts.ply");

    int iTopoCut = 0;
    for(const auto &topoCut : topoCuts) {

        printf("\n%d(%lu). compute cross cut\n", iTopoCut++, topoCuts.size());

        const auto &opaqueEdgeVec = topoCut;

        std::vector<std::vector<MyMesh::HalfedgeHandle>> thisCrossCuts = ComputeCrossCuts(mesh, opaqueEdgeVec);

        printf("%d. %lu cross cuts found\n", iTopoCut, thisCrossCuts.size());

//        getchar();

        for(auto &crossCut : thisCrossCuts) {
          printf(" %lu edges in cut\n", crossCut.size());
          crossCuts.emplace_back(std::move(crossCut));
        }

//        getchar();
    }

    // TODO : check number of pids 
    ComputeUniqueCuts(mesh, crossCuts);

//    SaveLoopsPLY(mesh, crossCuts, "crosscuts.ply");

    printf("%lu topo cross cuts found\n", crossCuts.size());
//    getchar();

//    ComputeCrossCuts(mesh, crossCuts);

    // add doubleCrossCuts to topoCuts

//    ComputeUniqueCuts(mesh, topoCuts);

    std::vector<std::vector<MyMesh::HalfedgeHandle>> doubleCrossCuts;


    int iLoop = 0;
    for(const auto &opaqueEdgeVec : crossCuts) {      

      printf("\nloop %d of %lu, %lu edges - cross cuts\n", iLoop++, crossCuts.size(), opaqueEdgeVec.size());
//      getchar();

      std::vector<std::vector<MyMesh::HalfedgeHandle>> thisCrossCuts = ComputeCrossCuts(mesh, opaqueEdgeVec);

//      getchar();

      for(auto &thisCrossCut : thisCrossCuts)
       doubleCrossCuts.emplace_back(std::move(thisCrossCut));
    }

    printf("%lu double cross cuts found\n", doubleCrossCuts.size());
//    SaveLoopsPLY(mesh, doubleCrossCuts, "doublecrosscuts.ply");

//    getchar();

// 2. filter out equivalent cuts
    crossCuts.insert(crossCuts.begin(), doubleCrossCuts.begin(), doubleCrossCuts.end());
    ComputeUniqueCuts(mesh, crossCuts);

//    getchar();

    for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
        mesh.property(opaqueEdges, *he_it) == -2;

    printf("%lu topo cuts, %lu cross cuts\n", topoCuts.size(), crossCuts.size());

    topoCuts.insert(topoCuts.end(), crossCuts.begin(), crossCuts.end());

//    SaveLoopsPLY(mesh, topoCuts, "totalcuts.ply");

    ComputeUniqueCuts(mesh, topoCuts);

    SaveLoopsPLY(mesh, topoCuts, "totalcuts.ply");

    printf("%lu unique total cuts\n", topoCuts.size());

    // TODO : extract this into ComputeAndSavePids
    int nOrigPids = 0;
    {     
        OpenMesh::FPropHandleT<int> pids;
        mesh.add_property(pids);

        for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
          mesh.property(opaqueEdges, *he_it) = -2;

        for(size_t iTopoCut = 0; iTopoCut < topoCuts.size(); iTopoCut++) {
            const auto &thisCut = topoCuts[iTopoCut];

            for(auto hedgeHandle : thisCut) {
              mesh.property(opaqueEdges, hedgeHandle) = 0;
              auto oppHedgeHandle = mesh.opposite_halfedge_handle(hedgeHandle);
              mesh.property(opaqueEdges, oppHedgeHandle) = 0;
            }

        }

//        for(size_t iCrossCut = 0; iCrossCut < crossCuts.size(); iCrossCut++) {
//            const auto &thisCut = crossCuts[iCrossCut];

//            for(auto hedgeHandle : thisCut) {
//              mesh.property(opaqueEdges, hedgeHandle) = 0;
//              auto oppHedgeHandle = mesh.opposite_halfedge_handle(hedgeHandle);
//              mesh.property(opaqueEdges, oppHedgeHandle) = 0;
//            }
//        }

        nOrigPids = FloodFill(mesh, opaqueEdges, 0, pids, true);
        printf("whole mesh with all the cuts: %d pids\n", nOrigPids);

        std::vector<std::set<int>> adjPids;

        for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it) {
          auto faceHandle = mesh.face_handle(*he_it);
          auto oppHandle = mesh.opposite_face_handle(*he_it);
 
          auto iPid = mesh.property(pids, faceHandle);
          auto jPid = mesh.property(pids, oppHandle);

          if(iPid != jPid) {
            if(adjPids.size() <= static_cast<size_t>(iPid))
              adjPids.resize(iPid+1);

            adjPids[iPid].insert(jPid);
          }
        }

        for(size_t iPid = 0; iPid < adjPids.size(); iPid++)
          printf("pid %lu : %lu adjacencies\n", iPid, adjPids[iPid].size());

        // generate pid colors
        std::vector<std::array<double, 3>> pidColors;

        for(int iPid = 0; iPid < nOrigPids; iPid++) {
          double r = 1.0*rand()/RAND_MAX;
          double g = 1.0*rand()/RAND_MAX;
          double b = 1.0*rand()/RAND_MAX;

          pidColors.push_back({r, g, b});
        }

        SavePLY(mesh, pids, "pidsfromtopocuts.ply");

        mesh.remove_property(pids);
    }

//    getchar();

    mesh.remove_property(opaqueEdges);


#if 0
    for(size_t iIncludeCut = 0; iIncludeCut < topoCuts.size(); iIncludeCut++)
    {     
        OpenMesh::FPropHandleT<int> pids;
        mesh.add_property(pids);

        const auto &thisCut = topoCuts[iIncludeCut];

        for(auto hedgeHandle : thisCut)
            mesh.property(opaqueEdges, hedgeHandle) = 0;

        auto nPids = FloodFill(mesh, opaqueEdges, 0, pids, true);
        printf("loop %lu : %d pids\n", iIncludeCut, nPids);

        if(nPids > nOrigPids) {
          printf("  skipping\n");
   
          for(auto hedgeHandle : thisCut)
            mesh.property(opaqueEdges, hedgeHandle) = -2;
        }
        else {
          double r = 1.0*rand()/RAND_MAX;
          double g = 1.0*rand()/RAND_MAX;
          double b = 1.0*rand()/RAND_MAX;

          for(auto &hedgeHandle : thisCut) {
            auto toVertex = mesh.to_vertex_handle(hedgeHandle);
            auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

            meshVertexPositions.push_back({fromVertex[0], fromVertex[1], fromVertex[2]});
            meshVertexColors.push_back({r, g, b});

            auto midPnt = (mesh.point(toVertex) + mesh.point(fromVertex))*0.5;

            meshVertexPositions.push_back({midPnt[0], midPnt[1], midPnt[2]});
            meshVertexColors.push_back({r, g, b});
          }
        }

        mesh.remove_property(pids);
    }
#endif

#if 0
    for(size_t iTry = 0; iTry <= topoCuts.size(); iTry++) {

        // mark opaque edges from all cuts except maybe one
        // mark all edges of both cuts as opaque

        for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
          mesh.property(opaqueEdges, *he_it) = -2;

        size_t iSkipCut = topoCuts.size()-iTry;
        for(size_t iTopoCut = 0; iTopoCut < topoCuts.size(); iTopoCut++) {
          if(iTopoCut != iSkipCut) {
            const auto &thisCut = topoCuts[iTopoCut];

            for(auto hedgeHandle : thisCut)
            mesh.property(opaqueEdges, hedgeHandle) = 0;
          }
        }

        OpenMesh::FPropHandleT<int> pids;
        mesh.add_property(pids);

        auto nPids = FloodFill(mesh, opaqueEdges, 0, pids, true);

        mesh.remove_property(pids);

        printf("%lu. skip cut %lu, %d pids\n", iTry, iSkipCut, nPids);
    }
#endif

    return topoCuts;
  };

#include <iostream>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <functional>
#include "happly.h"
#include "AABBTree.h"
#include <eigen3/Eigen/Core>
#include <list>
#include <map>

#include "common.h"

// ----------------------------------------------------------------------------
 
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;
 
int FloodFill(MyMesh &mesh, OpenMesh::HPropHandleT<int> &opaqueEdges, int iLoop, OpenMesh::FPropHandleT<int> &pids, bool shouldDebug)
{
  // flood fill through non-opaque mesh edges, group contiguous facets into pids 
  std::map<std::pair<int, int>, int> edgeCounter;

  for (MyMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it)
    mesh.property(pids, *f_it) = -1;

  int iPid = 0;
  for(MyMesh::FaceIter f_it=mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
    if( mesh.property(pids, *f_it) == -1 ) {

//      printf("seeding %d\n", f_it->idx());
      std::vector<MyMesh::FaceHandle> stack;

      mesh.property(pids, *f_it) = iPid;
      stack.push_back(*f_it);

      int nInPid = 0;
      while(stack.size()) {

        auto facetHandle = stack.back();

        stack.pop_back();
        nInPid++;

        // push neighbors onto stack
        for(MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(facetHandle); fh_it.is_valid(); ++fh_it) {

          if(mesh.property(opaqueEdges, *fh_it) != iLoop) {
            auto oppHedgeHandle = mesh.opposite_halfedge_handle(*fh_it);
            auto oppHandle = mesh.opposite_face_handle(*fh_it);

            if(!mesh.is_boundary(oppHedgeHandle)) {
              if(mesh.property(pids, oppHandle) == -1) {
                mesh.property(pids, oppHandle) = iPid;
                stack.push_back(oppHandle);
              }
              else if (mesh.property(pids, oppHandle) != iPid) {
                printf("!!! %d : %d in wrong pid %d !!!\n", facetHandle.idx(), oppHandle.idx(), mesh.property(pids, oppHandle));
                printf("%d in pid\n", nInPid);
              }
            }
          } 
        }
      }

      if(shouldDebug)
        printf("   pid %d has %d faces\n", iPid, nInPid);

      iPid++;
    }
  }

  return iPid;
}

void ComputeUniqueCuts(MyMesh &mesh, std::vector<std::vector<MyMesh::HalfedgeHandle>> &simpleCuts)
{
    // take input simpleCuts, remove redundant
    printf("computing unique cuts : starting with %lu\n", simpleCuts.size());

    OpenMesh::HPropHandleT<int> opaqueEdges;
    mesh.add_property(opaqueEdges);

    std::vector<std::pair<double, std::vector<MyMesh::HalfedgeHandle>>> topoCuts;

    int nTopoCuts = 0;

    while(simpleCuts.size()) {

      nTopoCuts++;
      // find shortest simple cut
      double minLen = -1.0;
      size_t iMinCut = simpleCuts.size();

      for(size_t iCut = 0; iCut < simpleCuts.size(); iCut++) {
        if(simpleCuts[iCut].size() == 0)
          continue;

        double cutLen = 0.0;

        for(const auto &hedgeHandle : simpleCuts[iCut]) {
          auto toVertex = mesh.to_vertex_handle(hedgeHandle);
          auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

          double edgeLen = (mesh.point(toVertex) - mesh.point(fromVertex)).norm();
          cutLen += edgeLen;
        }

        if(iMinCut == simpleCuts.size() || minLen > cutLen) {
          minLen = cutLen;
          iMinCut = iCut;
        }
      }

      if(iMinCut == simpleCuts.size())
        break;

      printf(" %d. cut has %lu edges, length %.6f\n", nTopoCuts-1, simpleCuts[iMinCut].size(), minLen);

      auto thisCut = std::move(simpleCuts[iMinCut]);

      if(thisCut.size() == 2)
        continue;

#if 0
      double r = 1.0*rand()/RAND_MAX;
      double g = 1.0*rand()/RAND_MAX;
      double b = 1.0*rand()/RAND_MAX;

      for(auto &hedgeHandle : thisCut) {
        auto toVertex = mesh.to_vertex_handle(hedgeHandle);
        auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

        auto midPnt = (mesh.point(toVertex) + mesh.point(fromVertex))*0.5;

        meshVertexPositions.push_back({midPnt[0], midPnt[1], midPnt[2]});
        meshVertexColors.push_back({r, g, b});
      }
#endif

      for(size_t iCut = 0; iCut < simpleCuts.size(); iCut++) {

        if(simpleCuts[iCut].size() == 0)
          continue;

        // mark opaque edges from both cuts
        // mark all edges of both cuts as opaque

        for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
          mesh.property(opaqueEdges, *he_it) = -2;

        for(auto hedgeHandle : thisCut) {
          mesh.property(opaqueEdges, hedgeHandle) = 0;
          mesh.property(opaqueEdges, mesh.opposite_halfedge_handle(hedgeHandle)) = 0;
        }

        for(auto hedgeHandle : simpleCuts[iCut]) {
          mesh.property(opaqueEdges, hedgeHandle) = 0;
          mesh.property(opaqueEdges, mesh.opposite_halfedge_handle(hedgeHandle)) = 0;
        }

        OpenMesh::FPropHandleT<int> pids;
        mesh.add_property(pids);

        // flood fill
        int nPids = FloodFill(mesh, opaqueEdges, 0, pids);

        bool isEquivalent = false;
        if( nPids > 1 )
        // then discard the other cut as it's equivalent to the min cut 
          isEquivalent = true;
        else {
          // check for coincident 

          std::vector<std::set<MyMesh::VertexHandle>> vertices;
          vertices.resize(2);

          for(auto &hedgeHandle : thisCut) {
            auto toVertex = mesh.to_vertex_handle(hedgeHandle);
            auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

            vertices[0].insert(toVertex);
            vertices[0].insert(fromVertex);
          }

          for(auto &hedgeHandle : simpleCuts[iCut]) {
            auto toVertex = mesh.to_vertex_handle(hedgeHandle);
            auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

            vertices[1].insert(toVertex);
            vertices[1].insert(fromVertex);
          }
          
          std::map<MyMesh::VertexHandle, int> vertexLoopCount;

          for(const auto &vertexSet : vertices)
            for(auto &vertexHandle : vertexSet)
              vertexLoopCount[vertexHandle]++;

          size_t nMultiple = 0;
          for(const auto &vertexRec : vertexLoopCount)
            if(vertexRec.second > 1)
              nMultiple++;

          if(nMultiple == vertexLoopCount.size())
            isEquivalent = true;
        }

        if(isEquivalent)
          simpleCuts[iCut].clear();

        mesh.remove_property(pids);
      }

      topoCuts.emplace_back(minLen, std::move(thisCut));
    }

    printf("%lu cuts, %lu unique\n", simpleCuts.size(), topoCuts.size());

    std::vector<std::vector<MyMesh::HalfedgeHandle>> uniqueCuts;

    for(auto &topoCut : topoCuts)
      uniqueCuts.emplace_back(std::move(topoCut.second));

    simpleCuts = std::move(uniqueCuts);

    mesh.remove_property(opaqueEdges);
}

bool IsCutSimple(MyMesh &mesh,
                 const std::vector<MyMesh::HalfedgeHandle> &opaqueEdgeVec)
{
  // check if this cut doesn't separate mesh into two contiguous pieces
  if(opaqueEdgeVec.size() <= 2)
    return false;

  // check if degenerate
  std::set<MyMesh::HalfedgeHandle> edges;

  for(auto hedgeHandle : opaqueEdgeVec) {
    edges.insert(hedgeHandle);
    edges.insert(mesh.opposite_halfedge_handle(hedgeHandle));
  }

  if(2*opaqueEdgeVec.size() - edges.size() >= 6) {
//    printf(" !!! susp loop : %lu loop edges, %lu closure !!!\n", opaqueEdgeVec.size(), edges.size());
    return false;
  }

  OpenMesh::HPropHandleT<int> opaqueEdges;
  mesh.add_property(opaqueEdges);

  for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
    mesh.property(opaqueEdges, *he_it) = -2;

  for(auto hedgeHandle : opaqueEdgeVec) {
    mesh.property(opaqueEdges, hedgeHandle) = 0;
    mesh.property(opaqueEdges, mesh.opposite_halfedge_handle(hedgeHandle)) = 0;
  }

  OpenMesh::FPropHandleT<int> pids;
  mesh.add_property(pids);

  // flood fill
  int nPids = FloodFill(mesh, opaqueEdges, 0, pids);

  mesh.remove_property(pids);
  mesh.remove_property(opaqueEdges);

  if(nPids != 1) {
//      printf(" !!! %lu edges cross cut not simple !!!\n", opaqueEdgeVec.size());
  }

  return nPids == 1;
}

std::vector<std::vector<MyMesh::HalfedgeHandle>> ComputeCrossCuts(MyMesh &mesh, 
                      const std::vector<MyMesh::HalfedgeHandle> &opaqueEdgeVec)
{
    // compute simple cuts "perpendicular" to given cut in opaqueEdgeVec
    bool shouldSave = true;

    OpenMesh::HPropHandleT<int> opaqueEdges;
    mesh.add_property(opaqueEdges);

    // find cross loop
    std::set<MyMesh::VertexHandle> opaqueVertices;

    std::vector<std::pair<double, std::vector<MyMesh::HalfedgeHandle>>> lenCrossCuts;

    // clear opaque edges
    for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it)
      mesh.property(opaqueEdges, *he_it) = -2;
      
    for(auto hedgeHandle : opaqueEdgeVec) {


              // set opaque edges
              mesh.property(opaqueEdges, hedgeHandle) = 0;
              mesh.property(opaqueEdges, mesh.opposite_halfedge_handle(hedgeHandle)) = 0;

              auto toVertex = mesh.to_vertex_handle(hedgeHandle);
              auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

              opaqueVertices.insert(toVertex);
              opaqueVertices.insert(fromVertex);

    }


#if 0

    for(auto hedgeHandle : opaqueEdgeVec) {

          auto toVertex = mesh.to_vertex_handle(hedgeHandle);
          auto fromVertex = mesh.from_vertex_handle(hedgeHandle);


          auto midPnt = (mesh.point(toVertex) + mesh.point(fromVertex))*0.5;

          meshVertexPositions.push_back({midPnt[0], midPnt[1], midPnt[2]});
          meshVertexColors.push_back({0.5, 0.5, 0.5});
    }


    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::array<double, 3>> meshVertexColors;
    std::vector<std::vector<size_t>> meshFaceIndices;

    happly::PLYData plyOut;

    plyOut.addVertexPositions(meshVertexPositions);
    plyOut.addFaceIndices(meshFaceIndices);
    plyOut.addVertexColors(meshVertexColors);

    // Write the object to file
    plyOut.write("thisloop.ply", happly::DataFormat::ASCII);

#endif

    for(size_t iEdge = 0; iEdge < opaqueEdgeVec.size(); iEdge++) {

              auto sideAVertex = mesh.from_vertex_handle(opaqueEdgeVec[iEdge]);
              auto sideBVertex = mesh.to_vertex_handle(opaqueEdgeVec[iEdge]);

              std::vector<MyMesh::VertexHandle> fromVertices;
              // iterate over sideAVertex adjacencies, find fromVertex that's not in opaqueVertices
              for(auto vh_iter = mesh.voh_iter(sideAVertex); vh_iter.is_valid(); ++vh_iter) {
                auto otherVertex = mesh.to_vertex_handle(*vh_iter);

                if(opaqueVertices.find(otherVertex) == opaqueVertices.end())
                  fromVertices.push_back(otherVertex);
              }

              std::vector<MyMesh::VertexHandle> toVertices;
              // iterate over sideBVertex adjacencies, find toVertex that's not in opaqueVertices
              for(auto vh_iter = mesh.voh_iter(sideBVertex); vh_iter.is_valid(); ++vh_iter) {
                auto otherVertex = mesh.to_vertex_handle(*vh_iter);

                if(opaqueVertices.find(otherVertex) == opaqueVertices.end())
                  toVertices.push_back(otherVertex);
              }

              if(toVertices.size() == 0 || fromVertices.size() == 0)
                continue;

              for(auto toVertex : toVertices)
                for(auto fromVertex : fromVertices)
                  if(toVertex != fromVertex) {

                    auto toPnt = mesh.point(toVertex);
                    auto fromPnt = mesh.point(fromVertex);

//                    printf(" from (%.6f, %.6f, %.6f) to (%.6f, %.6f, %.6f)\n", 
//                            fromPnt[0], fromPnt[1], fromPnt[2], toPnt[0], toPnt[1], toPnt[2]);

//                    getchar();

                    auto canWalkEdge = [&mesh, &opaqueEdges, &opaqueVertices, iLoop=0](MyMesh::HalfedgeHandle edgeHandle) {

                      if( mesh.property(opaqueEdges, edgeHandle) == iLoop )
                        return false;

                      // TODO : vertex property
                      auto sideAVertex = mesh.from_vertex_handle(edgeHandle);
                      auto sideBVertex = mesh.to_vertex_handle(edgeHandle);

                      if(opaqueVertices.find(sideAVertex) != opaqueVertices.end())
                        return false;

                      if(opaqueVertices.find(sideBVertex) != opaqueVertices.end())
                        return false;

                      return true;
                    };

                    auto isFinished = [&mesh, toVertex](MyMesh::VertexHandle vertexHandle) {
                      //return false;
                      return vertexHandle == toVertex;
                    };

                    std::map<MyMesh::VertexHandle, std::pair<double, int>> distMap;
                    int nComputed = DijxtraDistances(mesh, {fromVertex}, canWalkEdge, isFinished, distMap);

                    double loopLength = -1.0;

                    if( distMap.find(toVertex) != distMap.end() && 
                        (loopLength < 0.0 || distMap[toVertex].first < loopLength)) {
                      loopLength = distMap[toVertex].first;

                      std::vector<MyMesh::HalfedgeHandle> shortestPath = DijxtraPath(mesh, distMap, toVertex);

                      if(shortestPath.size() > 2) {

                        // shortest path from toVertex to fromVertex without opaqueEdges, to close the loop
                        std::map<MyMesh::VertexHandle, std::pair<double, int>> distMap;
                        int nComputed = DijxtraDistances(mesh, {toVertex}, 
                          [](MyMesh::HalfedgeHandle edgeHandle) {return true;}, 
                          [fromVertex](MyMesh::VertexHandle vertexHandle) { return vertexHandle == fromVertex; },
                          distMap);
        
                        if(nComputed > 0) {
                          std::vector<MyMesh::HalfedgeHandle> backPath = DijxtraPath(mesh, distMap, fromVertex);

//                          printf(" %lu + %lu cross cut path\n", shortestPath.size(), backPath.size());

                          for(const auto &bp: backPath)
                            shortestPath.push_back(bp);
                        }
            
                        // TODO : test if simple cut

                        if(!IsCutSimple(mesh, shortestPath)) {
                          shortestPath.clear();
                        }

//                        printf(" cross loop has %lu edges\n", shortestPath.size());
                      }

                      if(shortestPath.size() > 2) {
                        // test is cuts are equivalent

                        bool initOpaque = true;

                        OpenMesh::HPropHandleT<int> opaqueCrossEdges;
                        mesh.add_property(opaqueCrossEdges);

                        bool isEquivalent = false;
                        for(auto &crossCut : lenCrossCuts) {
                          // test is shortestPath is equivalent to crossCut.second

                          if(initOpaque) {
                            for(auto hedgeHandle : mesh.halfedges())
                              mesh.property(opaqueCrossEdges, hedgeHandle) = -2;

                            initOpaque = true;
                          }

                          for(auto hedgeHandle : shortestPath) {
                            mesh.property(opaqueCrossEdges, hedgeHandle) = 0;
                            mesh.property(opaqueCrossEdges, mesh.opposite_halfedge_handle(hedgeHandle)) = 0;
                          }

                          for(auto hedgeHandle : crossCut.second) {
                            mesh.property(opaqueCrossEdges, hedgeHandle) = 0;
                            mesh.property(opaqueCrossEdges, mesh.opposite_halfedge_handle(hedgeHandle)) = 0;
                          }

                          OpenMesh::FPropHandleT<int> pids;
                          mesh.add_property(pids);

                          auto nPids = FloodFill(mesh, opaqueCrossEdges, 0, pids);
      //                    printf("loop %lu : %d pids\n", iIncludeCut, nPids);

                          mesh.remove_property(pids);

                          for(auto hedgeHandle : crossCut.second) {
                            mesh.property(opaqueCrossEdges, hedgeHandle) = -2;
                            mesh.property(opaqueCrossEdges, mesh.opposite_halfedge_handle(hedgeHandle)) = -2;
                          }

                          if(nPids > 1) {
//                            printf("  equivalent to pre-existing\n");
                            isEquivalent = true;

                            // update pre-existing if shorter
                            if( loopLength < crossCut.first) {
//                              printf(" >> replacing %lu with %lu edges\n", crossCut.second.size(), shortestPath.size());
                              crossCut.first = loopLength;
                              crossCut.second = std::move(shortestPath);
                            }
        
                            break;
                          }

                          mesh.remove_property(pids);
                        }

                        mesh.remove_property(opaqueCrossEdges);

                        if(!isEquivalent) {
//                              printf(" >> adding %lu edges\n", shortestPath.size());
                              lenCrossCuts.emplace_back(loopLength, std::move(shortestPath));
                        }
                      }
                  }
            }
    }

    mesh.remove_property(opaqueEdges);

    std::vector<std::vector<MyMesh::HalfedgeHandle>> crossCuts;

    for(auto &crossCut : lenCrossCuts)
      crossCuts.emplace_back(std::move(crossCut.second));

    printf("%lu cross cuts found\n", crossCuts.size());

    return crossCuts;
}

void AssignPids(MyMesh &mesh, std::map<MyMesh::VertexHandle, std::pair<double, int>> &distMap, OpenMesh::FPropHandleT<int> &pids)
{
  // assign same pid to facets in same voronoi region defined in distMap 
      size_t nAssigned = 0;
      for(auto facetHandle : mesh.faces()) {
        mesh.property(pids, facetHandle) = -1;

        std::map<int, int> clusterCounts;

        // get facet vertices
        for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
          auto &vertexDist = distMap[*fv_it];
          clusterCounts[vertexDist.second]++;
        }

        for(const auto &clusterCount : clusterCounts)
          if(clusterCount.second == 3) {
            mesh.property(pids, facetHandle) = clusterCount.first;
            nAssigned++;
          }
      }    

      printf("vertices : %lu facets, %d assigned\n", mesh.n_faces(), nAssigned);

      while(nAssigned < mesh.n_faces()) {
      // assign pid based on neighbors
      for(auto facetHandle : mesh.faces())
        if(mesh.property(pids, facetHandle) == -1) {

            double bestEdgeLen = -1.0;
            int iBestPid = -1;

            for(MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(facetHandle); fh_it.is_valid(); ++fh_it) {

              // get opposite facet 

              auto oppHandle = mesh.opposite_face_handle(*fh_it);
              if(oppHandle.is_valid()) {

                int oppPid = mesh.property(pids, oppHandle);

                if(oppPid != -1) {
                  // get edge len
                  auto toVertexHandle = mesh.to_vertex_handle(*fh_it);
                  auto fromVertexHandle = mesh.from_vertex_handle(*fh_it);

                  double edgeLen = (mesh.point(toVertexHandle)-mesh.point(fromVertexHandle)).norm();

                  if(edgeLen > bestEdgeLen) {
                    bestEdgeLen = edgeLen;
                    iBestPid = oppPid;
                  }
                }
              }
            }

            if(iBestPid != -1) {
              nAssigned++;
              mesh.property(pids, facetHandle) = iBestPid;
            }
        }

      printf("neighbors : %d facets, %lu assigned\n", mesh.n_faces(), nAssigned);
      }

      for(auto facetHandle : mesh.faces())
        if(mesh.property(pids, facetHandle) == -1) {
          double minDist = -1.0;
          int iBestCluster = -1;
          for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {

            auto &vertexDist = distMap[*fv_it];
            if(minDist < 0.0 || vertexDist.first < minDist) {
              minDist = vertexDist.first;
              iBestCluster = vertexDist.second;
            }
          }

          mesh.property(pids, facetHandle) = iBestCluster;
        }
}

void FlipIsolatedFacets(MyMesh &mesh, OpenMesh::FPropHandleT<int> &pids)
{
  // assign isolated facet to neighboring pid
    int nFlip = -1;

    while(nFlip != 0) {

      nFlip = 0;
      for(auto facetHandle : mesh.faces()) {

        std::map<int, double> pidEdgeLen;

        auto facetPid = mesh.property(pids, facetHandle);
        bool pidContiguous = false;
        // iterate over halfedges
        for(MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(facetHandle); fh_it.is_valid(); ++fh_it) {

            // get opposite facet
            // check if the same pid
            // compute edgeLen

            auto oppHandle = mesh.opposite_face_handle(*fh_it);

            if(oppHandle.is_valid()) {
              auto oppPid = mesh.property(pids, oppHandle);

              if(oppPid == facetPid)
                pidContiguous = true;
              else {
                double edgeLen = 1.0; // TODO : compute edge len
                pidEdgeLen[oppPid] += edgeLen;
              }
            }
        }

        if(!pidContiguous) {

          int bestPid = -1;
          auto maxLen = -1.0;

          for(const auto &edgeLen : pidEdgeLen)
            if(edgeLen.second > maxLen) {
              maxLen = edgeLen.second;
              bestPid = edgeLen.first;
            }

          if(bestPid >= 0) {
            assert(bestPid != facetPid);
            mesh.property(pids, facetHandle) = bestPid;
            nFlip++;
          }
        }
      }

      printf("%d isolated pid facets flipped\n", nFlip);

    } while(nFlip);
}

void SaveLoopsPLY(MyMesh &mesh, const std::vector<std::vector<MyMesh::HalfedgeHandle>> &crossCuts, std::string fileName)
{
    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::array<double, 3>> meshVertexColors;
    std::vector<std::vector<size_t>> meshFaceIndices;

    happly::PLYData plyOut;


    for(const auto &crossCut : crossCuts) {

        double r = 1.0*rand()/RAND_MAX;
        double g = 1.0*rand()/RAND_MAX;
        double b = 1.0*rand()/RAND_MAX;

        for(auto hedgeHandle : crossCut) {

          auto toVertex = mesh.to_vertex_handle(hedgeHandle);
          auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

          auto midPnt = (mesh.point(toVertex) + mesh.point(fromVertex))*0.5;

          meshVertexPositions.push_back({midPnt[0], midPnt[1], midPnt[2]});
          meshVertexColors.push_back({r, g, b});
      }
    }

            plyOut.addVertexPositions(meshVertexPositions);
          plyOut.addFaceIndices(meshFaceIndices);
          plyOut.addVertexColors(meshVertexColors);

  plyOut.write(fileName.c_str(), happly::DataFormat::ASCII);
}


void SavePLY(MyMesh &mesh, const OpenMesh::FPropHandleT<int> &pids, std::string fileName)
{
  happly::PLYExport exportPly;

  std::vector<std::set<size_t>> adjPids;
//  adjPids.resize(maxPid+1);

  for(auto hedgeHandle : mesh.halfedges()) {
    // get two facets
    auto facetHandle = mesh.face_handle(hedgeHandle);
    auto oppHandle = mesh.opposite_face_handle(hedgeHandle);

    if(facetHandle.is_valid() && oppHandle.is_valid()) {
      int iPid = mesh.property(pids, facetHandle);
      int jPid = mesh.property(pids, oppHandle);

      if(iPid < jPid)
        std::swap(iPid, jPid);

      if(iPid >= static_cast<int>(adjPids.size()))
        adjPids.resize(iPid+1);

      // get two pids
      if(iPid != jPid) {
        adjPids[iPid].insert(jPid);
        adjPids[jPid].insert(iPid);
      }
    }
  }

  // TODO : print statistics
  // TODO : get pid contours, color

  auto pidLoops = ComputePatchBoundaries(mesh, pids);

//        std::vector<std::array<double, 3>> pidColors;

//       for(int iPid = 0; iPid < nOrigPids; iPid++) {
//          double r = 1.0*rand()/RAND_MAX;
//          double g = 1.0*rand()/RAND_MAX;
//          double b = 1.0*rand()/RAND_MAX;

//          pidColors.push_back({r, g, b});
//        }



  size_t nVertices = 0;

  for(int iPid = 0, maxPid = 0; iPid <= maxPid; iPid++) {

    double r = (0.3*rand())/RAND_MAX;
    double g = (0.3*rand())/RAND_MAX;
    double b = (0.3*rand())/RAND_MAX;

//    printf("%.6f, %.6f, %.6f\n", r, g, b);
    std::vector<MyMesh::VertexHandle> newSourceVertices;
    bool isNonManifold = IsPatchTopologyInvalid(mesh, pids, pidLoops[iPid], newSourceVertices);

    if(isNonManifold)
      r = std::min(1.0, r+0.8);
    
    if(!isNonManifold || pidLoops[iPid].size() > 1)
      b = std::min(1.0, b+0.8);

    std::map<MyMesh::VertexHandle, size_t> vertexMap;
    std::vector<MyMesh::FaceHandle> pidFacets;

    // add pid vertices
    for(auto facetHandle : mesh.faces()) {
      auto facetPid = mesh.property(pids, facetHandle);

      if(iPid == 0)
        maxPid = std::max(maxPid, facetPid);

      if(facetPid == iPid) {
        pidFacets.push_back(facetHandle);

        // iterate over vertices
        for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
          auto vertexHandle = *fv_it;
          auto insertIter = vertexMap.insert({vertexHandle, nVertices});
          if(insertIter.second) {
            exportPly.AddVertex(mesh.point(vertexHandle), r, g, b);
            nVertices++;
          }
        }
      }
    }

    // add pid facets

    for(const auto &facetHandle : pidFacets) {
      // iterate over vertices

      std::vector<size_t> vertices;
      for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
        vertices.push_back(vertexMap[*fv_it]);
      }

      if(vertices.size() >= 3)
        exportPly.AddFacet(vertices[0], vertices[1], vertices[2]);
    }
  }

  exportPly.WritePLY(fileName);
}

std::vector<std::list<std::vector<MyMesh::HalfedgeHandle>>> ComputePatchBoundaries(MyMesh &mesh, const OpenMesh::FPropHandleT<int> &pids)
{
  // edge chains corresponding to pid boundary

  // edge properties
  OpenMesh::HPropHandleT<int> pidHedges;

  mesh.add_property(pidHedges);
//  printf("ComputePatchBoundaries\n");

  // compute edges bordering on pid boundaries

  int maxPid = 0;

  for(auto hedgeHandle : mesh.halfedges()) {

    mesh.property(pidHedges, hedgeHandle) = -1;

    // get two facets
    auto facetHandle = mesh.face_handle(hedgeHandle);
    auto oppHandle = mesh.opposite_face_handle(hedgeHandle);

    maxPid = std::max(maxPid, mesh.property(pids, facetHandle));
    if (mesh.property(pids, facetHandle) >= 0 && mesh.property(pids, oppHandle) >= 0 && mesh.property(pids, facetHandle) != mesh.property(pids, oppHandle))
      mesh.property(pidHedges, hedgeHandle) = mesh.property(pids, facetHandle);  
  }

  // compute loops

  std::vector<std::list<std::vector<MyMesh::HalfedgeHandle>>> pidLoops;

  if(maxPid >= 0)
    pidLoops.resize(maxPid+1);

//  printf("%d max pid, %lu loops\n", maxPid, pidLoops.size());

  for(auto hedgeHandle : mesh.halfedges()) {

    int thisPid = mesh.property(pidHedges, hedgeHandle);

    if(thisPid >= 0) {

//      printf("%d, %d\n", hedgeHandle.idx(), thisPid);

      std::vector<MyMesh::HalfedgeHandle> loop;

      bool loopFinished = false;

      MyMesh::HalfedgeHandle nextHedge = hedgeHandle;

      // TODO : non-manifold loop vertices
      while(nextHedge.is_valid()) {
//        printf(" %d, %d\n", nextHedge.idx(), thisPid);


        loop.push_back(nextHedge);
        mesh.property(pidHedges, nextHedge) = -1;

        auto toVertex = mesh.to_vertex_handle(nextHedge);
        bool found = false;

        // loop at halfedges coming out of toVertex
        for(auto vh_iter = mesh.voh_iter(toVertex); vh_iter.is_valid(); ++vh_iter) {
            auto adjHedge = *vh_iter;

//            printf("   %d, %d\n", adjHedge.idx(), mesh.property(pidHedges, adjHedge) );
            if( adjHedge.is_valid() && mesh.property(pidHedges, adjHedge) == thisPid) {
              found = true;
              nextHedge = adjHedge;
//              printf("  >>\n");
              break;
            }
          }

        // TODO : check if loop is broken, not closed  
        if(!found)
          nextHedge = MyMesh::HalfedgeHandle();
      }

      if(thisPid >= 0 && pidLoops.size() > thisPid)
        pidLoops[thisPid].push_back(std::move(loop));
      else
        printf("oops\n");
    }

  }

  int nMult = 0;
  int iPid = 0;
  for(const auto &pidContours : pidLoops) {
//    printf("%d. %d loops\n", iPid++, pidContours.size());

    if(pidContours.size() > 1)
      nMult++;
  }
  printf("  %lu pids, %d with multiple contours\n", pidLoops.size(), nMult);

  mesh.remove_property(pidHedges);

  return pidLoops;
//  printf("~ComputePatchBoundaries\n");
}

bool IsPatchTopologyInvalid(MyMesh &mesh, const OpenMesh::FPropHandleT<int> &pids, 
                            const std::list<std::vector<MyMesh::HalfedgeHandle>> &pidBoundary, 
                            std::vector<MyMesh::VertexHandle> &newSourceVertices)
{
  // patch topology validity checker, must be contiguous, can't have two pieces of boundary with same neighboring pid

    bool invalidTopology = false;
    newSourceVertices.clear();

    if(pidBoundary.size() > 1) {
          printf("!!! %lu contours !!!\n", pidBoundary.size());

          for(const auto &pidLoop : pidBoundary) {
            auto he = pidLoop.front();

            newSourceVertices.push_back(mesh.from_vertex_handle(he));
          }
    }
    else if(pidBoundary.size() == 1) {
          std::vector<int> pidEdges;

          const auto &pidLoop = pidBoundary.front();          
          for(auto hedge : pidLoop) {

            auto oppFacet = mesh.opposite_face_handle(hedge);
            auto oppPid = mesh.property(pids, oppFacet);

            if(pidEdges.size() == 0 || pidEdges.back() != oppPid)
              pidEdges.push_back(oppPid);
          }

          if(pidEdges.size() >= 2 && pidEdges.front() == pidEdges.back())
            pidEdges.pop_back();

          if(pidEdges.size() <= 2)
            printf("!!! %lu total edges !!!\n", pidEdges.size());

          std::map<int, int> pidEdgeCount;

          for(auto iPid : pidEdges)
            pidEdgeCount[iPid]++;

          if(pidEdgeCount.size() < pidEdges.size()) {
            printf("!!! %lu double edges !!\n", pidEdges.size() - pidEdgeCount.size());

            for(auto pidCount : pidEdgeCount)
              if(pidCount.second > 1) {

                bool segStarted = false;

                int segCount = 0;

                for(auto hedge : pidLoop) {

                  auto oppFacet = mesh.opposite_face_handle(hedge);
                  auto oppPid = mesh.property(pids, oppFacet);

                  if(!segStarted && pidCount.first == oppPid) {
                    segStarted = true;
                    segCount++;
                  }
                  if(segStarted && pidCount.first != oppPid) {
                    segStarted = false;

                    newSourceVertices.push_back(mesh.from_vertex_handle(hedge));
                  }
                }
              }
          }
  }

  return invalidTopology || newSourceVertices.size() > 0;
}

#if 0
BuildBorderLoops
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

std::vector<std::vector<MyMesh::Point>> extract_border_contours(MyMesh& mesh) {
    // Ensure edge and halfedge status is available
    mesh.request_edge_status();
    mesh.request_halfedge_status();
    mesh.request_vertex_status();

    std::unordered_set<MyMesh::EdgeHandle> visited_edges;

    std::vector<std::vector<MyMesh::Point>> contours;

    for (auto eh : mesh.edges()) {
        if (mesh.is_boundary(eh) && visited_edges.find(eh) == visited_edges.end()) {
            // Start a new contour
            std::vector<MyMesh::Point> contour;

            // Start from one of the halfedges
            MyMesh::HalfedgeHandle heh = mesh.halfedge_handle(eh, 0);
            if (!mesh.is_boundary(heh)) {
                heh = mesh.opposite_halfedge_handle(heh);
            }

            // Traverse the boundary loop
            MyMesh::HalfedgeHandle start_heh = heh;
            do {
                visited_edges.insert(mesh.edge_handle(heh));
                contour.push_back(mesh.point(mesh.from_vertex_handle(heh)));

                // Find next boundary halfedge
                heh = mesh.next_halfedge_handle(heh);
                while (!mesh.is_boundary(heh)) {
                    heh = mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(heh));
                }
            } while (heh != start_heh);

            // Add last point to close the loop
            contour.push_back(mesh.point(mesh.from_vertex_handle(heh)));

            contours.push_back(contour);
        }
    }

    // Cleanup
    mesh.release_edge_status();
    mesh.release_halfedge_status();
    mesh.release_vertex_status();

    return contours;
}
#endif
#include "common.h"

#include <deque>
#include <cassert>
#include "MeshTools.h"
#include "happly.h"
#include "AABBTree.h"
#include "MeshSampling.h"
#include "ControlPrintf.h"
#include "MeshSampling.h"

// save as points into ply file, sample each edge as 11 points
void MeshSamplingSave(std::vector<Eigen::Vector3d> &meshPoints, std::vector<std::set<size_t>> &adjMatrix, std::string fileName) 
{
  happly::PLYExport exportPly;

  for (size_t iPnt = 0; iPnt < adjMatrix.size(); iPnt++) {

    printf("%lu. %lu adjacencies\n", iPnt, adjMatrix[iPnt].size());

    for (auto jPnt : adjMatrix[iPnt]) {

      int nPnts = 11;

      for (int jj = 1; jj < nPnts; jj += 2)
      {
        double param = 1.0 * jj / nPnts;
        auto midPnt = meshPoints[iPnt] * param + meshPoints[jPnt] * (1.0 - param);

        exportPly.AddVertex({midPnt[0], midPnt[1], midPnt[2]}, 0.8, 0.8, 0.5);
      }
    }
  }

  printf("write %s\n", fileName.c_str());
  exportPly.WritePLY(fileName);
}

// generate mesh sampling by flood filling mesh 
// popping a point from the stack and pushing 4 points formed as point + step*crossDir[ii : {0,3}}]
// can specify points
std::pair<std::vector<MyMesh::VertexHandle>, std::vector<Eigen::Vector3d>> MeshSamplingSampleUV(MyMesh &mesh, double step, int nIter, std::vector<Eigen::Vector3d> &initPoints, std::vector<Eigen::Vector3d> &crossField, std::vector<std::set<size_t>> &adjMatrix)
{
    AABBTree tree;
    tree.Build(mesh, 20);

    happly::PLYExport exportPly;

    std::vector<MyMesh::VertexHandle> sourceVertices;
    std::vector<Eigen::Vector3d> meshPoints;

    std::vector<std::pair<Eigen::Vector3d, bool>> stack;
    std::vector<std::pair<int, int>> stackEdges;
    std::vector<int> vertexIds;

    Eigen::Vector3d sourcePnt;
    
    if(initPoints.size()) {
      sourcePnt = initPoints.front();

      for(auto initPoint : initPoints)
        stack.push_back({initPoint, true});
    }
    else {
      sourcePnt = toVec(mesh.point(MyMesh::VertexHandle(0)));
      stack.push_back({sourcePnt, true});
    }

    vertexIds.resize(stack.size(),-1);

    int iIter = 0;

    while(stack.size()) {

      iIter++;

//      printf("%d. %d elements in stack\n", iIter, stack.size());

//      if(iIter >= 100)
//        break;

      // 1. find next point
      double maxDist = -1.0;

      size_t iBestPnt = std::numeric_limits<size_t>::max();
      bool oneFound = false;

      for(size_t iPnt = 0; iPnt < stack.size(); iPnt++) {
        const auto &pnt = stack[iPnt];
        if(pnt.second) {
          double dist = (pnt.first-sourcePnt).norm();

          if(maxDist < 0.0 || dist < maxDist) {
            maxDist = dist;
            iBestPnt = iPnt;
            oneFound = true;
          }
        }
      }

      if(!oneFound)
        break;

      // 2. TODO :  pop point from stack that's the farmost from source
      stack[iBestPnt].second = false;
      Eigen::Vector3d currPoint = stack[iBestPnt].first;

//      printf("curr pnt (%.6f, %.6f, %.6f)\n", currPoint[0], currPoint[1], currPoint[2]);

      auto np = tree.FindNearestPoint(currPoint);

//      printf("proj pnt (%.6f, %.6f, %.6f), dist %.6f\n", np.pnt[0], np.pnt[1], np.pnt[2], np.dist);

      MyMesh::VertexHandle closestVertex;
      if(np.facetIndex != -1) {
        MyMesh::FaceHandle closestFacetHandle(np.facetIndex);

        // TODO : proper cross field/normal barycentric interpolation
        double closestDist = -1.0;

        for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(closestFacetHandle); fv_it.is_valid(); ++fv_it) {
          MyMesh::VertexHandle facetVertex(fv_it->idx());

          double dist = (toVec(mesh.point(facetVertex))-currPoint).norm();

          if(closestDist < 0.0 || dist < closestDist) {
            closestDist = dist;
            closestVertex = facetVertex;
          }
        }
      }

      // 2. get cross field at the point
      std::vector<Eigen::Vector3d> uvPoints;

      if(closestVertex.is_valid()) {

        auto pnt = mesh.point(closestVertex);

//        printf("vertex pnt (%.6f, %.6f, %.6f)\n", pnt[0], pnt[1], pnt[2]);

        vertexIds[iBestPnt] = meshPoints.size();

        meshPoints.push_back(np.pnt);
        sourceVertices.push_back(closestVertex);

        auto crossVec = crossField[closestVertex.idx()];
        auto normal = toVec(mesh.normal(closestVertex));

        exportPly.AddVertex(np.pnt, 1.0, 1.0, 0.0);

        std::vector<MyMesh::HalfedgeHandle> neighborhood = MeshSamplingNeighborhood(mesh, closestVertex, 1.5*step);

//        printf("curr point %lu (%.6f, %.6f, %.6f)\n", iBestPnt, currPoint[0], currPoint[1], currPoint[2]);

        // make a step at each +-u/v direction
        for(int i = 0; i < 4; i++) {

          Eigen::Vector3d uvPnt = MeshSamplingNeighbor(mesh, currPoint, step, normal, crossVec, neighborhood);

//          printf("  %d.  (%.6f, %.6f, %.6f) at %.6f\n", i, uvPnt[0], uvPnt[1], uvPnt[2], (currPoint - uvPnt).norm());

//          uvPoints.push_back(uvPnt);
          auto np = tree.FindNearestPoint(currPoint + crossVec*step);

//          printf("  neighborX  (%.6f, %.6f, %.6f) at %.6f\n", uvPnt[0], uvPnt[1], uvPnt[2], (currPoint - uvPnt).norm());

//          printf("  projection (%.6f, %.6f, %.6f) at %.6f\n", np.pnt[0], np.pnt[1], np.pnt[2], (currPoint - np.pnt).norm());

          uvPoints.push_back(np.pnt);
//          uvPoints.push_back(uvPnt);

          crossVec = cross(crossVec, normal);
        }
      }

      int iPnt = 0;
      for(auto uvPoint : uvPoints) {

        printf("  %d.  (%.6f, %.6f, %.6f) at %.6f\n", iPnt++, uvPoint[0], uvPoint[1], uvPoint[2], (currPoint - uvPoint).norm());

        bool isClose = false;
        size_t iClosestPnt;
        double minDist = -1.0;

        // 5. put in the cluster with tol = 0.2*step
        for(size_t ii = 0; ii < stack.size(); ii++) 
          if( ii != iBestPnt ) {
            auto clusterPnt = stack[ii];
            double dist = (uvPoint-clusterPnt.first).norm();

            if(dist < 0.6*step) {
              isClose = true;
              if(minDist == -1.0 || minDist > dist) {
                iClosestPnt = ii;
                minDist = dist;
              }
            }
          }

        // 6. if new cluster then stack.push_back(crossPoint[i]);

        if(!isClose) {
          printf("   push at %lu\n", stack.size());
          exportPly.AddVertex(uvPoint, 0.5, 0.8, 0.5);

          stack.push_back({uvPoint, true});
          vertexIds.push_back(-1);
        }
        else {
          printf("   %.6f to %lu\n", minDist, iClosestPnt);
          stackEdges.push_back({iBestPnt, iClosestPnt});
        }
      }
    }

    printf("%d iters, %lu points, %lu edges\n", iIter, stack.size(), stackEdges.size());

///    exportPly.WritePLY("uvinitpoints.ply");


///    std::vector<std::set<size_t>> adjMatrix;

    printf("%lu %lu vertices\n", meshPoints.size(), sourceVertices.size());
    adjMatrix.resize(meshPoints.size());

    for(const auto &edge : stackEdges) {
      printf("edge %d-%d\n", vertexIds[edge.first], vertexIds[edge.second] );
      printf("edge %d-%d\n", vertexIds[edge.second], vertexIds[edge.first]);


      if(vertexIds[edge.first] >= 0 && vertexIds[edge.second] >= 0) {
        adjMatrix[vertexIds[edge.first]].insert(vertexIds[edge.second]);
        adjMatrix[vertexIds[edge.second]].insert(vertexIds[edge.first]);
      }
    }

    for(size_t iPnt = 0; iPnt < adjMatrix.size(); iPnt++) {

      printf("%lu. %lu adjacencies\n", iPnt, adjMatrix[iPnt].size());

      for(auto jPnt : adjMatrix[iPnt]) {

        int nPnts = 11;

        for(int jj = 1; jj < nPnts; jj+=2) {
            double param = 1.0*jj/nPnts;
            auto midPnt = meshPoints[iPnt]*param + meshPoints[jPnt]*(1.0-param);

            exportPly.AddVertex({midPnt[0], midPnt[1], midPnt[2]}, 0.8, 0.8, 0.5);
        }
      }
    }

    exportPly.WritePLY("uvpoints.ply");

    return {sourceVertices, meshPoints};
}

// chain of half edges bounding on rad wide neighborhood of sourceVertex
std::vector<MyMesh::HalfedgeHandle> MeshSamplingNeighborhood(MyMesh &mesh, MyMesh::VertexHandle &sourceVertex, double rad)
{
    FibonacciHeap<sVertexDistRec> prioRTQ;
 
    std::map<MyMesh::VertexHandle, std::pair<double, int>> distMap;

    int nComputed = DijxtraDistances(mesh, {sourceVertex},
                                    [](MyMesh::HalfedgeHandle edgeHandle) {return true;}, // can walk every edge
                                    [&distMap, rad](MyMesh::VertexHandle vertexHandle) {return distMap[vertexHandle].first > rad;}, // finished when all vertices with dist <= rad are extracted
                                    distMap, &prioRTQ);

    //  get all edges with both vertices in distMap
    
    std::vector<MyMesh::HalfedgeHandle> neighborhoodEdges;

    for(auto vmap : distMap) {

      auto currVertex = vmap.first;

      for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(currVertex); voh_it.is_valid(); ++voh_it) {
          auto testHedge = *voh_it;

          auto toVertex = mesh.to_vertex_handle(testHedge);

          if(distMap.find(toVertex) != distMap.end())
            neighborhoodEdges.push_back(testHedge);

        }
    }

    return neighborhoodEdges;
}

// make a step of length step from sourcePoint 
Eigen::Vector3d MeshSamplingNeighbor(MyMesh &mesh, Eigen::Vector3d sourcePnt, double rad, Eigen::Vector3d normal, Eigen::Vector3d dir, const std::vector<MyMesh::HalfedgeHandle> &neighborhood)
{
//  Eigen::Vector3d pnt = toVec(mesh.point(sourceVertex));

//  printf(" (%.6f, %.6f, %.6f) rad %.6f, %.6f normal, %.6f dir, %lu neighbors\n", sourcePnt[0], sourcePnt[1], sourcePnt[2], rad, normal.norm(), dir.norm(), neighborhood.size());


  // build a plane passing through point containing norm and dir 
  Eigen::Vector3d planeNormal = cross(normal, dir);
  Plane3d dirPlane(sourcePnt, planeNormal);

  double lowDist = -1.0, highDist = -1.0;
  Eigen::Vector3d lowPoint, highPoint;

  happly::PLYExport eply;

  for(auto hedgeHandle : neighborhood) {

    auto toPnt = toVec(mesh.point(mesh.to_vertex_handle(hedgeHandle)));
    auto fromPnt = toVec(mesh.point(mesh.from_vertex_handle(hedgeHandle)));

    eply.AddVertex(toPnt, 0.0, 0.5, 0.5);
    eply.AddVertex(fromPnt, 0.0, 0.5, 0.5);

    // intersect all edges with the plane
    Line3d edgeSeg(toPnt, fromPnt);

    Eigen::Vector3d xPoint;
    auto extype = dirPlane.Intersect(edgeSeg, xPoint);

    // find xpoints with dists to source points : max of those smaller than rad, min of those >= than rad

    if(extype == Plane3d::Point) { // TODO : handle Plane3d::Line
      double param = edgeSeg.ProjectParam(xPoint);

      double dirDist = dir.dot(xPoint-sourcePnt);
      if(dirDist > 0.0 && param >= 0.0 && param <= 1.0) {

        double sourceDist = (xPoint-sourcePnt).norm();

        double param = 0.5*sourceDist/rad;
        param = std::min(1.0, param);
        param = std::max(0.0, param);

        eply.AddVertex(xPoint, 1.0-param, 4.0*param*(1.0-param), param);

        if(sourceDist < rad) {
          if(lowDist < 0.0 || lowDist < sourceDist) {
            lowDist = sourceDist;
            lowPoint = xPoint;
          }
        }
        else if (sourceDist >= rad) {
          if(highDist < 0.0 || highDist > sourceDist) {
            highDist = sourceDist;
            highPoint = xPoint;
          }
        }
      }
    }
  }

  // interpolate and output the point

  Eigen::Vector3d npnt;

  if(lowDist == -1 && highDist == -1) {
    printf("!!! empty neighborhood !!!\n");
  }
  else if (lowDist == -1.0)
    npnt = highPoint;
  else if (highDist == -1.0)
    npnt = lowPoint;
  else {
    double paramDist = highDist - lowDist;

    if(fabs(paramDist) < eps)
      npnt = highPoint;
    else 
      npnt = lowPoint * (highDist-rad)/paramDist + highPoint * (rad-lowDist)/paramDist;
  }

  eply.AddVertex(npnt, 1.0, 0.0, 0.0);

//  eply.WritePLY("neighborx.ply");

  return npnt;
}

// sample mesh randomly placing vertices at step distance
std::vector<MyMesh::VertexHandle> MeshSamplingSampleRandom(MyMesh &mesh, double step)
{
  // random sampling
    std::vector<MyMesh::VertexHandle> sourceVertices;

    int vCount = mesh.n_vertices()/20+1;

    for(auto vertexHandle : mesh.vertices()) {

      auto pnt = mesh.point(vertexHandle);

      bool isClose = false;

      for(auto sourceVertex : sourceVertices) {
        auto sourcePnt = mesh.point(sourceVertex);

        double distpnt = (sourcePnt-pnt).norm();

        if(distpnt < step) {
          isClose = true;
          break;
        }
      }

      if(!isClose)
        sourceVertices.push_back(vertexHandle);
    }

//    sourceVertices.resize(4);
    printf("%lu vertices, %lu sources\n", mesh.n_vertices(), sourceVertices.size());

    return sourceVertices;
}

// compute centroidal voronoi using Lloyd clustering
std::pair<std::vector<MyMesh::VertexHandle>, std::vector<int>> MeshSamplingCluster(MyMesh &mesh, 
  std::vector<MyMesh::VertexHandle> &sourceVertices,  int nIter, std::vector<Eigen::Vector3d> &crossField)
{

//    OpenMesh::VPropHandleT<Eigen::Vector3d> dirField;
//    mesh.add_property(dirField);
//    std::vector<Eigen::Vector3d> crossField = ComputeCrossfield(mesh, dirField, 100);
//    mesh.remove_property(dirField);


    auto SaveSourcePoints = [&mesh](const std::vector<MyMesh::VertexHandle> &sourceVertices, std::string fileName) {

      happly::PLYExport exportPly;

//    std::vector<Eigen::Vector3d> points;
//    std::vector<Eigen::Vector3d> normals;
//    std::vector<Eigen::Vector3d> dirs;

      for(auto vertex : sourceVertices) {
        auto pnt = mesh.point(vertex);
    
 //       auto normal = mesh.normal(vertex);

//      points.push_back(toVec(pnt));
//      normals.push_back(toVec(normal));
//      dirs.push_back(crossField[vertex.idx()]);

        exportPly.AddVertex({pnt[0], pnt[1], pnt[2]}, 0.0, 1.0, 0.0);
      }

      exportPly.WritePLY(fileName);
    };

    SaveSourcePoints(sourceVertices, "initvoro.ply");

    AABBTree tree;
    tree.Build(mesh, 20);

    std::map<MyMesh::VertexHandle, std::pair<double, int>> distMap;

    auto SaveDistMap = [&mesh, &distMap, &sourceVertices](std::string voroDistMap) {
    happly::PLYExport exportPly;
    
    std::vector<std::array<double, 3>> colors(sourceVertices.size());

    for(size_t iPnt = 0; iPnt < sourceVertices.size(); iPnt++) {

      double r = 1.0*rand()/RAND_MAX;
      double g = 1.0*rand()/RAND_MAX;
      double b = 1.0*rand()/RAND_MAX; 

      double sum = r+g+b;

      if(sum > 0.0) {
        r /= sum; 
        g /= sum; 
        b /= sum;
      }

      colors[iPnt] = {r, g, b};
    }

    double maxDist = 0.0;
    for(const auto &dist : distMap)
      maxDist = std::max(maxDist, dist.second.first);

    for(auto vertexHandle : mesh.vertices()) {
      double thisDist = distMap[vertexHandle].first;
      int iPid = distMap[vertexHandle].second;

      double weight = (maxDist+thisDist)/(2.0*maxDist);

      exportPly.AddVertex(mesh.point(vertexHandle), colors[iPid][0]*weight, colors[iPid][1]*weight, colors[iPid][2]*weight);
    }

    for(const auto &facetHandle : mesh.faces()) {
      // iterate over vertices

      std::vector<size_t> vertices;
      for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
        vertices.push_back(fv_it->idx());
      }

      if(vertices.size() >= 3)
        exportPly.AddFacet(vertices[0], vertices[1], vertices[2]);
    }

    exportPly.WritePLY(voroDistMap);
    };


    OpenMesh::FPropHandleT<int> pids;
    mesh.add_property(pids);

    FibonacciHeap<sVertexDistRec> prioRTQ;
 
    // compute voronoi cells for sourceVertices, compute centroids, move sourceVertices to centroids
    for(int iIter = 0; iIter <= nIter; iIter++) {
      printf("\n%d. \n", iIter);

      distMap.clear();

      while(!prioRTQ.isEmpty())
        prioRTQ.removeMinimum();

      // compute distances

      std::vector<std::array<Eigen::Vector3d, 3>> localCoords;

      int nComputed = DijxtraDistances(mesh, sourceVertices, localCoords,
                                      [](MyMesh::HalfedgeHandle edgeHandle) {return true;}, // can walk every edge
                                      [](MyMesh::VertexHandle vertexHandle) {return false;}, // not finished till all vertices are extracted
                                      distMap, &prioRTQ);

      printf("  %lu sources, %d computed\n", sourceVertices.size(), nComputed);

      // compute patches

      // assign pid based on 2 or more vertices from same source
      AssignPids(mesh, distMap, pids);

      // check for isolated facets

      FlipIsolatedFacets(mesh, pids);

      // compute patch properties, like contiguity, contours
      auto pidBoundaries = ComputePatchBoundaries(mesh, pids);

//      SavePLY(mesh, pids, "voropids.ply");

      std::vector<MyMesh::VertexHandle> newSourceVertices;

      // compute edges, find broken topology

      if(nIter != 0) {
        for(const auto &pidBoundary : pidBoundaries) {
          std::vector<MyMesh::VertexHandle> thisSourceVertices;

          if( IsPatchTopologyInvalid(mesh, pids, pidBoundary, thisSourceVertices) )
            newSourceVertices.insert(newSourceVertices.begin(), thisSourceVertices.begin(), thisSourceVertices.end());
        }
      }

      if(newSourceVertices.size()) {
        printf("%lu new sources\n", newSourceVertices.size());

        {
          happly::PLYExport exportPly;

          for(auto vertex : newSourceVertices) {

            auto pnt = mesh.point(vertex);
            exportPly.AddVertex({pnt[0], pnt[1], pnt[2]}, 0.0, 1.0, 0.0);
          }

          exportPly.WritePLY("newsourcevoro.ply");
//          getchar();

        }
      }

      std::vector<int> clusterCounts;
      clusterCounts.resize(sourceVertices.size());

      for(auto vertexHandle : mesh.vertices()) {
        clusterCounts[distMap[vertexHandle].second]++;
      }

      // compute centroids
      std::vector<std::pair<Eigen::Vector3d, double>> clusterCentroids;
      clusterCentroids.resize(sourceVertices.size());

      double totalShift = 0.0;

        for(auto facetHandle : mesh.faces()) {
          Eigen::Vector3d centroid(0.0, 0.0, 0.0);

          std::vector<int> clusterIndices;
          std::vector<Eigen::Vector3d> facePoints;

          for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
              auto pnt = mesh.point(*fv_it);
              facePoints.push_back({pnt[0], pnt[1], pnt[2]});

              centroid += Eigen::Vector3d(pnt[0], pnt[1], pnt[2]);
              clusterIndices.push_back(distMap[*fv_it].second);
          }

          if(clusterIndices.size()) {
            // compute area from facePoints

            double facetArea = 0.0;
            for(size_t iTr = 0; iTr < facePoints.size()-2; iTr++) {
              Triangle3d triangle(facePoints[0], facePoints[iTr+1], facePoints[iTr+2]);

              facetArea += triangle.GetArea();
            }

            assert(facetArea > 0.0);

            centroid *= facetArea/clusterIndices.size();

            for(auto iCluster : clusterIndices) {
              clusterCentroids[iCluster].first += centroid;
              clusterCentroids[iCluster].second += facetArea;
            }
          }
        }


      // TODO : take multiple boundaries, put new points on multiple boundaries

      // if region is a cylinder, not a disk, then put more sources        
      // assign sources to cluster centroids

      int nSame = 0;
      for(size_t iC = 0; iC < clusterCentroids.size(); iC++) {
        
        auto &centroid = clusterCentroids[iC];
//        printf("   cluster %lu, source %d, area %.6f, count %d\n", iC, sourceVertices[iC].idx(), centroid.second, clusterCounts[iC]);
        if(centroid.second > 0.0) {
          centroid.first *= 1.0/centroid.second;

          // project centroid.first onto the mesh
          auto np = tree.FindNearestPoint(centroid.first);

          auto sourcePnt = mesh.point(sourceVertices[iC]);

//          printf("(%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f), (%.6f, %.6f, %.6f), dist %.6f\n",
//                  sourcePnt[0], sourcePnt[1], sourcePnt[2],
//                  centroid.first[0], centroid.first[1], centroid.first[2],
//                  np.pnt[0], np.pnt[1], np.pnt[2], np.dist);

          if(np.facetIndex != -1) {
            MyMesh::FaceHandle closestFacetHandle(np.facetIndex);

            MyMesh::VertexHandle closestVertexHandle;
            double closestDist = -1.0;

            for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(closestFacetHandle); fv_it.is_valid(); ++fv_it) {
              const auto &pnt = mesh.point(*fv_it);

              double currDist = (centroid.first - toVec(pnt)).norm();

              if(closestDist < 0.0 || currDist < closestDist) {
                closestVertexHandle = *fv_it;
                closestDist = currDist;
              }
            }

            if(closestVertexHandle.idx() != -1) {

              auto pnt1 = mesh.point(sourceVertices[iC]);
              auto pnt2 = mesh.point(closestVertexHandle);

              auto ddist = (pnt1-pnt2).dot(pnt1-pnt2);

              totalShift += ddist;
              if(sourceVertices[iC] == closestVertexHandle)
                nSame++;
              else
                sourceVertices[iC] = closestVertexHandle;
            }
          }
        }
      }


//      if(iIter == nIter-1) // nIter)
//        sourceVertices.insert(sourceVertices.end(), newSourceVertices.begin(), newSourceVertices.end());
      
      // compute duplicate sources
      std::set<MyMesh::VertexHandle> usedVertices;
      std::vector<size_t> dupVertices;

      happly::PLYExport dupExportPly;
      for(size_t iS = 0; iS < sourceVertices.size(); iS++) {
        auto insertIterr = usedVertices.insert(sourceVertices[iS]);
        auto pnt = mesh.point(sourceVertices[iS]);

        if(!insertIterr.second) {

          printf("  duplicate at (%.6f, %.6f, %.6f)\n", pnt[0], pnt[1], pnt[2]);

          dupVertices.push_back(iS);

          dupExportPly.AddVertex(pnt, 1.0, 0.0, 0.0);
        }
        else 
          dupExportPly.AddVertex(pnt, 0.0, 1.0, 0.0);
      }
      dupExportPly.WritePLY("dupsources.ply");
//      getchar();

      double maxDist = 0.0;

      for(const auto &distVert : distMap)
        maxDist = std::max(maxDist, distVert.second.first);

      printf("%d. %d placed, max dist %.6f, %d same, aver shift %.6f, %lu duplicate\n",
             iIter, nComputed, maxDist, nSame, nComputed ? totalShift/nComputed : 0.0, dupVertices.size());

      if(iIter == nIter)
        break;

      for(size_t iDup = 0; iDup < dupVertices.size(); iDup++) {

        double maxDist = -1.0;
        bool updated = false;

        for(const auto &dm : distMap) {

          double thisDist = dm.second.first;

          if (thisDist > maxDist) {
            auto fromPnt = mesh.point(dm.first);

            // stay away from vertices that have been just assigned
            for(auto usedVertex : usedVertices) {
              auto toPnt = mesh.point(usedVertex);
              double dist = (toPnt-fromPnt).norm();

              thisDist = std::min(dist, thisDist);
            }
          }

          if(thisDist > maxDist) {
            maxDist = thisDist;
            sourceVertices[dupVertices[iDup]] = dm.first;
            updated = true;
          }
        }

        if(updated) {
          usedVertices.insert(sourceVertices[dupVertices[iDup]]);
          auto pnt = mesh.point(sourceVertices[dupVertices[iDup]]);

          printf("  relocated to (%.6f, %.6f, %.6f), %lu relocated\n", pnt[0], pnt[1], pnt[2], usedVertices.size());
        }
      }
    }
    // TODO : print pids statistics

    SavePLY(mesh, pids, "voropids.ply");

    mesh.remove_property(pids);

    happly::PLYExport exportPly;

    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Eigen::Vector3d> dirs;

    for(auto vertex : sourceVertices) {
      auto pnt = mesh.point(vertex);
    
      auto normal = mesh.normal(vertex);

      points.push_back(toVec(pnt));
      normals.push_back(toVec(normal));
      dirs.push_back(crossField[vertex.idx()]);

      exportPly.AddVertex({pnt[0], pnt[1], pnt[2]}, 0.0, 1.0, 0.0);
    }

  // save source points crossfield at sample points only
  for(auto vertexHandle : sourceVertices) {
    auto vnormal = mesh.normal(vertexHandle);
    auto pnt = mesh.point(vertexHandle);
 
//    printf("   %.6f, %.6f, %.6f\n", vnormal[0], vnormal[1], vnormal[2]);

    Eigen::Vector3d vpoint(pnt[0], pnt[1], pnt[2]);

    auto cf = crossField[vertexHandle.idx()];

    Eigen::Matrix3d mv = Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d(vnormal[0], vnormal[1], vnormal[2])).toRotationMatrix();

    std::vector<Eigen::Vector3d> vertexRosy;

    vertexRosy.push_back(crossField[vertexHandle.idx()]);
    for(int iR = 1; iR < 4; iR++) {
      vertexRosy.push_back(mv*vertexRosy[iR-1]);
    }

    double weight = 0.5;//dotProds[vertexHandle.idx()];
//    weight = 1.0-5.0*(1.0-weight);

    double red = 1.0-weight;
    double blue = weight;
    double green = 4.0*red*blue;

    exportPly.AddVertex({vpoint[0], vpoint[1], vpoint[2]}, 0.5*red, 0.5*green, 0.5*blue);

    for(auto vr : vertexRosy) {
      Eigen::Vector3d vvr = vpoint + 0.1*vr;

      exportPly.AddVertex({vvr[0], vvr[1], vvr[2]}, red, green, blue);
    }
  }

  exportPly.WritePLY("centrovoro.ply");

  // save vertices with distance map showing how far away from the source are they

  SaveDistMap("vorodistmap.ply");

  // for each vertex, output the source vertex which voronoi cell does this vertex belong to
  std::vector<int> regionMap(mesh.n_vertices());

  for(auto vertexHandle : mesh.vertices())
    regionMap[vertexHandle.idx()] = distMap[vertexHandle].second;

  printf("%lu sources computed\n", sourceVertices.size());

  return std::make_pair(sourceVertices, regionMap);
} // sample

// given sample points, compute edges based on crossfield adjacency 
std::vector<std::set<size_t>> MeshSamplingCentroidVoronoiEdges(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<Eigen::Vector3d> &crossField)
{
  // start with points, use crossField to connect two points
    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Eigen::Vector3d> dirs;

    for(auto vertex : sourceVertices) {
      auto pnt = mesh.point(vertex);
    
      auto normal = mesh.normal(vertex);

      points.push_back(toVec(pnt));
      normals.push_back(toVec(normal));
//      dirs.push_back(crossField[vertex.idx()]);
    }

    happly::PLYExport exPly;

    std::vector<std::set<size_t>> adjMatrix;
    bool makeSymmetric = true;
  
    adjMatrix.resize(points.size());

    for(size_t iPntDir = 0; iPntDir < 4*points.size(); iPntDir++) {

      size_t iPnt = iPntDir/4;
      size_t iDir = iPntDir % 4;

      // walk the mesh from sourceVertices[iPnt] in each of 4 crossFields[iPnt] directions,
      MyMesh::VertexHandle vertexHandle = sourceVertices[iPnt];
      MyMesh::HalfedgeHandle hedgeHandle = *(mesh.voh_iter(vertexHandle));

      auto point = toVec(mesh.point(vertexHandle));
      auto normal = toVec(mesh.normal(vertexHandle));
      auto dirVec = crossField[vertexHandle.idx()];

      for(size_t iiDir = 0; iiDir < iDir; iiDir++)
        dirVec = dirVec.cross(normal);
      
//      printf("%lu. (%.6f, %.6f, %.6f)\n", iPnt, point[0], point[1], point[2]);

      bool reached = false;
      int iIter = 0;

      size_t iAdjPnt = iPnt;

      while(!reached) {

        // TODO : cluster check implementation

        exPly.AddVertex({point[0], point[1], point[2]}, 0.8, 0.2, 0.8);

//        printf("   %d. (%.6f, %.6f, %.6f) -> (%.6f, %.6f, %.6f)\n", iIter++, point[0], point[1], point[2], dirVec[0], dirVec[1], dirVec[2]);

        std::pair<Eigen::Vector3d, MyMesh::HalfedgeHandle> step = MakeStep(mesh, point, dirVec, hedgeHandle, vertexHandle);

        if(step.second.is_valid()) {
          // compute new direction

          MyMesh::VertexHandle toVertex = mesh.to_vertex_handle(step.second);
          MyMesh::VertexHandle fromVertex = mesh.from_vertex_handle(step.second);

          auto toPnt = mesh.point(toVertex);
          auto fromPnt = mesh.point(fromVertex);
          size_t iV = toVertex.idx();
          size_t jV = fromVertex.idx();

          point = step.first;
          hedgeHandle = step.second;
          
          Line3d edgeSeg(toVec(fromPnt), toVec(toPnt));
          double param = edgeSeg.ProjectParam(point);

          auto toNorm = mesh.normal(toVertex);
          auto fromNorm = mesh.normal(fromVertex);

          Eigen::Vector3d toNormal = toVec(toNorm); 
          Eigen::Vector3d fromNormal = toVec(fromNorm); 

          Eigen::Vector3d toDirVec = bestRosyFit(dirVec, toNormal, crossField[toVertex.idx()]);
          Eigen::Vector3d fromDirVec = bestRosyFit(dirVec, fromNormal, crossField[fromVertex.idx()]);

          Eigen::Vector3d nextDirVec = toDirVec*param + fromDirVec*(1.0-param);
          nextDirVec.stableNormalize();

          dirVec = nextDirVec;

          vertexHandle = MyMesh::VertexHandle();
          if(param < eps)
            vertexHandle = fromVertex;
          if(param > 1.0-eps)
            vertexHandle = toVertex;

          iAdjPnt = regionMap[mesh.to_vertex_handle(hedgeHandle).idx()];
          if(iAdjPnt != iPnt )
            reached = true;
            // TODO : add graph vertex here
        }
        else
          reached = true;


        // take current point, 
        // compute new cross field by edge interpolation

        // if edge then take the other facet

        // if vertes then circulate the vertex facets
        // project new dir onto facet plane
        // take e1.cross(dir), dir.cross(e2) where e1 and e2 are two edges coming out of the vertex
        // if cross products are coplanar then take that facet


        // move across the facet in dir till edge intersection

        // if selected facet belongs to another region then 
        if(iIter == 1000)
          reached = true;
      }

      if(iAdjPnt != iPnt) {

        double dist = (points[iPnt]-points[iAdjPnt]).norm();
        if(dist >= 6.0)
          printf("!!! long distance %.6f !!!\n", dist);

        else {
        adjMatrix[iPnt].insert(iAdjPnt);
        if(makeSymmetric)
          adjMatrix[iAdjPnt].insert(iPnt);
        }
      }
    }

    exPly.WritePLY("crosspath.ply");

    
    // TODO : smooth while constraining to surface

    return adjMatrix;
} // Edges

// given sample points, compute edges based on voronoi cell adjacency
std::vector<std::set<size_t>> MeshSamplingPidEdges(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<Eigen::Vector3d> &crossField)
{
    std::vector<std::set<size_t>> adjMatrix;
    bool makeSymmetric = true;
  
    adjMatrix.resize(sourceVertices.size());

    for(auto hedgeHandle : mesh.halfedges()) {

      // get two vertices
      auto toVertex = mesh.to_vertex_handle(hedgeHandle);
      auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

      int iRegion = regionMap[toVertex.idx()];
      int jRegion = regionMap[fromVertex.idx()];

      if(iRegion != jRegion) {
        adjMatrix[iRegion].insert(jRegion);
        adjMatrix[jRegion].insert(iRegion);
      }
    }

    return adjMatrix;
}

#if 0
std::vector<std::set<size_t>> MeshSamplingGeometricEdges(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, std::vector<Eigen::Vector3d> &points, std::vector<Eigen::Vector3d> &crossField)
{
  // look for geometric edges that fit the cross field node structure

    printf("MeshSamplingGeometricEdges\n");

    std::vector<Eigen::Vector3d> normals;
    std::vector<Eigen::Vector3d> dirs;

    for(auto vertex : sourceVertices) {    
      auto normal = mesh.normal(vertex);

      normals.push_back(toVec(normal));
      dirs.push_back(crossField[vertex.idx()]);
    }

    std::vector<std::set<size_t>> adjMatrix;
    bool makeSymmetric = true;
  
    adjMatrix.resize(points.size());

    double step = 2.0; // TODO : param

    for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {
      for(size_t jPnt = 0; jPnt < points.size(); jPnt++) {
        double pdist = sqrt(distSq(points[iPnt], points[jPnt]) );
        if(iPnt != jPnt && pdist > 0.75*step && pdist < 1.5*step) {

          Eigen::Vector3d edgeVec = points[jPnt]-points[iPnt];

          auto normal = normals[iPnt];
          auto fieldDir = dirs[iPnt];

          edgeVec -= normal*(edgeVec.dot(normal));
          edgeVec.stableNormalize();

          Eigen::Vector3d bestFit = bestRosyFit(edgeVec, normal, fieldDir);

          // test vector against cross field
          if(bestFit.dot(edgeVec) > 0.9) {
            adjMatrix[iPnt].insert(jPnt);
            adjMatrix[jPnt].insert(iPnt);
          }
        }
      }
    }

    happly::PLYExport eply;

    for(size_t iPnt = 0; iPnt < sourceVertices.size(); iPnt++) {
      // get all edges
      for(auto jPnt : adjMatrix[iPnt]) {

        // TODO : project onto the plane
        Eigen::Vector3d dirVec = points[iPnt]-points[jPnt];
        dirVec.stableNormalize();

        // iterate over directions, find the best dot prod

         // check if good match with dir
        // compute edge vec

        double bestDotProd = -1.0;

        auto normal = normals[iPnt];
        auto dir = dirs[iPnt];

        for(int ii = 0; ii < 4; ii++) {

          double dotProd = dir.dot(dirVec);

          bestDotProd = std::max(dotProd, bestDotProd);

          dir = cross(dir, normal);
        }

        // TODO : score edge based on dotprod

        bestDotProd = (bestDotProd -0.5)*2.0;
        bestDotProd = std::max(0.0, bestDotProd);
        bestDotProd = std::min(1.0, bestDotProd);
        double blue = bestDotProd;
        double red = 1.0-bestDotProd;
        double green = 4.0*blue*red;

        int nPnts = 11;


        for(int jj = 1; jj < nPnts; jj+=2) {

              double param = 1.0*jj/nPnts;
              auto midPnt = points[iPnt]*param + points[jPnt]*(1.0-param);

              eply.AddVertex({midPnt[0], midPnt[1], midPnt[2]}, red, green, blue);
        }
      }
    }

    eply.WritePLY("crossgeoedges.ply");

    return adjMatrix;
} // AllEdges
#endif

// quick function to compute edges satisfying sufficient condition for delaunay edge
std::vector<std::set<size_t>> MeshSamplingEdges(const std::vector<Eigen::Vector3d> &points)
{

  std::vector<std::set<size_t>> adjMatrix;
  adjMatrix.resize(points.size());

  for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

    // TODO : find all neighbor candidates efficiently
    std::vector<size_t> neighbors;
    for(size_t iNPnt = 0; iNPnt < points.size(); iNPnt++)
      if(iNPnt != iPnt)
        neighbors.push_back(iNPnt);

    // TODO : redo the algorithm
    // find n neighbor candidates
    // nxn matrix
    // for each neighbor, find those overshadowed by it
    // find neighbors are not overshadowed by anybody

    std::vector<bool> isSeen(points.size(), true);
    isSeen[iPnt] = false;

    for(size_t iNeighbor : neighbors) {
      Eigen::Vector3d normal = (points[iNeighbor] - points[iPnt]);
      normal.stableNormalize();
      Plane3d plane(points[iNeighbor], normal);
      for(size_t jNeighbor : neighbors) {
        if(iNeighbor != jNeighbor && plane.SignedDistance(points[jNeighbor]) >= 0.0)
          isSeen[jNeighbor] = false;
      }
    }

    int nSeen = 0;
    for(auto seen : isSeen)
      if(seen)
        nSeen++;

    std::set<size_t> addedNeighbors;

    while (neighbors.size()) {

//      printf("    >> ");
//      for(auto in : neighbors)
//        printf("%lu ", in);

//      printf("\n");

      // find closest neighbor such that
      // there is no add neighbor this neighbor is closer to than curr point
      // this is the quick sufficient condition for the edge (iPnt, neighbor)
      // to be in Delaunay triangulation

      double minDist = -1.0;
      size_t iToAddNeighbor = 0;
      std::vector<size_t> nextNeighbors;

      for(size_t ii = 0; ii < neighbors.size(); ii++) {

        size_t iNeighbor = neighbors[ii];

        double currDist = (points[iPnt] - points[iNeighbor]).norm();

        bool isDelaunay = true;
        bool alreadyAdded = false;

        for(auto addedNeighbor : addedNeighbors) {

          if(addedNeighbor == iNeighbor) {
            alreadyAdded = true;
            break;
          }
          else {
            Eigen::Vector3d planeNormal = points[addedNeighbor]-points[iPnt];
            planeNormal.stableNormalize();

            Plane3d plane(points[addedNeighbor], planeNormal);

            double neighborDist = plane.SignedDistance(points[iNeighbor]);
            if( neighborDist > 0.0 /*|| planeNormal.norm() < eps*/) {
              isDelaunay = false;
              break;
            }
          }
        }

        if(isDelaunay) {
          // put in next neighbors
          if(!alreadyAdded) {
            nextNeighbors.push_back(iNeighbor);

            // check if this should iToAddNeighbor
            if(minDist < 0.0 || minDist > currDist) {
              minDist = currDist;
              iToAddNeighbor = iNeighbor;
            }
          }
        }
        else if (!alreadyAdded && isSeen[iNeighbor])
          printf("!!! discarding visible neighbor %d !!!\n", iNeighbor);

      }

      if(minDist > -1.0) {
        printf("  add %d (%.6f, %.6f, %.6f) at dist %.6f\n", iToAddNeighbor, points[iToAddNeighbor][0], points[iToAddNeighbor][1], points[iToAddNeighbor][2], minDist);
        addedNeighbors.insert(iToAddNeighbor);
      }

      printf("  %lu before, %lu after\n", neighbors.size(), nextNeighbors.size());
      neighbors = std::move(nextNeighbors);
    }

    printf("%lu. (%.6f, %.6f, %.6f), seen %d, added %lu neighbors\n\n", iPnt, points[iPnt][0], points[iPnt][1], points[iPnt][2], nSeen, addedNeighbors.size());

    for(auto addedNeighbor : addedNeighbors)
//      if(isSeen[addedNeighbor])
        adjMatrix[iPnt].insert(addedNeighbor);

    happly::PLYExport eply;

    eply.AddVertex(points[iPnt], 1.0, 0.0, 0.0);

    for(auto in : adjMatrix[iPnt]) {

      eply.AddVertex(points[in], 0.0, 0.0, 1.0);

        int nPnts = 13;

        for(int jj = 1; jj < nPnts; jj+=2) {

              double param = 1.0*jj/nPnts;
              auto midPnt = points[iPnt]*param + points[in]*(1.0-param);

              eply.AddVertex({midPnt[0], midPnt[1], midPnt[2]}, 0.0, 1.0, 0.0);
        }

    }

    eply.WritePLY("delaunay.ply");

//    getchar();
  }

  // TODO : count edges that do not appear in both directions

  bool addSymmetric = false;

  if(addSymmetric) {
    for(size_t iPnt = 0; iPnt < adjMatrix.size(); iPnt++) {
      for(size_t jPnt : adjMatrix[iPnt])
        if(iPnt != jPnt) {
          adjMatrix[jPnt].insert(iPnt);
        }
    }
  }
  else {

    std::vector<std::pair<size_t, size_t>> edgesToRemove;

    for(size_t iPnt = 0; iPnt < adjMatrix.size(); iPnt++)
      for(int jPnt : adjMatrix[iPnt]) {
        if(iPnt == jPnt || adjMatrix[jPnt].find(iPnt) == adjMatrix[jPnt].end())
          edgesToRemove.push_back({iPnt, jPnt});
      }

    size_t nEdges = 0;

    for(const auto &adj : adjMatrix)
      nEdges += adj.size();

    printf("%lu vertices, %lu half edges, %lu to remove\n", adjMatrix.size(), nEdges, edgesToRemove.size());

    for(const auto &edgeToRemove : edgesToRemove) {
      adjMatrix[edgeToRemove.first].erase(edgeToRemove.second);
      adjMatrix[edgeToRemove.second].erase(edgeToRemove.first);
    }

  }

  return adjMatrix;
}

#if 0
void MeshSamplingQuadQualityEdges(MyMesh &mesh, std::vector<Eigen::Vector3d> &points, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<std::set<size_t>> &adjMatrix, std::vector<Eigen::Vector3d> &crossField)
{
  // show "quad quality" of each edge
    printf("MeshSamplingQuadEdges\n");

    std::vector<Eigen::Vector3d> normals;

    for(auto vertex : sourceVertices) {
      auto pnt = mesh.point(vertex);
      auto normal = mesh.normal(vertex);

      normals.push_back(toVec(normal));
    }

    happly::PLYExport eply;

    for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

//      printf("quad edges of %lu, %lu neighbors\n", iPnt, adjMatrix[iPnt].size());

      if(adjMatrix[iPnt].size() == 0)
        continue;

      auto point = points[iPnt];
      auto normal = normals[iPnt];

      std::vector<Eigen::Vector3d> projEdgeVecs;
      std::vector<Eigen::Vector3d> edgePoints;

      for(auto jPnt : adjMatrix[iPnt] ) {
          auto toPoint = points[jPnt];

          edgePoints.push_back(toPoint);
          Eigen::Vector3d edgeVec = toPoint - point;

          // project onto the plane
          edgeVec -= normal*(normal.dot(edgeVec));
          edgeVec.stableNormalize();

          projEdgeVecs.push_back(edgeVec);
      }

      Eigen::MatrixXd dotMatrix = Eigen::MatrixXd::Zero(projEdgeVecs.size(), projEdgeVecs.size());

      for(size_t iVec = 0; iVec < projEdgeVecs.size(); iVec++) {

        auto dirVec = projEdgeVecs[iVec];

        for(int ii = 0; ii < 4; ii++) {

        // compute angles between dir and edge vectors
          for(size_t jVec = 0; jVec < projEdgeVecs.size(); jVec++) {

            double dotProd = dirVec.dot(projEdgeVecs[jVec]);
            double currDotProd = dotMatrix(iVec, jVec);

            dotMatrix(iVec, jVec) = std::max(dotProd, currDotProd);            
          }

          dirVec = cross(dirVec, normal);
        }
      }

      for(size_t iVec = 0; iVec < projEdgeVecs.size(); iVec++) {

        double totalDot = 0.0;
        for(size_t jVec = 0; jVec < projEdgeVecs.size(); jVec++)
          totalDot += dotMatrix(iVec, jVec);

        totalDot /= projEdgeVecs.size();

        totalDot = (totalDot-0.5)*2.0;

        totalDot = std::max(totalDot, 0.0);
        totalDot = std::min(totalDot, 1.0);

        // color edges by totalDot

        int nPnts = 11;


        double red = 1.0-totalDot;
        double blue = totalDot;
        double green = 4.0*red*blue;


        for(int jj = 1; jj < nPnts; jj+=2) {

              double param = 1.0*jj/nPnts;
              auto midPnt = point*param + edgePoints[iVec]*(1.0-param);

              eply.AddVertex({midPnt[0], midPnt[1], midPnt[2]}, red, green, blue);
        }
      }

    }

    eply.WritePLY("crossdiredges.ply");


}
#endif

#if 0
void MeshSamplingCrossEdges(MyMesh &mesh, std::vector<Eigen::Vector3d> &points, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<std::set<size_t>> &adjMatrix, std::vector<Eigen::Vector3d> &crossField)
{
  // quad quality of edges based on crossfield alignment

    printf("MeshSamplingCrossEdges %lu\n", points.size());

    std::vector<Eigen::Vector3d> normals;
    std::vector<Eigen::Vector3d> dirs;

    for(auto vertex : sourceVertices) {
      auto pnt = mesh.point(vertex);
      auto normal = mesh.normal(vertex);

      normals.push_back(toVec(normal));
      dirs.push_back(crossField[vertex.idx()]);
    }

    int nCoinc = 0;

    struct sDirEdge {
      Eigen::Vector3d dir;

      int edgeI;
      double dotProd;
    };

    std::vector<std::vector<sDirEdge>> crossEdges;
    crossEdges.resize(adjMatrix.size());

    int nFourPoints = 0;
 
    struct sCalcDirEdge {
      Eigen::Vector3d dir;
      Eigen::Vector3d edgeDir;

      int edgeI;
      double dotProd;
    };

    for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

//      ControlPrintf printf;

      printf("smooth %lu, %lu neighbors\n", iPnt, adjMatrix[iPnt].size());

      if(adjMatrix[iPnt].size() == 0)
        continue;

      auto point = points[iPnt];
      auto normal = normals[iPnt];
      auto dir = dirs[iPnt];

      std::vector<Eigen::Vector3d> crossDirs;

      for(int ii = 0; ii < 4; ii++) {

        if(crossDirs.size() == 0)
          crossDirs.push_back(dir);
        else 
          crossDirs.push_back(crossDirs.back().cross(normal));
      }

      // 2. find edges corresponding to cross field directions
     
      for(int ii = 0; ii < 4; ii++) {

        std::vector<sCalcDirEdge> thisDirEdges;

        // compute angles between dir and edge vectors
        for(auto jPnt : adjMatrix[iPnt] ) {

          auto toPoint = points[jPnt];

          Eigen::Vector3d edgeVec = toPoint - point;

          // project onto the plane
          edgeVec -= normal*(normal.dot(edgeVec));
          edgeVec.stableNormalize();
    
          double dotProd = crossDirs[ii].dot(edgeVec);
          
          thisDirEdges.push_back({dir, edgeVec, jPnt, dotProd});
        }

        size_t iBest = 0;
        double bestDot = -1.0;

        // choose edge candidates
        for(size_t iEdge = 0; iEdge < thisDirEdges.size(); iEdge++) {
          if(bestDot < thisDirEdges[iEdge].dotProd) {
            iBest = iEdge;
            bestDot = thisDirEdges[iEdge].dotProd;
          }
        }

        // good fit, take the best
        if(bestDot > 0.98) {
          crossEdges[iPnt].push_back({thisDirEdges[iBest].dir, thisDirEdges[iBest].edgeI, thisDirEdges[iBest].dotProd});
        }
        // fair fit, take all the candidates
        else if (bestDot > 0.90) {

          for(size_t iEdge = 0; iEdge < thisDirEdges.size(); iEdge++)
            if( thisDirEdges[iEdge].dotProd > 0.90) {
              crossEdges[iPnt].push_back({thisDirEdges[iEdge].dir, thisDirEdges[iEdge].edgeI, thisDirEdges[iEdge].dotProd});
            }
        }
        else {
          crossEdges[iPnt].push_back({thisDirEdges[iBest].dir, thisDirEdges[iBest].edgeI, thisDirEdges[iBest].dotProd});
        } 
      }

      std::set<int> inCrossEdges;

      for(const auto &crossEdge : crossEdges[iPnt]) {
          inCrossEdges.insert(crossEdge.edgeI);
      }

      if(inCrossEdges.size() == crossEdges[iPnt].size())
        nFourPoints++;
    }

    // TODO : analyze all edges
    // if edge (i, j) found but crossEdges[i] to j && crossEdges[j] to i do not exist
    // then remove the edge

    printf("%d four points\n", nFourPoints);

    std::vector<std::pair<int, int>> edgesToRemove;

    for(int iPnt = 0; iPnt < points.size(); iPnt++) {

      // get vertex
      // get normal
      // get dir

      // iterate over edges
      auto point = points[iPnt];
      auto normal = normals[iPnt];
      auto dir = dirs[iPnt];

      for(int jPnt : adjMatrix[iPnt]) {

        bool thisFound = false;
        for(const auto &crossEdge : crossEdges[iPnt]) {
          if( crossEdge.edgeI == jPnt)
            thisFound = true;
        }

        bool otherFound = false;
        for(const auto &crossEdge : crossEdges[jPnt]) {
          if(crossEdge.edgeI == iPnt)
            otherFound = true;
        }

        if(!thisFound && !otherFound)
          edgesToRemove.push_back({iPnt, jPnt});
      }
    }


    size_t nCrossEdges = 0;
    size_t nEdges = 0;

    for(const auto &adj : adjMatrix)
      nEdges += adj.size();

    for(const auto &cross : crossEdges)
      nCrossEdges += cross.size();

    printf("%lu vertices, %lu edges, %lu cross, %lu to remove\n", adjMatrix.size(), nEdges/2, nCrossEdges/2, edgesToRemove.size());
//    getchar();

    for(const auto &edgeToRemove : edgesToRemove) {
      adjMatrix[edgeToRemove.first].erase(edgeToRemove.second);
      adjMatrix[edgeToRemove.second].erase(edgeToRemove.first);
    }

  // save in a ply file

    happly::PLYExport eply;
    for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

      // get vertex
      // get normal
      // get dir

      // iterate over edges
      auto point = points[iPnt];
      auto normal = normals[iPnt];
      auto dir = dirs[iPnt];

      eply.AddVertex(point, 0.0, 1.0, 0.0);
//      eply.AddVertex(points[iPnt], 0.1, 0.8, 0.1);

      for(auto jPnt : adjMatrix[iPnt] ){

        auto toPoint = points[jPnt];

        Eigen::Vector3d edgeVec = toPoint - point;

        // project onto the plane
        edgeVec -= normal*(normal.dot(edgeVec));
        edgeVec.stableNormalize();

        // TODO : fiond out what gets modified
        auto bestDirVec = bestRosyFit(edgeVec, normal, dir);
      
        double dotProd = bestDirVec.dot(edgeVec);
      
        double w = (dotProd-0.5)*2.0;

        Eigen::Vector3d midPoint = point*0.67 + toPoint*0.33;

        double red = 1.0-w;
        double blue = w;
        double green = 4.0*red*blue;

//        eply.AddVertex(midPoint, red, green, blue);
      }
    }

    eply.WritePLY("voroedgeffit.ply");
 
    // TODO : smooth while constraining to surface
} // smooth
#endif

#if 0
std::vector<Eigen::Vector3d> MeshSamplingSmooth(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<std::set<size_t>> &adjMatrix, std::vector<Eigen::Vector3d> &crossField)
{
    // smooth based on graph

    printf("MeshSamplingSmooth\n");

    // get points
    AABBTree tree;
    tree.Build(mesh, 20);

    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Eigen::Vector3d> dirs;

    for(auto vertex : sourceVertices) {
      auto pnt = mesh.point(vertex);
      auto normal = mesh.normal(vertex);

      points.push_back(toVec(pnt));
      normals.push_back(toVec(normal));
      dirs.push_back(crossField[vertex.idx()]);
    }

    for(int iIter = 0; iIter < 1; iIter++) {

      if(iIter > 0)
        for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {    
          auto vertex = sourceVertices[iPnt];
          auto normal = mesh.normal(vertex);

          normals[iPnt] = toVec(normal);
          dirs[iPnt] = crossField[vertex.idx()];
        }

    int nCoinc = 0;

    struct sDirEdge {
      Eigen::Vector3d dir;

      int edgeI;
      double dotProd;
    };

    std::vector<std::vector<sDirEdge>> crossEdges;
    crossEdges.resize(adjMatrix.size());

    int nFourPoints = 0;
 
    struct sCalcDirEdge {
      Eigen::Vector3d dir;
      Eigen::Vector3d edgeDir;

      int edgeI;
      double dotProd;
    };

    ControlPrintf printf;

    std::vector<Eigen::Vector3d> nextPoints;
    nextPoints.resize(points.size());

    for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

      printf("smooth %lu(%lu).\n", iPnt, points.size());

      auto point = points[iPnt];
      auto normal = normals[iPnt];
      auto dir = dirs[iPnt];

      std::vector<Eigen::Vector3d> crossDirs;

      for(int ii = 0; ii < 4; ii++) {

        if(crossDirs.size() == 0)
          crossDirs.push_back(dir);
        else 
          crossDirs.push_back(crossDirs.back().cross(normal));
      }

      nextPoints[iPnt] = Eigen::Vector3d(0.0,0.0,0.0);
      double totalWeight = 0.0;

      // 2. find edges corresponding to cross field directions
     
      for(int ii = 0; ii < 4; ii++) {

        std::vector<sCalcDirEdge> thisDirEdges;

        // compute angles between dir and edge vectors
        for(auto jPnt : adjMatrix[iPnt] ) {

          auto toPoint = points[jPnt];

          Eigen::Vector3d edgeVec = toPoint - point;

          // project onto the plane
          edgeVec -= normal*(normal.dot(edgeVec));
          edgeVec.stableNormalize();
    
          double dotProd = crossDirs[ii].dot(edgeVec);
          
          thisDirEdges.push_back({dir, edgeVec, jPnt, dotProd});
        }

        size_t iBest = 0;
        double bestDot = -1.0;

        // choose edge candidates
        for(size_t iEdge = 0; iEdge < thisDirEdges.size(); iEdge++) {
          if(bestDot < thisDirEdges[iEdge].dotProd) {
            iBest = iEdge;
            bestDot = thisDirEdges[iEdge].dotProd;
          }
        }

        // 3. smooth vertices based on cross edges

        // good fit, take the best
        if(bestDot > 0.98) {
          crossEdges[iPnt].push_back({thisDirEdges[iBest].dir, thisDirEdges[iBest].edgeI, thisDirEdges[iBest].dotProd});

          nextPoints[iPnt] += points[thisDirEdges[iBest].edgeI];
          totalWeight++;
        }
        // fair fit, take all the candidates
        else if (bestDot > 0.90) {

          for(size_t iEdge = 0; iEdge < thisDirEdges.size(); iEdge++)
            if( thisDirEdges[iEdge].dotProd > 0.90) {
              crossEdges[iPnt].push_back({thisDirEdges[iEdge].dir, thisDirEdges[iEdge].edgeI, thisDirEdges[iEdge].dotProd});

              nextPoints[iPnt] += points[thisDirEdges[iEdge].edgeI];
              totalWeight++;
            }
        }
        else {
          crossEdges[iPnt].push_back({thisDirEdges[iBest].dir, thisDirEdges[iBest].edgeI, thisDirEdges[iBest].dotProd});
          nextPoints[iPnt] += points[thisDirEdges[iBest].edgeI];
          totalWeight++;
        } 
      }

      if(totalWeight > 0.0)
        nextPoints[iPnt] *= 1.0/totalWeight;
      else
        nextPoints[iPnt] = points[iPnt];

      printf("   %.6f, (%.6f, %.6f, %.6f)\n", totalWeight, nextPoints[iPnt][0], nextPoints[iPnt][1], nextPoints[iPnt][2]);

      std::set<int> inCrossEdges;

      for(const auto &crossEdge : crossEdges[iPnt]) {
          inCrossEdges.insert(crossEdge.edgeI);
      }

      if(inCrossEdges.size())
        nFourPoints++;

    }

    //4. project onto mesh

    double totalShift = 0.0;

    for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

      auto np = tree.FindNearestPoint(nextPoints[iPnt]);

      totalShift += (points[iPnt]-np.pnt).norm();

      points[iPnt] = np.pnt;

      MyMesh::FaceHandle facetHandle(np.facetIndex);

      double closestDist = -1.0;

      MyMesh::VertexHandle closestVertex;

      for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {
        auto thisVertex = *fv_it;

        auto thisPnt = toVec(mesh.point(thisVertex));

        double thisDist = (points[iPnt]-thisPnt).norm();

        if(closestDist < 0.0 || closestDist > thisDist) {
          closestDist = thisDist;
          closestVertex = thisVertex;
        }
      }

      if(closestVertex.is_valid())
        sourceVertices[iPnt] = closestVertex;
    }

    printf("%lu vertices, %d four point, %.6f\n", points.size(), nFourPoints, totalShift);

    }

    // TODO : analyze all edges
    // if edge (i, j) found but crossEdges[i] to j && crossEdges[j] to i do not exist
    // then remove the edge

    std::vector<std::pair<int, int>> edgesToRemove;

#if 0    

    for(int iPnt = 0; iPnt < points.size(); iPnt++) {

      // get vertex
      // get normal
      // get dir

      // iterate over edges
      auto point = points[iPnt];
      auto normal = normals[iPnt];
      auto dir = dirs[iPnt];

      for(int jPnt : adjMatrix[iPnt]) {

        bool thisFound = false;
        for(const auto &crossEdge : crossEdges[iPnt]) {
          if( crossEdge.edgeI == jPnt)
            thisFound = true;
        }

        bool otherFound = false;
        for(const auto &crossEdge : crossEdges[jPnt]) {
          if(crossEdge.edgeI == iPnt)
            otherFound = true;
        }

        if(!thisFound && !otherFound)
          edgesToRemove.push_back({iPnt, jPnt});
      }
    }


    size_t nCrossEdges = 0;
    size_t nEdges = 0;

    for(const auto &adj : adjMatrix)
      nEdges += adj.size();

    for(const auto &cross : crossEdges)
      nCrossEdges += cross.size();

    printf("%lu vertices, %lu edges, %lu cross, %lu to remove\n", adjMatrix.size(), nEdges/2, nCrossEdges/2, edgesToRemove.size());
    getchar();
#endif

//    if(false)
//    for(const auto &edgeToRemove : edgesToRemove) {
//      adjMatrix[edgeToRemove.first].erase(edgeToRemove.second);
//      adjMatrix[edgeToRemove.second].erase(edgeToRemove.first);
//    }

  // save in a ply file

  happly::PLYExport eply;
  for(size_t iPnt = 0; iPnt < points.size(); iPnt++) {

      // get vertex
      // get normal
      // get dir

      // iterate over edges
      auto point = points[iPnt];
      auto normal = normals[iPnt];
      auto dir = dirs[iPnt];

      eply.AddVertex(point, 0.0, 1.0, 0.0);
//      eply.AddVertex(points[iPnt], 0.1, 0.8, 0.1);

      printf("%lu adns\n", adjMatrix[iPnt].size());
      for(auto jPnt : adjMatrix[iPnt] ){

        auto toPoint = points[jPnt];

        Eigen::Vector3d edgeVec = toPoint - point;

        // project onto the plane
        edgeVec -= normal*(normal.dot(edgeVec));
        edgeVec.stableNormalize();

        // TODO : fiond out what gets modified
        auto bestDirVec = bestRosyFit(edgeVec, normal, dir);
      
        double dotProd = bestDirVec.dot(edgeVec);
      
        double w = (dotProd-0.5)*2.0;

        Eigen::Vector3d midPoint = point*0.67 + toPoint*0.33;

        double red = 1.0-w;
        double blue = w;
        double green = 4.0*red*blue;

        eply.AddVertex(midPoint, red, green, blue);
      }
    }

    eply.WritePLY("voroedgeffit.ply");
 
    // TODO : smooth while constraining to surface

} // smooth
#endif

// quick function to look for triangles and quads to create mesh object for sourceVertices connected by edges
void MeshSamplingTriQuads(MyMesh &mesh, const std::vector<MyMesh::VertexHandle> &sourceVertices, std::vector<std::set<size_t>> &edges) 
{
    // look for 3 and 4 loops

    // 3. build loops

    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector3d> normals;
    std::vector<Eigen::Vector3d> dirs;

    for(auto vertex : sourceVertices) {
      auto pnt = mesh.point(vertex);
    
      auto normal = mesh.normal(vertex);

      points.push_back(toVec(pnt));
      normals.push_back(toVec(normal));
//      dirs.push_back(crossField[vertex.idx()]);
    std::vector<std::vector<size_t>> facets;

    std::map<std::pair<size_t, size_t>, size_t> edgeMap;

    int nMissed = 0;
    for(size_t iV = 0; iV < edges.size(); iV++) {
      // check if neighbor adjacency exists
      bool found = false;

      std::set<size_t> doubleEdges;

      for(size_t jV : edges[iV]) {
        if( edges[jV].find(iV) != edges[jV].end() ) {
          doubleEdges.insert(jV);
        }
      }

      nMissed += edges[iV].size()-doubleEdges.size();
      if(edges[iV].size() > doubleEdges.size())
        edges[iV] = std::move(doubleEdges);
    }

    int nTriangles = 0;
    int nQuads = 0;
    int averDegree = 0;

    for(size_t iV = 0; iV < edges.size(); iV++) {
          averDegree += edges[iV].size();

          for(size_t jV : edges[iV])
            for(size_t kV : edges[iV])
              if(jV != kV) {

                // check for triangle
                bool triangleFound = false;
                for(size_t kkV : edges[jV])
                  if(kkV == kV) {
                    nTriangles++;
                    triangleFound = true;

                    // TODO : filter out duplicate
                    facets.push_back({iV, jV, kV});
                  }

                // check for quad
                bool quadFound = false;
                if(!triangleFound )
                  for(size_t lV : edges[jV])
                    if(lV != iV)
                      for(size_t kkV : edges[lV])
                        if(kkV == kV && kkV != jV ) {
                          nQuads++;
                          quadFound = true;

                          // TODO : filter out duplicate
                          facets.push_back({iV, jV, lV, kV});
                        }
              }
    }

    printf("%lu vertices, %d missed edges, %.6f aver. degree, %d triangles, %d quads\n", edges.size(), nMissed, 1.0*averDegree/edges.size(), nTriangles/3, nQuads/4);

        int nPnts = 11;
        happly::PLYExport eply;
        for(size_t iV = 0; iV < edges.size(); iV++) {
          for(size_t jV : edges[iV]) {

            for(int iPnt = 1; iPnt < nPnts; iPnt+=2) {

              double param = 1.0*iPnt/nPnts;
              auto midPnt = points[iV]*param + points[jV]*(1.0-param);

              eply.AddVertex({midPnt[0], midPnt[1], midPnt[2]}, 0.0, 0.0, 1.0);
            }
          }
        }
        eply.WritePLY("edgessstriquads.ply");
   }
}

// function to detect simple non-separating cycles in graph specified by points and edges using DFS, create facet per each cycle
MyMesh MeshSamplingCycles(const std::vector<Eigen::Vector3d> &points, std::vector<std::set<size_t>> &edges)
{
  // detect cycles, create facets

    printf("MeshSamplingCycles\n");

    // 3. build loops

    int nEdges = 0;

    for(const auto &edge : edges)
      nEdges += edge.size();

    nEdges /= 2;

    printf("MeshSamplingCycles : %d edges\n", nEdges);

//    std::vector<std::vector<size_t>> facets;
//    std::map<std::pair<size_t, size_t>, size_t> edgeMap;

    int nMissed = 0;
    for(size_t iV = 0; iV < edges.size(); iV++) {
      // check if neighbor adjacency exists
      bool found = false;

      std::set<size_t> doubleEdges;

      for(size_t jV : edges[iV]) {
        if( edges[jV].find(iV) != edges[jV].end() ) {
          doubleEdges.insert(jV);
        }
      }

      nMissed += edges[iV].size()-doubleEdges.size();
      if(edges[iV].size() > doubleEdges.size()) {
        edges[iV] = std::move(doubleEdges);
        printf("remove not symmetric\n");
      }
    }

  #if 0  
    int nTriangles = 0;
    int nQuads = 0;
    int averDegree = 0;

    for(size_t iV = 0; iV < edges.size(); iV++) {
          averDegree += edges[iV].size();

          for(size_t jV : edges[iV])
            for(size_t kV : edges[iV])
              if(jV != kV) {

                // check for triangle
                bool triangleFound = false;
                for(size_t kkV : edges[jV])
                  if(kkV == kV) {
                    nTriangles++;
                    triangleFound = true;

                    // TODO : filter out duplicate
                    facets.push_back({iV, jV, kV});
                  }

                // check for quad
                bool quadFound = false;
                if(!triangleFound )
                  for(size_t lV : edges[jV])
                    if(lV != iV)
                      for(size_t kkV : edges[lV])
                        if(kkV == kV && kkV != jV ) {
                          nQuads++;
                          quadFound = true;

                          // TODO : filter out duplicate
                          facets.push_back({iV, jV, lV, kV});
                        }
              }
    }

    printf("%lu vertices, %d missed edges, %.6f aver. degree, %d triangles, %d quads\n", edges.size(), nMissed, 1.0*averDegree/edges.size(), nTriangles/3, nQuads/4);
    getchar();

#endif    

    // TODO : color by angle 
        int nPnts = 11;
        happly::PLYExport eply;
        for(size_t iV = 0; iV < edges.size(); iV++) {
          for(size_t jV : edges[iV]) {

            for(int iPnt = 1; iPnt < nPnts; iPnt+=2) {

              double param = 1.0*iPnt/nPnts;
              auto midPnt = points[iV]*param + points[jV]*(1.0-param);

              eply.AddVertex({midPnt[0], midPnt[1], midPnt[2]}, 0.0, 0.0, 1.0);
            }
          }
        }
        eply.WritePLY("edgesss.ply");

    happly::PLYExport eeply;

    for(auto pnt : points)
        eeply.AddVertex(pnt, 0.1, 0.8, 0.4);
        
    std::vector<std::vector<int>> loops;
    std::map<std::pair<int, int>, int> isEdgeInLoop;

    for(int ii = 0; ii < edges.size(); ii++)
      for(auto jj : edges[ii]) 
        isEdgeInLoop[{ii, jj}] = -1;

    int nVertices = static_cast<int>(edges.size());

    // TODO : extract this into a function, takes a vector of sourceIs, put them all in fifo in initial order
    // TODO : param to return the first loop it builds

    std::deque<std::pair<int, int>> edgeStack;
    // TODO : failure

    for(int iEdge = 0; iEdge < edges.size(); iEdge++)
        edgeStack.push_back({iEdge, -1});
    //edgeStack.push_back({0, -1}); edgeStack.push_back({1, -1}); edgeStack.push_back({2, -1});
    int iLoop = 0;

    while(edgeStack.size()) {

        auto thisEdge = edgeStack.front();
        edgeStack.pop_front();

        auto detectLoops = [&isEdgeInLoop, &edges, &eeply, &iLoop](int sourceI, int edgeI, bool returnOnSingleLoop, bool onlySourceLoops) { 

            bool shouldDebug = false;

            ControlPrintf printf;

            std::vector<std::vector<int>> loops;

            int iStep = 1;
            std::map<int, int> vertexStep;
            std::deque<int> fifo;

            fifo.push_back(sourceI);
            vertexStep[sourceI] = -1;

            if(edgeI != -1 && edgeI != sourceI) {
                fifo.push_back(edgeI);
                vertexStep[edgeI] = -1;
            }
      
            bool loopsDone = false;
            int nEdgesInLoop = 0;
            for(auto &el : isEdgeInLoop)
                if(el.second >= 0)
                    nEdgesInLoop++;

            printf("\n>>>>>> looping vertex %d to %d, %d edges in loops\n", sourceI, edgeI, nEdgesInLoop);

            while(fifo.size() && !loopsDone) {
                int currI = fifo.front();
                fifo.pop_front();
                vertexStep[currI] = iStep++;

                printf("\n >>> %d. vertex %d\n", iStep, currI);

            //      for(auto vs : vertexStep)
            //        printf("  %d (%d)\n", vs.first, vs.second);

                // iterate over vertex adjacencies
                for(int neighborI : edges[currI]) {
                    
                    printf("   edge to %d (%d)\n", neighborI, vertexStep[neighborI]);

                    if(vertexStep[neighborI] == 0) {
                    // don't put in fifo twice
                        vertexStep[neighborI] = -1;
                        fifo.push_back(neighborI);
                        printf("      push onto stack\n");
                    }
                    else if (vertexStep[neighborI] > 0) {
                      assert(vertexStep[neighborI] < vertexStep[currI]);

                      auto tracePath = [&edges, &vertexStep, &isEdgeInLoop](int currI, bool inReverse, int avoidI) {

//                        printf("trace_path_from %d(%d), reversed %d\n", currI, vertexStep[currI], inReverse);
                        std::vector<int> path;

                        int currStep = vertexStep[currI];

                        int iIter = 0;
                        while(currI != edges.size()) {
                        path.push_back(currI);

//                          printf(" %d. %d\n", iIter++, currI);

                        int bestI = edges.size();

                        for(auto nextI : edges[currI]) {

                            if(nextI == avoidI)
                              continue;

                            int loopId = inReverse ? isEdgeInLoop[{nextI, currI}] : isEdgeInLoop[{currI, nextI}];
//                            printf("   edge (%d - %d(%d)), loop %d\n", currI, nextI, vertexStep[nextI], loopId);

                            // walk on available edges towards the source 
                            if( (!inReverse && isEdgeInLoop[{currI, nextI}] == -1) || 
                                inReverse && isEdgeInLoop[{nextI, currI}] == -1) {

                            // by choosing a neighbor vertex stepped on the earliest
                            if(vertexStep[nextI] > 0 && vertexStep[nextI] < currStep) {
                                currStep = vertexStep[nextI];
                                bestI = nextI;
                            }
                            }
                        }
                        currI = bestI;
                        }

                        return path;
                    };

                    // trace path from neighborI back to sourceI by walking on revertsed edges that are !isEdgeInLoop
                    auto neighborPath = tracePath(neighborI, true, -1);

                    // detect path from currI back to sourceI by walking on edges that are !isEdgeInLoop

                    // TODO : make sure it doesn't go back through neighborI
                    auto currPath = tracePath(currI, false, neighborI);

                    std::map<int, std::pair<int, int>> vertexChainIdx;

                    printf("      back loop segment ");
                    for(int ii = 0; ii < currPath.size(); ii++) {
                        vertexChainIdx[currPath[ii]].first = ii+1;
                        printf("%d(%d) ", currPath[ii], vertexStep[currPath[ii]]);
                    }
                    printf("\n");

                    printf("      front loop segment ");
                    for(int jj = 0; jj < neighborPath.size(); jj++) {
                        printf("%d(%d) ", neighborPath[jj], vertexStep[neighborPath[jj]]);
                        vertexChainIdx[neighborPath[jj]].second = jj+1;
                    }

                    printf("\n");
                    int rootI = -1;
                    int rootIndex = 2*edges.size();

                    // find chain root vertex
                    for(auto vChain : vertexChainIdx)
                        if(vChain.second.first > 0 && vChain.second.second > 0) {

                        int vindex = vChain.second.first + vChain.second.second;

                        if(vindex < rootIndex) {
                            rootIndex = vindex;
                            rootI = vChain.first;
                        }
                    }

                    // check if paths form a loop
                    printf(" *** paths of length %lu, %lu, root %d\n", currPath.size(), neighborPath.size(), rootI);

                    std::vector<int> loop;

//                    printf(" %d, %d, %d\n", rootI, sourceI, !onlySourceLoops);

                    bool considerThisLoop = rootI == sourceI || !onlySourceLoops;

//                    printf(" consider this loop %d\n", considerThisLoop);

                    if(considerThisLoop) {

                        // TODO : check if the loop contains edge (sourceId, edgeId)
                        printf(" check for edge (%d, %d)\n", sourceI, edgeI);
                        bool edgeFound = edgeI == -1;
                        for(int ii = 0; !edgeFound && ii < currPath.size(); ii++) {
                            if( currPath[ii] == edgeI )
                                edgeFound = true;
                            if( currPath[ii] == rootI )
                                break;
                        }

                        for(int jj = 0; !edgeFound && jj < neighborPath.size(); jj++) {
                            if( neighborPath[jj] == edgeI )
                                edgeFound = true;
                            if( neighborPath[jj] == rootI )
                                break;
                        }

                        printf(" edge %sfound\n", edgeFound ? " " : "not ");
                        considerThisLoop = edgeFound;
                    }

                    if(considerThisLoop && rootI != -1 && rootI != currI && rootI != neighborI) {
                        printf("      ");

                        for(size_t ii = 0; ii < currPath.size(); ii++) {
                            if(currPath[ii] == rootI)
                                printf("|");
                
                            printf("%d", currPath[ii]);
                            loop.push_back(currPath[ii]);
                            if(currPath[ii] == rootI)
                                break;
                            printf(" ");
                        }

                        bool addThis = false;
                        for(size_t kk = 0; kk < neighborPath.size(); kk++) {
                                size_t jj = neighborPath.size()-kk-1;
                                if(addThis) {
                                    printf("%d ", neighborPath[jj]);
                                    loop.push_back(neighborPath[jj]);
                                }
                                else
                                    if(neighborPath[jj] == rootI) {
                                    addThis = true;
                                    printf("|");
                                }
                        }
                        printf("\n");

                        printf(" *** loop had %lu vertices\n", loop.size());

                        if(loop.size() >= 3) {
                            // make sure has at most one common edge with another loop

                            std::map<int, int> edgeLoopCount;

                            for(size_t ii = 0; ii < loop.size(); ii++) {
                                size_t nextI = (ii+1)%loop.size();

                                auto oppEdgeFound = isEdgeInLoop.find({loop[nextI], loop[ii]});

                                if(oppEdgeFound != isEdgeInLoop.end()) {
                                    int oppLoop = oppEdgeFound->second;
                                    printf("   edge %d, opposite loop %d\n", ii, oppLoop);
                                    if(oppLoop != -1)
                                        edgeLoopCount[oppLoop]++;
                                }
                            }

                            bool nmFound = false;
                            for(const auto &edgeCount : edgeLoopCount)
                                if(edgeCount.second > 1) {
                                    printf("  $$$ non-manifold $$$\n");
                                    nmFound = true; 
                                }

                            if(nmFound)
                                loop.clear();
                        }

                        if(loop.size() >= 3) {
                            // mark loop edges as used

                            eeply.StartFacet();

                            for(size_t ii = 0; ii < loop.size(); ii++) {
                                eeply.AddFacetVertex(loop[ii]);
                                size_t nextI = (ii+1)%loop.size();

                                auto edgeFound = isEdgeInLoop.find({loop[ii], loop[nextI]});

                                printf("  loop %lu: segment (%lu/%lu), add_edge (%d, %d) found %d, in loop %d\n", 
                                        loops.size(), ii, loop.size(), loop[ii], loop[nextI],
                                        edgeFound != isEdgeInLoop.end(), edgeFound != isEdgeInLoop.end() ? edgeFound->second: -1);

                //                assert(!isEdgeInLoop[{loop[ii], loop[nextI]}]);
                                if(isEdgeInLoop[{loop[ii], loop[nextI]}] >= 0)
                                    printf("!!!\n");
                                isEdgeInLoop[{loop[ii], loop[nextI]}] = iLoop;
                            }

                            eeply.EndFacet();

                            loops.push_back(std::move(loop));

                            iLoop++;
                            if(returnOnSingleLoop) {
                                loopsDone = true;
                                break;
                            }
                        }
                    }
                    printf("\n");
                };
            }
            printf("   %lu in stack\n", fifo.size());
            }

            return loops;
        };

        if(thisEdge.second == -1 || isEdgeInLoop[thisEdge] == -1) {

            auto sourceLoops = detectLoops(thisEdge.first, thisEdge.second, true, true);

//            printf("detectLoops : edge (%d, %d), %lu loops\n", thisEdge.first, thisEdge.second, sourceLoops.size());

            for(auto && sourceLoop : sourceLoops) {

                bool edgeFound = false;
                for(int ii = 0; ii < sourceLoop.size(); ii++) {
                    
//                    printf(" %d\n", sourceLoop[ii]);

                    int jj = (ii+1)%sourceLoop.size();

                    // get opposite edge
                    edgeStack.push_back({sourceLoop[jj], sourceLoop[ii]});

                    if(sourceLoop[ii] == thisEdge.first && sourceLoop[jj] == thisEdge.second) {
//                      printf("  edge\n");
                      edgeFound = true;
                    }
                    if(sourceLoop[jj] == thisEdge.first && sourceLoop[ii] == thisEdge.second) {
                      printf("  !!! back edge !!!\n");
                      edgeFound = true;
                    }
                }
            
                if(!edgeFound && sourceLoop.size())
                  printf("  !!! start edge not found in loop of %lu !!!\n", sourceLoop.size());
                loops.emplace_back(sourceLoop);
            }

//            if(loops.size() >= 50)
//                break;
        }
    }

    printf("%lu loops\n", loops.size());

    eeply.WritePLY("vorofaces.ply");

    // make a mesh, save

    MyMesh voroMesh;
    for(size_t iV = 0; iV < points.size(); iV++) {
      auto pnt = points[iV];

      voroMesh.add_vertex(MyMesh::Point(pnt[0], pnt[1], pnt[2]));
    }

    for(const auto &loop : loops) {
      std::vector<MyMesh::VertexHandle> handles;

      for(auto ii : loop)
        handles.push_back(MyMesh::VertexHandle(ii));

      voroMesh.add_face(handles);
    }

    printf("voro mesh has %d vertices, %d faces\n", voroMesh.n_vertices(), voroMesh.n_faces());

    return voroMesh;
}
  
// increase number of quad elements by finding triangles to pair up and remove common edge (quad diagonal)
void MeshSamplingRemoveDiagonals(MyMesh &voroMesh, const std::vector<Eigen::Vector3d> &points, MyMesh &mesh, const std::vector<Eigen::Vector3d> &crossField, const std::vector<std::set<size_t>> &uvAdjMatrix)
{
    // create quad dominant mesh by removind diagonals

    // TODO : the rest is a separate function QuadsByRemoveDiag
    // TODO : use crossField

    printf("%lu %lu\n", points.size(), uvAdjMatrix.size());

    AABBTree tree;
    tree.Build(mesh, 20);

    voroMesh.request_face_status();
    voroMesh.request_edge_status();
    voroMesh.request_vertex_status();

    auto numVertices = [&voroMesh](MyMesh::FaceHandle &facetHandle) {
        int nV = 0; 

        if(facetHandle.is_valid())
        for(MyMesh::FaceVertexIter fv_it = voroMesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it)
          nV++;
        
        return nV;
      };

    auto halfEdgeLen = [&voroMesh](MyMesh::HalfedgeHandle &hedgeHandle) {
        MyMesh::VertexHandle toVertex = voroMesh.to_vertex_handle(hedgeHandle);
        MyMesh::VertexHandle fromVertex = voroMesh.from_vertex_handle(hedgeHandle);

        auto toPnt = voroMesh.point(toVertex);
        auto fromPnt = voroMesh.point(fromVertex);

        return sqrt((toPnt-fromPnt).norm());
      };

    auto halfEdgeVec = [&voroMesh](MyMesh::HalfedgeHandle &hedgeHandle) {
        MyMesh::VertexHandle toVertex = voroMesh.to_vertex_handle(hedgeHandle);
        MyMesh::VertexHandle fromVertex = voroMesh.from_vertex_handle(hedgeHandle);

        auto toPnt = voroMesh.point(toVertex);
        auto fromPnt = voroMesh.point(fromVertex);

        return toVec(toPnt-fromPnt);
      };
    
    // TODO : iterate over edges of voro mesh
    // try to remove diagonal by testing the quad

    int nQuads = 0;

    std::vector<std::pair<double, MyMesh::HalfedgeHandle>> edgeLens;
    for(auto hedgeHandle : voroMesh.halfedges())
        edgeLens.push_back({halfEdgeLen(hedgeHandle), hedgeHandle});

    std::sort(edgeLens.begin(), edgeLens.end());
    std::reverse(edgeLens.begin(), edgeLens.end());

    OpenMesh::FPropHandleT<bool> wasMerged;
    voroMesh.add_property(wasMerged);

    for(auto facetHandle : voroMesh.faces())
      voroMesh.property(wasMerged, facetHandle) = false;

    happly::PLYExport eeeply;

    for(auto pnt : points)
        eeeply.AddVertex(pnt, 0.1, 0.8, 0.4);

    for(size_t iHedge = 0; iHedge < edgeLens.size(); iHedge++) {
      auto hedgeHandle = edgeLens[iHedge].second;

      if(!hedgeHandle.is_valid())
        continue;

      int iPnt = voroMesh.to_vertex_handle(hedgeHandle).idx();
      int jPnt = voroMesh.from_vertex_handle(hedgeHandle).idx();

//      if((iPnt == 1028 && jPnt == 1013 ) ||
//          (iPnt == 1013 && jPnt == 1028))
//        printf("opanki\n");

      if(iPnt < 0 || jPnt < 0 || iPnt >= uvAdjMatrix.size() || jPnt >= uvAdjMatrix.size())
        continue;

//      printf("eddge %d-%d\n", iPnt, jPnt);

      // check if diagonal doesn't appear in uvAdjMatrix
      if(uvAdjMatrix[iPnt].find(jPnt) != uvAdjMatrix[iPnt].end()) {
//        printf("eddge %d-%d\n", iPnt, jPnt);
        continue;
      }

//      auto currPoint = toVec(voroMesh.point(voroMesh.to_vertex_handle(hedgeHandle)));

//      auto np = tree.FindNearestPoint(currPoint);

      // TODO : extract into helper function
//      MyMesh::VertexHandle closestVertex;
//      if(np.facetIndex != -1) {
//        MyMesh::FaceHandle closestFacetHandle(np.facetIndex);

//        double closestDist = -1.0;

//        for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(closestFacetHandle); fv_it.is_valid(); ++fv_it) {
//          MyMesh::VertexHandle facetVertex(fv_it->idx());

//          double dist = (toVec(mesh.point(facetVertex))-currPoint).norm();

//          if(closestDist < 0.0 || dist < closestDist) {
//            closestDist = dist;
//            closestVertex = facetVertex;
//          }
//        }
//      }

//      if(closestVertex.is_valid()) {
//        auto normal = toVec(mesh.normal(closestVertex));
//        auto dirVec = halfEdgeVec(hedgeHandle);
//        dirVec.stableNormalize();

      // check against cross field
//        Eigen::Vector3d bestFitVec = bestRosyFit(dirVec, normal, crossField[closestVertex.idx()]);
//        if( dirVec.dot(bestFitVec) > 0.95 )
//          continue;
//      }

      auto facetHandle = voroMesh.face_handle(hedgeHandle);
      auto oppFacetHandle = voroMesh.opposite_face_handle(hedgeHandle);

      if(numVertices(facetHandle) == 3 && numVertices(oppFacetHandle) == 3 &&
        !voroMesh.property(wasMerged, facetHandle) && !voroMesh.property(wasMerged, oppFacetHandle)) {

        std::vector<MyMesh::HalfedgeHandle> hedgeHandles;

        std::vector<double> edgeLenss;

        for(MyMesh::FaceHalfedgeIter fh_it = voroMesh.fh_iter(facetHandle); fh_it.is_valid(); ++fh_it) {
          if(*fh_it != hedgeHandle) {
            auto hedgeHangle = *fh_it;
            hedgeHandles.push_back(hedgeHangle);
            double edgeLen = halfEdgeLen(hedgeHangle);
            edgeLenss.push_back(edgeLen);
          }
        }

        auto oppHedgeHandle = voroMesh.opposite_halfedge_handle(hedgeHandle);

        for(MyMesh::FaceHalfedgeIter fh_it = voroMesh.fh_iter(oppFacetHandle); fh_it.is_valid(); ++fh_it) {
          if(*fh_it != oppHedgeHandle) {
            auto hedgeHangle = *fh_it;
            hedgeHandles.push_back(hedgeHangle);
            double edgeLen = halfEdgeLen(hedgeHangle);
            edgeLenss.push_back(edgeLen);
          }
        }

        // check if all perspective quad edge lengths are roughly the same
        // TODO : check planarity
        // TODO : check diagonal lengths

        double averEdgeLen = 0;
        for(auto edgeLen : edgeLenss)
          averEdgeLen += edgeLen;

        averEdgeLen /= edgeLenss.size();

        bool oneOff = false;
//        for(auto edgeLen : edgeLenss) {

//          if( fabs(edgeLen - averEdgeLen) / averEdgeLen > 0.1)
//            oneOff = true;
//        }

        if(hedgeHandles.size() == 4 && !oneOff) {
          nQuads++;

          std::vector<MyMesh::VertexHandle> handles;

          // loop vertex handles from hedgeHandles

          MyMesh::HalfedgeHandle loopHandle = hedgeHandles.front();

          for(size_t iLoop = 0; iLoop < hedgeHandles.size(); iLoop++) {

            MyMesh::VertexHandle vertexHandle = voroMesh.from_vertex_handle(loopHandle);
            MyMesh::VertexHandle nextHandle = voroMesh.to_vertex_handle(loopHandle);

            handles.push_back(vertexHandle);

            for(auto nextHedgeHandle : hedgeHandles) {
              if(voroMesh.from_vertex_handle(nextHedgeHandle) == nextHandle) {
                loopHandle = nextHedgeHandle;
                break;
              }
            }
          }

          eeeply.StartFacet();

          for(size_t ii = 0; ii < handles.size(); ii++)
            eeeply.AddFacetVertex(handles[ii].idx());

          eeeply.EndFacet();

          voroMesh.property(wasMerged, facetHandle) = true;
          voroMesh.property(wasMerged, oppFacetHandle) = true;

          // delete two faces
          // add one quad

          if(facetHandle.is_valid())
            voroMesh.delete_face(facetHandle, true);
          
          if(oppFacetHandle.is_valid())
            voroMesh.delete_face(oppFacetHandle, true);

          voroMesh.add_face(handles);
        }
      }
    }

    for(auto facetHandle : voroMesh.faces())
      if( !voroMesh.property(wasMerged, facetHandle) ) {
          eeeply.StartFacet();

          for(MyMesh::FaceVertexIter fv_it = voroMesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it)
            eeeply.AddFacetVertex(fv_it->idx());

          eeeply.EndFacet();
      }

    eeeply.WritePLY("voroquads.ply");

    voroMesh.garbage_collection();
    OpenMesh::IO::write_mesh(voroMesh, "voroquadmesh.ply");

    voroMesh.remove_property(wasMerged);

    printf("%d quad diags to be removed\n", nQuads/2);
}

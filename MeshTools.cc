 #include <iostream>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <functional>
#include "happly.h"
#include "AABBTree.h"
#include "tools3d.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
// ----------------------------------------------------------------------------
 
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;
 
std::vector<double> ComputeEigenvector(MyMesh &mesh, const std::vector<std::vector<double>> &orthValsVec, int nIter) 
{
    // compute fiedler (associated with smallest positive eigenvalue) vector using iterative method
    // can be extended to compute n'th smallest positive eigenvalue vector if (n-1) previous eigenvectors are supplied
    std::vector<double> vals(mesh.n_vertices());
    std::vector<double> nextVals(mesh.n_vertices());

    double totalAver = 0.0;
    for(auto &val : vals) {
      val = (rand()-0.5)/RAND_MAX;
      totalAver += val;
    }

    totalAver /= vals.size();
    for(auto &val : vals) {
      val -= totalAver;
    }
  
    double prevNorm = mesh.n_vertices()/10.0;
    double devia = 0.0;

    for(int iIter = 0; iIter <= nIter; iIter++) {

      if(!(iIter %100) || iIter == nIter)
        printf("%d. %.6f ", iIter, devia);

      double totalAver = 0.0;
      int nVertices = 0;

      for (auto v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) {
        nVertices++;

        int valence = 0;
        
        int thisIdx = v_it->idx();

        double averVal = 0.0;

        for (auto vv_it=mesh.vv_iter( *v_it ); vv_it.is_valid(); ++vv_it) {
            int neighborIdx = vv_it->idx();

            averVal += vals[neighborIdx];

            ++valence;
        }

        if(valence)
          nextVals[thisIdx] = averVal/valence;

        totalAver += nextVals[thisIdx];
      }

      totalAver /= mesh.n_vertices();

      // move to the center
      for(auto &val : nextVals) {
        val -= totalAver;
      }

      // subtract other eigenvectors
      double normDot = 0.0;
      for(const auto &orthVals : orthValsVec) {
        double dotProd = 0.0;
        double orthNorm = 0.0;

        for(size_t ii = 0; ii < nextVals.size(); ii++) {
          dotProd += nextVals[ii]*orthVals[ii];
          orthNorm += orthVals[ii]*orthVals[ii];
        }

        if(fabs(dotProd) > 1e-10)
          for(size_t ii = 0; ii < nextVals.size(); ii++) {
            nextVals[ii] -= orthVals[ii]*(dotProd/orthNorm);
        }

        normDot += dotProd*dotProd;

        dotProd = 0.0;
        orthNorm = 0.0;

        for(size_t ii = 0; ii < nextVals.size(); ii++) {
          dotProd += nextVals[ii]*orthVals[ii];
          orthNorm += orthVals[ii]*orthVals[ii];
        }
      }

      // compute normal
      double vectorNorm = 0;
      for(auto &val : nextVals) {
        vectorNorm += val*val;
      }

      // normalize
      for(auto &val : nextVals ) {
        val *= sqrt(prevNorm/vectorNorm);
      }

      vectorNorm = sqrt(vectorNorm);
      if(!(iIter %100) || iIter == nIter) {
        if(orthValsVec.size())
          printf("%.6f, %.6f\n", vectorNorm, normDot);
        else
          printf("%.6f\n", vectorNorm);
      }

      devia = 0.0;
      for(size_t ii = 0; ii < vals.size(); ii++) {
        devia += (vals[ii] - nextVals[ii])*(vals[ii]-nextVals[ii]);
      }

      std::swap(vals, nextVals);
    }

    return vals;
}

std::vector<double> ComputeGradient(MyMesh &mesh, const std::vector<double> &scalar)
{
  // compute the absolute value of the gradient of scalar defined on mesh

  std::vector<double> gradientVal;

  for(auto vertexHandle : mesh.vertices()) {

    auto vpoint = toVec(mesh.point(vertexHandle));
    auto vnormal = toVec(mesh.normal(vertexHandle));

    // get edges    
    // 
    Eigen::Matrix3d mat;
    mat(0, 0) = 0.0;
    mat(0, 1) = 0.0;
    mat(0, 2) = 0.0;
    mat(1, 0) = 0.0;
    mat(1, 1) = 0.0;
    mat(1, 2) = 0.0;
    mat(2, 0) = 0.0;
    mat(2, 1) = 0.0;
    mat(2, 2) = 0.0;

    for(auto vh_iter = mesh.voh_iter(vertexHandle); vh_iter.is_valid(); ++vh_iter) {

      auto toVertex = mesh.to_vertex_handle(*vh_iter);

      auto vToPoint = toVec(mesh.point(toVertex));

      Eigen::Vector3d edgeVec = vToPoint-vpoint;

      double range = scalar[toVertex.idx()] - scalar[vertexHandle.idx()];

      // project onto the normal plane

      edgeVec -= vnormal*vnormal.dot(edgeVec);

      edgeVec *= range / edgeVec.dot(edgeVec);

      Eigen::Matrix3d outer = edgeVec*edgeVec.transpose();

      mat += outer;
    }

    // Create an EigenSolver object
    Eigen::EigenSolver<Eigen::Matrix3d> solver(mat);

    // Get the eigenvalues (as complex numbers)
    Eigen::VectorXcd eigenvalues = solver.eigenvalues();

//    std::cout << "Eigenvalues:\n" << eigenvalues << std::endl;

    auto grad = (eigenvalues[0]*eigenvalues[0] + eigenvalues[1]*eigenvalues[1] + eigenvalues[2]*eigenvalues[2]).real();

    gradientVal.push_back(grad);
  }

// TODO : extract into a function

  auto &vals = gradientVal;

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


  happly::PLYData plyOut;
  std::vector<std::array<double, 3>> meshVertexPositions;
  std::vector<std::array<double, 3>> meshVertexColors;
  std::vector<std::vector<size_t>> meshFaceIndices;

  meshVertexPositions.resize(mesh.n_vertices());
  meshVertexColors.resize(mesh.n_vertices());
  meshFaceIndices.resize(mesh.n_faces());

  for (auto v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it) {  
      int thisIdx = v_it->idx();

      auto pnt = mesh.point( *v_it );
      meshVertexPositions[thisIdx] = {pnt[0], pnt[1], pnt[2]};

      double weight = (vals[thisIdx]-0.5*averVal-minVal)/(averVal-minVal+1e-10);

      weight = std::max(0.0, weight);
      weight = std::min(1.0, weight);

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
  plyOut.write("gradient.ply", happly::DataFormat::ASCII);

  return gradientVal;
}

int ComputeOpaqueEdges(MyMesh &mesh, const std::vector<double> &vals, double val, OpenMesh::HPropHandleT<int> &opaqueEdges)
{
    // compute a closed chain of mesh edges best approximating the val isoline of vals defined on mesh vertices
    int nOpaque = 0;

    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::array<double, 3>> meshVertexColors;

    // Create an empty object

  // TODO : request_face_normals
  // TODO : request_vertex_colors, release_vertex_colors, has_vertex_colors

    // detect halfedges spanning over given value
    for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it) {

        auto toVertex = mesh.to_vertex_handle(*he_it);
        auto fromVertex = mesh.from_vertex_handle(*he_it);

        double valA = vals[toVertex.idx()]-val;
        double valB = vals[fromVertex.idx()]-val;

        if(valA*valB < 0.0 || fabs(valA) < 1e-10 || fabs(valB) < 1e-10) {
            mesh.property(opaqueEdges, *he_it) = -1;
            nOpaque++;
        }
        else
            mesh.property(opaqueEdges, *he_it) = -2;
    }

    printf("      %d opaque edges\n", nOpaque);

    // compute loops
    int nLoops = 0;

    for(auto he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it) {
      
      if(mesh.property(opaqueEdges, *he_it) == -1) {
        
        std::vector<MyMesh::HalfedgeHandle> stack;
        stack.push_back(*he_it);

        int nHedgesInLoop = 0;

        double r = 1.0*rand()/RAND_MAX;
        double g = 1.0*rand()/RAND_MAX;
        double b = 1.0*rand()/RAND_MAX;

        while(stack.size()) {

          auto hedgeHandle = stack.back();
          stack.pop_back();

          nHedgesInLoop++;
          mesh.property(opaqueEdges, hedgeHandle) = nLoops;
          mesh.property(opaqueEdges, mesh.opposite_halfedge_handle(hedgeHandle)) = nLoops;

          auto toVertex = mesh.to_vertex_handle(hedgeHandle);
          auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

          auto midPnt = (mesh.point(toVertex) + mesh.point(fromVertex))*0.5;
          meshVertexPositions.push_back({midPnt[0], midPnt[1], midPnt[2]});
          meshVertexColors.push_back({r, g, b});

          for(auto vh_iter = mesh.voh_iter(toVertex); vh_iter.is_valid(); ++vh_iter) {
            if(mesh.property(opaqueEdges, *vh_iter) == -1) {
              stack.push_back(*vh_iter);
            }
          }
        }

        printf("      loop %d has %d hedges\n", nLoops, nHedgesInLoop);
        nLoops++;
      } 
    }

 //   printf("%d loops found\n", nLoops);

///    plyOut.addVertexPositions(meshVertexPositions);
//    plyOut.addVertexColors(meshVertexColors);

    // Write the object to file
  //  plyOut.write("loops.ply", happly::DataFormat::ASCII);

    return nLoops;
}

std::vector<std::pair<MyMesh::VertexHandle, MyMesh::VertexHandle>> ComputeSharpEdges(MyMesh &mesh, OpenMesh::HPropHandleT<int> &curveIds, double dihedCos) 
{
  // detect sharp edges with dihedral angle cos below dihedCos, chain them into polylines (curves)
  std::vector<std::pair<MyMesh::VertexHandle, MyMesh::VertexHandle>> endpoints;

  int nCurves = 0;

  // detect dihedral edges
  for(auto hedgeHandle : mesh.halfedges()) {

    auto facetHandle = mesh.face_handle(hedgeHandle);
    auto oppHandle = mesh.opposite_face_handle(hedgeHandle);

    if(facetHandle.is_valid() && oppHandle.is_valid()) {

      auto norm1 = mesh.normal(facetHandle);
      auto norm2 = mesh.normal(oppHandle);

      double dihedDot = norm1.dot(norm2);

      if(dihedDot <= dihedCos + eps)
        mesh.property(curveIds, hedgeHandle) = -1;
      else
        mesh.property(curveIds, hedgeHandle) = -2;
    }
  }

  // group into curves
  for(auto hedgeHandle : mesh.halfedges()) {

    if(mesh.property(curveIds, hedgeHandle) == -1) {

      std::pair<MyMesh::VertexHandle, MyMesh::VertexHandle> curveEnds;

      std::vector<MyMesh::HalfedgeHandle> stack;

      mesh.property(curveIds, hedgeHandle) = nCurves;
      stack.push_back(hedgeHandle);

      auto oppHandle = mesh.opposite_halfedge_handle(hedgeHandle);
      mesh.property(curveIds, oppHandle) = nCurves;
      stack.push_back(oppHandle);

      while(stack.size()) {

        auto currHedge = stack.back();
        stack.pop_back();

        auto toVertex = mesh.to_vertex_handle(currHedge);

        int nSharpHedges = 0;
        std::vector<MyMesh::HalfedgeHandle> nextHedges;

        for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(toVertex); voh_it.is_valid(); ++voh_it) {
          auto nextHedge = *voh_it;

          int sharpId = mesh.property(curveIds, nextHedge);

          if(sharpId != -2) {
            nSharpHedges++;

            if(sharpId == -1)
              nextHedges.push_back(nextHedge);
          }

        }

        // count the number of sharp edges coming out of toVertex
        // if 2 and there is one unassigned to curve
        
        if(nSharpHedges == 2 && nextHedges.size() == 1) {
            auto nextHedge = nextHedges.back();
            mesh.property(curveIds, nextHedge) = nCurves;
            mesh.property(curveIds, mesh.opposite_halfedge_handle(nextHedge)) = nCurves;
            stack.push_back(nextHedge);
        }
        else if (nSharpHedges > 0) {
          // reached the end of a vertex
          if(curveEnds.first.is_valid())
            curveEnds.second = toVertex;
          else
            curveEnds.first = toVertex;
        }
      }

      if(!curveEnds.first.is_valid())
        printf("!!! no end point !!!\n");

      if(!curveEnds.second.is_valid()) // closed curve
        curveEnds.second = curveEnds.first;

      endpoints.push_back(curveEnds);
      nCurves++;
    }
  }

  printf("%d curves\n", nCurves);

  return endpoints;
}

std::vector<Eigen::Vector3d> SampleSharpEdges(MyMesh &mesh, OpenMesh::HPropHandleT<int> &curveIds, std::vector<std::pair<MyMesh::VertexHandle, MyMesh::VertexHandle>> &endpoints, double step) 
{
  // generate sampling points of the sharp edges, given step
  std::set<MyMesh::VertexHandle> endVertices;

  happly::PLYExport eply;

  std::vector<Eigen::Vector3d> pnts;
  // get curve edge chains
  size_t nCurves = endpoints.size();

  for(size_t iCurve = 0; iCurve < nCurves; iCurve++) {

    auto startVertex = endpoints[iCurve].first;
    auto endVertex = endpoints[iCurve].second;

    std::vector<MyMesh::HalfedgeHandle> hedgeChain;

    double chainLen = 0.0;

    MyMesh::VertexHandle currVertex;
    MyMesh::HalfedgeHandle currHedge;

    // trace chain from startVertex to endVertex

    while(currVertex != endVertex) {

      if(!currVertex.is_valid())
        currVertex = startVertex;

      MyMesh::HalfedgeHandle nextHedge;
      for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(currVertex); voh_it.is_valid(); ++voh_it) {
          auto testHedge = *voh_it;

          int sharpId = mesh.property(curveIds, testHedge);  

          // make sure it's not going back
          if(sharpId == static_cast<int>(iCurve) && (!currHedge.is_valid() || testHedge != mesh.opposite_halfedge_handle(currHedge))) {
            nextHedge = testHedge;
          }
      }

      if(nextHedge.is_valid()) {
        hedgeChain.push_back(nextHedge);
        currHedge = nextHedge;

        auto fromPnt = mesh.point(currVertex);

        currVertex = mesh.to_vertex_handle(currHedge);

        auto toPnt = mesh.point(currVertex);
        chainLen += (toPnt-fromPnt).norm();
      }
      else {
        printf("!!! broken chain !!!\n");
        break;
      }
    }

    printf("%lu. chain of %lu edges, %.6f length\n", iCurve, hedgeChain.size(), chainLen);

    if(chainLen > step) {
      int nInterPoints = static_cast<double>(chainLen/step);

      double chainStep = chainLen/nInterPoints;

      printf("  sample with %d pnts, %.6f step\n", nInterPoints, chainStep);

      double currLen = 0.0;
      double nextStepDist = chainStep;

      for(auto hedge : hedgeChain) {

        auto fromPnt = mesh.point(mesh.from_vertex_handle(hedge));
        auto toPnt = mesh.point(mesh.to_vertex_handle(hedge));

        auto edgeLen = (fromPnt-toPnt).norm();

        if(edgeLen > 0.0) {

          // check if intermediate point is between currLen and currLen + edgeLen

          bool tryOneStep = true;
          while(tryOneStep) {
            if(currLen + edgeLen > nextStepDist) {
              auto outputPoint = fromPnt + (toPnt-fromPnt)*(nextStepDist-currLen)/edgeLen;
              eply.AddVertex(outputPoint, 0.2, 0.4, 0.7);
              pnts.push_back(toVec(outputPoint));

              // TODO : outputPoint to ply file
              nextStepDist += chainStep;
            }
            else
              tryOneStep = false;
          }

          currLen += edgeLen;

        }
      }
    }

    // output vertex points, make sure each vertex is output exactly once

    auto insertIter = endVertices.insert(startVertex);
    if(insertIter.second) {
      eply.AddVertex(mesh.point(startVertex), 0.1, 0.5, 0.8);
      pnts.push_back(toVec(mesh.point(startVertex)));
    }

    auto insertIter1 = endVertices.insert(endVertex);
    if(insertIter1.second) {
      eply.AddVertex(mesh.point(endVertex), 0.1, 0.5, 0.8);
      pnts.push_back(toVec(mesh.point(endVertex)));
    }
  }

  eply.WritePLY("edgesamples.ply");

  return pnts;
}

void SetDirFieldFromSharpEdges(MyMesh &mesh, OpenMesh::HPropHandleT<int> &curveIds, OpenMesh::VPropHandleT<Eigen::Vector3d> &dirField)
{
  // set direction field to follow sharp edges
  for(auto vertexHandle : mesh.vertices()) {
    // get sharp edges of the vertex
    
    std::vector<MyMesh::HalfedgeHandle> sharpEdges;
    std::vector<Eigen::Vector3d> edgeVecs;

    for(auto vh_iter = mesh.voh_iter(vertexHandle); vh_iter.is_valid(); ++vh_iter) {
      if(mesh.property(curveIds, *vh_iter) >= 0) {
        sharpEdges.push_back(*vh_iter);

        Eigen::Vector3d edgeVec = toVec(mesh.point(mesh.to_vertex_handle(*vh_iter))-mesh.point(mesh.from_vertex_handle(*vh_iter)));
        edgeVec.stableNormalize();
        edgeVecs.push_back(edgeVec);
      }
    }

    // this is the precomputed field direction at mid-vertex of a chain of sharp edges
    if(edgeVecs.size() == 2) {
      Eigen::Vector3d vertexVec = edgeVecs[0]-edgeVecs[1];
      vertexVec.stableNormalize();
      mesh.property(dirField, vertexHandle) = vertexVec;
    }
  }
}

std::pair<Eigen::Vector3d, MyMesh::HalfedgeHandle> MakeStep(MyMesh &mesh, Eigen::Vector3d point, Eigen::Vector3d dirVec, MyMesh::HalfedgeHandle hedgeHandle, MyMesh::VertexHandle vertexHandle)
{
    // move along mesh surface from point to 
    Line3d shootSeg(point, point+dirVec);

    double exitParam = -1.0;
    Eigen::Vector3d exitPoint;
    MyMesh::HalfedgeHandle exitHedge;

    MyMesh::FaceHandle facetHandle;
    if(vertexHandle.is_valid()) {
        // circulate facets
        for(MyMesh::VertexIHalfedgeIter vh_iter = mesh.vih_iter(vertexHandle); vh_iter.is_valid(); ++vh_iter) {

          auto faceHandle = mesh.face_handle(*vh_iter);
          auto nextHe = mesh.next_halfedge_handle(*vh_iter);

          auto faceNormal = toVec(mesh.normal(faceHandle));

          auto prevPnt = toVec(mesh.point(mesh.from_vertex_handle(*vh_iter)));
          auto nextPnt = toVec(mesh.point(mesh.to_vertex_handle(nextHe)));
          auto currPnt = toVec(mesh.point(vertexHandle));

          auto toCross = (prevPnt-currPnt).cross(dirVec);
          auto fromCross = (currPnt-nextPnt).cross(dirVec);

//          printf("       facet %d\n", faceHandle.idx());
          if(toCross.dot(fromCross) > eps && faceNormal.dot(toCross) <- eps) {
            facetHandle = faceHandle;
//            printf("     >> take\n");
          }
        }
    }
    else {
      auto oppFaceHandle = mesh.opposite_face_handle(hedgeHandle);
      facetHandle = oppFaceHandle;
    }

    // face found

    // facet entrance currPnt
    // direction dirVec
    // find exit point

    // iterate over face half edges
    for(MyMesh::FaceHalfedgeIter fh_it = mesh.fh_iter(facetHandle); fh_it.is_valid(); ++fh_it) {
    // build line segment 

        auto toPnt = mesh.point(mesh.to_vertex_handle(*fh_it));
        auto fromPnt = mesh.point(mesh.from_vertex_handle(*fh_it));

        Line3d edgeSeg(toVec(toPnt), toVec(fromPnt));

        std::pair<double, double> params;
        bool isSkew = shootSeg.Skew(edgeSeg, &params);

//        printf("       (%.6f, %.6f, %.6f)-(%.6f, %.6f, %.6f)\n", toPnt[0], toPnt[1], toPnt[2], fromPnt[0], fromPnt[1], fromPnt[2]);

//        printf("        %d, %.6f, %.6f\n", isSkew, params.first, params.second);

        if(isSkew && params.first > eps && params.second >= 0.0 && params.second <= 1.0) {
            if(params.first > exitParam) {
                  exitParam = params.first;
                  exitPoint = edgeSeg.Eval(params.second);
                  exitHedge = *fh_it;
                  // TODO : check if there is exitVertex
              }
        }
    }

  return std::make_pair(exitPoint, exitHedge);
}

bool DoTrianglesIntersect(MyMesh &mesh, const std::vector<MyMesh::VertexHandle> &newTriangle, const std::vector <MyMesh::FaceHandle> & triangles)
{
  std::vector<std::array<double, 3>> triA;

  for(const auto &vertex : newTriangle) {
    auto point = mesh.point(vertex);
    triA.push_back({point[0], point[1], point[2]});
  }

  for(const auto &triangle : triangles) {

    std::vector<MyMesh::VertexHandle> vertices;
    std::vector<std::array<double, 3>> triB;

    for (auto fv_it = mesh.cfv_iter(triangle); fv_it.is_valid(); ++fv_it)
    {
      vertices.push_back(*fv_it);
      auto point = mesh.point(*fv_it);
      triB.push_back({point[0], point[1], point[2]});
    }

    // check if newTriangle and triangle have common vertices

    std::vector<MyMesh::VertexHandle> commonVertices;

    for(const auto vA : newTriangle)
      for(const auto vB : vertices)
        if(vA == vB)
          commonVertices.push_back(vA);

    Triangle3d triangleA(toVec(triA[0]), toVec(triA[1]), toVec(triA[2]));
    Triangle3d triangleB(toVec(triB[0]), toVec(triB[1]), toVec(triB[2]));

    // if has common edge
    if (commonVertices.size() == 2)
    {
      double dotProd = triangleA.GetNormal().dot(triangleB.GetNormal());

      if(dotProd < -0.999) {
        printf("!!! edge intersection !!!\n");
        return true;
      }

      return false;
    }
    //  if has common vertex
    else if (commonVertices.size() == 1)
    {
      auto commonVertex = commonVertices[0];
      // check if non-common edge intersects other facet


      // get non-common edge of newTriangle
      std::vector<MyMesh::VertexHandle> nonCommonA;
      for (const auto vA : newTriangle)
        if(vA != commonVertex)
          nonCommonA.push_back(vA);

      if (nonCommonA.size() != 2)
        printf("!!! wrong number of non common vertices !!!\n");

      // check if non-common edge intersects triangle
      {
      auto pntA = toVec(mesh.point(nonCommonA[0]));
      auto pntB = toVec(mesh.point(nonCommonA[1]));
      Line3d segA(pntA, pntB);

      std::pair<double, double> params = {-1, -1};
      int nX = triangleB.IntersectLine(segA, &params);

      std::vector<Eigen::Vector3d> testPoints;
      if (nX > 0 &&((0.0 <= params.first && params.first <= 1.0)))
        testPoints.push_back(segA.Eval(params.first));

      if (nX > 1 && (0.0 <= params.second && params.second <= 1.0))
        testPoints.push_back(segA.Eval(params.first));

      for (const auto &testPoint : testPoints)
        if (triangleA.IsPointInside(testPoint))
          {
            triangleA.Print();
            triangleB.Print();

            printf("common %d, non-common %d, %d\n", commonVertex.idx(), nonCommonA[0].idx(), nonCommonA[1].idx());

            printf("!!! vertex intersection !!!\n");
            getchar();
            return true;
          }
      }

      // get non-common edge of vertices
      std::vector<MyMesh::VertexHandle> nonCommonB;
      for (const auto vB : vertices)
        if (vB != commonVertex)
          nonCommonB.push_back(vB);

      if (nonCommonB.size() != 2)
        printf("!!! wrong number of non common vertices !!!\n");

        // TODO : a function to check intersection of triangle with segment
      {
      auto pntC = toVec(mesh.point(nonCommonB[0]));
      auto pntD = toVec(mesh.point(nonCommonB[1]));
      Line3d segB(pntC, pntD);

      std::pair<double, double> params = {-1, -1};
      int nX = triangleA.IntersectLine(segB, &params);

      std::vector<Eigen::Vector3d> testPoints;
      if (nX > 0 && ((0.0 <= params.first && params.first <= 1.0)))
        testPoints.push_back(segB.Eval(params.first));

      if( nX > 1 && (0.0 <= params.second && params.second <= 1.0))
        testPoints.push_back(segB.Eval(params.first));

      for(const auto &testPoint : testPoints)
        if (triangleA.IsPointInside(testPoint))
        {
          printf("%.6f %.6f %.6f 0 0 255\n", testPoint[0], testPoint[1], testPoint[2]);

          triangleA.Print();
          triangleB.Print();

          printf("common %d, non-common %d, %d\n", commonVertex.idx(), nonCommonB[0].idx(), nonCommonB[1].idx());

          printf("!!! vertex intersection !!!\n");
          getchar();
          return true;
        }
      }
      return false;
    }

    double V0[3] = { triA[0][0], triA[0][1], triA[0][2]};

    double V1[3] = { triA[1][0], triA[1][1], triA[1][2]};

    double V2[3] = { triA[2][0], triA[2][1], triA[2][2]};

    double U0[3] = { triB[0][0], triB[0][1], triB[0][2]};

    double U1[3] = { triB[1][0], triB[1][1], triB[1][2]};

    double U2[3] = { triB[2][0], triB[2][1], triB[2][2]};

    int doIntersect = tri_tri_intersect(V0, V1, V2, U0, U1, U2);

    if(doIntersect) {
      printf("!!! intersection !!!\n");
      return true;
    }

  }

  return false;
}

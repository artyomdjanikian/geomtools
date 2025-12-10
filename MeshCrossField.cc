#include "common.h"
#include "tools3d.h"
#include "happly.h"
#include "AABBTree.h"
#include "MeshTools.h"

std::vector<Eigen::Vector3d> ComputeCrossField(MyMesh &mesh, int nIter, std::vector<std::vector<MyMesh::HalfedgeHandle>> &topoCuts)
{
    // compute cross field, make it follow sharp edges and topoCuts
      
    // TODO : use const MyMesh iters

    mesh.request_vertex_normals();

    // build aabb tree
    AABBTree tree;
    tree.Build(mesh, 20);

    size_t nTopoCuts = topoCuts.size();
 
    printf("compute %lu smooth cuts\n", nTopoCuts);
    // smooth topo cuts

    auto __smoothLoops = [](MyMesh &mesh, const AABBTree &tree, std::vector<std::vector<MyMesh::HalfedgeHandle>> & topoCuts) {

      std::vector<std::vector<Eigen::Vector3d>> smoothLoops;

      happly::PLYExport epply;
      for (size_t iCut = 0, nTopoCuts = topoCuts.size(); iCut < nTopoCuts; iCut++) {

        std::map<MyMesh::VertexHandle, std::set<MyMesh::VertexHandle>> adjacencies;
        std::map<MyMesh::VertexHandle, Eigen::Vector3d> points;

        printf("%lu. %lu edges\n", iCut, topoCuts[iCut].size());

        for (auto &hedgeHandle : topoCuts[iCut])
        {
          auto toVertex = mesh.to_vertex_handle(hedgeHandle);
          auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

          points[toVertex] = toVec(mesh.point(toVertex));
          points[fromVertex] = toVec(mesh.point(fromVertex));
          adjacencies[toVertex].insert(fromVertex);
          adjacencies[fromVertex].insert(toVertex);
        }

        size_t nAdj = 0;
        for (const auto &adj : adjacencies)
          nAdj += adj.second.size();

        printf("%lu. %lu edges, %lu vertices, %lu adjacencies\n", iCut, topoCuts[iCut].size(), adjacencies.size(), nAdj);

        double r = 1.0 * rand() / RAND_MAX;
        double g = 1.0 * rand() / RAND_MAX;
        double b = 1.0 * rand() / RAND_MAX;

        double normalize = r + g + b;
        if (normalize > 0.0 && normalize < 1.0)
        {
          r /= normalize;
          g /= normalize;
          b /= normalize;
        }

        std::map<MyMesh::VertexHandle, Eigen::Vector3d> nextPoints;
        std::map<MyMesh::VertexHandle, int> projFacets;

        for (int iIter = 0; iIter < 1000; iIter++) {

          happly::PLYExport lepply;

          for (const auto &adj : adjacencies)
          {

            Eigen::Vector3d averPoint(0.0, 0.0, 0.0);

            for (const auto &neighbor : adj.second)
            {
              averPoint += points[neighbor];
            }

            if (adj.second.size())
            {
              averPoint *= 1.0 / adj.second.size();

              lepply.AddVertex(points[adj.first], 0.0, 0.0, 1.0);

              lepply.AddVertex(averPoint, 1.0, 0.0, 0.0);

              // project
              auto np = tree.FindNearestPoint(averPoint);

              lepply.AddVertex(np.pnt, 0.0, 1.0, 0.0);

              nextPoints[adj.first] = np.pnt;
              projFacets[adj.first] = np.facetIndex;
            }
          }

          //        if(iIter == 999) {
          //          printf("%g %g %g\n", r*255, g*255, b*255);
          //          lepply.WritePLY("lastiter.ply");
          //          getchar();
          //       }

          nextPoints.swap(points);
        }

        std::vector<Eigen::Vector3d> sortedLoop;

        std::set<MyMesh::VertexHandle> usedVertices;

        auto seed = adjacencies.begin();
        MyMesh::VertexHandle loopVertex = seed->first;

        usedVertices.insert(loopVertex);

        while (loopVertex.is_valid())
        {
          sortedLoop.push_back(points[loopVertex]);

          bool found = false;
          for (auto nextVertex : adjacencies[loopVertex])
          {
            auto insertIter = usedVertices.insert(nextVertex);
            if (insertIter.second)
            {
              loopVertex = nextVertex;
              found = true;
              break;
            }
          }

          if (!found)
            loopVertex = MyMesh::VertexHandle();
        }

        if (sortedLoop.size() != adjacencies.size())
          printf("!!! %lu %lu !!!\n", sortedLoop.size(), adjacencies.size());

        smoothLoops.push_back(std::move(sortedLoop));

        if (iCut < nTopoCuts)
        {
          b = 1.0;
          r = 0.0;
        }
        else
        {
          r = 1.0;
          b = 0.0;
        }

        for (const auto &pnt : points)
        {
          epply.AddVertex({pnt.second[0], pnt.second[1], pnt.second[2]}, r, g, b);
        }

        // TODO : iterate over edges
        // check two projFacets of two neighboring vertices, if facets share an edge
        // if not then insert another vertex in the middle

        // in the end, we have a contiguous band of facets
        // we can use the band to intersect two loops
        // TODO : establish loop direction at each vertex
      }

      epply.WritePLY("smoothloops.ply");

      return smoothLoops;
    };

    auto smoothLoops = __smoothLoops(mesh, tree, topoCuts);

    // TODO : mark vertices where smooth loops intersect with mesh

    OpenMesh::VPropHandleT<int> xLoops;
    mesh.add_property(xLoops);
    std::vector<std::vector<std::pair<int, Eigen::Vector3d>>> vertexXLoop;

    for(auto vertexHandle : mesh.vertices())
      mesh.property(xLoops, vertexHandle) = -1;

    // mesh vertex property xId, init to -1

    happly::PLYExport leply;

    for(size_t iLoop = 0; iLoop < smoothLoops.size(); iLoop++) {

      const auto &smoothLoop = smoothLoops[iLoop];
      // project points onto mesh using AABBtree
      // compute loop directions at each point

      std::vector<AABBTree::sRayX> projPoints(smoothLoop.size());
      std::vector<Eigen::Vector3d> loopDirs(smoothLoop.size());

      // iterate over proj points

      // TODO : #pragma omp parallel for
      for(int iPnt = 0, nPnts = smoothLoop.size(); iPnt < nPnts; iPnt++) {
        projPoints[iPnt] = tree.FindNearestPoint(smoothLoop[iPnt]);

        int iPrevPnt = (iPnt-1+nPnts)%nPnts;
        int iNextPnt = (iPnt+1)%nPnts;

        Eigen::Vector3d loopDir = smoothLoop[iPrevPnt]-smoothLoop[iNextPnt];
        loopDir.stableNormalize();

        MyMesh::FaceHandle facetHandle(projPoints[iPnt].facetIndex);

        // get vertices of projFacet
        for(MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle); fv_it.is_valid(); ++fv_it) {

          auto pnt = toVec(mesh.point(*fv_it));
          leply.AddVertex(pnt, 0.1, 0.5, 0.7);

          leply.AddVertex(pnt + loopDir*0.1);
          leply.AddVertex(pnt - loopDir*0.1);

          int xLoop = mesh.property(xLoops, *fv_it);

          if(xLoop == -1) {
            xLoop = vertexXLoop.size();
            vertexXLoop.resize(xLoop+1);

            mesh.property(xLoops, *fv_it) = xLoop;
          }

          vertexXLoop[xLoop].push_back(std::make_pair(iLoop, loopDir));
        }
      }
    }

    leply.WritePLY("loopvertices.ply");
    // TODO : detect loop intersections, angles

    int nLoopVertices = 0;
    int nXLoopVertices = 0;


    OpenMesh::VPropHandleT<Eigen::Vector3d> dirField;

    mesh.add_property(dirField);

    // TODO : init dirfield

    // TODO : function to set dirField from xLoops


    for(auto vertexHandle : mesh.vertices()) {
      mesh.property(dirField, vertexHandle) = Eigen::Vector3d(0.0, 0.0, 0.0);

      auto xLoop = mesh.property(xLoops, vertexHandle);

      if(xLoop >= 0) {
        nLoopVertices++;

        std::set<int> loopIds;

        bool useFirst = true;
        for(const auto &vxloop : vertexXLoop[xLoop]) {
          loopIds.insert(vxloop.first);
          if(useFirst) {

            auto dirVec = vxloop.second;
            Eigen::Vector3d vNormal = toVec(mesh.normal(vertexHandle));

            dirVec -= vNormal*vNormal.dot(dirVec);
            dirVec.stableNormalize();

            mesh.property(dirField, vertexHandle) = vxloop.second;
            useFirst = false;
          }
        }

        if(loopIds.size() > 1)
          nXLoopVertices++;
      }

    }

    printf("%d vertices on a loop, %d on loop intersection\n", nLoopVertices, nXLoopVertices);
    // TODO : perpendicularity test

    mesh.remove_property(xLoops);

    // TODO : return vertex directions, edge samples


    OpenMesh::HPropHandleT<int> curveIds;
    mesh.add_property(curveIds);

    auto endpoints = ComputeSharpEdges(mesh, curveIds, 0.5);
//    auto samplePnts = SampleSharpEdges(mesh, curveIds, endpoints, 2.0);

    // TODO : use samplePnts in MeshSamplingSampleUV

    SetDirFieldFromSharpEdges(mesh, curveIds, dirField);
    mesh.remove_property(curveIds);

    happly::PLYExport plyOut;

    for(size_t iIncludeCut = 0; iIncludeCut < topoCuts.size(); iIncludeCut++)
    {     
      double r = 1.0*rand()/RAND_MAX;
      double g = 1.0*rand()/RAND_MAX;
      double b = 1.0*rand()/RAND_MAX;

      double normalize = r+g+b;
      if(normalize > 0.0 && normalize < 1.0) {
        r /= normalize;
        g /= normalize;
        b /= normalize;
      }
  
      if(iIncludeCut < nTopoCuts) {
        b = 1.0;
        r = 0.0;
      }
      else {
        r = 1.0;
        b = 0.0;
      }

      for(auto &hedgeHandle : topoCuts[iIncludeCut]) {
            auto toVertex = mesh.to_vertex_handle(hedgeHandle);
            auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

            auto midPnt = (mesh.point(toVertex) + mesh.point(fromVertex))*0.5;

            plyOut.AddVertex(midPnt, r, g, b);
      }
    }

    plyOut.WritePLY("loops.ply");

    int nF = 0;
    for(auto vertexHandle : mesh.vertices()) {
      double norm = mesh.property(dirField,vertexHandle).norm();
      if(norm > 0.1)
        nF++;
    }

    printf("%d nonzero dirs\n", nF);

    std::vector<Eigen::Vector3d> crossField;
    
    SmoothCrossField(mesh, dirField, 200, 0, crossField);
    SmoothCrossField(mesh, dirField, 100, 1, crossField);
    for(int ii = 0; ii < 3 ;ii++) {
      SmoothCrossField(mesh, dirField, 1, 2, crossField);
      SmoothCrossField(mesh, dirField, 100, 1, crossField);
    }
  
//    getchar();

    mesh.remove_property(dirField);

    return crossField;
}

// TODO : input for constrained edges
void SmoothCrossField(MyMesh &mesh, OpenMesh::VPropHandleT<Eigen::Vector3d> &dirField, int nIter, int weightMethod, std::vector<Eigen::Vector3d> &crossField)
{
    int nF = 0;
    for (auto vertexHandle : mesh.vertices())
    {
        double norm = mesh.property(dirField, vertexHandle).norm();
        if (norm > 0.1)
            nF++;
    }

    printf("start SmoothCrossfield : %d nonzero dirs\n", nF);

    // compute vertex normals
    mesh.request_vertex_normals();

    bool shouldInit = crossField.size() == 0;
    crossField.resize(mesh.n_vertices());
    std::vector<Eigen::Vector3d> normals;

    // initialize cross field
#if 0  
  for(auto hedgeHandle : mesh.halfedges()) {
    auto toVertex = mesh.to_vertex_handle(hedgeHandle);
    auto fromVertex = mesh.from_vertex_handle(hedgeHandle);

    auto tpnt = mesh.point(toVertex); Eigen::Vector3d toPnt(tpnt[0], tpnt[1], tpnt[2]);
    auto fpnt = mesh.point(fromVertex); Eigen::Vector3d fromPnt(fpnt[0], fpnt[1], fpnt[2]);

    Eigen::Vector3d edgeVec = toPnt - fromPnt; 
    edgeVec.stableNormalize();

    if(edgeVec.dot(edgeVec) > eps) // TODO : check if constrained
      crossField[fromVertex.idx()] = edgeVec;num_
  }
#endif

    std::vector<bool> canMove(mesh.n_vertices(), true);
    int nFixed = 0;

    for (auto vertexHandle : mesh.vertices())
    {

        Eigen::Vector3d vec = mesh.property(dirField, vertexHandle);
        double len = vec.norm();

        if (len > 0.1)
        { // ititialize crossField
            nFixed++;
            canMove[vertexHandle.idx()] = false;
            if (shouldInit)
            {
                crossField[vertexHandle.idx()] = vec;
                crossField[vertexHandle.idx()].stableNormalize();
            }
        }
        else
        {
            auto vpoint = mesh.point(vertexHandle);
            auto vnormal = mesh.normal(vertexHandle);

            double vnlen = vnormal.norm();
            if (fabs(vnlen - 1.0) > 0.001)
                printf("!!!\n");

            // neighbor vertices

            double bestDot = -1.0;

            for (auto vh_iter = mesh.voh_iter(vertexHandle); vh_iter.is_valid(); ++vh_iter)
            {

                auto toVertex = mesh.to_vertex_handle(*vh_iter);
                auto topoint = mesh.point(toVertex);
                auto tonormal = mesh.normal(toVertex);

                double normDot = vnormal.dot(tonormal);

                if (shouldInit)
                {
                    if (normDot > bestDot)
                    {
                        bestDot = normDot;
                        crossField[vertexHandle.idx()] = toVec(vpoint - topoint);
                    }
                }
            }
        }
    }

    printf("%d fixed vertices\n", nFixed);
    //     mesh.update_vertex_normals();

    double totalOffset = 0.0;

    std::vector<Eigen::Vector3d> nextCrossField(mesh.n_vertices());
    std::vector<double> dotProds(mesh.n_vertices());

    for (int iIter = 0; iIter < nIter; iIter++)
    {
        // 1. zero out nextCrossField
        for (auto &c : nextCrossField)
            c = {0.0, 0.0, 0.0};

        // 2. average neighbor cross fields
        for (auto vertexHandle : mesh.vertices())
        {

            if (!canMove[vertexHandle.idx()])
            {
                nextCrossField[vertexHandle.idx()] = crossField[vertexHandle.idx()];
                continue;
            }

            auto vpoint = mesh.point(vertexHandle);
            auto vnormal = mesh.normal(vertexHandle);

            Eigen::Matrix3d mv = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d(vnormal[0], vnormal[1], vnormal[2])).toRotationMatrix();

            std::vector<Eigen::Vector3d> vertexRosy;

            vertexRosy.push_back(crossField[vertexHandle.idx()]);
            for (int iR = 1; iR < 4; iR++)
            {
                vertexRosy.push_back(mv * vertexRosy[iR - 1]);
            }

            dotProds[vertexHandle.idx()] = 0.0;
            double bestFitDot = -1.0;
            int nNeighbors = 0;

            // neighbor vertices
            for (auto vh_iter = mesh.voh_iter(vertexHandle); vh_iter.is_valid(); ++vh_iter)
            {

                auto toVertex = mesh.to_vertex_handle(*vh_iter);
                auto tonormal = mesh.normal(toVertex);

                Eigen::Matrix3d tomv = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d(tonormal[0], tonormal[1], tonormal[2])).toRotationMatrix();
                std::vector<Eigen::Vector3d> toVertexRosy;

                toVertexRosy.push_back(crossField[toVertex.idx()]);
                for (int iR = 1; iR < 4; iR++)
                {
                    toVertexRosy.push_back(tomv * toVertexRosy[iR - 1]);
                }

                double maxDot = -1.0;
                std::pair<size_t, size_t> bestRot;

                // get cross field of toVertex
                // find optimal rotation

                for (size_t ii = 0; ii < vertexRosy.size(); ii++)
                    for (size_t jj = 0; jj < toVertexRosy.size(); jj++)
                    {
                        double dotProd = vertexRosy[ii].dot(toVertexRosy[jj]);

                        if (dotProd > maxDot)
                        {
                            maxDot = dotProd;
                            bestRot = {ii, jj};
                        }
                    }

                dotProds[vertexHandle.idx()] += maxDot;
                nNeighbors++;

                auto &toCrossVec = toVertexRosy[bestRot.second];

                // rotate back
                if (bestRot.first > 0)
                {
                    if (bestRot.first <= 1)
                        toCrossVec = mv * toCrossVec;
                    if (bestRot.first <= 2)
                        toCrossVec = mv * toCrossVec;
                    if (bestRot.first <= 3)
                        toCrossVec = mv * toCrossVec;
                }

                if (maxDot < 0.0)
                    printf("!!!\n");

                // TODO : maybe average with weights, if poor angle then smaller weight

                double weight = 1.0;
                if (weightMethod == 1)
                    weight = maxDot;
                if (weightMethod == 2)
                {
                    if (maxDot > bestFitDot)
                    {
                        bestFitDot = maxDot;
                        nextCrossField[vertexHandle.idx()] = toCrossVec;
                    }
                }
                else
                    nextCrossField[vertexHandle.idx()] += toCrossVec * maxDot;
            }

            if (nNeighbors)
                dotProds[vertexHandle.idx()] /= nNeighbors;
        }

        // 3. project onto tangent plane
        double totalDot = 0.0;
        int nMoved = 0;
        //    dotProds.clear();

        for (auto vertexHandle : mesh.vertices())
            if (canMove[vertexHandle.idx()])
            {
                auto vnorm = mesh.normal(vertexHandle);
                auto pnt = mesh.point(vertexHandle);

                Eigen::Vector3d vnormal(vnorm[0], vnorm[1], vnorm[2]);

                double dotProd = vnormal.dot(nextCrossField[vertexHandle.idx()]);
                nextCrossField[vertexHandle.idx()] -= vnormal * dotProd;

                nextCrossField[vertexHandle.idx()].stableNormalize();

                if (nextCrossField[vertexHandle.idx()].norm() > eps)
                {
                    nMoved++;
                    //        dotProds.push_back(crossField[vertexHandle.idx()].dot(nextCrossField[vertexHandle.idx()]));
                    totalDot += crossField[vertexHandle.idx()].dot(nextCrossField[vertexHandle.idx()]);
                    crossField[vertexHandle.idx()] = nextCrossField[vertexHandle.idx()];
                }
            }

        if (iIter == 0 || iIter == nIter - 1 || !(iIter % 10))
            printf("%d. %.6f\n", iIter, nMoved ? totalDot / nMoved : 0.0);
    }

    happly::PLYExport eply;

    printf("crossfield computed\n");

    for (auto vertexHandle : mesh.vertices())
    {
        auto vnormal = mesh.normal(vertexHandle);
        auto pnt = mesh.point(vertexHandle);

        //    printf("   %.6f, %.6f, %.6f\n", vnormal[0], vnormal[1], vnormal[2]);

        Eigen::Vector3d vpoint(pnt[0], pnt[1], pnt[2]);

        auto cf = crossField[vertexHandle.idx()];

        Eigen::Matrix3d mv = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d(vnormal[0], vnormal[1], vnormal[2])).toRotationMatrix();

        std::vector<Eigen::Vector3d> vertexRosy;

        vertexRosy.push_back(crossField[vertexHandle.idx()]);
        for (int iR = 1; iR < 4; iR++)
        {
            vertexRosy.push_back(mv * vertexRosy[iR - 1]);
        }

        double weight = dotProds[vertexHandle.idx()];
        weight = 1.0 - 5.0 * (1.0 - weight);

        double red = 1.0 - weight;
        double blue = weight;
        double green = 4.0 * red * blue;

        if (!canMove[vertexHandle.idx()])
        {
            red = 0.5;
            blue = 0.5;
            green = 0.5;
        }

        eply.AddVertex({vpoint[0], vpoint[1], vpoint[2]}, 0.5 * red, 0.5 * green, 0.5 * blue);

        for (auto vr : vertexRosy)
        {
            Eigen::Vector3d vvr = vpoint + 0.1 * vr;

            eply.AddVertex({vvr[0], vvr[1], vvr[2]}, red, green, blue);
        }
    }

    eply.WritePLY("crossfield.ply");

    if (false)
    {
        // see how do mesh edges follow crossfield
        happly::PLYExport eply;

        for (auto vertexHandle : mesh.vertices())
        {

            auto point = toVec(mesh.point(vertexHandle));
            auto normal = toVec(mesh.normal(vertexHandle));

            for (auto vh_iter = mesh.voh_iter(vertexHandle); vh_iter.is_valid(); ++vh_iter)
            {

                auto toVertex = mesh.to_vertex_handle(*vh_iter);
                auto toPoint = toVec(mesh.point(toVertex));

                Eigen::Vector3d edgeVec = toPoint - point;
                edgeVec.stableNormalize();

                // TODO : fiond out what gets modified
                auto bestDirVec = bestRosyFit(edgeVec, normal, crossField[vertexHandle.idx()]);

                double dotProd = bestDirVec.dot(edgeVec);

                double w = (dotProd - 0.5) * 2.0;

                Eigen::Vector3d midPoint = point * 0.67 + toPoint * 0.33;

                double red = 1.0 - w;
                double blue = w;
                double green = 4.0 * red * blue;

                eply.AddVertex(midPoint, red, green, blue);
            }
        }

        eply.WritePLY("edgefits.ply");
    }
}

#include "common.h"
#include "tools3d.h"
#include "AABBTree.h"
#include "happly.h"
#include "MeshTools.h"
#include "MeshSampling.h"


/*    ../ build / Build / bin / SimpleCuts ~/ Documents / Models / topology - optimized - bottle - opener - 1.snapshot.2 / Bottle_Opener_remeshplus.ply */

void QuadRemeshTest(MyMesh &mesh, double step, int nIter)
{
      auto topoCuts = ComputeTopoCuts(mesh, nIter);
      auto crossField = ComputeCrossField(mesh, nIter, topoCuts);
      std::vector<std::set<size_t>> uvAdjMatrix;

      // TODO : SampleSharpEdges, use in MeshSamplingSampleUV

      OpenMesh::HPropHandleT<int> curveIds;
      mesh.add_property(curveIds);

      auto endpoints = ComputeSharpEdges(mesh, curveIds, 0.5);
      auto samplePnts = SampleSharpEdges(mesh, curveIds, endpoints, 2.0);

      mesh.remove_property(curveIds);

      auto uvpnts = MeshSamplingSampleUV(mesh, 2.0, 5, samplePnts, crossField, uvAdjMatrix);

      auto sourcePoints = uvpnts.first;
      auto points = uvpnts.second;

      auto ret = MeshSamplingCluster(mesh, sourcePoints, 0, crossField);

      auto sourceVertices = ret.first;
      auto regionMap = ret.second;

      //  MeshSamplingGeometricEdges(mesh, sourceVertices, points, crossField);

      auto edges = MeshSamplingPidEdges(mesh, sourceVertices, regionMap, crossField);

      MeshSamplingSave(points, edges, "uvvoropoints.ply");

      for (size_t iEdge = 0; iEdge < edges.size(); iEdge++)
      {

        std::set<size_t> result;

        std::set_intersection(uvAdjMatrix[iEdge].begin(), uvAdjMatrix[iEdge].end(),
                              edges[iEdge].begin(), edges[iEdge].end(),
                              std::inserter(result, result.begin()));

        uvAdjMatrix[iEdge].swap(result);
      }

      // TODO : intersection of edges[i] and uvAdjMatrix[i]
      //  auto edges = adjMatrix;

      // TODO : edge normals rotations best dot prod matrix, find 4 best fit cross edges
      //   MeshSamplingQuadQualityEdges(mesh, points, sourceVertices, regionMap, edges, crossField);
      //   MeshSamplingCrossEdges(mesh, points, sourceVertices, regionMap, edges, crossField);

      //  MeshSamplingTriQuads(mesh, sourceVertices, edges);
      auto voroMesh = MeshSamplingCycles(points, uvAdjMatrix);
      MeshSamplingRemoveDiagonals(voroMesh, points, mesh, crossField, uvAdjMatrix);
    }
 #include <iostream>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include "happly.h"
#include "fibonacci.hpp" 
#include "common.h"
#include "tools3d.h"

#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

// Type definitions
using Point3D = Eigen::Vector3d;
using DistanceMatrix = Eigen::MatrixXd;
using namespace boost;

// Define the graph using adjacency_list with edge weights
using Graph = adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double>>;
using Vertex = graph_traits<Graph>::vertex_descriptor;

// Dijxtra algorithm to find shortest path in a graph specified by vertex coordinates (points) and edge weights (dist_matrix)
std::vector<std::size_t> findShortestPath(
    const std::vector<Point3D>& points,
    const DistanceMatrix& dist_matrix,
    std::size_t source_index,
    std::size_t target_index) 
{
    std::size_t n = points.size();
    Graph g(n); // Create a graph with n vertices

    // Add edges based on the distance matrix
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            double dist = dist_matrix(i, j);
            if (std::isfinite(dist) && dist > 0.0) {
                add_edge(i, j, dist, g);
            }
        }
    }

    // Prepare containers for Dijkstra's output
    std::vector<Vertex> predecessors(n); // To store the shortest path tree
    std::vector<double> distances(n);    // To store distances from the source

    // Run Dijkstra's algorithm
    dijkstra_shortest_paths(g, source_index,
                            predecessor_map(&predecessors[0]).
                            distance_map(&distances[0]));

    // Reconstruct the path from source to target
    std::vector<std::size_t> path;
    for (Vertex v = target_index; v != source_index; v = predecessors[v]) {
        path.push_back(v);
        if (predecessors[v] == v) {
            std::cerr << "No path exists." << std::endl;
            return {};
        }
    }
    path.push_back(source_index);
    std::reverse(path.begin(), path.end());
    return path;
}

// ----------------------------------------------------------------------------
 
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

// Dijxtra algorith to find shortest distance along mesh edges 
// from source vertices to any other vertex
int DijxtraDistances(MyMesh &mesh,
                    const std::vector<MyMesh::VertexHandle> &sources, 
                    const std::function<bool(const MyMesh::HalfedgeHandle &)> &canWalkEdge,
                    const std::function<bool(const MyMesh::VertexHandle &)> &isFinished,
                    std::map<MyMesh::VertexHandle, std::pair<double, int>> &distMap,
                    FibonacciHeap<sVertexDistRec> * pPrioRTQ)
{
    FibonacciHeap<sVertexDistRec> localprioRTQ;
    FibonacciHeap<sVertexDistRec> &prioRTQ = pPrioRTQ == nullptr ? localprioRTQ : *pPrioRTQ;

    for(size_t iS = 0; iS < sources.size(); iS++)
        prioRTQ.decreaseKey({0.0, static_cast<int>(iS), sources[iS]});

    int iIter = 0;

    while(!prioRTQ.isEmpty()) {
        iIter++;

        // 1. extract vertex
        auto minVertex = prioRTQ.removeMinimum();
    
        int iSource = minVertex.iSource;
        double addDist = minVertex.dist;
        MyMesh::VertexHandle fromVertex = minVertex.vertexHandle;

        auto fromPos = mesh.point(fromVertex);

        distMap[fromVertex] = {addDist, iSource};

        if(isFinished(fromVertex)) {
//            printf("loop finished\n");
            break;
        }

        // 2.0 iterate over outgoing halfedges
        for(auto vh_iter = mesh.voh_iter(fromVertex); vh_iter.is_valid(); ++vh_iter) {
            if(canWalkEdge(*vh_iter)) {
                auto toVertex = mesh.to_vertex_handle(*vh_iter);
                // vertex wasn't extracted
                if( distMap.find(toVertex) == distMap.end() ) {
                    auto toPos = mesh.point(toVertex);

                    double edgeLength = 0.0;
                    for(int iC = 0; iC < 3; iC++)
                        edgeLength += (fromPos[iC]-toPos[iC])*(fromPos[iC]-toPos[iC]);

                    assert(edgeLength > 0.0);

                    edgeLength = sqrt(edgeLength);

                    prioRTQ.decreaseKey({addDist+edgeLength, iSource, toVertex});
                }
            }
        }
    }

//    printf("%d iters, %lu added\n", iIter, distMap.size());

    return iIter;
}

int DijxtraDistances(MyMesh &mesh,
                    const std::vector<MyMesh::VertexHandle> &sources,
                    const std::vector<std::array<Eigen::Vector3d, 3>> &localCoords,
                    const std::function<bool(const MyMesh::HalfedgeHandle &)> &canWalkEdge,
                    const std::function<bool(const MyMesh::VertexHandle &)> &isFinished,
                    std::map<MyMesh::VertexHandle, std::pair<double, int>> &distMap,
                    FibonacciHeap<sVertexDistRec> * pPrioRTQ)
{
    FibonacciHeap<sVertexDistRec> localprioRTQ;
    FibonacciHeap<sVertexDistRec> &prioRTQ = pPrioRTQ == nullptr ? localprioRTQ : *pPrioRTQ;

    for(size_t iS = 0; iS < sources.size(); iS++)
        prioRTQ.decreaseKey({0.0, static_cast<int>(iS), sources[iS]});

    int iIter = 0;

    while(!prioRTQ.isEmpty()) {
        iIter++;

        // 1. extract vertex
        auto minVertex = prioRTQ.removeMinimum();
    
        int iSource = minVertex.iSource;
        double addDist = minVertex.dist;
        MyMesh::VertexHandle fromVertex = minVertex.vertexHandle;

        auto fromPos = mesh.point(fromVertex);

        distMap[fromVertex] = {addDist, iSource};

        if(isFinished(fromVertex)) {
//            printf("loop finished\n");
            break;
        }

        // 2.0 iterate over outgoing halfedges
        for(auto vh_iter = mesh.voh_iter(fromVertex); vh_iter.is_valid(); ++vh_iter) {
            if(canWalkEdge(*vh_iter)) {
                auto toVertex = mesh.to_vertex_handle(*vh_iter);
                // vertex wasn't extracted
                if( distMap.find(toVertex) == distMap.end() ) {
                    auto toPos = mesh.point(toVertex);

                    double reachDist = 0.0;
                    if(localCoords.size()) {
                        Eigen::Vector3d uvw;

                        Eigen::Vector3d sourcePos = toVec(mesh.point(sources[iSource]));
                        Eigen::Vector3d toVertexVec = toVec(toPos) - sourcePos;

                        // find vertex local coords in the source coordinate system
                        for(size_t iC = 0; iC < 3; iC++)
                            uvw[iC] = fabs(localCoords[iSource][iC].dot(toVertexVec));

                        reachDist = std::max( std::max(uvw[0], uvw[1]), uvw[2]);
                    }
                    else {
                        double edgeLength = 0.0;
                        for(int iC = 0; iC < 3; iC++)
                            edgeLength += (fromPos[iC]-toPos[iC])*(fromPos[iC]-toPos[iC]);

                        assert(edgeLength > 0.0);

                        edgeLength = sqrt(edgeLength);

                        reachDist = addDist+edgeLength;
                    }

                    prioRTQ.decreaseKey({reachDist, iSource, toVertex});
                }
            }
        }
    }

//    printf("%d iters, %lu added\n", iIter, distMap.size());

    return iIter;
}

// compute shortest path from destination vertex using distance map computed by DijxtraDistances
std::vector<MyMesh::HalfedgeHandle> DijxtraPath(MyMesh &mesh, 
                                                const std::map<MyMesh::VertexHandle, std::pair<double, int>> &distMap,
                                                MyMesh::VertexHandle destVertex)
{
    std::vector<MyMesh::HalfedgeHandle> shortestPath;

    MyMesh::VertexHandle currVertex = destVertex;

    auto iter = distMap.find(currVertex);
    if(iter == distMap.end())
        return shortestPath;

    double dist = iter->second.first;

    while(dist > 0.0) {

        // find neighbor of curr vertex with the smallest distance to origin(s)
        MyMesh::VertexHandle nextVertex;
        MyMesh::HalfedgeHandle pathHedge;
        double nextVertexDist = dist;

        for(auto vh_iter = mesh.voh_iter(currVertex); vh_iter.is_valid(); ++vh_iter) {
            auto toVertex = mesh.to_vertex_handle(*vh_iter);
            auto iter = distMap.find(toVertex);
            if( iter != distMap.end() ) {
                double currDist = iter->second.first;

                if(currDist < nextVertexDist) {
                    nextVertexDist = currDist;
                    nextVertex = toVertex;
                    pathHedge = *vh_iter;
                }
            }
        }

        if(nextVertexDist < dist) {
            shortestPath.push_back(pathHedge);
            dist = nextVertexDist;
            currVertex = nextVertex;
        }
        else {
            // dead end, shouldn't happen
            printf("!!! dead end !!!\n");
            shortestPath.clear();
            break;
        }
    }

    return shortestPath;
}
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <functional>
#include <list>
#include <array>
#include <vector>
#include <map>

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
// -------------------- Eigen
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Dense>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

// Type definitions
using Point3D = Eigen::Vector3d;
using DistanceMatrix = Eigen::MatrixXd;

#include "fibonacci.hpp"

// TODO : find out why this doesn't work

//struct EigennTraits : OpenMesh::DefaultTraits {
//    using Point = Eigen::Vector3d;
//    using Normal = Eigen::Vector3d;
//    using TexCoord2D = Eigen::Vector2d;

// TODO : try these
//VertexAttributes(OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);
//FaceAttributes(OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);

//};

//using MyMesh = OpenMesh::PolyMesh_ArrayKernelT<EigennTraits>;

typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

int FloodFill(MyMesh &mesh, OpenMesh::HPropHandleT<int> &opaqueEdges, int iLoop, OpenMesh::FPropHandleT<int> &pids, bool shouldDebug = false);

std::vector<std::vector<MyMesh::HalfedgeHandle>> ComputeTopoCuts(MyMesh &mesh, int nIter);

void ComputeUniqueCuts(MyMesh &mesh, std::vector<std::vector<MyMesh::HalfedgeHandle>> &simpleCuts);
std::vector<std::vector<MyMesh::HalfedgeHandle>> ComputeCrossCuts(MyMesh &mesh, const std::vector<MyMesh::HalfedgeHandle> &opaqueEdgeVec);
bool IsCutSimple(MyMesh &mesh,
                 const std::vector<MyMesh::HalfedgeHandle> &opaqueEdgeVec);

void AssignPids(MyMesh &mesh, std::map<MyMesh::VertexHandle, std::pair<double, int>> &distMap, OpenMesh::FPropHandleT<int> &pids);

void FlipIsolatedFacets(MyMesh &mesh, OpenMesh::FPropHandleT<int> &pids);

std::vector<std::list<std::vector<MyMesh::HalfedgeHandle>>> ComputePatchBoundaries(MyMesh &mesh, const OpenMesh::FPropHandleT<int> &pids);

bool IsPatchTopologyInvalid(MyMesh &mesh, const OpenMesh::FPropHandleT<int> &pids,
                            const std::list<std::vector<MyMesh::HalfedgeHandle>> &pidBoundary,
                            std::vector<MyMesh::VertexHandle> &newSourceVertices);

void SavePLY(MyMesh &mesh, const OpenMesh::FPropHandleT<int> &pids, std::string fileName);
void SaveLoopsPLY(MyMesh &mesh, const std::vector<std::vector<MyMesh::HalfedgeHandle>> &crossCuts, std::string fileName);

struct sVertexDistRec {
        double dist;
        int    iSource;
        MyMesh::VertexHandle vertexHandle;

        bool operator==(const sVertexDistRec &rhs) const {
            return vertexHandle == rhs.vertexHandle;
        }

        bool operator <(const sVertexDistRec &rhs) const {
            return dist < rhs.dist;
        }

         bool operator >(const sVertexDistRec &rhs) const {
            return dist > rhs.dist;
        }       
};

int DijxtraDistances(MyMesh &mesh,
            const std::vector<MyMesh::VertexHandle> &sources, 
            const std::function<bool(const MyMesh::HalfedgeHandle &)> &canWalkEdge,
            const std::function<bool(const MyMesh::VertexHandle &)> &isFinished,
            std::map<MyMesh::VertexHandle, std::pair<double, int>> &distMap,
            FibonacciHeap<sVertexDistRec> * pPrioRTQ = nullptr);

int DijxtraDistances(MyMesh &mesh,
                    const std::vector<MyMesh::VertexHandle> &sources,
                    const std::vector<std::array<Eigen::Vector3d, 3>> &localCoords,
                    const std::function<bool(const MyMesh::HalfedgeHandle &)> &canWalkEdge,
                    const std::function<bool(const MyMesh::VertexHandle &)> &isFinished,
                    std::map<MyMesh::VertexHandle, std::pair<double, int>> &distMap,
                    FibonacciHeap<sVertexDistRec> * pPrioRTQ);

std::vector<std::size_t> findShortestPath(const std::vector<Point3D>& points, 
                                          const DistanceMatrix& dist_matrix, 
                                          std::size_t source_index,
                                          std::size_t target_index);

void SimulateCloth(MyMesh &mesh);

std::vector<MyMesh::HalfedgeHandle> DijxtraPath(MyMesh &mesh, 
                                                const std::map<MyMesh::VertexHandle, std::pair<double, int>> &distMap,
                                                MyMesh::VertexHandle destVertex);


std::vector<Eigen::Vector3d> ComputeCrossField(MyMesh &mesh, int nIter, std::vector<std::vector<MyMesh::HalfedgeHandle>> &topoCuts);

void SmoothCrossField(MyMesh &mesh, OpenMesh::VPropHandleT<Eigen::Vector3d> &dirField, int nIter, int weightMethod, std::vector<Eigen::Vector3d> &crossfield);

void PatchLoops(MyMesh &mesh, double loopStep);

void SimulateSandCastle(MyMesh &mesh, double gridStep, int nIter);

class AABBTree;

void ShapeDiameterFuncTest(MyMesh &mesh, AABBTree &tree, int nIter);

void QuadRemeshTest(MyMesh &mesh, double step, int nIter);

void OctreeTest(MyMesh &mesh, const AABBTree &tree, int nIter);

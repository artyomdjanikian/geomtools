#pragma once
  
#include "common.h"
#include <vector>

std::pair<std::vector<MyMesh::VertexHandle>, std::vector<Eigen::Vector3d>> MeshSamplingSampleUV(MyMesh &mesh, double step, int nIter, std::vector<Eigen::Vector3d> &initPoints, std::vector<Eigen::Vector3d> &crossField, std::vector<std::set<size_t>> &adjMatrix);

std::vector<MyMesh::HalfedgeHandle> MeshSamplingNeighborhood(MyMesh &mesh, MyMesh::VertexHandle &sourceVertex, double rad);
Eigen::Vector3d MeshSamplingNeighbor(MyMesh &mesh, Eigen::Vector3d sourcePnt, double rad, Eigen::Vector3d normal, Eigen::Vector3d dir, const std::vector<MyMesh::HalfedgeHandle> &neighborhood);

std::vector<MyMesh::VertexHandle> MeshSamplingSampleRandom(MyMesh &mesh, double step, int nIter, std::vector<Eigen::Vector3d> &crossField);

std::pair<std::vector<MyMesh::VertexHandle>, std::vector<int>> MeshSamplingCluster(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, int nIter, std::vector<Eigen::Vector3d> &crossField);
std::vector<std::set<size_t>> MeshSamplingEdges(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<Eigen::Vector3d> &crossField);
std::vector<std::set<size_t>> MeshSamplingPidEdges(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<Eigen::Vector3d> &crossField);
//std::vector<std::set<size_t>> MeshSamplingGeometricEdges(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, std::vector<Eigen::Vector3d> &points, std::vector<Eigen::Vector3d> &crossField);

std::vector<std::set<size_t>> MeshSamplingEdges(const std::vector<Eigen::Vector3d> &points);


//void MeshSamplingQuadQualityEdges(MyMesh &mesh, std::vector<Eigen::Vector3d> &points, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<std::set<size_t>> &adjMatrix, std::vector<Eigen::Vector3d> &crossField);
//void MeshSamplingCrossEdges(MyMesh &mesh, std::vector<Eigen::Vector3d> &points, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<std::set<size_t>> &adjMatrix, std::vector<Eigen::Vector3d> &crossField);

std::vector<Eigen::Vector3d> MeshSamplingSmooth(MyMesh &mesh, std::vector<MyMesh::VertexHandle> &sourceVertices, const std::vector<int> &regionMap, std::vector<std::set<size_t>> &adjMatrix, std::vector<Eigen::Vector3d> &crossField);

void MeshSamplingTriQuads(MyMesh &mesh, const std::vector<MyMesh::VertexHandle> &sourceVertices, std::vector<std::set<size_t>> &edges);
MyMesh MeshSamplingCycles(const std::vector<Eigen::Vector3d> &points, std::vector<std::set<size_t>> &edges);
void MeshSamplingRemoveDiagonals(MyMesh &voroMesh, const std::vector<Eigen::Vector3d> &points, MyMesh &mesh, const std::vector<Eigen::Vector3d> &crossField, const std::vector<std::set<size_t>> &uvAdjMatrix);
void MeshSamplingSave(std::vector<Eigen::Vector3d> &meshPoints, std::vector<std::set<size_t>> &adjMatrix, std::string fileName);

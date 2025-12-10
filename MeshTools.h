#pragma once
#include <vector>
#include "common.h"

std::vector<double> ComputeEigenvector(MyMesh &mesh, const std::vector<std::vector<double>> &orthValsVec, int nIter);
std::vector<double> ComputeGradient(MyMesh &mesh, const std::vector<double> &scalar);
std::pair<Eigen::Vector3d, MyMesh::HalfedgeHandle> MakeStep(MyMesh &mesh, Eigen::Vector3d point, Eigen::Vector3d dir, MyMesh::HalfedgeHandle halfedgeHandle, MyMesh::VertexHandle vertexHandle);
int ComputeOpaqueEdges(MyMesh &mesh, const std::vector<double> &vals, double val, OpenMesh::HPropHandleT<int> &opaqueEdges);
std::vector<std::pair<MyMesh::VertexHandle, MyMesh::VertexHandle>> ComputeSharpEdges(MyMesh &mesh, OpenMesh::HPropHandleT<int> &curveIds, double dihedCos);
std::vector<Eigen::Vector3d> SampleSharpEdges(MyMesh &mesh, OpenMesh::HPropHandleT<int> &curveIds, std::vector<std::pair<MyMesh::VertexHandle, MyMesh::VertexHandle>> &endpoints, double step);
void SetDirFieldFromSharpEdges(MyMesh &mesh, OpenMesh::HPropHandleT<int> &curveIds, OpenMesh::VPropHandleT<Eigen::Vector3d> &dirField);
// TODO : void SetDirFieldFromLoops(MyMesh &mesh, OpenMesh::VPropHandleT<int> &xLoops, OpenMesh::VPropHandleT<Eigen::Vector3d> &dirField);
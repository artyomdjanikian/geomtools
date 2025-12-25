// point set tools
// similar to point cloud but allowed to be sparce and/or have interior points

#pragma once
#include <vector>
#include <eigen3/Eigen/Core>

std::vector<Eigen::Vector3d> PointSetComputeNormals(const std::vector<Eigen::Vector3d> &points, size_t kNearest);
std::vector<double> PointSetGetExteriorWeights(const std::vector<Eigen::Vector3d> &points, size_t kNearest);
std::vector<Eigen::Vector3d> PointSetCluster(const std::vector<Eigen::Vector3d> &points, size_t nClusters, int nClusterIters);
void PointSetRemoveDuplicates(std::vector<Eigen::Vector3d> &points, double eps);
void PointSetSmooth(std::vector<Eigen::Vector3d> &points, int kNearest, int nIter);// auto weights = PointSetGetExteriorWeights(points, 20);
void PointSetSmooth(std::vector<Eigen::Vector3d> &points, const std::vector<bool> &pinnedPoints, int kNearest, int nIter); // auto weights = PointSetGetExteriorWeights(points, 20);

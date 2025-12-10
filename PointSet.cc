#include "PointSet.h"
#include "happly.h"
#include "tools3d.h"
#include <eigen3/Eigen/Dense>

// for each element:
// find nNearest nearest neighbors
// build covariance matrix as sum of (point[i]-point[neighborOfI]).transpose()*(point[i]-point[neighborOfI])
// find eigenvalues, ratio of smallest to largest eigenvalues tells how likely the point to lie
// on the point set boundary
std::vector<double> PointSetGetExteriorWeights(const std::vector<Eigen::Vector3d> &points, size_t kNearest)
{
    // min/max eigenvalue of covariance matrix
    std::vector<double> weights(points.size(), 0.0);

    // TODO : optimize with cSpatialHashGrid
    for (size_t iPnt = 0; iPnt < points.size(); iPnt++)
    {
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

        std::vector<std::pair<double, Eigen::Vector3d>> nearestNeighbors;

        for (size_t jPnt = 0; jPnt < points.size(); jPnt++) {
            if (iPnt != jPnt) {
                Eigen::Vector3d edgeVec = points[jPnt] - points[iPnt];
                double len = edgeVec.norm();
                nearestNeighbors.push_back({len, edgeVec});
            }
        }

        std::sort(nearestNeighbors.begin(), nearestNeighbors.end(),
                [](const auto &a, const auto &b) { return a.first < b.first; });

        if (kNearest > 0 && kNearest < nearestNeighbors.size())
            nearestNeighbors.resize(kNearest);

        happly::PLYExport ePly;

        for (size_t jPnt = 0; jPnt < nearestNeighbors.size(); jPnt++) {
            Eigen::Vector3d edgeVec = nearestNeighbors[jPnt].second;

            ePly.AddVertex(points[iPnt] + edgeVec, 0.0, 0.0, 1.0);

            Eigen::Matrix3d outer = edgeVec * edgeVec.transpose();
            mat += outer;
        }

        ePly.AddVertex(points[iPnt], 1.0, 0.0, 0.0);

        // Create an EigenSolver object
        Eigen::EigenSolver<Eigen::Matrix3d> solver(mat);

        // Get the eigenvalues (as complex numbers)
        Eigen::VectorXcd eigenvalues = solver.eigenvalues();

        Eigen::MatrixXcd eigenvectors = solver.eigenvectors();

        std::vector<Plane3d> planes;

        for (long ii = 0; ii < eigenvalues.size(); ii++) {
        //      printf("eigenvalue %.6f\n", eigenvalues[ii].real());
            auto evec = eigenvectors.col(ii);

            Eigen::Vector3d revec(evec[0].real(), evec[1].real(), evec[2].real());
            ePly.AddVertex(revec + points[iPnt], 0.0, 1.0, 0.0);

            planes.push_back(Plane3d(points[iPnt], revec));
        }

        Plane3d bestPlane = planes[0];
        double bestRatio = 0.0;

        for (const auto &plane : planes) {
            double totalDist = 0.0;
            double totalSignedDist = 0.0;

            for (size_t jPnt = 0; jPnt < nearestNeighbors.size(); jPnt++) {
                Eigen::Vector3d edgeVec = nearestNeighbors[jPnt].second;

                double dist = plane.SignedDistance(points[iPnt] + edgeVec);
                //        printf("dist = %.6f\n", dist);

                totalDist += fabs(dist);
                totalSignedDist += dist;
            }

            bool flipPlane = false;
            if (totalSignedDist < 0.0)
                flipPlane = true;

            totalSignedDist = fabs(totalSignedDist);
            double thisRatio = totalDist > eps ? totalSignedDist / totalDist : 0.0;

            //      printf("dist %.6f, signed %.6f, this ratio = %.6f\n", totalDist, totalSignedDist, thisRatio);

            if (bestRatio == 0.0 || thisRatio > bestRatio) {
                bestRatio = std::max(bestRatio, thisRatio);
                bestPlane = plane;
                if (flipPlane)
                    bestPlane.Flip();
            }
        }

        // find the halfspace containing only points[iPnt], use bestPlane as an estimate
        for (int iTry = 0; iTry < 10; iTry++) {

            int nNeg = 0;
            Eigen::Vector3d newNorm(0.0, 0.0, 0.0);
            double errDot = 0.0;

            double totalDist = 0.0;
            double totalSignedDist = 0.0;

            for (size_t jPnt = 0; jPnt < nearestNeighbors.size(); jPnt++)
            {
                Eigen::Vector3d edgeVec = nearestNeighbors[jPnt].second;

                auto testPnt = points[iPnt] + edgeVec;
                double dist = bestPlane.SignedDistance(testPnt);

                totalDist += fabs(dist);
                totalSignedDist += dist;

                if (dist < 0.0)
                {
                Eigen::Vector3d pntVec = -edgeVec;

                Eigen::Vector3d interNorm = -edgeVec.cross(bestPlane.GetNormal());

                Eigen::Vector3d thisNorm = edgeVec.cross(interNorm);

                thisNorm.stableNormalize();

                if (thisNorm.norm() > eps)
                {

                    double thisDot = thisNorm.dot(bestPlane.GetNormal());
                    printf("   (%.6f ,%.6f, %.6f), dist %.6f\n", testPnt[0], testPnt[1], testPnt[2], dist);
                    printf("    (%.6f, %.6f, %.6f), dot %.6f\n", thisNorm[0], thisNorm[1], thisNorm[2], thisDot);

                    nNeg++;
                    newNorm += thisNorm;
                    errDot += 1.0 - thisDot * thisDot;
                }
                }
            }

            double thisRatio = totalDist > eps ? totalSignedDist / totalDist : 0.0;

        //      printf("dist %.6f, signed %.6f, this ratio = %.6f\n", totalDist, totalSignedDist, thisRatio);

            if (thisRatio > bestRatio)
                bestRatio = std::max(bestRatio, thisRatio);

            auto plnormal = bestPlane.GetNormal();

            printf("%d. (%.6f, %.6f, %.6f), %d neg, %.6f errDot\n\n", iTry, plnormal[0], plnormal[1], plnormal[2], nNeg, errDot);
            if (nNeg == 0)
                break;

            newNorm.stableNormalize();
            bestPlane = Plane3d(points[iPnt], newNorm);

        //        printf("dist = %.6f\n", dist);
        }

        ePly.AddVertex(bestPlane.GetNormal() + points[iPnt], 1.0, 0.0, 1.0);

        //    ePly.WritePLY("knearest.ply");

        //    getchar();

        //    printf("best ratio = %.6f\n", bestRatio);

        //     getchar();

        //     double minVal = eigenvalues[0].real();
        //     double maxVal = eigenvalues[0].real();

        //     for(auto val : eigenvalues) {
        //         double realval = val.real();

        //         minVal = std::min(minVal, realval);
        //         maxVal = std::max(maxVal, realval);
        //     }

        weights[iPnt] = 1.0 - bestRatio; // fabs(minVal/maxVal);
    }

    return weights;
}

// cluster point set
std::vector<Eigen::Vector3d> PointSetCluster(const std::vector<Eigen::Vector3d> &points, size_t clusterSize, int nClusterIters)
{
    std::vector<Eigen::Vector3d> clusterPoints;

    // can't be less than 1 point per cluster, can't be less than 1 cluster
    clusterSize = std::min(clusterSize, points.size());
    if (clusterSize == 0)
        clusterSize = 1;

    for (size_t ii = 0; ii < points.size() / clusterSize; ii++)
        clusterPoints.push_back(points[ii]);

    printf("cluster %lu points\n", clusterPoints.size());

    // debug only - show points with exterior weights
    // the higher the weight, the more likely the point to lie in the point set interior
    {
        auto weights = PointSetGetExteriorWeights(points, 20);
        happly::PLYExport qPly;

        for (size_t iPnt = 0; iPnt < points.size(); iPnt++)
        {
            const auto &point = points[iPnt];
            qPly.AddVertex(point, weights[iPnt], 4.0 * weights[iPnt] * (1.0 - weights[iPnt]), (1.0 - weights[iPnt]));
        }

        qPly.WritePLY("toclusterpoints.ply");
    }

    // Lloyd's clustering algorithm
    for (int iIter = 0; iIter < nClusterIters; iIter++)
    {

        happly::PLYExport qPly;

        std::vector<std::pair<Eigen::Vector3d, int>> newClusterPoints(clusterPoints.size(), {Eigen::Vector3d::Zero(), 0});
        std::vector<double> clusterSizes(clusterPoints.size());

        std::vector<double> colorWeights;

        for (size_t iCluster = 0; iCluster < clusterPoints.size(); iCluster++)
        {
            colorWeights.push_back(1.0 * rand() / RAND_MAX);
        }

        // for each point, find the closest cluster

        for (const auto &point : points)
        {
            double minDist = -1.0;
            size_t iBestCluster = 0;

            for (size_t iCluster = 0; iCluster < clusterPoints.size(); iCluster++)
            {

                double clusterDist = (clusterPoints[iCluster] - point).norm();

                if (minDist < 0.0 || clusterDist < minDist)
                {
                    minDist = clusterDist;
                    iBestCluster = iCluster;
                }
            }

            double weight = colorWeights[iBestCluster];

            //        qPly.AddVertex(point, 0.75*weight, 2.0*(1.0-weight)*weight, 0.75*(1.0-weight));

            newClusterPoints[iBestCluster].first += point;
            newClusterPoints[iBestCluster].second++;
            clusterSizes[iBestCluster] += minDist;
        }

        for (size_t iCluster = 0; iCluster < clusterPoints.size(); iCluster++)
        {
            double weight = colorWeights[iCluster];

            qPly.AddVertex(clusterPoints[iCluster], weight, 0.5 * (1.0 - weight) * weight, (1.0 - weight));
        }

        // set the cluster point to the centroid of points in the cluster
        for (size_t iCluster = 0; iCluster < clusterPoints.size(); iCluster++)
        {
            if (newClusterPoints[iCluster].second)
                clusterPoints[iCluster] = newClusterPoints[iCluster].first * (1.0 / newClusterPoints[iCluster].second);
        }

        for (size_t iCluster = 0; iCluster < clusterPoints.size(); iCluster++)
        {
            double weight = colorWeights[iCluster];

            qPly.AddVertex(clusterPoints[iCluster], weight, 4.0 * (1.0 - weight) * weight, (1.0 - weight));
        }

        qPly.WritePLY("patchclusters.ply");

        ///      getchar();
    }

    return clusterPoints;
}

// remove points closer than eps
void PointSetRemoveDuplicates(std::vector<Eigen::Vector3d> &points, double eps)
{
    std::sort(points.begin(), points.end(), [](const auto &pnt1, const auto &pnt2) { 
      if(pnt1[0] < pnt2[0])
        return true;
      else if (pnt1[0] == pnt2[0]) {
        if(pnt1[1] < pnt2[1])
          return true;
        else if (pnt1[1] == pnt2[1] && pnt1[2] < pnt2[2])
          return true;
      }

      return false; });

    size_t jPnt = 0;
    for (size_t iPnt = 0; iPnt < points.size(); iPnt++)
    {
        bool coincident = false;
        if (iPnt > 0)
            coincident = (points[iPnt] - points[iPnt - 1]).norm() < 1e-5;

        if (!coincident)
        {
            if (jPnt < iPnt)
                points[jPnt] = points[iPnt];
            jPnt++;
        }
    }

    if (jPnt != points.size())
    {
        printf("%lu points now\n", jPnt);
        points.resize(jPnt);
    }
}

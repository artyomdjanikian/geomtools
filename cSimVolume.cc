#include <stdint.h>
#include <fstream>
#include <vector>
#include <array>
#include <math.h>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include "common.h"

#include "cSimVolume.h"
#include "tools3d.h"
#include "happly.h"
#include "AABBTree.h"

void cSimVolume::initAnisotropic(const Bounds3d &bounds, double _step)
{
  step = _step;
  Eigen::Vector3d span = bounds.maxBound - bounds.minBound;
  printf("initial span %.6f, %.6f, %.6f\n", span[0], span[1], span[2]);

  Bounds3d extendBounds = bounds;
  extendBounds += extendBounds.minBound-span*0.5;
  extendBounds += extendBounds.maxBound+span*0.5;

  auto origin = extendBounds.minBound;
  span = extendBounds.maxBound - extendBounds.minBound;

  Ni = static_cast<size_t>(ceil(span[0]/step));
  Nj = static_cast<size_t>(ceil(span[1]/step));
  Nk = static_cast<size_t>(ceil(span[2]/step));

  printf("extended span %.6f, %.6f, %.6f\n", span[0], span[1], span[2]);
  printf("%d * %d * %d simulation grid\n", Ni, Nj, Nk);

  simVol.resize(Ni*Nj*Nk);

  // fill out
  for(size_t i = 0; i < Ni; i++)
    for(size_t j = 0; j < Nj; j++)
      for(size_t k = 0; k < Nk; k++) {
          auto &mass = getMass(i, j, k);

          mass.ijk = {static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)};

          mass.pos = origin + Eigen::Vector3d(step*i, step*j, step*k);
          mass.origPos = mass.pos;
          mass.nextPos = mass.pos;
          mass.prevPos = mass.pos; 

          mass.invMass = 0.0;
          mass.force = {0.0, 0.0, 0.0};
      }
}

void cSimVolume::init(const Bounds3d &bounds, double extendFactor)
{
      simVol.resize(Ni*Nj*Nk);

      Eigen::Vector3d span = bounds.maxBound - bounds.minBound;
      auto extend = span * extendFactor;
      span += extend*2.0;
      auto origin = bounds.minBound - extend;

      step = std::max(span[0], span[1]);
      step = std::max(step, span[2]);
      step /= (N-1);

      // fill out
      for(size_t i = 0; i < Ni; i++)
        for(size_t j = 0; j < Nj; j++)
          for(size_t k = 0; k < Nk; k++) {
              auto &mass = getMass(i, j, k);

//              bool isFilled = (k>=N/2) &&(j < N/10 || k > N-N/10);
//              bool isFixed = isFilled && (j == 0 || k == N/2);
//              bool isForced = !isFixed; //j >= N/2;

              mass.ijk = {static_cast<int>(i), static_cast<int>(j), static_cast<int>(k)};

              //mass.pos = origin;
              mass.pos = origin + Eigen::Vector3d(step*i, step*j, step*k);
              mass.origPos = mass.pos;
              mass.nextPos = mass.pos;
              mass.prevPos = mass.pos; 

//              printf("%d, %d, %d => %.6f, %.6f, %.6f\n", i, j, k, mass.pos[0], mass.pos[1], mass.pos[2]);

//              if(isFilled)
//                mass.invMass = 1.0;
//              else
                mass.invMass = 0.0;

//              if(isFixed)
//                mass.invMass = 0.0;

//              if(isForced)
//                mass.force = {0.0, 0.0, -0.02};
//              else
                mass.force = {0.0, 0.0, 0.0};
          }

}

void cSimVolume::savePLY(std::string filename, int type)
{
  std::fstream plyfile;

  plyfile.open(filename.c_str(), std::ios::out);

  std::vector<double> weights;

  int nParticles = 0;
  for(size_t i = 0; i < Ni; i++)
    for(size_t j = 0; j < Nj; j++)
      for(size_t k = 0; k < Nk; k++) {
        auto &mass = getMass(i, j, k);

        if(mass.invMass >= 0.0) {

            // type 0 : set up
            // type 4/5 : get force
            // type 1 : vom mises stress
            // type 2 : get average strain
            // type 3:  get elastic energy

            double weight = -1.0;

            if(type == 0) {
              Eigen::Vector3d forceVec = getForce(mass);

              if(mass.invMass <= eps)
                weight = 0.0;
              else if (mass.force.norm() > eps)
                weight = 1.0;
              else
                weight = 0.5;

            }
            else if(type == 4 || type == 5) {
              Eigen::Vector3d forceVec = getForce(mass);

              if(type == 5)
                forceVec += mass.force;

              weight = forceVec.norm();
//              printf("%.6f\n", weight);
            }
            else if (type == 2) {
              weight = getStrain(mass);
            }
            else if (type == 3) {
              auto strensor = getStrainTensor(mass);
              auto stressor = getStressTensor(mass);

              for(int ii = 0; ii < 3; ii++)
                for(int jj = 0; jj < 3; jj++)
                  weight += strensor(ii, jj) * stressor(ii, jj);


//                std::vector<double> eivals = getStrainEigenValues(mass);
//                weight = eivals[0] + eivals[1] + eivals[2];
            }
            else if (type == 1) {
//              auto strensor = getStrainTensor(mass);
              auto stressor = getStressTensor(mass);

//              printf("\n");

//              printf("strain tensor\n");

              auto strensor = getStrainTensor(mass);
//              for(int ii = 0; ii < 3; ii++) {
//                printf("(%.6f, %.6f, %.6f)", strensor(ii, 0), strensor(ii, 1), strensor(ii, 2));
//                if(ii == 2)
//                  printf("\n");
//                else 
//                  printf(", ");
//             }

  
//              printf("stress tensor\n");
//              for(int ii = 0; ii < 3; ii++) {
//                printf("(%.6f, %.6f, %.6f)", stressor(ii, 0), stressor(ii, 1), stressor(ii, 2));
//                if(ii == 2)
//                  printf("\n");
//                else 
//                  printf(", ");
//              }

              std::vector<double> eivals;

              Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(stressor);
              for(int i = 0; i < 3; ++i) {
                auto E = eigensolver.eigenvalues().col(0)[i];

                eivals.push_back(E.real());
              }

//              printf("eigenvalues %.6f, %.6f, %.6f\n", eivals[0], eivals[1], eivals[2]);

              weight = 0.0;
              for(int ii = 0; ii < 3; ii++) {
                int jj = ii+1;
                if(jj == 3)
                  jj = 0;

                weight += (eivals[ii]-eivals[jj])*(eivals[ii]-eivals[jj]);
              }

              weight = sqrt(0.5*weight);
            }

            weights.push_back(weight);
            nParticles++;
          }
        }

  plyfile << "ply\n";
  plyfile << "format ascii 1.0\n";

  plyfile << "element vertex " << nParticles << "\n";
  plyfile << "property double x\n";
  plyfile << "property double y\n";
  plyfile << "property double z\n";
  plyfile << "property uint8 red\n";
  plyfile << "property uint8 green\n";
  plyfile << "property uint8 blue\n";

  plyfile << "element face 0\n";
  plyfile << "property list uint8 int32 vertex_indices\n";
  plyfile << "end_header\n";

  double maxWeight = 0.0;

  for(auto weight : weights)
    maxWeight = std::max(weight, maxWeight);
 
  if(maxWeight < 1e-10)
    maxWeight = 1.0;

//  printf("savePLY : max weight = %.6f\n", maxWeight);

  int iP = 0;

  for(size_t i = 0; i < Ni; i++)
    for(size_t j = 0; j < Nj; j++)
      for(size_t k = 0; k < Nk; k++) {
          auto &mass = getMass(i, j, k);

          if(mass.invMass >= 0.0) {

            weights[iP] /= maxWeight;

//            printf("%d.  %.6f\n", iP, weights[iP]);

            plyfile << mass.pos[0] << " " << mass.pos[1] << " " << mass.pos[2] << " ";

            uint8_t red = static_cast<uint8_t>(255*weights[iP]);  // mass.invMass == 0.0 ? 255 : 120;

            uint8_t green = static_cast<uint8_t>(255 * 4* (weights[iP])*(1.0-weights[iP]));

            uint8_t blue = static_cast<uint8_t>(255 * (1.0-weights[iP]));

            iP++;

            //if(fabs(mass.force[0])+fabs(mass.force[1])+fabs(mass.force[2]) > 0.0) {
            //    red = 250;
            //    green = 50;
            //    blue = 250;
           // }
            //else
            if (mass.invMass < 0.0) {
                red = 0;
                green = 0;
                blue = 0;
            }

            plyfile << std::to_string(red) << " " << std::to_string(green) << " " << std::to_string(blue) << "\n";
        }
      }

  plyfile << "\n";
}

Eigen::MatrixXd cSimVolume::getStrainTensor(const sMassPoint &mass) const
{
    int i = mass.ijk[0];
    int j = mass.ijk[1];
    int k = mass.ijk[2];

    Eigen::Vector3d dUdX(0.0, 0.0, 0.0);
    int ndUdX = 0;

    Eigen::Vector3d dUdY(0.0, 0.0, 0.0);
    int ndUdY = 0;

    Eigen::Vector3d dUdZ(0.0, 0.0, 0.0);
    int ndUdZ = 0;

    if(isValid(i+1, j, k)) {
      const auto &neighbor = getMass(i+1, j, k);

      dUdX += neighbor.pos-mass.pos;
      ndUdX++;
    }

    if(isValid(i-1, j, k)) {
      const auto &neighbor = getMass(i-1, j, k);

      dUdX -= neighbor.pos-mass.pos;
      ndUdX++;
    }

    if(isValid(i, j+1, k)) {
      const auto &neighbor = getMass(i, j+1, k);

      dUdY += neighbor.pos-mass.pos;
      ndUdY++;
    }

    if(isValid(i, j-1, k)) {
      const auto &neighbor = getMass(i, j-1, k);

      dUdY -= neighbor.pos-mass.pos;
      ndUdY++;
    }

    if(isValid(i, j, k+1)) {
      const auto &neighbor = getMass(i, j, k+1);

      dUdZ += neighbor.pos-mass.pos;
      ndUdZ++;
    }

    if(isValid(i, j, k-1)) {
      const auto &neighbor = getMass(i, j, k-1);

      dUdZ -= neighbor.pos-mass.pos;
      ndUdZ++;
    }

    dUdX *= 1.0/ndUdX;
    dUdY *= 1.0/ndUdY;
    dUdZ *= 1.0/ndUdZ;

    double strainTensor[3][3];

    strainTensor[0][0] = dUdX[0];
    strainTensor[1][1] = dUdY[1];
    strainTensor[2][2] = dUdZ[2];
    
    strainTensor[0][1] = 0.5*(dUdX[1] + dUdY[0]);
    strainTensor[1][0] = strainTensor[0][1];

    strainTensor[1][2] = 0.5*(dUdZ[1] + dUdY[2]);
    strainTensor[2][1] = strainTensor[1][2];

    strainTensor[2][0] = 0.5*(dUdX[2] + dUdZ[0]);
    strainTensor[0][2] = strainTensor[2][0];

    Eigen::MatrixXd strensor(3, 3);

    for(int ii = 0; ii < 3 ; ii++)
      for(int jj = 0; jj < 3; jj++)
        strensor(ii, jj) = strainTensor[ii][jj];

    return strensor;
}

Eigen::MatrixXd cSimVolume::getStressTensor(const sMassPoint &mass) const
{
    int i = mass.ijk[0];
    int j = mass.ijk[1];
    int k = mass.ijk[2];

    Eigen::MatrixXd strensor = getStrainTensor(mass);

    Eigen::MatrixXd strensor6(1, 6);

    strensor6(0, 0) = strensor(0, 0);
    strensor6(0, 1) = strensor(1, 1);
    strensor6(0, 2) = strensor(2, 2);
    strensor6(0, 3) = strensor(1, 2);
    strensor6(0, 4) = strensor(2, 0);
    strensor6(0, 5) = strensor(0, 1);

    // TODO : from strensor to strensor6
    
    double poissonR = 0.25;

    Eigen::MatrixXd elastensor(6, 6);

    // https://doc.comsol.com/5.5/doc/com.comsol.help.sme/sme_ug_theory.06.23.html
    // isotropic elasticity tensor

    //double EmodYoung = 1.0;
    // double factor = EmodYoung/((1.0+poissonR)*(1.0-2.0*poissonR));

    for(int ii = 0; ii < 6; ii++)
      for(int jj = 0; jj < 6; jj++) {
        elastensor(ii, jj) = 0.0;

        if( ii == jj ) {
          if(ii < 3)
            elastensor(ii, jj) = 1.0-poissonR;
          else
            elastensor(ii, jj) = (1.0-2.0*poissonR)/2.0;
        }
        else if (ii < 3 && jj < 3)
          elastensor(ii, jj) = poissonR;
        else
          elastensor(ii, jj) = 0.0;
      }

    auto stressor6 = elastensor*strensor6;

    Eigen::MatrixXd stressor(3, 3);

    // from stressor6 to stressor

    stressor(0, 0) = stressor6(0, 0);
    stressor(1, 1) = stressor6(0, 1);
    stressor(2, 2) = stressor6(0, 2);
    stressor(1, 2) = stressor6(0, 3);
    stressor(2, 1) = stressor6(0, 3);
    stressor(2, 0) = stressor6(0, 4);
    stressor(0, 2) = stressor6(0, 4);
    stressor(0, 1) = stressor6(0, 5);
    stressor(1, 0) = stressor6(0, 5);

    return stressor;
}

double cSimVolume::getStrain(const sMassPoint &mass) const
{
  double totalStrain = 0.0;

  auto strensor = getStrainTensor(mass);

  for(int ii = 0; ii < 3; ii++)
    for(int jj = 0; jj < 3; jj++)
      totalStrain+= strensor(ii, jj) * strensor(ii, jj);

  return totalStrain;
}

std::vector<double> cSimVolume::getStrainEigenValues(const sMassPoint &mass) const
{
  auto strensor = getStrainTensor(mass);

  std::vector<double> eivals;

  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(strensor);
  for(int i = 0; i < 3; ++i) {
    auto E = eigensolver.eigenvalues().col(0)[i];

    eivals.push_back(E.real());
//    printf("%.6f", E);
//   if(i == 2)  
//      printf("\n");
//    else  
//      printf(", ");
  }

  return eivals;
}

Eigen::Vector3d cSimVolume::getForce(const sMassPoint &mass) const
{
  bool shouldDebug = false;//mass.ijk[0] == 1 && mass.ijk[1] == 5 && mass.ijk[2] == 5;

  Eigen::Vector3d totalShift(0.0, 0.0, 0.0);

  if (mass.invMass > 0.0) {
    if(shouldDebug)
      printf("\ncurr pos (%.6f, %.6f, %.6f)\n", mass.pos[0], mass.pos[1], mass.pos[2]);

    int nNeighbors = 0;

    int i = mass.ijk[0];
    int j = mass.ijk[1];
    int k = mass.ijk[2];

    for (int ii = -1; ii <= 1; ii++)
      for (int jj = -1; jj <= 1; jj++)
        for (int kk = -1; kk <= 1; kk++)
        {
          int nZeros = 0;

          if (ii == 0)
            nZeros++;
          if (jj == 0)
            nZeros++;
          if (kk == 0)
            nZeros++;

          if (nZeros == 3) // || nZeros == 1)
            continue;

          if (isValid(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk))
          {
            auto &neighbor = getMass(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk);

            if(neighbor.invMass < 0.0)
              continue;

            if(shouldDebug)
              printf("(%d, %d, %d)+(%d, %d, %d)\n", i, j, k, ii, jj, kk);

            Eigen::Vector3d shiftVec = neighbor.pos - mass.pos;
            // compute length from

            double shiftLen = shiftVec.norm();

            if (shiftLen > eps)
            {
              double relaxLengthFactor = nZeros == 2 ? 1 : nZeros == 1 ? sqrt(2.0) : nZeros == 0 ? sqrt(3.0) : 0.0;
              double relaxLength = step*relaxLengthFactor;

              double stiffnessFactor = nZeros == 2 ? 1.0 : nZeros == 1 ? 0.5 : nZeros == 0 ? 0.1875 : 0.0;

              if(shouldDebug)
                printf("   (%d, %d, %d), relaxed len %.6f, curr len %.6f\n", ii, jj, kk, relaxLength, shiftLen);

              shiftVec *= stiffnessFactor * (shiftLen - relaxLength) / shiftLen;
              if(shouldDebug)
                printf("       force (%.6f, %.6f, %.6f)\n", shiftVec[0], shiftVec[1], shiftVec[2]);
              totalShift += shiftVec;
            }
          }
        }
  }

  return totalShift;//+ mass.force;
}

void cSimVolume::simulate(int nIter, double dT)
{
  // inverse mass proportional to number of occupied cells
//    Eigen::Vector3d pnt0(0.0, 0.0, 0.0);
//    Eigen::Vector3d pnt1(1.0, 1.0, 1.0);

//    Line3d line(pnt0, pnt1);
//    Eigen::Vector3d testPnt(0.6, 0.4, 0.5);

//    double ddd = line.SegmentDistanceSq(testPnt);
//    exit(0);

    for(size_t i = 0; i < Ni; i++)
      for(size_t j = 0; j < Nj; j++)
        for(size_t k = 0; k < Nk; k++) {

          auto &mass = getMass(i, j, k);
          mass.pos = mass.origPos;
          mass.nextPos = mass.origPos;
          mass.prevPos = mass.origPos;

          if(mass.invMass > 0.0) {
            int nNeighbors = 0;

            for (int ii = -1; ii <= 1; ii++)
                for (int jj = -1; jj <= 1; jj++)
                  for (int kk = -1; kk <= 1; kk++) {
                    int nZeros = 0;

                    if (ii == 0)
                      nZeros++;
                    if (jj == 0)
                      nZeros++;
                    if (kk == 0)
                      nZeros++;

                    if (nZeros != 0) {

                      if (isValid(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk)) {
                        auto &neighbor = getMass(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk);

                        if(neighbor.invMass >= 0.0) {

                          if( nZeros == 2 ) 
                            nNeighbors += 4;
                          else if (nZeros == 1)
                            nNeighbors += 2;
                          else if (nZeros == 0)
                            nNeighbors += 1;
                        }
                      }
                    }
                  }

            if(nNeighbors < 56)
              mass.invMass *= nNeighbors/56.0;

          }
        }

  std::vector<std::array<size_t, 3>> masses;

  double bottomZ = 0.0;
  bool oneFound = false;

  for(size_t i = 0; i < Ni; i++)
    for(size_t j = 0; j < Nj; j++)
      for(size_t k = 0; k < Nk; k++) {
//      printf("%d, %d, %d\n", i, j, k);

        auto &mass = getMass(i, j, k);

        if(mass.invMass > 0.0) {
          masses.push_back({i, j, k});
          if(!oneFound) {
            oneFound = true;
            bottomZ = mass.pos[2];
          }
          else
            bottomZ = std::min(bottomZ, mass.pos[2]);
        }
      }

  printf("bottom level at %.6f\n", bottomZ);

  printf("%d iterations, %lu particles, %lu movable\n", 
          nIter, N*N*N, masses.size());

  for(int iIter = 0; iIter < nIter; iIter++) {
#pragma omp parallel for
   for(const auto &ijk : masses) {
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      auto &mass = getMass(i, j, k);

      if(mass.invMass > 0.0) {

        Eigen::Vector3d totalForce = getForce(mass)+mass.force;

        double decr = 0.95;
        mass.nextPos = mass.pos + (mass.pos-mass.prevPos) * decr + (totalForce)*mass.invMass*dT*dT;

        if(mass.nextPos[2] < bottomZ)
          mass.nextPos[2] = bottomZ;

      }
    }

    double averShift = 0.0;
    int nMoved = 0;
#pragma omp parallel for reduction (+:averShift) reduction (+:nMoved)
  for(const auto &ijk : masses) {
      int i = ijk[0];
      int j = ijk[1];
      int k = ijk[2];

      auto &mass = getMass(i, j, k);

      if(mass.invMass > 0.0) {
        nMoved++;
        averShift += (mass.nextPos-mass.pos).norm();
        mass.prevPos = mass.pos;
        mass.pos = mass.nextPos;
      }
    }

    if(!(iIter % 100) || iIter == nIter-1)
      printf("%d. %d moved, aver shift %.6f\n", iIter, nMoved, averShift/nMoved);
  }
}

void cSimVolume::extrapolateStrain(int nIter)
{
  for(size_t iM = 0; iM < simVol.size(); iM++) {
    if(simVol[iM].invMass < 0.0)
      simVol[iM].pos = simVol[iM].origPos;
  }

  for(int iIter = 0; iIter < nIter; iIter++) {
    double totalShift = 0;
    int nMoved = 0;

#pragma omp parallel for reduction(+:nMoved) reduction(+:totalShift)   
    for(size_t i = 0; i < Ni; i++)
      for(size_t j = 0; j < Nj; j++)
        for(size_t k = 0; k < Nk; k++) {
  //      printf("%d, %d, %d\n", i, j, k);

          auto &mass = getMass(i, j, k);

          if(mass.invMass < 0.0) {
            int nNeighbors = 0;
            Eigen::Vector3d averPos(0.0, 0.0, 0.0);

            for (int ii = -1; ii <= 1; ii++)
                for (int jj = -1; jj <= 1; jj++)
                  for (int kk = -1; kk <= 1; kk++) {
                    int nZeros = 0;

                    if (ii == 0)
                      nZeros++;
                    if (jj == 0)
                      nZeros++;
                    if (kk == 0)
                      nZeros++;

                    if (nZeros == 2) {

                      if (isValid(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk)) {
                        auto &neighbor = getMass(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk);

                        nNeighbors++;

                        averPos += neighbor.pos;
                      }
                    }
                  }
//                      printf("%d\n", nNeighbors);
            if(nNeighbors == 6) {
              nMoved++;
              averPos *= 1.0/nNeighbors;
              totalShift += (mass.pos-averPos).norm();
              mass.nextPos = averPos;
            }
          }
        }
  
    if(iIter == nIter-1 || !(iIter%(nIter/50)))
      printf("%d. %d moved, aver shift %.6f\n", iIter, nMoved, totalShift/nMoved);

#pragma omp parallel for
    for(size_t iM = 0; iM < simVol.size(); iM++)
      simVol[iM].pos = simVol[iM].nextPos;

  }

  bool extendNotShrink = true;

  for(int iExtend = 0; iExtend < 0; iExtend++) {
    std::vector<std::array<int, 3>> extendVec;

    for(size_t i = 0; i < Ni; i++)
      for(size_t j = 0; j < Nj; j++)
        for(size_t k = 0; k < Nk; k++) {
  //      printf("%d, %d, %d\n", i, j, k);

          auto &mass = getMass(i, j, k);

          if(mass.invMass >= 0.0) {
            for (int ii = -1; ii <= 1; ii++)
                for (int jj = -1; jj <= 1; jj++)
                  for (int kk = -1; kk <= 1; kk++) {

                    if (isValid(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk)) {
                      auto &neighbor = getMass(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk);

                      if(neighbor.invMass < 0.0) {
                        if(extendNotShrink)
                          extendVec.push_back({mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk});
                        else
                          extendVec.push_back({mass.ijk[0], mass.ijk[1], mass.ijk[2]});
                      }
                    }
                  }
          }
        }

    printf("%lu extended\n", extendVec.size());

    for(auto extend : extendVec)
        getMass(extend[0], extend[1], extend[2]).invMass = extendNotShrink ? 2.0 : -2.0;

  }

  std::vector<Eigen::Vector3d> meax;
  std::vector<double> meaxShift;

  double maxShift = 0.0;
  std::array<int, 3> maxIjk;

  for(auto ijk : medianAxis) {
    const auto &mass = getMass(ijk[0], ijk[1], ijk[2]);

    meax.push_back(mass.pos);

    double shift = (mass.pos-mass.origPos).norm()*fabs(mass.sdf);
    meaxShift.push_back(shift);

    if(shift > maxShift) {
      maxShift = std::max(shift, maxShift);
      maxIjk = ijk;
    }
  }

  printf("max dist %.6f\n", maxShift);

  std::vector<Eigen::Vector3d> origPnts;

  for(int ii = -1; ii <=1 ; ii++)
    for(int jj = -1; jj <=1 ; jj++)
      for(int kk = -1; kk <=1 ; kk++)
        if(isValid(maxIjk[0] + ii, maxIjk[1] + jj, maxIjk[2] + kk)) {

          Eigen::Vector3d origPnt = getMass(maxIjk[0] + ii, maxIjk[1] + jj, maxIjk[2] + kk).nearestPos;
          origPnts.push_back(origPnt);
        }

  double maxDist = -1.0;
  std::pair<Eigen::Vector3d, Eigen::Vector3d> maxDistPair;

  for(size_t ii = 0; ii < origPnts.size(); ii++)
    for(size_t jj = 0; jj < origPnts.size(); jj++) {
      double pairDist = (origPnts[ii]-origPnts[jj]).norm();

      if(pairDist > maxDist) {
        maxDist = pairDist;
        maxDistPair = {origPnts[ii], origPnts[jj]};
      }
    }

  printf("%.6f distance\n", maxDist);
  printf("(%.6f, %.6f, %.6f) - (%.6f, %.6f, %.6f)\n", 
          maxDistPair.first[0], maxDistPair.first[1], maxDistPair.first[2],
          maxDistPair.second[0], maxDistPair.second[1], maxDistPair.second[2]);

  happly::PLYData plyOut;
  std::vector<std::array<double, 3>> meshVertexPositions;
  std::vector<std::array<double, 3>> meshVertexColors;

  for(size_t iV = 0; iV < meax.size(); iV++) {
    meshVertexPositions.push_back({meax[iV][0], meax[iV][1], meax[iV][2]});
    double weight = meaxShift[iV]/maxShift;
    meshVertexColors.push_back({weight, 4.0*weight*(1.0-weight), 1.0-weight});
  }

  plyOut.addVertexPositions(meshVertexPositions);
  plyOut.addVertexColors(meshVertexColors);

  // Write the object to file
  plyOut.write("mediaxdist.ply", happly::DataFormat::ASCII);    

  Line3d stiffRib(maxDistPair.first, maxDistPair.second);

  int nStiffRib = 0;
  for(size_t i = 0; i < Ni; i++)
    for(size_t j = 0; j < Nj; j++)
      for(size_t k = 0; k < Nk; k++) {
  //      printf("%d, %d, %d\n", i, j, k);

        auto &mass = getMass(i, j, k);

        if(mass.invMass < 0.0) {

          // distance to the cylinder
          double dist = 3.0*step-sqrt(stiffRib.SegmentDistanceSq(mass.origPos));
          double param = stiffRib.ProjectParam(mass.origPos);
          if(param < 0.0)
            param = 0.0;
          else if (param >= 1.0)
            param = 1.0;

          if(dist >= 0.0) {
            mass.invMass = 1.0;
            mass.sdf = dist;
            mass.nearestPos = stiffRib.Eval(param);

            nStiffRib++;
          }
          else if (dist > mass.sdf) { // distances are negative
            // here it's closer to the cylinder
            mass.sdf = dist;
            mass.nearestPos = stiffRib.Eval(param);
          }
        }
      }

  printf("%d points in stiffness rib\n", nStiffRib);
}

void cSimVolume::computeExternalMedialAxis()
{
    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::array<double, 3>> meshVertexColors;

    for(size_t i = 0; i < Ni; i++)
      for(size_t j = 0; j < Nj; j++)
        for(size_t k = 0; k < Nk; k++) {
  //      printf("%d, %d, %d\n", i, j, k);

          auto &mass = getMass(i, j, k);

          if(mass.invMass < 0.0) {
            int nNeighbors = 0;

            for (int ii = -1; ii <= 1; ii++)
                for (int jj = -1; jj <= 1; jj++)
                  for (int kk = -1; kk <= 1; kk++) {
                    int nZeros = 0;

                    if (ii == 0)
                      nZeros++;
                    if (jj == 0)
                      nZeros++;
                    if (kk == 0)
                      nZeros++;

                    if (nZeros == 2) {

                      if (isValid(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk)) {
                        auto &neighbor = getMass(mass.ijk[0] + ii, mass.ijk[1] + jj, mass.ijk[2] + kk);

                        if(neighbor.invMass < 0.0) {
                           if( (mass.nearestPos-neighbor.nearestPos).norm() > 3.0*step ) {
                              meshVertexPositions.push_back({mass.origPos[0], mass.origPos[1], mass.origPos[2]});
                              meshVertexColors.push_back({1.0, 0, 1.0});
                              medianAxis.push_back(mass.ijk);                          
                           }
                        }
                      }
                    }
                  }
          }
        }

    happly::PLYData plyOut;

    plyOut.addVertexPositions(meshVertexPositions);
    plyOut.addVertexColors(meshVertexColors);

    printf("%d medial axis points\n", medianAxis.size());

    // Write the object to file
    plyOut.write("mediax.ply", happly::DataFormat::ASCII);    
}

// apply loads, compute deformation, find external medial axis points that would be good candidates for placing stiffness ribs through
void SimVolumeTest(MyMesh &mesh, const AABBTree &tree, int nIter, int nCycles, Bounds3d bounds, double step, double dT, double ff)
{
  printf("%d iters, dT = %.6f\n", nIter, dT);
  printf("compute closest points using tree\n");

  happly::PLYData plyOut;

  std::vector<std::array<double, 3>> meshVertexPositions;
  std::vector<std::array<double, 3>> meshVertexColors;

  //  meshVertexPositions.resize(mesh.n_vertices());
  //  meshVertexColors.resize(mesh.n_vertices());

  mesh.request_face_normals();
  mesh.update_normals();

  cSimVolume simVolume;

  //    simVolume.init(bounds, 0.33);
  simVolume.initAnisotropic(bounds, step);

  printf("init, step %.6f\n", simVolume.getStep());

  int nFixedMasses = 0;
  int nForcedMasses = 0;

  for (size_t i = 0; i < simVolume.GetNi(); i++)
    for (size_t j = 0; j < simVolume.GetNj(); j++)
      for (size_t k = 0; k < simVolume.GetNk(); k++)
      {
        auto pnt = simVolume.GetPoint(i, j, k);

        //          if(j == 0 && k == 0)
        //            printf("test point (%6f, %.6f, %.6f)\n", pnt[0], pnt[1], pnt[2]);

        auto closestPnt = tree.FindNearestPoint(pnt);

        //          if(j == 0 && k == 0)
        //            printf("facet %d distance %.6f, (%.6f, %.6f, %.6f)\n", closestPnt.facetIndex, closestPnt.dist,
        //                  closestPnt.pnt[0], closestPnt.pnt[1], closestPnt.pnt[2]);

        //          if(pnt[2] < -3.9)
        //            isFixed = true;

        bool isInside = false;
        if (closestPnt.dist >= 0)
          isInside = true;

        bool isFixed = false;
        bool isForced = false;
        bool isForcedBack = false;

        if (isInside && closestPnt.dist < 1.5 * simVolume.getStep())
        {
          //            printf("(%.6f, %.6f, %.6f)\n", pnt[0], pnt[1], pnt[2]);
          //            printf("(%.6f, %.6f, %.6f), dist %.6f, facet %d\n", pnt[0], pnt[1], pnt[2], closestPnt.dist, closestPnt.facetIndex);

          MyMesh::FaceHandle facetHandle(closestPnt.facetIndex);

          int nFixed = 0;
          int nForced = 0;
          int nForcedBack = 0;

          MyMesh::FaceVertexIter fv_it = mesh.fv_iter(facetHandle);
          for (; fv_it.is_valid(); ++fv_it)
          {
            auto vcolor = mesh.color(*fv_it);

            //              printf("   %d, %d, %d\n", vcolor[0], vcolor[1], vcolor[2]);

            if (vcolor[0] == 255 && vcolor[1] == 0 && vcolor[2] == 0)
              nForced++;

            if (vcolor[2] == 255 && vcolor[1] == 0 && vcolor[0] == 0)
              nFixed++;

            if (vcolor[0] == 255 && vcolor[1] == 255 && vcolor[2] == 0)
              nForcedBack++;
          }

          if (nForced == 3)
            isForced = true;

          if (nForcedBack == 3)
          {
            isForcedBack = true;
            printf("*");
          }

          if (nFixed == 3)
            isFixed = true;

          //            if(isForced || isFixed)
          //              printf("   %d, %d\n", isForced, isFixed);
        }

        Eigen::Vector3d gForce(0.0, 0.0, 0.0);

        if (isForced)
          gForce[2] = ff;
        else if (isForcedBack)
          gForce[2] = -ff;

        if (isFixed)
          nFixedMasses++;

        if (isForced || isForcedBack)
          nForcedMasses++;

        double invMass = isInside ? isFixed ? 0.0 : 1.0 : -1.0;

        simVolume.SetProperties(i, j, k, gForce, invMass);

        //          if(invMass < 0.0)
        simVolume.SetExternalPoint(i, j, k, closestPnt.pnt, closestPnt.dist);
      }

  // Add mesh data (elements are created automatically)
  // plyOut.addVertexPositions(meshVertexPositions);
  // plyOut.addVertexColors(meshVertexColors);

  // Write the object to file
  // plyOut.write("pnts.ply", happly::DataFormat::ASCII);

  printf("%d fixed, %d forced masses\n", nFixedMasses, nForcedMasses);

  for (int iTry = 0; iTry < nCycles; iTry++)
  {

    simVolume.computeExternalMedialAxis();

    simVolume.savePLY("init.ply", 0);

    simVolume.simulate(nIter, dT);

    simVolume.extrapolateStrain(1000);

    simVolume.savePLY("simul.ply", 1);
  }

  //    simVolume.savePLY("simul1.ply", 1);

  int nGrid = 10;
  for (int i = 0; i <= nGrid; i++)
  {
    double x = (bounds.minBound[0] * (nGrid - i) + bounds.maxBound[0] * i) * 1.0 / nGrid;

    for (int j = 0; j <= nGrid; j++)
    {
      double y = (bounds.minBound[1] * (nGrid - j) + bounds.maxBound[1] * j) * 1.0 / nGrid;

      for (int k = 0; k <= nGrid; k++)
      {
        double z = (bounds.minBound[2] * (nGrid - k) + bounds.maxBound[2] * k) * 1.0 / nGrid;

        auto closestPnt = tree.FindNearestPoint({x, y, z});

        auto normal = mesh.normal(MyMesh::FaceHandle(closestPnt.facetIndex));
        //            printf("%.6f, %.6f, %.6f\n", normal[0], normal[1], normal[2]);

        meshVertexPositions.push_back({x, y, z});

        if (closestPnt.dist < 0.0)
          meshVertexColors.push_back({1.0, 0, 0});
        else
          meshVertexColors.push_back({0, 0, 1.0});

        meshVertexPositions.push_back({closestPnt.pnt[0], closestPnt.pnt[1], closestPnt.pnt[2]});
        meshVertexColors.push_back({0, 1.0, 0});
      }
    }
  }

  // Add mesh data (elements are created automatically)
  plyOut.addVertexPositions(meshVertexPositions);
  plyOut.addVertexColors(meshVertexColors);

  // Write the object to file
  plyOut.write("pnts.ply", happly::DataFormat::ASCII);
};

//  __SimVolumeTest(nIter, nCycles, bounds, step, dT, ff);

#pragma once

#include <stdint.h>
#include <fstream>
#include <vector>
#include <array>
#include <math.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

// mass and spring simulation class

struct Bounds3d;

struct sMassPoint
{
  std::array<int, 3> ijk;
  Eigen::Vector3d origPos;
  Eigen::Vector3d pos;
  Eigen::Vector3d nextPos;
  Eigen::Vector3d prevPos;
  double invMass; // zero mass - fixed, negative mass - mass is zero
  Eigen::Vector3d force;

  Eigen::Vector3d nearestPos;
  double sdf;

//  double strain;
//  bool mustBe;
};

class cSimVolume {

public:

  void init(const Bounds3d &, double extendFactor = 0.1);
  void initAnisotropic(const Bounds3d &bounds, double step);

  size_t GetN() const { return N; }
  size_t GetNi() const { return Ni; }
  size_t GetNj() const { return Nj; }
  size_t GetNk() const { return Nk; }

  Eigen::Vector3d GetPoint(int i, int j, int k) const { return getMass(i, j, k).pos; }
  void SetProperties(int i, int j, int k, Eigen::Vector3d force, double inverseMass) {
    auto &mass = getMass(i, j, k);

    mass.force = force;
    mass.invMass = inverseMass;
    // invMass < 0.0 - empty
    // invMass == 0.0 - fixed
    // invMass > 0.0 - regular simulation
  }

  void SetExternalPoint(int i, int j, int k, Eigen::Vector3d nearestPnt, double sdf) {
    auto &mass = getMass(i, j, k);

    mass.nearestPos = nearestPnt;
    mass.sdf = sdf;
    // invMass < 0.0 - empty
    // invMass == 0.0 - fixed
    // invMass > 0.0 - regular simulation
  }

  void simulate(int nIter, double dT);
  void extrapolateStrain(int nIter);
  void computeExternalMedialAxis();

//  void computeStrain();
  void savePLY(std::string filename, int type = 0);

  size_t getIndex(int i, int j, int k) const { return i + N*j + N*N*k;}
  double getStep() const { return step;}

private:

  const size_t N = 150;
  const double eps = 1e-10;

  bool isValid(int i, int j, int k) const { return i >= 0 && (size_t)i < Ni && j >= 0 && (size_t)j < Nj && k >= 0 && (size_t)k < Nk;}
  size_t index(int i, int j, int k) const { return i + Ni*j + Ni*Nj*k; }

  sMassPoint &getMass(int i, int j, int k) {
    size_t l = index(i, j, k);

    if(l >= simVol.size())
      printf("!!!! %lu : %d, %d, %d!!!\n", l, i, j, k);

    return simVol[l];
  }

  const sMassPoint &getMass(int i, int j, int k) const {
    return simVol[index(i, j, k)];
  }

  Eigen::Vector3d getForce(const sMassPoint &mass) const;

  Eigen::MatrixXd getStrainTensor(const sMassPoint &mass) const;
  Eigen::MatrixXd getStressTensor(const sMassPoint &mass) const;

  double getStrain(const sMassPoint &mass) const;
  std::vector<double> getStrainEigenValues(const sMassPoint &mass) const;

private:

  size_t Ni = N;
  size_t Nj = N;
  size_t Nk = N;
  
  double step;
  std::vector<sMassPoint> simVol;
  std::vector<std::array<int, 3>> medianAxis;
};

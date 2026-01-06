#include "tools3d.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include "happly.h"

Eigen::Vector3d cross(const Eigen::Vector3d &vec1, const Eigen::Vector3d &vec2)
 {
    //Eigen::Vector3d crossProd = {0.0, 0.0, 0.0};
    //auto cp = vec1.cross(vec2); for some reason, this is undefined for Vector3d, but defined for Vector3f

    return {vec1[1]*vec2[2] - vec1[2]*vec2[1],
            vec1[2]*vec2[0] - vec1[0]*vec2[2], 
            vec1[0]*vec2[1] - vec1[1]*vec2[0]};
}

// given dirVec and 4 vectors or rosy(cross) field 
// defined by fieldVec, fieldVec1 = cross(normal, fieldVec), fieldVec2 = cross(normal, fieldVec1), fieldVec3 = cross(normal, fieldVec2)
// return the cross field vector maximizing the dot product with dirVec
Eigen::Vector3d bestRosyFit(const Eigen::Vector3d &dirVec, const Eigen::Vector3d &normal, const Eigen::Vector3d &fieldVec)
{
    Eigen::Vector3d nextFieldVec = fieldVec;

    double bestDotProd = -1.0;
    Eigen::Vector3d bestFieldVec;

    for(int ii = 0; ii < 4; ii++) {

        double dotProd = dirVec.dot(nextFieldVec);

        if(dotProd > bestDotProd) {
            bestDotProd = dotProd;
            bestFieldVec = nextFieldVec;
        }

        nextFieldVec = nextFieldVec.cross(normal);
    }

    return bestFieldVec;
}

// sample segment by points equally spaced at max distance <= step
std::vector<Eigen::Vector3d> Sample(const Line3d &seg, double step)
{
  double segLen = seg.Length();

  std::vector<Eigen::Vector3d> pnts;

  if(segLen > step) {

    int nInternalPnts = static_cast<int>(round(segLen/step));

    for(int iPnt = 1; iPnt <= nInternalPnts; iPnt++) {
      double param = 1.0*iPnt/(nInternalPnts+1);

      auto pnt = seg.Eval(param);

      pnts.push_back(pnt);
    }
  }

  return pnts;
}

// TODO : Sample(const Triangle3d &triangle, double step);

bool DoIntersect(const Triangle3d &triA, const Triangle3d &triB)
{
  double V0[3] = { triA.Pnt(0)[0], triA.Pnt(0)[1], triA.Pnt(0)[2]};

  double V1[3] = { triA.Pnt(1)[0], triA.Pnt(1)[1], triA.Pnt(1)[2]};

  double V2[3] = { triA.Pnt(2)[0], triA.Pnt(2)[1], triA.Pnt(2)[2]};

  double U0[3] = { triB.Pnt(0)[0], triB.Pnt(0)[1], triB.Pnt(0)[2]};

  double U1[3] = { triB.Pnt(1)[0], triB.Pnt(1)[1], triB.Pnt(1)[2]};

  double U2[3] = { triB.Pnt(2)[0], triB.Pnt(2)[1], triB.Pnt(2)[2]};

  int doIntersect = tri_tri_intersect(V0, V1, V2, U0, U1, U2);

  return doIntersect;
}


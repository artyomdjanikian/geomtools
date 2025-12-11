#pragma once

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

static double eps = 1e-10;

//template <typename T> Eigen::Vector3d& operator=(Eigen::Vector3d &v, const T& t)
//{
//    v = Eigen::Vector3d(t[0], t[1], t[2]);
//    return v;
//}
// this isn't working, compiler complains that ‘Eigen::Vector3d& operator=(Eigen::Vector3d&, const T&)’ must be a non-static member function
// TODO : try w/o template

template <typename T> Eigen::Vector3d toVec(const T& t)
{
    return Eigen::Vector3d(t[0], t[1], t[2]);
}

namespace {
    double distSq(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
    {
        double distanceSq = 0.0;

        for (size_t iC = 0; iC < 3; iC++)
        {
            double thisDist = a[iC] - b[iC];
            distanceSq += thisDist * thisDist;
        }

        return distanceSq;
    }
}

// from tribox3.cc
int triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3]);

// from tritri_isectline.cc
int tri_tri_intersect_with_isectline(double V0[3],double V1[3],double V2[3],
				                    double U0[3],double U1[3],double U2[3],int *coplanar,
				                    double isectpt1[3],double isectpt2[3]);

// from tools3d.cc
Eigen::Vector3d cross(const Eigen::Vector3d &vec1, const Eigen::Vector3d &vec2);

Eigen::Vector3d bestRosyFit(const Eigen::Vector3d &dirVec, const Eigen::Vector3d &normal, const Eigen::Vector3d &nextDirVec );

struct Line3d {

    Line3d(Eigen::Vector3d pntA, Eigen::Vector3d pntB) { pnts[0] = pntA, pnts[1] = pntB;}

    bool IsDegen() const { return (pnts[1]-pnts[0]).norm() < eps;}

    Eigen::Vector3d Dir() const { Eigen::Vector3d dir = pnts[1] - pnts[0]; dir.stableNormalize(); return dir; }

    double Length() const { return(pnts[1] - pnts[0]).norm(); }

    Eigen::Vector3d Eval(double param) const { return pnts[0] * (1.0-param) + pnts[1]*param;}

    double ProjectParam( Eigen::Vector3d pnt ) const { auto span = pnts[1] - pnts[0];
                                                        double dotProd = (pnt - pnts[0]).dot(span);
                                                        double segLenSq = span.dot(span);
                                                        return dotProd/segLenSq;
                                                        }

    Eigen::Vector3d Project(Eigen::Vector3d pnt) const { auto segVec = pnts[1] - pnts[0];
                                                         double dotProd = (pnt - pnts[0]).dot(segVec);
                                                         double segLenSq = segVec.dot(segVec);
                                                         double factor = dotProd * eps > segLenSq ? dotProd/segLenSq : 0.5;
                                                         return pnts[0] + segVec*factor;
                                                       }

    double SegmentDistanceSq(Eigen::Vector3d pnt) const {
        auto span = pnts[1]-pnts[0];
        double dotProd = (pnt-pnts[0]).dot(span);
        double segLenSq = span.dot(span);

        if(segLenSq > eps*eps)
            dotProd /= segLenSq;  

        double pdistSq = -1.0;

        if(dotProd <= eps)
            pdistSq = distSq(pnts[0], pnt);
        else if(dotProd >= 1.0-eps)
            pdistSq = distSq(pnts[1], pnt);
        else
            pdistSq = distSq(Eval(dotProd), pnt);

        return pdistSq;
    }

    bool Skew(const Line3d &other, std::pair<double, double> *pParams) const {

        auto dir0 = pnts[1]-pnts[0];
        auto dir1 = other.pnts[1]-other.pnts[0];
        auto R = other.pnts[0]-pnts[0];

        auto crossDir = cross(dir0, dir1);
        double crossNorm = crossDir.norm();

        // parallel or coincident
        if(crossNorm < eps)
            return false;

        double t0 = cross(dir1, crossDir).dot(R)/(crossDir.dot(crossDir));
        double t1 = cross(dir0, crossDir).dot(R)/(crossDir.dot(crossDir));

        if(pParams != nullptr) {
            pParams->first = t0;
            pParams->second = t1;
        }

        return true;
    }

    Eigen::Vector3d pnts[2];
};

struct Plane3d {

    Plane3d(Eigen::Vector3d pntA, Eigen::Vector3d pntB, Eigen::Vector3d pntC) {
        double a1 = pntB[0] - pntA[0]; //x2 - x1;
        double b1 = pntB[1] - pntA[1]; //y2 - y1;
        double c1 = pntB[2] - pntA[2]; //z2 - z1;
        double a2 = pntC[0] - pntA[0]; //x3 - x1;
        double b2 = pntC[1] - pntA[1]; //y3 - y1;
        double c2 = pntC[2] - pntA[2]; //z3 - z1;

        a = b1 * c2 - b2 * c1;
        b = a2 * c1 - a1 * c2;
        c = a1 * b2 - b1 * a2;

        double totalDist = a*a + b*b + c*c;
        totalDist = sqrt(totalDist);
        if(totalDist > eps) {
            a/=totalDist;
            b/=totalDist;
            c/=totalDist;
        }

        d = (- a * pntA[0] - b * pntA[1] - c * pntA[2]);
    }

    Plane3d(Eigen::Vector3d pnt, Eigen::Vector3d normal) { a = normal[0], b = normal[1], c = normal[2], d = -normal.dot(pnt); }

    void Flip() { a *= -1.0; b *= -1.0; c *= -1.0; d *= -1.0; }
    bool IsDegen() const { return fabs(a) < eps && fabs(b) < eps && fabs(c) < eps; }

    Eigen::Vector3d GetNormal() const { return {a, b, c};}

    double SignedDistance(Eigen::Vector3d pnt) const { return pnt[0]*a + pnt[1]*b + pnt[2]*c + d; }

    Eigen::Vector3d Project(Eigen::Vector3d pnt) const { return pnt - GetNormal()*SignedDistance(pnt);}

    enum eXType {
        None=10,
        Point,
        Line
    };

    eXType Intersect(Line3d line, Eigen::Vector3d &xPnt) const {
        eXType xType = None;

        double distA = SignedDistance(line.pnts[0]);
        double distB = SignedDistance(line.pnts[1]);

        double distDiff = distA-distB;
        
        if( fabs(distDiff) > eps ) {
            xType = Point;
            xPnt = (line.pnts[1]*distA - line.pnts[0]*distB)/distDiff;
        }
        else if(fabs(distDiff) < eps)
            xType = Line;

        return xType;
    }

private:
    double a, b, c, d;
};

struct Bounds3d {
    Bounds3d & operator += (const Eigen::Vector3d &pnt) {

        if(!isInit) {
            minBound = pnt;
            maxBound = pnt;
            isInit = true;
        }
        else {
            for(size_t iC = 0; iC < 3; iC++) {
                minBound[iC] = std::min(minBound[iC], pnt[iC]);
                maxBound[iC] = std::max(maxBound[iC], pnt[iC]);
            }
        }

        return *this;
    }

    bool IsPointInside(const Eigen::Vector3d &pnt) const {
        for(size_t iC = 0; iC < 3; iC++)
            if(pnt[iC] < minBound[iC]-eps || pnt[iC] > maxBound[iC]+eps)
                return false;

        return true;
    }

    double PointDistanceSq(Eigen::Vector3d pnt) const {
        Eigen::Vector3d dists(0.0, 0.0, 0.0);

        for(size_t iC = 0; iC < 3; iC++)
            if(pnt[iC] > maxBound[iC])
                dists[iC] = pnt[iC] - maxBound[iC];
            else if (pnt[iC] < minBound[iC])
                dists[iC] = minBound[iC] - pnt[iC];


        return dists.dot(dists);
    }

    bool IntersectBounds(const Bounds3d &other) const {

        std::array<bool, 3> overlap = {false, false, false};

        for(size_t iC = 0; iC < 3; iC++) 
            if( minBound[iC] < other.maxBound[iC]+eps ||
                other.minBound[iC] < maxBound[iC]+eps)
                overlap[iC] = true;

        return overlap[0] && overlap[1] && overlap[2];
    }

    // TODO : IntersectTriangle

    bool IntersectLine(const Line3d &line, std::pair<double, double> * pParams = nullptr) const {

        std::array<bool, 6> hasInt = {false, false, false, false, false, false};
        double intParams[6];
        int nParams = 0;

//        std::vector<Eigen::Vector3d> xpnts;
//        xpnts.reserve(6);

        for(int iFace = 0; iFace < 6; iFace++) {

            Eigen::Vector3d normal;

            if(iFace == 0 || iFace == 3)
                normal = {0.0, 0.0, 1.0};
            else if (iFace == 1 || iFace == 4)
                normal = {0.0, 1.0, 0.0};
            else if (iFace == 2 || iFace == 5)
                normal = {1.0, 0.0, 0.0};

            Eigen::Vector3d pnt;

            if(iFace < 3)
                pnt = minBound;
            else
                pnt = maxBound;

            Plane3d plane(pnt, normal);

            double d0 = plane.SignedDistance(line.pnts[0]);
            double d1 = plane.SignedDistance(line.pnts[1]);

            if(fabs(d1-d0) > eps) {
                double tZeroDist = d0/(d0-d1);

                auto xPnt = line.Eval(tZeroDist);
//                xpnts.push_back(xPnt);

                if(IsPointInside(xPnt)) {
                    intParams[nParams++] = tZeroDist;
                    hasInt[iFace] = true;
                }
            }
        }

        if(nParams == 1) {
            printf("!!! wrong number of intersections !!!\n");

            if(pParams != nullptr)
                *pParams = {intParams[0], intParams[1]};

            return true;            
        }

        if(nParams >= 2) {
            double minT = intParams[0];
            double maxT = intParams[0];

            for(int iParam = 0; iParam < nParams; iParam++) {
                minT = std::min(minT, intParams[iParam]);
                maxT = std::max(maxT, intParams[iParam]);
            }

            if(pParams != nullptr)
                *pParams = {minT, maxT};

            return true;
        }

        return false;
    }

    bool isInit = false;
    Eigen::Vector3d minBound;
    Eigen::Vector3d maxBound;
};

struct Triangle3d {

    Triangle3d(Eigen::Vector3d pntA, Eigen::Vector3d pntB, Eigen::Vector3d pntC) { pnts[0] = pntA, pnts[1] = pntB, pnts[2] = pntC; }

    enum eTopoElem { Facet = 10, Edge01, Edge12, Edge20, Vertex0, Vertex1, Vertex2, Undef};

    const Eigen::Vector3d &Pnt(int i) const { return pnts[i];}

    Eigen::Vector3d GetAreaNormal() const {
        Eigen::Vector3d normal(0.0, 0.0, 0.0);

        for(size_t iV = 0; iV < 3; iV++) {
            size_t iNextV = (iV+1)%3;
            size_t iOppV = (iNextV +1)%3;

            for(size_t iC = 0; iC < 3; iC++) {
                size_t iC1 = (iC+1)%3;
                size_t iC2 = (iC1+1)%3;

                normal[iC] += pnts[iNextV][iC2]*pnts[iV][iC1] - pnts[iNextV][iC1]*pnts[iV][iC2];
            }
        }

        return normal*0.5;
    }

    Eigen::Vector3d GetNormal() const {
        auto normal = GetAreaNormal();
        normal.stableNormalize();

        return normal;
    }

    double GetArea() const {
        auto normal = GetAreaNormal();

        return normal.norm();
    }

    Plane3d GetPlane() const {
        return Plane3d(pnts[0], pnts[1], pnts[2]);
    }

    bool IsPointInside(Eigen::Vector3d pnt) const {

        auto normal = GetNormal();

        for(size_t iV = 0; iV < 3; iV++) {
            size_t iNextV = (iV+1)%3;
            size_t iOppV = (iNextV +1)%3;

            Eigen::Vector3d edgeVec = pnts[iNextV]-pnts[iV];
            edgeVec.stableNormalize();

            auto edgeNormal = cross(edgeVec, normal);
            double signDistPnt = edgeNormal.dot(pnt-pnts[iV]);
            double signDistOpp = edgeNormal.dot(pnts[iOppV]-pnts[iV]);

            if(signDistPnt * signDistOpp <= 0.0)
                return false;
        }

        return true;
    }

    // doesn't work if triangle or line are degen
 
    std::pair<Eigen::Vector3d, double> PointDistanceSq(Eigen::Vector3d pnt) const {

        auto plane = Plane3d(pnts[0], pnts[1], pnts[2]);

        double facetDist = -1.0;

        if(!plane.IsDegen()) {
            // project onto the plane

            auto projPnt = plane.Project(pnt);

            // check if projPnt is inside triangle

            bool isInside = IsPointInside(projPnt);

            if(isInside)
                return {projPnt, distSq(pnt, projPnt)};
        }

        double edgeDists[3];
        Eigen::Vector3d edgePnts[3];
        std::array<int, 3> vertexIds = {-1, -1, -1};

        for(size_t iEdge = 0; iEdge < 3; iEdge++) {

            size_t iNextEdge = (iEdge+1)%3;

            Line3d lineSeg(pnts[iEdge], pnts[iNextEdge]);

            if(!lineSeg.IsDegen()) {
                double param = lineSeg.ProjectParam(pnt);
                double segLen = lineSeg.Length();

//                printf("    edge %d(%.6f, %.6f, %.6f) - (%.6f, %.6f, %.6f), param %.6f\n",
//                         iEdge, pnts[iEdge][0], pnts[iEdge][1], pnts[iEdge][2], pnts[iNextEdge][0], pnts[iNextEdge][1], pnts[iNextEdge][2], param);

                Eigen::Vector3d segPnt(0.0, 0.0, 0.0);

                // snap to vertices
                if(param*segLen <= eps) {
                    segPnt = pnts[iEdge];
                    vertexIds[iEdge] = iEdge;
                }
                else if (param*segLen >= segLen-eps) {
                   segPnt = pnts[iNextEdge];
                    vertexIds[iEdge] = iNextEdge;
                }
                else // (param*segLen > eps && param*segLen < segLen-eps) - segment interior
                    segPnt = lineSeg.Eval(param);

                edgeDists[iEdge] = distSq(pnt, segPnt);
                edgePnts[iEdge] = segPnt;
            }
            else {
                auto segPnt = pnts[iEdge];
                vertexIds[iEdge] = iEdge;
                edgeDists[iEdge] = distSq(pnt, segPnt);
                edgePnts[iEdge] = segPnt;
            }
        }

        size_t iShortestDist = 0;

        for(size_t iC = 0; iC < 3; iC++)
            if(edgeDists[iC] < edgeDists[iShortestDist])
                iShortestDist = iC;

        return {edgePnts[iShortestDist], edgeDists[iShortestDist]};
    }

   int IntersectLine(const Line3d &line, std::pair<double, double> *pParams = nullptr) const {
        int nX = 0;

        auto plane = Plane3d(pnts[0], pnts[1], pnts[2]);

        double d0 = plane.SignedDistance(line.pnts[0]);
        double d1 = plane.SignedDistance(line.pnts[1]);

        std::vector<double> xSegParams;

        // line in plane
        if (fabs(d0) <= eps && fabs(d1) <= eps) {

            for(size_t iEdge = 0; iEdge < 3; iEdge++) {

                size_t iNextEdge = (iEdge+1)%3;

                if(iNextEdge >= 3)
                    iNextEdge = 0;

                Line3d edgeSeg(pnts[iEdge], pnts[iNextEdge]);

                std::pair<double, double> params;
                auto areSkew = line.Skew(edgeSeg, &params);

                if(areSkew) {
                    // TODO : snap to segment end points
                    if(params.second >= 0 && params.second <= 1.0) {
                        xSegParams.push_back(params.first);
                    }
                }
                else {
                     // check for coincident segments
                    auto projectPnt0 = line.Project(pnts[iEdge]);
                    auto projectPnt1 = line.Project(pnts[iNextEdge]);

                    if( (projectPnt0 - pnts[iEdge]).norm() < eps &&
                        (projectPnt1 - pnts[iNextEdge]).norm() < eps) {
                            xSegParams.push_back(line.ProjectParam(projectPnt0));
                            xSegParams.push_back(line.ProjectParam(projectPnt1));
                        }
                }
            }
        }
        else if (fabs(d0) <= eps) {
            // snap to end
            nX = 1;
            if(pParams != nullptr)
                pParams->first = 0.0;
        }
        else if (fabs(d1) <= eps) {
            // snap to end
            nX = 1;
            if(pParams != nullptr)
                pParams->first = 1.0;
        }
        // line intersects plane
        if(fabs(d0-d1) > eps) {
            double tZeroDist = d0/(d0-d1);

            auto xPnt = line.Eval(tZeroDist);

//            printf("      (%.6f, %.6f, %.6f)->%.6f\n", xPnt[0], xPnt[1], xPnt[2], plane.SignedDistance(xPnt));

            auto pntDist = PointDistanceSq(xPnt);

//            printf("         triangle dist %.6f (%.6f, %.6f, %.6f)\n", sqrt(pntDist.second), pntDist.first[0], pntDist.first[2], pntDist.first[2]);
            if(sqrt(pntDist.second) < eps) {
                nX++;

                if(pParams != nullptr)
                    pParams->first = tZeroDist;
            }
        }

        // TODO : fix nX and xSegParams
        if(nX + xSegParams.size() > 1) {
            if(nX)
                xSegParams.push_back(pParams->first);

            std::sort(xSegParams.begin(), xSegParams.end());

            if(pParams != nullptr) {
                pParams->first = xSegParams.front();
                pParams->second = xSegParams.back();

                nX = 2;
            }
        }

        return nX;
    }

    Eigen::Vector3d Barycentric(Eigen::Vector3d P) const
    {
        Eigen::Vector3d bary;

        auto normal = GetAreaNormal();

        const auto &a = pnts[0];
        const auto &b = pnts[1];
        const auto &c = pnts[2];

        double areaABC = normal.dot( cross( (b - a), (c - a) )  ) ;
        double areaPBC = normal.dot( cross( (b - P), (c - P) )  ) ;
        double areaPCA = normal.dot( cross( (c - P), (a - P) )  ) ;

        bary[0] = areaPBC / areaABC ; // alpha
        bary[1] = areaPCA / areaABC ; // beta
        bary[2] = 1.0f - bary[0] - bary[1] ; // gamma

        return bary;
    }

    Eigen::Vector3d EvalBarycentric(Eigen::Vector3d uvw) const { return pnts[0]*uvw[0] + pnts[1]*uvw[1] + pnts[2]*uvw[2]; }

    double AspectRatio() const {
        double a = (pnts[0]-pnts[1]).norm();
        double b = (pnts[1]-pnts[2]).norm(); 
        double c = (pnts[2]-pnts[0]).norm(); 

        double s = (a + b + c) / 2.0;
        double area = std::sqrt(s * (s - a) * (s - b) * (s - c));
        double circumradius = (a * b * c) / (4 * area);
        double inradius = area / s;
        return circumradius > eps ? (2 * inradius)/circumradius : 0.0;
    }

//    std::vector<Eigen::Vector3d> Sample(double step) const;
// TODO : implement

private:
    std::array<Eigen::Vector3d, 3> pnts;
};

std::vector<Eigen::Vector3d> Sample(const Line3d &seg, double step);


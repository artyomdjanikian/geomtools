#include <vector>
#include <list>
#include <tuple>
#include <eigen3/Eigen/Core>

// triangulate 2d polygon with otectomy (ear removal) algorithm

#pragma once

namespace polygon2d {

    // for numerical stability,
    // segments shorter than/areas smaller than
    // eps = 1e-10 are considered degenerate
    const double eps = 1e-10;

    // 2d line segment struct
    struct Line2d
    {
        Line2d(Eigen::Vector2d pntA, Eigen::Vector2d pntB) { pnts[0] = pntA, pnts[1] = pntB; }

        bool IsDegen() const { return (pnts[1] - pnts[0]).norm() < eps; }

        Eigen::Vector2d Eval(double param) const { return pnts[0] * (1.0 - param) + pnts[1] * param; }

        double SegmentDistance(Eigen::Vector2d pnt) const
        {
            auto span = pnts[1] - pnts[0];
            double dotProd = (pnt - pnts[0]).dot(span);
            double segLenSq = span.dot(span);

            dotProd /= segLenSq;

            double pointDist = -1.0;

            if (dotProd < 0.0)
                pointDist = (pnts[0] - pnt).norm();
            else if (dotProd > 1.0)
                pointDist = (pnts[1] - pnt).norm();
            else
                pointDist = (Eval(dotProd) - pnt).norm();

            return pointDist;
        }

        Eigen::Vector2d pnts[2];
    };

    class ConvexifyPolygon {

    public:

        // create input polygon
        void push_back(Eigen::Vector2d point) {
            points_.push_back(point);
        }

        // access polygon points
        const Eigen::Vector2d & operator[](size_t index) const {
            return points_[index];
        }

        // count polygon points
        size_t size() const {
            return points_.size();
        }
        
        // main convex partitioning function
        double computeConvexPartition() {

            // if loop has duplicate points at the ends, remove one
            if (points().size() && (points().front() - points().back()).norm() < eps)
                points().pop_back();

            // sanity check
            if(points().size() < 3)
                return 0.0;

            double totalArea = computeSignedArea();

            std::list<size_t> pointChain;

            for(size_t i = 0; i < points().size(); i++) {
                pointChain.push_back(i);
            }

            double area = totalArea;

            while(pointChain.size() > 0) {

                struct sEar
                {
                    size_t ii = std::numeric_limits<size_t>::max();
                    size_t jj = std::numeric_limits<size_t>::max();
                    size_t kk = std::numeric_limits<size_t>::max();
                    double area = 0.0;

                    bool isEmpty() const { return ii == SIZE_MAX; }
                };

                sEar bestEar;

                // loop over vertices, find convex triangle
                for(auto iter = pointChain.begin(); iter != pointChain.end(); ++iter) {

                    auto accessIter = iter;
                    size_t ii = *accessIter++;

                    if(accessIter == pointChain.end())
                        accessIter = pointChain.begin();

                    size_t jj = *accessIter++;

                    if(accessIter == pointChain.end())
                        accessIter = pointChain.begin();

                    size_t kk = *accessIter;

                    // compute area of this triangle
                    double triangleArea = triangleSignedArea(points()[ii], points()[jj], points()[kk]);

                    // degenerate chain, pick only degenerate triangle
                    if (fabs(area) <= eps) {
                        if (fabs(triangleArea) <= eps) {
                            bestEar = {ii, jj, kk, triangleArea};
                            break;
                        }
                    }
                    // regular chain, pick the smallest area triangle with orientation
                    // (clockwise or counterclockwise) matching the orientation of chain
                    // this is an ear
                    else if (fabs(triangleArea) > eps && triangleArea * area > 0.0 && fabs(area) >= fabs(triangleArea)-eps) {
                        // check if diagonal intersects polygon edges
                        Line2d diag(points()[ii], points()[kk]);

                        bool xFound = false;
                        for(auto segIter = pointChain.begin(); segIter != pointChain.end(); ++segIter) {
                        
                            auto accessIter = segIter;
                            size_t iiX = *accessIter++;

                            if(accessIter == pointChain.end())
                                accessIter = pointChain.begin();

                            size_t jjX = *accessIter++;

                            // those segments have an end point coinciding with a diagonal end points
                            // we don't need to check them
                            if(iiX == ii || jjX == ii || iiX == kk || jjX == kk)
                                continue;

                            xFound = doSegmentsIntersect(diag, Line2d(points()[iiX], points()[jjX]));

                            if(xFound)
                                break;
                        }

                        if (!xFound /*&& (bestEar.isEmpty() || fabs(bestEar.area) > fabs(triangleArea))*/ ) {
                            bestEar = {ii, jj, kk, triangleArea};
                            break;
                        }
                    }
                }

                // no ear found
                // this can only happen if the chain is self-intersecting
                // in this case, no valid convexification exists and area is undefined
                if (bestEar.isEmpty()) {
                    triangles_.clear();
                    return -1.0;
                }
                else {
                    // add ear to convex partition

                    triangles_.emplace_back(bestEar.ii, bestEar.jj, bestEar.kk);
                    area -= bestEar.area;

                    // last triangle added
                    if(pointChain.size() == 3)
                        pointChain.clear();
                    // remove "ear tip" vertex from the chain
                    else {
                        for (auto iter = pointChain.begin(); iter != pointChain.end(); ++iter) {
                            if(*iter == bestEar.jj) {
                                pointChain.erase(iter);
                                break;
                            }
                        }
                    }
                }
            }

            return fabs(totalArea);
        }

        // a triangle is represented by 3 vertex indices
        struct sTriangle
        {
            sTriangle(size_t _index0, size_t _index1, size_t _index2) : index0(_index0), index1(_index1), index2(_index2) {}
            
            size_t index0 = std::numeric_limits<size_t>::max();
            size_t index1 = std::numeric_limits<size_t>::max();
            size_t index2 = std::numeric_limits<size_t>::max();
        };

        // count output triangles
        size_t numTriangles() const {
            return triangles_.size();
        }

        // access triangles
        const sTriangle& getTriangle(size_t index) const {
            return triangles_[index];
        }

        void saveOFF(std::string filename) const {
            std::ofstream file(filename);

            file << "OFF" << std::endl;

            file << points().size() << " " << triangles().size() << std::endl;

            for (const auto &point : points())
                file << point[0] << " " << point[1] << " 0\n";

            for(const auto &tr : triangles())
                file << "3 " << tr.index0 << " " << tr.index1 << " " << tr.index2 << std::endl;

            file.close();
        }

    private:

        // signed area of the polygon, works for concave polygons
        double computeSignedArea() const {
            double area = 0.0;

            for(size_t i = 2; i < points().size(); i++) {
                size_t j = i-1;
                auto thisArea = triangleSignedArea(points()[0], points()[j], points()[i]);

                area += thisArea;
            }
            return area;
        }

        // signed area of a single triangle
        double triangleSignedArea(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1, const Eigen::Vector2d &p2) const {

            double x0 = p0[0], y0 = p0[1];
            double x1 = p1[0], y1 = p1[1];
            double x2 = p2[0], y2 = p2[1];

            double area = 0.5*( x0*(y1-y2) + x1*(y2-y0) + x2*(y0-y1) );

            return area;
        }

        // intersect 2 2d segments
        bool doSegmentsIntersect(const Line2d &seg, const Line2d &xSeg) const {

            // degenerate cases, for generality
            if(seg.IsDegen() && xSeg.IsDegen()) {
                return (seg.Eval(0.5) - xSeg.Eval(0.5)).norm() < eps;
            }
            else if (xSeg.IsDegen()) {
                return seg.SegmentDistance(xSeg.Eval(0.5)) < eps;
            }
            else if (seg.IsDegen()) {
                return xSeg.SegmentDistance(seg.Eval(0.5)) < eps;
            }

            // check for line intersection

//            double lenSeg = seg.Length();
//            double lenXSeg = xSeg.Length();

            double x00 = seg.pnts[0][0], y00 = seg.pnts[0][1];
            double x01 = seg.pnts[1][0], y01 = seg.pnts[1][1];
            double dx0 = x01 - x00, dy0 = y01 - y00;

//            double len0Sq = dx0 * dx0 + dy0 * dy0;

            double x10 = xSeg.pnts[0][0], y10 = xSeg.pnts[0][1];
            double x11 = xSeg.pnts[1][0], y11 = xSeg.pnts[1][1];
            double dx1 = x11 - x10, dy1 = y11 - y10;

//            double len1Sq = dx1 * dx1 + dy1 * dy1;

            double denom = dx0*dy1 - dy0*dx1;

            // segments are parallel
            if(fabs(denom) < eps) {

                // check for overlapping segments
                double dist00 = seg.SegmentDistance(xSeg.Eval(0.0));
                if(dist00 < eps)
                    return true;

                double dist01 = seg.SegmentDistance(xSeg.Eval(1.0));
                if(dist01 < eps)
                    return true;

                double dist10 = xSeg.SegmentDistance(seg.Eval(0.0));
                if(dist10 < eps)
                    return true;

                double dist11 = xSeg.SegmentDistance(seg.Eval(1.0));
                if(dist11 < eps)
                    return true;
            }
            else {
                // proper intersection of two lines
            
                double t = ((x10-x00)*dy1 - (y10-y00)*dx1)/denom;
                double u = ((x10-x00)*dy0 - (y10-y00)*dx0)/denom;

                // interior intersection
                if(0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0)
                    return true;
                else {
                    // clamp t, u to [0, 1]
                    // evaluate segments at t, u
                    // check if point distance is smaller than eps

                    t = std::max(0.0, std::min(1.0, t));
                    u = std::max(0.0, std::min(1.0, u));

                    auto pntA = seg.Eval(t);
                    auto pntB = xSeg.Eval(u);

                    if((pntA-pntB).norm() < eps)
                        return true;
                }
            }

            return false;
        }

        // access point vector
        const std::vector<Eigen::Vector2d> &points() const {
            return points_;
        }

        std::vector<Eigen::Vector2d> &points() {
            return points_;
        }

        // access triangle vector
        const std::vector<sTriangle> &triangles() const {
            return triangles_;
        }

        std::vector<sTriangle> &triangles() {
            return triangles_;
        }

        std::vector<Eigen::Vector2d> points_;

        std::vector<sTriangle> triangles_;
};
}; // namespace polygon2d
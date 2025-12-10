#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>

#include "ConvexifyPolygon.h"
#include "tools3d.h"

// unit tests

// a pentagon
TEST(ConvexifyPolyTest, Pentagon)
{
    polygon2d::ConvexifyPolygon convexPoly;

    convexPoly.push_back({1.000000, 0.000000});
    convexPoly.push_back({0.309017, 0.951057});
    convexPoly.push_back({-0.809017, 0.587785});
    convexPoly.push_back({-0.809017, -0.587785});
    convexPoly.push_back({0.309017, -0.951057});

    double area = convexPoly.computeConvexPartition();

    double checkArea = 2.377641895659;

    EXPECT_DOUBLE_EQ(area, checkArea);
}

// a triangle with vertices inserted in the middle of its edges
TEST(ConvexifyPolyTest, SubdivideTriangle)
{
    polygon2d::ConvexifyPolygon convexPoly;

    convexPoly.push_back({0, 0});
    convexPoly.push_back({0, 1});
    convexPoly.push_back({0, 2});
    convexPoly.push_back({1, 1});
    convexPoly.push_back({2, 0});
    convexPoly.push_back({1, 0});

    double area = convexPoly.computeConvexPartition();

    double checkArea = 2.0;

    EXPECT_DOUBLE_EQ(area, checkArea);
}

// invalid self-intersecting quad
TEST(ConvexifyPolyTest, BowtieQuad)
{
    polygon2d::ConvexifyPolygon convexPoly;

    convexPoly.push_back({0, 0});
    convexPoly.push_back({0, 1});
    convexPoly.push_back({1, 0});
    convexPoly.push_back({1, 1});

    double area = convexPoly.computeConvexPartition();

    double checkArea = -1.0;

    EXPECT_DOUBLE_EQ(area, checkArea);
}

// concave pentagon with diagonal intersecting edge
TEST(ConvexifyPolyTest, BootShapedConcave)
{
    polygon2d::ConvexifyPolygon convexPoly;

    convexPoly.push_back({0, 1});
    convexPoly.push_back({0, 0});
    convexPoly.push_back({1, 0});
    convexPoly.push_back({0.5, 0.1});
    convexPoly.push_back({0.5, 1});

    double area = convexPoly.computeConvexPartition();

    double checkArea = 0.525;

    EXPECT_DOUBLE_EQ(area, checkArea);
}

// triangle area and distance to a point
TEST(Tools3dTest, Triangle3dTest)
{
    Eigen::Vector3d pntA(0,0,0);
    Eigen::Vector3d pntB(1, 2, 3);
    Eigen::Vector3d pntC(4, -5, 6);
    Eigen::Vector3d pntD(2, 2, 2);

    Triangle3d triangle(pntA, pntB, pntC);

    auto closestPnt = triangle.PointDistanceSq(pntD);

    printf("%.6f\n", closestPnt.second);

    double checkDist = 1.2987012987012991;
    EXPECT_DOUBLE_EQ(closestPnt.second, checkDist);

    double checkArea = 30.561413579872251;
    double area = triangle.GetArea();

    EXPECT_DOUBLE_EQ(area, checkArea);
}

int main(int argc, char * argv[])
{
    testing::InitGoogleTest(&argc, argv);

    auto testStatus = RUN_ALL_TESTS();
    if (testStatus) {
        std::cerr << "Error: Unit test failure\n";
        exit(1);
    }

    exit(0);

    // TODO : plug here
    // SDFTest(nIter, mesh, tree);

    // OctreeTest(nIter, mesh, tree);

    // SimVolumeTest(nIter, nCycles, bounds, step, dT, ff);

    // SimulateCloth(mesh);

    // SimulateSandCastle(mesh, step, nIter);
}

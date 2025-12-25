#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

#include <gtest/gtest.h>
#include <eigen3/Eigen/Core>

#include "common.h"
#include "AABBTree.h"
#include "cGrid3d.h"
#include "cOctree.h"
#include "cSimVolume.h"
#include "MeshTools.h"
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
    Eigen::Vector3d pntA(0, 0, 3);
    Eigen::Vector3d pntB(0, 1, 3);
    Eigen::Vector3d pntC(1, 1, 3);
    Eigen::Vector3d pntD(0.5, 0.5, 2);

    Triangle3d triangle(pntA, pntB, pntC);

    auto closestPnt = triangle.PointDistanceSq(pntD);

    printf("%.6f\n", closestPnt.second);

    double checkDist = 1.0;
    EXPECT_DOUBLE_EQ(closestPnt.second, checkDist);

    double checkArea = 0.5;
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

    int nIter = argc > 2 ? atoi(argv[2]) : 100;
    int nCycles = argc > 3 ? atoi(argv[3]) : 3;
    double step = argc > 4 ? atof(argv[4]) : 5.0;
    double dT = argc > 5 ? atof(argv[5]) : 0.005;
    double ff = argc > 6 ? atof(argv[6]) : 1.0;

    MyMesh mesh;
    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.request_vertex_colors();
    // IO::Options ropt;

    // -------------------- read mesh

    // TODO : extract into OpenMeshMP

    printf("read mesh %s\n", argv[1]);

    OpenMesh::IO::Options ropt;
    ropt += OpenMesh::IO::Options::Binary;
    ropt += OpenMesh::IO::Options::VertexColor;

    if (!OpenMesh::IO::read_mesh(mesh, argv[1], ropt))
    {
        std::cerr << "Error loading mesh from file " << std::endl;
        return 1;
    }

    mesh.update_normals();

    AABBTree tree;
    tree.Build(mesh, 20);

    auto bounds = tree.GetBounds();

    ShapeDiameterFuncTest(mesh, tree, nIter);

    OctreeTest(mesh, tree, nIter);

    // SimVolumeTest(nIter, nCycles, bounds, step, dT, ff);

    // SimulateCloth(mesh);

    SimulateSandCastle(mesh, step, nIter);

    exit(0);
}

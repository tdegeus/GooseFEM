
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::MeshQuad4::Regular", "MeshQuad4.h")
{
    SECTION("Regular")
    {
        xt::xtensor<double, 2> coor_ = {{0., 0.}, {1., 0.}, {2., 0.}, {3., 0.}, {4., 0.}, {5., 0.},
                                        {0., 1.}, {1., 1.}, {2., 1.}, {3., 1.}, {4., 1.}, {5., 1.},
                                        {0., 2.}, {1., 2.}, {2., 2.}, {3., 2.}, {4., 2.}, {5., 2.},
                                        {0., 3.}, {1., 3.}, {2., 3.}, {3., 3.}, {4., 3.}, {5., 3.}};

        xt::xtensor<double, 2> conn_ = {{0, 1, 7, 6},
                                        {1, 2, 8, 7},
                                        {2, 3, 9, 8},
                                        {3, 4, 10, 9},
                                        {4, 5, 11, 10},
                                        {6, 7, 13, 12},
                                        {7, 8, 14, 13},
                                        {8, 9, 15, 14},
                                        {9, 10, 16, 15},
                                        {10, 11, 17, 16},
                                        {12, 13, 19, 18},
                                        {13, 14, 20, 19},
                                        {14, 15, 21, 20},
                                        {15, 16, 22, 21},
                                        {16, 17, 23, 22}};

        size_t nelem_ = 15;
        size_t nnode_ = 24;
        size_t nne_ = 4;
        size_t ndim_ = 2;
        size_t nnodePeriodic_ = 15;

        xt::xtensor<size_t, 1> nodesBottomEdge_ = {0, 1, 2, 3, 4, 5};
        xt::xtensor<size_t, 1> nodesTopEdge_ = {18, 19, 20, 21, 22, 23};
        xt::xtensor<size_t, 1> nodesLeftEdge_ = {0, 6, 12, 18};
        xt::xtensor<size_t, 1> nodesRightEdge_ = {5, 11, 17, 23};
        xt::xtensor<size_t, 1> nodesBottomOpenEdge_ = {1, 2, 3, 4};
        xt::xtensor<size_t, 1> nodesTopOpenEdge_ = {19, 20, 21, 22};
        xt::xtensor<size_t, 1> nodesLeftOpenEdge_ = {6, 12};
        xt::xtensor<size_t, 1> nodesRightOpenEdge_ = {11, 17};

        size_t nodesBottomLeftCorner_ = 0;
        size_t nodesBottomRightCorner_ = 5;
        size_t nodesTopLeftCorner_ = 18;
        size_t nodesTopRightCorner_ = 23;
        size_t nodesLeftBottomCorner_ = 0;
        size_t nodesLeftTopCorner_ = 18;
        size_t nodesRightBottomCorner_ = 5;
        size_t nodesRightTopCorner_ = 23;

        xt::xtensor<size_t, 2> dofs_ = {{0, 1},   {2, 3},   {4, 5},   {6, 7},   {8, 9},   {10, 11},
                                        {12, 13}, {14, 15}, {16, 17}, {18, 19}, {20, 21}, {22, 23},
                                        {24, 25}, {26, 27}, {28, 29}, {30, 31}, {32, 33}, {34, 35},
                                        {36, 37}, {38, 39}, {40, 41}, {42, 43}, {44, 45}, {46, 47}};

        xt::xtensor<size_t, 2> nodesPeriodic_ = {
            {0, 5}, {0, 23}, {0, 18}, {1, 19}, {2, 20}, {3, 21}, {4, 22}, {6, 11}, {12, 17}};

        size_t nodesOrigin_ = 0;

        xt::xtensor<size_t, 2> dofsPeriodic_ = {
            {0, 1},   {2, 3},   {4, 5},   {6, 7},   {8, 9},   {0, 1},   {10, 11}, {12, 13},
            {14, 15}, {16, 17}, {18, 19}, {10, 11}, {20, 21}, {22, 23}, {24, 25}, {26, 27},
            {28, 29}, {20, 21}, {0, 1},   {2, 3},   {4, 5},   {6, 7},   {8, 9},   {0, 1}};

        GooseFEM::Mesh::Quad4::Regular mesh(5, 3);

        xt::xtensor<double, 2> coor = mesh.coor();
        xt::xtensor<size_t, 2> conn = mesh.conn();

        size_t nelem = mesh.nelem();
        size_t nnode = mesh.nnode();
        size_t nne = mesh.nne();
        size_t ndim = mesh.ndim();
        size_t nnodePeriodic = mesh.nnode() - mesh.nodesPeriodic().shape(0);

        xt::xtensor<size_t, 1> nodesBottomEdge = mesh.nodesBottomEdge();
        xt::xtensor<size_t, 1> nodesTopEdge = mesh.nodesTopEdge();
        xt::xtensor<size_t, 1> nodesLeftEdge = mesh.nodesLeftEdge();
        xt::xtensor<size_t, 1> nodesRightEdge = mesh.nodesRightEdge();
        xt::xtensor<size_t, 1> nodesBottomOpenEdge = mesh.nodesBottomOpenEdge();
        xt::xtensor<size_t, 1> nodesTopOpenEdge = mesh.nodesTopOpenEdge();
        xt::xtensor<size_t, 1> nodesLeftOpenEdge = mesh.nodesLeftOpenEdge();
        xt::xtensor<size_t, 1> nodesRightOpenEdge = mesh.nodesRightOpenEdge();

        size_t nodesBottomLeftCorner = mesh.nodesBottomLeftCorner();
        size_t nodesBottomRightCorner = mesh.nodesBottomRightCorner();
        size_t nodesTopLeftCorner = mesh.nodesTopLeftCorner();
        size_t nodesTopRightCorner = mesh.nodesTopRightCorner();
        size_t nodesLeftBottomCorner = mesh.nodesLeftBottomCorner();
        size_t nodesLeftTopCorner = mesh.nodesLeftTopCorner();
        size_t nodesRightBottomCorner = mesh.nodesRightBottomCorner();
        size_t nodesRightTopCorner = mesh.nodesRightTopCorner();

        xt::xtensor<size_t, 2> dofs = mesh.dofs();

        xt::xtensor<size_t, 2> nodesPeriodic = mesh.nodesPeriodic();
        size_t nodesOrigin = mesh.nodesOrigin();
        xt::xtensor<size_t, 2> dofsPeriodic = mesh.dofsPeriodic();

        REQUIRE(nelem == nelem_);
        REQUIRE(nnode == nnode_);
        REQUIRE(nne == nne_);
        REQUIRE(ndim == ndim_);
        REQUIRE(nnodePeriodic == nnodePeriodic_);
        REQUIRE(nodesBottomLeftCorner == nodesBottomLeftCorner_);
        REQUIRE(nodesBottomRightCorner == nodesBottomRightCorner_);
        REQUIRE(nodesTopLeftCorner == nodesTopLeftCorner_);
        REQUIRE(nodesTopRightCorner == nodesTopRightCorner_);
        REQUIRE(nodesLeftBottomCorner == nodesLeftBottomCorner_);
        REQUIRE(nodesLeftTopCorner == nodesLeftTopCorner_);
        REQUIRE(nodesRightBottomCorner == nodesRightBottomCorner_);
        REQUIRE(nodesRightTopCorner == nodesRightTopCorner_);
        REQUIRE(nodesOrigin == nodesOrigin_);

        REQUIRE(xt::allclose(coor, coor_));

        REQUIRE(xt::all(xt::equal(conn, conn_)));
        REQUIRE(xt::all(xt::equal(nodesBottomEdge, nodesBottomEdge_)));
        REQUIRE(xt::all(xt::equal(nodesTopEdge, nodesTopEdge_)));
        REQUIRE(xt::all(xt::equal(nodesLeftEdge, nodesLeftEdge_)));
        REQUIRE(xt::all(xt::equal(nodesRightEdge, nodesRightEdge_)));
        REQUIRE(xt::all(xt::equal(nodesBottomOpenEdge, nodesBottomOpenEdge_)));
        REQUIRE(xt::all(xt::equal(nodesTopOpenEdge, nodesTopOpenEdge_)));
        REQUIRE(xt::all(xt::equal(nodesLeftOpenEdge, nodesLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(nodesRightOpenEdge, nodesRightOpenEdge_)));
        REQUIRE(xt::all(xt::equal(nodesPeriodic, nodesPeriodic_)));
        REQUIRE(xt::all(xt::equal(dofsPeriodic, dofsPeriodic_)));
    }
}

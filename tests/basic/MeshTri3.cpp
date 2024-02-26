#include <GooseFEM/MeshTri3.h>
#include <catch2/catch_all.hpp>

TEST_CASE("GooseFEM::MeshTri3", "MeshTri3.h")
{
    SECTION("Regular")
    {
        xt::xtensor<double, 2> coor_ = {{0., 0.}, {1., 0.}, {2., 0.}, {3., 0.}, {4., 0.}, {5., 0.},
                                        {0., 1.}, {1., 1.}, {2., 1.}, {3., 1.}, {4., 1.}, {5., 1.},
                                        {0., 2.}, {1., 2.}, {2., 2.}, {3., 2.}, {4., 2.}, {5., 2.},
                                        {0., 3.}, {1., 3.}, {2., 3.}, {3., 3.}, {4., 3.}, {5., 3.}};

        xt::xtensor<double, 2> conn_ = {
            {0, 1, 6},    {1, 7, 6},    {1, 2, 7},    {2, 8, 7},    {2, 3, 8},    {3, 9, 8},
            {3, 4, 9},    {4, 10, 9},   {4, 5, 10},   {5, 11, 10},  {6, 7, 12},   {7, 13, 12},
            {7, 8, 13},   {8, 14, 13},  {8, 9, 14},   {9, 15, 14},  {9, 10, 15},  {10, 16, 15},
            {10, 11, 16}, {11, 17, 16}, {12, 13, 18}, {13, 19, 18}, {13, 14, 19}, {14, 20, 19},
            {14, 15, 20}, {15, 21, 20}, {15, 16, 21}, {16, 22, 21}, {16, 17, 22}, {17, 23, 22}
        };

        size_t nelem_ = 15 * 2;
        size_t nnode_ = 24;
        size_t nne_ = 3;
        size_t ndim_ = 2;

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
            {0, 5}, {0, 23}, {0, 18}, {1, 19}, {2, 20}, {3, 21}, {4, 22}, {6, 11}, {12, 17}
        };

        size_t nodesOrigin_ = 0;

        xt::xtensor<size_t, 2> dofsPeriodic_ = {{0, 1},   {2, 3},   {4, 5},   {6, 7},   {8, 9},
                                                {0, 1},   {10, 11}, {12, 13}, {14, 15}, {16, 17},
                                                {18, 19}, {10, 11}, {20, 21}, {22, 23}, {24, 25},
                                                {26, 27}, {28, 29}, {20, 21}, {0, 1},   {2, 3},
                                                {4, 5},   {6, 7},   {8, 9},   {0, 1}};

        GooseFEM::Mesh::Tri3::Regular mesh(5, 3);

        REQUIRE(mesh.nelem() == nelem_);
        REQUIRE(mesh.nnode() == nnode_);
        REQUIRE(mesh.nne() == nne_);
        REQUIRE(mesh.ndim() == ndim_);
        REQUIRE(mesh.nodesBottomLeftCorner() == nodesBottomLeftCorner_);
        REQUIRE(mesh.nodesBottomRightCorner() == nodesBottomRightCorner_);
        REQUIRE(mesh.nodesTopLeftCorner() == nodesTopLeftCorner_);
        REQUIRE(mesh.nodesTopRightCorner() == nodesTopRightCorner_);
        REQUIRE(mesh.nodesLeftBottomCorner() == nodesLeftBottomCorner_);
        REQUIRE(mesh.nodesLeftTopCorner() == nodesLeftTopCorner_);
        REQUIRE(mesh.nodesRightBottomCorner() == nodesRightBottomCorner_);
        REQUIRE(mesh.nodesRightTopCorner() == nodesRightTopCorner_);
        REQUIRE(mesh.nodesOrigin() == nodesOrigin_);
        REQUIRE(xt::allclose(mesh.coor(), coor_));
        REQUIRE(xt::all(xt::equal(mesh.conn(), conn_)));
        REQUIRE(xt::all(xt::equal(mesh.dofs(), dofs_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomEdge(), nodesBottomEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopEdge(), nodesTopEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesLeftEdge(), nodesLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesRightEdge(), nodesRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomOpenEdge(), nodesBottomOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopOpenEdge(), nodesTopOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesLeftOpenEdge(), nodesLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesRightOpenEdge(), nodesRightOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesPeriodic(), nodesPeriodic_)));
        REQUIRE(xt::all(xt::equal(mesh.dofsPeriodic(), dofsPeriodic_)));
    }
}

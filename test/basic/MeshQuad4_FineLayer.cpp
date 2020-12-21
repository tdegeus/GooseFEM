
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::MeshQuad4::FineLayer", "MeshQuad4.h")
{
    SECTION("FineLayer")
    {
        xt::xtensor<double, 2> coor_ = {
            {0., 0.},  {3., 0.},  {6., 0.},  {9., 0.},  {0., 3.},  {3., 3.},  {6., 3.},  {9., 3.},
            {0., 6.},  {3., 6.},  {6., 6.},  {9., 6.},  {1., 7.},  {2., 7.},  {4., 7.},  {5., 7.},
            {7., 7.},  {8., 7.},  {0., 8.},  {1., 8.},  {2., 8.},  {3., 8.},  {4., 8.},  {5., 8.},
            {6., 8.},  {7., 8.},  {8., 8.},  {9., 8.},  {0., 9.},  {1., 9.},  {2., 9.},  {3., 9.},
            {4., 9.},  {5., 9.},  {6., 9.},  {7., 9.},  {8., 9.},  {9., 9.},  {0., 10.}, {1., 10.},
            {2., 10.}, {3., 10.}, {4., 10.}, {5., 10.}, {6., 10.}, {7., 10.}, {8., 10.}, {9., 10.},
            {0., 11.}, {1., 11.}, {2., 11.}, {3., 11.}, {4., 11.}, {5., 11.}, {6., 11.}, {7., 11.},
            {8., 11.}, {9., 11.}, {0., 12.}, {1., 12.}, {2., 12.}, {3., 12.}, {4., 12.}, {5., 12.},
            {6., 12.}, {7., 12.}, {8., 12.}, {9., 12.}, {0., 13.}, {1., 13.}, {2., 13.}, {3., 13.},
            {4., 13.}, {5., 13.}, {6., 13.}, {7., 13.}, {8., 13.}, {9., 13.}, {1., 14.}, {2., 14.},
            {4., 14.}, {5., 14.}, {7., 14.}, {8., 14.}, {0., 15.}, {3., 15.}, {6., 15.}, {9., 15.},
            {0., 18.}, {3., 18.}, {6., 18.}, {9., 18.}, {0., 21.}, {3., 21.}, {6., 21.}, {9., 21.}};

        xt::xtensor<double, 2> conn_ = {
            {0, 1, 5, 4},     {1, 2, 6, 5},     {2, 3, 7, 6},     {4, 5, 9, 8},
            {5, 6, 10, 9},    {6, 7, 11, 10},   {8, 9, 13, 12},   {9, 21, 20, 13},
            {12, 13, 20, 19}, {8, 12, 19, 18},  {9, 10, 15, 14},  {10, 24, 23, 15},
            {14, 15, 23, 22}, {9, 14, 22, 21},  {10, 11, 17, 16}, {11, 27, 26, 17},
            {16, 17, 26, 25}, {10, 16, 25, 24}, {18, 19, 29, 28}, {19, 20, 30, 29},
            {20, 21, 31, 30}, {21, 22, 32, 31}, {22, 23, 33, 32}, {23, 24, 34, 33},
            {24, 25, 35, 34}, {25, 26, 36, 35}, {26, 27, 37, 36}, {28, 29, 39, 38},
            {29, 30, 40, 39}, {30, 31, 41, 40}, {31, 32, 42, 41}, {32, 33, 43, 42},
            {33, 34, 44, 43}, {34, 35, 45, 44}, {35, 36, 46, 45}, {36, 37, 47, 46},
            {38, 39, 49, 48}, {39, 40, 50, 49}, {40, 41, 51, 50}, {41, 42, 52, 51},
            {42, 43, 53, 52}, {43, 44, 54, 53}, {44, 45, 55, 54}, {45, 46, 56, 55},
            {46, 47, 57, 56}, {48, 49, 59, 58}, {49, 50, 60, 59}, {50, 51, 61, 60},
            {51, 52, 62, 61}, {52, 53, 63, 62}, {53, 54, 64, 63}, {54, 55, 65, 64},
            {55, 56, 66, 65}, {56, 57, 67, 66}, {58, 59, 69, 68}, {59, 60, 70, 69},
            {60, 61, 71, 70}, {61, 62, 72, 71}, {62, 63, 73, 72}, {63, 64, 74, 73},
            {64, 65, 75, 74}, {65, 66, 76, 75}, {66, 67, 77, 76}, {68, 69, 78, 84},
            {69, 70, 79, 78}, {70, 71, 85, 79}, {78, 79, 85, 84}, {71, 72, 80, 85},
            {72, 73, 81, 80}, {73, 74, 86, 81}, {80, 81, 86, 85}, {74, 75, 82, 86},
            {75, 76, 83, 82}, {76, 77, 87, 83}, {82, 83, 87, 86}, {84, 85, 89, 88},
            {85, 86, 90, 89}, {86, 87, 91, 90}, {88, 89, 93, 92}, {89, 90, 94, 93},
            {90, 91, 95, 94}};

        size_t nelem_ = 81;
        size_t nnode_ = 96;
        size_t nne_ = 4;
        size_t ndim_ = 2;
        size_t shape_x_ = 9;
        size_t shape_y_ = 21;

        xt::xtensor<size_t, 1> elementsMiddleLayer_ = {36, 37, 38, 39, 40, 41, 42, 43, 44};

        xt::xtensor<size_t, 1> nodesBottomEdge_ = {0, 1, 2, 3};
        xt::xtensor<size_t, 1> nodesTopEdge_ = {92, 93, 94, 95};
        xt::xtensor<size_t, 1> nodesLeftEdge_ = {0, 4, 8, 18, 28, 38, 48, 58, 68, 84, 88, 92};
        xt::xtensor<size_t, 1> nodesRightEdge_ = {3, 7, 11, 27, 37, 47, 57, 67, 77, 87, 91, 95};
        xt::xtensor<size_t, 1> nodesBottomOpenEdge_ = {1, 2};
        xt::xtensor<size_t, 1> nodesTopOpenEdge_ = {93, 94};
        xt::xtensor<size_t, 1> nodesLeftOpenEdge_ = {4, 8, 18, 28, 38, 48, 58, 68, 84, 88};
        xt::xtensor<size_t, 1> nodesRightOpenEdge_ = {7, 11, 27, 37, 47, 57, 67, 77, 87, 91};

        size_t nodesBottomLeftCorner_ = 0;
        size_t nodesBottomRightCorner_ = 3;
        size_t nodesTopLeftCorner_ = 92;
        size_t nodesTopRightCorner_ = 95;
        size_t nodesLeftBottomCorner_ = 0;
        size_t nodesLeftTopCorner_ = 92;
        size_t nodesRightBottomCorner_ = 3;
        size_t nodesRightTopCorner_ = 95;

        xt::xtensor<size_t, 2> dofs_ = {
            {0, 1},     {2, 3},     {4, 5},     {6, 7},     {8, 9},     {10, 11},   {12, 13},
            {14, 15},   {16, 17},   {18, 19},   {20, 21},   {22, 23},   {24, 25},   {26, 27},
            {28, 29},   {30, 31},   {32, 33},   {34, 35},   {36, 37},   {38, 39},   {40, 41},
            {42, 43},   {44, 45},   {46, 47},   {48, 49},   {50, 51},   {52, 53},   {54, 55},
            {56, 57},   {58, 59},   {60, 61},   {62, 63},   {64, 65},   {66, 67},   {68, 69},
            {70, 71},   {72, 73},   {74, 75},   {76, 77},   {78, 79},   {80, 81},   {82, 83},
            {84, 85},   {86, 87},   {88, 89},   {90, 91},   {92, 93},   {94, 95},   {96, 97},
            {98, 99},   {100, 101}, {102, 103}, {104, 105}, {106, 107}, {108, 109}, {110, 111},
            {112, 113}, {114, 115}, {116, 117}, {118, 119}, {120, 121}, {122, 123}, {124, 125},
            {126, 127}, {128, 129}, {130, 131}, {132, 133}, {134, 135}, {136, 137}, {138, 139},
            {140, 141}, {142, 143}, {144, 145}, {146, 147}, {148, 149}, {150, 151}, {152, 153},
            {154, 155}, {156, 157}, {158, 159}, {160, 161}, {162, 163}, {164, 165}, {166, 167},
            {168, 169}, {170, 171}, {172, 173}, {174, 175}, {176, 177}, {178, 179}, {180, 181},
            {182, 183}, {184, 185}, {186, 187}, {188, 189}, {190, 191}};

        xt::xtensor<size_t, 2> nodesPeriodic_ = {{0, 3},
                                                 {0, 95},
                                                 {0, 92},
                                                 {1, 93},
                                                 {2, 94},
                                                 {4, 7},
                                                 {8, 11},
                                                 {18, 27},
                                                 {28, 37},
                                                 {38, 47},
                                                 {48, 57},
                                                 {58, 67},
                                                 {68, 77},
                                                 {84, 87},
                                                 {88, 91}};

        size_t nodesOrigin_ = 0;

        xt::xtensor<size_t, 2> dofsPeriodic_ = {
            {0, 1},     {2, 3},     {4, 5},     {0, 1},     {6, 7},     {8, 9},     {10, 11},
            {6, 7},     {12, 13},   {14, 15},   {16, 17},   {12, 13},   {18, 19},   {20, 21},
            {22, 23},   {24, 25},   {26, 27},   {28, 29},   {30, 31},   {32, 33},   {34, 35},
            {36, 37},   {38, 39},   {40, 41},   {42, 43},   {44, 45},   {46, 47},   {30, 31},
            {48, 49},   {50, 51},   {52, 53},   {54, 55},   {56, 57},   {58, 59},   {60, 61},
            {62, 63},   {64, 65},   {48, 49},   {66, 67},   {68, 69},   {70, 71},   {72, 73},
            {74, 75},   {76, 77},   {78, 79},   {80, 81},   {82, 83},   {66, 67},   {84, 85},
            {86, 87},   {88, 89},   {90, 91},   {92, 93},   {94, 95},   {96, 97},   {98, 99},
            {100, 101}, {84, 85},   {102, 103}, {104, 105}, {106, 107}, {108, 109}, {110, 111},
            {112, 113}, {114, 115}, {116, 117}, {118, 119}, {102, 103}, {120, 121}, {122, 123},
            {124, 125}, {126, 127}, {128, 129}, {130, 131}, {132, 133}, {134, 135}, {136, 137},
            {120, 121}, {138, 139}, {140, 141}, {142, 143}, {144, 145}, {146, 147}, {148, 149},
            {150, 151}, {152, 153}, {154, 155}, {150, 151}, {156, 157}, {158, 159}, {160, 161},
            {156, 157}, {0, 1},     {2, 3},     {4, 5},     {0, 1}};

        GooseFEM::Mesh::Quad4::FineLayer mesh(9, 18);

        xt::xtensor<double, 2> coor = mesh.coor();
        xt::xtensor<size_t, 2> conn = mesh.conn();

        size_t nelem = mesh.nelem();
        size_t nnode = mesh.nnode();
        size_t nne = mesh.nne();
        size_t ndim = mesh.ndim();
        size_t shape_x = mesh.nelx();
        size_t shape_y = mesh.nely();

        xt::xtensor<size_t, 1> elementsMiddleLayer = mesh.elementsMiddleLayer();

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
        REQUIRE(shape_x == shape_x_);
        REQUIRE(shape_y == shape_y_);
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

        REQUIRE(xt::all(xt::equal(elementsMiddleLayer, elementsMiddleLayer_)));
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

    SECTION("FineLayer::elementgrid_ravel - uniform")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(5, 5);
        xt::xtensor<size_t, 1> a = {
             0,  1,  2,  3,  4,
             5,  6,  7,  8,  9,
            10, 11, 12, 13, 14,
            15, 16, 17, 18, 19,
            20, 21, 22, 23, 24};
        REQUIRE(xt::all(xt::equal(a, mesh.elementgrid_ravel({0, 5}, {0, 5}))));
        REQUIRE(xt::all(xt::equal(a, mesh.elementgrid_ravel({0, 5}, {}))));
        REQUIRE(xt::all(xt::equal(a, mesh.elementgrid_ravel({}, {0, 5}))));
        REQUIRE(xt::all(xt::equal(a, mesh.elementgrid_ravel({}, {}))));

        xt::xtensor<size_t, 1> b = {
             5,  6,  7,  8,  9,
            10, 11, 12, 13, 14,
            15, 16, 17, 18, 19};
        REQUIRE(xt::all(xt::equal(b, mesh.elementgrid_ravel({1, 4}, {0, 5}))));

        xt::xtensor<size_t, 1> c = {
             6,  7,  8,
            11, 12, 13,
            16, 17, 18};
        REQUIRE(xt::all(xt::equal(c, mesh.elementgrid_ravel({1, 4}, {1, 4}))));
    }

    SECTION("FineLayer::elementgrid_ravel - refined")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(6, 18);
        xt::xtensor<size_t, 1> a = {
            19, 20, 21, 22,
            25, 26, 27, 28,
            31, 32, 33, 34};
        REQUIRE(xt::all(xt::equal(a, mesh.elementgrid_ravel({9, 12}, {1, 5}))));

        xt::xtensor<size_t, 1> b = {
            0, 1,
            2, 3};
        REQUIRE(xt::all(xt::equal(b, mesh.elementgrid_ravel({0, 6}, {0, 6}))));

        xt::xtensor<size_t, 1> c = {
            50, 51,
            52, 53};
        REQUIRE(xt::all(xt::equal(c, mesh.elementgrid_ravel({15, 21}, {0, 6}))));

        xt::xtensor<size_t, 1> d = {
            4, 5, 6, 7, 8, 9, 10, 11};
        REQUIRE(xt::all(xt::equal(d, mesh.elementgrid_ravel({6, 8}, {0, 6}))));
    }

    SECTION("FineLayer::elementgrid_ravel - uniform")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(5, 5);
        xt::xtensor<size_t, 1> a = {
             0,  1,  2,  3,  4,
             5,  6,  7,  8,  9,
            10, 11, 12, 13, 14,
            15, 16, 17, 18, 19,
            20, 21, 22, 23, 24};

        REQUIRE(xt::all(xt::equal(a, xt::sort(mesh.elementgrid_around_ravel(10, 2)))));
        REQUIRE(xt::all(xt::equal(a, xt::sort(mesh.elementgrid_around_ravel(11, 2)))));
        REQUIRE(xt::all(xt::equal(a, xt::sort(mesh.elementgrid_around_ravel(12, 2)))));
        REQUIRE(xt::all(xt::equal(a, xt::sort(mesh.elementgrid_around_ravel(13, 2)))));
        REQUIRE(xt::all(xt::equal(a, xt::sort(mesh.elementgrid_around_ravel(14, 2)))));

        xt::xtensor<size_t, 1> b10 = {
             5,  6,  9,
            10, 11, 14,
            15, 16, 19};

        xt::xtensor<size_t, 1> b11 = {
             5,  6,  7,
            10, 11, 12,
            15, 16, 17};

        xt::xtensor<size_t, 1> b12 = {
             6,  7,  8,
            11, 12, 13,
            16, 17, 18};

        xt::xtensor<size_t, 1> b13 = {
             7,  8,  9,
            12, 13, 14,
            17, 18, 19};

        xt::xtensor<size_t, 1> b14 = {
             5,  8,  9,
            10, 13, 14,
            15, 18, 19};

        REQUIRE(xt::all(xt::equal(xt::sort(b10), xt::sort(mesh.elementgrid_around_ravel(10, 1)))));
        REQUIRE(xt::all(xt::equal(xt::sort(b11), xt::sort(mesh.elementgrid_around_ravel(11, 1)))));
        REQUIRE(xt::all(xt::equal(xt::sort(b12), xt::sort(mesh.elementgrid_around_ravel(12, 1)))));
        REQUIRE(xt::all(xt::equal(xt::sort(b13), xt::sort(mesh.elementgrid_around_ravel(13, 1)))));
        REQUIRE(xt::all(xt::equal(xt::sort(b14), xt::sort(mesh.elementgrid_around_ravel(14, 1)))));

        REQUIRE(xt::all(xt::equal(10, xt::sort(mesh.elementgrid_around_ravel(10, 0)))));
        REQUIRE(xt::all(xt::equal(11, xt::sort(mesh.elementgrid_around_ravel(11, 0)))));
        REQUIRE(xt::all(xt::equal(12, xt::sort(mesh.elementgrid_around_ravel(12, 0)))));
        REQUIRE(xt::all(xt::equal(13, xt::sort(mesh.elementgrid_around_ravel(13, 0)))));
        REQUIRE(xt::all(xt::equal(14, xt::sort(mesh.elementgrid_around_ravel(14, 0)))));
    }

    SECTION("FineLayer::elementgrid_around_ravel")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(12, 18);
        xt::xtensor<size_t, 1> r0 = {
            36, 37, 71,
            48, 49, 59,
            60, 61, 47,
        };
        xt::xtensor<size_t, 1> r1 = {
            36, 37, 38,
            48, 49, 50,
            60, 61, 62,
        };
        xt::xtensor<size_t, 1> r12 = {
            46, 47, 36,
            58, 59, 48,
            70, 71, 60,
        };

        REQUIRE(xt::all(xt::equal(xt::sort(r0), xt::sort(mesh.elementgrid_around_ravel(48, 1)))));
        for (size_t n = 0; n < 10; ++n) {
            REQUIRE(xt::all(xt::equal(xt::sort(r1) + n, xt::sort(mesh.elementgrid_around_ravel(49 + n, 1)))));
        }
        REQUIRE(xt::all(xt::equal(xt::sort(r12), xt::sort(mesh.elementgrid_around_ravel(59, 1)))));
    }

    SECTION("FineLayer::elementgrid_leftright")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(12, 18);
        xt::xtensor<size_t, 1> r0 = {
            48, 49, 59,
        };
        xt::xtensor<size_t, 1> r1 = {
            48, 49, 50,
        };
        xt::xtensor<size_t, 1> r12 = {
            58, 59, 48,
        };

        REQUIRE(xt::all(xt::equal(xt::sort(r0), xt::sort(mesh.elementgrid_leftright(48, 1, 1)))));
        for (size_t n = 0; n < 10; ++n) {
            REQUIRE(xt::all(xt::equal(xt::sort(r1) + n, xt::sort(mesh.elementgrid_leftright(49 + n, 1, 1)))));
        }
        REQUIRE(xt::all(xt::equal(xt::sort(r12), xt::sort(mesh.elementgrid_leftright(59, 1, 1)))));
    }

    SECTION("FineLayer - replica - trivial")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(1, 1);
        GooseFEM::Mesh::Quad4::FineLayer replica(mesh.coor(), mesh.conn());
        REQUIRE(xt::all(xt::equal(mesh.conn(), replica.conn())));
        REQUIRE(xt::allclose(mesh.coor(), replica.coor()));
    }

    SECTION("FineLayer - replica - equi-sized")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(4, 4);
        GooseFEM::Mesh::Quad4::FineLayer replica(mesh.coor(), mesh.conn());
        REQUIRE(xt::all(xt::equal(mesh.conn(), replica.conn())));
        REQUIRE(xt::allclose(mesh.coor(), replica.coor()));
    }

    SECTION("FineLayer - replica")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(27, 27);
        GooseFEM::Mesh::Quad4::FineLayer replica(mesh.coor(), mesh.conn());
        REQUIRE(xt::all(xt::equal(mesh.conn(), replica.conn())));
        REQUIRE(xt::allclose(mesh.coor(), replica.coor()));
    }

    SECTION("FineLayer - replica - one")
    {
        double h = xt::numeric_constants<double>::PI;
        GooseFEM::Mesh::Quad4::FineLayer mesh(1, 1, h);
        GooseFEM::Mesh::Quad4::FineLayer replica(mesh.coor(), mesh.conn());
        REQUIRE(xt::all(xt::equal(mesh.conn(), replica.conn())));
        REQUIRE(xt::allclose(mesh.coor(), replica.coor()));
    }

    SECTION("FineLayer - replica - large")
    {
        size_t N = 2 * std::pow(3, 6);
        double h = xt::numeric_constants<double>::PI;
        GooseFEM::Mesh::Quad4::FineLayer mesh(N, N, h);
        GooseFEM::Mesh::Quad4::FineLayer replica(mesh.coor(), mesh.conn());
        REQUIRE(xt::all(xt::equal(mesh.conn(), replica.conn())));
        REQUIRE(xt::allclose(mesh.coor(), replica.coor()));
    }

    SECTION("FineLayer - roll - trivial")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(3, 3);
        xt::xtensor<size_t, 1> m0 = {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8
        };
        xt::xtensor<size_t, 1> m1 = {
            2, 0, 1,
            5, 3, 4,
            8, 6, 7
        };
        xt::xtensor<size_t, 1> m2 = {
            1, 2, 0,
            4, 5, 3,
            7, 8, 6,
        };
        REQUIRE(xt::all(xt::equal(mesh.roll(0), m0)));
        REQUIRE(xt::all(xt::equal(mesh.roll(1), m1)));
        REQUIRE(xt::all(xt::equal(mesh.roll(2), m2)));
        REQUIRE(xt::all(xt::equal(mesh.roll(3), m0)));
        REQUIRE(xt::all(xt::equal(mesh.roll(4), m1)));
        REQUIRE(xt::all(xt::equal(mesh.roll(5), m2)));
    }

    SECTION("FineLayer - roll - a")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(6, 18);
        xt::xtensor<size_t, 1> m0 = {
            0, 1,
            2, 3,
            4, 5, 6, 7, 8, 9, 10, 11,
            12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29,
            30, 31, 32, 33, 34, 35,
            36, 37, 38, 39, 40, 41,
            42, 43, 44, 45, 46, 47, 48, 49,
            50, 51,
            52, 53
        };
        xt::xtensor<size_t, 1> m1 = {
            1, 0,
            3, 2,
            8, 9, 10, 11, 4, 5, 6, 7,
            15, 16, 17, 12, 13, 14,
            21, 22, 23, 18, 19, 20,
            27, 28, 29, 24, 25, 26,
            33, 34, 35, 30, 31, 32,
            39, 40, 41, 36, 37, 38,
            46, 47, 48, 49, 42, 43, 44, 45,
            51, 50,
            53, 52
        };
        REQUIRE(xt::all(xt::equal(mesh.roll(0), m0)));
        REQUIRE(xt::all(xt::equal(mesh.roll(1), m1)));
        REQUIRE(xt::all(xt::equal(mesh.roll(2), m0)));
        REQUIRE(xt::all(xt::equal(mesh.roll(3), m1)));
    }

    SECTION("FineLayer - roll - b")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(9, 18);
        xt::xtensor<size_t, 1> m0 = {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 22, 23, 24, 25, 26,
            27, 28, 29, 30, 31, 32, 33, 34, 35,
            36, 37, 38, 39, 40, 41, 42, 43, 44,
            45, 46, 47, 48, 49, 50, 51, 52, 53,
            54, 55, 56, 57, 58, 59, 60, 61, 62,
            63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
            75, 76, 77,
            78, 79, 80
        };
        xt::xtensor<size_t, 1> m1 = {
            2, 0, 1,
            5, 3, 4,
            14, 15, 16, 17, 6, 7, 8, 9, 10, 11, 12, 13,
            24, 25, 26, 18, 19, 20, 21, 22, 23,
            33, 34, 35, 27, 28, 29, 30, 31, 32,
            42, 43, 44, 36, 37, 38, 39, 40, 41,
            51, 52, 53, 45, 46, 47, 48, 49, 50,
            60, 61, 62, 54, 55, 56, 57, 58, 59,
            71, 72, 73, 74, 63, 64, 65, 66, 67, 68, 69, 70,
            77, 75, 76,
            80, 78, 79
        };
        xt::xtensor<size_t, 1> m2 = {
            1, 2, 0,
            4, 5, 3,
            10, 11, 12, 13, 14, 15, 16, 17, 6, 7, 8, 9,
            21, 22, 23, 24, 25, 26, 18, 19, 20,
            30, 31, 32, 33, 34, 35, 27, 28, 29,
            39, 40, 41, 42, 43, 44, 36, 37, 38,
            48, 49, 50, 51, 52, 53, 45, 46, 47,
            57, 58, 59, 60, 61, 62, 54, 55, 56,
            67, 68, 69, 70, 71, 72, 73, 74, 63, 64, 65, 66,
            76, 77, 75,
            79, 80, 78
        };
        REQUIRE(xt::all(xt::equal(mesh.roll(0), m0)));
        REQUIRE(xt::all(xt::equal(mesh.roll(1), m1)));
        REQUIRE(xt::all(xt::equal(mesh.roll(2), m2)));
        REQUIRE(xt::all(xt::equal(mesh.roll(3), m0)));
        REQUIRE(xt::all(xt::equal(mesh.roll(4), m1)));
        REQUIRE(xt::all(xt::equal(mesh.roll(5), m2)));
    }

    SECTION("Map::RefineRegular")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(5, 4);

        GooseFEM::Mesh::Quad4::Map::RefineRegular refine(mesh, 5, 3);

        xt::xtensor<double, 1> a = xt::random::rand<double>({mesh.nelem()});
        auto a_ = refine.mapToCoarse(refine.mapToFine(a));

        REQUIRE(xt::allclose(a, xt::mean(a_, {1})));

        xt::xtensor<double, 2> b =
            xt::random::rand<double>(std::array<size_t, 2>{mesh.nelem(), 4ul});
        auto b_ = refine.mapToCoarse(refine.mapToFine(b));

        REQUIRE(xt::allclose(xt::mean(b, {1}), xt::mean(b_, {1})));

        xt::xtensor<double, 4> c =
            xt::random::rand<double>(std::array<size_t, 4>{mesh.nelem(), 4ul, 3ul, 3ul});
        auto c_ = refine.mapToCoarse(refine.mapToFine(c));

        REQUIRE(xt::allclose(xt::mean(c, {1}), xt::mean(c_, {1})));
    }

    SECTION("Map::FineLayer2Regular")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(5, 5);

        GooseFEM::Mesh::Quad4::Map::FineLayer2Regular map(mesh);

        xt::xtensor<double, 1> a = xt::random::rand<double>({mesh.nelem()});
        auto a_ = map.mapToRegular(a);

        REQUIRE(xt::allclose(a, a_));

        xt::xtensor<double, 2> b =
            xt::random::rand<double>(std::array<size_t, 2>{mesh.nelem(), 4ul});
        auto b_ = map.mapToRegular(b);

        REQUIRE(xt::allclose(b, b_));

        xt::xtensor<double, 4> c =
            xt::random::rand<double>(std::array<size_t, 4>{mesh.nelem(), 4ul, 3ul, 3ul});
        auto c_ = map.mapToRegular(c);

        REQUIRE(xt::allclose(c, c_));
    }
}

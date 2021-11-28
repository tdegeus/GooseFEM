#define CATCH_CONFIG_MAIN // tells Catch to provide a main() - only do this in one cpp file
#include <GooseFEM/MeshHex8.h>
#include <catch2/catch.hpp>

TEST_CASE("GooseFEM::MeshHex8", "MeshHex8.h")
{
    SECTION("Regular")
    {
        xt::xtensor<double, 2> coor_ = {
            {0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}, {3., 0., 0.}, {0., 1., 0.}, {1., 1., 0.},
            {2., 1., 0.}, {3., 1., 0.}, {0., 2., 0.}, {1., 2., 0.}, {2., 2., 0.}, {3., 2., 0.},
            {0., 3., 0.}, {1., 3., 0.}, {2., 3., 0.}, {3., 3., 0.}, {0., 0., 1.}, {1., 0., 1.},
            {2., 0., 1.}, {3., 0., 1.}, {0., 1., 1.}, {1., 1., 1.}, {2., 1., 1.}, {3., 1., 1.},
            {0., 2., 1.}, {1., 2., 1.}, {2., 2., 1.}, {3., 2., 1.}, {0., 3., 1.}, {1., 3., 1.},
            {2., 3., 1.}, {3., 3., 1.}, {0., 0., 2.}, {1., 0., 2.}, {2., 0., 2.}, {3., 0., 2.},
            {0., 1., 2.}, {1., 1., 2.}, {2., 1., 2.}, {3., 1., 2.}, {0., 2., 2.}, {1., 2., 2.},
            {2., 2., 2.}, {3., 2., 2.}, {0., 3., 2.}, {1., 3., 2.}, {2., 3., 2.}, {3., 3., 2.},
            {0., 0., 3.}, {1., 0., 3.}, {2., 0., 3.}, {3., 0., 3.}, {0., 1., 3.}, {1., 1., 3.},
            {2., 1., 3.}, {3., 1., 3.}, {0., 2., 3.}, {1., 2., 3.}, {2., 2., 3.}, {3., 2., 3.},
            {0., 3., 3.}, {1., 3., 3.}, {2., 3., 3.}, {3., 3., 3.}};

        xt::xtensor<double, 2> conn_ = {
            {0, 1, 5, 4, 16, 17, 21, 20},     {1, 2, 6, 5, 17, 18, 22, 21},
            {2, 3, 7, 6, 18, 19, 23, 22},     {4, 5, 9, 8, 20, 21, 25, 24},
            {5, 6, 10, 9, 21, 22, 26, 25},    {6, 7, 11, 10, 22, 23, 27, 26},
            {8, 9, 13, 12, 24, 25, 29, 28},   {9, 10, 14, 13, 25, 26, 30, 29},
            {10, 11, 15, 14, 26, 27, 31, 30}, {16, 17, 21, 20, 32, 33, 37, 36},
            {17, 18, 22, 21, 33, 34, 38, 37}, {18, 19, 23, 22, 34, 35, 39, 38},
            {20, 21, 25, 24, 36, 37, 41, 40}, {21, 22, 26, 25, 37, 38, 42, 41},
            {22, 23, 27, 26, 38, 39, 43, 42}, {24, 25, 29, 28, 40, 41, 45, 44},
            {25, 26, 30, 29, 41, 42, 46, 45}, {26, 27, 31, 30, 42, 43, 47, 46},
            {32, 33, 37, 36, 48, 49, 53, 52}, {33, 34, 38, 37, 49, 50, 54, 53},
            {34, 35, 39, 38, 50, 51, 55, 54}, {36, 37, 41, 40, 52, 53, 57, 56},
            {37, 38, 42, 41, 53, 54, 58, 57}, {38, 39, 43, 42, 54, 55, 59, 58},
            {40, 41, 45, 44, 56, 57, 61, 60}, {41, 42, 46, 45, 57, 58, 62, 61},
            {42, 43, 47, 46, 58, 59, 63, 62}};

        xt::xtensor<size_t, 1> nodesBottom_ = {
            0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51};

        xt::xtensor<size_t, 1> nodesTop_ = {
            12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63};

        xt::xtensor<size_t, 1> nodesLeft_ = {
            0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60};

        xt::xtensor<size_t, 1> nodesRight_ = {
            3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63};

        xt::xtensor<size_t, 1> nodesFront_ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

        xt::xtensor<size_t, 1> nodesBack_ = {
            48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};

        xt::xtensor<size_t, 1> nodesFrontTopEdge_ = {12, 13, 14, 15};
        xt::xtensor<size_t, 1> nodesFrontLeftEdge_ = {0, 4, 8, 12};
        xt::xtensor<size_t, 1> nodesFrontRightEdge_ = {3, 7, 11, 15};
        xt::xtensor<size_t, 1> nodesBackBottomEdge_ = {48, 49, 50, 51};
        xt::xtensor<size_t, 1> nodesBackTopEdge_ = {60, 61, 62, 63};
        xt::xtensor<size_t, 1> nodesBackLeftEdge_ = {48, 52, 56, 60};
        xt::xtensor<size_t, 1> nodesBackRightEdge_ = {51, 55, 59, 63};
        xt::xtensor<size_t, 1> nodesBottomLeftEdge_ = {0, 16, 32, 48};
        xt::xtensor<size_t, 1> nodesBottomRightEdge_ = {3, 19, 35, 51};
        xt::xtensor<size_t, 1> nodesTopLeftEdge_ = {12, 28, 44, 60};
        xt::xtensor<size_t, 1> nodesTopRightEdge_ = {15, 31, 47, 63};
        xt::xtensor<size_t, 1> nodesFrontFace_ = {5, 6, 9, 10};
        xt::xtensor<size_t, 1> nodesBackFace_ = {53, 54, 57, 58};
        xt::xtensor<size_t, 1> nodesLeftFace_ = {20, 24, 36, 40};
        xt::xtensor<size_t, 1> nodesRightFace_ = {23, 27, 39, 43};
        xt::xtensor<size_t, 1> nodesBottomFace_ = {17, 18, 33, 34};
        xt::xtensor<size_t, 1> nodesTopFace_ = {29, 30, 45, 46};
        xt::xtensor<size_t, 1> nodesFrontBottomOpenEdge_ = {1, 2};
        xt::xtensor<size_t, 1> nodesFrontTopOpenEdge_ = {13, 14};
        xt::xtensor<size_t, 1> nodesFrontLeftOpenEdge_ = {4, 8};
        xt::xtensor<size_t, 1> nodesFrontRightOpenEdge_ = {7, 11};
        xt::xtensor<size_t, 1> nodesBackBottomOpenEdge_ = {49, 50};
        xt::xtensor<size_t, 1> nodesBackTopOpenEdge_ = {61, 62};
        xt::xtensor<size_t, 1> nodesBackLeftOpenEdge_ = {52, 56};
        xt::xtensor<size_t, 1> nodesBackRightOpenEdge_ = {55, 59};
        xt::xtensor<size_t, 1> nodesBottomLeftOpenEdge_ = {16, 32};
        xt::xtensor<size_t, 1> nodesBottomRightOpenEdge_ = {19, 35};
        xt::xtensor<size_t, 1> nodesTopLeftOpenEdge_ = {28, 44};
        xt::xtensor<size_t, 1> nodesTopRightOpenEdge_ = {31, 47};
        size_t nodesFrontBottomLeftCorner_ = 0;
        size_t nodesFrontBottomRightCorner_ = 3;
        size_t nodesFrontTopLeftCorner_ = 12;
        size_t nodesFrontTopRightCorner_ = 15;
        size_t nodesBackBottomLeftCorner_ = 48;
        size_t nodesBackBottomRightCorner_ = 51;
        size_t nodesBackTopLeftCorner_ = 60;
        size_t nodesBackTopRightCorner_ = 63;

        xt::xtensor<size_t, 2> dofs_ = {
            {0, 1, 2},       {3, 4, 5},       {6, 7, 8},       {9, 10, 11},     {12, 13, 14},
            {15, 16, 17},    {18, 19, 20},    {21, 22, 23},    {24, 25, 26},    {27, 28, 29},
            {30, 31, 32},    {33, 34, 35},    {36, 37, 38},    {39, 40, 41},    {42, 43, 44},
            {45, 46, 47},    {48, 49, 50},    {51, 52, 53},    {54, 55, 56},    {57, 58, 59},
            {60, 61, 62},    {63, 64, 65},    {66, 67, 68},    {69, 70, 71},    {72, 73, 74},
            {75, 76, 77},    {78, 79, 80},    {81, 82, 83},    {84, 85, 86},    {87, 88, 89},
            {90, 91, 92},    {93, 94, 95},    {96, 97, 98},    {99, 100, 101},  {102, 103, 104},
            {105, 106, 107}, {108, 109, 110}, {111, 112, 113}, {114, 115, 116}, {117, 118, 119},
            {120, 121, 122}, {123, 124, 125}, {126, 127, 128}, {129, 130, 131}, {132, 133, 134},
            {135, 136, 137}, {138, 139, 140}, {141, 142, 143}, {144, 145, 146}, {147, 148, 149},
            {150, 151, 152}, {153, 154, 155}, {156, 157, 158}, {159, 160, 161}, {162, 163, 164},
            {165, 166, 167}, {168, 169, 170}, {171, 172, 173}, {174, 175, 176}, {177, 178, 179},
            {180, 181, 182}, {183, 184, 185}, {186, 187, 188}, {189, 190, 191}};

        xt::xtensor<size_t, 2> dofsPeriodic_ = {
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {0, 1, 2},    {9, 10, 11},  {12, 13, 14},
            {15, 16, 17}, {9, 10, 11},  {18, 19, 20}, {21, 22, 23}, {24, 25, 26}, {18, 19, 20},
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {0, 1, 2},    {27, 28, 29}, {30, 31, 32},
            {33, 34, 35}, {27, 28, 29}, {36, 37, 38}, {39, 40, 41}, {42, 43, 44}, {36, 37, 38},
            {45, 46, 47}, {48, 49, 50}, {51, 52, 53}, {45, 46, 47}, {27, 28, 29}, {30, 31, 32},
            {33, 34, 35}, {27, 28, 29}, {54, 55, 56}, {57, 58, 59}, {60, 61, 62}, {54, 55, 56},
            {63, 64, 65}, {66, 67, 68}, {69, 70, 71}, {63, 64, 65}, {72, 73, 74}, {75, 76, 77},
            {78, 79, 80}, {72, 73, 74}, {54, 55, 56}, {57, 58, 59}, {60, 61, 62}, {54, 55, 56},
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {0, 1, 2},    {9, 10, 11},  {12, 13, 14},
            {15, 16, 17}, {9, 10, 11},  {18, 19, 20}, {21, 22, 23}, {24, 25, 26}, {18, 19, 20},
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {0, 1, 2}};

        xt::xtensor<size_t, 2> nodesPeriodic_ = {
            {0, 3},   {0, 51},  {0, 48},  {0, 12},  {0, 15},  {0, 63},  {0, 60},  {1, 49},
            {2, 50},  {1, 61},  {2, 62},  {1, 13},  {2, 14},  {16, 19}, {32, 35}, {16, 31},
            {32, 47}, {16, 28}, {32, 44}, {4, 7},   {8, 11},  {4, 55},  {8, 59},  {4, 52},
            {8, 56},  {5, 53},  {6, 54},  {9, 57},  {10, 58}, {20, 23}, {24, 27}, {36, 39},
            {40, 43}, {17, 29}, {18, 30}, {33, 45}, {34, 46}};

        size_t nodesOrigin_ = 0;

        GooseFEM::Mesh::Hex8::Regular mesh(3, 3, 3);

        REQUIRE(xt::allclose(mesh.coor(), coor_));
        REQUIRE(xt::all(xt::equal(mesh.conn(), conn_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottom(), nodesBottom_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTop(), nodesTop_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesLeft(), nodesLeft_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesRight(), nodesRight_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFront(), nodesFront_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBack(), nodesBack_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontTopEdge(), nodesFrontTopEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontLeftEdge(), nodesFrontLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontRightEdge(), nodesFrontRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackBottomEdge(), nodesBackBottomEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackTopEdge(), nodesBackTopEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackLeftEdge(), nodesBackLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackRightEdge(), nodesBackRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomLeftEdge(), nodesBottomLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomRightEdge(), nodesBottomRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopLeftEdge(), nodesTopLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopRightEdge(), nodesTopRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontFace(), nodesFrontFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackFace(), nodesBackFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesLeftFace(), nodesLeftFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesRightFace(), nodesRightFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomFace(), nodesBottomFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopFace(), nodesTopFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontBottomOpenEdge(), nodesFrontBottomOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontTopOpenEdge(), nodesFrontTopOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontLeftOpenEdge(), nodesFrontLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontRightOpenEdge(), nodesFrontRightOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackBottomOpenEdge(), nodesBackBottomOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackTopOpenEdge(), nodesBackTopOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackLeftOpenEdge(), nodesBackLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackRightOpenEdge(), nodesBackRightOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomLeftOpenEdge(), nodesBottomLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomRightOpenEdge(), nodesBottomRightOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopLeftOpenEdge(), nodesTopLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopRightOpenEdge(), nodesTopRightOpenEdge_)));
        REQUIRE(mesh.nodesFrontBottomLeftCorner() == nodesFrontBottomLeftCorner_);
        REQUIRE(mesh.nodesFrontBottomRightCorner() == nodesFrontBottomRightCorner_);
        REQUIRE(mesh.nodesFrontTopLeftCorner() == nodesFrontTopLeftCorner_);
        REQUIRE(mesh.nodesFrontTopRightCorner() == nodesFrontTopRightCorner_);
        REQUIRE(mesh.nodesBackBottomLeftCorner() == nodesBackBottomLeftCorner_);
        REQUIRE(mesh.nodesBackBottomRightCorner() == nodesBackBottomRightCorner_);
        REQUIRE(mesh.nodesBackTopLeftCorner() == nodesBackTopLeftCorner_);
        REQUIRE(mesh.nodesBackTopRightCorner() == nodesBackTopRightCorner_);
        REQUIRE(xt::all(xt::equal(mesh.dofs(), dofs_)));
        REQUIRE(xt::all(xt::equal(mesh.dofsPeriodic(), dofsPeriodic_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesPeriodic(), nodesPeriodic_)));
        REQUIRE(mesh.nodesOrigin() == nodesOrigin_);
    }

    SECTION("FineLayer - basic")
    {
        GooseFEM::Mesh::Hex8::FineLayer mesh(3, 3, 3);

        xt::xtensor<double, 2> coor_ = {
            {0., 0., 0.}, {1., 0., 0.}, {2., 0., 0.}, {3., 0., 0.}, {0., 0., 1.}, {1., 0., 1.},
            {2., 0., 1.}, {3., 0., 1.}, {0., 0., 2.}, {1., 0., 2.}, {2., 0., 2.}, {3., 0., 2.},
            {0., 0., 3.}, {1., 0., 3.}, {2., 0., 3.}, {3., 0., 3.}, {0., 1., 0.}, {1., 1., 0.},
            {2., 1., 0.}, {3., 1., 0.}, {0., 1., 1.}, {1., 1., 1.}, {2., 1., 1.}, {3., 1., 1.},
            {0., 1., 2.}, {1., 1., 2.}, {2., 1., 2.}, {3., 1., 2.}, {0., 1., 3.}, {1., 1., 3.},
            {2., 1., 3.}, {3., 1., 3.}, {0., 2., 0.}, {1., 2., 0.}, {2., 2., 0.}, {3., 2., 0.},
            {0., 2., 1.}, {1., 2., 1.}, {2., 2., 1.}, {3., 2., 1.}, {0., 2., 2.}, {1., 2., 2.},
            {2., 2., 2.}, {3., 2., 2.}, {0., 2., 3.}, {1., 2., 3.}, {2., 2., 3.}, {3., 2., 3.},
            {0., 3., 0.}, {1., 3., 0.}, {2., 3., 0.}, {3., 3., 0.}, {0., 3., 1.}, {1., 3., 1.},
            {2., 3., 1.}, {3., 3., 1.}, {0., 3., 2.}, {1., 3., 2.}, {2., 3., 2.}, {3., 3., 2.},
            {0., 3., 3.}, {1., 3., 3.}, {2., 3., 3.}, {3., 3., 3.}};

        xt::xtensor<double, 2> conn_ = {
            {0, 1, 17, 16, 4, 5, 21, 20},     {1, 2, 18, 17, 5, 6, 22, 21},
            {2, 3, 19, 18, 6, 7, 23, 22},     {4, 5, 21, 20, 8, 9, 25, 24},
            {5, 6, 22, 21, 9, 10, 26, 25},    {6, 7, 23, 22, 10, 11, 27, 26},
            {8, 9, 25, 24, 12, 13, 29, 28},   {9, 10, 26, 25, 13, 14, 30, 29},
            {10, 11, 27, 26, 14, 15, 31, 30}, {16, 17, 33, 32, 20, 21, 37, 36},
            {17, 18, 34, 33, 21, 22, 38, 37}, {18, 19, 35, 34, 22, 23, 39, 38},
            {20, 21, 37, 36, 24, 25, 41, 40}, {21, 22, 38, 37, 25, 26, 42, 41},
            {22, 23, 39, 38, 26, 27, 43, 42}, {24, 25, 41, 40, 28, 29, 45, 44},
            {25, 26, 42, 41, 29, 30, 46, 45}, {26, 27, 43, 42, 30, 31, 47, 46},
            {32, 33, 49, 48, 36, 37, 53, 52}, {33, 34, 50, 49, 37, 38, 54, 53},
            {34, 35, 51, 50, 38, 39, 55, 54}, {36, 37, 53, 52, 40, 41, 57, 56},
            {37, 38, 54, 53, 41, 42, 58, 57}, {38, 39, 55, 54, 42, 43, 59, 58},
            {40, 41, 57, 56, 44, 45, 61, 60}, {41, 42, 58, 57, 45, 46, 62, 61},
            {42, 43, 59, 58, 46, 47, 63, 62}};

        xt::xtensor<size_t, 1> nodesBottom_ = {
            0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15};

        xt::xtensor<size_t, 1> nodesTop_ = {
            48, 52, 56, 60, 49, 53, 57, 61, 50, 54, 58, 62, 51, 55, 59, 63};

        xt::xtensor<size_t, 1> nodesLeft_ = {
            0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60};

        xt::xtensor<size_t, 1> nodesRight_ = {
            3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63};

        xt::xtensor<size_t, 1> nodesFront_ = {
            0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51};

        xt::xtensor<size_t, 1> nodesBack_ = {
            12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63};

        xt::xtensor<size_t, 1> nodesFrontTopEdge_ = {48, 49, 50, 51};
        xt::xtensor<size_t, 1> nodesFrontLeftEdge_ = {0, 16, 32, 48};
        xt::xtensor<size_t, 1> nodesFrontRightEdge_ = {3, 19, 35, 51};
        xt::xtensor<size_t, 1> nodesBackBottomEdge_ = {12, 13, 14, 15};
        xt::xtensor<size_t, 1> nodesBackTopEdge_ = {60, 61, 62, 63};
        xt::xtensor<size_t, 1> nodesBackLeftEdge_ = {12, 28, 44, 60};
        xt::xtensor<size_t, 1> nodesBackRightEdge_ = {15, 31, 47, 63};
        xt::xtensor<size_t, 1> nodesBottomLeftEdge_ = {0, 4, 8, 12};
        xt::xtensor<size_t, 1> nodesBottomRightEdge_ = {3, 7, 11, 15};
        xt::xtensor<size_t, 1> nodesTopLeftEdge_ = {48, 52, 56, 60};
        xt::xtensor<size_t, 1> nodesTopRightEdge_ = {51, 55, 59, 63};
        xt::xtensor<size_t, 1> nodesFrontFace_ = {17, 18, 33, 34};
        xt::xtensor<size_t, 1> nodesBackFace_ = {29, 30, 45, 46};
        xt::xtensor<size_t, 1> nodesLeftFace_ = {20, 24, 36, 40};
        xt::xtensor<size_t, 1> nodesRightFace_ = {23, 27, 39, 43};
        xt::xtensor<size_t, 1> nodesBottomFace_ = {5, 9, 6, 10};
        xt::xtensor<size_t, 1> nodesTopFace_ = {53, 57, 54, 58};
        xt::xtensor<size_t, 1> nodesFrontBottomOpenEdge_ = {1, 2};
        xt::xtensor<size_t, 1> nodesFrontTopOpenEdge_ = {49, 50};
        xt::xtensor<size_t, 1> nodesFrontLeftOpenEdge_ = {16, 32};
        xt::xtensor<size_t, 1> nodesFrontRightOpenEdge_ = {19, 35};
        xt::xtensor<size_t, 1> nodesBackBottomOpenEdge_ = {13, 14};
        xt::xtensor<size_t, 1> nodesBackTopOpenEdge_ = {61, 62};
        xt::xtensor<size_t, 1> nodesBackLeftOpenEdge_ = {28, 44};
        xt::xtensor<size_t, 1> nodesBackRightOpenEdge_ = {31, 47};
        xt::xtensor<size_t, 1> nodesBottomLeftOpenEdge_ = {4, 8};
        xt::xtensor<size_t, 1> nodesBottomRightOpenEdge_ = {7, 11};
        xt::xtensor<size_t, 1> nodesTopLeftOpenEdge_ = {52, 56};
        xt::xtensor<size_t, 1> nodesTopRightOpenEdge_ = {55, 59};
        size_t nodesFrontBottomLeftCorner_ = 0;
        size_t nodesFrontBottomRightCorner_ = 3;
        size_t nodesFrontTopLeftCorner_ = 48;
        size_t nodesFrontTopRightCorner_ = 51;
        size_t nodesBackBottomLeftCorner_ = 12;
        size_t nodesBackBottomRightCorner_ = 15;
        size_t nodesBackTopLeftCorner_ = 60;
        size_t nodesBackTopRightCorner_ = 63;

        xt::xtensor<size_t, 2> dofs_ = {
            {0, 1, 2},       {3, 4, 5},       {6, 7, 8},       {9, 10, 11},     {12, 13, 14},
            {15, 16, 17},    {18, 19, 20},    {21, 22, 23},    {24, 25, 26},    {27, 28, 29},
            {30, 31, 32},    {33, 34, 35},    {36, 37, 38},    {39, 40, 41},    {42, 43, 44},
            {45, 46, 47},    {48, 49, 50},    {51, 52, 53},    {54, 55, 56},    {57, 58, 59},
            {60, 61, 62},    {63, 64, 65},    {66, 67, 68},    {69, 70, 71},    {72, 73, 74},
            {75, 76, 77},    {78, 79, 80},    {81, 82, 83},    {84, 85, 86},    {87, 88, 89},
            {90, 91, 92},    {93, 94, 95},    {96, 97, 98},    {99, 100, 101},  {102, 103, 104},
            {105, 106, 107}, {108, 109, 110}, {111, 112, 113}, {114, 115, 116}, {117, 118, 119},
            {120, 121, 122}, {123, 124, 125}, {126, 127, 128}, {129, 130, 131}, {132, 133, 134},
            {135, 136, 137}, {138, 139, 140}, {141, 142, 143}, {144, 145, 146}, {147, 148, 149},
            {150, 151, 152}, {153, 154, 155}, {156, 157, 158}, {159, 160, 161}, {162, 163, 164},
            {165, 166, 167}, {168, 169, 170}, {171, 172, 173}, {174, 175, 176}, {177, 178, 179},
            {180, 181, 182}, {183, 184, 185}, {186, 187, 188}, {189, 190, 191}};

        xt::xtensor<size_t, 2> dofsPeriodic_ = {
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {0, 1, 2},    {9, 10, 11},  {12, 13, 14},
            {15, 16, 17}, {9, 10, 11},  {18, 19, 20}, {21, 22, 23}, {24, 25, 26}, {18, 19, 20},
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {0, 1, 2},    {27, 28, 29}, {30, 31, 32},
            {33, 34, 35}, {27, 28, 29}, {36, 37, 38}, {39, 40, 41}, {42, 43, 44}, {36, 37, 38},
            {45, 46, 47}, {48, 49, 50}, {51, 52, 53}, {45, 46, 47}, {27, 28, 29}, {30, 31, 32},
            {33, 34, 35}, {27, 28, 29}, {54, 55, 56}, {57, 58, 59}, {60, 61, 62}, {54, 55, 56},
            {63, 64, 65}, {66, 67, 68}, {69, 70, 71}, {63, 64, 65}, {72, 73, 74}, {75, 76, 77},
            {78, 79, 80}, {72, 73, 74}, {54, 55, 56}, {57, 58, 59}, {60, 61, 62}, {54, 55, 56},
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {0, 1, 2},    {9, 10, 11},  {12, 13, 14},
            {15, 16, 17}, {9, 10, 11},  {18, 19, 20}, {21, 22, 23}, {24, 25, 26}, {18, 19, 20},
            {0, 1, 2},    {3, 4, 5},    {6, 7, 8},    {0, 1, 2}};

        xt::xtensor<size_t, 2> nodesPeriodic_ = {
            {0, 3},   {0, 15},  {0, 12},  {0, 48},  {0, 51},  {0, 63},  {0, 60},  {1, 13},
            {2, 14},  {1, 61},  {2, 62},  {1, 49},  {2, 50},  {4, 7},   {8, 11},  {4, 55},
            {8, 59},  {4, 52},  {8, 56},  {16, 19}, {32, 35}, {16, 31}, {32, 47}, {16, 28},
            {32, 44}, {17, 29}, {18, 30}, {33, 45}, {34, 46}, {20, 23}, {24, 27}, {36, 39},
            {40, 43}, {5, 53},  {9, 57},  {6, 54},  {10, 58}};

        size_t nodesOrigin_ = 0;

        REQUIRE(xt::allclose(mesh.coor(), coor_));
        REQUIRE(xt::all(xt::equal(mesh.conn(), conn_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottom(), nodesBottom_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTop(), nodesTop_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesLeft(), nodesLeft_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesRight(), nodesRight_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFront(), nodesFront_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBack(), nodesBack_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontTopEdge(), nodesFrontTopEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontLeftEdge(), nodesFrontLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontRightEdge(), nodesFrontRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackBottomEdge(), nodesBackBottomEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackTopEdge(), nodesBackTopEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackLeftEdge(), nodesBackLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackRightEdge(), nodesBackRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomLeftEdge(), nodesBottomLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomRightEdge(), nodesBottomRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopLeftEdge(), nodesTopLeftEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopRightEdge(), nodesTopRightEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontFace(), nodesFrontFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackFace(), nodesBackFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesLeftFace(), nodesLeftFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesRightFace(), nodesRightFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomFace(), nodesBottomFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopFace(), nodesTopFace_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontBottomOpenEdge(), nodesFrontBottomOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontTopOpenEdge(), nodesFrontTopOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontLeftOpenEdge(), nodesFrontLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesFrontRightOpenEdge(), nodesFrontRightOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackBottomOpenEdge(), nodesBackBottomOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackTopOpenEdge(), nodesBackTopOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackLeftOpenEdge(), nodesBackLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBackRightOpenEdge(), nodesBackRightOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomLeftOpenEdge(), nodesBottomLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesBottomRightOpenEdge(), nodesBottomRightOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopLeftOpenEdge(), nodesTopLeftOpenEdge_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesTopRightOpenEdge(), nodesTopRightOpenEdge_)));
        REQUIRE(mesh.nodesFrontBottomLeftCorner() == nodesFrontBottomLeftCorner_);
        REQUIRE(mesh.nodesFrontBottomRightCorner() == nodesFrontBottomRightCorner_);
        REQUIRE(mesh.nodesFrontTopLeftCorner() == nodesFrontTopLeftCorner_);
        REQUIRE(mesh.nodesFrontTopRightCorner() == nodesFrontTopRightCorner_);
        REQUIRE(mesh.nodesBackBottomLeftCorner() == nodesBackBottomLeftCorner_);
        REQUIRE(mesh.nodesBackBottomRightCorner() == nodesBackBottomRightCorner_);
        REQUIRE(mesh.nodesBackTopLeftCorner() == nodesBackTopLeftCorner_);
        REQUIRE(mesh.nodesBackTopRightCorner() == nodesBackTopRightCorner_);
        REQUIRE(xt::all(xt::equal(mesh.dofs(), dofs_)));
        REQUIRE(xt::all(xt::equal(mesh.dofsPeriodic(), dofsPeriodic_)));
        REQUIRE(xt::all(xt::equal(mesh.nodesPeriodic(), nodesPeriodic_)));
        REQUIRE(mesh.nodesOrigin() == nodesOrigin_);
    }
}

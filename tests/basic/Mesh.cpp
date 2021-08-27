
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Mesh", "Mesh.h")
{
    SECTION("overlapping")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(5, 5, 1.0);

        auto coor_a = mesh.coor();
        auto overlap_a = mesh.nodesTopEdge();

        auto coor_b = mesh.coor();
        auto overlap_b = mesh.nodesBottomEdge();
        xt::view(coor_b, xt::all(), 1) += 5.0;

        auto overlap = GooseFEM::Mesh::overlapping(coor_a, coor_b);

        REQUIRE(xt::all(xt::equal(xt::view(overlap, 0, xt::all()), overlap_a)));
        REQUIRE(xt::all(xt::equal(xt::view(overlap, 1, xt::all()), overlap_b)));
    }

    SECTION("ManualStitch")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(5, 5, 1.0);

        auto coor_a = mesh.coor();
        auto conn_a = mesh.conn();
        auto overlap_a = mesh.nodesTopEdge();

        auto coor_b = mesh.coor();
        auto conn_b = mesh.conn();
        auto overlap_b = mesh.nodesBottomEdge();
        xt::view(coor_b, xt::all(), 1) += 5.0;

        GooseFEM::Mesh::ManualStitch stitch(coor_a, conn_a, overlap_a, coor_b, conn_b, overlap_b);

        GooseFEM::Mesh::Quad4::Regular res(5, 10, 1.0);

        REQUIRE(xt::allclose(stitch.coor(), res.coor()));
        REQUIRE(xt::all(xt::equal(stitch.conn(), res.conn())));

        REQUIRE(stitch.nodemap().size() == 2);
        REQUIRE(stitch.elemmap().size() == 2);
        REQUIRE(xt::all(xt::equal(stitch.nodemap(0), xt::arange<size_t>(6 * 6))));
        REQUIRE(xt::all(xt::equal(stitch.nodemap(1), xt::arange<size_t>(6 * 6) + 5 * 6)));
        REQUIRE(xt::all(xt::equal(stitch.elemmap(0), xt::arange<size_t>(5 * 5))));
        REQUIRE(xt::all(xt::equal(stitch.elemmap(1), xt::arange<size_t>(5 * 5) + 5 * 5)));

        REQUIRE(xt::all(xt::equal(stitch.nodeset(xt::eval(xt::arange<size_t>(6 * 6)), 0), xt::arange<size_t>(6 * 6))));
        REQUIRE(xt::all(xt::equal(stitch.nodeset(xt::eval(xt::arange<size_t>(6 * 6)), 1), xt::arange<size_t>(6 * 6) + 5 * 6)));
        REQUIRE(xt::all(xt::equal(stitch.elemset(xt::eval(xt::arange<size_t>(5 * 5)), 0), xt::arange<size_t>(5 * 5))));
        REQUIRE(xt::all(xt::equal(stitch.elemset(xt::eval(xt::arange<size_t>(5 * 5)), 1), xt::arange<size_t>(5 * 5) + 5 * 5)));

        REQUIRE(xt::all(xt::equal(stitch.nodeset(overlap_a, 0), stitch.nodeset(overlap_b, 1))));
    }

    SECTION("Stitch")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(5, 5, 1.0);

        auto coor_a = mesh.coor();
        auto conn_a = mesh.conn();
        auto overlap_a = mesh.nodesTopEdge();
        auto nset = mesh.nodesLeftEdge();
        auto eset = xt::eval(xt::arange<size_t>(mesh.nelem()));

        auto coor_b = mesh.coor();
        auto conn_b = mesh.conn();
        auto overlap_b = mesh.nodesBottomEdge();
        xt::view(coor_b, xt::all(), 1) += 5.0;

        GooseFEM::Mesh::Stitch stitch;
        stitch.push_back(coor_a, conn_a);
        stitch.push_back(coor_b, conn_b);

        GooseFEM::Mesh::Quad4::Regular res(5, 10, 1.0);

        REQUIRE(xt::allclose(stitch.coor(), res.coor()));
        REQUIRE(xt::all(xt::equal(stitch.conn(), res.conn())));

        REQUIRE(stitch.nodemap().size() == 2);
        REQUIRE(stitch.elemmap().size() == 2);
        REQUIRE(xt::all(xt::equal(stitch.nodemap(0), xt::arange<size_t>(6 * 6))));
        REQUIRE(xt::all(xt::equal(stitch.nodemap(1), xt::arange<size_t>(6 * 6) + 5 * 6)));
        REQUIRE(xt::all(xt::equal(stitch.elemmap(0), xt::arange<size_t>(5 * 5))));
        REQUIRE(xt::all(xt::equal(stitch.elemmap(1), xt::arange<size_t>(5 * 5) + 5 * 5)));

        REQUIRE(xt::all(xt::equal(stitch.nodeset(xt::eval(xt::arange<size_t>(6 * 6)), 0), xt::arange<size_t>(6 * 6))));
        REQUIRE(xt::all(xt::equal(stitch.nodeset(xt::eval(xt::arange<size_t>(6 * 6)), 1), xt::arange<size_t>(6 * 6) + 5 * 6)));
        REQUIRE(xt::all(xt::equal(stitch.elemset(xt::eval(xt::arange<size_t>(5 * 5)), 0), xt::arange<size_t>(5 * 5))));
        REQUIRE(xt::all(xt::equal(stitch.elemset(xt::eval(xt::arange<size_t>(5 * 5)), 1), xt::arange<size_t>(5 * 5) + 5 * 5)));

        REQUIRE(xt::all(xt::equal(stitch.nodeset(overlap_a, 0), stitch.nodeset(overlap_b, 1))));

        REQUIRE(xt::all(xt::equal(stitch.nodeset({nset, nset}), xt::arange<size_t>(0, 6 * 6 + 5 * 6, 6))));
        REQUIRE(xt::all(xt::equal(stitch.elemset({eset, eset}), xt::arange<size_t>(2 * 5 * 5))));
    }

    SECTION("Vstack - equivalent test as 'Stitch'")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(5, 5, 1.0);

        GooseFEM::Mesh::Vstack stitch;
        stitch.push_back(mesh.coor(), mesh.conn(), mesh.nodesBottomEdge(), mesh.nodesTopEdge());
        stitch.push_back(mesh.coor(), mesh.conn(), mesh.nodesBottomEdge(), mesh.nodesTopEdge());

        auto nset = mesh.nodesLeftEdge();
        auto eset = xt::eval(xt::arange<size_t>(mesh.nelem()));

        GooseFEM::Mesh::Quad4::Regular res(5, 10, 1.0);

        REQUIRE(xt::allclose(stitch.coor(), res.coor()));
        REQUIRE(xt::all(xt::equal(stitch.conn(), res.conn())));

        REQUIRE(stitch.nodemap().size() == 2);
        REQUIRE(stitch.elemmap().size() == 2);
        REQUIRE(xt::all(xt::equal(stitch.nodemap(0), xt::arange<size_t>(6 * 6))));
        REQUIRE(xt::all(xt::equal(stitch.nodemap(1), xt::arange<size_t>(6 * 6) + 5 * 6)));
        REQUIRE(xt::all(xt::equal(stitch.elemmap(0), xt::arange<size_t>(5 * 5))));
        REQUIRE(xt::all(xt::equal(stitch.elemmap(1), xt::arange<size_t>(5 * 5) + 5 * 5)));

        REQUIRE(xt::all(xt::equal(stitch.nodeset(xt::eval(xt::arange<size_t>(6 * 6)), 0), xt::arange<size_t>(6 * 6))));
        REQUIRE(xt::all(xt::equal(stitch.nodeset(xt::eval(xt::arange<size_t>(6 * 6)), 1), xt::arange<size_t>(6 * 6) + 5 * 6)));
        REQUIRE(xt::all(xt::equal(stitch.elemset(xt::eval(xt::arange<size_t>(5 * 5)), 0), xt::arange<size_t>(5 * 5))));
        REQUIRE(xt::all(xt::equal(stitch.elemset(xt::eval(xt::arange<size_t>(5 * 5)), 1), xt::arange<size_t>(5 * 5) + 5 * 5)));

        REQUIRE(xt::all(xt::equal(stitch.nodeset({nset, nset}), xt::arange<size_t>(0, 6 * 6 + 5 * 6, 6))));
        REQUIRE(xt::all(xt::equal(stitch.elemset({eset, eset}), xt::arange<size_t>(2 * 5 * 5))));
    }

    SECTION("Vstack - several layers")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(5, 5, 1.0);

        GooseFEM::Mesh::Vstack stitch;
        stitch.push_back(mesh.coor(), mesh.conn(), mesh.nodesBottomEdge(), mesh.nodesTopEdge());
        stitch.push_back(mesh.coor(), mesh.conn(), mesh.nodesBottomEdge(), mesh.nodesTopEdge());
        stitch.push_back(mesh.coor(), mesh.conn(), mesh.nodesBottomEdge(), mesh.nodesTopEdge());
        stitch.push_back(mesh.coor(), mesh.conn(), mesh.nodesBottomEdge(), mesh.nodesTopEdge());
        stitch.push_back(mesh.coor(), mesh.conn(), mesh.nodesBottomEdge(), mesh.nodesTopEdge());

        GooseFEM::Mesh::Quad4::Regular res(5, 5 * 5, 1.0);

        REQUIRE(xt::allclose(stitch.coor(), res.coor()));
        REQUIRE(xt::all(xt::equal(stitch.conn(), res.conn())));

        REQUIRE(stitch.nodemap().size() == 5);
        REQUIRE(stitch.elemmap().size() == 5);
    }

    SECTION("edgesize")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2, 10.0);
        auto s = GooseFEM::Mesh::edgesize(mesh.coor(), mesh.conn());
        auto t = GooseFEM::Mesh::edgesize(mesh.coor(), mesh.conn(), mesh.getElementType());
        REQUIRE(xt::allclose(s, 10.0));
        REQUIRE(xt::allclose(t, 10.0));
    }

    SECTION("coordination")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        auto N = GooseFEM::Mesh::coordination(mesh.conn());
        xt::xtensor<size_t, 1> ret = {1, 2, 2, 1, 2, 4, 4, 2, 2, 4, 4, 2, 1, 2, 2, 1};
        REQUIRE(xt::all(xt::equal(N, ret)));
    }

    SECTION("centers")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2, 2.0);
        xt::xtensor<double, 2> c = {
            {1.0, 1.0},
            {3.0, 1.0},
            {1.0, 3.0},
            {3.0, 3.0}};

        REQUIRE(xt::allclose(GooseFEM::Mesh::centers(mesh.coor(), mesh.conn()), c));
    }

    SECTION("elem2node")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        auto tonode = GooseFEM::Mesh::elem2node(mesh.conn());
        REQUIRE(tonode.size() == 16);
        REQUIRE(tonode[0] == std::vector<size_t>{0});
        REQUIRE(tonode[1] == std::vector<size_t>{0, 1});
        REQUIRE(tonode[2] == std::vector<size_t>{1, 2});
        REQUIRE(tonode[3] == std::vector<size_t>{2});
        REQUIRE(tonode[4] == std::vector<size_t>{0, 3});
        REQUIRE(tonode[5] == std::vector<size_t>{0, 1, 3, 4});
        REQUIRE(tonode[6] == std::vector<size_t>{1, 2, 4, 5});
        REQUIRE(tonode[7] == std::vector<size_t>{2, 5});
        REQUIRE(tonode[8] == std::vector<size_t>{3, 6});
        REQUIRE(tonode[9] == std::vector<size_t>{3, 4, 6, 7});
        REQUIRE(tonode[10] == std::vector<size_t>{4, 5, 7, 8});
        REQUIRE(tonode[11] == std::vector<size_t>{5, 8});
        REQUIRE(tonode[12] == std::vector<size_t>{6});
        REQUIRE(tonode[13] == std::vector<size_t>{6, 7});
        REQUIRE(tonode[14] == std::vector<size_t>{7, 8});
        REQUIRE(tonode[15] == std::vector<size_t>{8});
    }

    SECTION("elemmap2nodemap")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);

        xt::xtensor<size_t, 1> elmap0 = {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8
        };
        xt::xtensor<size_t, 1> elmap1 = {
            2, 0, 1,
            5, 3, 4,
            8, 6, 7
        };
        xt::xtensor<size_t, 1> elmap2 = {
            1, 2, 0,
            4, 5, 3,
            7, 8, 6
        };

        xt::xtensor<size_t, 1> nodemap0 = {
             0,  1,  2,  3,
             4,  5,  6,  7,
             8,  9, 10, 11,
            12, 13, 14, 15
        };
        xt::xtensor<size_t, 1> nodemap1 = {
              2,  0,  1,  2,
              6,  4,  5,  6,
             10,  8,  9, 10,
             14, 15, 13, 14
        };
        xt::xtensor<size_t, 1> nodemap2 = {
             1,  2,  0,  1,
             5,  6,  4,  5,
             9, 10,  8,  9,
            13, 14, 15, 13
        };

        REQUIRE(xt::all(xt::equal(GooseFEM::Mesh::elemmap2nodemap(elmap0, mesh.coor(), mesh.conn()), nodemap0)));
        REQUIRE(xt::all(xt::equal(GooseFEM::Mesh::elemmap2nodemap(elmap1, mesh.coor(), mesh.conn()), nodemap1)));
        REQUIRE(xt::all(xt::equal(GooseFEM::Mesh::elemmap2nodemap(elmap2, mesh.coor(), mesh.conn()), nodemap2)));
    }

    SECTION("elemmap2nodemap - example 1")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(3, 3);

        xt::xtensor<int, 1> elemval = {
            1, 0, 0,
            1, 0, 0,
            1, 0, 0
        };

        xt::xtensor<int, 1> elemval_r1 = {
            0, 1, 0,
            0, 1, 0,
            0, 1, 0
        };

        xt::xtensor<int, 1> elemval_r2 = {
            0, 0, 1,
            0, 0, 1,
            0, 0, 1
        };

        xt::xtensor<int, 1> nodeval = {
            1, 0, 0, 1,
            1, 0, 0, 1,
            1, 0, 0, 1,
            1, 0, 0, 1
        };

        xt::xtensor<int, 1> nodeval_r1 = {
            0, 1, 0, 0,
            0, 1, 0, 0,
            0, 1, 0, 0,
            0, 1, 0, 0
        };

        xt::xtensor<int, 1> nodeval_r2 = {
            0, 0, 1, 0,
            0, 0, 1, 0,
            0, 0, 1, 0,
            0, 0, 1, 0
        };

        {
            auto elemmap = mesh.roll(0);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(1);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval_r1, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval_r1, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(2);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval_r2, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval_r2, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(3);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(4);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval_r1, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval_r1, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(5);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval_r2, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval_r2, xt::view(nodeval, xt::keep(nodemap)))));
        }
    }

    SECTION("elemmap2nodemap - example 2")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(3, 3);

        xt::xtensor<int, 1> elemval = {
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        };

        xt::xtensor<int, 1> elemval_r1 = {
            0, 1, 0,
            0, 0, 1,
            1, 0, 0
        };

        xt::xtensor<int, 1> elemval_r2 = {
            0, 0, 1,
            1, 0, 0,
            0, 1, 0
        };

        xt::xtensor<int, 1> nodeval = {
            1, 0, 0, 1,
            0, 1, 0, 0,
            0, 0, 1, 0,
            1, 0, 0, 1
        };

        xt::xtensor<int, 1> nodeval_r1 = {
            0, 1, 0, 0,
            0, 0, 1, 0,
            1, 0, 0, 1,
            0, 1, 0, 0
        };

        xt::xtensor<int, 1> nodeval_r2 = {
            0, 0, 1, 0,
            1, 0, 0, 1,
            0, 1, 0, 0,
            0, 0, 1, 0
        };

        {
            auto elemmap = mesh.roll(0);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(1);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval_r1, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval_r1, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(2);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval_r2, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval_r2, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(3);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(4);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval_r1, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval_r1, xt::view(nodeval, xt::keep(nodemap)))));
        }

        {
            auto elemmap = mesh.roll(5);
            auto nodemap = GooseFEM::Mesh::elemmap2nodemap(elemmap, mesh.coor(), mesh.conn());
            REQUIRE(xt::all(xt::equal(elemval_r2, xt::view(elemval, xt::keep(elemmap)))));
            REQUIRE(xt::all(xt::equal(nodeval_r2, xt::view(nodeval, xt::keep(nodemap)))));
        }
    }
}

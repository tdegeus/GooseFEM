#define CATCH_CONFIG_MAIN // tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::VectorPartitioned", "VectorPartitioned.h")
{

    SECTION("basic")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        auto iip = xt::eval(xt::arange<size_t>(4 * 2));
        GooseFEM::VectorPartitioned vector(mesh.conn(), mesh.dofs(), iip);

        xt::xtensor<bool, 2> is_p = {
            {true, true},
            {true, true},
            {true, true},
            {true, true},
            {false, false},
            {false, false},
            {false, false},
            {false, false},
            {false, false},
            {false, false},
            {false, false},
            {false, false},
            {false, false},
            {false, false},
            {false, false},
            {false, false}};

        REQUIRE(xt::all(xt::equal(vector.dofs_is_p(), is_p)));
        REQUIRE(xt::all(xt::not_equal(vector.dofs_is_u(), is_p)));
    }

    SECTION("asDofs")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(9, 9);

        auto dofs = mesh.dofs();
        auto nodesLeft = mesh.nodesLeftOpenEdge();
        auto nodesRight = mesh.nodesRightOpenEdge();
        auto nodesTop = mesh.nodesTopEdge();
        auto nodesBottom = mesh.nodesBottomEdge();

        xt::view(dofs, xt::keep(nodesRight)) = xt::view(dofs, xt::keep(nodesLeft));

        size_t ni = nodesBottom.size();
        xt::xtensor<size_t, 1> iip = xt::empty<size_t>({4 * ni});
        xt::view(iip, xt::range(0 * ni, 1 * ni)) = xt::view(dofs, xt::keep(nodesTop), 0);
        xt::view(iip, xt::range(1 * ni, 2 * ni)) = xt::view(dofs, xt::keep(nodesTop), 1);
        xt::view(iip, xt::range(2 * ni, 3 * ni)) = xt::view(dofs, xt::keep(nodesBottom), 0);
        xt::view(iip, xt::range(3 * ni, 4 * ni)) = xt::view(dofs, xt::keep(nodesBottom), 1);

        xt::xtensor<double, 2> u = xt::random::randn<double>({mesh.nnode(), mesh.ndim()});
        xt::view(u, xt::keep(nodesRight)) = xt::view(u, xt::keep(nodesLeft));

        GooseFEM::VectorPartitioned vector(mesh.conn(), dofs, iip);
        auto u_u = vector.AsDofs_u(u);
        auto u_p = vector.AsDofs_p(u);
        REQUIRE(xt::allclose(u, vector.NodeFromPartitioned(u_u, u_p)));
    }

    SECTION("copy_u, copy_p")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(9, 9);

        auto dofs = mesh.dofs();
        auto nodesLeft = mesh.nodesLeftOpenEdge();
        auto nodesRight = mesh.nodesRightOpenEdge();
        auto nodesTop = mesh.nodesTopEdge();
        auto nodesBottom = mesh.nodesBottomEdge();

        xt::view(dofs, xt::keep(nodesRight)) = xt::view(dofs, xt::keep(nodesLeft));

        size_t ni = nodesBottom.size();
        xt::xtensor<size_t, 1> iip = xt::empty<size_t>({4 * ni});
        xt::view(iip, xt::range(0 * ni, 1 * ni)) = xt::view(dofs, xt::keep(nodesTop), 0);
        xt::view(iip, xt::range(1 * ni, 2 * ni)) = xt::view(dofs, xt::keep(nodesTop), 1);
        xt::view(iip, xt::range(2 * ni, 3 * ni)) = xt::view(dofs, xt::keep(nodesBottom), 0);
        xt::view(iip, xt::range(3 * ni, 4 * ni)) = xt::view(dofs, xt::keep(nodesBottom), 1);

        xt::xtensor<double, 2> u = xt::random::randn<double>({mesh.nnode(), mesh.ndim()});
        xt::view(u, xt::keep(nodesRight)) = xt::view(u, xt::keep(nodesLeft));

        GooseFEM::VectorPartitioned vector(mesh.conn(), dofs, iip);
        xt::xtensor<double, 2> v = xt::empty<double>({mesh.nnode(), mesh.ndim()});
        vector.copy_u(u, v);
        vector.copy_p(u, v);
        REQUIRE(xt::allclose(u, v));
    }
}

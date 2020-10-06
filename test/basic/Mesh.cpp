
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Mesh", "Mesh.h")
{

    SECTION("edgesize")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2, 10.0);
        auto s = GooseFEM::Mesh::edgesize(mesh.coor(), mesh.conn());
        auto t = GooseFEM::Mesh::edgesize(mesh.coor(), mesh.conn(), mesh.getElementType());
        REQUIRE(xt::allclose(s, 10.0));
        REQUIRE(xt::allclose(t, 10.0));
    }
}

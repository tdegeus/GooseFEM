
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::MeshQuad4::Map", "MeshQuad4.h")
{
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

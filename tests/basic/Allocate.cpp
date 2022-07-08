#include <GooseFEM/GooseFEM.h>
#include <catch2/catch_all.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xrandom.hpp>

#define ISCLOSE(a, b) REQUIRE_THAT((a), Catch::Matchers::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Allocate", "Allocate.h")
{

    SECTION("asTensor - pre-allocated output")
    {
        xt::xtensor<int, 2> a0 = {{1, 2}, {3, 4}};

        xt::xtensor<int, 3> a1 = {{{1, 1, 1}, {2, 2, 2}}, {{3, 3, 3}, {4, 4, 4}}};

        xt::xtensor<int, 4> a2 = {
            {{{1, 1}, {1, 1}, {1, 1}}, {{2, 2}, {2, 2}, {2, 2}}},
            {{{3, 3}, {3, 3}, {3, 3}}, {{4, 4}, {4, 4}, {4, 4}}}};

        xt::xtensor<int, 3> b1 = xt::empty<int>({2, 2, 3});
        xt::xtensor<int, 4> b2 = xt::empty<int>({2, 2, 3, 2});

        GooseFEM::asTensor(a0, b1);
        GooseFEM::asTensor(a0, b2);

        REQUIRE(xt::all(xt::equal(a1, b1)));
        REQUIRE(xt::all(xt::equal(a2, b2)));
    }

    SECTION("asTensor - return")
    {
        xt::xtensor<int, 2> a0 = {{1, 2}, {3, 4}};

        xt::xtensor<int, 3> a1 = {{{1, 1, 1}, {2, 2, 2}}, {{3, 3, 3}, {4, 4, 4}}};

        xt::xtensor<int, 4> a2 = {
            {{{1, 1}, {1, 1}, {1, 1}}, {{2, 2}, {2, 2}, {2, 2}}},
            {{{3, 3}, {3, 3}, {3, 3}}, {{4, 4}, {4, 4}, {4, 4}}}};

        auto b1 = GooseFEM::AsTensor<1>(a0, 3);
        auto c1 = GooseFEM::AsTensor(1, a0, 3);
        auto b2 = GooseFEM::AsTensor(a0, {3, 2});

        static_assert(std::is_same<decltype(b1), decltype(a1)>::value, "X");
        static_assert(std::is_same<decltype(b2), decltype(a2)>::value, "X");
        static_assert(std::is_same<decltype(c1), xt::xarray<int>>::value, "X");

        REQUIRE(xt::all(xt::equal(a1, b1)));
        REQUIRE(xt::all(xt::equal(a1, c1)));
        REQUIRE(xt::all(xt::equal(a2, b2)));
    }

    SECTION("as3d")
    {
        xt::xtensor<int, 2> a1 = {{1}, {3}};

        xt::xtensor<int, 2> b1 = {{1, 0, 0}, {3, 0, 0}};

        xt::xtensor<int, 2> a2 = {{1, 2}, {3, 4}};

        xt::xtensor<int, 2> b2 = {{1, 2, 0}, {3, 4, 0}};

        xt::xtensor<int, 2> a3 = {{1, 2, 5}, {3, 4, 6}};

        xt::xtensor<int, 2> b3 = {{1, 2, 5}, {3, 4, 6}};

        REQUIRE(xt::all(xt::equal(GooseFEM::as3d(a1), b1)));
        REQUIRE(xt::all(xt::equal(GooseFEM::as3d(a2), b2)));
        REQUIRE(xt::all(xt::equal(GooseFEM::as3d(a3), b3)));
    }
}

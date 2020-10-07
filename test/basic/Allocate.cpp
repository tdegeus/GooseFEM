
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Allocate", "Allocate.h")
{

    SECTION("asTensor - 1")
    {
        xt::xtensor<double, 2> a0 = {
            {1.0, 2.0},
            {3.0, 4.0}};

        xt::xtensor<double, 3> a1 = {
            {{1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}},
            {{3.0, 3.0, 3.0}, {4.0, 4.0, 4.0}}};

        xt::xtensor<double, 4> a2 = {
            {{{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}}, {{2.0, 2.0}, {2.0, 2.0}, {2.0, 2.0}}},
            {{{3.0, 3.0}, {3.0, 3.0}, {3.0, 3.0}}, {{4.0, 4.0}, {4.0, 4.0}, {4.0, 4.0}}}};

        xt::xtensor<double, 3> b1 = xt::empty<double>({2, 2, 3});
        xt::xtensor<double, 4> b2 = xt::empty<double>({2, 2, 3, 2});

        GooseFEM::asTensor<2, 1>(a0, b1);
        GooseFEM::asTensor<2, 2>(a0, b2);

        REQUIRE(xt::allclose(a1, b1));
        REQUIRE(xt::allclose(a2, b2));
    }

    SECTION("asTensor - 2")
    {
        xt::xtensor<double, 2> a0 = {
            {1.0, 2.0},
            {3.0, 4.0}};

        xt::xtensor<double, 3> a1 = {
            {{1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}},
            {{3.0, 3.0, 3.0}, {4.0, 4.0, 4.0}}};

        xt::xtensor<double, 4> a2 = {
            {{{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}}, {{2.0, 2.0}, {2.0, 2.0}, {2.0, 2.0}}},
            {{{3.0, 3.0}, {3.0, 3.0}, {3.0, 3.0}}, {{4.0, 4.0}, {4.0, 4.0}, {4.0, 4.0}}}};

        auto b1 = GooseFEM::AsTensor<2, 1>(a0, 3);
        auto b2 = GooseFEM::AsTensor<2, 2>(a0, {3, 2});

        REQUIRE(xt::allclose(a1, b1));
        REQUIRE(xt::allclose(a2, b2));
    }

    SECTION("asTensor - 3")
    {
        xt::xtensor<double, 2> a0 = {
            {1.0, 2.0},
            {3.0, 4.0}};

        xt::xtensor<double, 3> a1 = {
            {{1.0, 1.0}, {2.0, 2.0}},
            {{3.0, 3.0}, {4.0, 4.0}}};

        xt::xtensor<double, 4> a2 = {
            {{{1.0, 1.0}, {1.0, 1.0}}, {{2.0, 2.0}, {2.0, 2.0}}},
            {{{3.0, 3.0}, {3.0, 3.0}}, {{4.0, 4.0}, {4.0, 4.0}}}};

        auto b1 = GooseFEM::AsTensor<2, 1>(a0, 2);
        auto b2 = GooseFEM::AsTensor<2, 2>(a0, 2);

        REQUIRE(xt::allclose(a1, b1));
        REQUIRE(xt::allclose(a2, b2));
    }

    SECTION("asTensor - 4")
    {
        xt::xtensor<double, 2> a0 = {
            {1.0, 2.0},
            {3.0, 4.0}};

        xt::xtensor<double, 3> a1 = {
            {{1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}},
            {{3.0, 3.0, 3.0}, {4.0, 4.0, 4.0}}};

        xt::xtensor<double, 4> a2 = {
            {{{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}}, {{2.0, 2.0}, {2.0, 2.0}, {2.0, 2.0}}},
            {{{3.0, 3.0}, {3.0, 3.0}, {3.0, 3.0}}, {{4.0, 4.0}, {4.0, 4.0}, {4.0, 4.0}}}};

        auto b1 = GooseFEM::AsTensor(1, a0, 3);
        auto b2 = GooseFEM::AsTensor(2, a0, {3, 2});

        REQUIRE(xt::allclose(a1, b1));
        REQUIRE(xt::allclose(a2, b2));
    }

    SECTION("asTensor - 5")
    {
        xt::xtensor<double, 2> a0 = {
            {1.0, 2.0},
            {3.0, 4.0}};

        xt::xtensor<double, 3> a1 = {
            {{1.0, 1.0}, {2.0, 2.0}},
            {{3.0, 3.0}, {4.0, 4.0}}};

        xt::xtensor<double, 4> a2 = {
            {{{1.0, 1.0}, {1.0, 1.0}}, {{2.0, 2.0}, {2.0, 2.0}}},
            {{{3.0, 3.0}, {3.0, 3.0}}, {{4.0, 4.0}, {4.0, 4.0}}}};

        auto b1 = GooseFEM::AsTensor(1, a0, 2);
        auto b2 = GooseFEM::AsTensor(2, a0, 2);

        REQUIRE(xt::allclose(a1, b1));
        REQUIRE(xt::allclose(a2, b2));
    }

}

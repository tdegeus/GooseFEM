
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Iterate", "Iterate.h")
{

    SECTION("StopList")
    {
        GooseFEM::Iterate::StopList stop(5);

        REQUIRE(stop.stop(5.e+0, 1.e-3) == false);
        REQUIRE(stop.stop(5.e+1, 1.e-3) == false);
        REQUIRE(stop.stop(5.e-1, 1.e-3) == false);
        REQUIRE(stop.stop(5.e-2, 1.e-3) == false);
        REQUIRE(stop.stop(5.e-3, 1.e-3) == false);
        REQUIRE(stop.stop(5.e-4, 1.e-3) == false);
        REQUIRE(stop.stop(5.e-4, 1.e-3) == false);
        REQUIRE(stop.stop(5.e-4, 1.e-3) == false);
        REQUIRE(stop.stop(5.e-4, 1.e-3) == false);
        REQUIRE(stop.stop(5.e-4, 1.e-3) == true);
    }
}

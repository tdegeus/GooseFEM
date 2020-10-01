
#include "support.h"

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

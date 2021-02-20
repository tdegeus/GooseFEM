
#include <catch2/catch.hpp>
#include <GooseFEM/GooseFEM.h>

TEST_CASE("GooseFEM::version", "version.h")
{

    SECTION("basic")
    {
        auto v = GooseFEM::version();
        auto d = GooseFEM::version_dependencies();
    }

}

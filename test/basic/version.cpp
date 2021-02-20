
#include <catch2/catch.hpp>
#include <GooseFEM/GooseFEM.h>

TEST_CASE("GooseFEM::Version", "Version.h")
{

    SECTION("basic")
    {
        auto g = GooseFEM::git();
        auto v = GooseFEM::version();
        auto d = GooseFEM::version_dependencies();
    }

}

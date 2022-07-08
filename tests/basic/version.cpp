#include <GooseFEM/GooseFEM.h>
#include <catch2/catch_all.hpp>
#include <iostream>

TEST_CASE("GooseFEM::version", "version.h")
{

    SECTION("basic")
    {
        std::cout << GooseFEM::version() << std::endl;

        auto deps = GooseFEM::version_dependencies();

        for (auto& i : deps) {
            std::cout << i << std::endl;
        }
    }
}

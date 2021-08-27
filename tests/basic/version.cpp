#define CATCH_CONFIG_MAIN // tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>
#include <GooseFEM/GooseFEM.h>
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

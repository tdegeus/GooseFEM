
#include <catch2/catch.hpp>
#include <GooseFEM/GooseFEM.h>
#include <iostream>

TEST_CASE("GooseFEM::version", "version.h")
{

    SECTION("basic")
    {
        std::cout << GooseFEM::version() << std::endl;
        std::cout << GooseFEM::version_dependencies() << std::endl;
    }

}

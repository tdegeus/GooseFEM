
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Iterate", "Iterate.h")
{
    SECTION("StopList::stop_simple - sorted input")
    {
        GooseFEM::Iterate::StopList stop(5);

        std::vector<double> res = {5e+0,  5e+1,  5e-1,  5e-2,  5e-3,  5e-4,  4e-4,  3e-4,  2e-4,  1e-4};
        std::vector<bool> conv =  {false, false, false, false, false, false, false, false, false, true};

        for (size_t i = 0; i < res.size(); ++i) {
            REQUIRE(stop.stop_simple(res[i], 1e-3) == conv[i]);
        }

        std::vector<double> m_res = { 0.05  ,  0.005 ,  0.0005,  0.0004,  0.0005};
    }

    SECTION("StopList::stop_simple - unsorted input")
    {
        GooseFEM::Iterate::StopList a(5);
        GooseFEM::Iterate::StopList b(5);

        std::vector<double> res = {5e+0,  5e+1,  5e-1,  5e-2,  5e-3,  5e-4,  4e-4,  5e-4,  4e-4,  4e-4};
        std::vector<bool> conv =  {false, false, false, false, false, false, false, false, false, true};

        for (size_t i = 0; i < res.size(); ++i) {
            REQUIRE(a.stop_simple(res[i], 1e-3) == conv[i]);
            REQUIRE(b.stop(res[i], 1e-3) == false);
        }
    }

    SECTION("StopList::stop - sorted input")
    {
        GooseFEM::Iterate::StopList stop(5);

        std::vector<double> res = {5e+0,  5e+1,  5e-1,  5e-2,  5e-3,  5e-4,  4e-4,  3e-4,  2e-4,  1e-4};
        std::vector<bool> conv =  {false, false, false, false, false, false, false, false, false, true};

        for (size_t i = 0; i < res.size(); ++i) {
            REQUIRE(stop.stop(res[i], 1e-3) == conv[i]);
        }
    }

    SECTION("StopList::stop - unsorted input")
    {
        GooseFEM::Iterate::StopList stop(5);

        std::vector<double> res = {5e+0,  5e+1,  5e-1,  5e-2,  5e-3,  5e-4,  4e-4,  3e-4,  2e-4,  3e-4,  2e-4,  1e-4,  9e-5,  8e-5};
        std::vector<bool> conv =  {false, false, false, false, false, false, false, false, false, false, false, false, false, true};

        for (size_t i = 0; i < res.size(); ++i) {
            REQUIRE(stop.stop(res[i], 1e-3) == conv[i]);
        }
    }
}

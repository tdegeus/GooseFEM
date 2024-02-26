#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>
#include <catch2/catch_all.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xrandom.hpp>

#define ISCLOSE(a, b) REQUIRE_THAT((a), Catch::Matchers::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Iterate", "Iterate.h")
{
    SECTION("StopList::all_less - sorted input")
    {
        GooseFEM::Iterate::StopList residuals(5);

        std::vector<double> res = {5e+0, 5e+1, 5e-1, 5e-2, 5e-3, 5e-4, 4e-4, 3e-4, 2e-4, 1e-4};
        std::vector<bool> conv = {
            false, false, false, false, false, false, false, false, false, true
        };

        for (size_t i = 0; i < res.size(); ++i) {
            residuals.roll_insert(res[i]);
            REQUIRE(residuals.all_less(1e-3) == conv[i]);
        }
    }

    SECTION("StopList::all_less - unsorted input")
    {
        GooseFEM::Iterate::StopList residuals(5);

        std::vector<double> res = {5e+0, 5e+1, 5e-1, 5e-2, 5e-3, 5e-4, 4e-4, 5e-4, 4e-4, 4e-4};
        std::vector<bool> conv = {
            false, false, false, false, false, false, false, false, false, true
        };

        for (size_t i = 0; i < res.size(); ++i) {
            residuals.roll_insert(res[i]);
            REQUIRE(residuals.all_less(1e-3) == conv[i]);
        }
    }

    SECTION("StopList::descending && all_less - sorted input")
    {
        GooseFEM::Iterate::StopList residuals(5);

        std::vector<double> res = {5e+0, 5e+1, 5e-1, 5e-2, 5e-3, 5e-4, 4e-4, 3e-4, 2e-4, 1e-4};
        std::vector<bool> conv = {
            false, false, false, false, false, false, false, false, false, true
        };

        for (size_t i = 0; i < res.size(); ++i) {
            residuals.roll_insert(res[i]);
            REQUIRE((residuals.descending() && residuals.all_less(1e-3)) == conv[i]);
        }
    }

    SECTION("StopList::descending && all_less - unsorted input")
    {
        GooseFEM::Iterate::StopList residuals(5);

        std::vector<double> res = {
            5e+0, 5e+1, 5e-1, 5e-2, 5e-3, 5e-4, 4e-4, 3e-4, 2e-4, 3e-4, 2e-4, 1e-4, 9e-5, 8e-5
        };
        std::vector<bool> conv = {
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            true
        };

        for (size_t i = 0; i < res.size(); ++i) {
            residuals.roll_insert(res[i]);
            REQUIRE((residuals.descending() && residuals.all_less(1e-3)) == conv[i]);
        }
    }
}

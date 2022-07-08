#include <GooseFEM/GooseFEM.h>
#include <catch2/catch_all.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xrandom.hpp>

#define ISCLOSE(a, b) REQUIRE_THAT((a), Catch::Matchers::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Element", "Element.h")
{

    SECTION("QuadratureBase - AsTensor - e.g. Quad4")
    {
        size_t nelem = 3;
        size_t nip = 4;
        size_t nne = 4;
        size_t ndim = 2;
        xt::xtensor<double, 3> elemvec = xt::empty<double>({nelem, nne, ndim});
        GooseFEM::Element::Quad4::Quadrature quad(elemvec);

        xt::xtensor<double, 2> qscalar = xt::random::rand<double>({nelem, nip});
        auto qtensor = quad.AsTensor<2>(qscalar);

        REQUIRE(xt::has_shape(quad.allocate_qscalar<double>(), qscalar.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qscalar<double>(0.0), qscalar.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qtensor<2, double>(), qtensor.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qtensor<2, double>(0.0), qtensor.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qtensor<double>(2), qtensor.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qtensor<double>(2, 0.0), qtensor.shape()));

        for (size_t e = 0; e < nelem; ++e) {
            for (size_t q = 0; q < nip; ++q) {
                REQUIRE(xt::allclose(xt::view(qtensor, xt::keep(e), xt::keep(q)), qscalar(e, q)));
            }
        }
    }

    SECTION("QuadratureBase - AsTensor - e.g. Quad4")
    {
        size_t nelem = 3;
        size_t nip = 4;
        size_t nne = 4;
        size_t ndim = 2;
        xt::xtensor<double, 3> elemvec = xt::empty<double>({nelem, nne, ndim});
        GooseFEM::Element::Quad4::Quadrature quad(elemvec);

        xt::xtensor<double, 2> qscalar = xt::random::rand<double>({nelem, nip});
        auto qtensor = quad.AsTensor<2>(qscalar);

        REQUIRE(xt::has_shape(quad.allocate_qscalar<double>(), qscalar.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qscalar<double>(0.0), qscalar.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qtensor<2, double>(), qtensor.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qtensor<2, double>(0.0), qtensor.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qtensor<double>(2), qtensor.shape()));
        REQUIRE(xt::has_shape(quad.allocate_qtensor<double>(2, 0.0), qtensor.shape()));

        for (size_t e = 0; e < nelem; ++e) {
            for (size_t q = 0; q < nip; ++q) {
                REQUIRE(xt::allclose(xt::view(qtensor, xt::keep(e), xt::keep(q)), qscalar(e, q)));
            }
        }
    }
}

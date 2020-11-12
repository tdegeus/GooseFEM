
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Matrix", "Matrix.h")
{
    SECTION("solve")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        size_t nne = mesh.nne();
        size_t ndim = mesh.ndim();
        size_t nelem = mesh.nelem();
        size_t nnode = mesh.nnode();

        xt::xtensor<double, 3> a = xt::empty<double>({nelem, nne * ndim, nne * ndim});
        xt::xtensor<double, 1> b = xt::random::rand<double>({nnode * ndim});

        for (size_t e = 0; e < nelem; ++e) {
            xt::xtensor<double, 2> ae = xt::random::rand<double>({nne * ndim, nne * ndim});
            ae = (ae + xt::transpose(ae)) / 2.0;
            xt::view(a, e, xt::all(), xt::all()) = ae;
        }

        GooseFEM::Matrix A(mesh.conn(), mesh.dofs());
        GooseFEM::MatrixSolver<> Solver;
        A.assemble(a);
        xt::xtensor<double, 1> C = A.Dot(b);
        xt::xtensor<double, 1> B = Solver.Solve(A, C);

        REQUIRE(B.size() == b.size());
        REQUIRE(xt::allclose(B, b));
    }

    SECTION("set/dot/solve")
    {
        xt::xtensor<double, 2> A = xt::random::rand<double>({5, 5});
        xt::xtensor<double, 1> x = xt::random::rand<double>({5});
        xt::xtensor<double, 1> b = xt::zeros<double>({5});

        for (size_t i = 0; i < A.shape(0); ++i) {
            for (size_t j = 0; j < A.shape(1); ++j) {
                b(i) += A(i, j) * x(j);
            }
        }

        xt::xtensor<size_t, 2> conn = xt::zeros<size_t>({1, 5});
        xt::xtensor<size_t, 2> dofs = xt::zeros<size_t>({5, 2});

        GooseFEM::Matrix K(conn, dofs);
        GooseFEM::MatrixSolver<> Solver;
        K.set(xt::arange<size_t>(5), xt::arange<size_t>(5), A);

        auto B = K.Dot(x);
        auto X = Solver.Solve(K, b);

        REQUIRE(xt::allclose(X, x));
        REQUIRE(xt::allclose(B, b));
    }
}

#define CATCH_CONFIG_MAIN // tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::MatrixPartitioned", "MatrixPartitioned.h")
{
    SECTION("solve")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        size_t nne = mesh.nne();
        size_t ndim = mesh.ndim();
        size_t nelem = mesh.nelem();
        size_t nnode = mesh.nnode();
        auto dofs = mesh.dofs();
        size_t nnp = xt::amax(dofs)();
        nnp = (nnp - nnp % 2) / 2;
        xt::xtensor<size_t, 1> iip = xt::arange<size_t<(nnp);

        xt::xtensor<double, 3> a = xt::empty<double>({nelem, nne * ndim, nne * ndim});
        xt::xtensor<double, 1> b = xt::random::rand<double>({nnode * ndim});

        for (size_t e = 0; e < nelem; ++e) {
            xt::xtensor<double, 2> ae = xt::random::rand<double>({nne * ndim, nne * ndim});
            ae = (ae + xt::transpose(ae)) / 2.0;
            xt::view(a, e, xt::all(), xt::all()) = ae;
        }

        GooseFEM::MatrixPartitioned A(mesh.conn(), dofs, iip);
        GooseFEM::MatrixPartitionedSolver<> Solver;
        A.assemble(a);
        xt::xtensor<double, 1> C = A.Dot(b);
        xt::xtensor<double, 1> B = Solver.Solve(A, C);

        REQUIRE(B.size() == b.size());
        REQUIRE(xt::allclose(B, b));

        // check that allocating a different Solver instance still works
        GooseFEM::MatrixPartitionedSolver<> NewSolver;
        xt::xtensor<double, 1> NB = NewSolver.Solve(A, C);

        REQUIRE(NB.size() == b.size());
        REQUIRE(xt::allclose(NB, b));
    }

    SECTION("set/add/dot/solve - dofval")
    {
        xt::xtensor<double, 2> a = xt::random::rand<double>({10, 10});
        xt::xtensor<double, 1> x = xt::random::rand<double>({10});
        xt::xtensor<double, 1> b = xt::zeros<double>({10});

        xt::xtensor<double, 2> A = a + xt::transpose(a);

        for (size_t i = 0; i < A.shape(0); ++i) {
            for (size_t j = 0; j < A.shape(1); ++j) {
                b(i) += A(i, j) * x(j);
            }
        }

        xt::xtensor<size_t, 2> conn = xt::zeros<size_t>({1, 5});
        xt::xtensor<size_t, 2> dofs = xt::arange<size_t>(10).reshape({5, 2});
        xt::xtensor<size_t, 1> iip = {0, 2, 4};

        GooseFEM::MatrixPartitioned K(conn, dofs, iip);
        GooseFEM::MatrixPartitionedSolver<> Solver;
        K.set(xt::arange<size_t>(10), xt::arange<size_t>(10), a);
        K.add(xt::arange<size_t>(10), xt::arange<size_t>(10), xt::transpose(a));

        REQUIRE(xt::allclose(A, K.Todense()));
        REQUIRE(xt::allclose(b, K.Dot(x)));
        REQUIRE(xt::allclose(x, Solver.Solve(K, b)));
    }

    SECTION("set/add/dot/solve - nodevec")
    {
        xt::xtensor<double, 2> a = xt::random::rand<double>({10, 10});
        xt::xtensor<double, 2> x = xt::random::rand<double>({5, 2});
        xt::xtensor<double, 2> b = xt::zeros<double>({5, 2});

        xt::xtensor<double, 2> A = a + xt::transpose(a);

        for (size_t m = 0; m < x.shape(0); ++m) {
            for (size_t n = 0; n < x.shape(0); ++n) {
                for (size_t i = 0; i < x.shape(1); ++i) {
                    for (size_t j = 0; j < x.shape(1); ++j) {
                        b(m, i) += A(m * x.shape(1) + i, n * x.shape(1) + j) * x(n, j);
                    }
                }
            }
        }

        xt::xtensor<size_t, 2> conn = xt::zeros<size_t>({1, 5});
        xt::xtensor<size_t, 2> dofs = xt::arange<size_t>(10).reshape({5, 2});
        xt::xtensor<size_t, 1> iip = {0, 2, 4};

        GooseFEM::MatrixPartitioned K(conn, dofs, iip);
        GooseFEM::MatrixPartitionedSolver<> Solver;
        K.set(xt::arange<size_t>(10), xt::arange<size_t>(10), a);
        K.add(xt::arange<size_t>(10), xt::arange<size_t>(10), xt::transpose(a));

        REQUIRE(xt::allclose(A, K.Todense()));
        REQUIRE(xt::allclose(b, K.Dot(x)));
        REQUIRE(xt::allclose(x, Solver.Solve(K, b)));
    }
}

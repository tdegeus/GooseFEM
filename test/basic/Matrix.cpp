
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
}

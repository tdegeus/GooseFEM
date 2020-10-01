
#include "support.h"

TEST_CASE("GooseFEM::MatrixDiagonal", "MatrixDiagonal.h")
{

    SECTION("dot")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        xt::xtensor<double, 1> a = xt::random::rand<double>({mesh.nnode() * mesh.ndim()});
        xt::xtensor<double, 1> b = xt::random::rand<double>({mesh.nnode() * mesh.ndim()});
        xt::xtensor<double, 1> c = a * b;

        GooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
        A.set(a);
        xt::xtensor<double, 1> C = A.Dot(b);

        REQUIRE(C.size() == c.size());
        REQUIRE(xt::allclose(C, c));
    }

    SECTION("solve")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        xt::xtensor<double, 1> a = xt::random::rand<double>({mesh.nnode() * mesh.ndim()});
        xt::xtensor<double, 1> b = xt::random::rand<double>({mesh.nnode() * mesh.ndim()});
        xt::xtensor<double, 1> c = a * b;

        GooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
        A.set(a);
        xt::xtensor<double, 1> C = A.Dot(b);
        xt::xtensor<double, 1> B = A.Solve(C);

        REQUIRE(B.size() == b.size());
        REQUIRE(xt::allclose(B, b));
    }
}

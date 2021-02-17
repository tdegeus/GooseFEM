
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::ElementHex8", "ElementHex8.h")
{

    SECTION("int_N_scalar_NT_dV")
    {
        GooseFEM::Mesh::Hex8::Regular mesh(3, 3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::MatrixDiagonal mat(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Hex8::Quadrature quad(
            vec.AsElement(mesh.coor()),
            GooseFEM::Element::Hex8::Nodal::xi(),
            GooseFEM::Element::Hex8::Nodal::w());

        xt::xtensor<double, 2> rho = xt::ones<double>({mesh.nelem(), quad.nip()});

        mat.assemble(quad.Int_N_scalar_NT_dV(rho));
        auto M = mat.Todiagonal();

        REQUIRE(M.size() == vec.ndof());
        REQUIRE(xt::allclose(M, 1.));
    }

    SECTION("symGradN_vector")
    {
        GooseFEM::Mesh::Hex8::FineLayer mesh(27, 27, 27);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Hex8::Quadrature quad(vec.AsElement(mesh.coor()));

        xt::xtensor<double, 2> F = xt::zeros<double>({3, 3});
        xt::xtensor<double, 2> EPS = xt::zeros<double>({3, 3});

        F(0, 1) = 0.1;
        EPS(0, 1) = 0.05;
        EPS(1, 0) = 0.05;

        auto coor = mesh.coor();
        xt::xtensor<double, 2> disp = xt::zeros<double>(coor.shape());

        for (size_t n = 0; n < mesh.nnode(); ++n) {
            for (size_t i = 0; i < F.shape()[0]; ++i) {
                for (size_t j = 0; j < F.shape()[1]; ++j) {
                    disp(n, i) += F(i, j) * coor(n, j);
                }
            }
        }

        auto eps = quad.SymGradN_vector(vec.AsElement(disp));
        auto dV = quad.AsTensor<2>(quad.dV());
        auto epsbar = xt::average(eps, dV, {0, 1});

        REQUIRE(eps.shape()[0] == mesh.nelem());
        REQUIRE(eps.shape()[1] == quad.nip());
        REQUIRE(eps.shape()[2] == mesh.ndim());
        REQUIRE(eps.shape()[3] == mesh.ndim());

        for (size_t e = 0; e < mesh.nelem(); ++e) {
            for (size_t q = 0; q < quad.nip(); ++q) {
                REQUIRE(xt::allclose(xt::view(eps, e, q), EPS));
            }
        }

        REQUIRE(xt::allclose(epsbar, EPS));
    }

    SECTION("symGradN_vector, int_gradN_dot_tensor2s_dV")
    {
        GooseFEM::Mesh::Hex8::FineLayer mesh(27, 27, 27);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Hex8::Quadrature quad(vec.AsElement(mesh.coor()));

        xt::xtensor<double, 2> F = xt::zeros<double>({3, 3});

        F(0, 1) = 0.1;

        auto coor = mesh.coor();
        xt::xtensor<double, 2> disp = xt::zeros<double>(coor.shape());

        for (size_t n = 0; n < mesh.nnode(); ++n) {
            for (size_t i = 0; i < F.shape()[0]; ++i) {
                for (size_t j = 0; j < F.shape()[1]; ++j) {
                    disp(n, i) += F(i, j) * coor(n, j);
                }
            }
        }

        auto eps = quad.SymGradN_vector(vec.AsElement(disp));
        auto Fi = vec.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(eps));

        REQUIRE(Fi.size() == vec.ndof());
        REQUIRE(xt::allclose(Fi, 0.));
    }
}

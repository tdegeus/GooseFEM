
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

TEST_CASE("GooseFEM::ElementQuad4", "ElementQuad4.h")
{
    SECTION("GradN")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));
        auto dNdx = quad.GradN();
        REQUIRE(xt::has_shape(dNdx, {mesh.nelem(), quad.nip(), mesh.nne(), mesh.ndim()}));
    }

    SECTION("dV")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));
        auto dV = quad.dV();
        REQUIRE(xt::allclose(dV, 0.5 * 0.5));
    }

    SECTION("interp_N_vector")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Quad4::Quadrature quad(vector.AsElement(mesh.coor()));

        auto u = vector.AllocateNodevec(1.0);
        auto ue = vector.AsElement(u);
        auto uq = quad.Interp_N_vector(ue);

        REQUIRE(xt::allclose(uq, 1.0));
    }

    SECTION("GradN_vector")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));

        xt::xtensor<double, 2> F = xt::zeros<double>({2, 2});
        xt::xtensor<double, 2> EPS = xt::zeros<double>({2, 2});

        F(0, 1) = 0.1;

        auto coor = mesh.coor();
        auto disp = xt::zeros_like(coor);

        for (size_t n = 0; n < mesh.nnode(); ++n) {
            for (size_t i = 0; i < F.shape()[0]; ++i) {
                for (size_t j = 0; j < F.shape()[1]; ++j) {
                    disp(n, i) += F(i, j) * coor(n, j);
                }
            }
        }

        auto f = quad.GradN_vector(vec.AsElement(disp));

        REQUIRE(xt::has_shape(f, {mesh.nelem(), quad.nip(), mesh.ndim(), mesh.ndim()}));

        for (size_t e = 0; e < mesh.nelem(); ++e) {
            for (size_t q = 0; q < quad.nip(); ++q) {
                REQUIRE(xt::allclose(xt::view(f, e, q), xt::transpose(F)));
            }
        }
    }

    SECTION("GradN_vector_T")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));

        xt::xtensor<double, 2> F = xt::zeros<double>({2, 2});
        xt::xtensor<double, 2> EPS = xt::zeros<double>({2, 2});

        F(0, 1) = 0.1;

        auto coor = mesh.coor();
        auto disp = xt::zeros_like(coor);

        for (size_t n = 0; n < mesh.nnode(); ++n) {
            for (size_t i = 0; i < F.shape()[0]; ++i) {
                for (size_t j = 0; j < F.shape()[1]; ++j) {
                    disp(n, i) += F(i, j) * coor(n, j);
                }
            }
        }

        auto f = quad.GradN_vector_T(vec.AsElement(disp));

        REQUIRE(xt::has_shape(f, {mesh.nelem(), quad.nip(), mesh.ndim(), mesh.ndim()}));

        for (size_t e = 0; e < mesh.nelem(); ++e) {
            for (size_t q = 0; q < quad.nip(); ++q) {
                REQUIRE(xt::allclose(xt::view(f, e, q), F));
            }
        }
    }

    SECTION("SymGradN_vector")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));

        xt::xtensor<double, 2> F = xt::zeros<double>({2, 2});
        xt::xtensor<double, 2> EPS = xt::zeros<double>({2, 2});

        F(0, 1) = 0.1;
        EPS(0, 1) = 0.05;
        EPS(1, 0) = 0.05;

        auto coor = mesh.coor();
        auto disp = xt::zeros_like(coor);

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

        REQUIRE(xt::has_shape(eps, {mesh.nelem(), quad.nip(), mesh.ndim(), mesh.ndim()}));
        REQUIRE(xt::has_shape(dV, {mesh.nelem(), quad.nip(), mesh.ndim(), mesh.ndim()}));
        REQUIRE(xt::has_shape(epsbar, {mesh.ndim(), mesh.ndim()}));

        for (size_t e = 0; e < mesh.nelem(); ++e) {
            for (size_t q = 0; q < quad.nip(); ++q) {
                REQUIRE(xt::allclose(xt::view(eps, e, q), EPS));
            }
        }

        REQUIRE(xt::allclose(epsbar, EPS));
    }

    SECTION("Int_N_scalar_NT_dV")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::MatrixDiagonal mat(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Quad4::Quadrature quad(
            vec.AsElement(mesh.coor()),
            GooseFEM::Element::Quad4::Nodal::xi(),
            GooseFEM::Element::Quad4::Nodal::w());

        xt::xtensor<double, 2> rho = xt::ones<double>({mesh.nelem(), quad.nip()});

        mat.assemble(quad.Int_N_scalar_NT_dV(rho));
        auto M = mat.Todiagonal();

        REQUIRE(M.size() == vec.ndof());
        REQUIRE(xt::allclose(M, 1.0));
    }

    SECTION("Int_gradN_dot_tensor2_dV")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));

        xt::xtensor<double, 2> F = xt::zeros<double>({2, 2});

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
        REQUIRE(xt::allclose(Fi, 0.0));
    }
}

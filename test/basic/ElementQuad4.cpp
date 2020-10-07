
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::ElementQuad4", "ElementQuad4.h")
{

    SECTION("dV - Gauss")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        auto coor = mesh.coor();
        auto top = mesh.nodesTopEdge();
        auto right = mesh.nodesRightEdge();

        xt::view(coor, xt::keep(top), xt::keep(1)) += 1.0;
        xt::view(coor, xt::keep(right), xt::keep(0)) += 1.0;

        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(coor));

        auto dV = quad.dV();
        auto dV_tensor = quad.AsTensor<2>(dV);

        REQUIRE(xt::allclose(xt::view(dV, xt::keep(0)), 1. / 4.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(1)), 2. / 4.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(2)), 2. / 4.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(3)), 4. / 4.));

        for (size_t e = 0; e < mesh.nelem(); ++e) {
            for (size_t q = 0; q < quad.nip(); ++q) {
                REQUIRE(xt::allclose(xt::view(dV_tensor, xt::keep(e), xt::keep(q)), dV(e, q)));
            }
        }
    }

    SECTION("dV - Nodal")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        auto coor = mesh.coor();
        auto top = mesh.nodesTopEdge();
        auto right = mesh.nodesRightEdge();

        xt::view(coor, xt::keep(top), xt::keep(1)) += 1.0;
        xt::view(coor, xt::keep(right), xt::keep(0)) += 1.0;

        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(
            vec.AsElement(coor),
            GooseFEM::Element::Quad4::Nodal::xi(),
            GooseFEM::Element::Quad4::Nodal::w());

        auto dV = quad.dV();
        auto dV_tensor = quad.AsTensor<2>(dV);

        REQUIRE(xt::allclose(xt::view(dV, xt::keep(0)), 1. / 4.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(1)), 2. / 4.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(2)), 2. / 4.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(3)), 4. / 4.));

        for (size_t e = 0; e < mesh.nelem(); ++e) {
            for (size_t q = 0; q < quad.nip(); ++q) {
                REQUIRE(xt::allclose(xt::view(dV_tensor, xt::keep(e), xt::keep(q)), dV(e, q)));
            }
        }
    }

    SECTION("int_N_scalar_NT_dV")
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
        auto M = mat.AsDiagonal();

        REQUIRE(M.size() == vec.ndof());
        REQUIRE(xt::allclose(M, 1.));
    }

    SECTION("symGradN_vector")
    {
        GooseFEM::Mesh::Quad4::FineLayer mesh(27, 27);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));

        xt::xtensor<double, 2> F = xt::zeros<double>({2, 2});
        xt::xtensor<double, 2> EPS = xt::zeros<double>({2, 2});

        F(0, 1) = 0.1;
        EPS(0, 1) = 0.05;
        EPS(1, 0) = 0.05;

        xt::xtensor<double, 2> coor = mesh.coor();
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
        GooseFEM::Mesh::Quad4::FineLayer mesh(27, 27);
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
        REQUIRE(xt::allclose(Fi, 0.));
    }
}

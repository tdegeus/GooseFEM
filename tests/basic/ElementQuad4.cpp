#include <GooseFEM/GooseFEM.h>
#include <catch2/catch_all.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xrandom.hpp>

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
        REQUIRE(xt::has_shape(dV, quad.shape_qtensor<0>()));
        REQUIRE(xt::has_shape(dV, quad.shape_qscalar()));
        REQUIRE(xt::allclose(dV, 0.5 * 0.5));
    }

    SECTION("InterpQuad_vector")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Quad4::Quadrature quad(vector.AsElement(mesh.coor()));

        auto u = vector.allocate_nodevec(1.0);
        auto ue = vector.AsElement(u);
        auto uq = quad.InterpQuad_vector(ue);

        REQUIRE(xt::allclose(uq, 1.0));
    }

    SECTION("GradN_vector")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));

        size_t td = mesh.ndim();
        xt::xtensor<double, 2> F = xt::zeros<double>({td, td});
        xt::xtensor<double, 2> EPS = xt::zeros<double>({td, td});

        F(0, 1) = 0.1;

        auto coor = mesh.coor();
        auto disp = xt::zeros_like(coor);

        for (size_t n = 0; n < mesh.nnode(); ++n) {
            for (size_t i = 0; i < mesh.ndim(); ++i) {
                for (size_t j = 0; j < mesh.ndim(); ++j) {
                    disp(n, i) += F(i, j) * coor(n, j);
                }
            }
        }

        auto f = quad.GradN_vector(vec.AsElement(disp));

        REQUIRE(xt::has_shape(f, {mesh.nelem(), quad.nip(), td, td}));
        REQUIRE(xt::has_shape(f, quad.shape_qtensor<2>()));

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

        size_t td = mesh.ndim();
        xt::xtensor<double, 2> F = xt::zeros<double>({td, td});
        xt::xtensor<double, 2> EPS = xt::zeros<double>({td, td});

        F(0, 1) = 0.1;

        auto coor = mesh.coor();
        auto disp = xt::zeros_like(coor);

        for (size_t n = 0; n < mesh.nnode(); ++n) {
            for (size_t i = 0; i < mesh.ndim(); ++i) {
                for (size_t j = 0; j < mesh.ndim(); ++j) {
                    disp(n, i) += F(i, j) * coor(n, j);
                }
            }
        }

        auto f = quad.GradN_vector_T(vec.AsElement(disp));

        REQUIRE(xt::has_shape(f, {mesh.nelem(), quad.nip(), td, td}));
        REQUIRE(xt::has_shape(f, quad.shape_qtensor<2>()));

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

        size_t td = mesh.ndim();
        xt::xtensor<double, 2> F = xt::zeros<double>({td, td});
        xt::xtensor<double, 2> EPS = xt::zeros<double>({td, td});

        F(0, 1) = 0.1;
        EPS(0, 1) = 0.05;
        EPS(1, 0) = 0.05;

        auto coor = mesh.coor();
        auto disp = xt::zeros_like(coor);

        for (size_t n = 0; n < mesh.nnode(); ++n) {
            for (size_t i = 0; i < mesh.ndim(); ++i) {
                for (size_t j = 0; j < mesh.ndim(); ++j) {
                    disp(n, i) += F(i, j) * coor(n, j);
                }
            }
        }

        auto eps = quad.SymGradN_vector(vec.AsElement(disp));
        auto dV = quad.AsTensor<2>(quad.dV());
        auto epsbar = xt::average(eps, dV, {0, 1});

        REQUIRE(xt::has_shape(eps, {mesh.nelem(), quad.nip(), td, td}));
        REQUIRE(xt::has_shape(eps, quad.shape_qtensor<2>()));
        REQUIRE(xt::has_shape(dV, {mesh.nelem(), quad.nip(), td, td}));
        REQUIRE(xt::has_shape(dV, quad.shape_qtensor<2>()));
        REQUIRE(xt::has_shape(epsbar, {td, td}));

        for (size_t e = 0; e < mesh.nelem(); ++e) {
            for (size_t q = 0; q < quad.nip(); ++q) {
                REQUIRE(xt::allclose(xt::view(eps, e, q), EPS));
            }
        }

        REQUIRE(xt::allclose(epsbar, EPS));
    }

    SECTION("Int_N_vector_dV")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Quad4::Quadrature quad(vector.AsElement(mesh.coor()));

        xt::xtensor<double, 3> qvector = xt::empty<double>({mesh.nelem(), quad.nip(), mesh.ndim()});
        xt::view(qvector, xt::all(), xt::all(), 0) = 2.0;
        xt::view(qvector, xt::all(), xt::all(), 1) = 5.0;

        auto elemvec = quad.Int_N_vector_dV(qvector);
        auto nodevec = vector.AssembleNode(elemvec);

        REQUIRE(xt::allclose(xt::view(nodevec, xt::all(), 0), 2.0));
        REQUIRE(xt::allclose(xt::view(nodevec, xt::all(), 1), 5.0));
    }

    SECTION("Int_N_scalar_NT_dV")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::MatrixDiagonal mat(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Quad4::Quadrature quad(
            vec.AsElement(mesh.coor()),
            GooseFEM::Element::Quad4::Nodal::xi(),
            GooseFEM::Element::Quad4::Nodal::w()
        );

        xt::xtensor<double, 2> rho = xt::ones<double>({mesh.nelem(), quad.nip()});

        mat.assemble(quad.Int_N_scalar_NT_dV(rho));

        REQUIRE(mat.data().size() == vec.ndof());
        REQUIRE(xt::allclose(mat.data(), 1.0));
    }

    SECTION("Int_gradN_dot_tensor2_dV")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(3, 3);
        GooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());
        GooseFEM::Element::Quad4::Quadrature quad(vec.AsElement(mesh.coor()));

        size_t td = mesh.ndim();
        xt::xtensor<double, 2> F = xt::zeros<double>({td, td});

        F(0, 1) = 0.1;

        auto coor = mesh.coor();
        xt::xtensor<double, 2> disp = xt::zeros<double>(coor.shape());

        for (size_t n = 0; n < mesh.nnode(); ++n) {
            for (size_t i = 0; i < mesh.ndim(); ++i) {
                for (size_t j = 0; j < mesh.ndim(); ++j) {
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

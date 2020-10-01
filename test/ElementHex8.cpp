
#include "support.h"

TEST_CASE("GooseFEM::ElementHex8", "ElementHex8.h")
{

    using T2 = xt::xtensor_fixed<double, xt::xshape<3, 3>>;

    SECTION("dV - Gauss")
    {
        GooseFEM::Mesh::Hex8::Regular mesh(2, 2, 2);

        xt::xtensor<double, 2> coor = mesh.coor();

        xt::xtensor<size_t, 1> top = mesh.nodesTop();
        xt::xtensor<size_t, 1> right = mesh.nodesRight();
        xt::xtensor<size_t, 1> back = mesh.nodesBack();

        xt::view(coor, xt::keep(back), xt::keep(2)) += 1.;
        xt::view(coor, xt::keep(top), xt::keep(1)) += 1.;
        xt::view(coor, xt::keep(right), xt::keep(0)) += 1.;

        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());

        GooseFEM::Element::Hex8::Quadrature quad(vec.AsElement(coor));

        xt::xtensor<double, 2> dV = quad.DV();

        xt::xtensor<double, 4> dV_tensor = quad.DV(2);

        REQUIRE(xt::allclose(xt::view(dV, xt::keep(0)), 1. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(1)), 2. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(2)), 2. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(3)), 4. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(4)), 2. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(5)), 4. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(6)), 4. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(7)), 8. / 8.));

        for (size_t e = 0; e < mesh.nelem(); ++e)
            for (size_t q = 0; q < quad.nip(); ++q)
                REQUIRE(xt::allclose(xt::view(dV_tensor, xt::keep(e), xt::keep(q)), dV(e, q)));
    }

    SECTION("dV - Nodal")
    {
        GooseFEM::Mesh::Hex8::Regular mesh(2, 2, 2);

        xt::xtensor<double, 2> coor = mesh.coor();

        xt::xtensor<size_t, 1> top = mesh.nodesTop();
        xt::xtensor<size_t, 1> right = mesh.nodesRight();
        xt::xtensor<size_t, 1> back = mesh.nodesBack();

        xt::view(coor, xt::keep(back), xt::keep(2)) += 1.;
        xt::view(coor, xt::keep(top), xt::keep(1)) += 1.;
        xt::view(coor, xt::keep(right), xt::keep(0)) += 1.;

        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());

        GooseFEM::Element::Hex8::Quadrature quad(
            vec.AsElement(coor),
            GooseFEM::Element::Hex8::Nodal::xi(),
            GooseFEM::Element::Hex8::Nodal::w());

        xt::xtensor<double, 2> dV = quad.DV();

        xt::xtensor<double, 4> dV_tensor = quad.DV(2);

        REQUIRE(xt::allclose(xt::view(dV, xt::keep(0)), 1. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(1)), 2. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(2)), 2. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(3)), 4. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(4)), 2. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(5)), 4. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(6)), 4. / 8.));
        REQUIRE(xt::allclose(xt::view(dV, xt::keep(7)), 8. / 8.));

        for (size_t e = 0; e < mesh.nelem(); ++e)
            for (size_t q = 0; q < quad.nip(); ++q)
                REQUIRE(xt::allclose(xt::view(dV_tensor, xt::keep(e), xt::keep(q)), dV(e, q)));
    }

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

        xt::xtensor<double, 1> M = mat.AsDiagonal();

        REQUIRE(M.size() == vec.ndof());

        REQUIRE(xt::allclose(M, 1.));
    }

    SECTION("symGradN_vector")
    {
        GooseFEM::Mesh::Hex8::FineLayer mesh(27, 27, 27);

        GooseFEM::Vector vec(mesh.conn(), mesh.dofs());

        GooseFEM::Element::Hex8::Quadrature quad(vec.AsElement(mesh.coor()));

        T2 F = xt::zeros<double>({3, 3});
        T2 EPS = xt::zeros<double>({3, 3});

        F(0, 1) = 0.1;
        EPS(0, 1) = 0.05;
        EPS(1, 0) = 0.05;

        xt::xtensor<double, 2> coor = mesh.coor();
        xt::xtensor<double, 2> disp = xt::zeros<double>(coor.shape());

        for (size_t n = 0; n < mesh.nnode(); ++n)
            for (size_t i = 0; i < F.shape()[0]; ++i)
                for (size_t j = 0; j < F.shape()[1]; ++j)
                    disp(n, i) += F(i, j) * coor(n, j);

        xt::xtensor<double, 4> eps = quad.SymGradN_vector(vec.AsElement(disp));

        xt::xtensor<double, 4> dV = quad.DV(2);

        auto epsbar = xt::average(eps, dV, {0, 1});

        REQUIRE(eps.shape()[0] == mesh.nelem());
        REQUIRE(eps.shape()[1] == quad.nip());
        REQUIRE(eps.shape()[2] == mesh.ndim());
        REQUIRE(eps.shape()[3] == mesh.ndim());

        for (size_t e = 0; e < mesh.nelem(); ++e)
            for (size_t q = 0; q < quad.nip(); ++q)
                REQUIRE(xt::allclose(xt::view(eps, e, q), EPS));

        REQUIRE(xt::allclose(epsbar, EPS));
    }

    SECTION("symGradN_vector, int_gradN_dot_tensor2s_dV")
    {
        GooseFEM::Mesh::Hex8::FineLayer mesh(27, 27, 27);

        GooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());

        GooseFEM::Element::Hex8::Quadrature quad(vec.AsElement(mesh.coor()));

        T2 F = xt::zeros<double>({3, 3});

        F(0, 1) = 0.1;

        xt::xtensor<double, 2> coor = mesh.coor();
        xt::xtensor<double, 2> disp = xt::zeros<double>(coor.shape());

        for (size_t n = 0; n < mesh.nnode(); ++n)
            for (size_t i = 0; i < F.shape()[0]; ++i)
                for (size_t j = 0; j < F.shape()[1]; ++j)
                    disp(n, i) += F(i, j) * coor(n, j);

        xt::xtensor<double, 4> eps = quad.SymGradN_vector(vec.AsElement(disp));

        xt::xtensor<double, 1> Fi = vec.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(eps));

        REQUIRE(Fi.size() == vec.ndof());

        REQUIRE(xt::allclose(Fi, 0.));
    }
}

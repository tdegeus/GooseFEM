
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>
#include <GMatElastoPlasticQPot/Cartesian2d.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.0e-12));

TEST_CASE("hybrid-elastoplasticqpot", "GooseFEM.h")
{

    SECTION("Vector/Matrix - GMatElastoPlasticQPot")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N, 2 * N);
        xt::xtensor<size_t, 1> elem_b = xt::concatenate(xt::xtuple(
            xt::arange<size_t>(N),
            xt::arange<size_t>(2 * N, nelem)));
        size_t nelem_a = elem_a.size();
        size_t nelem_b = elem_b.size();

        auto dofs = mesh.dofs();
        auto conn = mesh.conn();
        xt::xtensor<size_t, 2> conn_a = xt::view(conn, xt::keep(elem_a), xt::all());
        xt::xtensor<size_t, 2> conn_b = xt::view(conn, xt::keep(elem_b), xt::all());

        GooseFEM::Vector vector(conn, dofs);
        GooseFEM::Vector vector_a(conn_a, dofs);
        GooseFEM::Vector vector_b(conn_b, dofs);

        GooseFEM::Matrix K_b(conn_b, dofs);

        auto coor = mesh.coor();

        GooseFEM::Element::Quad4::Quadrature quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::Quadrature quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::Quadrature quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat({nelem, nip});
        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat_a({nelem_a, nip});
        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat_b({nelem_b, nip});

        xt::xtensor<double, 2> epsy_a = 0.2 * xt::random::rand<double>(std::array<size_t, 2>{nelem_a, 100});
        epsy_a = xt::cumsum(epsy_a, 1);

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::xtensor<size_t, 2> idx = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_a), xt::all()) = 1;
            xt::view(idx, xt::keep(elem_a), xt::all()) = xt::arange<size_t>(nelem_a).reshape({-1, 1});
            mat.setCusp(I, idx, 3.0 * xt::ones<double>({nelem_a}), 4.0 * xt::ones<double>({nelem_a}), epsy_a);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_b), xt::all()) = 1;
            mat.setElastic(I, 5.0, 6.0);
        }

        {
            xt::xtensor<size_t, 2> I = xt::ones<size_t>({nelem_a, nip});
            xt::xtensor<size_t, 2> idx = xt::zeros<size_t>({nelem_a, nip});
            xt::view(idx, xt::range(0, nelem_a), xt::all()) = xt::arange<size_t>(nelem_a).reshape({-1, 1});
            mat_a.setCusp(I, idx, 3.0 * xt::ones<double>({nelem_a}), 4.0 * xt::ones<double>({nelem_a}), epsy_a);
        }

        {
            xt::xtensor<size_t, 2> I = xt::ones<size_t>({nelem_b, nip});
            mat_b.setElastic(I, 5.0, 6.0);
        }

        mat.check();
        mat_a.check();
        mat_b.check();

        for (size_t iter = 0; iter < 10; ++iter) {

            xt::xtensor<double, 2> disp = xt::random::rand<double>(coor.shape());

            auto Eps = quad.SymGradN_vector(vector.AsElement(disp));
            auto Eps_a = quad_a.SymGradN_vector(vector_a.AsElement(disp));
            auto Eps_b = quad_b.SymGradN_vector(vector_b.AsElement(disp));

            mat.setStrain(Eps);
            mat_a.setStrain(Eps_a);
            mat_b.setStrain(Eps_b);

            auto Sig = mat.Stress();
            auto Sig_a = mat_a.Stress();
            auto Sig_b = mat_b.Stress();

            auto C_b = mat_b.Tangent();

            K_b.assemble(quad_b.Int_gradN_dot_tensor4_dot_gradNT_dV(C_b));

            auto f = vector.AssembleNode(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto f_a = vector_a.AssembleNode(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto f_b = vector_b.AssembleNode(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(f, f_a + f_b));
            REQUIRE(xt::allclose(f_b, K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_a + K_b.Dot(disp)));

            auto fd = vector.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto fd_a = vector_a.AssembleDofs(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto fd_b = vector_b.AssembleDofs(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(fd, fd_a + fd_b));
            REQUIRE(xt::allclose(fd_b, K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_a + K_b.Dot(vector.AsDofs(disp))));

            auto Sig_c = decltype(Sig)::from_shape(Sig.shape());
            xt::view(Sig_c, xt::keep(elem_a)) = Sig_a;
            xt::view(Sig_c, xt::keep(elem_b)) = Sig_b;

            REQUIRE(xt::allclose(Sig, Sig_c));
        }
    }

    SECTION("Vector/Matrix - GMatElastoPlasticQPot - Periodic")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N, 2 * N);
        xt::xtensor<size_t, 1> elem_b = xt::concatenate(xt::xtuple(
            xt::arange<size_t>(N),
            xt::arange<size_t>(2 * N, nelem)));
        size_t nelem_a = elem_a.size();
        size_t nelem_b = elem_b.size();

        auto dofs = mesh.dofsPeriodic();
        auto conn = mesh.conn();
        xt::xtensor<size_t, 2> conn_a = xt::view(conn, xt::keep(elem_a), xt::all());
        xt::xtensor<size_t, 2> conn_b = xt::view(conn, xt::keep(elem_b), xt::all());

        GooseFEM::Vector vector(conn, dofs);
        GooseFEM::Vector vector_a(conn_a, dofs);
        GooseFEM::Vector vector_b(conn_b, dofs);

        GooseFEM::Matrix K_b(conn_b, dofs);

        auto coor = mesh.coor();

        GooseFEM::Element::Quad4::Quadrature quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::Quadrature quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::Quadrature quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat({nelem, nip});
        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat_a({nelem_a, nip});
        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat_b({nelem_b, nip});

        xt::xtensor<double, 2> epsy_a = 0.2 * xt::random::rand<double>(std::array<size_t, 2>{nelem_a, 100});
        epsy_a = xt::cumsum(epsy_a, 1);

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::xtensor<size_t, 2> idx = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_a), xt::all()) = 1;
            xt::view(idx, xt::keep(elem_a), xt::all()) = xt::arange<size_t>(nelem_a).reshape({-1, 1});
            mat.setCusp(I, idx, 3.0 * xt::ones<double>({nelem_a}), 4.0 * xt::ones<double>({nelem_a}), epsy_a);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_b), xt::all()) = 1;
            mat.setElastic(I, 5.0, 6.0);
        }

        {
            xt::xtensor<size_t, 2> I = xt::ones<size_t>({nelem_a, nip});
            xt::xtensor<size_t, 2> idx = xt::zeros<size_t>({nelem_a, nip});
            xt::view(idx, xt::range(0, nelem_a), xt::all()) = xt::arange<size_t>(nelem_a).reshape({-1, 1});
            mat_a.setCusp(I, idx, 3.0 * xt::ones<double>({nelem_a}), 4.0 * xt::ones<double>({nelem_a}), epsy_a);
        }

        {
            xt::xtensor<size_t, 2> I = xt::ones<size_t>({nelem_b, nip});
            mat_b.setElastic(I, 5.0, 6.0);
        }

        mat.check();
        mat_a.check();
        mat_b.check();

        for (size_t iter = 0; iter < 10; ++iter) {

            xt::xtensor<double, 2> disp = xt::random::rand<double>(coor.shape());
            disp = vector.AsNode(vector.AsDofs(disp));

            auto Eps = quad.SymGradN_vector(vector.AsElement(disp));
            auto Eps_a = quad_a.SymGradN_vector(vector_a.AsElement(disp));
            auto Eps_b = quad_b.SymGradN_vector(vector_b.AsElement(disp));

            mat.setStrain(Eps);
            mat_a.setStrain(Eps_a);
            mat_b.setStrain(Eps_b);

            auto Sig = mat.Stress();
            auto Sig_a = mat_a.Stress();
            auto Sig_b = mat_b.Stress();

            auto C_b = mat_b.Tangent();

            K_b.assemble(quad_b.Int_gradN_dot_tensor4_dot_gradNT_dV(C_b));

            auto f = vector.AssembleNode(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto f_a = vector_a.AssembleNode(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto f_b = vector_b.AssembleNode(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(f, f_a + f_b));
            REQUIRE(xt::allclose(f_b, K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_a + K_b.Dot(disp)));

            auto fd = vector.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto fd_a = vector_a.AssembleDofs(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto fd_b = vector_b.AssembleDofs(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(fd, fd_a + fd_b));
            REQUIRE(xt::allclose(fd_b, K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_a + K_b.Dot(vector.AsDofs(disp))));

            auto Sig_c = decltype(Sig)::from_shape(Sig.shape());
            xt::view(Sig_c, xt::keep(elem_a)) = Sig_a;
            xt::view(Sig_c, xt::keep(elem_b)) = Sig_b;

            REQUIRE(xt::allclose(Sig, Sig_c));
        }
    }

    SECTION("VectorPartitioned/Matrix - GMatElastoPlasticQPot - Semi-Periodic")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        size_t ndim = mesh.ndim();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N, 2 * N);
        xt::xtensor<size_t, 1> elem_b = xt::concatenate(xt::xtuple(
            xt::arange<size_t>(N),
            xt::arange<size_t>(2 * N, nelem)));
        size_t nelem_a = elem_a.size();
        size_t nelem_b = elem_b.size();

        auto dofs = mesh.dofs();
        auto conn = mesh.conn();
        xt::xtensor<size_t, 2> conn_a = xt::view(conn, xt::keep(elem_a), xt::all());
        xt::xtensor<size_t, 2> conn_b = xt::view(conn, xt::keep(elem_b), xt::all());

        auto left = mesh.nodesLeftOpenEdge();
        auto right = mesh.nodesRightOpenEdge();
        xt::view(dofs, xt::keep(right), 0) = xt::view(dofs, xt::keep(left), 0);
        xt::view(dofs, xt::keep(right), 1) = xt::view(dofs, xt::keep(left), 1);

        auto top = mesh.nodesTopEdge();
        auto bottom = mesh.nodesBottomEdge();
        size_t nfix = top.size();
        xt::xtensor<size_t, 1> iip = xt::empty<size_t>({2 * ndim * nfix});
        xt::view(iip, xt::range(0 * nfix, 1 * nfix)) = xt::view(dofs, xt::keep(bottom), 0);
        xt::view(iip, xt::range(1 * nfix, 2 * nfix)) = xt::view(dofs, xt::keep(bottom), 1);
        xt::view(iip, xt::range(2 * nfix, 3 * nfix)) = xt::view(dofs, xt::keep(top), 0);
        xt::view(iip, xt::range(3 * nfix, 4 * nfix)) = xt::view(dofs, xt::keep(top), 1);

        GooseFEM::VectorPartitioned vector(conn, dofs, iip);
        GooseFEM::VectorPartitioned vector_a(conn_a, dofs, iip);
        GooseFEM::VectorPartitioned vector_b(conn_b, dofs, iip);

        GooseFEM::Matrix K_b(conn_b, dofs);

        auto coor = mesh.coor();

        GooseFEM::Element::Quad4::Quadrature quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::Quadrature quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::Quadrature quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat({nelem, nip});
        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat_a({nelem_a, nip});
        GMatElastoPlasticQPot::Cartesian2d::Array<2> mat_b({nelem_b, nip});

        xt::xtensor<double, 2> epsy_a = 0.2 * xt::random::rand<double>(std::array<size_t, 2>{nelem_a, 100});
        epsy_a = xt::cumsum(epsy_a, 1);

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::xtensor<size_t, 2> idx = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_a), xt::all()) = 1;
            xt::view(idx, xt::keep(elem_a), xt::all()) = xt::arange<size_t>(nelem_a).reshape({-1, 1});
            mat.setCusp(I, idx, 3.0 * xt::ones<double>({nelem_a}), 4.0 * xt::ones<double>({nelem_a}), epsy_a);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_b), xt::all()) = 1;
            mat.setElastic(I, 5.0, 6.0);
        }

        {
            xt::xtensor<size_t, 2> I = xt::ones<size_t>({nelem_a, nip});
            xt::xtensor<size_t, 2> idx = xt::zeros<size_t>({nelem_a, nip});
            xt::view(idx, xt::range(0, nelem_a), xt::all()) = xt::arange<size_t>(nelem_a).reshape({-1, 1});
            mat_a.setCusp(I, idx, 3.0 * xt::ones<double>({nelem_a}), 4.0 * xt::ones<double>({nelem_a}), epsy_a);
        }

        {
            xt::xtensor<size_t, 2> I = xt::ones<size_t>({nelem_b, nip});
            mat_b.setElastic(I, 5.0, 6.0);
        }

        mat.check();
        mat_a.check();
        mat_b.check();

        for (size_t iter = 0; iter < 10; ++iter) {

            xt::xtensor<double, 2> disp = xt::random::rand<double>(coor.shape());
            disp = vector.AsNode(vector.AsDofs(disp));

            auto Eps = quad.SymGradN_vector(vector.AsElement(disp));
            auto Eps_a = quad_a.SymGradN_vector(vector_a.AsElement(disp));
            auto Eps_b = quad_b.SymGradN_vector(vector_b.AsElement(disp));

            mat.setStrain(Eps);
            mat_a.setStrain(Eps_a);
            mat_b.setStrain(Eps_b);

            auto Sig = mat.Stress();
            auto Sig_a = mat_a.Stress();
            auto Sig_b = mat_b.Stress();

            auto C_b = mat_b.Tangent();

            K_b.assemble(quad_b.Int_gradN_dot_tensor4_dot_gradNT_dV(C_b));

            auto f = vector.AssembleNode(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto f_a = vector_a.AssembleNode(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto f_b = vector_b.AssembleNode(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(f, f_a + f_b));
            REQUIRE(xt::allclose(f_b, K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_a + K_b.Dot(disp)));

            auto fd = vector.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto fd_a = vector_a.AssembleDofs(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto fd_b = vector_b.AssembleDofs(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(fd, fd_a + fd_b));
            REQUIRE(xt::allclose(fd_b, K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_a + K_b.Dot(vector.AsDofs(disp))));

            auto Sig_c = decltype(Sig)::from_shape(Sig.shape());
            xt::view(Sig_c, xt::keep(elem_a)) = Sig_a;
            xt::view(Sig_c, xt::keep(elem_b)) = Sig_b;

            REQUIRE(xt::allclose(Sig, Sig_c));
        }
    }

}

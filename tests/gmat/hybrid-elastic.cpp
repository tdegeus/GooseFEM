#define CATCH_CONFIG_MAIN // tells Catch to provide a main() - only do this in one cpp file
#include <Eigen/Eigen>
#include <GMatElastic/Cartesian3d.h>
#include <GooseFEM/GooseFEM.h>
#include <catch2/catch.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xrandom.hpp>

#define ISCLOSE(a, b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.0e-12));

TEST_CASE("hybrid-elastic", "GooseFEM.h")
{

    SECTION("Vector/Matrix - GMatElastic")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N, 2 * N);
        xt::xtensor<size_t, 1> elem_b =
            xt::concatenate(xt::xtuple(xt::arange<size_t>(N), xt::arange<size_t>(2 * N, nelem)));
        size_t nelem_a = elem_a.size();
        size_t nelem_b = elem_b.size();

        auto dofs = mesh.dofs();
        auto conn = mesh.conn();
        xt::xtensor<size_t, 2> conn_a = xt::view(conn, xt::keep(elem_a), xt::all());
        xt::xtensor<size_t, 2> conn_b = xt::view(conn, xt::keep(elem_b), xt::all());

        GooseFEM::Vector vector(conn, dofs);
        GooseFEM::Vector vector_a(conn_a, dofs);
        GooseFEM::Vector vector_b(conn_b, dofs);

        GooseFEM::Matrix K(conn, dofs);
        GooseFEM::Matrix K_a(conn_a, dofs);
        GooseFEM::Matrix K_b(conn_b, dofs);

        auto coor = mesh.coor();

        GooseFEM::Element::Quad4::QuadraturePlanar quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastic::Cartesian3d::Array<2> mat({nelem, nip});
        GMatElastic::Cartesian3d::Array<2> mat_a({nelem_a, nip}, 3.0, 4.0);
        GMatElastic::Cartesian3d::Array<2> mat_b({nelem_b, nip}, 5.0, 6.0);

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_a), xt::all()) = 1;
            mat.setElastic(I, 3.0, 4.0);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_b), xt::all()) = 1;
            mat.setElastic(I, 5.0, 6.0);
        }

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

            auto C = mat.Tangent();
            auto C_a = mat_a.Tangent();
            auto C_b = mat_b.Tangent();

            K.assemble(quad.Int_gradN_dot_tensor4_dot_gradNT_dV(C));
            K_a.assemble(quad_a.Int_gradN_dot_tensor4_dot_gradNT_dV(C_a));
            K_b.assemble(quad_b.Int_gradN_dot_tensor4_dot_gradNT_dV(C_b));

            auto f = vector.AssembleNode(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto f_a = vector_a.AssembleNode(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto f_b = vector_b.AssembleNode(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(f, f_a + f_b));
            REQUIRE(xt::allclose(f, K.Dot(disp)));
            REQUIRE(xt::allclose(f_a, K_a.Dot(disp)));
            REQUIRE(xt::allclose(f_b, K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_a + K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_b + K_a.Dot(disp)));

            auto fd = vector.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto fd_a = vector_a.AssembleDofs(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto fd_b = vector_b.AssembleDofs(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(fd, fd_a + fd_b));
            REQUIRE(xt::allclose(fd, K.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd_a, K_a.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd_b, K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_a + K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_b + K_a.Dot(vector.AsDofs(disp))));

            auto Sig_c = decltype(Sig)::from_shape(Sig.shape());
            xt::view(Sig_c, xt::keep(elem_a)) = Sig_a;
            xt::view(Sig_c, xt::keep(elem_b)) = Sig_b;

            REQUIRE(xt::allclose(Sig, Sig_c));
        }
    }

    SECTION("Vector/Matrix - GMatElastic - Periodic")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N, 2 * N);
        xt::xtensor<size_t, 1> elem_b =
            xt::concatenate(xt::xtuple(xt::arange<size_t>(N), xt::arange<size_t>(2 * N, nelem)));
        size_t nelem_a = elem_a.size();
        size_t nelem_b = elem_b.size();

        auto dofs = mesh.dofsPeriodic();
        auto conn = mesh.conn();
        xt::xtensor<size_t, 2> conn_a = xt::view(conn, xt::keep(elem_a), xt::all());
        xt::xtensor<size_t, 2> conn_b = xt::view(conn, xt::keep(elem_b), xt::all());

        GooseFEM::Vector vector(conn, dofs);
        GooseFEM::Vector vector_a(conn_a, dofs);
        GooseFEM::Vector vector_b(conn_b, dofs);

        GooseFEM::Matrix K(conn, dofs);
        GooseFEM::Matrix K_a(conn_a, dofs);
        GooseFEM::Matrix K_b(conn_b, dofs);

        auto coor = mesh.coor();

        GooseFEM::Element::Quad4::QuadraturePlanar quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastic::Cartesian3d::Array<2> mat({nelem, nip});
        GMatElastic::Cartesian3d::Array<2> mat_a({nelem_a, nip}, 3.0, 4.0);
        GMatElastic::Cartesian3d::Array<2> mat_b({nelem_b, nip}, 5.0, 6.0);

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_a), xt::all()) = 1;
            mat.setElastic(I, 3.0, 4.0);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_b), xt::all()) = 1;
            mat.setElastic(I, 5.0, 6.0);
        }

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

            auto C = mat.Tangent();
            auto C_a = mat_a.Tangent();
            auto C_b = mat_b.Tangent();

            K.assemble(quad.Int_gradN_dot_tensor4_dot_gradNT_dV(C));
            K_a.assemble(quad_a.Int_gradN_dot_tensor4_dot_gradNT_dV(C_a));
            K_b.assemble(quad_b.Int_gradN_dot_tensor4_dot_gradNT_dV(C_b));

            auto f = vector.AssembleNode(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto f_a = vector_a.AssembleNode(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto f_b = vector_b.AssembleNode(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(f, f_a + f_b));
            REQUIRE(xt::allclose(f, K.Dot(disp)));
            REQUIRE(xt::allclose(f_a, K_a.Dot(disp)));
            REQUIRE(xt::allclose(f_b, K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_a + K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_b + K_a.Dot(disp)));

            auto fd = vector.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto fd_a = vector_a.AssembleDofs(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto fd_b = vector_b.AssembleDofs(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(fd, fd_a + fd_b));
            REQUIRE(xt::allclose(fd, K.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd_a, K_a.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd_b, K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_a + K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_b + K_a.Dot(vector.AsDofs(disp))));

            auto Sig_c = decltype(Sig)::from_shape(Sig.shape());
            xt::view(Sig_c, xt::keep(elem_a)) = Sig_a;
            xt::view(Sig_c, xt::keep(elem_b)) = Sig_b;

            REQUIRE(xt::allclose(Sig, Sig_c));
        }
    }

    SECTION("VectorPartitioned/Matrix - GMatElastic - Semi-Periodic")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        size_t ndim = mesh.ndim();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N, 2 * N);
        xt::xtensor<size_t, 1> elem_b =
            xt::concatenate(xt::xtuple(xt::arange<size_t>(N), xt::arange<size_t>(2 * N, nelem)));
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

        GooseFEM::Matrix K(conn, dofs);
        GooseFEM::Matrix K_a(conn_a, dofs);
        GooseFEM::Matrix K_b(conn_b, dofs);

        auto coor = mesh.coor();

        GooseFEM::Element::Quad4::QuadraturePlanar quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastic::Cartesian3d::Array<2> mat({nelem, nip});
        GMatElastic::Cartesian3d::Array<2> mat_a({nelem_a, nip}, 3.0, 4.0);
        GMatElastic::Cartesian3d::Array<2> mat_b({nelem_b, nip}, 5.0, 6.0);

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_a), xt::all()) = 1;
            mat.setElastic(I, 3.0, 4.0);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, xt::keep(elem_b), xt::all()) = 1;
            mat.setElastic(I, 5.0, 6.0);
        }

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

            auto C = mat.Tangent();
            auto C_a = mat_a.Tangent();
            auto C_b = mat_b.Tangent();

            K.assemble(quad.Int_gradN_dot_tensor4_dot_gradNT_dV(C));
            K_a.assemble(quad_a.Int_gradN_dot_tensor4_dot_gradNT_dV(C_a));
            K_b.assemble(quad_b.Int_gradN_dot_tensor4_dot_gradNT_dV(C_b));

            auto f = vector.AssembleNode(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto f_a = vector_a.AssembleNode(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto f_b = vector_b.AssembleNode(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(f, f_a + f_b));
            REQUIRE(xt::allclose(f, K.Dot(disp)));
            REQUIRE(xt::allclose(f_a, K_a.Dot(disp)));
            REQUIRE(xt::allclose(f_b, K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_a + K_b.Dot(disp)));
            REQUIRE(xt::allclose(f, f_b + K_a.Dot(disp)));

            auto fd = vector.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(Sig));
            auto fd_a = vector_a.AssembleDofs(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
            auto fd_b = vector_b.AssembleDofs(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

            REQUIRE(xt::allclose(fd, fd_a + fd_b));
            REQUIRE(xt::allclose(fd, K.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd_a, K_a.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd_b, K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_a + K_b.Dot(vector.AsDofs(disp))));
            REQUIRE(xt::allclose(fd, fd_b + K_a.Dot(vector.AsDofs(disp))));

            auto Sig_c = decltype(Sig)::from_shape(Sig.shape());
            xt::view(Sig_c, xt::keep(elem_a)) = Sig_a;
            xt::view(Sig_c, xt::keep(elem_b)) = Sig_b;

            REQUIRE(xt::allclose(Sig, Sig_c));
        }
    }
}


#include "support.h"
#include <xtensor/xrandom.hpp>
#include <GMatElastic/Cartesian3d.h>

TEST_CASE("Example_hybrid-material", "GooseFEM.h")
{

    SECTION("Vector/Matrix - GMatElastic")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N);
        xt::xtensor<size_t, 1> elem_b = xt::arange<size_t>(N, nelem);
        size_t nelem_a = elem_a.size();
        size_t nelem_b = elem_b.size();

        xt::xtensor<size_t, 2> conn = mesh.conn();
        xt::xtensor<size_t, 2> conn_a = xt::view(conn, xt::keep(elem_a), xt::all());
        xt::xtensor<size_t, 2> conn_b = xt::view(conn, xt::keep(elem_b), xt::all());

        GooseFEM::Vector vector(conn, mesh.dofs());
        GooseFEM::Vector vector_a(conn_a, mesh.dofs());
        GooseFEM::Vector vector_b(conn_b, mesh.dofs());

        GooseFEM::Matrix K(conn, mesh.dofs());
        GooseFEM::Matrix K_a(conn_a, mesh.dofs());
        GooseFEM::Matrix K_b(conn_b, mesh.dofs());

        xt::xtensor<double, 2> coor = mesh.coor();
        xt::xtensor<double, 2> disp = xt::random::rand<double>(coor.shape());

        GooseFEM::Element::Quad4::QuadraturePlanar quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastic::Cartesian3d::Matrix mat(nelem, nip);
        GMatElastic::Cartesian3d::Matrix mat_a(nelem_a, nip, 3.0, 4.0);
        GMatElastic::Cartesian3d::Matrix mat_b(nelem_b, nip, 5.0, 6.0);

        xt::xtensor<size_t, 2> I = xt::empty<size_t>({nelem, nip});
        I.fill(0);
        xt::view(I, xt::keep(elem_a), xt::all()) = 1;
        mat.setElastic(I, 3.0, 4.0);
        I.fill(0);
        xt::view(I, xt::keep(elem_b), xt::all()) = 1;
        mat.setElastic(I, 5.0, 6.0);

        mat.check();
        mat_a.check();
        mat_b.check();

        xt::xtensor<double, 4> Eps = quad.SymGradN_vector(vector.AsElement(disp));
        xt::xtensor<double, 4> Eps_a = quad_a.SymGradN_vector(vector_a.AsElement(disp));
        xt::xtensor<double, 4> Eps_b = quad_b.SymGradN_vector(vector_b.AsElement(disp));
        xt::xtensor<double, 4> Sig = mat.Stress(Eps);
        xt::xtensor<double, 4> Sig_a = mat_a.Stress(Eps_a);
        xt::xtensor<double, 4> Sig_b = mat_b.Stress(Eps_b);
        xt::xtensor<double, 6> C;
        xt::xtensor<double, 6> C_a;
        xt::xtensor<double, 6> C_b;
        std::tie(Sig, C) = mat.Tangent(Eps);
        std::tie(Sig_a, C_a) = mat_a.Tangent(Eps_a);
        std::tie(Sig_b, C_b) = mat_b.Tangent(Eps_b);

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

        auto fd = vector.AssembleDofs(quad.Int_gradN_dot_tensor2_dV(Sig));
        auto fd_a = vector_a.AssembleDofs(quad_a.Int_gradN_dot_tensor2_dV(Sig_a));
        auto fd_b = vector_b.AssembleDofs(quad_b.Int_gradN_dot_tensor2_dV(Sig_b));

        REQUIRE(xt::allclose(fd, fd_a + fd_b));
        REQUIRE(xt::allclose(fd, K.Dot(vector.AsDofs(disp))));
        REQUIRE(xt::allclose(fd_a, K_a.Dot(vector.AsDofs(disp))));
        REQUIRE(xt::allclose(fd_b, K_b.Dot(vector.AsDofs(disp))));
    }

    SECTION("Vector/Matrix - GMatElastic - Periodic")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N);
        xt::xtensor<size_t, 1> elem_b = xt::arange<size_t>(N, nelem);
        size_t nelem_a = elem_a.size();
        size_t nelem_b = elem_b.size();

        xt::xtensor<size_t, 2> conn = mesh.conn();
        xt::xtensor<size_t, 2> conn_a = xt::view(conn, xt::keep(elem_a), xt::all());
        xt::xtensor<size_t, 2> conn_b = xt::view(conn, xt::keep(elem_b), xt::all());

        GooseFEM::Vector vector(conn, mesh.dofsPeriodic());
        GooseFEM::Vector vector_a(conn_a, mesh.dofsPeriodic());
        GooseFEM::Vector vector_b(conn_b, mesh.dofsPeriodic());

        GooseFEM::Matrix K(conn, mesh.dofsPeriodic());
        GooseFEM::Matrix K_a(conn_a, mesh.dofsPeriodic());
        GooseFEM::Matrix K_b(conn_b, mesh.dofsPeriodic());

        xt::xtensor<double, 2> coor = mesh.coor();
        xt::xtensor<double, 2> disp = xt::random::rand<double>(coor.shape());
        disp = vector.AsNode(vector.AsDofs(disp));

        GooseFEM::Element::Quad4::QuadraturePlanar quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastic::Cartesian3d::Matrix mat(nelem, nip);
        GMatElastic::Cartesian3d::Matrix mat_a(nelem_a, nip, 3.0, 4.0);
        GMatElastic::Cartesian3d::Matrix mat_b(nelem_b, nip, 5.0, 6.0);

        xt::xtensor<size_t, 2> I = xt::empty<size_t>({nelem, nip});
        I.fill(0);
        xt::view(I, xt::keep(elem_a), xt::all()) = 1;
        mat.setElastic(I, 3.0, 4.0);
        I.fill(0);
        xt::view(I, xt::keep(elem_b), xt::all()) = 1;
        mat.setElastic(I, 5.0, 6.0);

        mat.check();
        mat_a.check();
        mat_b.check();

        xt::xtensor<double, 4> Eps = quad.SymGradN_vector(vector.AsElement(disp));
        xt::xtensor<double, 4> Eps_a = quad_a.SymGradN_vector(vector_a.AsElement(disp));
        xt::xtensor<double, 4> Eps_b = quad_b.SymGradN_vector(vector_b.AsElement(disp));
        xt::xtensor<double, 4> Sig = mat.Stress(Eps);
        xt::xtensor<double, 4> Sig_a = mat_a.Stress(Eps_a);
        xt::xtensor<double, 4> Sig_b = mat_b.Stress(Eps_b);
        xt::xtensor<double, 6> C;
        xt::xtensor<double, 6> C_a;
        xt::xtensor<double, 6> C_b;
        std::tie(Sig, C) = mat.Tangent(Eps);
        std::tie(Sig_a, C_a) = mat_a.Tangent(Eps_a);
        std::tie(Sig_b, C_b) = mat_b.Tangent(Eps_b);

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
    }

    SECTION("VectorPartitioned/Matrix - GMatElastic - Periodic")
    {
        size_t N = 5;
        GooseFEM::Mesh::Quad4::Regular mesh(N, N);
        size_t nelem = mesh.nelem();
        size_t ndim = mesh.ndim();
        xt::xtensor<size_t, 1> elem_a = xt::arange<size_t>(N);
        xt::xtensor<size_t, 1> elem_b = xt::arange<size_t>(N, nelem);
        size_t nelem_a = elem_a.size();
        size_t nelem_b = elem_b.size();
        auto dofs = mesh.dofs();

        // periodicity in horizontal direction : eliminate 'dependent' DOFs
        auto left = mesh.nodesLeftOpenEdge();
        auto right = mesh.nodesRightOpenEdge();
        xt::view(dofs, xt::keep(right), 0) = xt::view(dofs, xt::keep(left), 0);
        xt::view(dofs, xt::keep(right), 1) = xt::view(dofs, xt::keep(left), 1);

        // fixed top and bottom
        auto top = mesh.nodesTopEdge();
        auto bottom = mesh.nodesBottomEdge();
        size_t nfix = top.size();
        xt::xtensor<size_t, 1> iip = xt::empty<size_t>({2 * ndim * nfix});
        xt::view(iip, xt::range(0 * nfix, 1 * nfix)) = xt::view(dofs, xt::keep(bottom), 0);
        xt::view(iip, xt::range(1 * nfix, 2 * nfix)) = xt::view(dofs, xt::keep(bottom), 1);
        xt::view(iip, xt::range(2 * nfix, 3 * nfix)) = xt::view(dofs, xt::keep(top), 0);
        xt::view(iip, xt::range(3 * nfix, 4 * nfix)) = xt::view(dofs, xt::keep(top), 1);

        xt::xtensor<size_t, 2> conn = mesh.conn();
        xt::xtensor<size_t, 2> conn_a = xt::view(conn, xt::keep(elem_a), xt::all());
        xt::xtensor<size_t, 2> conn_b = xt::view(conn, xt::keep(elem_b), xt::all());

        GooseFEM::VectorPartitioned vector(conn, dofs, iip);
        GooseFEM::VectorPartitioned vector_a(conn_a, dofs, iip);
        GooseFEM::VectorPartitioned vector_b(conn_b, dofs, iip);

        GooseFEM::Matrix K(conn, dofs);
        GooseFEM::Matrix K_a(conn_a, dofs);
        GooseFEM::Matrix K_b(conn_b, dofs);

        xt::xtensor<double, 2> coor = mesh.coor();
        xt::xtensor<double, 2> disp = xt::random::rand<double>(coor.shape());
        disp = vector.AsNode(vector.AsDofs(disp));

        GooseFEM::Element::Quad4::QuadraturePlanar quad(vector.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_a(vector_a.AsElement(coor));
        GooseFEM::Element::Quad4::QuadraturePlanar quad_b(vector_b.AsElement(coor));
        size_t nip = quad.nip();

        GMatElastic::Cartesian3d::Matrix mat(nelem, nip);
        GMatElastic::Cartesian3d::Matrix mat_a(nelem_a, nip, 3.0, 4.0);
        GMatElastic::Cartesian3d::Matrix mat_b(nelem_b, nip, 5.0, 6.0);

        xt::xtensor<size_t, 2> I = xt::empty<size_t>({nelem, nip});
        I.fill(0);
        xt::view(I, xt::keep(elem_a), xt::all()) = 1;
        mat.setElastic(I, 3.0, 4.0);
        I.fill(0);
        xt::view(I, xt::keep(elem_b), xt::all()) = 1;
        mat.setElastic(I, 5.0, 6.0);

        mat.check();
        mat_a.check();
        mat_b.check();

        xt::xtensor<double, 4> Eps = quad.SymGradN_vector(vector.AsElement(disp));
        xt::xtensor<double, 4> Eps_a = quad_a.SymGradN_vector(vector_a.AsElement(disp));
        xt::xtensor<double, 4> Eps_b = quad_b.SymGradN_vector(vector_b.AsElement(disp));
        xt::xtensor<double, 4> Sig = mat.Stress(Eps);
        xt::xtensor<double, 4> Sig_a = mat_a.Stress(Eps_a);
        xt::xtensor<double, 4> Sig_b = mat_b.Stress(Eps_b);
        xt::xtensor<double, 6> C;
        xt::xtensor<double, 6> C_a;
        xt::xtensor<double, 6> C_b;
        std::tie(Sig, C) = mat.Tangent(Eps);
        std::tie(Sig_a, C_a) = mat_a.Tangent(Eps_a);
        std::tie(Sig_b, C_b) = mat_b.Tangent(Eps_b);

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
    }
}

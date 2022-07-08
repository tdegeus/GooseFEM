#include <GooseFEM/GooseFEM.h>
#include <catch2/catch_all.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xrandom.hpp>

#define ISCLOSE(a, b) REQUIRE_THAT((a), Catch::Matchers::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Vector", "Vector.h")
{

    SECTION("asDofs - nodevec")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

        xt::xtensor<double, 2> v = xt::empty<double>({mesh.nnode(), std::size_t(2)});
        v(0, 0) = 1.0;
        v(0, 1) = 0.0;
        v(1, 0) = 1.0;
        v(1, 1) = 0.0;
        v(2, 0) = 1.0;
        v(2, 1) = 0.0;
        v(3, 0) = 1.5;
        v(3, 1) = 0.0;
        v(4, 0) = 1.5;
        v(4, 1) = 0.0;
        v(5, 0) = 1.5;
        v(5, 1) = 0.0;
        v(6, 0) = 1.0;
        v(6, 1) = 0.0;
        v(7, 0) = 1.0;
        v(7, 1) = 0.0;
        v(8, 0) = 1.0;
        v(8, 1) = 0.0;

        auto V = vector.AsDofs(v);

        REQUIRE(xt::has_shape(vector.allocate_dofval(), V.shape()));
        REQUIRE(xt::has_shape(vector.allocate_nodevec(), v.shape()));
        REQUIRE(xt::has_shape(vector.allocate_dofval(), vector.AsDofs(v).shape()));
        REQUIRE(xt::has_shape(vector.allocate_elemvec(), vector.AsElement(v).shape()));
        REQUIRE(xt::has_shape(vector.allocate_nodevec(), vector.AsNode(V).shape()));
        REQUIRE(xt::has_shape(vector.allocate_elemvec(), vector.AsElement(V).shape()));
        REQUIRE(V.size() == (mesh.nnode() - mesh.nodesPeriodic().shape(0)) * mesh.ndim());
        ISCLOSE(V(0), v(0, 0));
        ISCLOSE(V(1), v(0, 1));
        ISCLOSE(V(2), v(1, 0));
        ISCLOSE(V(3), v(1, 1));
        ISCLOSE(V(4), v(3, 0));
        ISCLOSE(V(5), v(3, 1));
        ISCLOSE(V(6), v(4, 0));
        ISCLOSE(V(7), v(4, 1));
    }

    SECTION("asDofs - elemvec")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

        xt::xtensor<double, 2> v = xt::empty<double>({mesh.nnode(), std::size_t(2)});
        v(0, 0) = 1.0;
        v(0, 1) = 0.0;
        v(1, 0) = 1.0;
        v(1, 1) = 0.0;
        v(2, 0) = 1.0;
        v(2, 1) = 0.0;
        v(3, 0) = 1.5;
        v(3, 1) = 0.0;
        v(4, 0) = 1.5;
        v(4, 1) = 0.0;
        v(5, 0) = 1.5;
        v(5, 1) = 0.0;
        v(6, 0) = 1.0;
        v(6, 1) = 0.0;
        v(7, 0) = 1.0;
        v(7, 1) = 0.0;
        v(8, 0) = 1.0;
        v(8, 1) = 0.0;

        auto V = vector.AsDofs(vector.AsElement(vector.AsDofs(v)));

        REQUIRE(V.size() == (mesh.nnode() - mesh.nodesPeriodic().shape(0)) * mesh.ndim());
        ISCLOSE(V(0), v(0, 0));
        ISCLOSE(V(1), v(0, 1));
        ISCLOSE(V(2), v(1, 0));
        ISCLOSE(V(3), v(1, 1));
        ISCLOSE(V(4), v(3, 0));
        ISCLOSE(V(5), v(3, 1));
        ISCLOSE(V(6), v(4, 0));
        ISCLOSE(V(7), v(4, 1));
    }

    SECTION("asDofs - assembleDofs")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

        xt::xtensor<double, 2> f = xt::empty<double>({mesh.nnode(), std::size_t(2)});
        f(0, 0) = -1.0;
        f(0, 1) = -1.0;
        f(1, 0) = 0.0;
        f(1, 1) = -1.0;
        f(2, 0) = 1.0;
        f(2, 1) = -1.0;
        f(3, 0) = -1.0;
        f(3, 1) = 0.0;
        f(4, 0) = 0.0;
        f(4, 1) = 0.0;
        f(5, 0) = 1.0;
        f(5, 1) = 0.0;
        f(6, 0) = -1.0;
        f(6, 1) = 1.0;
        f(7, 0) = 0.0;
        f(7, 1) = 1.0;
        f(8, 0) = 1.0;
        f(8, 1) = 1.0;

        xt::xtensor<double, 1> F = vector.AssembleDofs(f);

        REQUIRE(F.size() == (mesh.nnode() - mesh.nodesPeriodic().shape(0)) * mesh.ndim());
        ISCLOSE(F(0), 0);
        ISCLOSE(F(1), 0);
        ISCLOSE(F(2), 0);
        ISCLOSE(F(3), 0);
        ISCLOSE(F(4), 0);
        ISCLOSE(F(5), 0);
        ISCLOSE(F(6), 0);
        ISCLOSE(F(7), 0);
    }

    SECTION("asDofs - assembleNode")
    {
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

        xt::xtensor<double, 2> f = xt::empty<double>({mesh.nnode(), std::size_t(2)});
        f(0, 0) = -1.0;
        f(0, 1) = -1.0;
        f(1, 0) = 0.0;
        f(1, 1) = -1.0;
        f(2, 0) = 1.0;
        f(2, 1) = -1.0;
        f(3, 0) = -1.0;
        f(3, 1) = 0.0;
        f(4, 0) = 0.0;
        f(4, 1) = 0.0;
        f(5, 0) = 1.0;
        f(5, 1) = 0.0;
        f(6, 0) = -1.0;
        f(6, 1) = 1.0;
        f(7, 0) = 0.0;
        f(7, 1) = 1.0;
        f(8, 0) = 1.0;
        f(8, 1) = 1.0;

        xt::xtensor<double, 1> F = vector.AssembleDofs(vector.AsElement(f));

        REQUIRE(F.size() == (mesh.nnode() - mesh.nodesPeriodic().shape(0)) * mesh.ndim());
        ISCLOSE(F(0), 0);
        ISCLOSE(F(1), 0);
        ISCLOSE(F(2), 0);
        ISCLOSE(F(3), 0);
        ISCLOSE(F(4), 0);
        ISCLOSE(F(5), 0);
        ISCLOSE(F(6), 0);
        ISCLOSE(F(7), 0);
    }
}

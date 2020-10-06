
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor/xmath.hpp>
#include <GooseFEM/GooseFEM.h>

#define ISCLOSE(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

TEST_CASE("GooseFEM::Vector", "Vector.h")
{

    SECTION("asDofs - nodevec")
    {
        // mesh
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        // vector-definition
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

        // velocity field
        // - allocate
        xt::xtensor<double, 2> v = xt::empty<double>({mesh.nnode(), std::size_t(2)});
        // - set periodic
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

        // convert to DOFs
        xt::xtensor<double, 1> V = vector.AsDofs(v);

        // check
        // - size
        REQUIRE(V.size() == (mesh.nnode() - mesh.nodesPeriodic().shape(0)) * mesh.ndim());
        // - individual entries
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
        // mesh
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        // vector-definition
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

        // velocity field
        // - allocate
        xt::xtensor<double, 2> v = xt::empty<double>({mesh.nnode(), std::size_t(2)});
        // - set periodic
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

        // convert to DOFs - element - DOFs
        xt::xtensor<double, 1> V = vector.AsDofs(vector.AsElement(vector.AsDofs(v)));

        // check
        // - size
        REQUIRE(V.size() == (mesh.nnode() - mesh.nodesPeriodic().shape(0)) * mesh.ndim());
        // - individual entries
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
        // mesh
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        // vector-definition
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

        // force field
        // - allocate
        xt::xtensor<double, 2> f = xt::empty<double>({mesh.nnode(), std::size_t(2)});
        // - set periodic
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

        // assemble as DOFs
        xt::xtensor<double, 1> F = vector.AssembleDofs(f);

        // check
        // - size
        REQUIRE(F.size() == (mesh.nnode() - mesh.nodesPeriodic().shape(0)) * mesh.ndim());
        // - 'analytical' result
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
        // mesh
        GooseFEM::Mesh::Quad4::Regular mesh(2, 2);

        // vector-definition
        GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

        // force field
        // - allocate
        xt::xtensor<double, 2> f = xt::empty<double>({mesh.nnode(), std::size_t(2)});
        // - set periodic
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

        // convert to element, assemble as DOFs
        xt::xtensor<double, 1> F = vector.AssembleDofs(vector.AsElement(f));

        // check
        // - size
        REQUIRE(F.size() == (mesh.nnode() - mesh.nodesPeriodic().shape(0)) * mesh.ndim());
        // - 'analytical' result
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

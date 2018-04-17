
#include <catch/catch.hpp>

#define CPPMAT_NOCONVERT
#include <cppmat/cppmat.h>

#include <Eigen/Eigen>

#include <GooseFEM/GooseFEM.h>

// =================================================================================================

TEST_CASE("GooseFEM::Vector", "Vector.h")
{

// =================================================================================================

SECTION( "asDofs - nodevec" )
{
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);
  GooseFEM::MatD v(mesh.nnode(),2);

  v(0,0) = 1.0;  v(0,1) = 0.0;
  v(1,0) = 1.0;  v(1,1) = 0.0;
  v(2,0) = 1.0;  v(2,1) = 0.0;
  v(3,0) = 1.5;  v(3,1) = 0.0;
  v(4,0) = 1.5;  v(4,1) = 0.0;
  v(5,0) = 1.5;  v(5,1) = 0.0;
  v(6,0) = 1.0;  v(6,1) = 0.0;
  v(7,0) = 1.0;  v(7,1) = 0.0;
  v(8,0) = 1.0;  v(8,1) = 0.0;

  GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

  GooseFEM::ColD V = vector.asDofs(v);

  REQUIRE( V.size() == 4*2 );

  REQUIRE( V(0) == v(0,0) );
  REQUIRE( V(1) == v(0,1) );
  REQUIRE( V(2) == v(1,0) );
  REQUIRE( V(3) == v(1,1) );
  REQUIRE( V(4) == v(3,0) );
  REQUIRE( V(5) == v(3,1) );
  REQUIRE( V(6) == v(4,0) );
  REQUIRE( V(7) == v(4,1) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "asDofs - elemvec" )
{
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);
  GooseFEM::MatD v(mesh.nnode(),2);

  v(0,0) = 1.0;  v(0,1) = 0.0;
  v(1,0) = 1.0;  v(1,1) = 0.0;
  v(2,0) = 1.0;  v(2,1) = 0.0;
  v(3,0) = 1.5;  v(3,1) = 0.0;
  v(4,0) = 1.5;  v(4,1) = 0.0;
  v(5,0) = 1.5;  v(5,1) = 0.0;
  v(6,0) = 1.0;  v(6,1) = 0.0;
  v(7,0) = 1.0;  v(7,1) = 0.0;
  v(8,0) = 1.0;  v(8,1) = 0.0;

  GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

  GooseFEM::ColD V = vector.asDofs(v);

  V = vector.asDofs(vector.asElement(V));

  REQUIRE( V.size() == 4*2 );

  REQUIRE( V(0) == v(0,0) );
  REQUIRE( V(1) == v(0,1) );
  REQUIRE( V(2) == v(1,0) );
  REQUIRE( V(3) == v(1,1) );
  REQUIRE( V(4) == v(3,0) );
  REQUIRE( V(5) == v(3,1) );
  REQUIRE( V(6) == v(4,0) );
  REQUIRE( V(7) == v(4,1) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "asDofs - assembleDofs" )
{
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);
  GooseFEM::MatD f(mesh.nnode(),2);

  f(0,0) = -1.0;  f(0,1) = -1.0;
  f(1,0) =  0.0;  f(1,1) = -1.0;
  f(2,0) =  1.0;  f(2,1) = -1.0;
  f(3,0) = -1.0;  f(3,1) =  0.0;
  f(4,0) =  0.0;  f(4,1) =  0.0;
  f(5,0) =  1.0;  f(5,1) =  0.0;
  f(6,0) = -1.0;  f(6,1) =  1.0;
  f(7,0) =  0.0;  f(7,1) =  1.0;
  f(8,0) =  1.0;  f(8,1) =  1.0;

  GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

  GooseFEM::ColD F = vector.assembleDofs(f);

  REQUIRE( F.size() == 4*2 );

  REQUIRE( F(0) == 0.0 );
  REQUIRE( F(1) == 0.0 );
  REQUIRE( F(2) == 0.0 );
  REQUIRE( F(3) == 0.0 );
  REQUIRE( F(4) == 0.0 );
  REQUIRE( F(5) == 0.0 );
  REQUIRE( F(6) == 0.0 );
  REQUIRE( F(7) == 0.0 );
}

// -------------------------------------------------------------------------------------------------

SECTION( "asDofs - assembleNode" )
{
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);
  GooseFEM::MatD f(mesh.nnode(),2);

  f(0,0) = -1.0;  f(0,1) = -1.0;
  f(1,0) =  0.0;  f(1,1) = -1.0;
  f(2,0) =  1.0;  f(2,1) = -1.0;
  f(3,0) = -1.0;  f(3,1) =  0.0;
  f(4,0) =  0.0;  f(4,1) =  0.0;
  f(5,0) =  1.0;  f(5,1) =  0.0;
  f(6,0) = -1.0;  f(6,1) =  1.0;
  f(7,0) =  0.0;  f(7,1) =  1.0;
  f(8,0) =  1.0;  f(8,1) =  1.0;

  GooseFEM::Vector vector(mesh.conn(), mesh.dofsPeriodic());

  GooseFEM::ColD F = vector.assembleDofs( vector.asElement(f) );

  REQUIRE( F.size() == 4*2 );

  REQUIRE( F(0) == 0.0 );
  REQUIRE( F(1) == 0.0 );
  REQUIRE( F(2) == 0.0 );
  REQUIRE( F(3) == 0.0 );
  REQUIRE( F(4) == 0.0 );
  REQUIRE( F(5) == 0.0 );
  REQUIRE( F(6) == 0.0 );
  REQUIRE( F(7) == 0.0 );
}

// -------------------------------------------------------------------------------------------------

// =================================================================================================

}

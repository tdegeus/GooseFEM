
#include <catch/catch.hpp>

#define CPPMAT_NOCONVERT
#include <cppmat/cppmat.h>

#include <Eigen/Eigen>

#include <GooseFEM/GooseFEM.h>

// =================================================================================================

TEST_CASE("GooseFEM::MatrixDiagonal", "MatrixDiagonal.h")
{

// =================================================================================================

SECTION( "dot" )
{
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);

  GooseFEM::MatD a = GooseFEM::MatD::Random(9*2,9*2);
  GooseFEM::ColD b = GooseFEM::ColD::Random(9*2);
  GooseFEM::ColD c;

  for ( auto i = 0 ; i < a.rows() ; ++i ) {
    for ( auto j = 0 ; j < a.cols() ; ++j ) {
      if ( i != j ) a(i,j)  = 0.0;
      else          a(i,j) += 1.0;
    }
  }

  c = a * b;

  GooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
  GooseFEM::ColD C;

  for ( auto i = 0 ; i < a.rows() ; ++i )
    A(i,i) = a(i,i);

  C = A * b;

  REQUIRE( C.size() == c.size() );

  for ( auto i = 0 ; i < c.rows() ; ++i )
    REQUIRE( C(i) == c(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "solve" )
{
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);

  GooseFEM::MatD a = GooseFEM::MatD::Random(9*2,9*2);
  GooseFEM::ColD b = GooseFEM::ColD::Random(9*2);
  GooseFEM::ColD c;

  for ( auto i = 0 ; i < a.rows() ; ++i ) {
    for ( auto j = 0 ; j < a.cols() ; ++j ) {
      if ( i != j ) a(i,j)  = 0.0;
      else          a(i,j) += 1.0;
    }
  }

  c = a * b;

  GooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
  GooseFEM::ColD B, C;

  for ( auto i = 0 ; i < a.rows() ; ++i )
    A(i,i) = a(i,i);

  C = A * b;

  B = A.solve(C);

  REQUIRE( B.size() == b.size() );

  for ( auto i = 0 ; i < b.rows() ; ++i )
    REQUIRE( std::abs(B(i)-b(i)) < 1.e-12 );
}

// =================================================================================================

}

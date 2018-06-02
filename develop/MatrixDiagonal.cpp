
#include "support.h"

// =================================================================================================

TEST_CASE("GooseFEM::MatrixDiagonal", "MatrixDiagonal.h")
{

// =================================================================================================

SECTION( "dot" )
{
  // mesh
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);

  // random matrix and column
  GooseFEM::MatD a = GooseFEM::MatD::Random(mesh.nnode()*mesh.ndim(),mesh.nnode()*mesh.ndim());
  GooseFEM::ColD b = GooseFEM::ColD::Random(mesh.nnode()*mesh.ndim());
  GooseFEM::ColD c;

  // set diagonal
  for ( auto i = 0 ; i < a.rows() ; ++i ) {
    for ( auto j = 0 ; j < a.cols() ; ++j ) {
      if ( i != j ) a(i,j)  = 0.0;
      else          a(i,j) += 1.0;
    }
  }

  // compute product using Eigen
  c = a * b;

  // convert to GooseFEM
  // - allocate
  GooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
  GooseFEM::ColD C;
  // - set
  for ( auto i = 0 ; i < a.rows() ; ++i )
    A(i,i) = a(i,i);

  // compute product
  C = A * b;

  // check
  // - size
  REQUIRE( C.size() == c.size() );
  // - components
  for ( auto i = 0 ; i < c.rows() ; ++i )
    EQ( C(i), c(i) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "solve" )
{
  // mesh
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);

  // random matrix and column
  GooseFEM::MatD a = GooseFEM::MatD::Random(mesh.nnode()*mesh.ndim(),mesh.nnode()*mesh.ndim());
  GooseFEM::ColD b = GooseFEM::ColD::Random(mesh.nnode()*mesh.ndim());
  GooseFEM::ColD c;

  // set diagonal
  for ( auto i = 0 ; i < a.rows() ; ++i ) {
    for ( auto j = 0 ; j < a.cols() ; ++j ) {
      if ( i != j ) a(i,j)  = 0.0;
      else          a(i,j) += 1.0;
    }
  }

  // compute product using Eigen
  c = a * b;

  // convert to GooseFEM
  // - allocate
  GooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
  GooseFEM::ColD B, C;
  // - set
  for ( auto i = 0 ; i < a.rows() ; ++i )
    A(i,i) = a(i,i);

  // compute product
  C = A * b;

  // solve
  B = A.solve(C);

  // check
  // - size
  REQUIRE( B.size() == b.size() );
  // - components
  for ( auto i = 0 ; i < b.rows() ; ++i )
    EQ( B(i), b(i) );
}

// =================================================================================================

}

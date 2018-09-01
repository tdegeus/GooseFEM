
#include "support.h"

// =================================================================================================

TEST_CASE("xGooseFEM::MatrixDiagonal", "MatrixDiagonal.h")
{

// =================================================================================================

SECTION( "dot" )
{
  // mesh
  xGooseFEM::Mesh::Quad4::Regular mesh(2,2);

  // random matrix and column
  xt::xtensor<double,1> a = xt::random::rand<double>({mesh.nnode()*mesh.ndim()});
  xt::xtensor<double,1> b = xt::random::rand<double>({mesh.nnode()*mesh.ndim()});
  xt::xtensor<double,1> c;

  // compute product
  c = a * b;

  // convert to xGooseFEM
  // - allocate
  xGooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
  xt::xtensor<double,1> C;
  // - set
  for ( size_t i = 0 ; i < a.shape()[0] ; ++i )
    A(i,i) = a(i);

  // compute product
  C = A.dot(b);

  // check
  // - size
  REQUIRE( C.size() == c.size() );
  // - components
  REQUIRE( xt::allclose(C, c) );
}

// -------------------------------------------------------------------------------------------------

SECTION( "solve" )
{
  // mesh
  xGooseFEM::Mesh::Quad4::Regular mesh(2,2);

  // random matrix and column
  xt::xtensor<double,1> a = xt::random::rand<double>({mesh.nnode()*mesh.ndim()});
  xt::xtensor<double,1> b = xt::random::rand<double>({mesh.nnode()*mesh.ndim()});
  xt::xtensor<double,1> c;

  // compute product
  c = a * b;

  // convert to xGooseFEM
  // - allocate
  xGooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
  xt::xtensor<double,1> B, C;
  // - set
  for ( size_t i = 0 ; i < a.shape()[0] ; ++i )
    A(i,i) = a(i);

  // compute product
  C = A.dot(b);

  // solve
  B = A.solve(C);

  // check
  // - size
  REQUIRE( B.size() == b.size() );
  // - components
  REQUIRE( xt::allclose(B, b) );
}

// // =================================================================================================

}

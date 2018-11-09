
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
  xt::xtensor<double,1> a = xt::random::rand<double>({mesh.nnode()*mesh.ndim()});
  xt::xtensor<double,1> b = xt::random::rand<double>({mesh.nnode()*mesh.ndim()});
  xt::xtensor<double,1> c;

  // compute product
  c = a * b;

  // convert to GooseFEM
  // - allocate
  GooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
  xt::xtensor<double,1> C;
  // - set
  A.set(a);

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
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);

  // random matrix and column
  xt::xtensor<double,1> a = xt::random::rand<double>({mesh.nnode()*mesh.ndim()});
  xt::xtensor<double,1> b = xt::random::rand<double>({mesh.nnode()*mesh.ndim()});
  xt::xtensor<double,1> c;

  // compute product
  c = a * b;

  // convert to GooseFEM
  // - allocate
  GooseFEM::MatrixDiagonal A(mesh.conn(), mesh.dofs());
  xt::xtensor<double,1> B, C;
  // - set
  A.set(a);

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

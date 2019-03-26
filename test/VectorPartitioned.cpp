
#include "support.h"

// =================================================================================================

TEST_CASE("GooseFEM::VectorPartitioned", "VectorPartitioned.h")
{

// =================================================================================================

SECTION("asDofs")
{
  // mesh

  GooseFEM::Mesh::Quad4::Regular mesh(9,9);

  // prescribed DOFs

  xt::xtensor<size_t,2> dofs = mesh.dofs();

  xt::xtensor<size_t,1> nodesLeft   = mesh.nodesLeftOpenEdge();
  xt::xtensor<size_t,1> nodesRight  = mesh.nodesRightOpenEdge();
  xt::xtensor<size_t,1> nodesTop    = mesh.nodesTopEdge();
  xt::xtensor<size_t,1> nodesBottom = mesh.nodesBottomEdge();

  for ( size_t i = 0 ; i < nodesLeft.size() ; ++i )
    for ( size_t j = 0 ; j < mesh.ndim() ; ++j )
      dofs(nodesRight(i),j) = dofs(nodesLeft(i),j);

  xt::xtensor<size_t,1> iip = xt::empty<size_t>({4*nodesBottom.size()});

  size_t i = 0;
  for ( auto &n : nodesTop    ) { iip(i) = dofs(n,0); ++i; }
  for ( auto &n : nodesTop    ) { iip(i) = dofs(n,1); ++i; }
  for ( auto &n : nodesBottom ) { iip(i) = dofs(n,0); ++i; }
  for ( auto &n : nodesBottom ) { iip(i) = dofs(n,1); ++i; }

  // random displacement

  xt::xtensor<double,2> u = xt::random::randn<double>({mesh.nnode(), mesh.ndim()});

  for ( size_t i = 0 ; i < nodesLeft.size() ; ++i )
    for ( size_t j = 0 ; j < mesh.ndim() ; ++j )
      u(nodesRight(i),j) = u(nodesLeft(i),j);

  // vector-definition

  GooseFEM::VectorPartitioned vector(mesh.conn(), dofs, iip);

  // convert

  xt::xtensor<double,1> u_u = vector.AsDofs_u(u);
  xt::xtensor<double,1> u_p = vector.AsDofs_p(u);

  REQUIRE( xt::allclose(u, vector.AsNode(u_u, u_p)) );

}

// =================================================================================================

SECTION("copy_u, copy_p")
{
  // mesh

  GooseFEM::Mesh::Quad4::Regular mesh(9,9);

  // prescribed DOFs

  xt::xtensor<size_t,2> dofs = mesh.dofs();

  xt::xtensor<size_t,1> nodesLeft   = mesh.nodesLeftOpenEdge();
  xt::xtensor<size_t,1> nodesRight  = mesh.nodesRightOpenEdge();
  xt::xtensor<size_t,1> nodesTop    = mesh.nodesTopEdge();
  xt::xtensor<size_t,1> nodesBottom = mesh.nodesBottomEdge();

  for ( size_t i = 0 ; i < nodesLeft.size() ; ++i )
    for ( size_t j = 0 ; j < mesh.ndim() ; ++j )
      dofs(nodesRight(i),j) = dofs(nodesLeft(i),j);

  xt::xtensor<size_t,1> iip = xt::empty<size_t>({4*nodesBottom.size()});

  size_t i = 0;
  for ( auto &n : nodesTop    ) { iip(i) = dofs(n,0); ++i; }
  for ( auto &n : nodesTop    ) { iip(i) = dofs(n,1); ++i; }
  for ( auto &n : nodesBottom ) { iip(i) = dofs(n,0); ++i; }
  for ( auto &n : nodesBottom ) { iip(i) = dofs(n,1); ++i; }

  // random displacement

  xt::xtensor<double,2> u = xt::random::randn<double>({mesh.nnode(), mesh.ndim()});

  for ( size_t i = 0 ; i < nodesLeft.size() ; ++i )
    for ( size_t j = 0 ; j < mesh.ndim() ; ++j )
      u(nodesRight(i),j) = u(nodesLeft(i),j);

  // vector-definition

  GooseFEM::VectorPartitioned vector(mesh.conn(), dofs, iip);

  // convert

  xt::xtensor<double,2> v = xt::empty<double>({mesh.nnode(), mesh.ndim()});

  vector.copy_u(u, v);
  vector.copy_p(u, v);

  REQUIRE( xt::allclose(u, v) );

}

// =================================================================================================

}

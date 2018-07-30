
#include "support.h"

// =================================================================================================

TEST_CASE("GooseFEM::ElementQuad4", "ElementQuad4.h")
{

using T2 = xt::xtensor_fixed<double, xt::xshape<2,2>>;

// =================================================================================================

SECTION( "int_N_scalar_NT_dV" )
{
  // mesh
  GooseFEM::Mesh::Quad4::Regular mesh(3,3);

  // vector-definition, and a diagonal matrix
  GooseFEM::Vector         vec(mesh.conn(), mesh.dofsPeriodic());
  GooseFEM::MatrixDiagonal mat(mesh.conn(), mesh.dofsPeriodic());

  // element definition, with nodal quadrature
  GooseFEM::Element::Quad4::Quadrature quad(
    vec.asElement(mesh.coor()),
    GooseFEM::Element::Quad4::Nodal::xi(),
    GooseFEM::Element::Quad4::Nodal::w()
  );

  // scalar per quadrature point (e.g. mass-density "rho")
  xt::xtensor<double,2> rho = xt::ones<double>({mesh.nelem(), quad.nip()});

  // evaluate integral and assemble diagonal matrix (e.g. mass matrix)
  mat.assemble(quad.int_N_scalar_NT_dV(rho));

  // check matrix
  // - get the matrix
  xt::xtensor<double,1> M = mat.asDiagonal();
  // - check the size
  REQUIRE( M.size() == vec.ndof() );
  // - check each component
  for ( size_t i = 0 ; i < M.size() ; ++i )
    EQ( M(i), 1 );
}

// =================================================================================================

SECTION( "symGradN_vector" )
{
  // mesh
  GooseFEM::Mesh::Quad4::FineLayer mesh(27,27);

  // vector-definition
  GooseFEM::Vector vec(mesh.conn(), mesh.dofs());

  // element definition, with Gauss quadrature
  GooseFEM::Element::Quad4::Quadrature quad( vec.asElement(mesh.coor()) );

  // macroscopic deformation gradient and strain
  // - zero-initialize
  T2 F   = xt::zeros<double>({2,2});
  T2 EPS = xt::zeros<double>({2,2});
  // - set non-zero components
  F  (0,1) = 0.1;
  EPS(0,1) = 0.05;
  EPS(1,0) = 0.05;

  // nodal coordinates and displacement
  xt::xtensor<double,2> coor = mesh.coor();;
  xt::xtensor<double,2> disp = xt::zeros<double>(coor.shape());

  // apply macroscopic deformation gradient
  for ( size_t n = 0 ; n < mesh.nnode() ; ++n )
    for ( size_t i = 0 ; i < F.shape()[0] ; ++i )
      for ( size_t j = 0 ; j < F.shape()[1] ; ++j )
        disp(n,i) += F(i,j) * coor(n,j);

  // compute quadrature point tensors
  xt::xtensor<double,4> eps = quad.symGradN_vector(vec.asElement(disp));

  // integration point volume
  xt::xtensor<double,4> dV = eps;
  quad.dV(dV);

  auto epsbar = xt::sum(dV*eps, {0,1}) / xt::sum(dV, {0,1});

  // check
  // - check sizes
  REQUIRE( eps.shape()[0] == mesh.nelem() );
  REQUIRE( eps.shape()[1] == quad.nip()   );
  REQUIRE( eps.shape()[2] == mesh.ndim()   );
  REQUIRE( eps.shape()[3] == mesh.ndim()   );
  // - check all components
  for ( size_t e = 0 ; e < mesh.nelem() ; ++e ) {
    for ( size_t k = 0 ; k < quad.nip() ; ++k ) {
      auto Eps = xt::view(eps, e, k, xt::all(), xt::all());
      for ( size_t i = 0 ; i < Eps.shape()[0] ; ++i )
        for ( size_t j = 0 ; j < Eps.shape()[1] ; ++j )
          EQ( Eps(i,j), EPS(i,j) );
    }
  }

  // check macroscopic tensor
  for ( size_t i = 0 ; i < epsbar.shape()[0] ; ++i )
    for ( size_t j = 0 ; j < epsbar.shape()[1] ; ++j )
      EQ( epsbar(i,j), EPS(i,j) );
}

// =================================================================================================

SECTION( "symGradN_vector, int_gradN_dot_tensor2s_dV" )
{
  // mesh
  GooseFEM::Mesh::Quad4::FineLayer mesh(27,27);

  // vector-definition
  GooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());

  // element definition, with Gauss quadrature
  GooseFEM::Element::Quad4::Quadrature quad( vec.asElement(mesh.coor()) );

  // macroscopic deformation gradient and strain
  // - zero-initialize
  T2 F = xt::zeros<double>({2,2});
  // - set non-zero components
  F(0,1) = 0.1;

  // nodal coordinates and displacement
  xt::xtensor<double,2> coor = mesh.coor();;
  xt::xtensor<double,2> disp = xt::zeros<double>(coor.shape());

  // apply macroscopic deformation gradient
  for ( size_t n = 0 ; n < mesh.nnode() ; ++n )
    for ( size_t i = 0 ; i < F.shape()[0] ; ++i )
      for ( size_t j = 0 ; j < F.shape()[1] ; ++j )
        disp(n,i) += F(i,j) * coor(n,j);

  // compute quadrature point tensors
  xt::xtensor<double,4> eps = quad.symGradN_vector(vec.asElement(disp));

  // nodal force vector (should be zero, as it is only sensitive to periodic fluctuations)
  xt::xtensor<double,1> Fi = vec.assembleDofs(quad.int_gradN_dot_tensor2_dV(eps));

  // check
  // - size
  REQUIRE( Fi.size() == vec.ndof() );
  // - check all components
  for ( size_t i = 0 ; i < vec.ndof() ; ++i )
    EQ( Fi(i), 0 );

}

// =================================================================================================

}

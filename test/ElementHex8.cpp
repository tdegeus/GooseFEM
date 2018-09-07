
#include "support.h"

// =================================================================================================

TEST_CASE("xGooseFEM::ElementHex8", "ElementHex8.h")
{

using T2 = xt::xtensor_fixed<double, xt::xshape<3,3>>;

// =================================================================================================

SECTION( "int_N_scalar_NT_dV" )
{
  // mesh
  xGooseFEM::Mesh::Hex8::Regular mesh(3,3,3);

  // vector-definition, and a diagonal matrix
  xGooseFEM::Vector         vec(mesh.conn(), mesh.dofsPeriodic());
  xGooseFEM::MatrixDiagonal mat(mesh.conn(), mesh.dofsPeriodic());

  // element definition, with nodal quadrature
  xGooseFEM::Element::Hex8::Quadrature quad(
    vec.asElement(mesh.coor()),
    xGooseFEM::Element::Hex8::Nodal::xi(),
    xGooseFEM::Element::Hex8::Nodal::w()
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
  REQUIRE( xt::allclose(M, 1.) );
}

// =================================================================================================

SECTION( "symGradN_vector" )
{
  // mesh
  xGooseFEM::Mesh::Hex8::FineLayer mesh(27,27,27);

  // vector-definition
  xGooseFEM::Vector vec(mesh.conn(), mesh.dofs());

  // element definition, with Gauss quadrature
  xGooseFEM::Element::Hex8::Quadrature quad( vec.asElement(mesh.coor()) );

  // macroscopic deformation gradient and strain
  // - zero-initialize
  T2 F   = xt::zeros<double>({3,3});
  T2 EPS = xt::zeros<double>({3,3});
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

  // volume averaged strain tensor
  auto epsbar = xt::average(eps, dV, {0,1});

  // check local strain tensors
  // - check sizes
  REQUIRE( eps.shape()[0] == mesh.nelem() );
  REQUIRE( eps.shape()[1] == quad.nip()   );
  REQUIRE( eps.shape()[2] == mesh.ndim()  );
  REQUIRE( eps.shape()[3] == mesh.ndim()  );
  // - check all components
  for ( size_t e = 0 ; e < mesh.nelem() ; ++e ) {
    for ( size_t k = 0 ; k < quad.nip() ; ++k ) {
      auto Eps = xt::view(eps, e, k);
      REQUIRE( xt::allclose(Eps, EPS));
    }
  }

  // check macroscopic tensor
  REQUIRE( xt::allclose(epsbar, EPS));
}

// =================================================================================================

SECTION( "symGradN_vector, int_gradN_dot_tensor2s_dV" )
{
  // mesh
  xGooseFEM::Mesh::Hex8::FineLayer mesh(27,27,27);

  // vector-definition
  xGooseFEM::Vector vec(mesh.conn(), mesh.dofsPeriodic());

  // element definition, with Gauss quadrature
  xGooseFEM::Element::Hex8::Quadrature quad( vec.asElement(mesh.coor()) );

  // macroscopic deformation gradient and strain
  // - zero-initialize
  T2 F = xt::zeros<double>({3,3});
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
  // - check each component
  REQUIRE( xt::allclose(Fi, 0.) );
}

// =================================================================================================

}

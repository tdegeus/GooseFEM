
#include <catch/catch.hpp>

#define CPPMAT_NOCONVERT
#include <cppmat/cppmat.h>

#include <Eigen/Eigen>

#include <GooseFEM/GooseFEM.h>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

using T2  = cppmat::tiny::cartesian::tensor2 <double,2>;
using T2s = cppmat::tiny::cartesian::tensor2s<double,2>;

// =================================================================================================

TEST_CASE("GooseFEM::ElementQuad4", "ElementQuad4.h")
{

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
  GooseFEM::ArrD rho({mesh.nelem(), quad.nip()});
  rho.setConstant(1.);

  // evaluate integral and assemble diagonal matrix (e.g. mass matrix)
  mat.assemble(quad.int_N_scalar_NT_dV(rho));

  // check matrix
  // - get the matrix
  GooseFEM::ColD M = mat.asDiagonal();
  // - check the size
  REQUIRE( M.size() == vec.ndof() );
  // - check each component
  for ( auto i = 0 ; i < M.size() ; ++i )
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

  // macroscopic deformation gradient
  // - allocate
  T2 F;
  // - zero-initialize
  F.setZero();
  // - set non-zero components
  F(0,1) = 0.1;

  // convert the macroscopic strain tensor
  T2 EPS = .5 * ( F + F.T() );

  // nodal coordinates
  GooseFEM::MatD coor = mesh.coor();;

  // nodal displacement
  // - allocate
  GooseFEM::MatD disp(mesh.nnode(), mesh.ndim());
  // - zero-initialize
  disp.setZero();

  // apply macroscopic deformation gradient
  for ( size_t n = 0 ; n < mesh.nnode() ; ++n )
    for ( size_t i = 0 ; i < F.ndim() ; ++i )
      for ( size_t j = 0 ; j < F.ndim() ; ++j )
        disp(n,i) += F(i,j) * coor(n,j);

  // compute quadrature point tensors
  GooseFEM::ArrD eps = quad.symGradN_vector(vec.asElement(disp));

  // compute volume averaged tensor
  GooseFEM::ArrD epsbar = eps.average(quad.dV(eps.shape(-1)), {0,1});

  // check
  // - temporary tensor, to view the tensors
  cppmat::view::cartesian::tensor2s<double,2> Eps;
  // - check sizes
  REQUIRE( eps.shape(0) == mesh.nelem() );
  REQUIRE( eps.shape(1) == quad.nip()   );
  REQUIRE( eps.shape(2) == Eps.size()   );
  // - check all components
  for ( size_t e = 0 ; e < mesh.nelem() ; ++e ) {
    for ( size_t k = 0 ; k < quad.nip() ; ++k ) {
      Eps.setMap(&eps(e,k));
      for ( size_t i = 0 ; i < Eps.ndim() ; ++i )
        for ( size_t j = 0 ; j < Eps.ndim() ; ++j )
          EQ( Eps(i,j), EPS(i,j) );
    }
  }

  // check macroscopic tensor
  // - convert to tensor object
  T2s Epsbar = T2s::Copy(epsbar.begin(), epsbar.end());
  // - check all components
  for ( size_t i = 0 ; i < Epsbar.ndim() ; ++i )
    for ( size_t j = 0 ; j < Epsbar.ndim() ; ++j )
      EQ( Epsbar(i,j), EPS(i,j) );
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

  // macroscopic deformation gradient
  // - allocate
  T2 F;
  // - zero-initialize
  F.setZero();
  // - set non-zero components
  F(0,1) = 0.1;

  // nodal coordinates
  GooseFEM::MatD coor = mesh.coor();;

  // nodal displacement
  // - allocate
  GooseFEM::MatD disp(mesh.nnode(), mesh.ndim());
  // - zero-initialize
  disp.setZero();

  // apply macroscopic deformation gradient
  for ( size_t n = 0 ; n < mesh.nnode() ; ++n )
    for ( size_t i = 0 ; i < F.ndim() ; ++i )
      for ( size_t j = 0 ; j < F.ndim() ; ++j )
        disp(n,i) += F(i,j) * coor(n,j);

  // compute quadrature point tensors
  GooseFEM::ArrD eps = quad.symGradN_vector(vec.asElement(disp));

  // nodal force vector (should be zero, as it is only sensitive to periodic fluctuations)
  GooseFEM::ColD Fi = vec.assembleDofs(quad.int_gradN_dot_tensor2s_dV(eps));

  // check
  // - size
  REQUIRE( Fi.size() == vec.ndof() );
  // - check all components
  for ( size_t i = 0 ; i < vec.ndof() ; ++i )
    EQ( Fi(i), 0 );

}

// =================================================================================================

}

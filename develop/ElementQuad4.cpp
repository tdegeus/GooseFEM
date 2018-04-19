
#include <catch/catch.hpp>

#define CPPMAT_NOCONVERT
#include <cppmat/cppmat.h>

#include <Eigen/Eigen>

#include <GooseFEM/GooseFEM.h>

// =================================================================================================

TEST_CASE("GooseFEM::ElementQuad4", "ElementQuad4.h")
{

// =================================================================================================

SECTION( "int_N_scalar_NT_dV" )
{
  GooseFEM::Mesh::Quad4::Regular mesh(3,3);

  GooseFEM::Vector         vec(mesh.conn(), mesh.dofsPeriodic());
  GooseFEM::MatrixDiagonal mat(mesh.conn(), mesh.dofsPeriodic());

  GooseFEM::Element::Quad4::Quadrature quad(
    vec.asElement(mesh.coor()),
    GooseFEM::Element::Quad4::Nodal::coordinates(),
    GooseFEM::Element::Quad4::Nodal::weights()
  );

  GooseFEM::ArrD rho({mesh.nelem(), quad.nip()});
  rho.setConstant(1.);

  mat.assemble(quad.int_N_scalar_NT_dV(rho));

  GooseFEM::ColD M = mat.asDiagonal();

  REQUIRE( M.size() == 9*2 );

  for ( auto i = 0 ; i < M.size() ; ++i )
    REQUIRE( M(i) == 1 );
}

// =================================================================================================

SECTION( "symGradN_vector" )
{
  GooseFEM::Mesh::Quad4::FineLayer mesh(27,27);

  GooseFEM::Vector         vec(mesh.conn(), mesh.dofs());
  GooseFEM::MatrixDiagonal mat(mesh.conn(), mesh.dofs());

  GooseFEM::Element::Quad4::Quadrature quad( vec.asElement(mesh.coor()) );

  cppmat::cartesian2d::tensor2<double> F;

  F.setZero();
  F(0,1) = 0.1;

  cppmat::cartesian2d::tensor2<double> EPS = .5 * ( F + F.T() );

  GooseFEM::MatD coor = mesh.coor();;

  GooseFEM::MatD disp(mesh.nnode(), mesh.ndim());
  disp.setZero();

  for ( size_t n = 0 ; n < mesh.nnode() ; ++n )
    for ( size_t i = 0 ; i < F.ndim() ; ++i )
      for ( size_t j = 0 ; j < F.ndim() ; ++j )
        disp(n,i) += F(i,j) * coor(n,j);

  GooseFEM::ArrD eps = quad.symGradN_vector(vec.asElement(disp));

  cppmat::view::cartesian2d::tensor2s<double> Eps;

  REQUIRE( eps.shape(0) == mesh.nelem() );
  REQUIRE( eps.shape(1) == quad.nip()   );
  REQUIRE( eps.shape(2) == Eps.size()   );

  for ( size_t e = 0 ; e < mesh.nelem() ; ++e ) {
    for ( size_t k = 0 ; k < quad.nip() ; ++k ) {
      Eps.map(&eps(e,k));
      for ( size_t i = 0 ; i < Eps.ndim() ; ++i )
        for ( size_t j = 0 ; j < Eps.ndim() ; ++j )
          REQUIRE( std::abs( Eps(i,j) - EPS(i,j) ) < 1.e-12 );
    }
  }
}

// =================================================================================================

SECTION( "symGradN_vector" )
{
  GooseFEM::Mesh::Quad4::FineLayer mesh(27,27);

  GooseFEM::Vector         vec(mesh.conn(), mesh.dofsPeriodic());
  GooseFEM::MatrixDiagonal mat(mesh.conn(), mesh.dofsPeriodic());

  GooseFEM::Element::Quad4::Quadrature quad( vec.asElement(mesh.coor()) );

  cppmat::cartesian2d::tensor2<double> F;

  F.setZero();
  F(0,1) = 0.1;

  cppmat::cartesian2d::tensor2<double> EPS = .5 * ( F + F.T() );

  GooseFEM::MatD coor = mesh.coor();;

  GooseFEM::MatD disp(mesh.nnode(), mesh.ndim());
  disp.setZero();

  for ( size_t n = 0 ; n < mesh.nnode() ; ++n )
    for ( size_t i = 0 ; i < F.ndim() ; ++i )
      for ( size_t j = 0 ; j < F.ndim() ; ++j )
        disp(n,i) += F(i,j) * coor(n,j);

  GooseFEM::ArrD eps = quad.symGradN_vector(vec.asElement(disp));

  GooseFEM::ArrD epsbar = eps.average(quad.dV(3), -2);

  std::cout << epsbar << std::endl;

  GooseFEM::ColD Fi = vec.assembleDofs(quad.int_gradN_dot_tensor2s_dV(eps));

  REQUIRE( Fi.size() == vec.ndof() );

  for ( size_t i = 0 ; i < vec.ndof() ; ++i )
    REQUIRE( std::abs( Fi(i) ) < 1.e-12 );

}

// =================================================================================================

}

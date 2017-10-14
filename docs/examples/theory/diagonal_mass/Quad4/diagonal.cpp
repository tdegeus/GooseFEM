/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <iostream>

#include <Eigen/Eigen>

#include <GooseFEM/GooseFEM.h>

// alias relevant types for <Eigen/Eigen>
using MatD = GooseFEM::MatD;
using MatS = GooseFEM::MatS;
using ColD = GooseFEM::ColD;
using ColS = GooseFEM::ColS;

// =================================================================================================

int main()
{

  // --------
  // geometry
  // --------

  // create a mesh
  GooseFEM::Mesh::Quad4::Regular mesh ( 2 , 2 );

  // element type for extrapolation and quadrature
  GooseFEM::Quad4 el;

  // extract relevant mesh data
  size_t nnode = mesh.nnode(); // number of nodes
  size_t nelem = mesh.nelem(); // number of elements
  size_t ndim  = mesh.ndim (); // number of dimensions        (== 2)
  size_t nne   = mesh.nne  (); // number of nodes per element (== 3)
  MatS   conn  = mesh.conn (); // connectivity (nodes of each element) : one row per element
  MatD   x0    = mesh.coor (); // nodal position in x- and y-direction : one row per node
  size_t ndof  = nnode * ndim; // total number of DOFs

  // DOF-numbers for each vector component of each node (sequential)
  MatS dofs = GooseFEM::Mesh::dofs(nnode,ndim);

  // perturb position for testing purpose
  x0(0,0) -= .5 ; x0(0,1) -= .2 ;
  x0(1,0) -= .1 ; x0(1,1) -= .1 ;
  x0(2,0) += .2 ; x0(2,1) -= .15;
  x0(3,0) += .1 ; x0(3,1) += .2 ;
  x0(4,0) += .1 ; x0(4,1) += .1 ;
  x0(5,0) -= .2 ; x0(5,1) += .15;
  x0(6,0) -= .5 ; x0(6,1) += .4 ;
  x0(7,0) -= .1 ; x0(7,1) += .2 ;
  x0(8,0) += .2 ; x0(8,1) += .3 ;

  // -------------------------------------
  // element & quadrature-point - allocate
  // -------------------------------------

  // element arrays
  ColS dof_e(nne*ndim     );  // array with DOF numbers of each row/column of "m_e"
  MatD x0_e (nne     ,ndim);  // nodal positions and displacements (column of vectors)

	// quadrature-point position and weight factor (this depends on the element type!)
  ColD xi_k(ndim);

  // ------------------------
  // global system - allocate
  // ------------------------

  // mass matrix: diagonal only
  ColD M(ndof);

  // density (constant in space)
  double rho = 1.0;

  // ------------------------
  // global system - assemble
  // ------------------------

  // zero-initialize global system
  M.setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - DOF numbers and nodal positions for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) x0_e .row(m)               = x0  .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) dof_e.segment(m*ndim,ndim) = dofs.row(conn(e,m));

    // - loop to integrate and assemble (N.B. quadrature-points coincide with the nodes)
    for ( size_t m = 0 ; m < nne ; ++m )
    {
      // -- define position of the quadrature-point in area coordinates (coincides with node "m")
      if      ( m == 0 ) { xi_k(0) = -1.; xi_k(1) = -1.; }
      else if ( m == 1 ) { xi_k(0) = +1.; xi_k(1) = -1.; }
      else if ( m == 2 ) { xi_k(0) = +1.; xi_k(1) = +1.; }
      else if ( m == 3 ) { xi_k(0) = -1.; xi_k(1) = +1.; }
      // -- evaluate the (gradient of the) shape functions
      el.eval( x0_e, xi_k, 1. );
      // -- assemble element mass matrix to global system
      for ( size_t i = 0; i < ndim ; ++i )
        M( dof_e(m*ndim+i) ) += rho * el.V; // N.B. "el.N(m) == 1", so omitted here
    }
  }

  // ---------------
  // print to screen
  // ---------------

  std::cout << M << std::endl;

  return 0;
}

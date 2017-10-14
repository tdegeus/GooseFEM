/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <iostream>
#include <Eigen/Eigen>
#include <cppmat/cppmat.h>
#include <GooseFEM/GooseFEM.h>
#include <GooseMaterial/Metal/LinearStrain/Elastic.h>

// alias relevant types for cppmat - exclusively 3-D
using T4  = cppmat::cartesian3d::tensor4 <double>;
using T2  = cppmat::cartesian3d::tensor2 <double>;
using T2s = cppmat::cartesian3d::tensor2s<double>;

// alias relevant types for Eigen
using MatD = GooseFEM::MatD;
using MatS = GooseFEM::MatS;
using ColD = GooseFEM::ColD;
using ColS = GooseFEM::ColS;

// alias
namespace GooseElas = GooseMaterial::Metal::LinearStrain::Elastic;

// =================================================================================================

int main()
{

  // --------
  // geometry
  // --------

  // create a mesh
  GooseFEM::Mesh::Quad4::Regular mesh(3,3);

  // element type for extrapolation and quadrature
  GooseFEM::Quad4 el;

  // extract relevant mesh data
  size_t nnode = mesh.nnode(); // number of nodes
  size_t nelem = mesh.nelem(); // number of elements
  size_t ndim  = mesh.ndim (); // number of dimensions        (== 2)
  size_t nne   = mesh.nne  (); // number of nodes per element (== 4)
  MatS   conn  = mesh.conn (); // connectivity (nodes of each element) : one row per element
  MatD   x0    = mesh.coor (); // nodal position in x- and y-direction : one row per node
  size_t ndof  = nnode * ndim; // total number of DOFs

  // nodal displacements and boundary tractions
  // - allocate as a column of vectors; like "x0"
  MatD u(nnode,ndim);
  MatD t(nnode,ndim);
  // - zero-initialize
  u.setZero();
  t.setZero();

  // DOF-numbers for each vector component of each node (sequential)
  MatS dofs = GooseFEM::Mesh::dofs(nnode,ndim);

  // -------------------
  // boundary conditions
  // -------------------

  // allocate size of the partitioned parts and displacement update (only things we need below)
  size_t nnp,nnu;
  ColD   DU_p;

  // set the prescribed displacements, and renumber DOFs in the order [ iiu , iip ]
  {
    // - get boundary nodes
    ColS   nodesRight  = mesh.nodesRight();
    ColS   nodesLeft   = mesh.nodesLeft();
    ColS   nodesBottom = mesh.nodesBottom();
    // - get sizes
    size_t nR = nodesRight .size();
    size_t nL = nodesLeft  .size();
    size_t nB = nodesBottom.size();
    // - set size of each partition
    nnp = nR + nL + nB;
    nnu = ndof - nnp;
    // - allocate prescribed displacement update
    DU_p.resize(nnp);
    // - allocate lists with prescribed DOFs
    ColS iip(nnp);
    // - fill
    size_t j = 0;
    for ( size_t i = 0 ; i < nR ; ++i ) { iip(j) = dofs(nodesRight (i), 0); DU_p(j) = 0.01; ++j; }
    for ( size_t i = 0 ; i < nL ; ++i ) { iip(j) = dofs(nodesLeft  (i), 0); DU_p(j) = 0.  ; ++j; }
    for ( size_t i = 0 ; i < nB ; ++i ) { iip(j) = dofs(nodesBottom(i), 1); DU_p(j) = 0.  ; ++j; }
    // - renumber DOFs in the following order [ iiu , iip ]; with iiu = setdiff ( dofs , iip )
    dofs = GooseFEM::Mesh::renumber(dofs,iip,"end");
  }

  // --------
  // material
  // --------

  // set elastic parameters
  double E  = 1.0;
  double nu = 0.3;
  double kappa,G;

  // convert elastic parameters to the pair required by "GooseElas::Cartesian3d::Material"
  std::tie(kappa,G) = GooseElas::ConvertParameters("E,nu",E,nu,"K,G");

  // linear elastic material
  // N.B. one definition is sufficient here because the example is homogeneous, and the material is
  // history independent; generally one definition per quadrature-point is needed
  GooseElas::Cartesian3d::Material mat = GooseElas::Cartesian3d::Material(kappa,G);

  // --------------------------------------
  // element & quadrature-point - variables
  // --------------------------------------

  // quadrature-point tensors
  // - displacement gradient (2-D)
  MatD gradu(ndim,ndim);
  // - strain, stress, and stiffness (3-D; strain has zero z-components)
  T2s eps(0.0), sig;
  T4  K4;

  // element arrays
  MatD k_e  (nne*ndim,nne*ndim); // element stiffness (matrix)
  ColD f_e  (nne*ndim         ); // element internal force (column)
  ColS dof_e(nne*ndim         ); // array with DOF numbers of each row/column of "f_e" and "k_e"
  MatD x0_e (nne     ,    ndim); // nodal positions     (column of vectors)
  MatD u_e  (nne     ,    ndim); // nodal displacements (column of vectors)

  // ------------------------
  // global system - allocate
  // ------------------------

  // partitioned system
  MatD K_uu(nnu,nnu); // stiffness           (matrix)
  MatD K_up(nnu,nnp); // "
  MatD K_pu(nnp,nnu); // "
  MatD K_pp(nnp,nnp); // "
  ColD F_u (nnu    ); // internal force      (column)
  ColD F_p (nnp    ); // "
  ColD T_u (nnu    ); // boundary traction   (column)
  ColD DU_u(nnu    ); // displacement update (column)

  // ------------------------
  // global system - assemble
  // ------------------------

  // zero-initialize global system
  K_uu.setZero();
  K_up.setZero();
  K_pu.setZero();
  K_pp.setZero();
  F_u .setZero();
  F_p .setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - DOF numbers; nodal positions and displacements for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) x0_e .row(m)               = x0  .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) u_e  .row(m)               = u   .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) dof_e.segment(m*ndim,ndim) = dofs.row(conn(e,m));

    // - zero initialize element stiffness and internal force
    k_e.setZero();
    f_e.setZero();

    // - loop over the quadrature-points
    for ( size_t k = 0 ; k < el.QuadGaussNumPoints() ; ++k )
    {
      // -- evaluate the (gradient of the) shape functions
      el.eval( x0_e, k );

      // -- local displacement gradient
      gradu = el.dNdx.transpose() * u_e;

      // -- local strain tensor (stored symmetrically, 2-D plane strain: all z-components zero)
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = i ; j < ndim ; ++j )
          eps(i,j) = 0.5 * ( gradu(i,j) + gradu(j,i) );

      // -- local constitutive response: stiffness and stress tensors
      std::tie(K4,sig) = mat.tangent_stress(eps);

      // -- add stress and stiffness tensors to element internal force column and stiffness matrix
      f_e += el.gradN_tensor2(sig);
      k_e += el.gradN_tensor4_gradNT(K4);
    }

    // - assemble element internal force to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) {
      if ( dof_e(i) < nnu ) F_u ( dof_e(i)     ) += f_e ( i );
      else                  F_p ( dof_e(i)-nnu ) += f_e ( i );
    }

    // - assemble element stiffness to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) {
      for ( size_t j = 0 ; j < nne*ndim ; ++j ) {
        if      ( dof_e(i) < nnu and dof_e(j) < nnu ) K_uu(dof_e(i)    , dof_e(j)    ) += k_e(i,j);
        else if ( dof_e(i) < nnu                    ) K_up(dof_e(i)    , dof_e(j)-nnu) += k_e(i,j);
        else if (                    dof_e(j) < nnu ) K_pu(dof_e(i)-nnu, dof_e(j)    ) += k_e(i,j);
        else                                          K_pp(dof_e(i)-nnu, dof_e(j)-nnu) += k_e(i,j);
      }
    }
  }

  // assemble boundary tractions (only relevant components)
  for ( size_t i = 0 ; i < ndof ; ++i )
    if ( dofs(i) < nnu )
      T_u ( dofs(i) ) = t ( i );

  // -----
  // solve
  // -----

  // solve for the unknown DOFs
  DU_u = K_uu.ldlt().solve( T_u - F_u - K_up * DU_p );

  // add displacement update to global system
  for ( size_t i = 0 ; i < ndof ; ++i ) {
    if ( dofs(i) < nnu ) u ( i ) += DU_u ( dofs(i)       );
    else                 u ( i ) += DU_p ( dofs(i) - nnu );
  }

  // print the result to the screen
  std::cout << "new nodal coordinates = " << std::endl << x0 + u << std::endl << std::endl;

  // compute and print reaction forces
  // - compute (only because the problem is linear, otherwise equilibrate to new "F_p")
  ColD T_p = K_pu * DU_u + K_pp * DU_p;
  // - place to global system
  for ( size_t i = 0 ; i < ndof ; ++i )
    if ( dofs(i) >= nnu )
      t ( i ) = T_p ( dofs(i) - nnu );
  // - print
  std::cout << "boundary traction = " << std::endl << t << std::endl << std::endl;

  return 0;
}

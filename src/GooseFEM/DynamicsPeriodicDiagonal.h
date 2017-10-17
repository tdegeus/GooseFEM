/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_DYNAMICS_PERIODIC_DIAGONAL_H
#define GOOSEFEM_DYNAMICS_PERIODIC_DIAGONAL_H

#include "Macros.h"
#include <cppmat/cppmat.h>

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Dynamics {
namespace Periodic {

// =========================================== OVERVIEW ============================================

template <class Element>
class DiagonalMass
{
public:

  // variables
  // ---------

  // element/quadrature/material definition (should match mesh)
  std::shared_ptr<Element> elem;

  // mesh : nodal quantities and connectivity
  MatS   dofs;    // DOF-numbers       of each node    (column of 'vectors')
  MatS   conn;    // node numbers      of each element (column of elements )
  MatD   x0;      // initial positions of each node    (column of vectors)
  MatD   u;       // displacements     of each node    (column of vectors)
  MatD   v;       // velocities        of each node    (column of vectors)
  MatD   v_n;     // velocities        of each node    (column of vectors, last increment)
  MatD   dv;      // velocities        of each node    (column of vectors, update)
  MatD   a;       // accelerations     of each node    (column of vectors)
  size_t nnode;   // number of nodes
  size_t nelem;   // number of elements
  size_t nne;     // number of nodes-per-element
  size_t ndim;    // number of dimensions
  size_t ndof;    // number of DOFs (after eliminating periodic dependencies)
  size_t refNode; // reference node (whose position will stay fixed at all times)

  // linear system
  ColD Minv;      // inverse of mass matrix (constructed diagonal -> inverse == diagonal)
  ColD Fu;        // force that depends on displacement (column)
  ColD Fv;        // force that depends on velocity (and displacement) (column)
  ColD A;         // acceleration (column)
  ColD A_n;       // acceleration (column, last increment)

  // time integration
  double dt;      // time step
  double t = 0.0; // current time
  double alpha;   // non-Galilean damping coefficient

  // constructor
  // -----------

  DiagonalMass(std::shared_ptr<Element> elem, const MatD &x0,const MatS &conn,const MatS &dofs,
    size_t refNode, double dt, double alpha=0.0 );

  // functions
  // ---------

  void increment();
  void removeRigidBodyMotion();
  void computeInverseMassMatrix();
  void computeInternalForce();
  void computeDampingForce();

};

// ========================================== SOURCE CODE ==========================================

template<class Element>
DiagonalMass<Element>::DiagonalMass(
  std::shared_ptr<Element> _elem, const MatD &_x0, const MatS &_conn, const MatS &_dofs,
  size_t _refNode, double _dt, double _alpha
)
{
  // problem specific settings
  elem    = _elem;

  // mesh definition
  x0      = _x0;
  conn    = _conn;
  dofs    = _dofs;
  refNode = _refNode;

  // time integration
  dt      = _dt;
  alpha   = _alpha;

  // extract mesh size from definition
  nnode   = static_cast<size_t>(x0  .rows());
  ndim    = static_cast<size_t>(x0  .cols());
  nelem   = static_cast<size_t>(conn.rows());
  nne     = static_cast<size_t>(conn.cols());
  ndof    = dofs.maxCoeff()+1;

  // nodal quantities
  u  .conservativeResize(nnode,ndim); u  .setZero();
  v  .conservativeResize(nnode,ndim); v  .setZero();
  v_n.conservativeResize(nnode,ndim); v_n.setZero();
  dv .conservativeResize(nnode,ndim); dv .setZero();
  a  .conservativeResize(nnode,ndim); a  .setZero();

  // linear system
  Minv .conservativeResize(ndof); Minv .setZero();
  Fu   .conservativeResize(ndof); Fu   .setZero();
  Fv   .conservativeResize(ndof); Fv   .setZero();
  A    .conservativeResize(ndof); A    .setZero();
  A_n  .conservativeResize(ndof); A_n  .setZero();

  // compute inverse mass matrix
  computeInverseMassMatrix();
}

// =================================================================================================

template<class Element>
void DiagonalMass<Element>::increment()
{
  // (1) new nodal positions (displacements)
  // - apply update (nodal) : x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  u += dt * v + ( .5 * std::pow(dt,2.) ) * a;
  // - remove rigid body motion ( u(refNode,:) == 0. )
  removeRigidBodyMotion();
  // - compute forces that are dependent on the displacement, but not on the velocity
  computeInternalForce();

  // (2a) estimate nodal velocities
  // - update nodal velocities
  v = v_n + dt * a;
  // - compute forces that are dependent on the velocity
  computeDampingForce();
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - Fu - Fv );
  // - update nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<ndof; ++i ) v(i) = v_n(i) + ( .5 * dt ) * ( A(dofs(i)) + A_n(dofs(i)) );
  // - compute forces that are dependent on the velocity
  computeDampingForce();

  // (2b) new nodal velocities
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - Fu - Fv );
  // - update nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<ndof; ++i ) v(i) = v_n(i) + ( .5 * dt ) * ( A(dofs(i)) + A_n(dofs(i)) );
  // - compute forces that are dependent on the velocity
  computeDampingForce();

  // (3) new nodal accelerations
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - Fu - Fv );
  // - set nodal accelerations (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<ndof; ++i ) a(i) = A(dofs(i));

  // store history
  A_n = A;  // accelerations (DOFs)
  v_n = v;  // nodal velocities
  t  += dt; // current time instance

  // N.B. at this point:
  // "a", "A", "A_n"  ->  new nodal accelerations, their DOF equivalents, and a 'back-up'
  // "v",      "v_n"  ->  new nodal velocities,                           and a 'back-up'
  // "u"              ->  new nodal displacements
  // The forces "Fu" and "Fv" correspond to this state of the system
}

// =================================================================================================

template<class Element>
void DiagonalMass<Element>::removeRigidBodyMotion()
{
  // remove the amount of deformation applied the reference node
  for ( size_t i = 0 ; i < nnode ; ++i )
    for ( size_t j = 0 ; j < ndim ; ++j )
      if ( i != refNode )
        u(i,j) -= u(refNode,j);

  // remove the deformation of the reference node
  for ( size_t j = 0 ; j < ndim ; ++j )
    u(refNode,j) = 0.;
}

// =================================================================================================

template<class Element>
void DiagonalMass<Element>::computeInverseMassMatrix()
{
  // array with DOF numbers of each node
  cppmat::matrix2<size_t> dofe(nne,ndim);

  // zero-initialize mass matrix (only the inverse is stored)
  ColD M(ndof); M.setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - nodal positions, displacement, and DOF numbers for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) {
      for ( size_t i = 0 ; i < ndim ; ++i ) {
        elem->xe(m,i) = x0  (conn(e,m),i);
        elem->ue(m,i) = u   (conn(e,m),i);
        dofe    (m,i) = dofs(conn(e,m),i);
      }
    }

    // - compute element mass matrix
    elem->computeMassMatrix(e);

    // - assemble element mass matrix : take only the diagonal (approximation!!)
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) M(dofe[i]) += elem->M(i,i);
  }

  // compute inverse of the mass matrix
  Minv = M.cwiseInverse();
}

// =================================================================================================

template<class Element>
void DiagonalMass<Element>::computeInternalForce()
{
  // array with DOF numbers of each node
  cppmat::matrix2<size_t> dofe(nne,ndim);

  // zero-initialize internal force
  Fu.setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - nodal positions, displacement, and DOF numbers for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) {
      for ( size_t i = 0 ; i < ndim ; ++i ) {
        elem->xe(m,i) = x0  (conn(e,m),i);
        elem->ue(m,i) = u   (conn(e,m),i);
        dofe    (m,i) = dofs(conn(e,m),i);
      }
    }

    // - compute element mass matrix
    elem->computeInternalForce(e);

    // - assemble internal force to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) Fu(dofe[i]) += elem->fu(i);
  }
}

// =================================================================================================

template<class Element>
void DiagonalMass<Element>::computeDampingForce()
{
  // array with DOF numbers of each node
  cppmat::matrix2<size_t> dofe(nne,ndim);

  // zero-initialize internal force
  Fv.setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - nodal positions, displacement, velocity, and DOF numbers for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) {
      for ( size_t i = 0 ; i < ndim ; ++i ) {
        elem->xe(m,i) = x0  (conn(e,m),i);
        elem->ue(m,i) = u   (conn(e,m),i);
        elem->ve(m,i) = v   (conn(e,m),i);
        dofe    (m,i) = dofs(conn(e,m),i);
      }
    }

    // - compute element mass matrix
    elem->computeDampingForce(e);

    // - assemble internal force to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) Fv(dofe[i]) += elem->fv(i);
  }
}

// =================================================================================================

} // namespace Periodic
} // namespace Dynamics
} // namespace GooseFEM

#endif

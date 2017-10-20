/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_DYNAMICS_DIAGONAL_PERIODIC_H
#define GOOSEFEM_DYNAMICS_DIAGONAL_PERIODIC_H

#include "Macros.h"
#include <cppmat/cppmat.h>

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Dynamics {
namespace Diagonal {

// =========================================== OVERVIEW ============================================

template <class Element>
class Periodic
{
public:

  // variables
  // ---------

  // element/quadrature/material definition (should match mesh)
  std::shared_ptr<Element> elem;

  // mesh : nodal quantities and connectivity
  MatS   dofs;    // DOF-numbers       of each node    (column of 'vectors')
  MatS   conn;    // node numbers      of each element (column of elements)
  MatD   x0;      // initial positions of each node    (column of vectors)
  MatD   u;       // displacements     of each node    (column of vectors)
  MatD   v;       // velocities        of each node    (column of vectors)
  MatD   a;       // accelerations     of each node    (column of vectors)
  size_t nnode;   // number of nodes
  size_t nelem;   // number of elements
  size_t nne;     // number of nodes-per-element
  size_t ndim;    // number of dimensions
  size_t ndof;    // number of DOFs (after eliminating periodic dependencies)

  // linear system
  ColD Minv;      // inverse of mass matrix (constructed diagonal -> inverse == diagonal)
  ColD Fu;        // force that depends on displacement (column)
  ColD Fv;        // force that depends on velocity (and displacement) (column)
  ColD V;         // velocity (column)
  ColD V_n;       // velocity (column, last increment)
  ColD A;         // acceleration (column)
  ColD A_n;       // acceleration (column, last increment)

  // time integration
  double dt;      // time step
  double t = 0.0; // current time
  double alpha;   // non-Galilean damping coefficient

  // constructor
  // -----------

  Periodic(std::shared_ptr<Element> elem, const MatD &x0, const MatS &conn, const MatS &dofs,
    double dt, double alpha=0.0 );

  // functions
  // ---------

  void velocityVerlet();  // one time step of time integrator
  void Verlet();          // one time step of time integrator (Fv == 0)
  void computeMinv();     // element loop: inverse mass matrix            <- "elem->computeM"
  void computeFu();       // element loop: displacement dependent forces  <- "elem->computeFu"
  void computeFv();       // element loop: velocity dependent forces      <- "elem->computeFv"
  void post();            // element loop: post-process                   <- "elem->post"

};

// ========================================== SOURCE CODE ==========================================

template<class Element>
Periodic<Element>::Periodic(
  std::shared_ptr<Element> _elem, const MatD &_x0, const MatS &_conn, const MatS &_dofs,
  double _dt, double _alpha
)
{
  // problem specific settings
  elem  = _elem;

  // mesh definition
  x0    = _x0;
  conn  = _conn;
  dofs  = _dofs;

  // time integration
  dt    = _dt;
  alpha = _alpha;

  // extract mesh size from definition
  nnode = static_cast<size_t>(x0  .rows());
  ndim  = static_cast<size_t>(x0  .cols());
  nelem = static_cast<size_t>(conn.rows());
  nne   = static_cast<size_t>(conn.cols());
  ndof  = dofs.maxCoeff()+1;

  // basic check (mostly the user is 'trusted')
  assert( dofs.size() == nnode * ndim );
  assert( ndof        <  nnode * ndim );

  // nodal quantities
  u.conservativeResize(nnode,ndim); u.setZero();
  v.conservativeResize(nnode,ndim); v.setZero();
  a.conservativeResize(nnode,ndim); a.setZero();

  // linear system (DOFs)
  Minv.conservativeResize(ndof); Minv.setZero();
  Fu  .conservativeResize(ndof); Fu  .setZero();
  Fv  .conservativeResize(ndof); Fv  .setZero();
  V   .conservativeResize(ndof); V   .setZero();
  V_n .conservativeResize(ndof); V_n .setZero();
  A   .conservativeResize(ndof); A   .setZero();
  A_n .conservativeResize(ndof); A_n .setZero();

  // compute inverse mass matrix : assumed constant in time
  computeMinv();
}

// =================================================================================================

template<class Element>
void Periodic<Element>::velocityVerlet()
{
  // (1) new nodal positions (displacements)
  // - apply update (nodal) : x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  u += dt * v + ( .5 * std::pow(dt,2.) ) * a;
  // - compute forces that are dependent on the displacement, but not on the velocity
  computeFu();

  // (2a) estimate nodal velocities
  // - update velocities (DOFs)
  V.noalias() = V_n + dt * A;
  // - convert to nodal velocity (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - compute forces that are dependent on the velocity
  computeFv();
  // - solve for accelerations (DOFs)
  A.noalias() = Minv.cwiseProduct( - Fu - Fv - alpha*V );
  // - update velocities (DOFs)
  V.noalias() = V_n + ( .5 * dt ) * ( A_n + A );
  // - convert to nodal velocity (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - compute forces that are dependent on the velocity
  computeFv();

  // (2b) new nodal velocities
  // - solve for accelerations (DOFs)
  A.noalias() = Minv.cwiseProduct( - Fu - Fv - alpha*V );
  // - update velocities (DOFs)
  V.noalias() = V_n + ( .5 * dt ) * ( A_n + A );
  // - convert to nodal velocity (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - compute forces that are dependent on the velocity
  computeFv();

  // (3) new nodal accelerations
  // - solve for accelerations (DOFs)
  A.noalias() = Minv.cwiseProduct( - Fu - Fv - alpha*V );
  // - convert to nodal acceleration (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) a(i) = A(dofs(i));

  // store history
  A_n = A;  // accelerations (DOFs)
  V_n = V;  // nodal velocities
  t  += dt; // current time instance

  // N.B. at this point:
  // "a" == "A" == "A_n"  ->  new nodal accelerations, their DOF equivalents, and a 'back-up'
  // "v" ==        "V_n"  ->  new nodal velocities,                           and a 'back-up'
  // "u"                  ->  new nodal displacements
  // The forces "Fu" and "Fv" correspond to this state of the system
}

// =================================================================================================

template<class Element>
void Periodic<Element>::Verlet()
{
  // (1) new nodal positions (displacements)
  // - apply update (nodal) : x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  u += dt * v + ( .5 * std::pow(dt,2.) ) * a;
  // - compute forces that are dependent on the displacement, but not on the velocity
  computeFu();
  // - solve for accelerations (DOFs)
  A.noalias() = Minv.cwiseProduct( - Fu );
  // - convert to nodal acceleration (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) a(i) = A(dofs(i));

  // (2) propagate velocities
  // - update velocities (DOFs)
  V.noalias() = V_n + ( .5 * dt ) * ( A_n + A );
  // - convert to nodal velocity (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));

  // store history
  A_n = A;  // accelerations (DOFs)
  V_n = V;  // nodal velocities
  t  += dt; // current time instance

  // N.B. at this point:
  // "a" == "A" == "A_n"  ->  new nodal accelerations, their DOF equivalents, and a 'back-up'
  // "v" ==        "V_n"  ->  new nodal velocities,                           and a 'back-up'
  // "u"                  ->  new nodal displacements
  // The forces "Fu" correspond to this state of the system
}

// =================================================================================================

template<class Element>
void Periodic<Element>::computeMinv()
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

    // - compute element mass matrix using problem specific "Element" class
    elem->computeM(e);

    // - assemble element mass matrix : take only the diagonal
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) M(dofe[i]) += elem->M(i,i);

    // - check that the user provided a diagonal mass matrix
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < nne*ndim ; ++i ) {
        for ( size_t j = 0 ; j < nne*ndim ; ++j ) {
          if ( i != j ) assert( ! elem->M(i,j) );
          else          assert(   elem->M(i,i) );
        }
      }
    #endif
  }

  // compute inverse of the mass matrix
  Minv = M.cwiseInverse();
}

// =================================================================================================

template<class Element>
void Periodic<Element>::computeFu()
{
  // array with DOF numbers of each node
  cppmat::matrix2<size_t> dofe(nne,ndim);

  // zero-initialize displacement dependent force
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

    // - compute element force using problem specific "Element" class
    elem->computeFu(e);

    // - assemble force to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) Fu(dofe[i]) += elem->fu(i);
  }
}

// =================================================================================================

template<class Element>
void Periodic<Element>::computeFv()
{
  // array with DOF numbers of each node
  cppmat::matrix2<size_t> dofe(nne,ndim);

  // zero-initialize velocity dependent force
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

    // - compute element force using problem specific "Element" class
    elem->computeFv(e);

    // - assemble force to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) Fv(dofe[i]) += elem->fv(i);
  }
}

// =================================================================================================

template<class Element>
void Periodic<Element>::post()
{
  // array with DOF numbers of each node
  cppmat::matrix2<size_t> dofe(nne,ndim);

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

    // - run post-process function of the problem specific "Element" class (output is handled there)
    elem->post(e);
  }
}

// =================================================================================================

}}} // namespace GooseFEM::Dynamics::Diagonal

#endif

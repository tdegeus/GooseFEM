/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_DYNAMICS_DIAGONAL_SEMIPERIODIC_H
#define GOOSEFEM_DYNAMICS_DIAGONAL_SEMIPERIODIC_H

#include "Macros.h"
#include <cppmat/cppmat.h>

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Dynamics {
namespace Diagonal {

// =========================================== OVERVIEW ============================================

template <class Element>
class SemiPeriodic
{
public:

  // variables
  // ---------

  // element/quadrature/material definition
  std::unique_ptr<Element> elem;

  // mesh: nodal position/displacement/velocity/acceleration, DOF-numbers, connectivity, dimensions
  size_t nnode, nelem, nne, ndim, ndof;
  MatS   conn, dofs;
  MatD   x, u, v, a;

  // fixed DOFs: prescribed velocity/acceleration, DOF-numbers, dimensions
  size_t nfixed;
  ColS   fixedDofs;
  ColD   fixedV, fixedA;

  // linear system: columns (also "M" and "D" which are diagonal)
  ColD   M, Minv, D, F, V, V_n, A, A_n;

  // time integration
  double dt, t=0.0;

  // constructor
  // -----------

  SemiPeriodic(){};

  SemiPeriodic(
    std::unique_ptr<Element> elem, const MatD &x, const MatS &conn, const MatS &dofs,
    const ColS &fixedDofs, double dt
  );

  // functions
  // ---------

  void velocityVerlet();           // one time step of time integrator
  void Verlet();                   // one time step of time integrator (no velocity dependence)
  void updated_x();                // process update in "x"
  void updated_u(bool init=false); // process update in "u", if init all possible updates are made
  void updated_v(bool init=false); // process update in "v", if init all possible updates are made
  void assemble_M();               // assemble the mass matrix from the element mass matrices
  void assemble_D();               // assemble the damping matrix from the element damping matrices
  void assemble_F();               // assemble the force from the element forces
};

// ========================================== SOURCE CODE ==========================================

template<class Element>
SemiPeriodic<Element>::SemiPeriodic(
  std::unique_ptr<Element> _elem, const MatD &_x, const MatS &_conn, const MatS &_dofs,
  const ColS &_fixedDofs, double _dt
)
{
  // copy input
  elem      = std::move(_elem);
  x         = _x;
  conn      = _conn;
  dofs      = _dofs;
  dt        = _dt;
  fixedDofs = _fixedDofs;

  // compute sizes
  nfixed    = static_cast<size_t>(fixedDofs.size());
  nnode     = static_cast<size_t>(x.rows());
  ndim      = static_cast<size_t>(x.cols());
  nelem     = static_cast<size_t>(conn.rows());
  nne       = static_cast<size_t>(conn.cols());
  ndof      = dofs.maxCoeff()+1;

  // basic checks (mostly the user is 'trusted')
  assert( dofs.size() == nnode * ndim );
  assert( ndof        <  nnode * ndim );

  #ifndef NDEBUG
  for ( size_t i = 0 ; i < nfixed ; ++i ) assert( fixedDofs(i) < ndof );
  #endif

  // allocate and zero-initialize nodal quantities
  u.conservativeResize(nnode,ndim); u.setZero();
  v.conservativeResize(nnode,ndim); v.setZero();
  a.conservativeResize(nnode,ndim); a.setZero();

  // allocate and zero-initialize linear system (DOFs)
  Minv.conservativeResize(ndof);
  M   .conservativeResize(ndof);
  D   .conservativeResize(ndof);
  F   .conservativeResize(ndof);
  V   .conservativeResize(ndof);
  V_n .conservativeResize(ndof); V_n .setZero();
  A   .conservativeResize(ndof);
  A_n .conservativeResize(ndof); A_n .setZero();

  // fixed DOFs : default zero velocity and acceleration
  fixedV.conservativeResize(nfixed); fixedV.setZero();
  fixedA.conservativeResize(nfixed); fixedA.setZero();

  // initialize all fields
  updated_x();
  updated_u(true);
  updated_v(true);
}

// =================================================================================================

template<class Element>
void SemiPeriodic<Element>::velocityVerlet()
{
  // (1) new positions (displacements)
  // - apply update (nodal) : x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  u += dt * v + ( .5 * std::pow(dt,2.) ) * a;
  // - process update in displacements
  updated_u();

  // (2a) estimate new velocities
  // - update velocities (DOFs)
  V.noalias() = V_n + dt * A;
  // - apply the fixed velocities
  for ( size_t i=0; i<nfixed; ++i ) V(fixedDofs(i)) = fixedV(i);
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - process update in velocities
  updated_v();
  // - solve for accelerations (DOFs)
  A.noalias() = Minv.cwiseProduct( - F - D.cwiseProduct(V) );
  // - update velocities (DOFs)
  V.noalias() = V_n + ( .5 * dt ) * ( A_n + A );
  // - apply the fixed velocities
  for ( size_t i=0; i<nfixed; ++i ) V(fixedDofs(i)) = fixedV(i);
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - process update in velocities
  updated_v();

  // (2b) new velocities
  // - solve for accelerations (DOFs)
  A.noalias() = Minv.cwiseProduct( - F - D.cwiseProduct(V) );
  // - update velocities (DOFs)
  V.noalias() = V_n + ( .5 * dt ) * ( A_n + A );
  // - apply the fixed velocities
  for ( size_t i=0; i<nfixed; ++i ) V(fixedDofs(i)) = fixedV(i);
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - process update in velocities
  updated_v();

  // (3) new accelerations
  // - solve for accelerations (DOFs)
  A.noalias() = Minv.cwiseProduct( - F - D.cwiseProduct(V) );
  // - apply the fixed accelerations
  for ( size_t i=0; i<nfixed; ++i ) A(fixedDofs(i)) = fixedA(i);
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
  // The forces "F" correspond to this state of the system
}

// =================================================================================================

template<class Element>
void SemiPeriodic<Element>::Verlet()
{
  // (1) new nodal positions (displacements)
  // - apply update (nodal) : x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  u += dt * v + ( .5 * std::pow(dt,2.) ) * a;
  // - process update in displacements
  updated_u();
  // - solve for accelerations (DOFs)
  A.noalias() = Minv.cwiseProduct( - F );
  // - apply the fixed accelerations
  for ( size_t i=0; i<nfixed; ++i ) A(fixedDofs(i)) = fixedA(i);
  // - convert to nodal acceleration (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) a(i) = A(dofs(i));

  // (2) propagate velocities
  // - update velocities (DOFs)
  V.noalias() = V_n + ( .5 * dt ) * ( A_n + A );
  // - apply the fixed velocities
  for ( size_t i=0; i<nfixed; ++i ) V(fixedDofs(i)) = fixedV(i);
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));

  // store history
  A_n = A;  // accelerations (DOFs)
  V_n = V;  // nodal velocities
  t  += dt; // current time instance

  // N.B. at this point:
  // "a" == "A" == "A_n"  ->  new nodal accelerations, their DOF equivalents, and a 'back-up'
  // "v" ==        "V_n"  ->  new nodal velocities,                           and a 'back-up'
  // "u"                  ->  new nodal displacements
  // The forces "F" correspond to this state of the system
}

// =================================================================================================

template<class Element>
void SemiPeriodic<Element>::updated_x()
{
  // set the nodal positions of all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < nelem ; ++e )
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        elem->x(e,m,i) = x(conn(e,m),i);

  // signal update
  elem->updated_x();
}

// =================================================================================================

template<class Element>
void SemiPeriodic<Element>::updated_u(bool init)
{
  // set the nodal displacements of all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < nelem ; ++e )
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        elem->u(e,m,i) = u(conn(e,m),i);

  // signal update
  elem->updated_u();

  // update
  if ( elem->changed_M or init ) assemble_M();
  if ( elem->changed_D or init ) assemble_D();
  if ( elem->changed_f or init ) assemble_F();
}

// =================================================================================================

template<class Element>
void SemiPeriodic<Element>::updated_v(bool init)
{
  // set the nodal displacements of all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < nelem ; ++e )
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        elem->v(e,m,i) = v(conn(e,m),i);

  // signal update
  elem->updated_v();

  // update
  if ( elem->changed_M or init ) assemble_M();
  if ( elem->changed_D or init ) assemble_D();
  if ( elem->changed_f or init ) assemble_F();
}

// =================================================================================================

template<class Element>
void SemiPeriodic<Element>::assemble_M()
{
  // zero-initialize
  M.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);
  // assemble
  #pragma omp parallel
  {
    // - mass matrix, per thread
    ColD M_(ndof);
    M_.setZero();

    // - assemble diagonal mass matrix, per thread
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
      for ( size_t m = 0 ; m < nne ; ++m )
        for ( size_t i = 0 ; i < ndim ; ++i )
          M_(dofs(conn(e,m),i)) += elem->M(e,m*ndim+i,m*ndim+i);

    // - reduce "M_" per thread to total "M"
    #pragma omp critical
      M += M_;
  }
  // automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  // compute inverse of the mass matrix
  Minv = M.cwiseInverse();
}

// =================================================================================================

template<class Element>
void SemiPeriodic<Element>::assemble_D()
{
  // zero-initialize
  D.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);
  // assemble
  #pragma omp parallel
  {
    // - damping matrix, per thread
    ColD D_(ndof);
    D_.setZero();

    // - assemble diagonal damping matrix, per thread
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
      for ( size_t m = 0 ; m < nne ; ++m )
        for ( size_t i = 0 ; i < ndim ; ++i )
          D_(dofs(conn(e,m),i)) += elem->D(e,m*ndim+i,m*ndim+i);

    // - reduce "D_" per thread to total "D"
    #pragma omp critical
      D += D_;
  }
  // automatic parallelization by Eigen
  Eigen::setNbThreads(0);
}

// =================================================================================================

template<class Element>
void SemiPeriodic<Element>::assemble_F()
{
  // zero-initialize
  F.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);
  // assemble
  #pragma omp parallel
  {
    // - force, per thread
    ColD F_(ndof);
    F_.setZero();

    // - assemble force, per thread
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
      for ( size_t m = 0 ; m < nne ; ++m )
        for ( size_t i = 0 ; i < ndim ; ++i )
          F_(dofs(conn(e,m),i)) += elem->f(e,m,i);

    // - reduce "F_" per thread to total "F"
    #pragma omp critical
      F += F_;
  }
  // automatic parallelization by Eigen
  Eigen::setNbThreads(0);
}

// =================================================================================================

}}} // namespace GooseFEM::Dynamics::Diagonal

#endif

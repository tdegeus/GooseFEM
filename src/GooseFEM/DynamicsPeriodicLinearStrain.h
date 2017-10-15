/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_DYNAMICS_PERIODIC_LINEARSTRAIN_H
#define GOOSEFEM_DYNAMICS_PERIODIC_LINEARSTRAIN_H

#include "Macros.h"

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Dynamics {
namespace Periodic {
namespace LinearStrain {
namespace DiagonalMass {

// =========================================== OVERVIEW ============================================

template<class Material, class Element, class Tensor>
class Simulation
{
public:

  // variables
  // ---------

  // problem specific settings
  std::shared_ptr<Material> mat;  // material definition per integration point per element
  std::shared_ptr<Element>  elem; // element definition (should match mesh)

  // mesh : nodal quantities and connectivity
  MatS   dofs;    // DOF-numbers       of each node     (column of 'vectors')
  MatS   conn;    // node number       of each element  (column of elements )
  MatD   x0;      // initial positions of each node     (column of vectors  )
  MatD   u;       // displacements     of each node     (column of vectors  )
  MatD   v;       // velocities        of each node     (column of vectors  )
  MatD   v_n;     // velocities        of each node     (column of vectors  ), last increment
  MatD   dv;      // velocities        of each node     (column of vectors  ), update
  MatD   a;       // accelerations     of each node     (column of vectors  )
  size_t nnode;   // number of nodes
  size_t nelem;   // number of elements
  size_t nne;     // number of nodes-per-element  (==4)
  size_t ndim;    // number of dimensions         (==2)
  size_t ndof;    // number of DOFs (after eliminating periodic dependencies)

  // linear system
  ColD Minv;      // inverse of mass matrix (constructed diagonal -> inverse diagonal)
  ColD Fint;      // internal forces at the DOFs
  ColD Fdamp;     // damping  forces at the DOFs
  ColD A;         // accelerations   of the DOFs
  ColD A_n;       // accelerations   of the DOFs, last increment

  // time integration
  double dt;      // time step
  double t = 0.0; // current time

  // constructor
  // -----------

  Simulation(
    std::shared_ptr<Material> mat,
    std::shared_ptr<Element>  el,
    const MatD &x0,
    const MatS &conn,
    const MatS &dofs,
    double      dt
  );

  // functions
  // ---------

  void increment();
  void suppressRigidBodyMotion();
  void computeInverseMassMatrix();
  void computeInternalForce();
  void computeDampingForce();

};

// ========================================== SOURCE CODE ==========================================

template<class Material, class Element, class Tensor>
Simulation<Material,Element,Tensor>::Simulation(
  std::shared_ptr<Material> _mat,
  std::shared_ptr<Element>  _el,
  const MatD               &_x0,
  const MatS               &_conn,
  const MatS               &_dofs,
  double                    _dt
)
{
  // problem specific settings
  mat   = _mat;
  elem  = _el;

  // mesh definition
  x0    = _x0;
  conn  = _conn;
  dofs  = _dofs;

  // time integration
  dt    = _dt;

  // extract mesh size from definition
  nnode = static_cast<size_t>(x0  .rows());
  ndim  = static_cast<size_t>(x0  .cols());
  nelem = static_cast<size_t>(conn.rows());
  nne   = static_cast<size_t>(conn.cols());
  ndof  = dofs.maxCoeff()+1;

  // nodal quantities
  u  .conservativeResize(nnode,ndim); u  .setZero();
  v  .conservativeResize(nnode,ndim); v  .setZero();
  v_n.conservativeResize(nnode,ndim); v_n.setZero();
  dv .conservativeResize(nnode,ndim); dv .setZero();
  a  .conservativeResize(nnode,ndim); a  .setZero();

  // linear system
  Minv .conservativeResize(ndof); Minv .setZero();
  Fint .conservativeResize(ndof); Fint .setZero();
  Fdamp.conservativeResize(ndof); Fdamp.setZero();
  A    .conservativeResize(ndof); A    .setZero();
  A_n  .conservativeResize(ndof); A_n  .setZero();
}

// =================================================================================================

template<class Material, class Element, class Tensor>
void Simulation<Material,Element,Tensor>::increment()
{
  // (1) new nodal positions (displacements)
  // - apply update (nodal) : x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
  u += dt * v + ( .5 * std::pow(dt,2.) ) * a;
  // - compute forces that are dependent on the displacement, but not on the velocity
  computeInternalForce();

  // (2a) estimate nodal velocities
  // - initial estimate the nodal velocity update
  dv = dt * a;
  // - suppress rigid body motion: mean(dv)==0.
  suppressRigidBodyMotion();
  // - update nodal velocities
  v = v_n + dv;
  // - compute forces that are dependent on the velocity
  computeDampingForce();
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - Fint - Fdamp );
  // - set nodal velocity update (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<ndof; ++i ) dv(i) = ( .5 * dt ) * ( A(dofs(i)) + A_n(dofs(i)) );
  // - suppress rigid body motion: mean(dv)==0.
  suppressRigidBodyMotion();
  // - update nodal velocities
  v = v_n + dv;
  // - compute forces that are dependent on the velocity
  computeDampingForce();

  // (2b) new nodal velocities
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - Fint - Fdamp );
  // - set nodal velocity update (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<ndof; ++i ) dv(i) = ( .5 * dt ) * ( A(dofs(i)) + A_n(dofs(i)) );
  // - suppress rigid body motion: mean(dv)==0.
  suppressRigidBodyMotion();
  // - update nodal velocities
  v = v_n + dv;
  // - compute forces that are dependent on the velocity
  computeDampingForce();

  // (3) new nodal accelerations
  // - solve for accelerations (DOFs)
  A = Minv.cwiseProduct( - Fint - Fdamp );
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
  // The forces "Fint" and "Fdamp" correspond to this state of the system
}

// =================================================================================================

template<class Material, class Element, class Tensor>
void Simulation<Material,Element,Tensor>::suppressRigidBodyMotion()
{
  // zero-initialize mean velocity
  ColD dvmean(ndim); dvmean.setZero();

  // compute mean velocity == rigid body motions
  // - sum of all nodal velocities
  for ( size_t i = 0 ; i < nnode ; ++i )
    for ( size_t j = 0 ; j < ndim ; ++j )
      dvmean(j) += dv(i,j);
  // - normalize (sum -> mean)
  dvmean /= static_cast<double>(nnode);

  // remove mean to avoid rigid body motion
  for ( size_t i = 0 ; i < nnode ; ++i )
    for ( size_t j = 0 ; j < ndim ; ++j )
      dv(i,j) -= dvmean(j);
}

// =================================================================================================

template<class Material, class Element, class Tensor>
void Simulation<Material,Element,Tensor>::computeInverseMassMatrix()
{
  // element arrays
  MatD x0_e (nne,ndim); // nodal positions (column of vectors)
  ColS dof_e(nne*ndim); // array with DOF numbers of each row/column of "M"

  // zero-initialize mass matrix (only the inverse is stored)
  ColD M(ndof);
  M.setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - DOF numbers; nodal positions for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) x0_e .row(m)               = x0  .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) dof_e.segment(m*ndim,ndim) = dofs.row(conn(e,m));

    // - loop to integrate and assemble (N.B. quadrature-points coincide with the nodes)
    for ( size_t m = 0 ; m < nne ; ++m )
    {
      // -- evaluate the shape functions
      elem->eval( x0_e, elem->QuadNodesPosition(m), elem->QuadNodesWeight(m) );

      // -- assemble element mass matrix to global system
      for ( size_t i = 0 ; i < ndim; ++i )
        M( dof_e(m*ndim+i) ) += mat->density(e,m) * elem->V; // N.B. "N(m) == 1", so omitted here
    }
  }

  // compute inverse of the mass matrix
  Minv = M.cwiseInverse();
}

// =================================================================================================

template<class Material, class Element, class Tensor>
void Simulation<Material,Element,Tensor>::computeInternalForce()
{
  // quadrature-point tensors
  // - displacement gradient
  MatD gradu(ndim,ndim);
  // - strain / stress
  Tensor eps, sig;
  // - zero-initialize (in case that the number of dimensions is higher)
  eps.zeros();
  sig.zeros();

  // element arrays
  MatD x0_e (nne,ndim); // nodal positions     (column of vectors)
  MatD u_e  (nne,ndim); // nodal displacements (column of vectors)
  ColS dof_e(nne*ndim); // array with DOF numbers of each row of "f_e"
  ColD f_e  (nne*ndim); // element internal force (column)

  // zero-initialize internal force
  Fint.setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - DOF numbers; nodal positions and displacements for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) x0_e .row(m)               = x0  .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) u_e  .row(m)               = u   .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) dof_e.segment(m*ndim,ndim) = dofs.row(conn(e,m));

    // - zero initialize element internal force
    f_e.setZero();

    // - loop over the quadrature-points
    for ( size_t k = 0 ; k < elem->QuadGaussNumPoints() ; ++k )
    {
      // -- evaluate the (gradient of the) shape functions
      elem->evalGradN( x0_e + u_e , k );

      // -- local displacement gradient
      gradu = elem->dNdx.transpose() * u_e;

      // -- local strain tensor
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = i ; j < ndim ; ++j )
          eps(i,j) = 0.5 * ( gradu(i,j) + gradu(j,i) );

      // -- local constitutive response: stiffness and stress tensors
      sig = mat->stressStrain(e,k,eps,elem->V);

      // -- add stress to element internal force column
      f_e += elem->gradN_tensor2(sig);
    }

    // - assemble internal force to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) Fint( dof_e(i) ) += f_e(i);
  }
}

// =================================================================================================

template<class Material, class Element, class Tensor>
void Simulation<Material,Element,Tensor>::computeDampingForce()
{
  // quadrature-point tensors
  // - displacement gradient
  MatD gradv(ndim,ndim);
  // - strain-rate / stress
  Tensor epsdot, sig;
  // - zero-initialize (in case that the number of dimensions is higher)
  epsdot.zeros();
  sig   .zeros();

  // element arrays
  MatD x0_e (nne,ndim); // nodal positions     (column of vectors)
  MatD u_e  (nne,ndim); // nodal displacements (column of vectors)
  MatD v_e  (nne,ndim); // nodal velocities    (column of vectors)
  ColS dof_e(nne*ndim); // array with DOF numbers of each row of "f_e"
  ColD f_e  (nne*ndim); // element internal force (column)

  // zero-initialize internal force
  Fdamp.setZero();

  // loop over all elements
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // - DOF numbers; nodal positions and displacements for element "e"
    for ( size_t m = 0 ; m < nne ; ++m ) x0_e .row(m)               = x0  .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) u_e  .row(m)               = u   .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) v_e  .row(m)               = v   .row(conn(e,m));
    for ( size_t m = 0 ; m < nne ; ++m ) dof_e.segment(m*ndim,ndim) = dofs.row(conn(e,m));

    // - zero initialize element internal force
    f_e.setZero();

    // - loop over the quadrature-points
    for ( size_t k = 0 ; k < elem->QuadGaussNumPoints() ; ++k )
    {
      // -- evaluate the (gradient of the) shape functions
      elem->evalGradN( x0_e + u_e , k );

      // -- local displacement gradient
      gradv = elem->dNdx.transpose() * v_e;

      // -- local strain tensor
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = i ; j < ndim ; ++j )
          epsdot(i,j) = 0.5 * ( gradv(i,j) + gradv(j,i) );

      // -- local constitutive response: stiffness and stress tensors
      sig = mat->stressStrainRate(e,k,epsdot,elem->V);

      // -- add stress to element internal force column
      f_e += elem->gradN_tensor2(sig);
    }

    // - assemble internal force to global system
    for ( size_t i = 0 ; i < nne*ndim ; ++i ) Fdamp( dof_e(i) ) += f_e(i);
  }
}

// =================================================================================================

} // namespace DiagonalMass
} // namespace LinearStrain
} // namespace Periodic
} // namespace Dynamics
} // namespace GooseFEM

#endif

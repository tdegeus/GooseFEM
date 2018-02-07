/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_OVERDAMPEDDYNAMICS_DIAGONAL_CPP
#define GOOSEFEM_OVERDAMPEDDYNAMICS_DIAGONAL_CPP

// -------------------------------------------------------------------------------------------------

#include "OverdampedDynamicsDiagonal.h"

// ============================ GooseFEM::OverdampedDynamics::Diagonal =============================

namespace GooseFEM {
namespace OverdampedDynamics {
namespace Diagonal {

// ===================================== Simulation - Periodic =====================================

template<class Element>
inline Periodic<Element>::Periodic(
  std::unique_ptr<Element> _elem, const MatD &_x, const MatS &_conn, const MatS &_dofs, double _dt
)
: elem(std::move(_elem)), conn(_conn), dofs(_dofs), x(_x), dt(_dt)
{
  // compute sizes
  nnode = static_cast<size_t>(x.rows());
  ndim  = static_cast<size_t>(x.cols());
  nelem = static_cast<size_t>(conn.rows());
  nne   = static_cast<size_t>(conn.cols());
  ndof  = dofs.maxCoeff()+1;

  // basic checks (mostly the user is 'trusted')
  assert( static_cast<size_t>(dofs.size()) == nnode * ndim );
  assert(                     ndof         <  nnode * ndim );

  // allocate and zero-initialize nodal quantities
  u  .conservativeResize(nnode,ndim); u  .setZero();
  u_n.conservativeResize(nnode,ndim); u_n.setZero();
  v  .conservativeResize(nnode,ndim); v  .setZero();

  // allocate and zero-initialize linear system (DOFs)
  D   .conservativeResize(ndof);
  Dinv.conservativeResize(ndof);
  F   .conservativeResize(ndof);
  V   .conservativeResize(ndof);

  // initialize all fields
  updated_x();
  updated_u(true);
  updated_v(true);
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void Periodic<Element>::forwardEuler()
{
  // (1) compute the velocities
  // - solve for velocities (DOFs)
  V.noalias() = Dinv.cwiseProduct( - F );
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));

  // (2) update positions
  u += dt * v;

  // process update in displacements and velocities
  updated_u();
  updated_v();

  // update time
  t += dt;
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void Periodic<Element>::midpoint()
{
  // back-up history
  u_n = u;

  // (1) compute trial state
  // - solve for trial velocities (DOFs)
  V.noalias() = Dinv.cwiseProduct( - F );
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - update to trial positions
  u = u_n + 0.5 * dt * v;
  // - process update in displacements and velocities
  updated_u();
  updated_v();

  // (2) compute actual state
  // - solve for velocities (DOFs)
  V.noalias() = Dinv.cwiseProduct( - F );
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - update to trial positions
  u = u_n + dt * v;
  // - process update in displacements and velocities
  updated_u();
  updated_v();

  // update time
  t += dt;
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void Periodic<Element>::updated_x()
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

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void Periodic<Element>::updated_u(bool init)
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
  if ( elem->changed_D or init ) assemble_D();
  if ( elem->changed_f or init ) assemble_F();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void Periodic<Element>::updated_v(bool init)
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
  if ( elem->changed_D or init ) assemble_D();
  if ( elem->changed_f or init ) assemble_F();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void Periodic<Element>::assemble_D()
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

  // compute inverse
  Dinv = D.cwiseInverse();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void Periodic<Element>::assemble_F()
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

// =================================== Simulation - SemiPeriodic ===================================

template<class Element>
inline SemiPeriodic<Element>::SemiPeriodic(
  std::unique_ptr<Element> _elem, const MatD &_x, const MatS &_conn, const MatS &_dofs,
  const ColS &_fixedDofs, double _dt
)
: elem(std::move(_elem)), conn(_conn), dofs(_dofs), x(_x), fixedDofs(_fixedDofs), dt(_dt)
{
  // compute sizes
  nfixed    = static_cast<size_t>(fixedDofs.size());
  nnode     = static_cast<size_t>(x.rows());
  ndim      = static_cast<size_t>(x.cols());
  nelem     = static_cast<size_t>(conn.rows());
  nne       = static_cast<size_t>(conn.cols());
  ndof      = dofs.maxCoeff()+1;

  // basic checks (mostly the user is 'trusted')
  assert( static_cast<size_t>(dofs.size()) == nnode * ndim );
  assert(                     ndof         <  nnode * ndim );

  #ifndef NDEBUG
  for ( size_t i = 0 ; i < nfixed ; ++i ) assert( fixedDofs(i) < ndof );
  #endif

  // allocate and zero-initialize nodal quantities
  u  .conservativeResize(nnode,ndim); u  .setZero();
  u_n.conservativeResize(nnode,ndim); u_n.setZero();
  v  .conservativeResize(nnode,ndim); v  .setZero();

  // allocate and zero-initialize linear system (DOFs)
  D   .conservativeResize(ndof);
  Dinv.conservativeResize(ndof);
  F   .conservativeResize(ndof);
  V   .conservativeResize(ndof);

  // fixed DOFs : default zero velocity
  fixedV.conservativeResize(nfixed); fixedV.setZero();

  // initialize all fields
  updated_x();
  updated_u(true);
  updated_v(true);
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void SemiPeriodic<Element>::forwardEuler()
{
  // (1) compute the velocities
  // - solve for velocities (DOFs)
  V.noalias() = Dinv.cwiseProduct( - F );
  // - apply the fixed velocities
  for ( size_t i=0; i<nfixed; ++i ) V(fixedDofs(i)) = fixedV(i);
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));

  // (2) update positions
  u += dt * v;

  // process update in displacements and velocities
  updated_u();
  updated_v();

  // update time
  t += dt;
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void SemiPeriodic<Element>::midpoint()
{
  // back-up history
  u_n = u;

  // (1) compute trial state
  // - solve for trial velocities (DOFs)
  V.noalias() = Dinv.cwiseProduct( - F );
  // - apply the fixed velocities
  for ( size_t i=0; i<nfixed; ++i ) V(fixedDofs(i)) = fixedV(i);
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - update to trial positions
  u = u_n + 0.5 * dt * v;
  // - process update in displacements and velocities
  updated_u();
  updated_v();

  // (2) compute actual state
  // - solve for velocities (DOFs)
  V.noalias() = Dinv.cwiseProduct( - F );
  // - apply the fixed velocities
  for ( size_t i=0; i<nfixed; ++i ) V(fixedDofs(i)) = fixedV(i);
  // - convert to nodal velocities (periodicity implies that several nodes depend on the same DOF)
  for ( size_t i=0; i<nnode*ndim; ++i ) v(i) = V(dofs(i));
  // - update to trial positions
  u = u_n + dt * v;
  // - process update in displacements and velocities
  updated_u();
  updated_v();

  // update time
  t += dt;
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void SemiPeriodic<Element>::updated_x()
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

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void SemiPeriodic<Element>::updated_u(bool init)
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
  if ( elem->changed_D or init ) assemble_D();
  if ( elem->changed_f or init ) assemble_F();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void SemiPeriodic<Element>::updated_v(bool init)
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
  if ( elem->changed_D or init ) assemble_D();
  if ( elem->changed_f or init ) assemble_F();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void SemiPeriodic<Element>::assemble_D()
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

  // compute inverse
  Dinv = D.cwiseInverse();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
inline void SemiPeriodic<Element>::assemble_F()
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

// ======================================== Element - Quad4 ========================================

namespace SmallStrain {

// -------------------------------------------------------------------------------------------------

template<class Material>
inline Quad4<Material>::Quad4(std::unique_ptr<Material> _mat, size_t _nelem)
{
  // copy from input
  nelem = _nelem;
  mat   = std::move(_mat);

  // allocate matrices
  // -
  x     .resize({nelem,nne,ndim});
  u     .resize({nelem,nne,ndim});
  v     .resize({nelem,nne,ndim});
  f     .resize({nelem,nne,ndim});
  // -
  D     .resize({nelem,nne*ndim,nne*ndim});
  // -
  dNx   .resize({nelem,nip,nne,ndim});
  // -
  vol   .resize({nelem,nip});
  vol_n .resize({nelem,nne});
  // -
  dNxi  .resize({nip,nne,ndim});
  dNxi_n.resize({nne,nne,ndim});
  // -
  w     .resize({nip});
  w_n   .resize({nne});

  // zero initialize matrices (only diagonal is written, no new zero-initialization necessary)
  D.setZero();

  // shape function gradient at all Gauss points, in local coordinates
  // - k == 0
  dNxi(0,0,0) = -.25*(1.+1./std::sqrt(3.));     dNxi(0,0,1) = -.25*(1.+1./std::sqrt(3.));
  dNxi(0,1,0) = +.25*(1.+1./std::sqrt(3.));     dNxi(0,1,1) = -.25*(1.-1./std::sqrt(3.));
  dNxi(0,2,0) = +.25*(1.-1./std::sqrt(3.));     dNxi(0,2,1) = +.25*(1.-1./std::sqrt(3.));
  dNxi(0,3,0) = -.25*(1.-1./std::sqrt(3.));     dNxi(0,3,1) = +.25*(1.+1./std::sqrt(3.));
  // - k == 1
  dNxi(1,0,0) = -.25*(1.+1./std::sqrt(3.));     dNxi(1,0,1) = -.25*(1.-1./std::sqrt(3.));
  dNxi(1,1,0) = +.25*(1.+1./std::sqrt(3.));     dNxi(1,1,1) = -.25*(1.+1./std::sqrt(3.));
  dNxi(1,2,0) = +.25*(1.-1./std::sqrt(3.));     dNxi(1,2,1) = +.25*(1.+1./std::sqrt(3.));
  dNxi(1,3,0) = -.25*(1.-1./std::sqrt(3.));     dNxi(1,3,1) = +.25*(1.-1./std::sqrt(3.));
  // - k == 2
  dNxi(2,0,0) = -.25*(1.-1./std::sqrt(3.));     dNxi(2,0,1) = -.25*(1.-1./std::sqrt(3.));
  dNxi(2,1,0) = +.25*(1.-1./std::sqrt(3.));     dNxi(2,1,1) = -.25*(1.+1./std::sqrt(3.));
  dNxi(2,2,0) = +.25*(1.+1./std::sqrt(3.));     dNxi(2,2,1) = +.25*(1.+1./std::sqrt(3.));
  dNxi(2,3,0) = -.25*(1.+1./std::sqrt(3.));     dNxi(2,3,1) = +.25*(1.-1./std::sqrt(3.));
  // - k == 3
  dNxi(3,0,0) = -.25*(1.-1./std::sqrt(3.));     dNxi(3,0,1) = -.25*(1.+1./std::sqrt(3.));
  dNxi(3,1,0) = +.25*(1.-1./std::sqrt(3.));     dNxi(3,1,1) = -.25*(1.-1./std::sqrt(3.));
  dNxi(3,2,0) = +.25*(1.+1./std::sqrt(3.));     dNxi(3,2,1) = +.25*(1.-1./std::sqrt(3.));
  dNxi(3,3,0) = -.25*(1.+1./std::sqrt(3.));     dNxi(3,3,1) = +.25*(1.+1./std::sqrt(3.));

  // integration point weight at all Gauss points
  w(0) = 1.;
  w(1) = 1.;
  w(2) = 1.;
  w(3) = 1.;

  // shape function gradient at all nodes, in local coordinates
  // - k == 0
  dNxi_n(0,0,0) = -0.5;     dNxi_n(0,0,1) = -0.5;
  dNxi_n(0,1,0) = +0.5;     dNxi_n(0,1,1) =  0.0;
  dNxi_n(0,2,0) =  0.0;     dNxi_n(0,2,1) =  0.0;
  dNxi_n(0,3,0) =  0.0;     dNxi_n(0,3,1) = +0.5;
  // - k == 1
  dNxi_n(1,0,0) = -0.5;     dNxi_n(1,0,1) =  0.0;
  dNxi_n(1,1,0) = +0.5;     dNxi_n(1,1,1) = -0.5;
  dNxi_n(1,2,0) =  0.0;     dNxi_n(1,2,1) = +0.5;
  dNxi_n(1,3,0) =  0.0;     dNxi_n(1,3,1) =  0.0;
  // - k == 2
  dNxi_n(2,0,0) =  0.0;     dNxi_n(2,0,1) =  0.0;
  dNxi_n(2,1,0) =  0.0;     dNxi_n(2,1,1) = -0.5;
  dNxi_n(2,2,0) = +0.5;     dNxi_n(2,2,1) = +0.5;
  dNxi_n(2,3,0) = -0.5;     dNxi_n(2,3,1) =  0.0;
  // - k == 3
  dNxi_n(3,0,0) =  0.0;     dNxi_n(3,0,1) = -0.5;
  dNxi_n(3,1,0) =  0.0;     dNxi_n(3,1,1) =  0.0;
  dNxi_n(3,2,0) = +0.5;     dNxi_n(3,2,1) =  0.0;
  dNxi_n(3,3,0) = -0.5;     dNxi_n(3,3,1) = +0.5;

  // integration point weight at all nodes
  w_n(0) = 1.;
  w_n(1) = 1.;
  w_n(2) = 1.;
  w_n(3) = 1.;

  // Note, the above is a specialization of the following:
  //
  // - Shape function gradients
  //
  //    dNxi(0,0) = -.25*(1.-xi(k,1)); dNxi(0,1) = -.25*(1.-xi(k,0));
  //    dNxi(1,0) = +.25*(1.-xi(k,1)); dNxi(1,1) = -.25*(1.+xi(k,0));
  //    dNxi(2,0) = +.25*(1.+xi(k,1)); dNxi(2,1) = +.25*(1.+xi(k,0));
  //    dNxi(3,0) = -.25*(1.+xi(k,1)); dNxi(3,1) = +.25*(1.-xi(k,0));
  //
  // - Gauss point coordinates and weights
  //
  //    xi(0,0) = -1./std::sqrt(3.); xi(0,1) = -1./std::sqrt(3.); w(0) = 1.;
  //    xi(1,0) = +1./std::sqrt(3.); xi(1,1) = -1./std::sqrt(3.); w(1) = 1.;
  //    xi(2,0) = +1./std::sqrt(3.); xi(2,1) = +1./std::sqrt(3.); w(2) = 1.;
  //    xi(3,0) = -1./std::sqrt(3.); xi(3,1) = +1./std::sqrt(3.); w(3) = 1.;
  //
  // - Nodal coordinates and weights
  //
  //    xi(0,0) = -1.; xi(0,1) = -1.; w(0) = 1.;
  //    xi(1,0) = +1.; xi(1,1) = -1.; w(1) = 1.;
  //    xi(2,0) = +1.; xi(2,1) = +1.; w(2) = 1.;
  //    xi(3,0) = -1.; xi(3,1) = +1.; w(3) = 1.;

}

// =================================================================================================

template<class Material>
inline void Quad4<Material>::updated_x()
{
#pragma omp parallel
{
  // intermediate quantities
  cppmat::cartesian2d::tensor2<double> J_, Jinv_;
  double Jdet_;
  // local views of the global arrays (speeds up indexing, and increases readability)
  cppmat::tiny::matrix2<double,8,8> M_, D_;
  cppmat::tiny::matrix2<double,4,2> dNxi_, dNx_, x_;
  cppmat::tiny::vector <double,4>   w_, vol_;

  // loop over all elements (in parallel)
  #pragma omp for
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // pointer to element positions
    x_.map(&x(e));

    // nodal quadrature
    // ----------------

    // pointer to element damping matrix, nodal volume, and integration weight
    D_  .map(&D    (e));
    vol_.map(&vol_n(e));
    w_  .map(&w_n  (0));

    // loop over nodes, i.e. the integration points
    for ( size_t k = 0 ; k < nne ; ++k )
    {
      // - pointer to the shape function gradients (local coordinates)
      dNxi_.map(&dNxi_n(k));

      // - Jacobian
      //   J(i,j) += dNxi(m,i) * xe(m,j)
      J_(0,0) = dNxi_(0,0)*x_(0,0) + dNxi_(1,0)*x_(1,0) + dNxi_(2,0)*x_(2,0) + dNxi_(3,0)*x_(3,0);
      J_(0,1) = dNxi_(0,0)*x_(0,1) + dNxi_(1,0)*x_(1,1) + dNxi_(2,0)*x_(2,1) + dNxi_(3,0)*x_(3,1);
      J_(1,0) = dNxi_(0,1)*x_(0,0) + dNxi_(1,1)*x_(1,0) + dNxi_(2,1)*x_(2,0) + dNxi_(3,1)*x_(3,0);
      J_(1,1) = dNxi_(0,1)*x_(0,1) + dNxi_(1,1)*x_(1,1) + dNxi_(2,1)*x_(2,1) + dNxi_(3,1)*x_(3,1);

      // - determinant of the Jacobian
      Jdet_ = J_.det();

      // - integration point volume
      vol_(k) = w_(k) * Jdet_;

      // - assemble element non-Galilean damping matrix
      //   D(m+i,n+i) += N(m) * alpha * vol * N(n);
      D_(k*2  ,k*2  ) = mat->alpha(e,k) * vol_(k);
      D_(k*2+1,k*2+1) = mat->alpha(e,k) * vol_(k);
    }

    // Gaussian quadrature
    // -------------------

    // pointer to element integration volume and weight
    vol_.map(&vol(e));
    w_  .map(&w  (0));

    // loop over Gauss points
    for ( size_t k = 0 ; k < nip ; ++k )
    {
      // - pointer to the shape function gradients (local/global coordinates)
      dNxi_.map(&dNxi(  k));
      dNx_ .map(&dNx (e,k));

      // - Jacobian
      //   J(i,j) += dNxi(m,i) * xe(m,j)
      J_(0,0) = dNxi_(0,0)*x_(0,0) + dNxi_(1,0)*x_(1,0) + dNxi_(2,0)*x_(2,0) + dNxi_(3,0)*x_(3,0);
      J_(0,1) = dNxi_(0,0)*x_(0,1) + dNxi_(1,0)*x_(1,1) + dNxi_(2,0)*x_(2,1) + dNxi_(3,0)*x_(3,1);
      J_(1,0) = dNxi_(0,1)*x_(0,0) + dNxi_(1,1)*x_(1,0) + dNxi_(2,1)*x_(2,0) + dNxi_(3,1)*x_(3,0);
      J_(1,1) = dNxi_(0,1)*x_(0,1) + dNxi_(1,1)*x_(1,1) + dNxi_(2,1)*x_(2,1) + dNxi_(3,1)*x_(3,1);

      // - determinant and inverse of the Jacobian
      Jdet_ = J_.det();
      Jinv_ = J_.inv();

      // - integration point volume
      vol_(k) = w_(k) * Jdet_;

      // - shape function gradients (global coordinates)
      //   dNx(m,i) += Jinv(i,j) * dNxi(m,j)
      for ( size_t m = 0 ; m < nne ; ++m )
      {
        dNx_(m,0) = Jinv_(0,0) * dNxi_(m,0) + Jinv_(0,1) * dNxi_(m,1);
        dNx_(m,1) = Jinv_(1,0) * dNxi_(m,0) + Jinv_(1,1) * dNxi_(m,1);
      }
    }
  }

  // set signals
  changed_f = false;
  changed_D = true;

} // #pragma omp parallel
}

// =================================================================================================

template<class Material>
inline void Quad4<Material>::updated_u()
{
#pragma omp parallel
{
  // intermediate quantities
  cppmat::cartesian2d::tensor2 <double> gradu_;
  cppmat::cartesian2d::tensor2s<double> eps_, sig_;
  // local views of the global arrays (speeds up indexing, and increases readability)
  cppmat::tiny::matrix2<double,4,2> dNx_, u_, f_;
  cppmat::tiny::vector <double,4>   vol_;

  // loop over all elements (in parallel)
  #pragma omp for
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // pointer to element forces, displacements, and integration volume
    f_  .map(&f  (e));
    u_  .map(&u  (e));
    vol_.map(&vol(e));

    // zero initialize forces
    f_.setZero();

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < nip ; ++k )
    {
      // - pointer to the shape function gradients, strain and stress tensor (stored symmetric)
      dNx_.map(&dNx     (e,k));
      eps_.map(&mat->eps(e,k));
      sig_.map(&mat->sig(e,k));

      // - displacement gradient
      //   gradu_(i,j) += dNx(m,i) * ue(m,j)
      gradu_(0,0) = dNx_(0,0)*u_(0,0) + dNx_(1,0)*u_(1,0) + dNx_(2,0)*u_(2,0) + dNx_(3,0)*u_(3,0);
      gradu_(0,1) = dNx_(0,0)*u_(0,1) + dNx_(1,0)*u_(1,1) + dNx_(2,0)*u_(2,1) + dNx_(3,0)*u_(3,1);
      gradu_(1,0) = dNx_(0,1)*u_(0,0) + dNx_(1,1)*u_(1,0) + dNx_(2,1)*u_(2,0) + dNx_(3,1)*u_(3,0);
      gradu_(1,1) = dNx_(0,1)*u_(0,1) + dNx_(1,1)*u_(1,1) + dNx_(2,1)*u_(2,1) + dNx_(3,1)*u_(3,1);

      // - strain (stored symmetric)
      //   eps(i,j) = .5 * ( gradu_(i,j) + gradu_(j,i) )
      eps_(0,0) =        gradu_(0,0);
      eps_(0,1) = .5 * ( gradu_(0,1) + gradu_(1,0) );
      eps_(1,1) =        gradu_(1,1);

      // - constitutive response
      mat->updated_eps(e,k);

      // - assemble to element force
      //   f(m,j) += dNx(m,i) * sig(i,j) * vol;
      for ( size_t m = 0 ; m < nne ; ++m )
      {
        f_(m,0) += dNx_(m,0) * sig_(0,0) * vol_(k) + dNx_(m,1) * sig_(1,0) * vol_(k);
        f_(m,1) += dNx_(m,0) * sig_(0,1) * vol_(k) + dNx_(m,1) * sig_(1,1) * vol_(k);
      }
    }
  }

  // set signals
  changed_f = true;
  changed_D = false;
}
}

// =================================================================================================

template<class Material>
inline void Quad4<Material>::updated_v()
{
  changed_f = true;
  changed_D = false;
}

// =================================================================================================

template<class Material>
inline cppmat::cartesian2d::tensor2s<double> Quad4<Material>::mean_eps(size_t e)
{
  cppmat::cartesian2d::tensor2s<double> eps_, tot_eps_(0.0);
  cppmat::tiny::vector <double,4> vol_;
  double tot_vol_ = 0.0;

  // pointer to integration volume
  vol_.map(&vol(e));

  // loop over all integration points in element "e"
  for ( size_t k = 0 ; k < nip ; ++k )
  {
    // - pointer to strain tensor (stored symmetric)
    eps_.map(&mat->eps(e,k));

    // - add to average
    tot_eps_ += vol_(k) * eps_;
    tot_vol_ += vol_(k);
  }

  // return volume average
  return ( tot_eps_ / tot_vol_ );
}

// =================================================================================================

template<class Material>
inline cppmat::cartesian2d::tensor2s<double> Quad4<Material>::mean_sig(size_t e)
{
  cppmat::cartesian2d::tensor2s<double> sig_, tot_sig_(0.0);
  cppmat::tiny::vector <double,4> vol_;
  double tot_vol_ = 0.0;

  // pointer to integration volume
  vol_.map(&vol(e));

  // loop over all integration points in element "e"
  for ( size_t k = 0 ; k < nip ; ++k )
  {
    // - pointer to strain tensor (stored symmetric)
    sig_.map(&mat->sig(e,k));

    // - add to average
    tot_sig_ += vol_(k) * sig_;
    tot_vol_ += vol_(k);
  }

  // return volume average
  return ( tot_sig_ / tot_vol_ );
}

// =================================================================================================

template<class Material>
inline cppmat::cartesian2d::tensor2s<double> Quad4<Material>::mean_eps()
{
  cppmat::cartesian2d::tensor2s<double> eps_, tot_eps_(0.0);
  cppmat::tiny::vector <double,4> vol_;
  double tot_vol_ = 0.0;

  // loop over all elements (in parallel)
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // pointer to integration volume
    vol_.map(&vol(e));

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < nip ; ++k )
    {
      // - pointer to strain tensor (stored symmetric)
      eps_.map(&mat->eps(e,k));

      // - add to average
      tot_eps_ += vol_(k) * eps_;
      tot_vol_ += vol_(k);
    }
  }

  // return volume average
  return ( tot_eps_ / tot_vol_ );
}

// =================================================================================================

template<class Material>
inline cppmat::cartesian2d::tensor2s<double> Quad4<Material>::mean_sig()
{
  cppmat::cartesian2d::tensor2s<double> sig_, tot_sig_(0.0);
  cppmat::tiny::vector <double,4> vol_;
  double tot_vol_ = 0.0;

  // loop over all elements (in parallel)
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // pointer to integration volume
    vol_.map(&vol(e));

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < nip ; ++k )
    {
      // - pointer to strain tensor (stored symmetric)
      sig_.map(&mat->sig(e,k));

      // - add to average
      tot_sig_ += vol_(k) * sig_;
      tot_vol_ += vol_(k);
    }
  }

  // return volume average
  return ( tot_sig_ / tot_vol_ );
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

}}} // namespace ...

// =================================================================================================

#endif

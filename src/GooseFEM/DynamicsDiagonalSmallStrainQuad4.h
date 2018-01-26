/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_DYNAMICS_DIAGONAL_SMALLSTRAIN_QUAD4_H
#define GOOSEFEM_DYNAMICS_DIAGONAL_SMALLSTRAIN_QUAD4_H

// -------------------------------------------------------------------------------------------------

#include "Macros.h"
#include <cppmat/cppmat.h>

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Element {
namespace Diagonal {
namespace SmallStrain {

// ========================== N.B. most loops are unrolled for efficiency ==========================

template<class Material>
class Quad4
{
public:

  // class variables
  // ---------------

  // constitutive response per integration point, per element
  std::unique_ptr<Material> mat;

  // matrices to store the element data
  cppmat::matrix<double> x, u, v, f, M, D, dNx, dNxi, w, vol, dNxi_n, w_n, vol_n;

  // dimensions
  size_t nelem, nip=4, nne=4, ndim=2;

  // signals for solvers
  bool changed_f=false, changed_M=true, changed_D=true;

  // class functions
  // ---------------

  // constructor
  Quad4(std::unique_ptr<Material> mat, size_t nelem);

  // recompute relevant quantities after "x", "u", or "v" have been externally updated
  void updated_x();
  void updated_u();
  void updated_v();

  // return volume average stress and strain
  cppmat::cartesian2d::tensor2s<double> mean_eps(size_t e);   // of element "e"
  cppmat::cartesian2d::tensor2s<double> mean_sig(size_t e);   // of element "e"
  cppmat::cartesian2d::tensor2s<double> mean_eps();           // of all elements
  cppmat::cartesian2d::tensor2s<double> mean_sig();           // of all elements
};

// =================================================================================================

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
  M     .resize({nelem,nne*ndim,nne*ndim});
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
  M.setZero();
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

    // pointer to element mass/damping matrix, nodal volume, and integration weight
    M_  .map(&M    (e));
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

      // - assemble element mass matrix
      //   M(m+i,n+i) += N(m) * rho * vol * N(n);
      M_(k*2  ,k*2  ) = mat->rho(e,k) * vol_(k);
      M_(k*2+1,k*2+1) = mat->rho(e,k) * vol_(k);

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
  changed_M = true;
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
  changed_M = false;
  changed_D = false;
}
}

// =================================================================================================

template<class Material>
inline void Quad4<Material>::updated_v()
{
#pragma omp parallel
{
  // intermediate quantities
  cppmat::cartesian2d::tensor2 <double> gradv_;
  cppmat::cartesian2d::tensor2s<double> epsdot_, sig_;
  // local views of the global arrays (speeds up indexing, and increases readability)
  cppmat::tiny::matrix2<double,4,2> dNx_, v_, f_;
  cppmat::tiny::vector <double,4>   vol_;

  // loop over all elements (in parallel)
  #pragma omp for
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // pointer to element forces, displacements, and integration volume
    f_  .map(&f  (e));
    v_  .map(&v  (e));
    vol_.map(&vol(e));

    // zero initialize forces
    f_.setZero();

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < nip ; ++k )
    {
      // - pointer to the shape function gradients, strain-rate and stress tensor (stored symmetric)
      dNx_   .map(&dNx        (e,k));
      epsdot_.map(&mat->epsdot(e,k));
      sig_   .map(&mat->sig   (e,k));

      // - displacement gradient
      //   gradv_(i,j) += dNx(m,i) * ue(m,j)
      gradv_(0,0) = dNx_(0,0)*v_(0,0) + dNx_(1,0)*v_(1,0) + dNx_(2,0)*v_(2,0) + dNx_(3,0)*v_(3,0);
      gradv_(0,1) = dNx_(0,0)*v_(0,1) + dNx_(1,0)*v_(1,1) + dNx_(2,0)*v_(2,1) + dNx_(3,0)*v_(3,1);
      gradv_(1,0) = dNx_(0,1)*v_(0,0) + dNx_(1,1)*v_(1,0) + dNx_(2,1)*v_(2,0) + dNx_(3,1)*v_(3,0);
      gradv_(1,1) = dNx_(0,1)*v_(0,1) + dNx_(1,1)*v_(1,1) + dNx_(2,1)*v_(2,1) + dNx_(3,1)*v_(3,1);

      // - strain (stored symmetric)
      //   epsdot(i,j) = .5 * ( gradv_(i,j) + gradv_(j,i) )
      epsdot_(0,0) =        gradv_(0,0);
      epsdot_(0,1) = .5 * ( gradv_(0,1) + gradv_(1,0) );
      epsdot_(1,1) =        gradv_(1,1);

      // - constitutive response
      mat->updated_epsdot(e,k);

      // - assemble to element force
      //   f(m,j) += dNdx(m,i) * sig(i,j) * vol;
      for ( size_t m = 0 ; m < nne ; ++m )
      {
        f_(m,0) += dNx_(m,0) * sig_(0,0) * vol_(k) + dNx_(m,1) * sig_(1,0) * vol_(k);
        f_(m,1) += dNx_(m,0) * sig_(0,1) * vol_(k) + dNx_(m,1) * sig_(1,1) * vol_(k);
      }
    }
  }

  // set signals
  changed_f = true;
  changed_M = false;
  changed_D = false;
}
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

// =================================================================================================

}}}} // namespace ...

#endif

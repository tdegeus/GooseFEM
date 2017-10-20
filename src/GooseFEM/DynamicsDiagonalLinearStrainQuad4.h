/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_DYNAMICS_DIAGONAL_QUAD4_H
#define GOOSEFEM_DYNAMICS_DIAGONAL_QUAD4_H

#include "Macros.h"
#include <cppmat/cppmat.h>

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Dynamics {
namespace Diagonal {
namespace LinearStrain {

// -------------------------------------------------------------------------------------------------

using vec = cppmat::cartesian2d::vector  <double>;
using T2  = cppmat::cartesian2d::tensor2 <double>;
using T2s = cppmat::cartesian2d::tensor2s<double>;
using T2d = cppmat::cartesian2d::tensor2d<double>;

// ========================== N.B. most loops are unrolled for efficiency ==========================

template <class QuadraturePoint>
class Quad4
{
public:

  // variables
  // ---------

  // arrays / matrices
  cppmat::tiny::matrix2<double,4,2> xe, ue, ve, xi, xi_n, dNdxi, dNdx;
  cppmat::tiny::vector <double,4>   w, w_n;
  cppmat::tiny::matrix2<double,8,8> M;
  cppmat::tiny::vector <double,8>   fu, fv;
  // tensors
  cppmat::cartesian2d::tensor2<double> J, Jinv, gradu, gradv;
  // scalars
  double Jdet, V;
  // sizes (nodes per element, dimensions, quadrature points)
  size_t nne=4, ndim=2, nk=4;
  // quadrature point : provides the constitutive response
  std::shared_ptr<QuadraturePoint> quad;

  // constructor
  // -----------

  Quad4(std::shared_ptr<QuadraturePoint> quad);

  // functions
  // ---------

  void computeM (size_t elem); // mass matrix                     <- quad->density
  void computeFu(size_t elem); // displacement dependent forces   <- quad->stressStrain
  void computeFv(size_t elem); // displacement dependent forces   <- quad->stressStrainRate
  void post     (size_t elem); // post-process                    <- quad->stressStrain(Rate)

};

// =================================================================================================

template <class QuadraturePoint>
Quad4<QuadraturePoint>::Quad4(std::shared_ptr<QuadraturePoint> _quad)
{
  // quadrature point routines
  quad = _quad;

  // integration point coordinates/weights: normal Gauss integration
  xi(0,0) = -1./std::sqrt(3.); xi(0,1) = -1./std::sqrt(3.); w(0) = 1.;
  xi(1,0) = +1./std::sqrt(3.); xi(1,1) = -1./std::sqrt(3.); w(1) = 1.;
  xi(2,0) = +1./std::sqrt(3.); xi(2,1) = +1./std::sqrt(3.); w(2) = 1.;
  xi(3,0) = -1./std::sqrt(3.); xi(3,1) = +1./std::sqrt(3.); w(3) = 1.;

  // integration point coordinates/weights: integration at the nodes
  xi_n(0,0) = -1.; xi_n(0,1) = -1.; w_n(0) = 1.;
  xi_n(1,0) = +1.; xi_n(1,1) = -1.; w_n(1) = 1.;
  xi_n(2,0) = +1.; xi_n(2,1) = +1.; w_n(2) = 1.;
  xi_n(3,0) = -1.; xi_n(3,1) = +1.; w_n(3) = 1.;
}

// =================================================================================================

template <class QuadraturePoint>
void Quad4<QuadraturePoint>::computeM(size_t elem)
{
  // zero-initialize element mass matrix
  M.zeros();

  // loop over integration points (coincide with the nodes to get a diagonal mass matrix)
  for ( size_t k = 0 ; k < nne ; ++k )
  {
    // - shape function gradients (local coordinates)
    dNdxi(0,0) = -.25*(1.-xi_n(k,1)); dNdxi(0,1) = -.25*(1.-xi_n(k,0));
    dNdxi(1,0) = +.25*(1.-xi_n(k,1)); dNdxi(1,1) = -.25*(1.+xi_n(k,0));
    dNdxi(2,0) = +.25*(1.+xi_n(k,1)); dNdxi(2,1) = +.25*(1.+xi_n(k,0));
    dNdxi(3,0) = -.25*(1.+xi_n(k,1)); dNdxi(3,1) = +.25*(1.-xi_n(k,0));

    // - Jacobian
    //   J(i,j) += dNdxi(m,i) * xe(m,j)
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m ) {
      J(0,0) += dNdxi(m,0) * xe(m,0);
      J(0,1) += dNdxi(m,0) * xe(m,1);
      J(1,0) += dNdxi(m,1) * xe(m,0);
      J(1,1) += dNdxi(m,1) * xe(m,1);
    }

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume (== part of the element volume associated with the node)
    V = w_n(k) * Jdet;

    // - assemble to element mass matrix (use the delta properties of the shape functions)
    //   M(m+i,n+i) = N(m) * rho * V * N(n);
    M(k*2  ,k*2  ) = quad->density(elem,k,V) * V;
    M(k*2+1,k*2+1) = quad->density(elem,k,V) * V;
  }
}

// =================================================================================================

template <class QuadraturePoint>
void Quad4<QuadraturePoint>::computeFu(size_t elem)
{
  // zero-initialize element force
  fu.zeros();

  // loop over integration points
  for ( size_t k = 0 ; k < nk ; ++k )
  {
    // - shape function gradients (local coordinates)
    dNdxi(0,0) = -.25*(1.-xi(k,1)); dNdxi(0,1) = -.25*(1.-xi(k,0));
    dNdxi(1,0) = +.25*(1.-xi(k,1)); dNdxi(1,1) = -.25*(1.+xi(k,0));
    dNdxi(2,0) = +.25*(1.+xi(k,1)); dNdxi(2,1) = +.25*(1.+xi(k,0));
    dNdxi(3,0) = -.25*(1.+xi(k,1)); dNdxi(3,1) = +.25*(1.-xi(k,0));

    // - Jacobian
    //   J(i,j) += dNdxi(m,i) * xe(m,j)
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m ) {
      J(0,0) += dNdxi(m,0) * xe(m,0);
      J(0,1) += dNdxi(m,0) * xe(m,1);
      J(1,0) += dNdxi(m,1) * xe(m,0);
      J(1,1) += dNdxi(m,1) * xe(m,1);
    }

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;

    // - shape function gradients (global coordinates)
    //   dNdx(m,i) += Jinv(i,j) * dNdxi(m,j)
    for ( size_t m = 0 ; m < nne ; ++m ) {
      dNdx(m,0) = Jinv(0,0) * dNdxi(m,0) + Jinv(0,1) * dNdxi(m,1);
      dNdx(m,1) = Jinv(1,0) * dNdxi(m,0) + Jinv(1,1) * dNdxi(m,1);
    }

    // - displacement gradient
    //   gradu(i,j) += dNdx(m,i) * ue(m,j)
    gradu.zeros();
    for ( size_t m = 0 ; m < nne ; ++m ) {
      gradu(0,0) += dNdx(m,0) * ue(m,0);
      gradu(0,1) += dNdx(m,0) * ue(m,1);
      gradu(1,0) += dNdx(m,1) * ue(m,0);
      gradu(1,1) += dNdx(m,1) * ue(m,1);
    }

    // - strain tensor (symmetric part of "gradu")
    //   quad->eps(i,j) = .5 * ( gradu(i,j) + gradu(j,i) )
    quad->eps(0,0) = gradu(0,0);
    quad->eps(0,1) = .5 * ( gradu(0,1) + gradu(1,0) );
    quad->eps(1,0) = quad->eps(0,1);
    quad->eps(1,1) = gradu(1,1);

    // - constitutive response
    quad->stressStrain(elem,k,V);

    // - assemble to element force
    //   fu(m*ndim+j) += dNdx(m,i) * quad->sig(i,j) * V;
    for ( size_t m = 0 ; m < nne ; ++m ) {
      fu(m*ndim+0) += dNdx(m,0) * quad->sig(0,0) * V + dNdx(m,1) * quad->sig(1,0) * V;
      fu(m*ndim+1) += dNdx(m,0) * quad->sig(0,1) * V + dNdx(m,1) * quad->sig(1,1) * V;
    }
  }
}

// =================================================================================================

template <class QuadraturePoint>
void Quad4<QuadraturePoint>::computeFv(size_t elem)
{
  // zero-initialize element force
  fv.zeros();

  // loop over integration points
  for ( size_t k = 0 ; k < nk ; ++k )
  {
    // - shape function gradients (local coordinates)
    dNdxi(0,0) = -.25*(1.-xi(k,1)); dNdxi(0,1) = -.25*(1.-xi(k,0));
    dNdxi(1,0) = +.25*(1.-xi(k,1)); dNdxi(1,1) = -.25*(1.+xi(k,0));
    dNdxi(2,0) = +.25*(1.+xi(k,1)); dNdxi(2,1) = +.25*(1.+xi(k,0));
    dNdxi(3,0) = -.25*(1.+xi(k,1)); dNdxi(3,1) = +.25*(1.-xi(k,0));

    // - Jacobian
    //   J(i,j) += dNdxi(m,i) * xe(m,j)
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m ) {
      J(0,0) += dNdxi(m,0) * xe(m,0);
      J(0,1) += dNdxi(m,0) * xe(m,1);
      J(1,0) += dNdxi(m,1) * xe(m,0);
      J(1,1) += dNdxi(m,1) * xe(m,1);
    }

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;

    // - shape function gradients (global coordinates)
		//   dNdx(m,i) += Jinv(i,j) * dNdxi(m,j)
    for ( size_t m = 0 ; m < nne ; ++m ) {
      dNdx(m,0) = Jinv(0,0) * dNdxi(m,0) + Jinv(0,1) * dNdxi(m,1);
      dNdx(m,1) = Jinv(1,0) * dNdxi(m,0) + Jinv(1,1) * dNdxi(m,1);
    }

    // - velocity gradient
		//   gradv(i,j) += dNdx(m,i) * ve(m,j)
    gradv.zeros();
    for ( size_t m = 0 ; m < nne ; ++m ) {
      gradv(0,0) += dNdx(m,0) * ve(m,0);
      gradv(0,1) += dNdx(m,0) * ve(m,1);
      gradv(1,0) += dNdx(m,1) * ve(m,0);
      gradv(1,1) += dNdx(m,1) * ve(m,1);
    }

    // - strain-rate tensor (symmetric part of "gradv")
		//   quad->epsdot(i,j) = .5 * ( gradv(i,j) + gradv(j,i) )
    quad->epsdot(0,0) = gradv(0,0);
    quad->epsdot(0,1) = .5 * ( gradv(0,1) + gradv(1,0) );
    quad->epsdot(1,0) = quad->epsdot(0,1);
    quad->epsdot(1,1) = gradv(1,1);

    // - constitutive response
    quad->stressStrainRate(elem,k,V);

    // - assemble to element force
		//   fv(m*ndim+j) += dNdx(m,i) * quad->sig(i,j) * V;
    for ( size_t m = 0 ; m < nne ; ++m ) {
      fv(m*2+0) += dNdx(m,0) * quad->sig(0,0) * V + dNdx(m,1) * quad->sig(1,0) * V;
      fv(m*2+1) += dNdx(m,0) * quad->sig(0,1) * V + dNdx(m,1) * quad->sig(1,1) * V;
    }
  }
}

// =================================================================================================

template <class QuadraturePoint>
void Quad4<QuadraturePoint>::post(size_t elem)
{
  // loop over integration points
  for ( size_t k = 0 ; k < nk ; ++k )
  {
    // - shape function gradients (local coordinates)
    dNdxi(0,0) = -.25*(1.-xi(k,1)); dNdxi(0,1) = -.25*(1.-xi(k,0));
    dNdxi(1,0) = +.25*(1.-xi(k,1)); dNdxi(1,1) = -.25*(1.+xi(k,0));
    dNdxi(2,0) = +.25*(1.+xi(k,1)); dNdxi(2,1) = +.25*(1.+xi(k,0));
    dNdxi(3,0) = -.25*(1.+xi(k,1)); dNdxi(3,1) = +.25*(1.-xi(k,0));

    // - Jacobian
    //   J(i,j) += dNdxi(m,i) * xe(m,j)
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m ) {
      J(0,0) += dNdxi(m,0) * xe(m,0);
      J(0,1) += dNdxi(m,0) * xe(m,1);
      J(1,0) += dNdxi(m,1) * xe(m,0);
      J(1,1) += dNdxi(m,1) * xe(m,1);
    }

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;

    // - shape function gradients (global coordinates)
    //   dNdx(m,i) += Jinv(i,j) * dNdxi(m,j)
    for ( size_t m = 0 ; m < nne ; ++m ) {
      dNdx(m,0) = Jinv(0,0) * dNdxi(m,0) + Jinv(0,1) * dNdxi(m,1);
      dNdx(m,1) = Jinv(1,0) * dNdxi(m,0) + Jinv(1,1) * dNdxi(m,1);
    }

    // - displacement gradient
    //   gradu(i,j) += dNdx(m,i) * ue(m,j)
    gradu.zeros();
    for ( size_t m = 0 ; m < nne ; ++m ) {
      gradu(0,0) += dNdx(m,0) * ue(m,0);
      gradu(0,1) += dNdx(m,0) * ue(m,1);
      gradu(1,0) += dNdx(m,1) * ue(m,0);
      gradu(1,1) += dNdx(m,1) * ue(m,1);
    }

    // - strain tensor (symmetric part of "gradu")
    //   quad->eps(i,j) = .5 * ( gradu(i,j) + gradu(j,i) )
    quad->eps(0,0) = gradu(0,0);
    quad->eps(0,1) = .5 * ( gradu(0,1) + gradu(1,0) );
    quad->eps(1,0) = quad->eps(0,1);
    quad->eps(1,1) = gradu(1,1);

    // - velocity gradient
    //   gradv(i,j) += dNdx(m,i) * ve(m,j)
    gradv.zeros();
    for ( size_t m = 0 ; m < nne ; ++m ) {
      gradv(0,0) += dNdx(m,0) * ve(m,0);
      gradv(0,1) += dNdx(m,0) * ve(m,1);
      gradv(1,0) += dNdx(m,1) * ve(m,0);
      gradv(1,1) += dNdx(m,1) * ve(m,1);
    }

    // - strain-rate tensor (symmetric part of "gradv")
    //   quad->epsdot(i,j) = .5 * ( gradv(i,j) + gradv(j,i) )
    quad->epsdot(0,0) = gradv(0,0);
    quad->epsdot(0,1) = .5 * ( gradv(0,1) + gradv(1,0) );
    quad->epsdot(1,0) = quad->epsdot(0,1);
    quad->epsdot(1,1) = gradv(1,1);

    // - constitutive response
    quad->stressStrainPost    (elem,k,V);
    quad->stressStrainRatePost(elem,k,V);
  }
}

// =================================================================================================

}}}} // namespace GooseFEM::Dynamics::Diagonal::LinearStrain

#endif

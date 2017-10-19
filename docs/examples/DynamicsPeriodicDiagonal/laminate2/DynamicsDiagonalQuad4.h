/* ========================================== DESCRIPTION ==========================================

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

// -------------------------------------------------------------------------------------------------

using vec = cppmat::cartesian2d::vector  <double>;
using T2  = cppmat::cartesian2d::tensor2 <double>;
using T2s = cppmat::cartesian2d::tensor2s<double>;
using T2d = cppmat::cartesian2d::tensor2d<double>;

// =================================================================================================

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
  // quadrature point
  std::shared_ptr<QuadraturePoint> quad;

  // constructor
  // -----------

  Quad4(std::shared_ptr<QuadraturePoint> quad);

  // functions
  // ---------

  void computeMassMatrix   (size_t elem);
  void computeInternalForce(size_t elem);
  void computeDampingForce (size_t elem);
  void postProcess         (size_t elem);

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
void Quad4<QuadraturePoint>::computeMassMatrix(size_t elem)
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
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume (== volume associated with the node, in this element)
    V = w_n(k) * Jdet;

    // - assemble to element mass matrix (use the delta properties of the shape functions)
    for ( size_t i = 0 ; i < ndim ; ++i )
      M(i,i) += quad->density(elem,k,V) * V;
  }
}

// =================================================================================================

template <class QuadraturePoint>
void Quad4<QuadraturePoint>::computeInternalForce(size_t elem)
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
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;

    // - shape function gradients (global coordinates)
    dNdx.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          dNdx(m,i) += Jinv(i,j) * dNdxi(m,j);

    // - displacement gradient
    gradu.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradu(i,j) += dNdx(m,i) * ue(m,j);

    // - strain tensor (symmetric part of "gradu")
    for ( size_t i = 0 ; i < ndim ; ++i ) {
      for ( size_t j = i ; j < ndim ; ++j ) {
        quad->eps(i,j) = .5 * ( gradu(i,j) + gradu(j,i) );
        quad->eps(j,i) = quad->eps(i,j);
      }
    }

    // - constitutive response
    quad->stressStrain(elem,k,V)

    // - assemble to element force
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          fu(m*ndim+j) += dNdx(m,i) * quad->sig(i,j) * V;
  }
}

// =================================================================================================

template <class QuadraturePoint>
void Quad4<QuadraturePoint>::computeDampingForce(size_t elem)
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
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;

    // - shape function gradients (global coordinates)
    dNdx.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          dNdx(m,i) += Jinv(i,j) * dNdxi(m,j);

    // - velocity gradient
    gradv.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradv(i,j) += dNdx(m,i) * ve(m,j);

    // - strain-rate tensor (symmetric part of "gradu")
    for ( size_t i = 0 ; i < ndim ; ++i ) {
      for ( size_t j = i ; j < ndim ; ++j ) {
        quad->epsdot(i,j) = .5 * ( gradu(i,j) + gradu(j,i) );
        quad->epsdot(j,i) = quad->epsdot(i,j);
      }
    }

    // - constitutive response
    quad->stressStrainRate(elem,k,V)

    // - assemble to element force
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          fv(m*ndim+j) += dNdx(m,i) * quad->sig(i,j) * V;
  }
}

// =================================================================================================

template <class QuadraturePoint>
void Quad4<QuadraturePoint>::postProcess(size_t elem)
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
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;
    // - add to total volume
    Vbar += V;

    // - shape function gradients (global coordinates)
    dNdx.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          dNdx(m,i) += Jinv(i,j) * dNdxi(m,j);

    // - displacement gradient
    gradu.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradu(i,j) += dNdx(m,i) * ue(m,j);

    // - strain tensor (symmetric part of "gradu")
    for ( size_t i = 0 ; i < ndim ; ++i ) {
      for ( size_t j = i ; j < ndim ; ++j ) {
        quad->eps(i,j) = .5 * ( gradu(i,j) + gradu(j,i) );
        quad->eps(j,i) = quad->eps(i,j);
      }
    }

    // - velocity gradient
    gradv.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradv(i,j) += dNdx(m,i) * ve(m,j);

    // - strain-rate tensor (symmetric part of "gradu")
    for ( size_t i = 0 ; i < ndim ; ++i ) {
      for ( size_t j = i ; j < ndim ; ++j ) {
        quad->epsdot(i,j) = .5 * ( gradu(i,j) + gradu(j,i) );
        quad->epsdot(j,i) = quad->epsdot(i,j);
      }
    }

    // - constitutive response
    quad->stressStrainPost(elem,k,V)
    quad->stressStrainRatePost(elem,k,V)
  }
}

// =================================================================================================

}}} // namespace GooseFEM::Dynamics::Diagonal

#endif

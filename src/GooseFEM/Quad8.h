/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_QUAD8_H
#define GOOSEFEM_QUAD8_H

#include "Macros.h"

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {

// ======================================= ELEMENT DEFINITION ======================================

class Quad8
{
public:
  size_t   nne   = 8; // number of nodes per element
  size_t   ndim  = 2; // number of dimensions
  double   thick    ; // thickness (modify to your needs)
  double   V        ; // volume
  ColD     x        ; // global coordinates                            [ndim     ]
  MatD     J        ; // Jacobian                                      [ndim,ndim]
  ColD     N        ; // shape functions                               [nne      ]
  MatD     dNdxi    ; // shape function gradients (local  coordinates) [nne ,ndim]
  MatD     dNdx     ; // shape function gradients (global coordinates) [nne ,ndim]

  // constructor
  Quad8(double thickness=1.);

  // integration point definitions: Gauss quadrature
  size_t QuadGaussNumPoints(        ); // number   of integration points
  ColD   QuadGaussPosition (size_t k); // position of integration point "k"
  double QuadGaussWeight   (size_t k); // weight   of integration point "k"

  // evaluate the shape function (gradient), at a specific local coordinate:
  // - "xe" : nodal coordinates of the nodes
  // - "xi" : integration point coordinates (in local coordinate-frame)
  // - "w"  : weight factor of the integration point -> sets "V" correctly
  void eval     (const MatD &xe, const ColD &xi, double w=1.);
  void evalGradN(const MatD &xe, const ColD &xi, double w=1.); // do not evaluate "x" and "N"

  // evaluate the shape function (gradient), at one of the default Gauss points for this element
  // (corresponds to the definitions in "QuadGaussPosition" and "QuadGaussWeight")
  void eval     (const MatD &xe, size_t k);
  void evalGradN(const MatD &xe, size_t k); // do not evaluate "x" and "N" (for speed)

  // evaluate the product 'column = ( gradN . tensor2 ) * V'
  // in index notation : column( m*ndim + j ) = dNdx(m,i) * tensor(i,j) * V
  ColD gradN_tensor2(const cppmat::cartesian  ::tensor2 <double> &tensor);
  ColD gradN_tensor2(const cppmat::cartesian  ::tensor2s<double> &tensor);
  ColD gradN_tensor2(const cppmat::cartesian2d::tensor2 <double> &tensor);
  ColD gradN_tensor2(const cppmat::cartesian2d::tensor2s<double> &tensor);
  ColD gradN_tensor2(const cppmat::cartesian3d::tensor2 <double> &tensor);
  ColD gradN_tensor2(const cppmat::cartesian3d::tensor2s<double> &tensor);

  // evaluate the product 'matrix = ( N * scalar * N^T ) * V'
  // in index notation : matrix( m*ndim+i , n*ndim+i ) = N(m) * scalar * N(n) * V
  MatD N_scalar_NT(const double &scalar);

  // evaluate the product 'matrix = ( gradN . tensor4 . gradN^T ) * V'
  // in index notation : matrix( m*ndim+j, n*ndim+k ) = dNdx(m,i) * tensor(i,j,k,l) * dNdx(n,l) * V
  MatD gradN_tensor4_gradNT(const cppmat::cartesian  ::tensor4<double> &tensor);
  MatD gradN_tensor4_gradNT(const cppmat::cartesian2d::tensor4<double> &tensor);
  MatD gradN_tensor4_gradNT(const cppmat::cartesian3d::tensor4<double> &tensor);

  // support functions -> avoid slow a slow iterative scheme
  double det (const MatD &A); // determinant
  MatD   invT(const MatD &A); // transpose of the inverse
};

// ========================================== SOURCE CODE ==========================================

Quad8::Quad8(double thickness)
{
  thick = thickness;

  x    .conservativeResize(nne      );
  J    .conservativeResize(ndim,ndim);
  N    .conservativeResize(nne      );
  dNdxi.conservativeResize(nne ,ndim);
  dNdx .conservativeResize(nne ,ndim);
}

// -------------------------------------------------------------------------------------------------

size_t Quad8::QuadGaussNumPoints()
{
  return 9;
}

// -------------------------------------------------------------------------------------------------

ColD Quad8::QuadGaussPosition(size_t k)
{
  ColD xi(ndim);

  switch ( k ) {
    case 0: xi(0) = -std::sqrt(3./5.); xi(1) = -std::sqrt(3./5.); return xi;
    case 1: xi(0) = +std::sqrt(3./5.); xi(1) = -std::sqrt(3./5.); return xi;
    case 2: xi(0) = +std::sqrt(3./5.); xi(1) = +std::sqrt(3./5.); return xi;
    case 3: xi(0) = -std::sqrt(3./5.); xi(1) = +std::sqrt(3./5.); return xi;
    case 4: xi(0) =  0.              ; xi(1) = -std::sqrt(3./5.); return xi;
    case 5: xi(0) = +std::sqrt(3./5.); xi(1) =  0.              ; return xi;
    case 6: xi(0) =  0.              ; xi(1) = +std::sqrt(3./5.); return xi;
    case 7: xi(0) = -std::sqrt(3./5.); xi(1) =  0.              ; return xi;
    case 8: xi(0) =  0.              ; xi(1) =  0.              ; return xi;
    default: throw std::runtime_error("Unknown input index");
  }
}

// -------------------------------------------------------------------------------------------------

double Quad8::QuadGaussWeight(size_t k)
{
  switch ( k ) {
    case 0: return 25./81.;
    case 1: return 25./81.;
    case 2: return 25./81.;
    case 3: return 25./81.;
    case 4: return 40./81.;
    case 5: return 40./81.;
    case 6: return 40./81.;
    case 7: return 40./81.;
    case 8: return 64./81.;
    default: throw std::runtime_error("Unknown input index");
  }
}

// -------------------------------------------------------------------------------------------------

void Quad8::eval(const MatD &xe, const ColD &xi, double w)
{
  // shape functions at the integration point "xi"
  N(0) = -0.25 * (1.-         xi(0)    ) * (1.-         xi(1)    ) * (1.+xi(0)+xi(1));
  N(1) = -0.25 * (1.+         xi(0)    ) * (1.-         xi(1)    ) * (1.-xi(0)+xi(1));
  N(2) = -0.25 * (1.+         xi(0)    ) * (1.+         xi(1)    ) * (1.-xi(0)-xi(1));
  N(3) = -0.25 * (1.-         xi(0)    ) * (1.+         xi(1)    ) * (1.+xi(0)-xi(1));
  N(4) = +0.5  * (1.-std::pow(xi(0),2.)) * (1.-         xi(1)    )                   ;
  N(5) = +0.5  * (1.+         xi(0)    ) * (1.-std::pow(xi(1),2.))                   ;
  N(6) = +0.5  * (1.-std::pow(xi(0),2.)) * (1.+         xi(1)    )                   ;
  N(7) = +0.5  * (1.-         xi(0)    ) * (1.-std::pow(xi(1),2.))                   ;

  // shape function gradients at the integration point "xi"
  // NB: this is still the gradient w.r.t. the local coordinate "xi" (conversion below)
  dNdxi(0,0) = +0.25*(1.-xi(1))*(2.*xi(0)+xi(1)); dNdxi(0,1) = +0.25*(1-xi(0))*(2.*xi(1)+xi(0));
  dNdxi(1,0) = +0.25*(1.-xi(1))*(2.*xi(0)-xi(1)); dNdxi(1,1) = +0.25*(1+xi(0))*(2.*xi(1)-xi(0));
  dNdxi(2,0) = +0.25*(1.+xi(1))*(2.*xi(0)+xi(1)); dNdxi(2,1) = +0.25*(1+xi(0))*(2.*xi(1)+xi(0));
  dNdxi(3,0) = +0.25*(1.+xi(1))*(2.*xi(0)-xi(1)); dNdxi(3,1) = +0.25*(1-xi(0))*(2.*xi(1)-xi(0));
  dNdxi(4,0) = -1.  *(1.-xi(1))*    xi(0)       ; dNdxi(4,1) = -0.5 *(1+xi(0))*(1.-xi(0))      ;
  dNdxi(5,0) = +0.5 *(1.+xi(1))*(1.-xi(1))      ; dNdxi(5,1) = -1.  *(1+xi(0))*    xi(1)       ;
  dNdxi(6,0) = -1.  *(1.+xi(1))*    xi(0)       ; dNdxi(6,1) = +0.5 *(1+xi(0))*(1.-xi(0))      ;
  dNdxi(7,0) = -0.5 *(1.+xi(1))*(1.-xi(1))      ; dNdxi(7,1) = -1.  *(1-xi(0))*    xi(1)       ;

  // Jacobian at the integration point
  J    = dNdxi.transpose() * xe;
  // integration point volume: 'weight factor' * 'determinant of the Jacobian'
  V    = thick * w * det(J);
  // global coordinates of the integration point
  x    = N.transpose() * xe;
  // shape function gradients at the integration point w.r.t. the global coordinates (i.e. at "x")
  dNdx = dNdxi * invT(J);
}

// -------------------------------------------------------------------------------------------------

void Quad8::evalGradN(const MatD &xe, const ColD &xi, double w)
{
  // shape function gradients at the integration point "xi"
  // NB: this is still the gradient w.r.t. the local coordinate "xi" (conversion below)
  dNdxi(0,0) = +0.25*(1.-xi(1))*(2.*xi(0)+xi(1)); dNdxi(0,1) = +0.25*(1-xi(0))*(2.*xi(1)+xi(0));
  dNdxi(1,0) = +0.25*(1.-xi(1))*(2.*xi(0)-xi(1)); dNdxi(1,1) = +0.25*(1+xi(0))*(2.*xi(1)-xi(0));
  dNdxi(2,0) = +0.25*(1.+xi(1))*(2.*xi(0)+xi(1)); dNdxi(2,1) = +0.25*(1+xi(0))*(2.*xi(1)+xi(0));
  dNdxi(3,0) = +0.25*(1.+xi(1))*(2.*xi(0)-xi(1)); dNdxi(3,1) = +0.25*(1-xi(0))*(2.*xi(1)-xi(0));
  dNdxi(4,0) = -1.  *(1.-xi(1))*    xi(0)       ; dNdxi(4,1) = -0.5 *(1+xi(0))*(1.-xi(0))      ;
  dNdxi(5,0) = +0.5 *(1.+xi(1))*(1.-xi(1))      ; dNdxi(5,1) = -1.  *(1+xi(0))*    xi(1)       ;
  dNdxi(6,0) = -1.  *(1.+xi(1))*    xi(0)       ; dNdxi(6,1) = +0.5 *(1+xi(0))*(1.-xi(0))      ;
  dNdxi(7,0) = -0.5 *(1.+xi(1))*(1.-xi(1))      ; dNdxi(7,1) = -1.  *(1-xi(0))*    xi(1)       ;

  // Jacobian at the integration point
  J    = dNdxi.transpose() * xe;
  // integration point volume: 'weight factor' * 'determinant of the Jacobian'
  V    = thick * w * det(J);
  // shape function gradients at the integration point w.r.t. the global coordinates (i.e. at "x")
  dNdx = dNdxi * invT(J);
}

// -------------------------------------------------------------------------------------------------

void Quad8::eval(const MatD &xe, size_t k)
{
  // integration point coordinates
  // - allocate
	ColD xi(ndim);
  // - set
  switch ( k ) {
    case 0: xi(0) = -std::sqrt(3./5.); xi(1) = -std::sqrt(3./5.); break;
    case 1: xi(0) = +std::sqrt(3./5.); xi(1) = -std::sqrt(3./5.); break;
    case 2: xi(0) = +std::sqrt(3./5.); xi(1) = +std::sqrt(3./5.); break;
    case 3: xi(0) = -std::sqrt(3./5.); xi(1) = +std::sqrt(3./5.); break;
    case 4: xi(0) =  0.              ; xi(1) = -std::sqrt(3./5.); break;
    case 5: xi(0) = +std::sqrt(3./5.); xi(1) =  0.              ; break;
    case 6: xi(0) =  0.              ; xi(1) = +std::sqrt(3./5.); break;
    case 7: xi(0) = -std::sqrt(3./5.); xi(1) =  0.              ; break;
    case 8: xi(0) =  0.              ; xi(1) =  0.              ; break;
    default: throw std::runtime_error("Unknown input index");
  }

  // integration point weight factor
  // - allocate
  double w;
  // - set
  switch ( k ) {
    case 0:  w = 25./81.; break;
    case 1:  w = 25./81.; break;
    case 2:  w = 25./81.; break;
    case 3:  w = 25./81.; break;
    case 4:  w = 40./81.; break;
    case 5:  w = 40./81.; break;
    case 6:  w = 40./81.; break;
    case 7:  w = 40./81.; break;
    case 8:  w = 64./81.; break;
    default: throw std::runtime_error("Unknown input index");
  }

  // shape functions at the integration point "xi"
  N(0) = -0.25 * (1.-         xi(0)    ) * (1.-         xi(1)    ) * (1.+xi(0)+xi(1));
  N(1) = -0.25 * (1.+         xi(0)    ) * (1.-         xi(1)    ) * (1.-xi(0)+xi(1));
  N(2) = -0.25 * (1.+         xi(0)    ) * (1.+         xi(1)    ) * (1.-xi(0)-xi(1));
  N(3) = -0.25 * (1.-         xi(0)    ) * (1.+         xi(1)    ) * (1.+xi(0)-xi(1));
  N(4) = +0.5  * (1.-std::pow(xi(0),2.)) * (1.-         xi(1)    )                   ;
  N(5) = +0.5  * (1.+         xi(0)    ) * (1.-std::pow(xi(1),2.))                   ;
  N(6) = +0.5  * (1.-std::pow(xi(0),2.)) * (1.+         xi(1)    )                   ;
  N(7) = +0.5  * (1.-         xi(0)    ) * (1.-std::pow(xi(1),2.))                   ;

  // shape function gradients at the integration point "xi"
  // NB: this is still the gradient w.r.t. the local coordinate "xi" (conversion below)
  dNdxi(0,0) = +0.25*(1.-xi(1))*(2.*xi(0)+xi(1)); dNdxi(0,1) = +0.25*(1-xi(0))*(2.*xi(1)+xi(0));
  dNdxi(1,0) = +0.25*(1.-xi(1))*(2.*xi(0)-xi(1)); dNdxi(1,1) = +0.25*(1+xi(0))*(2.*xi(1)-xi(0));
  dNdxi(2,0) = +0.25*(1.+xi(1))*(2.*xi(0)+xi(1)); dNdxi(2,1) = +0.25*(1+xi(0))*(2.*xi(1)+xi(0));
  dNdxi(3,0) = +0.25*(1.+xi(1))*(2.*xi(0)-xi(1)); dNdxi(3,1) = +0.25*(1-xi(0))*(2.*xi(1)-xi(0));
  dNdxi(4,0) = -1.  *(1.-xi(1))*    xi(0)       ; dNdxi(4,1) = -0.5 *(1+xi(0))*(1.-xi(0))      ;
  dNdxi(5,0) = +0.5 *(1.+xi(1))*(1.-xi(1))      ; dNdxi(5,1) = -1.  *(1+xi(0))*    xi(1)       ;
  dNdxi(6,0) = -1.  *(1.+xi(1))*    xi(0)       ; dNdxi(6,1) = +0.5 *(1+xi(0))*(1.-xi(0))      ;
  dNdxi(7,0) = -0.5 *(1.+xi(1))*(1.-xi(1))      ; dNdxi(7,1) = -1.  *(1-xi(0))*    xi(1)       ;

  // Jacobian at the integration point
  J    = dNdxi.transpose() * xe;
  // integration point volume: 'weight factor' * 'determinant of the Jacobian'
  V    = thick * w * det(J);
  // global coordinates of the integration point
  x    = N.transpose() * xe;
  // shape function gradients at the integration point w.r.t. the global coordinates (i.e. at "x")
  dNdx = dNdxi * invT(J);
}

// -------------------------------------------------------------------------------------------------

void Quad8::evalGradN(const MatD &xe, size_t k)
{
  // integration point coordinates
  // - allocate
  ColD xi(ndim);
  // - set
  switch ( k ) {
    case 0: xi(0) = -std::sqrt(3./5.); xi(1) = -std::sqrt(3./5.); break;
    case 1: xi(0) = +std::sqrt(3./5.); xi(1) = -std::sqrt(3./5.); break;
    case 2: xi(0) = +std::sqrt(3./5.); xi(1) = +std::sqrt(3./5.); break;
    case 3: xi(0) = -std::sqrt(3./5.); xi(1) = +std::sqrt(3./5.); break;
    case 4: xi(0) =  0.              ; xi(1) = -std::sqrt(3./5.); break;
    case 5: xi(0) = +std::sqrt(3./5.); xi(1) =  0.              ; break;
    case 6: xi(0) =  0.              ; xi(1) = +std::sqrt(3./5.); break;
    case 7: xi(0) = -std::sqrt(3./5.); xi(1) =  0.              ; break;
    case 8: xi(0) =  0.              ; xi(1) =  0.              ; break;
    default: throw std::runtime_error("Unknown input index");
  }

  // integration point weight factor
  // - allocate
  double w;
  // - set
  switch ( k ) {
    case 0:  w = 25./81.; break;
    case 1:  w = 25./81.; break;
    case 2:  w = 25./81.; break;
    case 3:  w = 25./81.; break;
    case 4:  w = 40./81.; break;
    case 5:  w = 40./81.; break;
    case 6:  w = 40./81.; break;
    case 7:  w = 40./81.; break;
    case 8:  w = 64./81.; break;
    default: throw std::runtime_error("Unknown input index");
  }

  // shape function gradients at the integration point "xi"
  // NB: this is still the gradient w.r.t. the local coordinate "xi" (conversion below)
  dNdxi(0,0) = +0.25*(1.-xi(1))*(2.*xi(0)+xi(1)); dNdxi(0,1) = +0.25*(1-xi(0))*(2.*xi(1)+xi(0));
  dNdxi(1,0) = +0.25*(1.-xi(1))*(2.*xi(0)-xi(1)); dNdxi(1,1) = +0.25*(1+xi(0))*(2.*xi(1)-xi(0));
  dNdxi(2,0) = +0.25*(1.+xi(1))*(2.*xi(0)+xi(1)); dNdxi(2,1) = +0.25*(1+xi(0))*(2.*xi(1)+xi(0));
  dNdxi(3,0) = +0.25*(1.+xi(1))*(2.*xi(0)-xi(1)); dNdxi(3,1) = +0.25*(1-xi(0))*(2.*xi(1)-xi(0));
  dNdxi(4,0) = -1.  *(1.-xi(1))*    xi(0)       ; dNdxi(4,1) = -0.5 *(1+xi(0))*(1.-xi(0))      ;
  dNdxi(5,0) = +0.5 *(1.+xi(1))*(1.-xi(1))      ; dNdxi(5,1) = -1.  *(1+xi(0))*    xi(1)       ;
  dNdxi(6,0) = -1.  *(1.+xi(1))*    xi(0)       ; dNdxi(6,1) = +0.5 *(1+xi(0))*(1.-xi(0))      ;
  dNdxi(7,0) = -0.5 *(1.+xi(1))*(1.-xi(1))      ; dNdxi(7,1) = -1.  *(1-xi(0))*    xi(1)       ;

  // Jacobian at the integration point
  J    = dNdxi.transpose() * xe;
  // integration point volume: 'weight factor' * 'determinant of the Jacobian'
  V    = thick * w * det(J);
  // shape function gradients at the integration point w.r.t. the global coordinates (i.e. at "x")
  dNdx = dNdxi * invT(J);
}

// -------------------------------------------------------------------------------------------------

ColD Quad8::gradN_tensor2(const cppmat::cartesian::tensor2<double> &tensor)
{
  ColD col(nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    col( m*2     )  = dNdx(m,0) * tensor(0,0) * V;
    col( m*2     ) += dNdx(m,1) * tensor(1,0) * V;
    col( m*2 + 1 )  = dNdx(m,0) * tensor(0,1) * V;
    col( m*2 + 1 ) += dNdx(m,1) * tensor(1,1) * V;
  }

  return col;
}

// -------------------------------------------------------------------------------------------------

ColD Quad8::gradN_tensor2(const cppmat::cartesian::tensor2s<double> &tensor)
{
  ColD col(nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    col( m*2     )  = dNdx(m,0) * tensor(0,0) * V;
    col( m*2     ) += dNdx(m,1) * tensor(1,0) * V;
    col( m*2 + 1 )  = dNdx(m,0) * tensor(0,1) * V;
    col( m*2 + 1 ) += dNdx(m,1) * tensor(1,1) * V;
  }

  return col;
}

// -------------------------------------------------------------------------------------------------

ColD Quad8::gradN_tensor2(const cppmat::cartesian2d::tensor2<double> &tensor)
{
  ColD col(nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    col( m*2     )  = dNdx(m,0) * tensor(0,0) * V;
    col( m*2     ) += dNdx(m,1) * tensor(1,0) * V;
    col( m*2 + 1 )  = dNdx(m,0) * tensor(0,1) * V;
    col( m*2 + 1 ) += dNdx(m,1) * tensor(1,1) * V;
  }

  return col;
}

// -------------------------------------------------------------------------------------------------

ColD Quad8::gradN_tensor2(const cppmat::cartesian2d::tensor2s<double> &tensor)
{
  ColD col(nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    col( m*2     )  = dNdx(m,0) * tensor(0,0) * V;
    col( m*2     ) += dNdx(m,1) * tensor(1,0) * V;
    col( m*2 + 1 )  = dNdx(m,0) * tensor(0,1) * V;
    col( m*2 + 1 ) += dNdx(m,1) * tensor(1,1) * V;
  }

  return col;
}

// -------------------------------------------------------------------------------------------------

ColD Quad8::gradN_tensor2(const cppmat::cartesian3d::tensor2<double> &tensor)
{
  ColD col(nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    col( m*2     )  = dNdx(m,0) * tensor(0,0) * V;
    col( m*2     ) += dNdx(m,1) * tensor(1,0) * V;
    col( m*2 + 1 )  = dNdx(m,0) * tensor(0,1) * V;
    col( m*2 + 1 ) += dNdx(m,1) * tensor(1,1) * V;
  }

  return col;
}

// -------------------------------------------------------------------------------------------------

ColD Quad8::gradN_tensor2(const cppmat::cartesian3d::tensor2s<double> &tensor)
{
  ColD col(nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    col( m*2     )  = dNdx(m,0) * tensor(0,0) * V;
    col( m*2     ) += dNdx(m,1) * tensor(1,0) * V;
    col( m*2 + 1 )  = dNdx(m,0) * tensor(0,1) * V;
    col( m*2 + 1 ) += dNdx(m,1) * tensor(1,1) * V;
  }

  return col;
}

// -------------------------------------------------------------------------------------------------

MatD Quad8::N_scalar_NT(const double &scalar)
{
  MatD mat(nne*ndim,nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    for ( size_t n = 0 ; n < nne ; ++n ) {
      mat( m*2   , n*2   ) = N(m) * scalar * N(n) * V;
      mat( m*2+1 , n*2+1 ) = N(m) * scalar * N(n) * V;
      mat( m*2   , n*2+1 ) = 0.;
      mat( m*2+1 , n*2   ) = 0.;
   }
 }

  return mat;
}

// -------------------------------------------------------------------------------------------------

MatD Quad8::gradN_tensor4_gradNT(const cppmat::cartesian::tensor4<double> &tensor)
{
  MatD mat(nne*ndim,nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    for ( size_t n = 0 ; n < nne ; ++n ) {
      mat( m*2   , n*2   )  = dNdx(m,0) * tensor(0,0,0,0) * dNdx(n,0) * V;
      mat( m*2   , n*2   ) += dNdx(m,0) * tensor(0,0,0,1) * dNdx(n,1) * V;
      mat( m*2   , n*2+1 )  = dNdx(m,0) * tensor(0,0,1,0) * dNdx(n,0) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,0) * tensor(0,0,1,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2   )  = dNdx(m,0) * tensor(0,1,0,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,0) * tensor(0,1,0,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2+1 )  = dNdx(m,0) * tensor(0,1,1,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,0) * tensor(0,1,1,1) * dNdx(n,1) * V;
      mat( m*2   , n*2   ) += dNdx(m,1) * tensor(1,0,0,0) * dNdx(n,0) * V;
      mat( m*2   , n*2   ) += dNdx(m,1) * tensor(1,0,0,1) * dNdx(n,1) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,1) * tensor(1,0,1,0) * dNdx(n,0) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,1) * tensor(1,0,1,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,1) * tensor(1,1,0,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,1) * tensor(1,1,0,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,1) * tensor(1,1,1,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,1) * tensor(1,1,1,1) * dNdx(n,1) * V;
    }
  }

  return mat;
}

// -------------------------------------------------------------------------------------------------

MatD Quad8::gradN_tensor4_gradNT(const cppmat::cartesian2d::tensor4<double> &tensor)
{
  MatD mat(nne*ndim,nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    for ( size_t n = 0 ; n < nne ; ++n ) {
      mat( m*2   , n*2   )  = dNdx(m,0) * tensor(0,0,0,0) * dNdx(n,0) * V;
      mat( m*2   , n*2   ) += dNdx(m,0) * tensor(0,0,0,1) * dNdx(n,1) * V;
      mat( m*2   , n*2+1 )  = dNdx(m,0) * tensor(0,0,1,0) * dNdx(n,0) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,0) * tensor(0,0,1,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2   )  = dNdx(m,0) * tensor(0,1,0,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,0) * tensor(0,1,0,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2+1 )  = dNdx(m,0) * tensor(0,1,1,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,0) * tensor(0,1,1,1) * dNdx(n,1) * V;
      mat( m*2   , n*2   ) += dNdx(m,1) * tensor(1,0,0,0) * dNdx(n,0) * V;
      mat( m*2   , n*2   ) += dNdx(m,1) * tensor(1,0,0,1) * dNdx(n,1) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,1) * tensor(1,0,1,0) * dNdx(n,0) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,1) * tensor(1,0,1,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,1) * tensor(1,1,0,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,1) * tensor(1,1,0,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,1) * tensor(1,1,1,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,1) * tensor(1,1,1,1) * dNdx(n,1) * V;
    }
  }

  return mat;
}

// -------------------------------------------------------------------------------------------------

MatD Quad8::gradN_tensor4_gradNT(const cppmat::cartesian3d::tensor4<double> &tensor)
{
  MatD mat(nne*ndim,nne*ndim);

  for ( size_t m = 0 ; m < nne ; ++m ) {
    for ( size_t n = 0 ; n < nne ; ++n ) {
      mat( m*2   , n*2   )  = dNdx(m,0) * tensor(0,0,0,0) * dNdx(n,0) * V;
      mat( m*2   , n*2   ) += dNdx(m,0) * tensor(0,0,0,1) * dNdx(n,1) * V;
      mat( m*2   , n*2+1 )  = dNdx(m,0) * tensor(0,0,1,0) * dNdx(n,0) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,0) * tensor(0,0,1,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2   )  = dNdx(m,0) * tensor(0,1,0,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,0) * tensor(0,1,0,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2+1 )  = dNdx(m,0) * tensor(0,1,1,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,0) * tensor(0,1,1,1) * dNdx(n,1) * V;
      mat( m*2   , n*2   ) += dNdx(m,1) * tensor(1,0,0,0) * dNdx(n,0) * V;
      mat( m*2   , n*2   ) += dNdx(m,1) * tensor(1,0,0,1) * dNdx(n,1) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,1) * tensor(1,0,1,0) * dNdx(n,0) * V;
      mat( m*2   , n*2+1 ) += dNdx(m,1) * tensor(1,0,1,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,1) * tensor(1,1,0,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2   ) += dNdx(m,1) * tensor(1,1,0,1) * dNdx(n,1) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,1) * tensor(1,1,1,0) * dNdx(n,0) * V;
      mat( m*2+1 , n*2+1 ) += dNdx(m,1) * tensor(1,1,1,1) * dNdx(n,1) * V;
    }
  }

  return mat;
}

// =================================================================================================

double Quad8::det(const MatD &A)
{
  return A(0,0) * A(1,1) - A(0,1) * A(1,0);
}

// =================================================================================================

MatD Quad8::invT(const MatD &A)
{
  MatD   C(2,2);
  double D = det(A);

  C(0,0) =       A(1,1) / D;
  C(1,0) = -1. * A(0,1) / D;
  C(0,1) = -1. * A(1,0) / D;
  C(1,1) =       A(0,0) / D;

  return C;
}

// =================================================================================================

} // namespace GooseFEM

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENTQUAD4_H
#define GOOSEFEM_ELEMENTQUAD4_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// =================================== GooseFEM::Element::Quad4 ====================================

namespace GooseFEM {
namespace Element {
namespace Quad4 {

// ================================ GooseFEM::Element::Quad4::Gauss ================================

namespace Gauss {
inline size_t nip(); // number of integration points
inline ArrD   xi();  // integration point coordinates (local coordinates)
inline ArrD   w();   // integration point weights
}

// ================================ GooseFEM::Element::Quad4::Nodal ================================

namespace Nodal {
inline size_t nip(); // number of integration points
inline ArrD   xi();  // integration point coordinates (local coordinates)
inline ArrD   w();   // integration point weights
}

// =================================================================================================

// ------------------------------------------ quadrature -------------------------------------------

class Quadrature
{
private:

  // dimensions (flexible)
  size_t m_nelem; // number of elements
  size_t m_nip;   // number of integration points

  // dimensions (fixed for this element type)
  static const size_t m_nne=4;  // number of nodes per element
  static const size_t m_ndim=2; // number of dimensions

  // data arrays
  ArrD m_x;    // nodal positions stored per element [nelem, nne, ndim]
  ArrD m_w;    // weight of each integration point [nip]
  ArrD m_xi;   // local coordinate of each integration point [nip, ndim]
  ArrD m_N;    // shape functions [nip, nne]
  ArrD m_dNxi; // shape function gradients w.r.t. local  coordinate [nip, nne, ndim]
  ArrD m_dNx;  // shape function gradients w.r.t. global coordinate [nelem, nip, nne, ndim]
  ArrD m_vol;  // integration point volume [nelem, nip]

private:

  // compute "vol" and "dNdx" based on current "x"
  void compute_dN();

public:

  // convention:
  //    "elemmat"  -  matrices stored per element       -  ArrD  -  [nelem, nne*ndim, nne*ndim]
  //    "elemvec"  -  nodal vectors stored per element  -  ArrD  -  [nelem, nne, ndim]
  //    "qtensor"  -  integration point tensor          -  ArrD  -  [nelem, nip, #tensor-components]
  //    "qscalar"  -  integration point scalar          -  ArrD  -  [nelem, nip]
  //
  // alias:
  //    T2   = cppmat::cartesian2d::tensor2<double>         -  #tensor-components = 4
  //    T2s  = cppmat::cartesian2d::tensor2s<double>        -  #tensor-components = 3
  //    vT2  = cppmat::view::cartesian2d::tensor2<double>   -  #tensor-components = 4
  //    vT2s = cppmat::view::cartesian2d::tensor2s<double>  -  #tensor-components = 3

  // constructor: integration point coordinates and weights are optional (default: Gauss)
  Quadrature(){};
  Quadrature(const ArrD &x, const ArrD &xi=ArrD(), const ArrD &w=ArrD());

  // update the nodal positions (shape of "x" should match the earlier definition)
  void update_x(const ArrD &x);

  // return dimensions
  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t ndim()  const; // number of dimension
  size_t nip()   const; // number of integration points

  // return integration volume
  ArrD dV() const;              // returns: qscalar
  ArrD dV(size_t ncomp) const;  // returns: qtensor (same volume per tensor-component)

  // dyadic product "qtensor(i,j) += dNdx(m,i) * elemvec(m,j)", its transpose and its symmetric part
  // - allow template (e.g. T2/T2s, or higher dimensional tensors)
  template<class T> ArrD gradN_vector   (const ArrD &elemvec) const; // returns: qtensor
  template<class T> ArrD gradN_vector_T (const ArrD &elemvec) const; // returns: qtensor
  template<class T> ArrD symGradN_vector(const ArrD &elemvec) const; // returns: qtensor
  // - default template
  ArrD gradN_vector   (const ArrD &elemvec) const; // template: T2
  ArrD gradN_vector_T (const ArrD &elemvec) const; // template: T2
  ArrD symGradN_vector(const ArrD &elemvec) const; // template: T2s

  // integral of the scalar product "elemmat(m*ndim+i,n*ndim+i) += N(m) * qscalar * N(n) * dV"
  ArrD int_N_scalar_NT_dV(const ArrD &qscalar) const; // returns: elemmat

  // integral of the dot product "elemvec(m,j) += dNdx(m,i) * qtensor(i,j) * dV"
  // - allow template (e.g. vT2/vT2s, or higher dimensional tensors)
  template<class T> ArrD int_gradN_dot_tensor2_dV(const ArrD &qtensor) const; // returns: elemvec
  // - default template
  ArrD int_gradN_dot_tensor2_dV (const ArrD &qtensor) const; // template: vT2/vT2s (auto-select)
  ArrD int_gradN_dot_tensor2s_dV(const ArrD &qtensor) const; // template: vT2s

};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

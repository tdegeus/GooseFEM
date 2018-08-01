/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENTHEX8_H
#define GOOSEFEM_ELEMENTHEX8_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ==================================== GooseFEM::Element::Hex8 ====================================

namespace GooseFEM {
namespace Element {
namespace Hex8 {

// ======================================== tensor algebra =========================================

static const size_t ndim = 3;

using T2 = cppmat::tiny::cartesian::tensor2<double,ndim>;

inline double inv(const T2 &A, T2 &Ainv);

// ================================ GooseFEM::Element::Hex8::Gauss =================================

namespace Gauss {
inline size_t nip(); // number of integration points
inline ArrD   xi();  // integration point coordinates (local coordinates)
inline ArrD   w();   // integration point weights
}

// ================================ GooseFEM::Element::Hex8::Nodal =================================

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
  static const size_t m_nne=8;  // number of nodes per element
  static const size_t m_ndim=3; // number of dimensions

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
  //    T2   = cppmat::tiny::cartesian::tensor2<double,3>         -  #tensor-components = 9
  //    T2s  = cppmat::tiny::cartesian::tensor2s<double,3>        -  #tensor-components = 6

  // constructor: integration point coordinates and weights are optional (default: Gauss)
  Quadrature() = default;
  Quadrature(const ArrD &x);
  Quadrature(const ArrD &x, const ArrD &xi, const ArrD &w);

  // update the nodal positions (shape of "x" should match the earlier definition)
  void update_x(const ArrD &x);

  // return dimensions
  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t ndim()  const; // number of dimension
  size_t nip()   const; // number of integration points

  // return integration volume
  ArrD dV(size_t ncomp=0) const;  // returns: qscalar/qtensor (same volume per tensor-component)

  // dyadic product "qtensor(i,j) += dNdx(m,i) * elemvec(m,j)", its transpose and its symmetric part
  // - allow template (e.g. T2/T2s, or higher dimensional tensors)
  template<class T> void gradN_vector   (const ArrD &elemvec, ArrD &qtensor) const;
  template<class T> void gradN_vector_T (const ArrD &elemvec, ArrD &qtensor) const;
  template<class T> void symGradN_vector(const ArrD &elemvec, ArrD &qtensor) const;
  // -
  template<class T> ArrD gradN_vector   (const ArrD &elemvec) const; // returns: qtensor
  template<class T> ArrD gradN_vector_T (const ArrD &elemvec) const; // returns: qtensor
  template<class T> ArrD symGradN_vector(const ArrD &elemvec) const; // returns: qtensor
  // - default template
  void gradN_vector   (const ArrD &elemvec, ArrD &qtensor) const; // template: T2
  void gradN_vector_T (const ArrD &elemvec, ArrD &qtensor) const; // template: T2
  void symGradN_vector(const ArrD &elemvec, ArrD &qtensor) const; // template: T2s
  // -
  ArrD gradN_vector   (const ArrD &elemvec) const; // template: T2
  ArrD gradN_vector_T (const ArrD &elemvec) const; // template: T2
  ArrD symGradN_vector(const ArrD &elemvec) const; // template: T2s

  // integral of the scalar product "elemmat(m*ndim+i,n*ndim+i) += N(m) * qscalar * N(n) * dV"
  void int_N_scalar_NT_dV(const ArrD &qscalar, ArrD &elemmat) const; // returns: elemmat
  // -
  ArrD int_N_scalar_NT_dV(const ArrD &qscalar) const; // returns: elemmat

  // integral of the dot product "elemvec(m,j) += dNdx(m,i) * qtensor(i,j) * dV"
  // - allow template (e.g. T2/T2s, or higher dimensional tensors)
  template<class T> void int_gradN_dot_tensor2_dV(const ArrD &qtensor, ArrD &elemvec) const; // returns: elemvec
  //
  template<class T> ArrD int_gradN_dot_tensor2_dV(const ArrD &qtensor) const; // returns: elemvec
  // - default template
  void int_gradN_dot_tensor2_dV (const ArrD &qtensor, ArrD &elemvec) const; // template: T2/T2s (auto-select)
  void int_gradN_dot_tensor2s_dV(const ArrD &qtensor, ArrD &elemvec) const; // template: T2s
  // -
  ArrD int_gradN_dot_tensor2_dV (const ArrD &qtensor) const; // template: T2/T2s (auto-select)
  ArrD int_gradN_dot_tensor2s_dV(const ArrD &qtensor) const; // template: T2s

};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

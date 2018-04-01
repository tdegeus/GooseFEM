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
size_t nip();         // number of integration points
ArrD   coordinates(); // integration point coordinates (local coordinates)
ArrD   weights();     // integration point weights
}

// ================================ GooseFEM::Element::Quad4::Nodal ================================

namespace Nodal {
size_t nip();         // number of integration points
ArrD   coordinates(); // integration point coordinates (local coordinates)
ArrD   weights();     // integration point weights
}

// =================================================================================================

// ------------------------------------------ quadrature -------------------------------------------

class Quadrature
{
private:

  // dimensions
  size_t m_nelem;               // number of elements
  size_t m_nip;                 // number of integration positions
  static const size_t m_nne=4;  // number of nodes per element
  static const size_t m_ndim=2; // number of dimensions
  // data arrays
  ArrD   m_x;    // element vector with nodal positions [nelem, nne, ndim]
  ArrD   m_w;    // weight of each integration point [nip]
  ArrD   m_xi;   // local coordinate of each integration point [nip, ndim]
  ArrD   m_N;    // shape functions w.r.t. local coordinate [nip, nne]
  ArrD   m_dNxi; // shape function gradients w.r.t. local coordinate [nip, nne, ndim]
  ArrD   m_dNx;  // shape function gradients w.r.t. global coordinate [nelem, nip, nne, ndim]
  ArrD   m_vol;  // integration point volume [nelem, nip]

private:

  // compute "vol" and "dNdx" based on current "x"
  void compute_dN();

public:

  // constructor, integration point coordinates and weights are optional (default: Gauss)
  Quadrature(const ArrD &x, const ArrD &xi=ArrD(), const ArrD &w=ArrD());

  // update the nodal positions (shape of "x" should match the earlier definition)
  void update_x(const ArrD &x);

  // return dimensions
  size_t nelem();
  size_t nne();
  size_t ndim();
  size_t nip();

  // dyadic product "qtensor(i,j) += dNdx(m,i) * elemvec(m,j)", its transpose and its symmetric part
  //
  // input : element vector            -  [nelem, nne, ndim]
  // output: integration point tensor  -  [nelem, nip, #tensor-components]
  //
  // - allow template (e.g. 'cppmat::cartesian2d::tensor2<double>')
  template<class T> ArrD gradN_vector   (const ArrD &elemvec);
  template<class T> ArrD gradN_vector_T (const ArrD &elemvec);
  template<class T> ArrD symGradN_vector(const ArrD &elemvec);
  // - default template with cppmat::cartesian2d::...<double>
  ArrD gradN_vector   (const ArrD &elemvec); // tensor2  : #tensor-components = ndim*ndim
  ArrD gradN_vector_T (const ArrD &elemvec); // tensor2  : #tensor-components = ndim*ndim
  ArrD symGradN_vector(const ArrD &elemvec); // tensor2s : #tensor-components = (ndim+1)*ndim/2

  // integral of the scalar product "elemmat(m*ndim+i,n*ndim+i) += N(m) * qscalar * N(n) * dV"
  //
  // input : integration point scalar  -  [nelem, nip]
  // output: element matrix            -  [nelem, nne*ndim, nne*ndim]
  //
  ArrD int_N_scalar_NT_dV(const ArrD &qscalar);

  // integral of the dot product "elemvec(m,j) += dNdx(m,i) * qtensor(i,j) * dV"
  //
  // input : integration point tensor  -  [nelem, nip, #tensor-components]
  // output: element vector            -  [nelem, nne, ndim]
  //
  // - allow template (e.g. 'cppmat::cartesian2d::tensor2<double>')
  template<class T> ArrD int_gradN_dot_tensor2_dV(const ArrD &qtensor);
  // - default template with cppmat::cartesian2d::...<double>
  ArrD int_gradN_dot_tensor2_dV (const ArrD &qtensor); // tensor2 / tensor2s (automatic selection)
  ArrD int_gradN_dot_tensor2s_dV(const ArrD &qtensor); // tensor2s

  // integral of a tensor "tensor(i,j) += qtensor(i,j) * dV" (a.k.a. volume average)
  //
  // input : integration point tensor  -  [nelem, nip, #tensor-components]
  // output: (element) tensor          -  [#tensor-components]
  //
  // - allow template (e.g. 'cppmat::cartesian2d::tensor2<double>')
  template<class T> T int_tensor2_dV(const ArrD &qtensor);
  template<class T> T int_tensor2_dV(const ArrD &qtensor, size_t e);
  // - default template with cppmat::cartesian2d::...<double>
  cppmat::cartesian2d::tensor2 <double> int_tensor2_dV (const ArrD &qtensor);           // tensor2
  cppmat::cartesian2d::tensor2s<double> int_tensor2s_dV(const ArrD &qtensor);           // tensor2s
  cppmat::cartesian2d::tensor2 <double> int_tensor2_dV (const ArrD &qtensor, size_t e); // tensor2
  cppmat::cartesian2d::tensor2s<double> int_tensor2s_dV(const ArrD &qtensor, size_t e); // tensor2s

};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

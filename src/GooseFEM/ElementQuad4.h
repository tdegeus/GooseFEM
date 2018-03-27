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

// --------------------------------------- Gauss integration ---------------------------------------

class Gauss
{
private:

  // dimensions
  size_t m_nelem;               // number of elements
  static const size_t m_nne=4;  // number of nodes per element
  static const size_t m_ndim=2; // number of dimensions
  static const size_t m_nip=4;  // number of integration positions

  // data arrays
  ArrD   m_x;    // element vector with nodal positions [nelem, nne, ndim]
  ArrD   m_w;    // weight of each integration point [nip]
  ArrD   m_xi;   // local coordinate of each integration point [nip, ndim]
  ArrD   m_dNxi; // shape function gradients wrt local coordinate [nip, nne, ndim]
  ArrD   m_dNx;  // shape function gradients wrt global coordinate [nelem, nip, nne, ndim]
  ArrD   m_vol;  // integration point volume [nelem, nip]

private:

  // compute "vol" and "dNdx" based on current "x"
  void compute_dN();

public:

  // constructor
  Gauss(const ArrD &x);

  // update the element vectors with nodal positions
  // (the shape of "x" should match the earlier definition)
  void update_x(const ArrD &x);

  // return dimensions
  size_t nelem();
  size_t nne();
  size_t ndim();
  size_t nip();

  // dyadic product "out(i,j) = dNdx(m,i) * inp(m,j)"
  // - allow template (allows modification of tensor storage)
  template<class T>
  ArrD gradN_vector(const ArrD &vector);
  // - default template: "cppmat::cartesian2d::tensor2<double>"
  ArrD gradN_vector(const ArrD &vector);

  // dyadic product "out(j,i) = dNdx(m,i) * inp(m,j)"
  // - allow template (allows modification of tensor storage)
  template<class T>
  ArrD vector_GradN(const ArrD &vector);
  // - default template: "cppmat::cartesian2d::tensor2<double>"
  ArrD vector_GradN(const ArrD &vector);

  // symmetric dyadic product "out(i,j) = dNdx(m,i) * inp(m,j)"
  // (the output is symmetrized)
  // - allow template (allows modification of tensor storage)
  template<class T>
  ArrD symGradN_vector(const ArrD &vector);
  // - default template: "cppmat::cartesian2d::tensor2s<double>"
  ArrD symGradN_vector(const ArrD &vector);

  // dot product "out(m,j) = dNdx(m,i) * inp(i,j)"
  // - allow template (allows modification of tensor storage)
  template<class T>
  ArrD gradN_dot_tensor2(const ArrD &tensor);
  // - infer template from input, may be:
  //   "cppmat::cartesian2d::tensor2<double>" or "cppmat::cartesian2d::tensor2s<double>"
  ArrD gradN_dot_tensor2 (const ArrD &tensor);
  // - default template: "cppmat::cartesian2d::tensor2s<double>"
  ArrD gradN_dot_tensor2s(const ArrD &tensor);

  // dot product "out(m,j) = dNdx(m,i) * inp(i,j) * dV"
  // - allow template (allows modification of tensor storage)
  template<class T>
  ArrD gradN_dot_tensor2_dV(const ArrD &tensor);
  // - infer template from input, may be:
  //   "cppmat::cartesian2d::tensor2<double>" or "cppmat::cartesian2d::tensor2s<double>"
  ArrD gradN_dot_tensor2_dV(const ArrD &tensor);

  // compute volume averaged tensor of an element
  // (has to be templated with the proper cppmat-tensor)
  template<class T>
  T average_tensor2(const ArrD &inp, size_t e);

  cppmat::cartesian2d::tensor2 <double> average_tensor2 (const ArrD &inp, size_t e);
  cppmat::cartesian2d::tensor2s<double> average_tensor2s(const ArrD &inp, size_t e);

  // compute volume averaged tensor of all elements
  // (has to be templated with the proper cppmat-tensor)
  template<class T>
  T average_tensor2(const ArrD &inp);

  cppmat::cartesian2d::tensor2 <double> average_tensor2 (const ArrD &inp);
  cppmat::cartesian2d::tensor2s<double> average_tensor2s(const ArrD &inp);
};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

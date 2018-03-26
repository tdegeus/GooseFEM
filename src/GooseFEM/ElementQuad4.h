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

  size_t m_nelem;         // number of elements
  static const size_t m_nne=4;   // number of nodes per element
  static const size_t m_ndim=2;  // number of dimensions
  static const size_t m_nip=4;   // number of integration positions
  ArrD   m_x;             // element vector with nodal positions [nelem, nne, ndim]
  ArrD   m_w;             // weight of each integration point [nip]
  ArrD   m_xi;            // local coordinate of each integration point [nip, ndim]
  ArrD   m_dNxi;          // shape function gradients wrt local coordinate [nip, nne, ndim]
  ArrD   m_dNx;           // shape function gradients wrt global coordinate [nelem, nip, nne, ndim]
  ArrD   m_vol;           // integration point volume [nelem, nip]

private:

  // compute "vol" and "dNdx" based on current "x"
  void compute_dN();

public:

  // constructor
  Gauss(const ArrD &x);

  void update_x(const ArrD &x);

  // number of integration points
  size_t nip();

  // dyadic product of the shape function gradient and an element vector [nelem, nne, ndim]
  ArrD gradN_vector(const ArrD &vector);

  ArrD vector_GradN(const ArrD &vector);

  ArrD symGradN_vector(const ArrD &vector);

  ArrD gradN_dot_tensor(const ArrD &tensor);

  template<class T>
  ArrD gradN_dot_tensor(const ArrD &tensor);

  // ArrD dN_vector(const ArrD &vector);
  // void dN_vector(const ArrD &vector, const ArrD &tensor);

  // ArrD sym_dN_vector(const ArrD &vector);
  // void sym_dN_vector(const ArrD &vector, const ArrD &tensor);

  // ArrD dN_tensor(const ArrD &tensor);
  // void dN_tensor(const ArrD &tensor, const ArrD &vector);

  // ArrD vol_dN_tensor(const ArrD &tensor);
  // void vol_dN_tensor(const ArrD &tensor, const ArrD &vector);


};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

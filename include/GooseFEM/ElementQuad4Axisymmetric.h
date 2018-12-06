/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENTQUAD4AXISYMMETRIC_H
#define GOOSEFEM_ELEMENTQUAD4AXISYMMETRIC_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// =================================================================================================

namespace GooseFEM {
namespace Element {
namespace Quad4 {

// -------------------------------------------------------------------------------------------------

class QuadratureAxisymmetric
{
public:

  // fixed dimensions:
  //    ndim = 2   -  number of dimensions
  //    nne  = 4   -  number of nodes per element
  //    tdim = 3   -  number of dimensions of tensors
  //
  // naming convention:
  //    "elemmat"  -  matrices stored per element       -  [nelem, nne*ndim, nne*ndim]
  //    "elemvec"  -  nodal vectors stored per element  -  [nelem, nne, ndim]
  //    "qtensor"  -  integration point tensor          -  [nelem, nip, tdim, tdim]
  //    "qscalar"  -  integration point scalar          -  [nelem, nip]

  // constructor: integration point coordinates and weights are optional (default: Gauss)

  QuadratureAxisymmetric() = default;

  QuadratureAxisymmetric(const xt::xtensor<double,3> &x);

  QuadratureAxisymmetric(const xt::xtensor<double,3> &x, const xt::xtensor<double,2> &xi, const xt::xtensor<double,1> &w);

  // update the nodal positions (shape of "x" should match the earlier definition)

  void update_x(const xt::xtensor<double,3> &x);

  // return dimensions

  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t ndim()  const; // number of dimension
  size_t nip()   const; // number of integration points

  // return integration volume

  void dV(xt::xtensor<double,2> &qscalar) const;

  void dV(xt::xtensor<double,4> &qtensor) const; // same volume for all tensor components

  // dyadic product "qtensor(i,j) += B(m,i,j,k) * elemvec(m,k)", its transpose and its symmetric part

  void gradN_vector(const xt::xtensor<double,3> &elemvec,
    xt::xtensor<double,4> &qtensor) const;

  void gradN_vector_T(const xt::xtensor<double,3> &elemvec,
    xt::xtensor<double,4> &qtensor) const;

  void symGradN_vector(const xt::xtensor<double,3> &elemvec,
    xt::xtensor<double,4> &qtensor) const;

  // integral of the scalar product "elemmat(m*ndim+i,n*ndim+i) += N(m) * qscalar * N(n) * dV"

  void int_N_scalar_NT_dV(const xt::xtensor<double,2> &qscalar,
    xt::xtensor<double,3> &elemmat) const;

  // integral of the assembled product "fm = ( Bm^T : qtensor ) dV"

  void int_gradN_dot_tensor2_dV(const xt::xtensor<double,4> &qtensor,
    xt::xtensor<double,3> &elemvec) const;

  // integral of the assembled product "Kmn = ( Bm^T : qtensor : Bn ) dV

  void int_gradN_dot_tensor4_dot_gradNT_dV(const xt::xtensor<double,6> &qtensor,
    xt::xtensor<double,3> &elemmat) const;

  // auto allocation of the functions above

  xt::xtensor<double,2> dV() const;

  xt::xtensor<double,4> dVtensor() const;

  xt::xtensor<double,4> gradN_vector(const xt::xtensor<double,3> &elemvec) const;

  xt::xtensor<double,4> gradN_vector_T(const xt::xtensor<double,3> &elemvec) const;

  xt::xtensor<double,4> symGradN_vector(const xt::xtensor<double,3> &elemvec) const;

  xt::xtensor<double,3> int_N_scalar_NT_dV(const xt::xtensor<double,2> &qscalar) const;

  xt::xtensor<double,3> int_gradN_dot_tensor2_dV(const xt::xtensor<double,4> &qtensor) const;

  xt::xtensor<double,3> int_gradN_dot_tensor4_dot_gradNT_dV(const xt::xtensor<double,6> &qtensor) const;

private:

  // compute "vol" and "B" based on current "x"
  void compute_dN();

private:

  // dimensions (flexible)
  size_t m_nelem; // number of elements
  size_t m_nip;   // number of integration points

  // dimensions (fixed for this element type)
  static const size_t m_nne=4;  // number of nodes per element
  static const size_t m_ndim=2; // number of dimensions
  static const size_t m_tdim=3; // number of dimensions of tensors

  // data arrays
  xt::xtensor<double,3> m_x;    // nodal positions stored per element [nelem, nne, ndim]
  xt::xtensor<double,1> m_w;    // weight of each integration point [nip]
  xt::xtensor<double,2> m_xi;   // local coordinate of each integration point [nip, ndim]
  xt::xtensor<double,2> m_N;    // shape functions [nip, nne]
  xt::xtensor<double,3> m_dNxi; // shape function gradients w.r.t. local  coordinate [nip, nne, ndim]
  xt::xtensor<double,6> m_B;    // B-matrix [nelem, nne, tdim, tdim, tdim]
  xt::xtensor<double,2> m_vol;  // integration point volume [nelem, nip]

};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#include "ElementQuad4Axisymmetric.hpp"

// =================================================================================================

#endif

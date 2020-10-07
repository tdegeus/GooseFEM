/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_ELEMENTQUAD4AXISYMMETRIC_H
#define GOOSEFEM_ELEMENTQUAD4AXISYMMETRIC_H

#include "config.h"

namespace GooseFEM {
namespace Element {
namespace Quad4 {

class QuadratureAxisymmetric {
public:
    // Fixed dimensions:
    //    ndim = 2   -  number of dimensions
    //    nne  = 4   -  number of nodes per element
    //    tdim = 3   -  number of dimensions of tensors
    //
    // Naming convention:
    //    "elemmat"  -  matrices stored per element       -  [nelem, nne*ndim, nne*ndim]
    //    "elemvec"  -  nodal vectors stored per element  -  [nelem, nne, ndim]
    //    "qtensor"  -  integration point tensor          -  [nelem, nip, tdim, tdim]
    //    "qscalar"  -  integration point scalar          -  [nelem, nip]

    // Constructor: integration point coordinates and weights are optional (default: Gauss)
    QuadratureAxisymmetric() = default;

    QuadratureAxisymmetric(const xt::xtensor<double, 3>& x);

    QuadratureAxisymmetric(
        const xt::xtensor<double, 3>& x,
        const xt::xtensor<double, 2>& xi,
        const xt::xtensor<double, 1>& w);

    // Update the nodal positions (shape of "x" should match the earlier definition)
    void update_x(const xt::xtensor<double, 3>& x);

    // Return dimensions
    size_t nelem() const; // number of elements
    size_t nne() const;   // number of nodes per element
    size_t ndim() const;  // number of dimension
    size_t nip() const;   // number of integration points

    xt::xtensor<double, 4> GradN() const;

    // Convert "qscalar" to "qtensor" of certain rank
    template <size_t rank = 0>
    void asTensor(const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 2 + rank>& qtensor) const;

    // Return integration volume
    xt::xtensor<double, 2> dV() const;

    // Dyadic product (and its transpose and symmetric part)
    // qtensor(i,j) += B(m,i,j,k) * elemvec(m,k)
    void gradN_vector(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;
    void gradN_vector_T(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;
    void symGradN_vector(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;

    // Integral of the scalar product
    // elemmat(m*ndim+i,n*ndim+i) += N(m) * qscalar * N(n) * dV
    void int_N_scalar_NT_dV(
        const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 3>& elemmat) const;

    // Integral of the assembled product
    // fm = ( Bm^T : qtensor ) dV
    void int_gradN_dot_tensor2_dV(
        const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const;

    // Integral of the assembled product
    // Kmn = ( Bm^T : qtensor : Bn ) dV
    void int_gradN_dot_tensor4_dot_gradNT_dV(
        const xt::xtensor<double, 6>& qtensor, xt::xtensor<double, 3>& elemmat) const;

    // Auto-allocation of the functions above
    xt::xtensor<double, 4> GradN_vector(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 4> GradN_vector_T(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 4> SymGradN_vector(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 3> Int_N_scalar_NT_dV(const xt::xtensor<double, 2>& qscalar) const;
    xt::xtensor<double, 3> Int_gradN_dot_tensor2_dV(const xt::xtensor<double, 4>& qtensor) const;
    xt::xtensor<double, 3> Int_gradN_dot_tensor4_dot_gradNT_dV(const xt::xtensor<double, 6>& qtensor) const;

    // Convert "qscalar" to "qtensor" of certain rank
    template <size_t rank = 0>
    xt::xtensor<double, 2 + rank> AsTensor(const xt::xtensor<double, 2>& qscalar) const;

    xt::xarray<double> AsTensor(size_t rank, const xt::xtensor<double, 2>& qscalar) const;

    template <size_t rank = 0>
    xt::xtensor<double, rank + 2> AllocateQtensor() const;

    template <size_t rank = 0>
    xt::xtensor<double, rank + 2> AllocateQtensor(double val) const;

    xt::xarray<double> AllocateQtensor(size_t rank) const;
    xt::xarray<double> AllocateQtensor(size_t rank, double val) const;

    xt::xtensor<double, 2> AllocateQscalar() const;
    xt::xtensor<double, 2> AllocateQscalar(double val) const;

private:
    // Compute "vol" and "B" based on current "x"
    void compute_dN();

private:
    // Dimensions (flexible)
    size_t m_nelem; // number of elements
    size_t m_nip;   // number of integration points

    // Dimensions (fixed for this element type)
    static const size_t m_nne = 4;  // number of nodes per element
    static const size_t m_ndim = 2; // number of dimensions
    static const size_t m_tdim = 3; // number of dimensions of tensors

    // Data arrays
    xt::xtensor<double, 3> m_x;    // nodal positions stored per element [nelem, nne, ndim]
    xt::xtensor<double, 1> m_w;    // weight of each integration point [nip]
    xt::xtensor<double, 2> m_xi;   // local coordinate of each integration point [nip, ndim]
    xt::xtensor<double, 2> m_N;    // shape functions [nip, nne]
    xt::xtensor<double, 3> m_dNxi; // shape function grad. wrt local  coor. [nip, nne, ndim]
    xt::xtensor<double, 6> m_B;    // B-matrix [nelem, nne, tdim, tdim, tdim]
    xt::xtensor<double, 2> m_vol;  // integration point volume [nelem, nip]
};

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#include "ElementQuad4Axisymmetric.hpp"

#endif

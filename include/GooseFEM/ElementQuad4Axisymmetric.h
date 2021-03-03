/**
Quadrature for 4-noded quadrilateral element in 2d (GooseFEM::Mesh::ElementType::Quad4),
in an axisymmetric coordinated system.

\file ElementQuad4Axisymmetric.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENTQUAD4AXISYMMETRIC_H
#define GOOSEFEM_ELEMENTQUAD4AXISYMMETRIC_H

#include "config.h"

namespace GooseFEM {
namespace Element {
namespace Quad4 {

class QuadratureAxisymmetric : public GooseFEM::Element::QuadratureBase<4, 2, 3> {
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

    xt::xtensor<double, 4> GradN() const;

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

private:
    // Compute "vol" and "B" based on current "x"
    void compute_dN();

private:
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

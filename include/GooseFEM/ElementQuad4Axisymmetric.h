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

/**
Interpolation and quadrature.

Fixed dimensions:
-   ``ndim = 2``: number of dimensions.
-   ``tdim = 3``: number of dimensions or tensor.
-   ``nne = 4``: number of nodes per element.

Naming convention:
-    ``elemmat``:  matrices stored per element, [#nelem, #nne * #ndim, #nne * #ndim]
-    ``elemvec``:  nodal vectors stored per element, [#nelem, #nne, #ndim]
-    ``qtensor``:  integration point tensor, [#nelem, #nip, #tdim, #tdim]
-    ``qscalar``:  integration point scalar, [#nelem, #nip]
*/
class QuadratureAxisymmetric : public GooseFEM::Element::QuadratureBaseCartesian<4, 2, 3> {
public:

    QuadratureAxisymmetric() = default;

    QuadratureAxisymmetric(const xt::xtensor<double, 3>& x);

    QuadratureAxisymmetric(
        const xt::xtensor<double, 3>& x,
        const xt::xtensor<double, 2>& xi,
        const xt::xtensor<double, 1>& w);

    /**
    Get the B-matrix (shape function gradients) (in global coordinates).
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \return ``B`` matrix stored per element, per integration point [#nelem, #nne, #tdim, #tdim, #tdim]
    */
    xt::xtensor<double, 6> B() const;

    // qtensor(e, q, i, j) += B(e, q, m, i, j, k) * elemvec(e, q, m, k)
    void gradN_vector(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const override;

    void gradN_vector_T(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const override;

    void symGradN_vector(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const override;

    // elemmat(e, q, m * ndim + i, n * ndim + i) +=
    //     N(e, q, m) * qscalar(e, q) * N(e, q, n) * dV(e, q)
    void int_N_scalar_NT_dV(
        const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 3>& elemmat) const override;

    // fm = ( Bm^T : qtensor ) dV
    void int_gradN_dot_tensor2_dV(
        const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const override;

    // Kmn = ( Bm^T : qtensor : Bn ) dV
    void int_gradN_dot_tensor4_dot_gradNT_dV(
        const xt::xtensor<double, 6>& qtensor, xt::xtensor<double, 3>& elemmat) const override;

protected:
    void compute_dN() override;

private:
    xt::xtensor<double, 4> GradN() const override;

private:
    xt::xtensor<double, 6> m_B; ///< B-matrix [#nelem, #nne, #tdim, #tdim, #tdim]
};

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#include "ElementQuad4Axisymmetric.hpp"

#endif

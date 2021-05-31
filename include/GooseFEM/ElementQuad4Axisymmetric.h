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
class QuadratureAxisymmetric : public QuadratureBaseCartesian<QuadratureAxisymmetric> {
public:

    QuadratureAxisymmetric() = default;

    /**
    Constructor: use default Gauss integration.
    The following is pre-computed during construction:
    -   the shape functions,
    -   the shape function gradients (in local and global) coordinates,
    -   the integration points volumes.
    They can be reused without any cost.
    They only have to be recomputed when the nodal position changes
    (note that they are assumed to be constant under a small-strain assumption).
    In that case use update_x() to update the nodal positions and
    to recompute the above listed quantities.

    \param x nodal coordinates (``elemvec``).
    */
    template <class T>
    QuadratureAxisymmetric(const T& x);

    /**
    Constructor with custom integration.
    The following is pre-computed during construction:
    -   the shape functions,
    -   the shape function gradients (in local and global) coordinates,
    -   the integration points volumes.
    They can be reused without any cost.
    They only have to be recomputed when the nodal position changes
    (note that they are assumed to be constant under a small-strain assumption).
    In that case use update_x() to update the nodal positions and
    to recompute the above listed quantities.

    \param x nodal coordinates (``elemvec``).
    \param xi Integration point coordinates (local coordinates) [#nip].
    \param w Integration point weights [#nip].
    */
    template <class T, class X, class W>
    QuadratureAxisymmetric(const T& x, const X& xi, const W& w);

    /**
    Get the B-matrix (shape function gradients) (in global coordinates).
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \return ``B`` matrix stored per element, per integration point [#nelem, #nne, #tdim, #tdim, #tdim]
    */
    xt::xtensor<double, 6> B() const;

private:

    friend QuadratureBase<QuadratureAxisymmetric>;
    friend QuadratureBaseCartesian<QuadratureAxisymmetric>;

    // qtensor(e, q, i, j) += B(e, q, m, i, j, k) * elemvec(e, q, m, k)
    template <class T, class R>
    void gradN_vector_impl(const T& elemvec, R& qtensor) const;

    template <class T, class R>
    void gradN_vector_T_impl(const T& elemvec, R& qtensor) const;

    template <class T, class R>
    void symGradN_vector_impl(const T& elemvec, R& qtensor) const;

    template <class T, class R>
    void int_N_vector_dV_impl(const T& qvector, R& elemvec) const;

    // elemmat(e, q, m * ndim + i, n * ndim + i) +=
    //     N(e, q, m) * qscalar(e, q) * N(e, q, n) * dV(e, q)
    template <class T, class R>
    void int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const;

    // fm = ( Bm^T : qtensor ) dV
    template <class T, class R>
    void int_gradN_dot_tensor2_dV_impl(const T& qtensor, R& elemvec) const;

    // Kmn = ( Bm^T : qtensor : Bn ) dV
    template <class T, class R>
    void int_gradN_dot_tensor4_dot_gradNT_dV_impl(const T& qtensor, R& elemmat) const;

    void compute_dN_impl();

    constexpr static size_t s_nne = 4;  ///< Number of nodes per element.
    constexpr static size_t s_ndim = 2; ///< Number of dimensions for nodal vectors.
    constexpr static size_t s_tdim = 3; ///< Number of dimensions for tensors.
    size_t m_tdim = 3; ///< Dynamic alias of s_tdim (remove in C++17)
    size_t m_nelem; ///< Number of elements.
    size_t m_nip;   ///< Number of integration points per element.
    xt::xtensor<double, 3> m_x;    ///< nodal positions stored per element [#nelem, #nne, #ndim]
    xt::xtensor<double, 1> m_w;    ///< weight of each integration point [nip]
    xt::xtensor<double, 2> m_xi;   ///< local coordinate of each integration point [#nip, #ndim]
    xt::xtensor<double, 2> m_N;    ///< shape functions [#nip, #nne]
    xt::xtensor<double, 3> m_dNxi; ///< shape function grad. wrt local  coor. [#nip, #nne, #ndim]
    xt::xtensor<double, 2> m_vol;  ///< integration point volume [#nelem, #nip]
    xt::xtensor<double, 6> m_B; ///< B-matrix [#nelem, #nne, #tdim, #tdim, #tdim]
};

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#include "ElementQuad4Axisymmetric.hpp"

#endif

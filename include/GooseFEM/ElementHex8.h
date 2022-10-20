/**
 * Quadrature for 8-noded hexahedral element in 3d (GooseFEM::Mesh::ElementType::Hex8),
 * in a Cartesian coordinate system.
 *
 * @file ElementHex8.h
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#ifndef GOOSEFEM_ELEMENTHEX8_H
#define GOOSEFEM_ELEMENTHEX8_H

#include "config.h"

namespace GooseFEM {
namespace Element {

/**
 * 8-noded hexahedral element in 3d (GooseFEM::Mesh::ElementType::Hex8).
 */
namespace Hex8 {

/**
 * gauss quadrature: quadrature points such that integration is exact for these bi-linear elements::
 */
namespace Gauss {

/**
 * Number of integration points:
 *
 *      nip = nne = 8
 *
 * @return unsigned int
 */
inline size_t nip()
{
    return 8;
}

/**
 * Integration point coordinates (local coordinates).
 *
 * @return Coordinates [#nip, ndim], with `ndim = 3`.
 */
inline array_type::tensor<double, 2> xi()
{
    size_t nip = 8;
    size_t ndim = 3;

    array_type::tensor<double, 2> xi = xt::empty<double>({nip, ndim});

    xi(0, 0) = -1.0 / std::sqrt(3.0);
    xi(0, 1) = -1.0 / std::sqrt(3.0);
    xi(0, 2) = -1.0 / std::sqrt(3.0);

    xi(1, 0) = +1.0 / std::sqrt(3.0);
    xi(1, 1) = -1.0 / std::sqrt(3.0);
    xi(1, 2) = -1.0 / std::sqrt(3.0);

    xi(2, 0) = +1.0 / std::sqrt(3.0);
    xi(2, 1) = +1.0 / std::sqrt(3.0);
    xi(2, 2) = -1.0 / std::sqrt(3.0);

    xi(3, 0) = -1.0 / std::sqrt(3.0);
    xi(3, 1) = +1.0 / std::sqrt(3.0);
    xi(3, 2) = -1.0 / std::sqrt(3.0);

    xi(4, 0) = -1.0 / std::sqrt(3.0);
    xi(4, 1) = -1.0 / std::sqrt(3.0);
    xi(4, 2) = +1.0 / std::sqrt(3.0);

    xi(5, 0) = +1.0 / std::sqrt(3.0);
    xi(5, 1) = -1.0 / std::sqrt(3.0);
    xi(5, 2) = +1.0 / std::sqrt(3.0);

    xi(6, 0) = +1.0 / std::sqrt(3.0);
    xi(6, 1) = +1.0 / std::sqrt(3.0);
    xi(6, 2) = +1.0 / std::sqrt(3.0);

    xi(7, 0) = -1.0 / std::sqrt(3.0);
    xi(7, 1) = +1.0 / std::sqrt(3.0);
    xi(7, 2) = +1.0 / std::sqrt(3.0);

    return xi;
}

/**
 * Integration point weights.
 *
 * @return Coordinates [#nip].
 */
inline array_type::tensor<double, 1> w()
{
    size_t nip = 8;

    array_type::tensor<double, 1> w = xt::empty<double>({nip});

    w(0) = 1.0;
    w(1) = 1.0;
    w(2) = 1.0;
    w(3) = 1.0;
    w(4) = 1.0;
    w(5) = 1.0;
    w(6) = 1.0;
    w(7) = 1.0;

    return w;
}

} // namespace Gauss

/**
 * nodal quadrature: quadrature points coincide with the nodes.
 * The order is the same as in the connectivity.
 */
namespace Nodal {

/**
 * Number of integration points:
 *
 *      nip = nne = 8
 *
 * @return unsigned int
 */
inline size_t nip()
{
    return 8;
}

/**
 * Integration point coordinates (local coordinates).
 *
 * @return Coordinates [#nip, `ndim`], with ``ndim = 3``.
 */
inline array_type::tensor<double, 2> xi()
{
    size_t nip = 8;
    size_t ndim = 3;

    array_type::tensor<double, 2> xi = xt::empty<double>({nip, ndim});

    xi(0, 0) = -1.0;
    xi(0, 1) = -1.0;
    xi(0, 2) = -1.0;

    xi(1, 0) = +1.0;
    xi(1, 1) = -1.0;
    xi(1, 2) = -1.0;

    xi(2, 0) = +1.0;
    xi(2, 1) = +1.0;
    xi(2, 2) = -1.0;

    xi(3, 0) = -1.0;
    xi(3, 1) = +1.0;
    xi(3, 2) = -1.0;

    xi(4, 0) = -1.0;
    xi(4, 1) = -1.0;
    xi(4, 2) = +1.0;

    xi(5, 0) = +1.0;
    xi(5, 1) = -1.0;
    xi(5, 2) = +1.0;

    xi(6, 0) = +1.0;
    xi(6, 1) = +1.0;
    xi(6, 2) = +1.0;

    xi(7, 0) = -1.0;
    xi(7, 1) = +1.0;
    xi(7, 2) = +1.0;

    return xi;
}

/**
 * Integration point weights.
 *
 * @return Coordinates [#nip].
 */
inline array_type::tensor<double, 1> w()
{
    size_t nip = 8;

    array_type::tensor<double, 1> w = xt::empty<double>({nip});

    w(0) = 1.0;
    w(1) = 1.0;
    w(2) = 1.0;
    w(3) = 1.0;
    w(4) = 1.0;
    w(5) = 1.0;
    w(6) = 1.0;
    w(7) = 1.0;

    return w;
}

} // namespace Nodal

/**
 * Interpolation and quadrature.
 *
 * Fixed dimensions:
 * -   `ndim = 3`: number of dimensions.
 * -   `nne = 8`: number of nodes per element.
 *
 * Naming convention:
 * -    `elemmat`:  matrices stored per element, [#nelem, #nne * #ndim, #nne * #ndim]
 * -    `elemvec`:  nodal vectors stored per element, [#nelem, #nne, #ndim]
 * -    `qtensor`:  integration point tensor, [#nelem, #nip, #ndim, #ndim]
 * -    `qscalar`:  integration point scalar, [#nelem, #nip]
 */
class Quadrature : public QuadratureBaseCartesian<Quadrature> {
public:
    Quadrature() = default;

    /**
     * Constructor: use default Gauss integration.
     * The following is pre-computed during construction:
     * -   the shape functions,
     * -   the shape function gradients (in local and global) coordinates,
     * -   the integration points volumes.
     * They can be reused without any cost.
     * They only have to be recomputed when the nodal position changes
     * (note that they are assumed to be constant under a small-strain assumption).
     * In that case use update_x() to update the nodal positions and
     * to recompute the above listed quantities.
     *
     * @param x nodal coordinates (``elemvec``).
     */
    template <class T>
    Quadrature(const T& x) : Quadrature(x, Gauss::xi(), Gauss::w())
    {
    }

    /**
     * Constructor with custom integration.
     * The following is pre-computed during construction:
     * -   the shape functions,
     * -   the shape function gradients (in local and global) coordinates,
     * -   the integration points volumes.
     * They can be reused without any cost.
     * They only have to be recomputed when the nodal position changes
     * (note that they are assumed to be constant under a small-strain assumption).
     * In that case use update_x() to update the nodal positions and
     * to recompute the above listed quantities.
     *
     * @param x nodal coordinates (``elemvec``).
     * @param xi Integration point coordinates (local coordinates) [#nip].
     * @param w Integration point weights [#nip].
     */
    template <class T, class X, class W>
    Quadrature(const T& x, const X& xi, const W& w)
    {
        m_x = x;
        m_w = w;
        m_xi = xi;
        m_nip = w.size();
        m_nelem = m_x.shape(0);
        m_N = xt::empty<double>({m_nip, s_nne});
        m_dNxi = xt::empty<double>({m_nip, s_nne, s_ndim});

        for (size_t q = 0; q < m_nip; ++q) {
            m_N(q, 0) = 0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1)) * (1.0 - xi(q, 2));
            m_N(q, 1) = 0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1)) * (1.0 - xi(q, 2));
            m_N(q, 2) = 0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1)) * (1.0 - xi(q, 2));
            m_N(q, 3) = 0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1)) * (1.0 - xi(q, 2));
            m_N(q, 4) = 0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1)) * (1.0 + xi(q, 2));
            m_N(q, 5) = 0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1)) * (1.0 + xi(q, 2));
            m_N(q, 6) = 0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1)) * (1.0 + xi(q, 2));
            m_N(q, 7) = 0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1)) * (1.0 + xi(q, 2));
        }

        for (size_t q = 0; q < m_nip; ++q) {
            // - dN / dxi_0
            m_dNxi(q, 0, 0) = -0.125 * (1.0 - xi(q, 1)) * (1.0 - xi(q, 2));
            m_dNxi(q, 1, 0) = +0.125 * (1.0 - xi(q, 1)) * (1.0 - xi(q, 2));
            m_dNxi(q, 2, 0) = +0.125 * (1.0 + xi(q, 1)) * (1.0 - xi(q, 2));
            m_dNxi(q, 3, 0) = -0.125 * (1.0 + xi(q, 1)) * (1.0 - xi(q, 2));
            m_dNxi(q, 4, 0) = -0.125 * (1.0 - xi(q, 1)) * (1.0 + xi(q, 2));
            m_dNxi(q, 5, 0) = +0.125 * (1.0 - xi(q, 1)) * (1.0 + xi(q, 2));
            m_dNxi(q, 6, 0) = +0.125 * (1.0 + xi(q, 1)) * (1.0 + xi(q, 2));
            m_dNxi(q, 7, 0) = -0.125 * (1.0 + xi(q, 1)) * (1.0 + xi(q, 2));
            // - dN / dxi_1
            m_dNxi(q, 0, 1) = -0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 2));
            m_dNxi(q, 1, 1) = -0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 2));
            m_dNxi(q, 2, 1) = +0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 2));
            m_dNxi(q, 3, 1) = +0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 2));
            m_dNxi(q, 4, 1) = -0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 2));
            m_dNxi(q, 5, 1) = -0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 2));
            m_dNxi(q, 6, 1) = +0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 2));
            m_dNxi(q, 7, 1) = +0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 2));
            // - dN / dxi_2
            m_dNxi(q, 0, 2) = -0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1));
            m_dNxi(q, 1, 2) = -0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1));
            m_dNxi(q, 2, 2) = -0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1));
            m_dNxi(q, 3, 2) = -0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1));
            m_dNxi(q, 4, 2) = +0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1));
            m_dNxi(q, 5, 2) = +0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1));
            m_dNxi(q, 6, 2) = +0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1));
            m_dNxi(q, 7, 2) = +0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1));
        }

        GOOSEFEM_ASSERT(m_x.shape(1) == s_nne);
        GOOSEFEM_ASSERT(m_x.shape(2) == s_ndim);
        GOOSEFEM_ASSERT(xt::has_shape(m_xi, {m_nip, s_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(m_w, {m_nip}));
        GOOSEFEM_ASSERT(xt::has_shape(m_N, {m_nip, s_nne}));
        GOOSEFEM_ASSERT(xt::has_shape(m_dNxi, {m_nip, s_nne, s_ndim}));

        m_dNx = xt::empty<double>({m_nelem, m_nip, s_nne, s_ndim});
        m_vol = xt::empty<double>(this->shape_qscalar());

        this->compute_dN();
    }

private:
    friend QuadratureBase<Quadrature>;
    friend QuadratureBaseCartesian<Quadrature>;

    template <class T, class R>
    void int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(qscalar, this->shape_qscalar()));
        GOOSEFEM_ASSERT(xt::has_shape(elemmat, this->shape_elemmat()));

        elemmat.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto M = xt::adapt(&elemmat(e, 0, 0), xt::xshape<s_nne * s_ndim, s_nne * s_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto N = xt::adapt(&m_N(q, 0), xt::xshape<s_nne>());
                auto& vol = m_vol(e, q);
                auto& rho = qscalar(e, q);

                // M(m * ndim + i, n * ndim + i) += N(m) * scalar * N(n) * dV
                for (size_t m = 0; m < s_nne; ++m) {
                    for (size_t n = 0; n < s_nne; ++n) {
                        M(m * s_ndim + 0, n * s_ndim + 0) += N(m) * rho * N(n) * vol;
                        M(m * s_ndim + 1, n * s_ndim + 1) += N(m) * rho * N(n) * vol;
                        M(m * s_ndim + 2, n * s_ndim + 2) += N(m) * rho * N(n) * vol;
                    }
                }
            }
        }
    }

    template <class T, class R>
    void int_gradN_dot_tensor2_dV_impl(const T& qtensor, R& elemvec) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));

        elemvec.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto f = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<s_nne, s_ndim>());
                auto sig = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_ndim, s_ndim>());
                auto& v = m_vol(e, q);

                for (size_t m = 0; m < s_nne; ++m) {
                    f(m, 0) +=
                        (dNx(m, 0) * sig(0, 0) + dNx(m, 1) * sig(1, 0) + dNx(m, 2) * sig(2, 0)) * v;
                    f(m, 1) +=
                        (dNx(m, 0) * sig(0, 1) + dNx(m, 1) * sig(1, 1) + dNx(m, 2) * sig(2, 1)) * v;
                    f(m, 2) +=
                        (dNx(m, 0) * sig(0, 2) + dNx(m, 1) * sig(1, 2) + dNx(m, 2) * sig(2, 2)) * v;
                }
            }
        }
    }

    constexpr static size_t s_nne = 8; ///< Number of nodes per element.
    constexpr static size_t s_ndim = 3; ///< Number of dimensions for nodal vectors.
    constexpr static size_t s_tdim = 3; ///< Number of dimensions for tensors.
    size_t m_tdim = 3; ///< Dynamic alias of s_tdim (remove in C++17)
    size_t m_nelem; ///< Number of elements.
    size_t m_nip; ///< Number of integration points per element.
    array_type::tensor<double, 3> m_x; ///< nodal positions stored per element [#nelem, #nne, #ndim]
    array_type::tensor<double, 1> m_w; ///< weight of each integration point [nip]
    array_type::tensor<double, 2> m_xi; ///< local coordinate per integration point [#nip, #ndim]
    array_type::tensor<double, 2> m_N; ///< shape functions [#nip, #nne]
    array_type::tensor<double, 3> m_dNxi; ///< local shape func grad [#nip, #nne, #ndim]
    array_type::tensor<double, 4> m_dNx; ///< global shape func grad  [#nelem, #nip, #nne, #ndim]
    array_type::tensor<double, 2> m_vol; ///< integration point volume [#nelem, #nip]
};

} // namespace Hex8
} // namespace Element
} // namespace GooseFEM

#endif

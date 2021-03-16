/**
Implementation of ElementHex8.h

\file ElementHex8.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENTHEX8_HPP
#define GOOSEFEM_ELEMENTHEX8_HPP

#include "ElementHex8.h"
#include "detail.hpp"

namespace GooseFEM {
namespace Element {
namespace Hex8 {

namespace Gauss {

inline size_t nip()
{
    return 8;
}

inline xt::xtensor<double, 2> xi()
{
    size_t nip = 8;
    size_t ndim = 3;

    xt::xtensor<double, 2> xi = xt::empty<double>({nip, ndim});

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

inline xt::xtensor<double, 1> w()
{
    size_t nip = 8;

    xt::xtensor<double, 1> w = xt::empty<double>({nip});

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

namespace Nodal {

inline size_t nip()
{
    return 8;
}

inline xt::xtensor<double, 2> xi()
{
    size_t nip = 8;
    size_t ndim = 3;

    xt::xtensor<double, 2> xi = xt::empty<double>({nip, ndim});

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

inline xt::xtensor<double, 1> w()
{
    size_t nip = 8;

    xt::xtensor<double, 1> w = xt::empty<double>({nip});

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

inline Quadrature::Quadrature(const xt::xtensor<double, 3>& x)
    : Quadrature(x, Gauss::xi(), Gauss::w())
{
}

inline Quadrature::Quadrature(
    const xt::xtensor<double, 3>& x,
    const xt::xtensor<double, 2>& xi,
    const xt::xtensor<double, 1>& w)
{
    size_t nip = w.size();
    xt::xtensor<double, 2> N = xt::empty<double>({nip, m_nne});
    xt::xtensor<double, 3> dNxi = xt::empty<double>({nip, m_nne, m_ndim});

    for (size_t q = 0; q < nip; ++q) {
        N(q, 0) = 0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1)) * (1.0 - xi(q, 2));
        N(q, 1) = 0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1)) * (1.0 - xi(q, 2));
        N(q, 2) = 0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1)) * (1.0 - xi(q, 2));
        N(q, 3) = 0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1)) * (1.0 - xi(q, 2));
        N(q, 4) = 0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1)) * (1.0 + xi(q, 2));
        N(q, 5) = 0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1)) * (1.0 + xi(q, 2));
        N(q, 6) = 0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1)) * (1.0 + xi(q, 2));
        N(q, 7) = 0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1)) * (1.0 + xi(q, 2));
    }

    for (size_t q = 0; q < nip; ++q) {
        // - dN / dxi_0
        dNxi(q, 0, 0) = -0.125 * (1.0 - xi(q, 1)) * (1.0 - xi(q, 2));
        dNxi(q, 1, 0) = +0.125 * (1.0 - xi(q, 1)) * (1.0 - xi(q, 2));
        dNxi(q, 2, 0) = +0.125 * (1.0 + xi(q, 1)) * (1.0 - xi(q, 2));
        dNxi(q, 3, 0) = -0.125 * (1.0 + xi(q, 1)) * (1.0 - xi(q, 2));
        dNxi(q, 4, 0) = -0.125 * (1.0 - xi(q, 1)) * (1.0 + xi(q, 2));
        dNxi(q, 5, 0) = +0.125 * (1.0 - xi(q, 1)) * (1.0 + xi(q, 2));
        dNxi(q, 6, 0) = +0.125 * (1.0 + xi(q, 1)) * (1.0 + xi(q, 2));
        dNxi(q, 7, 0) = -0.125 * (1.0 + xi(q, 1)) * (1.0 + xi(q, 2));
        // - dN / dxi_1
        dNxi(q, 0, 1) = -0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 2));
        dNxi(q, 1, 1) = -0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 2));
        dNxi(q, 2, 1) = +0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 2));
        dNxi(q, 3, 1) = +0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 2));
        dNxi(q, 4, 1) = -0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 2));
        dNxi(q, 5, 1) = -0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 2));
        dNxi(q, 6, 1) = +0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 2));
        dNxi(q, 7, 1) = +0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 2));
        // - dN / dxi_2
        dNxi(q, 0, 2) = -0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1));
        dNxi(q, 1, 2) = -0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1));
        dNxi(q, 2, 2) = -0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1));
        dNxi(q, 3, 2) = -0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1));
        dNxi(q, 4, 2) = +0.125 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1));
        dNxi(q, 5, 2) = +0.125 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1));
        dNxi(q, 6, 2) = +0.125 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1));
        dNxi(q, 7, 2) = +0.125 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1));
    }

    this->initQuadratureBaseCartesian(x, xi, w, N, dNxi);
}

inline void Quadrature::int_N_scalar_NT_dV(
    const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 3>& elemmat) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qscalar, {m_nelem, m_nip}));
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));

    elemmat.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto M = xt::adapt(&elemmat(e, 0, 0), xt::xshape<m_nne * m_ndim, m_nne * m_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto N = xt::adapt(&m_N(q, 0), xt::xshape<m_nne>());
            auto& vol = m_vol(e, q);
            auto& rho = qscalar(e, q);

            // M(m * ndim + i, n * ndim + i) += N(m) * scalar * N(n) * dV
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t n = 0; n < m_nne; ++n) {
                    M(m * m_ndim + 0, n * m_ndim + 0) += N(m) * rho * N(n) * vol;
                    M(m * m_ndim + 1, n * m_ndim + 1) += N(m) * rho * N(n) * vol;
                    M(m * m_ndim + 2, n * m_ndim + 2) += N(m) * rho * N(n) * vol;
                }
            }
        }
    }
}

inline void Quadrature::int_gradN_dot_tensor2_dV(
    const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, {m_nelem, m_nip, m_ndim, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));

    elemvec.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto f = xt::adapt(&elemvec(e, 0, 0), xt::xshape<m_nne, m_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<m_nne, m_ndim>());
            auto sig = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<m_ndim, m_ndim>());
            auto& vol = m_vol(e, q);

            for (size_t m = 0; m < m_nne; ++m) {
                f(m, 0) +=
                    (dNx(m, 0) * sig(0, 0) + dNx(m, 1) * sig(1, 0) + dNx(m, 2) * sig(2, 0)) * vol;
                f(m, 1) +=
                    (dNx(m, 0) * sig(0, 1) + dNx(m, 1) * sig(1, 1) + dNx(m, 2) * sig(2, 1)) * vol;
                f(m, 2) +=
                    (dNx(m, 0) * sig(0, 2) + dNx(m, 1) * sig(1, 2) + dNx(m, 2) * sig(2, 2)) * vol;
            }
        }
    }
}

} // namespace Hex8
} // namespace Element
} // namespace GooseFEM

#endif

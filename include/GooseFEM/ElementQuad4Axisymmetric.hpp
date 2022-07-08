/**
Implementation of ElementQuad4Axisymmetric.h

\file ElementQuad4Axisymmetric.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENTQUADAXISYMMETRIC_HPP
#define GOOSEFEM_ELEMENTQUADAXISYMMETRIC_HPP

#include "ElementQuad4Axisymmetric.h"
#include "detail.hpp"

namespace GooseFEM {
namespace Element {
namespace Quad4 {

template <class T>
inline QuadratureAxisymmetric::QuadratureAxisymmetric(const T& x)
    : QuadratureAxisymmetric(x, Gauss::xi(), Gauss::w())
{
}

template <class T, class X, class W>
inline QuadratureAxisymmetric::QuadratureAxisymmetric(const T& x, const X& xi, const W& w)
{
    m_x = x;
    m_w = w;
    m_xi = xi;
    m_nip = w.size();
    m_nelem = m_x.shape(0);

    GOOSEFEM_ASSERT(m_x.shape(1) == s_nne);
    GOOSEFEM_ASSERT(m_x.shape(2) == s_ndim);
    GOOSEFEM_ASSERT(xt::has_shape(m_xi, {m_nip, s_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(m_w, {m_nip}));

    m_N = xt::empty<double>({m_nip, s_nne});
    m_dNxi = xt::empty<double>({m_nip, s_nne, s_ndim});
    m_B = xt::empty<double>({m_nelem, m_nip, s_nne, s_tdim, s_tdim, s_tdim});
    m_vol = xt::empty<double>(this->shape_qscalar());

    for (size_t q = 0; q < m_nip; ++q) {
        m_N(q, 0) = 0.25 * (1.0 - m_xi(q, 0)) * (1.0 - m_xi(q, 1));
        m_N(q, 1) = 0.25 * (1.0 + m_xi(q, 0)) * (1.0 - m_xi(q, 1));
        m_N(q, 2) = 0.25 * (1.0 + m_xi(q, 0)) * (1.0 + m_xi(q, 1));
        m_N(q, 3) = 0.25 * (1.0 - m_xi(q, 0)) * (1.0 + m_xi(q, 1));
    }

    for (size_t q = 0; q < m_nip; ++q) {
        // - dN / dxi_0
        m_dNxi(q, 0, 0) = -0.25 * (1.0 - m_xi(q, 1));
        m_dNxi(q, 1, 0) = +0.25 * (1.0 - m_xi(q, 1));
        m_dNxi(q, 2, 0) = +0.25 * (1.0 + m_xi(q, 1));
        m_dNxi(q, 3, 0) = -0.25 * (1.0 + m_xi(q, 1));
        // - dN / dxi_1
        m_dNxi(q, 0, 1) = -0.25 * (1.0 - m_xi(q, 0));
        m_dNxi(q, 1, 1) = -0.25 * (1.0 + m_xi(q, 0));
        m_dNxi(q, 2, 1) = +0.25 * (1.0 + m_xi(q, 0));
        m_dNxi(q, 3, 1) = +0.25 * (1.0 - m_xi(q, 0));
    }

    this->compute_dN_impl();
}

inline void QuadratureAxisymmetric::compute_dN_impl()
{
    // most components remain zero, and are not written
    m_B.fill(0.0);

#pragma omp parallel
    {
        array_type::tensor<double, 2> J = xt::empty<double>({2, 2});
        array_type::tensor<double, 2> Jinv = xt::empty<double>({2, 2});

#pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto x = xt::adapt(&m_x(e, 0, 0), xt::xshape<s_nne, s_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto dNxi = xt::adapt(&m_dNxi(q, 0, 0), xt::xshape<s_nne, s_ndim>());
                auto B =
                    xt::adapt(&m_B(e, q, 0, 0, 0, 0), xt::xshape<s_nne, s_tdim, s_tdim, s_tdim>());
                auto N = xt::adapt(&m_N(q, 0), xt::xshape<s_nne>());

                // J(i,j) += dNxi(m,i) * x(m,j);
                J(0, 0) = dNxi(0, 0) * x(0, 0) + dNxi(1, 0) * x(1, 0) + dNxi(2, 0) * x(2, 0) +
                          dNxi(3, 0) * x(3, 0);
                J(0, 1) = dNxi(0, 0) * x(0, 1) + dNxi(1, 0) * x(1, 1) + dNxi(2, 0) * x(2, 1) +
                          dNxi(3, 0) * x(3, 1);
                J(1, 0) = dNxi(0, 1) * x(0, 0) + dNxi(1, 1) * x(1, 0) + dNxi(2, 1) * x(2, 0) +
                          dNxi(3, 1) * x(3, 0);
                J(1, 1) = dNxi(0, 1) * x(0, 1) + dNxi(1, 1) * x(1, 1) + dNxi(2, 1) * x(2, 1) +
                          dNxi(3, 1) * x(3, 1);

                double Jdet = detail::tensor<2>::inv(J, Jinv);

                // radius for computation of volume
                double rq = N(0) * x(0, 1) + N(1) * x(1, 1) + N(2) * x(2, 1) + N(3) * x(3, 1);

                // dNx(m,i) += Jinv(i,j) * dNxi(m,j)
                for (size_t m = 0; m < s_nne; ++m) {
                    // B(m, r, r, r) = dNdx(m,1)
                    B(m, 0, 0, 0) = Jinv(1, 0) * dNxi(m, 0) + Jinv(1, 1) * dNxi(m, 1);
                    // B(m, r, z, z) = dNdx(m,1)
                    B(m, 0, 2, 2) = Jinv(1, 0) * dNxi(m, 0) + Jinv(1, 1) * dNxi(m, 1);
                    // B(m, t, t, r)
                    B(m, 1, 1, 0) = 1.0 / rq * N(m);
                    // B(m, z, r, r) = dNdx(m,0)
                    B(m, 2, 0, 0) = Jinv(0, 0) * dNxi(m, 0) + Jinv(0, 1) * dNxi(m, 1);
                    // B(m, z, z, z) = dNdx(m,0)
                    B(m, 2, 2, 2) = Jinv(0, 0) * dNxi(m, 0) + Jinv(0, 1) * dNxi(m, 1);
                }

                m_vol(e, q) = m_w(q) * Jdet * 2.0 * M_PI * rq;
            }
        }
    }
}

inline array_type::tensor<double, 6> QuadratureAxisymmetric::B() const
{
    return m_B;
}

template <class T, class R>
inline void QuadratureAxisymmetric::gradN_vector_impl(const T& elemvec, R& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));

    qtensor.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto B = xt::adapt(&m_B(e, q, 0, 0, 0, 0), xt::xshape<s_nne, s_tdim, s_tdim, s_tdim>());
            auto gradu = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_tdim, s_tdim>());

            // gradu(i,j) += B(m,i,j,k) * u(m,perm(k))
            // (where perm(0) = 1, perm(2) = 0)
            gradu(0, 0) = B(0, 0, 0, 0) * u(0, 1) + B(1, 0, 0, 0) * u(1, 1) +
                          B(2, 0, 0, 0) * u(2, 1) + B(3, 0, 0, 0) * u(3, 1);
            gradu(1, 1) = B(0, 1, 1, 0) * u(0, 1) + B(1, 1, 1, 0) * u(1, 1) +
                          B(2, 1, 1, 0) * u(2, 1) + B(3, 1, 1, 0) * u(3, 1);
            gradu(2, 2) = B(0, 2, 2, 2) * u(0, 0) + B(1, 2, 2, 2) * u(1, 0) +
                          B(2, 2, 2, 2) * u(2, 0) + B(3, 2, 2, 2) * u(3, 0);
            gradu(0, 2) = B(0, 0, 2, 2) * u(0, 0) + B(1, 0, 2, 2) * u(1, 0) +
                          B(2, 0, 2, 2) * u(2, 0) + B(3, 0, 2, 2) * u(3, 0);
            gradu(2, 0) = B(0, 2, 0, 0) * u(0, 1) + B(1, 2, 0, 0) * u(1, 1) +
                          B(2, 2, 0, 0) * u(2, 1) + B(3, 2, 0, 0) * u(3, 1);
        }
    }
}

template <class T, class R>
inline void QuadratureAxisymmetric::gradN_vector_T_impl(const T& elemvec, R& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));

    qtensor.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto B = xt::adapt(&m_B(e, q, 0, 0, 0, 0), xt::xshape<s_nne, s_tdim, s_tdim, s_tdim>());
            auto gradu = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_tdim, s_tdim>());

            // gradu(j,i) += B(m,i,j,k) * u(m,perm(k))
            // (where perm(0) = 1, perm(2) = 0)
            gradu(0, 0) = B(0, 0, 0, 0) * u(0, 1) + B(1, 0, 0, 0) * u(1, 1) +
                          B(2, 0, 0, 0) * u(2, 1) + B(3, 0, 0, 0) * u(3, 1);
            gradu(1, 1) = B(0, 1, 1, 0) * u(0, 1) + B(1, 1, 1, 0) * u(1, 1) +
                          B(2, 1, 1, 0) * u(2, 1) + B(3, 1, 1, 0) * u(3, 1);
            gradu(2, 2) = B(0, 2, 2, 2) * u(0, 0) + B(1, 2, 2, 2) * u(1, 0) +
                          B(2, 2, 2, 2) * u(2, 0) + B(3, 2, 2, 2) * u(3, 0);
            gradu(2, 0) = B(0, 0, 2, 2) * u(0, 0) + B(1, 0, 2, 2) * u(1, 0) +
                          B(2, 0, 2, 2) * u(2, 0) + B(3, 0, 2, 2) * u(3, 0);
            gradu(0, 2) = B(0, 2, 0, 0) * u(0, 1) + B(1, 2, 0, 0) * u(1, 1) +
                          B(2, 2, 0, 0) * u(2, 1) + B(3, 2, 0, 0) * u(3, 1);
        }
    }
}

template <class T, class R>
inline void QuadratureAxisymmetric::symGradN_vector_impl(const T& elemvec, R& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));

    qtensor.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto B = xt::adapt(&m_B(e, q, 0, 0, 0, 0), xt::xshape<s_nne, s_tdim, s_tdim, s_tdim>());
            auto eps = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_tdim, s_tdim>());

            // gradu(j,i) += B(m,i,j,k) * u(m,perm(k))
            // eps(j,i) = 0.5 * (gradu(i,j) + gradu(j,i))
            // (where perm(0) = 1, perm(2) = 0)
            eps(0, 0) = B(0, 0, 0, 0) * u(0, 1) + B(1, 0, 0, 0) * u(1, 1) +
                        B(2, 0, 0, 0) * u(2, 1) + B(3, 0, 0, 0) * u(3, 1);
            eps(1, 1) = B(0, 1, 1, 0) * u(0, 1) + B(1, 1, 1, 0) * u(1, 1) +
                        B(2, 1, 1, 0) * u(2, 1) + B(3, 1, 1, 0) * u(3, 1);
            eps(2, 2) = B(0, 2, 2, 2) * u(0, 0) + B(1, 2, 2, 2) * u(1, 0) +
                        B(2, 2, 2, 2) * u(2, 0) + B(3, 2, 2, 2) * u(3, 0);
            eps(2, 0) =
                0.5 * (B(0, 0, 2, 2) * u(0, 0) + B(1, 0, 2, 2) * u(1, 0) + B(2, 0, 2, 2) * u(2, 0) +
                       B(3, 0, 2, 2) * u(3, 0) + B(0, 2, 0, 0) * u(0, 1) + B(1, 2, 0, 0) * u(1, 1) +
                       B(2, 2, 0, 0) * u(2, 1) + B(3, 2, 0, 0) * u(3, 1));
            eps(0, 2) = eps(2, 0);
        }
    }
}

template <class T, class R>
inline void QuadratureAxisymmetric::int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const
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

            // M(m*ndim+i,n*ndim+i) += N(m) * scalar * N(n) * dV
            for (size_t m = 0; m < s_nne; ++m) {
                for (size_t n = 0; n < s_nne; ++n) {
                    M(m * s_ndim + 0, n * s_ndim + 0) += N(m) * rho * N(n) * vol;
                    M(m * s_ndim + 1, n * s_ndim + 1) += N(m) * rho * N(n) * vol;
                }
            }
        }
    }
}

template <class T, class R>
inline void
QuadratureAxisymmetric::int_gradN_dot_tensor2_dV_impl(const T& qtensor, R& elemvec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));

    elemvec.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto f = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto B = xt::adapt(&m_B(e, q, 0, 0, 0, 0), xt::xshape<s_nne, s_tdim, s_tdim, s_tdim>());
            auto sig = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_tdim, s_tdim>());
            auto& vol = m_vol(e, q);

            // f(m,i) += B(m,i,j,perm(k)) * sig(i,j) * dV
            // (where perm(0) = 1, perm(2) = 0)
            for (size_t m = 0; m < s_nne; ++m) {
                f(m, 0) += vol * (B(m, 2, 2, 2) * sig(2, 2) + B(m, 0, 2, 2) * sig(0, 2));
                f(m, 1) += vol * (B(m, 0, 0, 0) * sig(0, 0) + B(m, 1, 1, 0) * sig(1, 1) +
                                  B(m, 2, 0, 0) * sig(2, 0));
            }
        }
    }
}

template <class T, class R>
inline void
QuadratureAxisymmetric::int_gradN_dot_tensor4_dot_gradNT_dV_impl(const T& qtensor, R& elemmat) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<4>()));
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, this->shape_elemmat()));

    elemmat.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto K = xt::adapt(&elemmat(e, 0, 0), xt::xshape<s_nne * s_ndim, s_nne * s_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto B = xt::adapt(&m_B(e, q, 0, 0, 0, 0), xt::xshape<s_nne, s_tdim, s_tdim, s_tdim>());
            auto C =
                xt::adapt(&qtensor(e, q, 0, 0, 0, 0), xt::xshape<s_tdim, s_tdim, s_tdim, s_tdim>());
            auto& vol = m_vol(e, q);

            // K(m*s_ndim+perm(c), n*s_ndim+perm(f)) = B(m,a,b,c) * C(a,b,d,e) * B(n,e,d,f) * vol;
            // (where perm(0) = 1, perm(2) = 0)
            for (size_t m = 0; m < s_nne; ++m) {
                for (size_t n = 0; n < s_nne; ++n) {
                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 0, 0, 0) * C(0, 0, 0, 0) * B(n, 0, 0, 0) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 0, 0, 0) * C(0, 0, 1, 1) * B(n, 1, 1, 0) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 0) +=
                        B(m, 0, 0, 0) * C(0, 0, 2, 2) * B(n, 2, 2, 2) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 0) +=
                        B(m, 0, 0, 0) * C(0, 0, 2, 0) * B(n, 0, 2, 2) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 0, 0, 0) * C(0, 0, 0, 2) * B(n, 2, 0, 0) * vol;

                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 1, 1, 0) * C(1, 1, 0, 0) * B(n, 0, 0, 0) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 1, 1, 0) * C(1, 1, 1, 1) * B(n, 1, 1, 0) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 0) +=
                        B(m, 1, 1, 0) * C(1, 1, 2, 2) * B(n, 2, 2, 2) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 0) +=
                        B(m, 1, 1, 0) * C(1, 1, 2, 0) * B(n, 0, 2, 2) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 1, 1, 0) * C(1, 1, 0, 2) * B(n, 2, 0, 0) * vol;

                    K(m * s_ndim + 0, n * s_ndim + 1) +=
                        B(m, 2, 2, 2) * C(2, 2, 0, 0) * B(n, 0, 0, 0) * vol;
                    K(m * s_ndim + 0, n * s_ndim + 1) +=
                        B(m, 2, 2, 2) * C(2, 2, 1, 1) * B(n, 1, 1, 0) * vol;
                    K(m * s_ndim + 0, n * s_ndim + 0) +=
                        B(m, 2, 2, 2) * C(2, 2, 2, 2) * B(n, 2, 2, 2) * vol;
                    K(m * s_ndim + 0, n * s_ndim + 0) +=
                        B(m, 2, 2, 2) * C(2, 2, 2, 0) * B(n, 0, 2, 2) * vol;
                    K(m * s_ndim + 0, n * s_ndim + 1) +=
                        B(m, 2, 2, 2) * C(2, 2, 0, 2) * B(n, 2, 0, 0) * vol;

                    K(m * s_ndim + 0, n * s_ndim + 1) +=
                        B(m, 0, 2, 2) * C(0, 2, 0, 0) * B(n, 0, 0, 0) * vol;
                    K(m * s_ndim + 0, n * s_ndim + 1) +=
                        B(m, 0, 2, 2) * C(0, 2, 1, 1) * B(n, 1, 1, 0) * vol;
                    K(m * s_ndim + 0, n * s_ndim + 0) +=
                        B(m, 0, 2, 2) * C(0, 2, 2, 2) * B(n, 2, 2, 2) * vol;
                    K(m * s_ndim + 0, n * s_ndim + 0) +=
                        B(m, 0, 2, 2) * C(0, 2, 2, 0) * B(n, 0, 2, 2) * vol;
                    K(m * s_ndim + 0, n * s_ndim + 1) +=
                        B(m, 0, 2, 2) * C(0, 2, 0, 2) * B(n, 2, 0, 0) * vol;

                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 2, 0, 0) * C(2, 0, 0, 0) * B(n, 0, 0, 0) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 2, 0, 0) * C(2, 0, 1, 1) * B(n, 1, 1, 0) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 0) +=
                        B(m, 2, 0, 0) * C(2, 0, 2, 2) * B(n, 2, 2, 2) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 0) +=
                        B(m, 2, 0, 0) * C(2, 0, 2, 0) * B(n, 0, 2, 2) * vol;
                    K(m * s_ndim + 1, n * s_ndim + 1) +=
                        B(m, 2, 0, 0) * C(2, 0, 0, 2) * B(n, 2, 0, 0) * vol;
                }
            }
        }
    }
}

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#endif

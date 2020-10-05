/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_ELEMENTQUAD4PLANAR_HPP
#define GOOSEFEM_ELEMENTQUAD4PLANAR_HPP

#include "ElementQuad4Planar.h"

namespace GooseFEM {
namespace Element {
namespace Quad4 {

inline QuadraturePlanar::QuadraturePlanar(const xt::xtensor<double, 3>& x, double thick)
    : QuadraturePlanar(x, Gauss::xi(), Gauss::w(), thick)
{
}

inline QuadraturePlanar::QuadraturePlanar(
    const xt::xtensor<double, 3>& x,
    const xt::xtensor<double, 2>& xi,
    const xt::xtensor<double, 1>& w,
    double thick)
    : m_x(x), m_w(w), m_xi(xi), m_thick(thick)
{
    GOOSEFEM_ASSERT(m_x.shape(1) == m_nne);
    GOOSEFEM_ASSERT(m_x.shape(2) == m_ndim);

    m_nelem = m_x.shape(0);
    m_nip = m_w.size();

    GOOSEFEM_ASSERT(m_xi.shape(0) == m_nip);
    GOOSEFEM_ASSERT(m_xi.shape(1) == m_ndim);
    GOOSEFEM_ASSERT(m_w.size() == m_nip);

    m_N = xt::empty<double>({m_nip, m_nne});
    m_dNxi = xt::empty<double>({m_nip, m_nne, m_ndim});
    m_dNx = xt::empty<double>({m_nelem, m_nip, m_nne, m_ndim});
    m_vol = xt::empty<double>({m_nelem, m_nip});

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

    compute_dN();
}

inline size_t QuadraturePlanar::nelem() const
{
    return m_nelem;
}

inline size_t QuadraturePlanar::nne() const
{
    return m_nne;
}

inline size_t QuadraturePlanar::ndim() const
{
    return m_ndim;
}

inline size_t QuadraturePlanar::nip() const
{
    return m_nip;
}

inline xt::xtensor<double, 4> QuadraturePlanar::GradN() const
{
    return m_dNx;
}

inline void QuadraturePlanar::dV(xt::xtensor<double, 2>& qscalar) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qscalar, {m_nelem, m_nip}));

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t q = 0; q < m_nip; ++q) {
            qscalar(e, q) = m_vol(e, q);
        }
    }
}

inline void QuadraturePlanar::dV(xt::xtensor<double, 4>& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, {m_nelem, m_nip, m_tdim, m_tdim}));

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t q = 0; q < m_nip; ++q) {
            for (size_t i = 0; i < m_tdim; ++i) {
                for (size_t j = 0; j < m_tdim; ++j) {
                    qtensor(e, q, i, j) = m_vol(e, q);
                }
            }
        }
    }
}

inline void QuadraturePlanar::dV(xt::xarray<double>& qtensor) const
{
    GOOSEFEM_ASSERT(qtensor.shape(0) == m_nelem);
    GOOSEFEM_ASSERT(qtensor.shape(1) == m_nip);

    xt::dynamic_shape<ptrdiff_t> strides = {static_cast<ptrdiff_t>(m_vol.strides()[0]),
                                            static_cast<ptrdiff_t>(m_vol.strides()[1])};

    for (size_t i = 2; i < qtensor.shape().size(); ++i) {
        strides.push_back(0);
    }

    qtensor =
        xt::strided_view(m_vol, qtensor.shape(), std::move(strides), 0ul, xt::layout_type::dynamic);
}

inline void QuadraturePlanar::update_x(const xt::xtensor<double, 3>& x)
{
    GOOSEFEM_ASSERT(x.shape() == m_x.shape());
    xt::noalias(m_x) = x;
    compute_dN();
}

inline void QuadraturePlanar::compute_dN()
{
    #pragma omp parallel
    {
        xt::xtensor<double, 2> J = xt::empty<double>({2, 2});
        xt::xtensor<double, 2> Jinv = xt::empty<double>({2, 2});

        #pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto x = xt::adapt(&m_x(e, 0, 0), xt::xshape<m_nne, m_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto dNxi = xt::adapt(&m_dNxi(q, 0, 0), xt::xshape<m_nne, m_ndim>());
                auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<m_nne, m_ndim>());

                // J(i,j) += dNxi(m,i) * x(m,j);
                J(0, 0) = dNxi(0, 0) * x(0, 0) + dNxi(1, 0) * x(1, 0) + dNxi(2, 0) * x(2, 0) +
                          dNxi(3, 0) * x(3, 0);
                J(0, 1) = dNxi(0, 0) * x(0, 1) + dNxi(1, 0) * x(1, 1) + dNxi(2, 0) * x(2, 1) +
                          dNxi(3, 0) * x(3, 1);
                J(1, 0) = dNxi(0, 1) * x(0, 0) + dNxi(1, 1) * x(1, 0) + dNxi(2, 1) * x(2, 0) +
                          dNxi(3, 1) * x(3, 0);
                J(1, 1) = dNxi(0, 1) * x(0, 1) + dNxi(1, 1) * x(1, 1) + dNxi(2, 1) * x(2, 1) +
                          dNxi(3, 1) * x(3, 1);

                double Jdet = inv(J, Jinv);

                // dNx(m,i) += Jinv(i,j) * dNxi(m,j);
                for (size_t m = 0; m < m_nne; ++m) {
                    dNx(m, 0) = Jinv(0, 0) * dNxi(m, 0) + Jinv(0, 1) * dNxi(m, 1);
                    dNx(m, 1) = Jinv(1, 0) * dNxi(m, 0) + Jinv(1, 1) * dNxi(m, 1);
                }

                m_vol(e, q) = m_w(q) * Jdet * m_thick;
            }
        }
    }
}

inline void QuadraturePlanar::gradN_vector(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, {m_nelem, m_nip, m_tdim, m_tdim}));

    qtensor.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<m_nne, m_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<m_nne, m_ndim>());
            auto gradu = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<m_tdim, m_tdim>());

            // gradu(i,j) += dNx(m,i) * u(m,j)
            gradu(0, 0) = dNx(0, 0) * u(0, 0) + dNx(1, 0) * u(1, 0) + dNx(2, 0) * u(2, 0) +
                          dNx(3, 0) * u(3, 0);
            gradu(0, 1) = dNx(0, 0) * u(0, 1) + dNx(1, 0) * u(1, 1) + dNx(2, 0) * u(2, 1) +
                          dNx(3, 0) * u(3, 1);
            gradu(1, 0) = dNx(0, 1) * u(0, 0) + dNx(1, 1) * u(1, 0) + dNx(2, 1) * u(2, 0) +
                          dNx(3, 1) * u(3, 0);
            gradu(1, 1) = dNx(0, 1) * u(0, 1) + dNx(1, 1) * u(1, 1) + dNx(2, 1) * u(2, 1) +
                          dNx(3, 1) * u(3, 1);
        }
    }
}

inline void QuadraturePlanar::gradN_vector_T(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, {m_nelem, m_nip, m_tdim, m_tdim}));

    qtensor.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<m_nne, m_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<m_nne, m_ndim>());
            auto gradu = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<m_tdim, m_tdim>());

            // gradu(j,i) += dNx(m,i) * u(m,j)
            gradu(0, 0) = dNx(0, 0) * u(0, 0) + dNx(1, 0) * u(1, 0) + dNx(2, 0) * u(2, 0) +
                          dNx(3, 0) * u(3, 0);
            gradu(1, 0) = dNx(0, 0) * u(0, 1) + dNx(1, 0) * u(1, 1) + dNx(2, 0) * u(2, 1) +
                          dNx(3, 0) * u(3, 1);
            gradu(0, 1) = dNx(0, 1) * u(0, 0) + dNx(1, 1) * u(1, 0) + dNx(2, 1) * u(2, 0) +
                          dNx(3, 1) * u(3, 0);
            gradu(1, 1) = dNx(0, 1) * u(0, 1) + dNx(1, 1) * u(1, 1) + dNx(2, 1) * u(2, 1) +
                          dNx(3, 1) * u(3, 1);
        }
    }
}

inline void QuadraturePlanar::symGradN_vector(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, {m_nelem, m_nip, m_tdim, m_tdim}));

    qtensor.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<m_nne, m_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<m_nne, m_ndim>());
            auto eps = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<m_tdim, m_tdim>());

            // gradu(i,j) += dNx(m,i) * u(m,j)
            // eps(j,i) = 0.5 * (gradu(i,j) + gradu(j,i))
            eps(0, 0) = dNx(0, 0) * u(0, 0) + dNx(1, 0) * u(1, 0) + dNx(2, 0) * u(2, 0) +
                        dNx(3, 0) * u(3, 0);
            eps(1, 1) = dNx(0, 1) * u(0, 1) + dNx(1, 1) * u(1, 1) + dNx(2, 1) * u(2, 1) +
                        dNx(3, 1) * u(3, 1);
            eps(0, 1) = 0.5 * (dNx(0, 0) * u(0, 1) + dNx(1, 0) * u(1, 1) + dNx(2, 0) * u(2, 1) +
                               dNx(3, 0) * u(3, 1) + dNx(0, 1) * u(0, 0) + dNx(1, 1) * u(1, 0) +
                               dNx(2, 1) * u(2, 0) + dNx(3, 1) * u(3, 0));
            eps(1, 0) = eps(0, 1);
        }
    }
}

inline void QuadraturePlanar::int_N_scalar_NT_dV(
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

            // M(m*ndim+i,n*ndim+i) += N(m) * scalar * N(n) * dV
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t n = 0; n < m_nne; ++n) {
                    M(m * m_ndim + 0, n * m_ndim + 0) += N(m) * rho * N(n) * vol;
                    M(m * m_ndim + 1, n * m_ndim + 1) += N(m) * rho * N(n) * vol;
                }
            }
        }
    }
}

inline void QuadraturePlanar::int_gradN_dot_tensor2_dV(
    const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, {m_nelem, m_nip, m_tdim, m_tdim}));
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));

    elemvec.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto f = xt::adapt(&elemvec(e, 0, 0), xt::xshape<m_nne, m_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<m_nne, m_ndim>());
            auto sig = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<m_tdim, m_tdim>());
            auto& vol = m_vol(e, q);

            for (size_t m = 0; m < m_nne; ++m) {
                f(m, 0) += (dNx(m, 0) * sig(0, 0) + dNx(m, 1) * sig(1, 0)) * vol;
                f(m, 1) += (dNx(m, 0) * sig(0, 1) + dNx(m, 1) * sig(1, 1)) * vol;
            }
        }
    }
}

inline void QuadraturePlanar::int_gradN_dot_tensor4_dot_gradNT_dV(
    const xt::xtensor<double, 6>& qtensor, xt::xtensor<double, 3>& elemmat) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, {m_nelem, m_nip, m_tdim, m_tdim, m_tdim, m_tdim}));
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));

    elemmat.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto K = xt::adapt(&elemmat(e, 0, 0), xt::xshape<m_nne * m_ndim, m_nne * m_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<m_nne, m_ndim>());
            auto C = xt::adapt(&qtensor(e, q, 0, 0, 0, 0), xt::xshape<m_tdim, m_tdim, m_tdim, m_tdim>());
            auto& vol = m_vol(e, q);

            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t n = 0; n < m_nne; ++n) {
                    for (size_t i = 0; i < m_ndim; ++i) {
                        for (size_t j = 0; j < m_ndim; ++j) {
                            for (size_t k = 0; k < m_ndim; ++k) {
                                for (size_t l = 0; l < m_ndim; ++l) {
                                    K(m * m_ndim + j, n * m_ndim + k) +=
                                        dNx(m, i) * C(i, j, k, l) * dNx(n, l) * vol;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

inline xt::xtensor<double, 2> QuadraturePlanar::DV() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({m_nelem, m_nip});
    this->dV(ret);
    return ret;
}

inline xt::xarray<double> QuadraturePlanar::DV(size_t rank) const
{
    std::vector<size_t> shape = {m_nelem, m_nip};

    for (size_t i = 0; i < rank; ++i) {
        shape.push_back(static_cast<size_t>(m_tdim));
    }

    xt::xarray<double> ret = xt::empty<double>(shape);
    this->dV(ret);
    return ret;
}

inline xt::xtensor<double, 4>
QuadraturePlanar::GradN_vector(const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 4> qtensor = xt::empty<double>({m_nelem, m_nip, m_tdim, m_tdim});
    this->gradN_vector(elemvec, qtensor);
    return qtensor;
}

inline xt::xtensor<double, 4>
QuadraturePlanar::GradN_vector_T(const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 4> qtensor = xt::empty<double>({m_nelem, m_nip, m_tdim, m_tdim});
    this->gradN_vector_T(elemvec, qtensor);
    return qtensor;
}

inline xt::xtensor<double, 4>
QuadraturePlanar::SymGradN_vector(const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 4> qtensor = xt::empty<double>({m_nelem, m_nip, m_tdim, m_tdim});
    this->symGradN_vector(elemvec, qtensor);
    return qtensor;
}

inline xt::xtensor<double, 3>
QuadraturePlanar::Int_N_scalar_NT_dV(const xt::xtensor<double, 2>& qscalar) const
{
    xt::xtensor<double, 3> elemmat = xt::empty<double>({m_nelem, m_nne * m_ndim, m_nne * m_ndim});
    this->int_N_scalar_NT_dV(qscalar, elemmat);
    return elemmat;
}

inline xt::xtensor<double, 3>
QuadraturePlanar::Int_gradN_dot_tensor2_dV(const xt::xtensor<double, 4>& qtensor) const
{
    xt::xtensor<double, 3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});
    this->int_gradN_dot_tensor2_dV(qtensor, elemvec);
    return elemvec;
}

inline xt::xtensor<double, 3>
QuadraturePlanar::Int_gradN_dot_tensor4_dot_gradNT_dV(const xt::xtensor<double, 6>& qtensor) const
{
    xt::xtensor<double, 3> elemmat = xt::empty<double>({m_nelem, m_ndim * m_nne, m_ndim * m_nne});
    this->int_gradN_dot_tensor4_dot_gradNT_dV(qtensor, elemmat);
    return elemmat;
}

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#endif

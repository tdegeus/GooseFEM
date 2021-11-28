/**
Implementation of Element.h

\file Element.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENT_HPP
#define GOOSEFEM_ELEMENT_HPP

#include "Element.h"
#include "detail.hpp"

namespace GooseFEM {
namespace Element {

inline xt::xtensor<double, 3>
asElementVector(const xt::xtensor<size_t, 2>& conn, const xt::xtensor<double, 2>& nodevec)
{
    size_t nelem = conn.shape(0);
    size_t nne = conn.shape(1);
    size_t ndim = nodevec.shape(1);

    xt::xtensor<double, 3> elemvec = xt::empty<double>({nelem, nne, ndim});

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {
        for (size_t m = 0; m < nne; ++m) {
            for (size_t i = 0; i < ndim; ++i) {
                elemvec(e, m, i) = nodevec(conn(e, m), i);
            }
        }
    }

    return elemvec;
}

inline xt::xtensor<double, 2>
assembleNodeVector(const xt::xtensor<size_t, 2>& conn, const xt::xtensor<double, 3>& elemvec)
{
    size_t nelem = conn.shape(0);
    size_t nne = conn.shape(1);
    size_t ndim = elemvec.shape(2);
    size_t nnode = xt::amax(conn)() + 1;

    GOOSEFEM_ASSERT(elemvec.shape(0) == nelem);
    GOOSEFEM_ASSERT(elemvec.shape(1) == nne);

    xt::xtensor<double, 2> nodevec = xt::zeros<double>({nnode, ndim});

    for (size_t e = 0; e < nelem; ++e) {
        for (size_t m = 0; m < nne; ++m) {
            for (size_t i = 0; i < ndim; ++i) {
                nodevec(conn(e, m), i) += elemvec(e, m, i);
            }
        }
    }

    return nodevec;
}

template <class E>
inline bool isSequential(const E& dofs)
{
    size_t ndof = xt::amax(dofs)() + 1;

    xt::xtensor<int, 1> exists = xt::zeros<int>({ndof});

    for (auto& i : dofs) {
        exists[i]++;
    }

    for (auto& i : dofs) {
        if (exists[i] == 0) {
            return false;
        }
    }

    return true;
}

inline bool isDiagonal(const xt::xtensor<double, 3>& elemmat)
{
    GOOSEFEM_ASSERT(elemmat.shape(1) == elemmat.shape(2));

    size_t nelem = elemmat.shape(0);
    size_t N = elemmat.shape(1);

    double eps = std::numeric_limits<double>::epsilon();

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                if (i != j) {
                    if (std::abs(elemmat(e, i, j)) > eps) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

template <class D>
inline auto QuadratureBase<D>::nelem() const
{
    return derived_cast().m_nelem;
}

template <class D>
inline auto QuadratureBase<D>::nne() const
{
    return D::s_nne;
}

template <class D>
inline auto QuadratureBase<D>::ndim() const
{
    return D::s_ndim;
}

template <class D>
inline auto QuadratureBase<D>::tdim() const
{
    return D::s_tdim;
}

template <class D>
inline auto QuadratureBase<D>::nip() const
{
    return derived_cast().m_nip;
}

template <class D>
template <class T, class R>
inline void QuadratureBase<D>::asTensor(const T& arg, R& ret) const
{
    GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_qscalar()));
    GooseFEM::asTensor(arg, ret);
}

template <class D>
template <size_t rank, class T>
inline auto QuadratureBase<D>::AsTensor(const T& arg) const
{
    return GooseFEM::AsTensor<rank>(arg, derived_cast().m_tdim);
}

template <class D>
template <class T>
inline auto QuadratureBase<D>::AsTensor(size_t rank, const T& arg) const
{
    return GooseFEM::AsTensor(rank, arg, derived_cast().m_tdim);
}

template <class D>
inline auto QuadratureBase<D>::shape_elemvec() const -> std::array<size_t, 3>
{
    return std::array<size_t, 3>{derived_cast().m_nelem, D::s_nne, D::s_ndim};
}

template <class D>
inline auto QuadratureBase<D>::shape_elemvec(size_t arg) const -> std::array<size_t, 3>
{
    return std::array<size_t, 3>{derived_cast().m_nelem, D::s_nne, arg};
}

template <class D>
inline auto QuadratureBase<D>::shape_elemmat() const -> std::array<size_t, 3>
{
    return std::array<size_t, 3>{
        derived_cast().m_nelem, D::s_nne * D::s_ndim, D::s_nne * D::s_ndim};
}

template <class D>
template <size_t rank>
inline auto QuadratureBase<D>::shape_qtensor() const -> std::array<size_t, 2 + rank>
{
    std::array<size_t, 2 + rank> shape;
    shape[0] = derived_cast().m_nelem;
    shape[1] = derived_cast().m_nip;
    std::fill(shape.begin() + 2, shape.end(), derived_cast().m_tdim);
    return shape;
}

template <class D>
inline auto QuadratureBase<D>::shape_qtensor(size_t rank) const -> std::vector<size_t>
{
    std::vector<size_t> shape(2 + rank);
    shape[0] = derived_cast().m_nelem;
    shape[1] = derived_cast().m_nip;
    std::fill(shape.begin() + 2, shape.end(), derived_cast().m_tdim);
    return shape;
}

template <class D>
template <size_t trank>
inline auto QuadratureBase<D>::shape_qtensor(size_t rank, size_t arg) const
    -> std::array<size_t, 2 + trank>
{
    GOOSEFEM_ASSERT(trank == rank);
    std::array<size_t, 2 + trank> shape;
    shape[0] = derived_cast().m_nelem;
    shape[1] = derived_cast().m_nip;
    std::fill(shape.begin() + 2, shape.end(), arg);
    return shape;
}

template <class D>
inline auto QuadratureBase<D>::shape_qtensor(size_t rank, size_t arg) const -> std::vector<size_t>
{
    std::vector<size_t> shape(2 + rank);
    shape[0] = derived_cast().m_nelem;
    shape[1] = derived_cast().m_nip;
    std::fill(shape.begin() + 2, shape.end(), arg);
    return shape;
}

template <class D>
inline auto QuadratureBase<D>::shape_qscalar() const -> std::array<size_t, 2>
{
    return std::array<size_t, 2>{derived_cast().m_nelem, derived_cast().m_nip};
}

template <class D>
inline auto QuadratureBase<D>::shape_qvector() const -> std::array<size_t, 3>
{
    return std::array<size_t, 3>{derived_cast().m_nelem, derived_cast().m_nip, D::s_tdim};
}

template <class D>
inline auto QuadratureBase<D>::shape_qvector(size_t arg) const -> std::array<size_t, 3>
{
    return std::array<size_t, 3>{derived_cast().m_nelem, derived_cast().m_nip, arg};
}

template <class D>
template <class R>
inline auto QuadratureBase<D>::allocate_elemvec() const
{
    return xt::xtensor<R, 3>::from_shape(this->shape_elemvec());
}

template <class D>
template <class R>
inline auto QuadratureBase<D>::allocate_elemvec(R val) const
{
    auto ret = xt::xtensor<R, 3>::from_shape(this->shape_elemvec());
    ret.fill(val);
    return ret;
}

template <class D>
template <class R>
inline auto QuadratureBase<D>::allocate_elemmat() const
{
    return xt::xtensor<R, 3>::from_shape(this->shape_elemmat());
}

template <class D>
template <class R>
inline auto QuadratureBase<D>::allocate_elemmat(R val) const
{
    auto ret = xt::xtensor<R, 3>::from_shape(this->shape_elemmat());
    ret.fill(val);
    return ret;
}

template <class D>
template <size_t rank, class R>
inline auto QuadratureBase<D>::allocate_qtensor() const
{
    return xt::xtensor<R, 2 + rank>::from_shape(this->shape_qtensor<rank>());
}

template <class D>
template <size_t rank, class R>
inline auto QuadratureBase<D>::allocate_qtensor(R val) const
{
    auto ret = xt::xtensor<R, 2 + rank>::from_shape(this->shape_qtensor<rank>());
    ret.fill(val);
    return ret;
}

template <class D>
template <class R>
inline auto QuadratureBase<D>::allocate_qtensor(size_t rank) const
{
    return xt::xarray<R>::from_shape(this->shape_qtensor(rank));
}

template <class D>
template <class R>
inline auto QuadratureBase<D>::allocate_qtensor(size_t rank, R val) const
{
    auto ret = xt::xarray<R>::from_shape(this->shape_qtensor(rank));
    ret.fill(val);
    return ret;
}

template <class D>
template <class R>
inline auto QuadratureBase<D>::allocate_qscalar() const
{
    return this->allocate_qtensor<0, R>();
}

template <class D>
template <class R>
inline auto QuadratureBase<D>::allocate_qscalar(R val) const
{
    return this->allocate_qtensor<0, R>(val);
}

template <class D>
inline auto QuadratureBase<D>::derived_cast() -> derived_type&
{
    return *static_cast<derived_type*>(this);
}

template <class D>
inline auto QuadratureBase<D>::derived_cast() const -> const derived_type&
{
    return *static_cast<const derived_type*>(this);
}

template <class D>
inline void QuadratureBaseCartesian<D>::compute_dN()
{
    derived_cast().compute_dN_impl();
}

template <class D>
inline void QuadratureBaseCartesian<D>::compute_dN_impl()
{
    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& vol = derived_cast().m_vol;
    auto& w = derived_cast().m_w;
    auto& dNxi = derived_cast().m_dNxi;
    auto& dNx = derived_cast().m_dNx;
    auto& x = derived_cast().m_x;

    dNx.fill(0.0);

#pragma omp parallel
    {
        auto J = xt::xtensor<double, 2>::from_shape({D::s_ndim, D::s_ndim});
        auto Jinv = xt::xtensor<double, 2>::from_shape({D::s_ndim, D::s_ndim});

#pragma omp for
        for (size_t e = 0; e < nelem; ++e) {

            auto xe = xt::adapt(&x(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

            for (size_t q = 0; q < nip; ++q) {

                auto dNxiq = xt::adapt(&dNxi(q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
                auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

                J.fill(0.0);

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t i = 0; i < D::s_ndim; ++i) {
                        for (size_t j = 0; j < D::s_ndim; ++j) {
                            J(i, j) += dNxiq(m, i) * xe(m, j);
                        }
                    }
                }

                double Jdet = detail::tensor<D::s_ndim>::inv(J, Jinv);

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t i = 0; i < D::s_ndim; ++i) {
                        for (size_t j = 0; j < D::s_ndim; ++j) {
                            dNxq(m, i) += Jinv(i, j) * dNxiq(m, i);
                        }
                    }
                }

                vol(e, q) = w(q) * Jdet;
            }
        }
    }
}

template <class D>
inline auto QuadratureBaseCartesian<D>::GradN() const -> xt::xtensor<double, 4>
{
    return derived_cast().m_dNx;
}

template <class D>
inline auto QuadratureBaseCartesian<D>::dV() const -> xt::xtensor<double, 2>
{
    return derived_cast().m_vol;
}

template <class D>
template <class T>
inline void QuadratureBaseCartesian<D>::update_x(const T& x)
{
    GOOSEFEM_ASSERT(xt::has_shape(x, derived_cast().m_x.shape()));
    xt::noalias(derived_cast().m_x) = x;
    derived_cast().compute_dN_impl();
}

template <class D>
template <class T>
inline auto QuadratureBaseCartesian<D>::InterpQuad_vector(const T& elemvec) const
    -> xt::xtensor<double, 3>
{
    size_t n = elemvec.shape(2);
    auto qvector = xt::xtensor<double, 3>::from_shape(this->shape_qvector(n));
    derived_cast().interpQuad_vector_impl(elemvec, qvector);
    return qvector;
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::interpQuad_vector(const T& elemvec, R& qvector) const
{
    derived_cast().interpQuad_vector_impl(elemvec, qvector);
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::interpQuad_vector_impl(const T& elemvec, R& qvector) const
{
    size_t n = elemvec.shape(2);
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec(n)));
    GOOSEFEM_ASSERT(xt::has_shape(qvector, this->shape_qvector(n)));

    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& N = derived_cast().m_N;

    qvector.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {

        auto fq = &elemvec(e, 0, 0);

        for (size_t q = 0; q < nip; ++q) {

            auto Nq = &N(q, 0);
            auto tq = &qvector(e, q, 0);

            for (size_t m = 0; m < D::s_nne; ++m) {
                for (size_t i = 0; i < n; ++i) {
                    tq[i] += Nq[m] * fq[m * n + i];
                }
            }
        }
    }
}

template <class D>
template <class T>
inline auto QuadratureBaseCartesian<D>::GradN_vector(const T& elemvec) const
    -> xt::xtensor<double, 4>
{
    auto qtensor = xt::xtensor<double, 4>::from_shape(this->template shape_qtensor<2>());
    derived_cast().gradN_vector_impl(elemvec, qtensor);
    return qtensor;
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::gradN_vector(const T& elemvec, R& qtensor) const
{
    derived_cast().gradN_vector_impl(elemvec, qtensor);
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::gradN_vector_impl(const T& elemvec, R& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<2>()));

    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& dNx = derived_cast().m_dNx;

    qtensor.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {

        auto ue = xt::adapt(&elemvec(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

        for (size_t q = 0; q < nip; ++q) {

            auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
            auto graduq = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<D::s_tdim, D::s_tdim>());

            for (size_t m = 0; m < D::s_nne; ++m) {
                for (size_t i = 0; i < D::s_ndim; ++i) {
                    for (size_t j = 0; j < D::s_ndim; ++j) {
                        graduq(i, j) += dNxq(m, i) * ue(m, j);
                    }
                }
            }
        }
    }
}

template <class D>
template <class T>
inline auto QuadratureBaseCartesian<D>::GradN_vector_T(const T& elemvec) const
    -> xt::xtensor<double, 4>
{
    auto qtensor = xt::xtensor<double, 4>::from_shape(this->template shape_qtensor<2>());
    derived_cast().gradN_vector_T_impl(elemvec, qtensor);
    return qtensor;
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::gradN_vector_T(const T& elemvec, R& qtensor) const
{
    derived_cast().gradN_vector_T_impl(elemvec, qtensor);
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::gradN_vector_T_impl(const T& elemvec, R& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<2>()));

    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& dNx = derived_cast().m_dNx;

    qtensor.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {

        auto ue = xt::adapt(&elemvec(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

        for (size_t q = 0; q < nip; ++q) {

            auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
            auto graduq = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<D::s_tdim, D::s_tdim>());

            for (size_t m = 0; m < D::s_nne; ++m) {
                for (size_t i = 0; i < D::s_ndim; ++i) {
                    for (size_t j = 0; j < D::s_ndim; ++j) {
                        graduq(j, i) += dNxq(m, i) * ue(m, j);
                    }
                }
            }
        }
    }
}

template <class D>
template <class T>
inline auto QuadratureBaseCartesian<D>::SymGradN_vector(const T& elemvec) const
    -> xt::xtensor<double, 4>
{
    auto qtensor = xt::xtensor<double, 4>::from_shape(this->template shape_qtensor<2>());
    derived_cast().symGradN_vector_impl(elemvec, qtensor);
    return qtensor;
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::symGradN_vector(const T& elemvec, R& qtensor) const
{
    derived_cast().symGradN_vector_impl(elemvec, qtensor);
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::symGradN_vector_impl(const T& elemvec, R& qtensor) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<2>()));

    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& dNx = derived_cast().m_dNx;

    qtensor.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {

        auto ue = xt::adapt(&elemvec(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

        for (size_t q = 0; q < nip; ++q) {

            auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
            auto epsq = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<D::s_tdim, D::s_tdim>());

            for (size_t m = 0; m < D::s_nne; ++m) {
                for (size_t i = 0; i < D::s_ndim; ++i) {
                    for (size_t j = 0; j < D::s_ndim; ++j) {
                        epsq(i, j) += 0.5 * dNxq(m, i) * ue(m, j);
                        epsq(j, i) += 0.5 * dNxq(m, i) * ue(m, j);
                    }
                }
            }
        }
    }
}

template <class D>
template <class T>
inline auto QuadratureBaseCartesian<D>::Int_N_vector_dV(const T& qvector) const
    -> xt::xtensor<double, 3>
{
    size_t n = qvector.shape(2);
    auto elemvec = xt::xtensor<double, 3>::from_shape(this->shape_elemvec(n));
    derived_cast().int_N_vector_dV_impl(qvector, elemvec);
    return elemvec;
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::int_N_vector_dV(const T& qvector, R& elemvec) const
{
    derived_cast().int_N_vector_dV_impl(qvector, elemvec);
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::int_N_vector_dV_impl(const T& qvector, R& elemvec) const
{
    size_t n = qvector.shape(2);
    GOOSEFEM_ASSERT(xt::has_shape(qvector, this->shape_qvector(n)));
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec(n)));

    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& N = derived_cast().m_N;
    auto& vol = derived_cast().m_vol;

    elemvec.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {

        auto f = &elemvec(e, 0, 0);

        for (size_t q = 0; q < nip; ++q) {

            auto Ne = &N(q, 0);
            auto tq = &qvector(e, q, 0);
            auto& volq = vol(e, q);

            for (size_t m = 0; m < D::s_nne; ++m) {
                for (size_t i = 0; i < n; ++i) {
                    f[m * n + i] += Ne[m] * tq[i] * volq;
                }
            }
        }
    }
}

template <class D>
template <class T>
inline auto QuadratureBaseCartesian<D>::Int_N_scalar_NT_dV(const T& qscalar) const
    -> xt::xtensor<double, 3>
{
    auto elemmat = xt::xtensor<double, 3>::from_shape(this->shape_elemmat());
    derived_cast().int_N_scalar_NT_dV_impl(qscalar, elemmat);
    return elemmat;
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::int_N_scalar_NT_dV(const T& qscalar, R& elemmat) const
{
    derived_cast().int_N_scalar_NT_dV_impl(qscalar, elemmat);
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qscalar, this->shape_qscalar()));
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, this->shape_elemmat()));

    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& N = derived_cast().m_N;
    auto& vol = derived_cast().m_vol;

    elemmat.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {

        auto Me =
            xt::adapt(&elemmat(e, 0, 0), xt::xshape<D::s_nne * D::s_ndim, D::s_nne * D::s_ndim>());

        for (size_t q = 0; q < nip; ++q) {

            auto Ne = xt::adapt(&N(q, 0), xt::xshape<D::s_nne>());
            auto& volq = vol(e, q);
            auto& rho = qscalar(e, q);

            // M(m * D::s_ndim + i, n * D::s_ndim + i) += N(m) * scalar * N(n) * dV
            for (size_t m = 0; m < D::s_nne; ++m) {
                for (size_t n = 0; n < D::s_nne; ++n) {
                    for (size_t i = 0; i < D::s_ndim; ++i) {
                        Me(m * D::s_ndim + i, n * D::s_ndim + i) += Ne(m) * rho * Ne(n) * volq;
                    }
                }
            }
        }
    }
}

template <class D>
template <class T>
inline auto QuadratureBaseCartesian<D>::Int_gradN_dot_tensor2_dV(const T& qtensor) const
    -> xt::xtensor<double, 3>
{
    auto elemvec = xt::xtensor<double, 3>::from_shape(this->shape_elemvec());
    derived_cast().int_gradN_dot_tensor2_dV_impl(qtensor, elemvec);
    return elemvec;
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::int_gradN_dot_tensor2_dV(const T& qtensor, R& elemvec) const
{
    derived_cast().int_gradN_dot_tensor2_dV_impl(qtensor, elemvec);
}

template <class D>
template <class T, class R>
inline void
QuadratureBaseCartesian<D>::int_gradN_dot_tensor2_dV_impl(const T& qtensor, R& elemvec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<2>()));
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));

    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& dNx = derived_cast().m_m_dNx;
    auto& vol = derived_cast().m_vol;

    elemvec.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {

        auto fe = xt::adapt(&elemvec(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

        for (size_t q = 0; q < nip; ++q) {

            auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
            auto sigq = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<D::s_tdim, D::s_tdim>());
            auto& volq = vol(e, q);

            for (size_t m = 0; m < D::s_nne; ++m) {
                for (size_t i = 0; i < D::s_ndim; ++i) {
                    for (size_t j = 0; j < D::s_ndim; ++j) {
                        fe(m, j) += dNxq(m, i) * sigq(i, j) * volq;
                    }
                }
            }
        }
    }
}

template <class D>
template <class T>
inline auto QuadratureBaseCartesian<D>::Int_gradN_dot_tensor4_dot_gradNT_dV(const T& qtensor) const
    -> xt::xtensor<double, 3>
{
    auto elemmat = xt::xtensor<double, 3>::from_shape(this->shape_elemmat());
    derived_cast().int_gradN_dot_tensor4_dot_gradNT_dV_impl(qtensor, elemmat);
    return elemmat;
}

template <class D>
template <class T, class R>
inline void
QuadratureBaseCartesian<D>::int_gradN_dot_tensor4_dot_gradNT_dV(const T& qtensor, R& elemmat) const
{
    derived_cast().int_gradN_dot_tensor4_dot_gradNT_dV_impl(qtensor, elemmat);
}

template <class D>
template <class T, class R>
inline void QuadratureBaseCartesian<D>::int_gradN_dot_tensor4_dot_gradNT_dV_impl(
    const T& qtensor,
    R& elemmat) const
{
    GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<4>()));
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, this->shape_elemmat()));

    auto nelem = derived_cast().m_nelem;
    auto nip = derived_cast().m_nip;
    auto& dNx = derived_cast().m_dNx;
    auto& vol = derived_cast().m_vol;

    elemmat.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {

        auto K =
            xt::adapt(&elemmat(e, 0, 0), xt::xshape<D::s_nne * D::s_ndim, D::s_nne * D::s_ndim>());

        for (size_t q = 0; q < nip; ++q) {

            auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
            auto Cq = xt::adapt(
                &qtensor(e, q, 0, 0, 0, 0),
                xt::xshape<D::s_tdim, D::s_tdim, D::s_tdim, D::s_tdim>());
            auto& volq = vol(e, q);

            for (size_t m = 0; m < D::s_nne; ++m) {
                for (size_t n = 0; n < D::s_nne; ++n) {
                    for (size_t i = 0; i < D::s_ndim; ++i) {
                        for (size_t j = 0; j < D::s_ndim; ++j) {
                            for (size_t k = 0; k < D::s_ndim; ++k) {
                                for (size_t l = 0; l < D::s_ndim; ++l) {
                                    K(m * D::s_ndim + j, n * D::s_ndim + k) +=
                                        dNxq(m, i) * Cq(i, j, k, l) * dNxq(n, l) * volq;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <class D>
inline auto QuadratureBaseCartesian<D>::derived_cast() -> derived_type&
{
    return *static_cast<derived_type*>(this);
}

template <class D>
inline auto QuadratureBaseCartesian<D>::derived_cast() const -> const derived_type&
{
    return *static_cast<const derived_type*>(this);
}

} // namespace Element
} // namespace GooseFEM

#endif

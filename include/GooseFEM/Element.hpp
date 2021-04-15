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

inline xt::xtensor<double, 3> asElementVector(
    const xt::xtensor<size_t, 2>& conn, const xt::xtensor<double, 2>& nodevec)
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

inline xt::xtensor<double, 2> assembleNodeVector(
    const xt::xtensor<size_t, 2>& conn, const xt::xtensor<double, 3>& elemvec)
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

template <size_t ne, size_t nd, size_t td>
inline QuadratureBase<ne, nd, td>::QuadratureBase(size_t nelem, size_t nip)
{
    this->initQuadratureBase(nelem, nip);
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBase<ne, nd, td>::initQuadratureBase(size_t nelem, size_t nip)
{
    m_nelem = nelem;
    m_nip = nip;
}

template <size_t ne, size_t nd, size_t td>
inline size_t QuadratureBase<ne, nd, td>::nelem() const
{
    return m_nelem;
}

template <size_t ne, size_t nd, size_t td>
inline size_t QuadratureBase<ne, nd, td>::nne() const
{
    return m_nne;
}

template <size_t ne, size_t nd, size_t td>
inline size_t QuadratureBase<ne, nd, td>::ndim() const
{
    return m_ndim;
}

template <size_t ne, size_t nd, size_t td>
inline size_t QuadratureBase<ne, nd, td>::tdim() const
{
    return m_tdim;
}

template <size_t ne, size_t nd, size_t td>
inline size_t QuadratureBase<ne, nd, td>::nip() const
{
    return m_nip;
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank, class T>
inline void QuadratureBase<ne, nd, td>::asTensor(
    const xt::xtensor<T, 2>& arg,
    xt::xtensor<T, 2 + rank>& ret) const
{
    GOOSEFEM_ASSERT(xt::has_shape(arg, {m_nelem, m_nne}));
    GooseFEM::asTensor<2, rank>(arg, ret);
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank, class T>
inline xt::xtensor<T, 2 + rank> QuadratureBase<ne, nd, td>::AsTensor(
    const xt::xtensor<T, 2>& qscalar) const
{
    return GooseFEM::AsTensor<2, rank>(qscalar, m_tdim);
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xarray<T> QuadratureBase<ne, nd, td>::AsTensor(
    size_t rank,
    const xt::xtensor<T, 2>& qscalar) const
{
    return GooseFEM::AsTensor(rank, qscalar, m_tdim);
}

template <size_t ne, size_t nd, size_t td>
inline std::array<size_t, 3> QuadratureBase<ne, nd, td>::shape_elemvec() const
{
    return std::array<size_t, 3>{m_nelem, m_nne, m_ndim};
}

template <size_t ne, size_t nd, size_t td>
inline std::array<size_t, 3> QuadratureBase<ne, nd, td>::shape_elemmat() const
{
    return std::array<size_t, 3>{m_nelem, m_nne * m_ndim, m_nne * m_ndim};
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank>
inline std::array<size_t, 2 + rank> QuadratureBase<ne, nd, td>::shape_qtensor() const
{
    std::array<size_t, 2 + rank> shape;
    shape[0] = m_nelem;
    shape[1] = m_nip;
    std::fill(shape.begin() + 2, shape.end(), td);
    return shape;
}

template <size_t ne, size_t nd, size_t td>
inline std::vector<size_t> QuadratureBase<ne, nd, td>::shape_qtensor(size_t rank) const
{
    std::vector<size_t> shape(2 + rank);
    shape[0] = m_nelem;
    shape[1] = m_nip;
    std::fill(shape.begin() + 2, shape.end(), td);
    return shape;
}

template <size_t ne, size_t nd, size_t td>
inline std::vector<size_t> QuadratureBase<ne, nd, td>::shape_qscalar() const
{
    std::vector<size_t> shape(2);
    shape[0] = m_nelem;
    shape[1] = m_nip;
    return shape;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 3> QuadratureBase<ne, nd, td>::allocate_elemvec() const
{
    xt::xtensor<T, 3> ret = xt::empty<T>(this->shape_elemvec());
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 3> QuadratureBase<ne, nd, td>::allocate_elemvec(T val) const
{
    xt::xtensor<T, 3> ret = xt::empty<T>(this->shape_elemvec());
    ret.fill(val);
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 3> QuadratureBase<ne, nd, td>::allocate_elemmat() const
{
    xt::xtensor<T, 3> ret = xt::empty<T>(this->shape_elemmat());
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 3> QuadratureBase<ne, nd, td>::allocate_elemmat(T val) const
{
    xt::xtensor<T, 3> ret = xt::empty<T>(this->shape_elemmat());
    ret.fill(val);
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank, class T>
inline xt::xtensor<T, 2 + rank> QuadratureBase<ne, nd, td>::allocate_qtensor() const
{
    xt::xtensor<T, 2 + rank> ret = xt::empty<T>(this->shape_qtensor<rank>());
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank, class T>
inline xt::xtensor<T, 2 + rank> QuadratureBase<ne, nd, td>::allocate_qtensor(T val) const
{
    xt::xtensor<T, 2 + rank> ret = xt::empty<T>(this->shape_qtensor<rank>());
    ret.fill(val);
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xarray<T> QuadratureBase<ne, nd, td>::allocate_qtensor(size_t rank) const
{
    xt::xarray<T> ret = xt::empty<T>(this->shape_qtensor(rank));
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xarray<T> QuadratureBase<ne, nd, td>::allocate_qtensor(size_t rank, T val) const
{
    xt::xarray<T> ret = xt::empty<T>(this->shape_qtensor(rank));
    ret.fill(val);
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 2> QuadratureBase<ne, nd, td>::allocate_qscalar() const
{
    return this->allocate_qtensor<0, T>();
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 2> QuadratureBase<ne, nd, td>::allocate_qscalar(T val) const
{
    return this->allocate_qtensor<0, T>(val);
}

/**
\cond
*/

template <size_t ne, size_t nd, size_t td>
template <size_t rank>
inline std::array<size_t, 2 + rank> QuadratureBase<ne, nd, td>::ShapeQtensor() const
{
    GOOSEFEM_WARNING("Deprecation warning: ShapeQtensor<rank> -> shape_qtensor<rank>");
    return this->shape_qtensor<rank>();
}

template <size_t ne, size_t nd, size_t td>
inline std::vector<size_t> QuadratureBase<ne, nd, td>::ShapeQtensor(size_t rank) const
{
    GOOSEFEM_WARNING("Deprecation warning: ShapeQtensor(rank) -> shape_qtensor(rank)");
    return this->shape_qtensor(rank);
}

template <size_t ne, size_t nd, size_t td>
inline std::vector<size_t> QuadratureBase<ne, nd, td>::ShapeQscalar() const
{
    GOOSEFEM_WARNING("Deprecation warning: ShapeQscalar -> shape_qscalar");
    return this->shape_qscalar();
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank, class T>
inline xt::xtensor<T, 2 + rank> QuadratureBase<ne, nd, td>::AllocateQtensor() const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateQtensor<rank, T> -> allocate_qtensor<rank, T>");
    return this->allocate_qtensor<rank, T>();
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank, class T>
inline xt::xtensor<T, 2 + rank> QuadratureBase<ne, nd, td>::AllocateQtensor(T val) const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateQtensor<rank, T> -> allocate_qtensor<rank, T>");
    return this->allocate_qtensor<rank, T>(val);
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xarray<T> QuadratureBase<ne, nd, td>::AllocateQtensor(size_t rank) const
{
    GOOSEFEM_WARNING_PYTHON("Deprecation warning: use np.empty(this.shape_qtensor(rank))")
    GOOSEFEM_WARNING("Deprecation warning: AllocateQtensor(rank) -> allocate_qtensor(rank)");
    return this->allocate_qtensor<T>(rank);
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xarray<T> QuadratureBase<ne, nd, td>::AllocateQtensor(size_t rank, T val) const
{
    GOOSEFEM_WARNING_PYTHON("Deprecation warning: use val * np.ones(this.shape_qtensor(rank))")
    GOOSEFEM_WARNING("Deprecation warning: AllocateQtensor(rank) -> allocate_qtensor(rank)");
    return this->allocate_qtensor<T>(rank, val);
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 2> QuadratureBase<ne, nd, td>::AllocateQscalar() const
{
    GOOSEFEM_WARNING_PYTHON("Deprecation warning: use np.empty(this.shape_qscalar())")
    GOOSEFEM_WARNING("Deprecation warning: AllocateQscalar -> allocate_qscalar");
    return this->allocate_qtensor<0, T>();
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 2> QuadratureBase<ne, nd, td>::AllocateQscalar(T val) const
{
    GOOSEFEM_WARNING_PYTHON("Deprecation warning: use np.empty(this.shape_qscalar())")
    GOOSEFEM_WARNING("Deprecation warning: AllocateQscalar -> allocate_qscalar");
    return this->allocate_qtensor<0, T>(val);
}

/**
\endcond
*/

template <size_t ne, size_t nd, size_t td>
inline QuadratureBaseCartesian<ne, nd, td>::QuadratureBaseCartesian(
    const xt::xtensor<double, 3>& x,
    const xt::xtensor<double, 2>& xi,
    const xt::xtensor<double, 1>& w,
    const xt::xtensor<double, 2>& N,
    const xt::xtensor<double, 3>& dNdxi)
{
    this->initQuadratureBaseCartesian(x, xi, w, N, dNdxi);
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::initQuadratureBaseCartesian(
    const xt::xtensor<double, 3>& x,
    const xt::xtensor<double, 2>& xi,
    const xt::xtensor<double, 1>& w,
    const xt::xtensor<double, 2>& N,
    const xt::xtensor<double, 3>& dNdxi)
{
    m_x = x;
    m_w = w;
    m_xi = xi;
    m_N = N;
    m_dNxi = dNdxi;

    this->initQuadratureBase(m_x.shape(0), m_w.size());

    GOOSEFEM_ASSERT(m_x.shape(1) == m_nne);
    GOOSEFEM_ASSERT(m_x.shape(2) == m_ndim);
    GOOSEFEM_ASSERT(xt::has_shape(m_xi, {m_nip, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(m_w, {m_nip}));
    GOOSEFEM_ASSERT(xt::has_shape(m_N, {m_nip, m_nne}));
    GOOSEFEM_ASSERT(xt::has_shape(m_dNxi, {m_nip, m_nne, m_ndim}));

    m_dNx = xt::empty<double>({m_nelem, m_nip, m_nne, m_ndim});
    m_vol = xt::empty<double>(this->shape_qscalar());

    this->compute_dN();
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::compute_dN()
{
    #pragma omp parallel
    {
        xt::xtensor<double, 2> J = xt::empty<double>({m_ndim, m_ndim});
        xt::xtensor<double, 2> Jinv = xt::empty<double>({m_ndim, m_ndim});
        m_dNx.fill(0.0);

        #pragma omp for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto x = xt::adapt(&m_x(e, 0, 0), xt::xshape<m_nne, m_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto dNxi = xt::adapt(&m_dNxi(q, 0, 0), xt::xshape<m_nne, m_ndim>());
                auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<m_nne, m_ndim>());

                J.fill(0.0);

                for (size_t m = 0; m < m_nne; ++m) {
                    for (size_t i = 0; i < m_ndim; ++i) {
                        for (size_t j = 0; j < m_ndim; ++j) {
                            J(i, j) += dNxi(m, i) * x(m, j);
                        }
                    }
                }

                double Jdet = detail::tensor<m_ndim>::inv(J, Jinv);

                for (size_t m = 0; m < m_nne; ++m) {
                    for (size_t i = 0; i < m_ndim; ++i) {
                        for (size_t j = 0; j < m_ndim; ++j) {
                            dNx(m, i) += Jinv(i, j) * dNxi(m, i);
                        }
                    }
                }

                m_vol(e, q) = m_w(q) * Jdet;
            }
        }
    }
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 4> QuadratureBaseCartesian<ne, nd, td>::GradN() const
{
    return m_dNx;
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 2> QuadratureBaseCartesian<ne, nd, td>::dV() const
{
    return m_vol;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::update_x(const xt::xtensor<double, 3>& x)
{
    GOOSEFEM_ASSERT(x.shape() == m_x.shape());
    xt::noalias(m_x) = x;
    this->compute_dN();
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 3> QuadratureBaseCartesian<ne, nd, td>::InterpQuad_vector(
    const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 3> qvector = xt::empty<double>({m_nelem, m_nip, m_ndim});
    this->interpQuad_vector(elemvec, qvector);
    return qvector;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::interpQuad_vector(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 3>& qvector) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(qvector, {m_nelem, m_nip, m_ndim}));

    qvector.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<m_nne, m_ndim>());

        for (size_t q = 0; q < m_nip; ++q) {

            auto N = xt::adapt(&m_N(q, 0), xt::xshape<m_nne>());
            auto ui = xt::adapt(&qvector(e, q, 0), xt::xshape<m_ndim>());

            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    ui(i) += N(m) * u(m, i);
                }
            }
        }
    }
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 3> QuadratureBaseCartesian<ne, nd, td>::Interp_N_vector(
    const xt::xtensor<double, 3>& elemvec) const
{
    GOOSEFEM_WARNING("Deprecation warning: Interp_N_vector -> InterpQuad_vector");
    GOOSEFEM_WARNING_PYTHON("Deprecation warning: Interp_N_vector -> InterpQuad_vector")
    return this->InterpQuad_vector(elemvec);
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::interp_N_vector(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 3>& qvector) const
{
    GOOSEFEM_WARNING("Deprecation warning: interp_N_vector -> interpQuad_vector");
    GOOSEFEM_WARNING_PYTHON("Deprecation warning: interp_N_vector -> interpQuad_vector")
    this->interpQuad_vector(elemvec, qvector);
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 4> QuadratureBaseCartesian<ne, nd, td>::GradN_vector(
    const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 4> qtensor = xt::empty<double>({m_nelem, m_nip, m_tdim, m_tdim});
    this->gradN_vector(elemvec, qtensor);
    return qtensor;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::gradN_vector(
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

            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    for (size_t j = 0; j < m_ndim; ++j) {
                        gradu(i, j) += dNx(m, i) * u(m, j);
                    }
                }
            }
        }
    }
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 4> QuadratureBaseCartesian<ne, nd, td>::GradN_vector_T(
    const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 4> qtensor = xt::empty<double>({m_nelem, m_nip, m_tdim, m_tdim});
    this->gradN_vector_T(elemvec, qtensor);
    return qtensor;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::gradN_vector_T(
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

            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    for (size_t j = 0; j < m_ndim; ++j) {
                        gradu(j, i) += dNx(m, i) * u(m, j);
                    }
                }
            }
        }
    }
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 4> QuadratureBaseCartesian<ne, nd, td>::SymGradN_vector(
    const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 4> qtensor = xt::empty<double>({m_nelem, m_nip, m_tdim, m_tdim});
    this->symGradN_vector(elemvec, qtensor);
    return qtensor;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::symGradN_vector(
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

            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    for (size_t j = 0; j < m_ndim; ++j) {
                        eps(i, j) += 0.5 * dNx(m, i) * u(m, j);
                        eps(j, i) += 0.5 * dNx(m, i) * u(m, j);
                    }
                }
            }
        }
    }
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 3> QuadratureBaseCartesian<ne, nd, td>::Int_N_vector_dV(
    const xt::xtensor<double, 3>& qvector) const
{
    size_t n = qvector.shape(2);
    xt::xtensor<double, 3> elemvec = xt::empty<double>({m_nelem, m_nne, n});
    this->int_N_vector_dV(qvector, elemvec);
    return elemvec;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::int_N_vector_dV(
    const xt::xtensor<double, 3>& qvector, xt::xtensor<double, 3>& elemvec) const
{
    size_t n = qvector.shape(2);
    GOOSEFEM_ASSERT(xt::has_shape(qvector, {m_nelem, m_nip, n}));
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, n}));

    elemvec.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {

        auto f = &elemvec(e, 0, 0);

        for (size_t q = 0; q < m_nip; ++q) {

            auto N = &m_N(q, 0);
            auto t = &qvector(e, q, 0);
            auto& vol = m_vol(e, q);

            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < n; ++i) {
                    f[m * n + i] += N[m] * t[i] * vol;
                }
            }
        }
    }
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 3> QuadratureBaseCartesian<ne, nd, td>::Int_N_scalar_NT_dV(
    const xt::xtensor<double, 2>& qscalar) const
{
    xt::xtensor<double, 3> elemmat = xt::empty<double>({m_nelem, m_nne * m_ndim, m_nne * m_ndim});
    this->int_N_scalar_NT_dV(qscalar, elemmat);
    return elemmat;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::int_N_scalar_NT_dV(
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
                    for (size_t i = 0; i < m_ndim; ++i) {
                        M(m * m_ndim + i, n * m_ndim + i) += N(m) * rho * N(n) * vol;
                    }
                }
            }
        }
    }
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 3> QuadratureBaseCartesian<ne, nd, td>::Int_gradN_dot_tensor2_dV(
    const xt::xtensor<double, 4>& qtensor) const
{
    xt::xtensor<double, 3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});
    this->int_gradN_dot_tensor2_dV(qtensor, elemvec);
    return elemvec;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::int_gradN_dot_tensor2_dV(
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
                for (size_t i = 0; i < m_ndim; ++i) {
                    for (size_t j = 0; j < m_ndim; ++j) {
                        f(m, j) += dNx(m, i) * sig(i, j) * vol;
                    }
                }
            }
        }
    }
}

template <size_t ne, size_t nd, size_t td>
inline xt::xtensor<double, 3> QuadratureBaseCartesian<ne, nd, td>::Int_gradN_dot_tensor4_dot_gradNT_dV(
    const xt::xtensor<double, 6>& qtensor) const
{
    xt::xtensor<double, 3> elemmat = xt::empty<double>({m_nelem, m_ndim * m_nne, m_ndim * m_nne});
    this->int_gradN_dot_tensor4_dot_gradNT_dV(qtensor, elemmat);
    return elemmat;
}

template <size_t ne, size_t nd, size_t td>
inline void QuadratureBaseCartesian<ne, nd, td>::int_gradN_dot_tensor4_dot_gradNT_dV(
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


} // namespace Element
} // namespace GooseFEM

#endif

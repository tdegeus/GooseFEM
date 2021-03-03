/**
Implementation of Element.h

\file Element.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENT_HPP
#define GOOSEFEM_ELEMENT_HPP

#include "Element.h"

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
template <size_t rank>
inline std::array<size_t, 2 + rank> QuadratureBase<ne, nd, td>::ShapeQtensor() const
{
    std::array<size_t, 2 + rank> shape;
    shape[0] = m_nelem;
    shape[1] = m_nip;
    std::fill(shape.begin() + 2, shape.end(), td);
    return shape;
}

template <size_t ne, size_t nd, size_t td>
inline std::vector<size_t> QuadratureBase<ne, nd, td>::ShapeQtensor(size_t rank) const
{
    std::vector<size_t> shape(2 + rank);
    shape[0] = m_nelem;
    shape[1] = m_nip;
    std::fill(shape.begin() + 2, shape.end(), td);
    return shape;
}

template <size_t ne, size_t nd, size_t td>
inline std::vector<size_t> QuadratureBase<ne, nd, td>::ShapeQscalar() const
{
    std::vector<size_t> shape(2);
    shape[0] = m_nelem;
    shape[1] = m_nip;
    return shape;
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank, class T>
inline xt::xtensor<T, 2 + rank> QuadratureBase<ne, nd, td>::AllocateQtensor() const
{
    xt::xtensor<T, 2 + rank> ret = xt::empty<T>(this->ShapeQtensor<rank>());
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <size_t rank, class T>
inline xt::xtensor<T, 2 + rank> QuadratureBase<ne, nd, td>::AllocateQtensor(T val) const
{
    xt::xtensor<T, 2 + rank> ret = xt::empty<T>(this->ShapeQtensor<rank>());
    ret.fill(val);
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xarray<T> QuadratureBase<ne, nd, td>::AllocateQtensor(size_t rank) const
{
    xt::xarray<T> ret = xt::empty<T>(this->ShapeQtensor(rank));
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xarray<T> QuadratureBase<ne, nd, td>::AllocateQtensor(size_t rank, T val) const
{
    xt::xarray<T> ret = xt::empty<T>(this->ShapeQtensor(rank));
    ret.fill(val);
    return ret;
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 2> QuadratureBase<ne, nd, td>::AllocateQscalar() const
{
    return this->AllocateQtensor<0, T>();
}

template <size_t ne, size_t nd, size_t td>
template <class T>
inline xt::xtensor<T, 2> QuadratureBase<ne, nd, td>::AllocateQscalar(T val) const
{
    return this->AllocateQtensor<0, T>(val);
}

} // namespace Element
} // namespace GooseFEM

#endif

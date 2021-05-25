/**
Implementation of Vector.h

\file Vector.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VECTOR_HPP
#define GOOSEFEM_VECTOR_HPP

#include "Vector.h"

namespace GooseFEM {

template <class S, class T>
inline Vector::Vector(const S& conn, const T& dofs)
{
    GOOSEFEM_ASSERT(conn.dimension() == 2);
    GOOSEFEM_ASSERT(dofs.dimension() == 2);

    m_conn = conn;
    m_dofs = dofs;

    m_nelem = m_conn.shape(0);
    m_nne = m_conn.shape(1);
    m_nnode = m_dofs.shape(0);
    m_ndim = m_dofs.shape(1);
    m_ndof = xt::amax(m_dofs)() + 1;

    GOOSEFEM_ASSERT(xt::amax(m_conn)() + 1 <= m_nnode);
    GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
}

inline size_t Vector::nelem() const
{
    return m_nelem;
}

inline size_t Vector::nne() const
{
    return m_nne;
}

inline size_t Vector::nnode() const
{
    return m_nnode;
}

inline size_t Vector::ndim() const
{
    return m_ndim;
}

inline size_t Vector::ndof() const
{
    return m_ndof;
}

inline xt::xtensor<size_t, 2> Vector::conn() const
{
    return m_conn;
}

inline xt::xtensor<size_t, 2> Vector::dofs() const
{
    return m_dofs;
}

template <class T>
inline void Vector::copy(const T& nodevec_src, T& nodevec_dest) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec_src, this->shape_nodevec()));
    GOOSEFEM_ASSERT(xt::has_shape(nodevec_dest, this->shape_nodevec()));

    xt::noalias(nodevec_dest) = nodevec_src;
}

template <class T, class R>
void Vector::asDofs(const T& arg, R& dofval) const
{
    this->asDofs_impl(arg, dofval);
}

// asDofs : distribution

template <class T, class R>
inline
typename std::enable_if_t<std::equal<xt::get_rank<T>::value, 1>>
Vector::asDofs_impl(const T& arg, R& dofval) const
{
    dofval = arg;
}

template <class T, class R>
inline
typename std::enable_if_t<std::equal<xt::get_rank<T>::value, 2>>
Vector::asDofs_impl(const T& arg, R& dofval) const
{
    this->asDofs_impl_nodevec(arg, dofval);
}

template <class T, class R>
inline
typename std::enable_if_t<std::equal<xt::get_rank<T>::value, 3>>
Vector::asDofs_impl(const T& arg, R& dofval) const
{
    this->asDofs_impl_elemvec(arg, dofval);
}

template <class T, class R>
inline
typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value>
Vector::asDofs_impl(const T& arg, R& dofval) const
{
    if (arg.dimension() == 1) {
        dofval = arg;
    }
    else if (arg.dimension() == 2) {
        this->asDofs_impl_nodevec(arg, dofval);
    }
    else if (arg.dimension() == 3) {
        this->asDofs_impl_elemvec(arg, dofval);
    }
    else {
        throw std::runtime_error("Vector::asDofs unknown dimension first argument");
    }
}

template <class T, class R>
inline void Vector::asNode(const T& arg, R& nodevec) const
{
    if (arg.dimension() == 1) {
        this->asNode_dofval(arg, nodevec);
    }
    else if (arg.dimension() == 3) {
        this->asNode_elemvec(arg, nodevec);
    }
    else if (arg.dimension() == 2) {
        nodevec = arg;
    }
    else {
        throw std::runtime_error("Vector::asNode unknown dimension first argument");
    }
}

template <class T, class R>
inline void Vector::asElement(const T& arg, R& elemvec) const
{
    if (arg.dimension() == 1) {
        this->asElement_dofval(arg, elemvec);
    }
    else if (arg.dimension() == 2) {
        this->asElement_nodevec(arg, elemvec);
    }
    else if (arg.dimension() == 3) {
        elemvec = arg;
    }
    else {
        throw std::runtime_error("Vector::asElement unknown dimension first argument");
    }
}

template <class T, class R>
inline void Vector::assembleDofs(const T& arg, R& dofval) const
{
    if (arg.dimension() == 2) {
        this->assembleDofs_nodevec(arg, dofval);
    }
    else if (arg.dimension() == 3) {
        this->assembleDofs_elemvec(arg, dofval);
    }
    else if (arg.dimension() == 1) {
        dofval = arg;
    }
    else {
        throw std::runtime_error("Vector::assembleDofs unknown dimension first argument");
    }
}

template <class T, class R>
inline void Vector::assembleNode(const T& arg, R& nodevec) const
{
    if (arg.dimension() == 3) {
        this->assembleNode_elemvec(arg, nodevec);
    }
    else if (arg.dimension() == 2) {
        nodevec = arg;
    }
    else {
        throw std::runtime_error("Vector::assembleNode unknown dimension first argument");
    }
}

// asDofs : implementation

template <class T, class R>
inline void Vector::asDofs_impl_dofval(const T& arg, R& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_dofval()));
    GOOSEFEM_ASSERT(xt::has_shape(dofval, this->shape_dofval()));

    dofval = arg;
}

template <class T, class R>
inline void Vector::asDofs_impl_nodevec(const T& arg, R& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_nodevec()));
    GOOSEFEM_ASSERT(xt::has_shape(dofval, this->shape_dofval()));

    dofval.fill(0.0);

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            dofval(m_dofs(m, i)) = arg(m, i);
        }
    }
}

template <class T, class R>
inline void Vector::asDofs_impl_elemvec(const T& arg, R& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(dofval, this->shape_dofval()));

    dofval.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                dofval(m_dofs(m_conn(e, m), i)) = arg(e, m, i);
            }
        }
    }
}

template <class T, class R>
inline void Vector::asNode_dofval(const T& dofval, R& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(dofval, this->shape_dofval()));
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            nodevec(m, i) = dofval(m_dofs(m, i));
        }
    }
}

template <class T, class R>
inline void Vector::asNode_elemvec(const T& elemvec, R& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));

    nodevec.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                nodevec(m_conn(e, m), i) = elemvec(e, m, i);
            }
        }
    }
}

template <class T, class R>
inline void Vector::asElement_dofval(const T& dofval, R& elemvec) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                elemvec(e, m, i) = dofval(m_dofs(m_conn(e, m), i));
            }
        }
    }
}

template <class T, class R>
inline void Vector::asElement_nodevec(const T& nodevec,R& elemvec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                elemvec(e, m, i) = nodevec(m_conn(e, m), i);
            }
        }
    }
}

template <class T, class R>
inline void Vector::assembleDofs_nodevec(const T& nodevec, R& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));
    GOOSEFEM_ASSERT(xt::has_shape(dofval, this->shape_dofval()));

    dofval.fill(0.0);

    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            dofval(m_dofs(m, i)) += nodevec(m, i);
        }
    }
}

template <class T, class R>
inline void Vector::assembleDofs_elemvec(const T& elemvec, R& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(dofval, this->shape_dofval()));

    dofval.fill(0.0);

    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                dofval(m_dofs(m_conn(e, m), i)) += elemvec(e, m, i);
            }
        }
    }
}

template <class T, class R>
inline void Vector::assembleNode_elemvec(const T& elemvec, R& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));

    xt::xtensor<double, 1> dofval = this->AssembleDofs(elemvec);
    this->asNode(dofval, nodevec);
}

template <class T>
inline xt::xtensor<double, 1> Vector::AsDofs(const T& arg) const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
    this->asDofs(arg, dofval);
    return dofval;
}

template <class T>
inline xt::xtensor<double, 2> Vector::AsNode(const T& arg) const
{
    xt::xtensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
    this->asNode(arg, nodevec);
    return nodevec;
}

template <class T>
inline xt::xtensor<double, 3> Vector::AsElement(const T& arg) const
{
    xt::xtensor<double, 3> elemvec = xt::empty<double>(this->shape_elemvec());
    this->asElement(arg, elemvec);
    return elemvec;
}

template <class T>
inline xt::xtensor<double, 1> Vector::AssembleDofs(const T& arg) const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
    this->assembleDofs(arg, dofval);
    return dofval;
}

template <class T>
inline xt::xtensor<double, 2> Vector::AssembleNode(const T& arg) const
{
    xt::xtensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
    this->assembleNode(arg, nodevec);
    return nodevec;
}

template <class T>
inline T Vector::Copy(const T& nodevec_src, const T& nodevec_dest) const
{
    T ret = T::from_shape(nodevec_dest.shape());
    this->copy(nodevec_src, ret);
    return ret;
}

inline std::array<size_t, 1> Vector::shape_dofval() const
{
    std::array<size_t, 1> shape;
    shape[0] = m_ndof;
    return shape;
}

inline std::array<size_t, 2> Vector::shape_nodevec() const
{
    std::array<size_t, 2> shape;
    shape[0] = m_nnode;
    shape[1] = m_ndim;
    return shape;
}

inline std::array<size_t, 3> Vector::shape_elemvec() const
{
    std::array<size_t, 3> shape;
    shape[0] = m_nelem;
    shape[1] = m_nne;
    shape[2] = m_ndim;
    return shape;
}

inline std::array<size_t, 3> Vector::shape_elemmat() const
{
    std::array<size_t, 3> shape;
    shape[0] = m_nelem;
    shape[1] = m_nne * m_ndim;
    shape[2] = m_nne * m_ndim;
    return shape;
}

inline xt::xtensor<double, 1> Vector::allocate_dofval() const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
    return dofval;
}

inline xt::xtensor<double, 2> Vector::allocate_nodevec() const
{
    xt::xtensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
    return nodevec;
}

inline xt::xtensor<double, 3> Vector::allocate_elemvec() const
{
    xt::xtensor<double, 3> elemvec = xt::empty<double>(this->shape_elemvec());
    return elemvec;
}

inline xt::xtensor<double, 3> Vector::allocate_elemmat() const
{
    xt::xtensor<double, 3> elemmat = xt::empty<double>(this->shape_elemmat());
    return elemmat;
}

inline xt::xtensor<double, 1> Vector::allocate_dofval(double val) const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
    dofval.fill(val);
    return dofval;
}

inline xt::xtensor<double, 2> Vector::allocate_nodevec(double val) const
{
    xt::xtensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
    nodevec.fill(val);
    return nodevec;
}

inline xt::xtensor<double, 3> Vector::allocate_elemvec(double val) const
{
    xt::xtensor<double, 3> elemvec = xt::empty<double>(this->shape_elemvec());
    elemvec.fill(val);
    return elemvec;
}

inline xt::xtensor<double, 3> Vector::allocate_elemmat(double val) const
{
    xt::xtensor<double, 3> elemmat = xt::empty<double>(this->shape_elemmat());
    elemmat.fill(val);
    return elemmat;
}

} // namespace GooseFEM

#endif

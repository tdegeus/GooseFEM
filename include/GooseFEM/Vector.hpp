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

inline Vector::Vector(const xt::xtensor<size_t, 2>& conn, const xt::xtensor<size_t, 2>& dofs)
    : m_conn(conn), m_dofs(dofs)
{
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

inline void Vector::copy(
    const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec_src, this->shape_nodevec()));
    GOOSEFEM_ASSERT(xt::has_shape(nodevec_dest, this->shape_nodevec()));

    xt::noalias(nodevec_dest) = nodevec_src;
}

inline void Vector::asDofs(
    const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    dofval.fill(0.0);

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            dofval(m_dofs(m, i)) = nodevec(m, i);
        }
    }
}

inline void Vector::asDofs(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    dofval.fill(0.0);

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                dofval(m_dofs(m_conn(e, m), i)) = elemvec(e, m, i);
            }
        }
    }
}

inline void Vector::asNode(
    const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            nodevec(m, i) = dofval(m_dofs(m, i));
        }
    }
}

inline void Vector::asNode(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 2>& nodevec) const
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

inline void Vector::asElement(
    const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 3>& elemvec) const
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

inline void Vector::asElement(
    const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 3>& elemvec) const
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

inline void Vector::assembleDofs(
    const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    dofval.fill(0.0);

    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            dofval(m_dofs(m, i)) += nodevec(m, i);
        }
    }
}

inline void Vector::assembleDofs(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    dofval.fill(0.0);

    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                dofval(m_dofs(m_conn(e, m), i)) += elemvec(e, m, i);
            }
        }
    }
}

inline void Vector::assembleNode(
    const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, this->shape_nodevec()));

    xt::xtensor<double, 1> dofval = this->AssembleDofs(elemvec);
    this->asNode(dofval, nodevec);
}

inline xt::xtensor<double, 1> Vector::AsDofs(const xt::xtensor<double, 2>& nodevec) const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
    this->asDofs(nodevec, dofval);
    return dofval;
}

inline xt::xtensor<double, 1> Vector::AsDofs(const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
    this->asDofs(elemvec, dofval);
    return dofval;
}

inline xt::xtensor<double, 2> Vector::AsNode(const xt::xtensor<double, 1>& dofval) const
{
    xt::xtensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
    this->asNode(dofval, nodevec);
    return nodevec;
}

inline xt::xtensor<double, 2> Vector::AsNode(const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
    this->asNode(elemvec, nodevec);
    return nodevec;
}

inline xt::xtensor<double, 3> Vector::AsElement(const xt::xtensor<double, 1>& dofval) const
{
    xt::xtensor<double, 3> elemvec = xt::empty<double>(this->shape_elemvec());
    this->asElement(dofval, elemvec);
    return elemvec;
}

inline xt::xtensor<double, 3> Vector::AsElement(const xt::xtensor<double, 2>& nodevec) const
{
    xt::xtensor<double, 3> elemvec = xt::empty<double>(this->shape_elemvec());
    this->asElement(nodevec, elemvec);
    return elemvec;
}

inline xt::xtensor<double, 1> Vector::AssembleDofs(const xt::xtensor<double, 2>& nodevec) const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
    this->assembleDofs(nodevec, dofval);
    return dofval;
}

inline xt::xtensor<double, 1> Vector::AssembleDofs(const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
    this->assembleDofs(elemvec, dofval);
    return dofval;
}

inline xt::xtensor<double, 2> Vector::AssembleNode(const xt::xtensor<double, 3>& elemvec) const
{
    xt::xtensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
    this->assembleNode(elemvec, nodevec);
    return nodevec;
}

inline xt::xtensor<double, 2> Vector::Copy(
    const xt::xtensor<double, 2>& nodevec_src, const xt::xtensor<double, 2>& nodevec_dest) const
{
    xt::xtensor<double, 2> ret = nodevec_dest;
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

/**
\cond
*/

inline std::array<size_t, 1> Vector::ShapeDofval() const
{
    GOOSEFEM_WARNING("Deprecation warning: ShapeDofval -> shape_dofval");
    return this->shape_dofval();
}

inline std::array<size_t, 2> Vector::ShapeNodevec() const
{
    GOOSEFEM_WARNING("Deprecation warning: ShapeNodevec -> shape_nodevec");
    return this->shape_nodevec();
}

inline std::array<size_t, 3> Vector::ShapeElemvec() const
{
    GOOSEFEM_WARNING("Deprecation warning: ShapeElemvec -> shape_elemvec");
    return this->shape_elemvec();
}

inline std::array<size_t, 3> Vector::ShapeElemmat() const
{
    GOOSEFEM_WARNING("Deprecation warning: ShapeElemmat -> shape_elemmat");
    return this->shape_elemmat();
}

inline xt::xtensor<double, 1> Vector::AllocateDofval() const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateDofval -> allocate_dofval");
    GOOSEFEM_WARNING_PYTHON("Deprecation warnings: using np.empty(this.allocate_dofval())")
    return this->allocate_dofval();
}

inline xt::xtensor<double, 2> Vector::AllocateNodevec() const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateNodevec -> allocate_nodevec");
    GOOSEFEM_WARNING_PYTHON("Deprecation warnings: using np.empty(this.allocate_nodevec())")
    return this->allocate_nodevec();
}

inline xt::xtensor<double, 3> Vector::AllocateElemvec() const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateElemvec -> allocate_elemvec");
    GOOSEFEM_WARNING_PYTHON("Deprecation warnings: using np.empty(this.allocate_elemvec())")
    return this->allocate_elemvec();
}

inline xt::xtensor<double, 3> Vector::AllocateElemmat() const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateElemmat -> allocate_elemmat");
    GOOSEFEM_WARNING_PYTHON("Deprecation warnings: using np.empty(this.allocate_elemmat())")
    return this->allocate_elemmat();
}

inline xt::xtensor<double, 1> Vector::AllocateDofval(double val) const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateDofval -> allocate_dofval");
    GOOSEFEM_WARNING_PYTHON("Deprecation warnings: using val * np.ones(this.allocate_dofval())")
    return this->allocate_dofval(val);
}

inline xt::xtensor<double, 2> Vector::AllocateNodevec(double val) const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateNodevec -> allocate_nodevec");
    GOOSEFEM_WARNING_PYTHON("Deprecation warnings: using val * np.ones(this.allocate_nodevec())")
    return this->allocate_nodevec(val);
}

inline xt::xtensor<double, 3> Vector::AllocateElemvec(double val) const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateElemvec -> allocate_elemvec");
    GOOSEFEM_WARNING_PYTHON("Deprecation warnings: using val * np.ones(this.allocate_elemvec())")
    return this->allocate_elemvec(val);
}

inline xt::xtensor<double, 3> Vector::AllocateElemmat(double val) const
{
    GOOSEFEM_WARNING("Deprecation warning: AllocateElemmat -> allocate_elemmat");
    GOOSEFEM_WARNING_PYTHON("Deprecation warnings: using val * np.ones(this.allocate_elemmat())")
    return this->allocate_elemmat(val);
}

/**
\endcond
*/

} // namespace GooseFEM

#endif

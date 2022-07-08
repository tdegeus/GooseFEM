/**
Implementation of VectorPartitioned.h

\file VectorPartitioned.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VECTORPARTITIONED_HPP
#define GOOSEFEM_VECTORPARTITIONED_HPP

#include "Mesh.h"
#include "VectorPartitioned.h"
#include "assertions.h"

namespace GooseFEM {

inline VectorPartitioned::VectorPartitioned(
    const array_type::tensor<size_t, 2>& conn,
    const array_type::tensor<size_t, 2>& dofs,
    const array_type::tensor<size_t, 1>& iip)
    : Vector(conn, dofs), m_iip(iip)
{
    GOOSEFEM_ASSERT(is_unique(iip));

    m_iiu = xt::setdiff1d(m_dofs, m_iip);
    m_nnp = m_iip.size();
    m_nnu = m_iiu.size();
    m_part = Mesh::Reorder({m_iiu, m_iip}).apply(m_dofs);

    GOOSEFEM_ASSERT(xt::amax(m_iip)() <= xt::amax(m_dofs)());
}

inline size_t VectorPartitioned::nnu() const
{
    return m_nnu;
}

inline size_t VectorPartitioned::nnp() const
{
    return m_nnp;
}

inline array_type::tensor<size_t, 1> VectorPartitioned::iiu() const
{
    return m_iiu;
}

inline array_type::tensor<size_t, 1> VectorPartitioned::iip() const
{
    return m_iip;
}

inline array_type::tensor<bool, 2> VectorPartitioned::dofs_is_u() const
{
    array_type::tensor<bool, 2> ret = xt::zeros<bool>(this->shape_nodevec());

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) < m_nnu) {
                ret(m, i) = true;
            }
        }
    }

    return ret;
}

inline array_type::tensor<bool, 2> VectorPartitioned::dofs_is_p() const
{
    array_type::tensor<bool, 2> ret = xt::zeros<bool>(this->shape_nodevec());

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) >= m_nnu) {
                ret(m, i) = true;
            }
        }
    }

    return ret;
}

inline void VectorPartitioned::copy_u(
    const array_type::tensor<double, 2>& nodevec_src,
    array_type::tensor<double, 2>& nodevec_dest) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec_src, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(nodevec_dest, {m_nnode, m_ndim}));

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) < m_nnu) {
                nodevec_dest(m, i) = nodevec_src(m, i);
            }
        }
    }
}

inline void VectorPartitioned::copy_p(
    const array_type::tensor<double, 2>& nodevec_src,
    array_type::tensor<double, 2>& nodevec_dest) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec_src, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(nodevec_dest, {m_nnode, m_ndim}));

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) >= m_nnu) {
                nodevec_dest(m, i) = nodevec_src(m, i);
            }
        }
    }
}

inline void VectorPartitioned::dofsFromParitioned(
    const array_type::tensor<double, 1>& dofval_u,
    const array_type::tensor<double, 1>& dofval_p,
    array_type::tensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);
    GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    dofval.fill(0.0);

#pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        dofval(m_iiu(d)) = dofval_u(d);
    }

#pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        dofval(m_iip(d)) = dofval_p(d);
    }
}

inline void VectorPartitioned::asDofs_u(
    const array_type::tensor<double, 1>& dofval,
    array_type::tensor<double, 1>& dofval_u) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);
    GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);

#pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        dofval_u(d) = dofval(m_iiu(d));
    }
}

inline void VectorPartitioned::asDofs_u(
    const array_type::tensor<double, 2>& nodevec,
    array_type::tensor<double, 1>& dofval_u) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);

    dofval_u.fill(0.0);

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) < m_nnu) {
                dofval_u(m_part(m, i)) = nodevec(m, i);
            }
        }
    }
}

inline void VectorPartitioned::asDofs_p(
    const array_type::tensor<double, 1>& dofval,
    array_type::tensor<double, 1>& dofval_p) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);
    GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);

#pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        dofval_p(d) = dofval(m_iip(d));
    }
}

inline void VectorPartitioned::asDofs_p(
    const array_type::tensor<double, 2>& nodevec,
    array_type::tensor<double, 1>& dofval_p) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);

    dofval_p.fill(0.0);

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) >= m_nnu) {
                dofval_p(m_part(m, i) - m_nnu) = nodevec(m, i);
            }
        }
    }
}

inline void VectorPartitioned::asDofs_u(
    const array_type::tensor<double, 3>& elemvec,
    array_type::tensor<double, 1>& dofval_u) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));
    GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);

    dofval_u.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m_conn(e, m), i) < m_nnu) {
                    dofval_u(m_part(m_conn(e, m), i)) = elemvec(e, m, i);
                }
            }
        }
    }
}

inline void VectorPartitioned::asDofs_p(
    const array_type::tensor<double, 3>& elemvec,
    array_type::tensor<double, 1>& dofval_p) const
{
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));
    GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);

    dofval_p.fill(0.0);

#pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m_conn(e, m), i) >= m_nnu) {
                    dofval_p(m_part(m_conn(e, m), i) - m_nnu) = elemvec(e, m, i);
                }
            }
        }
    }
}

inline void VectorPartitioned::nodeFromPartitioned(
    const array_type::tensor<double, 1>& dofval_u,
    const array_type::tensor<double, 1>& dofval_p,
    array_type::tensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);
    GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) < m_nnu) {
                nodevec(m, i) = dofval_u(m_part(m, i));
            }
            else {
                nodevec(m, i) = dofval_p(m_part(m, i) - m_nnu);
            }
        }
    }
}

inline void VectorPartitioned::elementFromPartitioned(
    const array_type::tensor<double, 1>& dofval_u,
    const array_type::tensor<double, 1>& dofval_p,
    array_type::tensor<double, 3>& elemvec) const
{
    GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);
    GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);
    GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));

#pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m_conn(e, m), i) < m_nnu) {
                    elemvec(e, m, i) = dofval_u(m_part(m_conn(e, m), i));
                }
                else {
                    elemvec(e, m, i) = dofval_p(m_part(m_conn(e, m), i) - m_nnu);
                }
            }
        }
    }
}

inline array_type::tensor<double, 1> VectorPartitioned::DofsFromParitioned(
    const array_type::tensor<double, 1>& dofval_u,
    const array_type::tensor<double, 1>& dofval_p) const
{
    array_type::tensor<double, 1> dofval = xt::empty<double>({m_ndof});
    this->dofsFromParitioned(dofval_u, dofval_p, dofval);
    return dofval;
}

inline array_type::tensor<double, 1>
VectorPartitioned::AsDofs_u(const array_type::tensor<double, 1>& dofval) const
{
    array_type::tensor<double, 1> dofval_u = xt::empty<double>({m_nnu});
    this->asDofs_u(dofval, dofval_u);
    return dofval_u;
}

inline array_type::tensor<double, 1>
VectorPartitioned::AsDofs_u(const array_type::tensor<double, 2>& nodevec) const
{
    array_type::tensor<double, 1> dofval_u = xt::empty<double>({m_nnu});
    this->asDofs_u(nodevec, dofval_u);
    return dofval_u;
}

inline array_type::tensor<double, 1>
VectorPartitioned::AsDofs_p(const array_type::tensor<double, 1>& dofval) const
{
    array_type::tensor<double, 1> dofval_p = xt::empty<double>({m_nnp});
    this->asDofs_p(dofval, dofval_p);
    return dofval_p;
}

inline array_type::tensor<double, 1>
VectorPartitioned::AsDofs_p(const array_type::tensor<double, 2>& nodevec) const
{
    array_type::tensor<double, 1> dofval_p = xt::empty<double>({m_nnp});
    this->asDofs_p(nodevec, dofval_p);
    return dofval_p;
}

inline array_type::tensor<double, 1>
VectorPartitioned::AsDofs_u(const array_type::tensor<double, 3>& elemvec) const
{
    array_type::tensor<double, 1> dofval_u = xt::empty<double>({m_nnu});
    this->asDofs_u(elemvec, dofval_u);
    return dofval_u;
}

inline array_type::tensor<double, 1>
VectorPartitioned::AsDofs_p(const array_type::tensor<double, 3>& elemvec) const
{
    array_type::tensor<double, 1> dofval_p = xt::empty<double>({m_nnp});
    this->asDofs_p(elemvec, dofval_p);
    return dofval_p;
}

inline array_type::tensor<double, 2> VectorPartitioned::NodeFromPartitioned(
    const array_type::tensor<double, 1>& dofval_u,
    const array_type::tensor<double, 1>& dofval_p) const
{
    array_type::tensor<double, 2> nodevec = xt::empty<double>({m_nnode, m_ndim});
    this->nodeFromPartitioned(dofval_u, dofval_p, nodevec);
    return nodevec;
}

inline array_type::tensor<double, 3> VectorPartitioned::ElementFromPartitioned(
    const array_type::tensor<double, 1>& dofval_u,
    const array_type::tensor<double, 1>& dofval_p) const
{
    array_type::tensor<double, 3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});
    this->elementFromPartitioned(dofval_u, dofval_p, elemvec);
    return elemvec;
}

inline array_type::tensor<double, 2> VectorPartitioned::Copy_u(
    const array_type::tensor<double, 2>& nodevec_src,
    const array_type::tensor<double, 2>& nodevec_dest) const
{
    array_type::tensor<double, 2> ret = nodevec_dest;
    this->copy_u(nodevec_src, ret);
    return ret;
}

inline array_type::tensor<double, 2> VectorPartitioned::Copy_p(
    const array_type::tensor<double, 2>& nodevec_src,
    const array_type::tensor<double, 2>& nodevec_dest) const
{
    array_type::tensor<double, 2> ret = nodevec_dest;
    this->copy_p(nodevec_src, ret);
    return ret;
}

} // namespace GooseFEM

#endif

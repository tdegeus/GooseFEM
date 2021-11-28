/**
Implementation of VectorPartitionedTyings.h

\file VectorPartitionedTyings.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VECTORPARTITIONEDTYINGS_HPP
#define GOOSEFEM_VECTORPARTITIONEDTYINGS_HPP

#include "VectorPartitionedTyings.h"

namespace GooseFEM {

template <class E, class M>
inline VectorPartitionedTyings::VectorPartitionedTyings(
    const E& conn,
    const E& dofs,
    const M& Cdu,
    const M& Cdp,
    const M& Cdi)
    : Vector(conn, dofs), m_Cdu(Cdu), m_Cdp(Cdp), m_Cdi(Cdi)
{
    GOOSEFEM_ASSERT(Cdu.rows() == Cdp.rows());
    GOOSEFEM_ASSERT(Cdi.rows() == Cdp.rows());

    m_nnu = static_cast<size_t>(m_Cdu.cols());
    m_nnp = static_cast<size_t>(m_Cdp.cols());
    m_nnd = static_cast<size_t>(m_Cdp.rows());
    m_nni = m_nnu + m_nnp;
    m_iiu = xt::arange<size_t>(m_nnu);
    m_iip = xt::arange<size_t>(m_nnu, m_nnu + m_nnp);
    m_iid = xt::arange<size_t>(m_nni, m_nni + m_nnd);
    m_Cud = m_Cdu.transpose();
    m_Cpd = m_Cdp.transpose();
    m_Cid = m_Cdi.transpose();

    GOOSEFEM_ASSERT(static_cast<size_t>(m_Cdi.cols()) == m_nni);
    GOOSEFEM_ASSERT(m_ndof == xt::amax(m_dofs)() + 1);
}

inline size_t VectorPartitionedTyings::nnu() const
{
    return m_nnu;
}

inline size_t VectorPartitionedTyings::nnp() const
{
    return m_nnp;
}

inline size_t VectorPartitionedTyings::nni() const
{
    return m_nni;
}

inline size_t VectorPartitionedTyings::nnd() const
{
    return m_nnd;
}

inline xt::xtensor<size_t, 1> VectorPartitionedTyings::iiu() const
{
    return m_iiu;
}

inline xt::xtensor<size_t, 1> VectorPartitionedTyings::iip() const
{
    return m_iip;
}

inline xt::xtensor<size_t, 1> VectorPartitionedTyings::iii() const
{
    return xt::arange<size_t>(m_nni);
}

inline xt::xtensor<size_t, 1> VectorPartitionedTyings::iid() const
{
    return m_iid;
}

template <class T>
inline void VectorPartitionedTyings::copy_p(const T& dofval_src, T& dofval_dest) const
{
    GOOSEFEM_ASSERT(dofval_src.dimension() == 1);
    GOOSEFEM_ASSERT(dofval_dest.dimension() == 1);
    GOOSEFEM_ASSERT(dofval_src.size() == m_ndof || dofval_src.size() == m_nni);
    GOOSEFEM_ASSERT(dofval_dest.size() == m_ndof || dofval_dest.size() == m_nni);

#pragma omp parallel for
    for (size_t i = m_nnu; i < m_nni; ++i) {
        dofval_dest(i) = dofval_src(i);
    }
}

template <class T, class R>
inline void
VectorPartitionedTyings::asDofs_i(const T& nodevec, R& dofval_i, bool apply_tyings) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(dofval_i.size() == m_nni);

    dofval_i.fill(0.0);

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_dofs(m, i) < m_nni) {
                dofval_i(m_dofs(m, i)) = nodevec(m, i);
            }
        }
    }

    if (!apply_tyings) {
        return;
    }

    Eigen::VectorXd Dofval_d = this->Eigen_asDofs_d(nodevec);
    Eigen::VectorXd Dofval_i = m_Cid * Dofval_d;

#pragma omp parallel for
    for (size_t i = 0; i < m_nni; ++i) {
        dofval_i(i) += Dofval_i(i);
    }
}

template <class T>
inline xt::xtensor<double, 1> VectorPartitionedTyings::AsDofs_i(const T& nodevec) const
{
    xt::xtensor<double, 1> dofval = xt::empty<double>({m_nni});
    this->asDofs_i(nodevec, dofval);
    return dofval;
}

template <class T>
inline Eigen::VectorXd VectorPartitionedTyings::Eigen_asDofs_d(const T& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

    Eigen::VectorXd dofval_d(m_nnd, 1);

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_dofs(m, i) >= m_nni) {
                dofval_d(m_dofs(m, i) - m_nni) = nodevec(m, i);
            }
        }
    }

    return dofval_d;
}

} // namespace GooseFEM

#endif

/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_VECTORPARTITIONEDTYINGS_HPP
#define GOOSEFEM_VECTORPARTITIONEDTYINGS_HPP

#include "VectorPartitionedTyings.h"

namespace GooseFEM {

inline VectorPartitionedTyings::VectorPartitionedTyings(
    const xt::xtensor<size_t,2>& conn,
    const xt::xtensor<size_t,2>& dofs,
    const Eigen::SparseMatrix<double>& Cdu,
    const Eigen::SparseMatrix<double>& Cdp,
    const Eigen::SparseMatrix<double>& Cdi)
    : m_conn(conn), m_dofs(dofs), m_Cdu(Cdu), m_Cdp(Cdp), m_Cdi(Cdi)
{
    GOOSEFEM_ASSERT(Cdu.rows() == Cdp.rows());
    GOOSEFEM_ASSERT(Cdi.rows() == Cdp.rows());

    m_nnu = static_cast<size_t>(m_Cdu.cols());
    m_nnp = static_cast<size_t>(m_Cdp.cols());
    m_nnd = static_cast<size_t>(m_Cdp.rows());
    m_nni = m_nnu + m_nnp;
    m_ndof = m_nni + m_nnd;
    m_iiu = xt::arange<size_t>(m_nnu);
    m_iip = xt::arange<size_t>(m_nnu, m_nnu + m_nnp);
    m_iid = xt::arange<size_t>(m_nni, m_nni + m_nnd);
    m_nelem = m_conn.shape(0);
    m_nne = m_conn.shape(1);
    m_nnode = m_dofs.shape(0);
    m_ndim = m_dofs.shape(1);
    m_Cud = m_Cdu.transpose();
    m_Cpd = m_Cdp.transpose();
    m_Cid = m_Cdi.transpose();

    GOOSEFEM_ASSERT(static_cast<size_t>(m_Cdi.cols()) == m_nni);
    GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
    GOOSEFEM_ASSERT(m_ndof == xt::amax(m_dofs)[0] + 1);
}

inline size_t VectorPartitionedTyings::nelem() const
{
    return m_nelem;
}

inline size_t VectorPartitionedTyings::nne() const
{
    return m_nne;
}

inline size_t VectorPartitionedTyings::nnode() const
{
    return m_nnode;
}

inline size_t VectorPartitionedTyings::ndim() const
{
    return m_ndim;
}

inline size_t VectorPartitionedTyings::ndof() const
{
    return m_ndof;
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

inline xt::xtensor<size_t,2> VectorPartitionedTyings::dofs() const
{
    return m_dofs;
}

inline xt::xtensor<size_t,1> VectorPartitionedTyings::iiu() const
{
    return m_iiu;
}

inline xt::xtensor<size_t,1> VectorPartitionedTyings::iip() const
{
    return m_iip;
}

inline xt::xtensor<size_t,1> VectorPartitionedTyings::iii() const
{
    return xt::arange<size_t>(m_nni);
}

inline xt::xtensor<size_t,1> VectorPartitionedTyings::iid() const
{
    return m_iid;
}

inline void VectorPartitionedTyings::copy_p(
    const xt::xtensor<double,1>& dofval_src, xt::xtensor<double,1>& dofval_dest) const
{
    GOOSEFEM_ASSERT(dofval_src.size() == m_ndof || dofval_src.size() == m_nni);
    GOOSEFEM_ASSERT(dofval_dest.size() == m_ndof || dofval_dest.size() == m_nni);

    #pragma omp parallel for
    for (size_t i = m_nnu; i < m_nni; ++i) {
        dofval_dest(i) = dofval_src(i);
    }
}

inline void VectorPartitionedTyings::asDofs_i(
    const xt::xtensor<double,2>& nodevec,
    xt::xtensor<double,1>& dofval_i,
    bool apply_tyings) const
{
    GOOSEFEM_ASSERT(
        nodevec.shape() == std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(dofval_i.size() == m_nni);

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_dofs(m, i) < m_nni) {
                dofval_i(m_dofs(m, i)) = nodevec(m, i);
            }
        }
    }

    if (!apply_tyings)
        return;

    Eigen::VectorXd Dofval_d = this->Eigen_asDofs_d(nodevec);
    Eigen::VectorXd Dofval_i = m_Cid * Dofval_d;

    #pragma omp parallel for
    for (size_t i = 0; i < m_nni; ++i) {
        dofval_i(i) += Dofval_i(i);
    }
}

inline void VectorPartitionedTyings::asNode(
    const xt::xtensor<double,1>& dofval, xt::xtensor<double,2>& nodevec) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);
    GOOSEFEM_ASSERT(
        nodevec.shape() == std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            nodevec(m, i) = dofval(m_dofs(m, i));
        }
    }
}

inline void VectorPartitionedTyings::asElement(
    const xt::xtensor<double,2>& nodevec, xt::xtensor<double,3>& elemvec) const
{
    GOOSEFEM_ASSERT(
        nodevec.shape() == std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(
        elemvec.shape() == std::decay_t<decltype(elemvec)>::shape_type({m_nelem, m_nne, m_ndim}));

    #pragma omp parallel for
    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                elemvec(e, m, i) = nodevec(m_conn(e, m), i);
            }
        }
    }
}

inline void VectorPartitionedTyings::assembleDofs(
    const xt::xtensor<double,3>& elemvec, xt::xtensor<double,1>& dofval) const
{
    GOOSEFEM_ASSERT(
        elemvec.shape() == std::decay_t<decltype(elemvec)>::shape_type({m_nelem, m_nne, m_ndim}));
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

inline void VectorPartitionedTyings::assembleNode(
    const xt::xtensor<double,3>& elemvec, xt::xtensor<double,2>& nodevec) const
{
    GOOSEFEM_ASSERT(
        elemvec.shape() == std::decay_t<decltype(elemvec)>::shape_type({m_nelem, m_nne, m_ndim}));
    GOOSEFEM_ASSERT(
        nodevec.shape() == std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

    xt::xtensor<double,1> dofval = this->AssembleDofs(elemvec);
    this->asNode(dofval, nodevec);
}

inline xt::xtensor<double,1>
VectorPartitionedTyings::AsDofs_i(const xt::xtensor<double,2>& nodevec) const
{
    xt::xtensor<double,1> dofval = xt::empty<double>({m_nni});
    this->asDofs_i(nodevec, dofval);
    return dofval;
}

inline xt::xtensor<double,2>
VectorPartitionedTyings::AsNode(const xt::xtensor<double,1>& dofval) const
{
    xt::xtensor<double,2> nodevec = xt::empty<double>({m_nnode, m_ndim});
    this->asNode(dofval, nodevec);
    return nodevec;
}

inline xt::xtensor<double,3>
VectorPartitionedTyings::AsElement(const xt::xtensor<double,2>& nodevec) const
{
    xt::xtensor<double,3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});
    this->asElement(nodevec, elemvec);
    return elemvec;
}

inline xt::xtensor<double,1>
VectorPartitionedTyings::AssembleDofs(const xt::xtensor<double,3>& elemvec) const
{
    xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});
    this->assembleDofs(elemvec, dofval);
    return dofval;
}

inline xt::xtensor<double,2>
VectorPartitionedTyings::AssembleNode(const xt::xtensor<double,3>& elemvec) const
{
    xt::xtensor<double,2> nodevec = xt::empty<double>({m_nnode, m_ndim});
    this->assembleNode(elemvec, nodevec);
    return nodevec;
}

inline Eigen::VectorXd
VectorPartitionedTyings::Eigen_asDofs_d(const xt::xtensor<double,2>& nodevec) const
{
    GOOSEFEM_ASSERT(
        nodevec.shape() == std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

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

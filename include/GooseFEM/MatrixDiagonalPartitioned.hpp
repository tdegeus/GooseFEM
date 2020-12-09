/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_MATRIXDIAGONALPARTITIONED_HPP
#define GOOSEFEM_MATRIXDIAGONALPARTITIONED_HPP

#include "MatrixDiagonalPartitioned.h"
#include "Mesh.h"

namespace GooseFEM {

inline MatrixDiagonalPartitioned::MatrixDiagonalPartitioned(
    const xt::xtensor<size_t, 2>& conn,
    const xt::xtensor<size_t, 2>& dofs,
    const xt::xtensor<size_t, 1>& iip)
    : m_conn(conn), m_dofs(dofs), m_iip(iip)
{
    m_nelem = m_conn.shape(0);
    m_nne = m_conn.shape(1);
    m_nnode = m_dofs.shape(0);
    m_ndim = m_dofs.shape(1);
    m_iiu = xt::setdiff1d(dofs, iip);
    m_ndof = xt::amax(m_dofs)() + 1;
    m_nnp = m_iip.size();
    m_nnu = m_iiu.size();
    m_part = Mesh::Reorder({m_iiu, m_iip}).get(m_dofs);
    m_Auu = xt::empty<double>({m_nnu});
    m_App = xt::empty<double>({m_nnp});
    m_inv_uu = xt::empty<double>({m_nnu});

    GOOSEFEM_ASSERT(xt::amax(m_conn)() + 1 <= m_nnode);
    GOOSEFEM_ASSERT(xt::amax(m_iip)() <= xt::amax(m_dofs)());
    GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
}

inline size_t MatrixDiagonalPartitioned::nelem() const
{
    return m_nelem;
}

inline size_t MatrixDiagonalPartitioned::nne() const
{
    return m_nne;
}

inline size_t MatrixDiagonalPartitioned::nnode() const
{
    return m_nnode;
}

inline size_t MatrixDiagonalPartitioned::ndim() const
{
    return m_ndim;
}

inline size_t MatrixDiagonalPartitioned::ndof() const
{
    return m_ndof;
}

inline size_t MatrixDiagonalPartitioned::nnu() const
{
    return m_nnu;
}

inline size_t MatrixDiagonalPartitioned::nnp() const
{
    return m_nnp;
}

inline xt::xtensor<size_t, 2> MatrixDiagonalPartitioned::dofs() const
{
    return m_dofs;
}

inline xt::xtensor<size_t, 1> MatrixDiagonalPartitioned::iiu() const
{
    return m_iiu;
}

inline xt::xtensor<size_t, 1> MatrixDiagonalPartitioned::iip() const
{
    return m_iip;
}

inline void MatrixDiagonalPartitioned::factorize()
{
    if (!m_factor) {
        return;
    }

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        m_inv_uu(d) = 1.0 / m_Auu(d);
    }

    m_factor = false;
}

inline void MatrixDiagonalPartitioned::assemble(const xt::xtensor<double, 3>& elemmat)
{
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));
    GOOSEFEM_ASSERT(Element::isDiagonal(elemmat));

    m_Auu.fill(0.0);
    m_App.fill(0.0);

    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {

                size_t d = m_part(m_conn(e, m), i);

                if (d < m_nnu) {
                    m_Auu(d) += elemmat(e, m * m_ndim + i, m * m_ndim + i);
                }
                else {
                    m_App(d - m_nnu) += elemmat(e, m * m_ndim + i, m * m_ndim + i);
                }
            }
        }
    }

    m_factor = true;
}

inline void
MatrixDiagonalPartitioned::dot(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const
{
    GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {

            size_t d = m_part(m, i);

            if (d < m_nnu) {
                b(m, i) = m_Auu(d) * x(m, i);
            }
            else {
                b(m, i) = m_App(d - m_nnu) * x(m, i);
            }
        }
    }
}

inline void
MatrixDiagonalPartitioned::dot(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const
{
    GOOSEFEM_ASSERT(x.size() == m_ndof);
    GOOSEFEM_ASSERT(b.size() == m_ndof);

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        b(m_iiu(d)) = m_Auu(d) * x(m_iiu(d));
    }

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        b(m_iip(d)) = m_App(d) * x(m_iip(d));
    }
}

inline void MatrixDiagonalPartitioned::dot_u(
    const xt::xtensor<double, 1>& x_u,
    const xt::xtensor<double, 1>& x_p,
    xt::xtensor<double, 1>& b_u) const
{
    UNUSED(x_p);

    GOOSEFEM_ASSERT(x_u.size() == m_nnu);
    GOOSEFEM_ASSERT(x_p.size() == m_nnp);
    GOOSEFEM_ASSERT(b_u.size() == m_nnu);

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        b_u(d) = m_Auu(d) * x_u(d);
    }
}

inline void MatrixDiagonalPartitioned::dot_p(
    const xt::xtensor<double, 1>& x_u,
    const xt::xtensor<double, 1>& x_p,
    xt::xtensor<double, 1>& b_p) const
{
    UNUSED(x_u);

    GOOSEFEM_ASSERT(x_u.size() == m_nnu);
    GOOSEFEM_ASSERT(x_p.size() == m_nnp);
    GOOSEFEM_ASSERT(b_p.size() == m_nnp);

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        b_p(d) = m_App(d) * x_p(d);
    }
}

inline void
MatrixDiagonalPartitioned::solve(const xt::xtensor<double, 2>& b, xt::xtensor<double, 2>& x)
{
    GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));

    this->factorize();

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) < m_nnu) {
                x(m, i) = m_inv_uu(m_part(m, i)) * b(m, i);
            }
        }
    }
}

inline void
MatrixDiagonalPartitioned::solve(const xt::xtensor<double, 1>& b, xt::xtensor<double, 1>& x)
{
    GOOSEFEM_ASSERT(b.size() == m_ndof);
    GOOSEFEM_ASSERT(x.size() == m_ndof);

    this->factorize();

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        x(m_iiu(d)) = m_inv_uu(d) * b(m_iiu(d));
    }
}

inline void MatrixDiagonalPartitioned::solve_u(
    const xt::xtensor<double, 1>& b_u,
    const xt::xtensor<double, 1>& x_p,
    xt::xtensor<double, 1>& x_u)
{
    UNUSED(x_p);

    GOOSEFEM_ASSERT(b_u.size() == m_nnu);
    GOOSEFEM_ASSERT(x_p.size() == m_nnp);
    GOOSEFEM_ASSERT(x_u.size() == m_nnu);

    this->factorize();

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        x_u(d) = m_inv_uu(d) * b_u(d);
    }
}

inline void MatrixDiagonalPartitioned::reaction(
    const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const
{
    GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) >= m_nnu) {
                b(m, i) = m_App(m_part(m, i) - m_nnu) * x(m, i);
            }
        }
    }
}

inline void MatrixDiagonalPartitioned::reaction(
    const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const
{
    GOOSEFEM_ASSERT(x.size() == m_ndof);
    GOOSEFEM_ASSERT(b.size() == m_ndof);

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        b(m_iip(d)) = m_App(d) * x(m_iip(d));
    }
}

inline void MatrixDiagonalPartitioned::reaction_p(
    const xt::xtensor<double, 1>& x_u,
    const xt::xtensor<double, 1>& x_p,
    xt::xtensor<double, 1>& b_p) const
{
    UNUSED(x_u);

    GOOSEFEM_ASSERT(x_u.size() == m_nnu);
    GOOSEFEM_ASSERT(x_p.size() == m_nnp);
    GOOSEFEM_ASSERT(b_p.size() == m_nnp);

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        b_p(d) = m_App(d) * x_p(d);
    }
}

inline xt::xtensor<double, 1> MatrixDiagonalPartitioned::Todiagonal() const
{
    xt::xtensor<double, 1> ret = xt::zeros<double>({m_ndof});

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        ret(m_iiu(d)) = m_Auu(d);
    }

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        ret(m_iip(d)) = m_App(d);
    }

    return ret;
}

inline xt::xtensor<double, 2> MatrixDiagonalPartitioned::Dot(const xt::xtensor<double, 2>& x) const
{
    xt::xtensor<double, 2> b = xt::empty<double>({m_nnode, m_ndim});
    this->dot(x, b);
    return b;
}

inline xt::xtensor<double, 1> MatrixDiagonalPartitioned::Dot(const xt::xtensor<double, 1>& x) const
{
    xt::xtensor<double, 1> b = xt::empty<double>({m_ndof});
    this->dot(x, b);
    return b;
}

inline xt::xtensor<double, 1> MatrixDiagonalPartitioned::Dot_u(
    const xt::xtensor<double, 1>& x_u, const xt::xtensor<double, 1>& x_p) const
{
    xt::xtensor<double, 1> b_u = xt::empty<double>({m_nnu});
    this->dot_u(x_u, x_p, b_u);
    return b_u;
}

inline xt::xtensor<double, 1> MatrixDiagonalPartitioned::Dot_p(
    const xt::xtensor<double, 1>& x_u, const xt::xtensor<double, 1>& x_p) const
{
    xt::xtensor<double, 1> b_p = xt::empty<double>({m_nnp});
    this->dot_p(x_u, x_p, b_p);
    return b_p;
}

inline xt::xtensor<double, 2>
MatrixDiagonalPartitioned::Solve(const xt::xtensor<double, 2>& b, const xt::xtensor<double, 2>& x)
{
    xt::xtensor<double, 2> ret = x;
    this->solve(b, ret);
    return ret;
}

inline xt::xtensor<double, 1>
MatrixDiagonalPartitioned::Solve(const xt::xtensor<double, 1>& b, const xt::xtensor<double, 1>& x)
{
    xt::xtensor<double, 1> ret = x;
    this->solve(b, ret);
    return ret;
}

inline xt::xtensor<double, 1> MatrixDiagonalPartitioned::Solve_u(
    const xt::xtensor<double, 1>& b_u, const xt::xtensor<double, 1>& x_p)
{
    xt::xtensor<double, 1> x_u = xt::empty<double>({m_nnu});
    this->solve_u(b_u, x_p, x_u);
    return x_u;
}

inline xt::xtensor<double, 2> MatrixDiagonalPartitioned::Reaction(
    const xt::xtensor<double, 2>& x, const xt::xtensor<double, 2>& b) const
{
    xt::xtensor<double, 2> ret = b;
    this->reaction(x, ret);
    return ret;
}

inline xt::xtensor<double, 1> MatrixDiagonalPartitioned::Reaction(
    const xt::xtensor<double, 1>& x, const xt::xtensor<double, 1>& b) const
{
    xt::xtensor<double, 1> ret = b;
    this->reaction(x, ret);
    return ret;
}

inline xt::xtensor<double, 1> MatrixDiagonalPartitioned::Reaction_p(
    const xt::xtensor<double, 1>& x_u, const xt::xtensor<double, 1>& x_p) const
{
    xt::xtensor<double, 1> b_p = xt::empty<double>({m_nnp});
    this->reaction_p(x_u, x_p, b_p);
    return b_p;
}

} // namespace GooseFEM

#endif

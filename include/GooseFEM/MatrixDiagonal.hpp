/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_MATRIXDIAGONAL_HPP
#define GOOSEFEM_MATRIXDIAGONAL_HPP

#include "MatrixDiagonal.h"

namespace GooseFEM {

inline MatrixDiagonal::MatrixDiagonal(
    const xt::xtensor<size_t, 2>& conn, const xt::xtensor<size_t, 2>& dofs)
    : m_conn(conn), m_dofs(dofs)
{
    m_nelem = m_conn.shape(0);
    m_nne = m_conn.shape(1);
    m_nnode = m_dofs.shape(0);
    m_ndim = m_dofs.shape(1);
    m_ndof = xt::amax(m_dofs)[0] + 1;
    m_A = xt::empty<double>({m_ndof});
    m_inv = xt::empty<double>({m_ndof});

    GOOSEFEM_ASSERT(xt::amax(m_conn)[0] + 1 == m_nnode);
    GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
}

inline size_t MatrixDiagonal::nelem() const
{
    return m_nelem;
}

inline size_t MatrixDiagonal::nne() const
{
    return m_nne;
}

inline size_t MatrixDiagonal::nnode() const
{
    return m_nnode;
}

inline size_t MatrixDiagonal::ndim() const
{
    return m_ndim;
}

inline size_t MatrixDiagonal::ndof() const
{
    return m_ndof;
}

inline xt::xtensor<size_t, 2> MatrixDiagonal::dofs() const
{
    return m_dofs;
}

inline void MatrixDiagonal::factorize()
{
    if (!m_factor) {
        return;
    }

    #pragma omp parallel for
    for (size_t d = 0; d < m_ndof; ++d)
        m_inv(d) = 1.0 / m_A(d);

    m_factor = false;
}

inline void MatrixDiagonal::set(const xt::xtensor<double, 1>& A)
{
    GOOSEFEM_ASSERT(A.size() == m_ndof);
    std::copy(A.begin(), A.end(), m_A.begin());
    m_factor = true;
}

inline void MatrixDiagonal::assemble(const xt::xtensor<double, 3>& elemmat)
{
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));
    GOOSEFEM_ASSERT(Element::isDiagonal(elemmat));

    m_A.fill(0.0);

    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                m_A(m_dofs(m_conn(e, m), i)) += elemmat(e, m * m_ndim + i, m * m_ndim + i);
            }
        }
    }

    m_factor = true;
}

inline void MatrixDiagonal::dot(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const
{
    GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            b(m, i) = m_A(m_dofs(m, i)) * x(m, i);
        }
    }
}

inline void MatrixDiagonal::dot(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const
{
    GOOSEFEM_ASSERT(x.size() == m_ndof);
    GOOSEFEM_ASSERT(b.size() == m_ndof);

    xt::noalias(b) = m_A * x;
}

inline void MatrixDiagonal::solve(const xt::xtensor<double, 2>& b, xt::xtensor<double, 2>& x)
{
    GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));

    this->factorize();

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            x(m, i) = m_inv(m_dofs(m, i)) * b(m, i);
        }
    }
}

inline void MatrixDiagonal::solve(const xt::xtensor<double, 1>& b, xt::xtensor<double, 1>& x)
{
    GOOSEFEM_ASSERT(b.size() == m_ndof);
    GOOSEFEM_ASSERT(x.size() == m_ndof);
    this->factorize();
    xt::noalias(x) = m_inv * b;
}

inline xt::xtensor<double, 1> MatrixDiagonal::AsDiagonal() const
{
    return m_A;
}

inline xt::xtensor<double, 2> MatrixDiagonal::Dot(const xt::xtensor<double, 2>& x) const
{
    xt::xtensor<double, 2> b = xt::empty<double>({m_nnode, m_ndim});
    this->dot(x, b);
    return b;
}

inline xt::xtensor<double, 1> MatrixDiagonal::Dot(const xt::xtensor<double, 1>& x) const
{
    xt::xtensor<double, 1> b = xt::empty<double>({m_ndof});
    this->dot(x, b);
    return b;
}

inline xt::xtensor<double, 2> MatrixDiagonal::Solve(const xt::xtensor<double, 2>& b)
{
    xt::xtensor<double, 2> x = xt::empty<double>({m_nnode, m_ndim});
    this->solve(b, x);
    return x;
}

inline xt::xtensor<double, 1> MatrixDiagonal::Solve(const xt::xtensor<double, 1>& b)
{
    xt::xtensor<double, 1> x = xt::empty<double>({m_ndof});
    this->solve(b, x);
    return x;
}

} // namespace GooseFEM

#endif

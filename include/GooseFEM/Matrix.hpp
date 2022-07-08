/**
Implementation of Matrix.h

\file Matrix.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIX_HPP
#define GOOSEFEM_MATRIX_HPP

#include "Matrix.h"

namespace GooseFEM {

inline Matrix::Matrix(
    const array_type::tensor<size_t, 2>& conn,
    const array_type::tensor<size_t, 2>& dofs)
{
    m_conn = conn;
    m_dofs = dofs;
    m_nelem = m_conn.shape(0);
    m_nne = m_conn.shape(1);
    m_nnode = m_dofs.shape(0);
    m_ndim = m_dofs.shape(1);
    m_ndof = xt::amax(m_dofs)() + 1;
    m_T.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_A.resize(m_ndof, m_ndof);

    GOOSEFEM_ASSERT(xt::amax(m_conn)() + 1 <= m_nnode);
    GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
}

inline size_t Matrix::nelem() const
{
    return m_nelem;
}

inline size_t Matrix::nne() const
{
    return m_nne;
}

inline size_t Matrix::nnode() const
{
    return m_nnode;
}

inline size_t Matrix::ndim() const
{
    return m_ndim;
}

inline size_t Matrix::ndof() const
{
    return m_ndof;
}

inline array_type::tensor<size_t, 2> Matrix::dofs() const
{
    return m_dofs;
}

inline void Matrix::assemble(const array_type::tensor<double, 3>& elemmat)
{
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));

    m_T.clear();

    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                for (size_t n = 0; n < m_nne; ++n) {
                    for (size_t j = 0; j < m_ndim; ++j) {
                        m_T.push_back(Eigen::Triplet<double>(
                            m_dofs(m_conn(e, m), i),
                            m_dofs(m_conn(e, n), j),
                            elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                    }
                }
            }
        }
    }

    m_A.setFromTriplets(m_T.begin(), m_T.end());
    m_changed = true;
}

inline void Matrix::set(
    const array_type::tensor<size_t, 1>& rows,
    const array_type::tensor<size_t, 1>& cols,
    const array_type::tensor<double, 2>& matrix)
{
    GOOSEFEM_ASSERT(rows.size() == matrix.shape(0));
    GOOSEFEM_ASSERT(cols.size() == matrix.shape(1));
    GOOSEFEM_ASSERT(xt::amax(cols)() < m_ndof);
    GOOSEFEM_ASSERT(xt::amax(rows)() < m_ndof);

    std::vector<Eigen::Triplet<double>> T;

    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = 0; j < cols.size(); ++j) {
            T.push_back(Eigen::Triplet<double>(rows(i), cols(j), matrix(i, j)));
        }
    }

    m_A.setFromTriplets(T.begin(), T.end());
    m_changed = true;
}

inline void Matrix::add(
    const array_type::tensor<size_t, 1>& rows,
    const array_type::tensor<size_t, 1>& cols,
    const array_type::tensor<double, 2>& matrix)
{
    GOOSEFEM_ASSERT(rows.size() == matrix.shape(0));
    GOOSEFEM_ASSERT(cols.size() == matrix.shape(1));
    GOOSEFEM_ASSERT(xt::amax(cols)() < m_ndof);
    GOOSEFEM_ASSERT(xt::amax(rows)() < m_ndof);

    std::vector<Eigen::Triplet<double>> T;

    Eigen::SparseMatrix<double> A(m_ndof, m_ndof);

    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = 0; j < cols.size(); ++j) {
            T.push_back(Eigen::Triplet<double>(rows(i), cols(j), matrix(i, j)));
        }
    }

    A.setFromTriplets(T.begin(), T.end());
    m_A += A;
    m_changed = true;
}

inline void Matrix::todense(array_type::tensor<double, 2>& ret) const
{
    GOOSEFEM_ASSERT(xt::has_shape(ret, {m_ndof, m_ndof}));

    ret.fill(0.0);

    for (int k = 0; k < m_A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_A, k); it; ++it) {
            ret(it.row(), it.col()) = it.value();
        }
    }
}

inline array_type::tensor<double, 2> Matrix::Todense() const
{
    array_type::tensor<double, 2> ret = xt::empty<double>({m_ndof, m_ndof});
    this->todense(ret);
    return ret;
}

inline void
Matrix::dot(const array_type::tensor<double, 2>& x, array_type::tensor<double, 2>& b) const
{
    GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
    this->asNode(m_A * this->AsDofs(x), b);
}

inline void
Matrix::dot(const array_type::tensor<double, 1>& x, array_type::tensor<double, 1>& b) const
{
    GOOSEFEM_ASSERT(b.size() == m_ndof);
    GOOSEFEM_ASSERT(x.size() == m_ndof);

    Eigen::Map<Eigen::VectorXd>(b.data(), b.size()).noalias() =
        m_A * Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
}

inline array_type::tensor<double, 2> Matrix::Dot(const array_type::tensor<double, 2>& x) const
{
    array_type::tensor<double, 2> b = xt::empty<double>({m_nnode, m_ndim});
    this->dot(x, b);
    return b;
}

inline array_type::tensor<double, 1> Matrix::Dot(const array_type::tensor<double, 1>& x) const
{
    array_type::tensor<double, 1> b = xt::empty<double>({m_ndof});
    this->dot(x, b);
    return b;
}

inline Eigen::VectorXd Matrix::AsDofs(const array_type::tensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

    Eigen::VectorXd dofval = Eigen::VectorXd::Zero(m_ndof, 1);

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            dofval(m_dofs(m, i)) = nodevec(m, i);
        }
    }

    return dofval;
}

inline void
Matrix::asNode(const Eigen::VectorXd& dofval, array_type::tensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(static_cast<size_t>(dofval.size()) == m_ndof);
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            nodevec(m, i) = dofval(m_dofs(m, i));
        }
    }
}

template <class Solver>
inline void MatrixSolver<Solver>::factorize(Matrix& matrix)
{
    if (!matrix.m_changed && !m_factor) {
        return;
    }
    m_solver.compute(matrix.m_A);
    m_factor = false;
    matrix.m_changed = false;
}

template <class Solver>
inline void MatrixSolver<Solver>::solve(
    Matrix& matrix,
    const array_type::tensor<double, 2>& b,
    array_type::tensor<double, 2>& x)
{
    GOOSEFEM_ASSERT(xt::has_shape(b, {matrix.m_nnode, matrix.m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(x, {matrix.m_nnode, matrix.m_ndim}));
    this->factorize(matrix);
    Eigen::VectorXd X = m_solver.solve(matrix.AsDofs(b));
    matrix.asNode(X, x);
}

template <class Solver>
inline void MatrixSolver<Solver>::solve(
    Matrix& matrix,
    const array_type::tensor<double, 1>& b,
    array_type::tensor<double, 1>& x)
{
    GOOSEFEM_ASSERT(b.size() == matrix.m_ndof);
    GOOSEFEM_ASSERT(x.size() == matrix.m_ndof);
    this->factorize(matrix);
    Eigen::Map<Eigen::VectorXd>(x.data(), x.size()).noalias() =
        m_solver.solve(Eigen::Map<const Eigen::VectorXd>(b.data(), matrix.m_ndof));
}

template <class Solver>
inline array_type::tensor<double, 2>
MatrixSolver<Solver>::Solve(Matrix& matrix, const array_type::tensor<double, 2>& b)
{
    array_type::tensor<double, 2> x = xt::empty<double>({matrix.m_nnode, matrix.m_ndim});
    this->solve(matrix, b, x);
    return x;
}

template <class Solver>
inline array_type::tensor<double, 1>
MatrixSolver<Solver>::Solve(Matrix& matrix, const array_type::tensor<double, 1>& b)
{
    array_type::tensor<double, 1> x = xt::empty<double>({matrix.m_ndof});
    this->solve(matrix, b, x);
    return x;
}

} // namespace GooseFEM

#endif

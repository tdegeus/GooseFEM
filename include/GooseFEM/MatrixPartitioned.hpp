/**
Implementation of MatrixPartitioned.h

\file MatrixPartitioned.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXPARTITIONED_HPP
#define GOOSEFEM_MATRIXPARTITIONED_HPP

#include "MatrixPartitioned.h"
#include "Mesh.h"

namespace GooseFEM {

inline MatrixPartitioned::MatrixPartitioned(
    const xt::xtensor<size_t, 2>& conn,
    const xt::xtensor<size_t, 2>& dofs,
    const xt::xtensor<size_t, 1>& iip)
{
    m_conn = conn;
    m_dofs = dofs;
    m_iip = iip;
    m_nelem = m_conn.shape(0);
    m_nne = m_conn.shape(1);
    m_nnode = m_dofs.shape(0);
    m_ndim = m_dofs.shape(1);
    m_ndof = xt::amax(m_dofs)() + 1;

    GOOSEFEM_ASSERT(is_unique(iip));
    GOOSEFEM_ASSERT(xt::amax(conn)() + 1 <= m_nnode);
    GOOSEFEM_ASSERT(xt::amax(iip)() <= xt::amax(dofs)());
    GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);

    m_iiu = xt::setdiff1d(dofs, iip);
    m_nnp = m_iip.size();
    m_nnu = m_iiu.size();

    m_part = Mesh::Reorder({m_iiu, m_iip}).apply(m_dofs);
    m_Tuu.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tup.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tpu.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tpp.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Auu.resize(m_nnu, m_nnu);
    m_Aup.resize(m_nnu, m_nnp);
    m_Apu.resize(m_nnp, m_nnu);
    m_App.resize(m_nnp, m_nnp);
}

inline size_t MatrixPartitioned::nnu() const
{
    return m_nnu;
}

inline size_t MatrixPartitioned::nnp() const
{
    return m_nnp;
}

inline xt::xtensor<size_t, 1> MatrixPartitioned::iiu() const
{
    return m_iiu;
}

inline xt::xtensor<size_t, 1> MatrixPartitioned::iip() const
{
    return m_iip;
}

inline void MatrixPartitioned::assemble(const xt::xtensor<double, 3>& elemmat)
{
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));

    m_Tuu.clear();
    m_Tup.clear();
    m_Tpu.clear();
    m_Tpp.clear();

    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {

                size_t di = m_part(m_conn(e, m), i);

                for (size_t n = 0; n < m_nne; ++n) {
                    for (size_t j = 0; j < m_ndim; ++j) {

                        size_t dj = m_part(m_conn(e, n), j);

                        if (di < m_nnu && dj < m_nnu) {
                            m_Tuu.push_back(Eigen::Triplet<double>(
                                di, dj, elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (di < m_nnu) {
                            m_Tup.push_back(Eigen::Triplet<double>(
                                di, dj - m_nnu, elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (dj < m_nnu) {
                            m_Tpu.push_back(Eigen::Triplet<double>(
                                di - m_nnu, dj, elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else {
                            m_Tpp.push_back(Eigen::Triplet<double>(
                                di - m_nnu,
                                dj - m_nnu,
                                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                    }
                }
            }
        }
    }

    m_Auu.setFromTriplets(m_Tuu.begin(), m_Tuu.end());
    m_Aup.setFromTriplets(m_Tup.begin(), m_Tup.end());
    m_Apu.setFromTriplets(m_Tpu.begin(), m_Tpu.end());
    m_App.setFromTriplets(m_Tpp.begin(), m_Tpp.end());
    m_changed = true;
}

inline void MatrixPartitioned::set(
    const xt::xtensor<size_t, 1>& rows,
    const xt::xtensor<size_t, 1>& cols,
    const xt::xtensor<double, 2>& matrix)
{
    GOOSEFEM_ASSERT(rows.size() == matrix.shape(0));
    GOOSEFEM_ASSERT(cols.size() == matrix.shape(1));
    GOOSEFEM_ASSERT(xt::amax(cols)() < m_ndof);
    GOOSEFEM_ASSERT(xt::amax(rows)() < m_ndof);

    std::vector<Eigen::Triplet<double>> Tuu;
    std::vector<Eigen::Triplet<double>> Tup;
    std::vector<Eigen::Triplet<double>> Tpu;
    std::vector<Eigen::Triplet<double>> Tpp;

    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = 0; j < cols.size(); ++j) {
            size_t di = rows(i);
            size_t dj = cols(j);
            double v = matrix(i, j);

            if (di < m_nnu && dj < m_nnu) {
                Tuu.push_back(Eigen::Triplet<double>(di, dj, v));
            }
            else if (di < m_nnu) {
                Tup.push_back(Eigen::Triplet<double>(di, dj - m_nnu, v));
            }
            else if (dj < m_nnu) {
                Tpu.push_back(Eigen::Triplet<double>(di - m_nnu, dj, v));
            }
            else {
                Tpp.push_back(Eigen::Triplet<double>(di - m_nnu, dj - m_nnu, v));
            }
        }
    }

    m_Auu.setFromTriplets(Tuu.begin(), Tuu.end());
    m_Aup.setFromTriplets(Tup.begin(), Tup.end());
    m_Apu.setFromTriplets(Tpu.begin(), Tpu.end());
    m_App.setFromTriplets(Tpp.begin(), Tpp.end());
    m_changed = true;
}

inline void MatrixPartitioned::add(
    const xt::xtensor<size_t, 1>& rows,
    const xt::xtensor<size_t, 1>& cols,
    const xt::xtensor<double, 2>& matrix)
{
    GOOSEFEM_ASSERT(rows.size() == matrix.shape(0));
    GOOSEFEM_ASSERT(cols.size() == matrix.shape(1));
    GOOSEFEM_ASSERT(xt::amax(cols)() < m_ndof);
    GOOSEFEM_ASSERT(xt::amax(rows)() < m_ndof);

    std::vector<Eigen::Triplet<double>> Tuu;
    std::vector<Eigen::Triplet<double>> Tup;
    std::vector<Eigen::Triplet<double>> Tpu;
    std::vector<Eigen::Triplet<double>> Tpp;

    Eigen::SparseMatrix<double> Auu(m_nnu, m_nnu);
    Eigen::SparseMatrix<double> Aup(m_nnu, m_nnp);
    Eigen::SparseMatrix<double> Apu(m_nnp, m_nnu);
    Eigen::SparseMatrix<double> App(m_nnp, m_nnp);

    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = 0; j < cols.size(); ++j) {
            size_t di = rows(i);
            size_t dj = cols(j);
            double v = matrix(i, j);

            if (di < m_nnu && dj < m_nnu) {
                Tuu.push_back(Eigen::Triplet<double>(di, dj, v));
            }
            else if (di < m_nnu) {
                Tup.push_back(Eigen::Triplet<double>(di, dj - m_nnu, v));
            }
            else if (dj < m_nnu) {
                Tpu.push_back(Eigen::Triplet<double>(di - m_nnu, dj, v));
            }
            else {
                Tpp.push_back(Eigen::Triplet<double>(di - m_nnu, dj - m_nnu, v));
            }
        }
    }

    Auu.setFromTriplets(Tuu.begin(), Tuu.end());
    Aup.setFromTriplets(Tup.begin(), Tup.end());
    Apu.setFromTriplets(Tpu.begin(), Tpu.end());
    App.setFromTriplets(Tpp.begin(), Tpp.end());
    m_Auu += Auu;
    m_Aup += Aup;
    m_Apu += Apu;
    m_App += App;
    m_changed = true;
}

inline void MatrixPartitioned::todense(xt::xtensor<double, 2>& ret) const
{
    GOOSEFEM_ASSERT(xt::has_shape(ret, {m_ndof, m_ndof}));

    ret.fill(0.0);

    for (int k = 0; k < m_Auu.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_Auu, k); it; ++it) {
            ret(it.row(), it.col()) = it.value();
        }
    }

    for (int k = 0; k < m_Aup.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_Aup, k); it; ++it) {
            ret(it.row(), it.col() + m_nnu) = it.value();
        }
    }

    for (int k = 0; k < m_Apu.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_Apu, k); it; ++it) {
            ret(it.row() + m_nnu, it.col()) = it.value();
        }
    }

    for (int k = 0; k < m_App.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_App, k); it; ++it) {
            ret(it.row() + m_nnu, it.col() + m_nnu) = it.value();
        }
    }
}

inline void MatrixPartitioned::dot(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const
{
    GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));

    Eigen::VectorXd X_u = this->AsDofs_u(x);
    Eigen::VectorXd X_p = this->AsDofs_p(x);
    Eigen::VectorXd B_u = m_Auu * X_u + m_Aup * X_p;
    Eigen::VectorXd B_p = m_Apu * X_u + m_App * X_p;

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) < m_nnu) {
                b(m, i) = B_u(m_part(m, i));
            }
            else{
                b(m, i) = B_p(m_part(m, i) - m_nnu);
            }
        }
    }
}

inline void MatrixPartitioned::dot(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const
{
    GOOSEFEM_ASSERT(b.size() == m_ndof);
    GOOSEFEM_ASSERT(x.size() == m_ndof);

    Eigen::VectorXd X_u = this->AsDofs_u(x);
    Eigen::VectorXd X_p = this->AsDofs_p(x);

    Eigen::Map<Eigen::VectorXd>(b.data(), m_nnu).noalias() = m_Auu * X_u + m_Aup * X_p;
    Eigen::Map<Eigen::VectorXd>(b.data() + m_nnu, m_ndof).noalias() = m_Auu * X_u + m_Aup * X_p;
}

inline void
MatrixPartitioned::reaction(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const
{
    GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));

    Eigen::VectorXd X_u = this->AsDofs_u(x);
    Eigen::VectorXd X_p = this->AsDofs_p(x);
    Eigen::VectorXd B_p = m_Apu * X_u + m_App * X_p;

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) >= m_nnu) {
                b(m, i) = B_p(m_part(m, i) - m_nnu);
            }
        }
    }
}

inline void
MatrixPartitioned::reaction(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const
{
    GOOSEFEM_ASSERT(x.size() == m_ndof);
    GOOSEFEM_ASSERT(b.size() == m_ndof);

    Eigen::VectorXd X_u = this->AsDofs_u(x);
    Eigen::VectorXd X_p = this->AsDofs_p(x);
    Eigen::VectorXd B_p = m_Apu * X_u + m_App * X_p;

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        b(m_iip(d)) = B_p(d);
    }
}

inline void MatrixPartitioned::reaction_p(
    const xt::xtensor<double, 1>& x_u,
    const xt::xtensor<double, 1>& x_p,
    xt::xtensor<double, 1>& b_p) const
{
    GOOSEFEM_ASSERT(x_u.size() == m_nnu);
    GOOSEFEM_ASSERT(x_p.size() == m_nnp);
    GOOSEFEM_ASSERT(b_p.size() == m_nnp);

    Eigen::Map<Eigen::VectorXd>(b_p.data(), b_p.size()).noalias() =
        m_Apu * Eigen::Map<const Eigen::VectorXd>(x_u.data(), x_u.size()) +
        m_App * Eigen::Map<const Eigen::VectorXd>(x_p.data(), x_p.size());
}

inline xt::xtensor<double, 2>
MatrixPartitioned::Reaction(const xt::xtensor<double, 2>& x, const xt::xtensor<double, 2>& b) const
{
    xt::xtensor<double, 2> ret = b;
    this->reaction(x, ret);
    return ret;
}

inline xt::xtensor<double, 1>
MatrixPartitioned::Reaction(const xt::xtensor<double, 1>& x, const xt::xtensor<double, 1>& b) const
{
    xt::xtensor<double, 1> ret = b;
    this->reaction(x, ret);
    return ret;
}

inline xt::xtensor<double, 1> MatrixPartitioned::Reaction_p(
    const xt::xtensor<double, 1>& x_u, const xt::xtensor<double, 1>& x_p) const
{
    xt::xtensor<double, 1> b_p = xt::empty<double>({m_nnp});
    this->reaction_p(x_u, x_p, b_p);
    return b_p;
}

inline Eigen::VectorXd MatrixPartitioned::AsDofs_u(const xt::xtensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    Eigen::VectorXd dofval_u(m_nnu, 1);

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        dofval_u(d) = dofval(m_iiu(d));
    }

    return dofval_u;
}

inline Eigen::VectorXd MatrixPartitioned::AsDofs_u(const xt::xtensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

    Eigen::VectorXd dofval_u = Eigen::VectorXd::Zero(m_nnu, 1);

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) < m_nnu) {
                dofval_u(m_part(m, i)) = nodevec(m, i);
            }
        }
    }

    return dofval_u;
}

inline Eigen::VectorXd MatrixPartitioned::AsDofs_p(const xt::xtensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    Eigen::VectorXd dofval_p(m_nnp, 1);

    #pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        dofval_p(d) = dofval(m_iip(d));
    }

    return dofval_p;
}

inline Eigen::VectorXd MatrixPartitioned::AsDofs_p(const xt::xtensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

    Eigen::VectorXd dofval_p = Eigen::VectorXd::Zero(m_nnp, 1);

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_part(m, i) >= m_nnu) {
                dofval_p(m_part(m, i) - m_nnu) = nodevec(m, i);
            }
        }
    }

    return dofval_p;
}

template <class Solver>
inline void MatrixPartitionedSolver<Solver>::factorize(MatrixPartitioned& matrix)
{
    if (!matrix.m_changed && !m_factor) {
        return;
    }
    m_solver.compute(matrix.m_Auu);
    m_factor = false;
    matrix.m_changed = false;
}

template <class Solver>
inline void MatrixPartitionedSolver<Solver>::solve(
    MatrixPartitioned& matrix, const xt::xtensor<double, 2>& b, xt::xtensor<double, 2>& x)
{
    GOOSEFEM_ASSERT(xt::has_shape(b, {matrix.m_nnode, matrix.m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(x, {matrix.m_nnode, matrix.m_ndim}));

    this->factorize(matrix);
    Eigen::VectorXd B_u = matrix.AsDofs_u(b);
    Eigen::VectorXd X_p = matrix.AsDofs_p(x);
    Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - matrix.m_Aup * X_p));

    #pragma omp parallel for
    for (size_t m = 0; m < matrix.m_nnode; ++m) {
        for (size_t i = 0; i < matrix.m_ndim; ++i) {
            if (matrix.m_part(m, i) < matrix.m_nnu) {
                x(m, i) = X_u(matrix.m_part(m, i));
            }
        }
    }
}

template <class Solver>
inline void MatrixPartitionedSolver<Solver>::solve(
    MatrixPartitioned& matrix, const xt::xtensor<double, 1>& b, xt::xtensor<double, 1>& x)
{
    GOOSEFEM_ASSERT(b.size() == matrix.m_ndof);
    GOOSEFEM_ASSERT(x.size() == matrix.m_ndof);

    this->factorize(matrix);
    Eigen::VectorXd B_u = matrix.AsDofs_u(b);
    Eigen::VectorXd X_p = matrix.AsDofs_p(x);
    Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - matrix.m_Aup * X_p));

    #pragma omp parallel for
    for (size_t d = 0; d < matrix.m_nnu; ++d) {
        x(matrix.m_iiu(d)) = X_u(d);
    }
}

template <class Solver>
inline void MatrixPartitionedSolver<Solver>::solve_u(
    MatrixPartitioned& matrix,
    const xt::xtensor<double, 1>& b_u,
    const xt::xtensor<double, 1>& x_p,
    xt::xtensor<double, 1>& x_u)
{
    GOOSEFEM_ASSERT(b_u.size() == matrix.m_nnu);
    GOOSEFEM_ASSERT(x_p.size() == matrix.m_nnp);
    GOOSEFEM_ASSERT(x_u.size() == matrix.m_nnu);

    this->factorize(matrix);

    Eigen::Map<Eigen::VectorXd>(x_u.data(), x_u.size()).noalias() = m_solver.solve(Eigen::VectorXd(
        Eigen::Map<const Eigen::VectorXd>(b_u.data(), b_u.size()) -
        matrix.m_Aup * Eigen::Map<const Eigen::VectorXd>(x_p.data(), x_p.size())));
}

template <class Solver>
inline xt::xtensor<double, 2> MatrixPartitionedSolver<Solver>::Solve(
    MatrixPartitioned& matrix, const xt::xtensor<double, 2>& b, const xt::xtensor<double, 2>& x)
{
    xt::xtensor<double, 2> ret = x;
    this->solve(matrix, b, ret);
    return ret;
}

template <class Solver>
inline xt::xtensor<double, 1> MatrixPartitionedSolver<Solver>::Solve(
    MatrixPartitioned& matrix, const xt::xtensor<double, 1>& b, const xt::xtensor<double, 1>& x)
{
    xt::xtensor<double, 1> ret = x;
    this->solve(matrix, b, ret);
    return ret;
}

template <class Solver>
inline xt::xtensor<double, 1> MatrixPartitionedSolver<Solver>::Solve_u(
    MatrixPartitioned& matrix, const xt::xtensor<double, 1>& b_u, const xt::xtensor<double, 1>& x_p)
{
    xt::xtensor<double, 1> x_u = xt::empty<double>({matrix.m_nnu});
    this->solve_u(matrix, b_u, x_p, x_u);
    return x_u;
}

} // namespace GooseFEM

#endif

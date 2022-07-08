/**
Implementation of MatrixPartitionedTyings.h

\file MatrixPartitionedTyings.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXPARTITIONEDTYINGS_HPP
#define GOOSEFEM_MATRIXPARTITIONEDTYINGS_HPP

#include "MatrixPartitionedTyings.h"

namespace GooseFEM {

inline MatrixPartitionedTyings::MatrixPartitionedTyings(
    const array_type::tensor<size_t, 2>& conn,
    const array_type::tensor<size_t, 2>& dofs,
    const Eigen::SparseMatrix<double>& Cdu,
    const Eigen::SparseMatrix<double>& Cdp)
{
    GOOSEFEM_ASSERT(Cdu.rows() == Cdp.rows());

    m_conn = conn;
    m_dofs = dofs;
    m_Cdu = Cdu;
    m_Cdp = Cdp;
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
    m_Tuu.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tup.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tpu.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tpp.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tud.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tpd.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tdu.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tdp.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Tdd.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_Auu.resize(m_nnu, m_nnu);
    m_Aup.resize(m_nnu, m_nnp);
    m_Apu.resize(m_nnp, m_nnu);
    m_App.resize(m_nnp, m_nnp);
    m_Aud.resize(m_nnu, m_nnd);
    m_Apd.resize(m_nnp, m_nnd);
    m_Adu.resize(m_nnd, m_nnu);
    m_Adp.resize(m_nnd, m_nnp);
    m_Add.resize(m_nnd, m_nnd);

    GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
    GOOSEFEM_ASSERT(m_ndof == xt::amax(m_dofs)() + 1);
}

inline size_t MatrixPartitionedTyings::nnu() const
{
    return m_nnu;
}

inline size_t MatrixPartitionedTyings::nnp() const
{
    return m_nnp;
}

inline size_t MatrixPartitionedTyings::nni() const
{
    return m_nni;
}

inline size_t MatrixPartitionedTyings::nnd() const
{
    return m_nnd;
}

inline array_type::tensor<size_t, 1> MatrixPartitionedTyings::iiu() const
{
    return m_iiu;
}

inline array_type::tensor<size_t, 1> MatrixPartitionedTyings::iip() const
{
    return m_iip;
}

inline array_type::tensor<size_t, 1> MatrixPartitionedTyings::iii() const
{
    return xt::arange<size_t>(m_nni);
}

inline array_type::tensor<size_t, 1> MatrixPartitionedTyings::iid() const
{
    return m_iid;
}

inline void MatrixPartitionedTyings::assemble(const array_type::tensor<double, 3>& elemmat)
{
    GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));

    m_Tuu.clear();
    m_Tup.clear();
    m_Tpu.clear();
    m_Tpp.clear();
    m_Tud.clear();
    m_Tpd.clear();
    m_Tdu.clear();
    m_Tdp.clear();
    m_Tdd.clear();

    for (size_t e = 0; e < m_nelem; ++e) {
        for (size_t m = 0; m < m_nne; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {

                size_t di = m_dofs(m_conn(e, m), i);

                for (size_t n = 0; n < m_nne; ++n) {
                    for (size_t j = 0; j < m_ndim; ++j) {

                        size_t dj = m_dofs(m_conn(e, n), j);

                        if (di < m_nnu && dj < m_nnu) {
                            m_Tuu.push_back(Eigen::Triplet<double>(
                                di, dj, elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (di < m_nnu && dj < m_nni) {
                            m_Tup.push_back(Eigen::Triplet<double>(
                                di, dj - m_nnu, elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (di < m_nnu) {
                            m_Tud.push_back(Eigen::Triplet<double>(
                                di, dj - m_nni, elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (di < m_nni && dj < m_nnu) {
                            m_Tpu.push_back(Eigen::Triplet<double>(
                                di - m_nnu, dj, elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (di < m_nni && dj < m_nni) {
                            m_Tpp.push_back(Eigen::Triplet<double>(
                                di - m_nnu,
                                dj - m_nnu,
                                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (di < m_nni) {
                            m_Tpd.push_back(Eigen::Triplet<double>(
                                di - m_nnu,
                                dj - m_nni,
                                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (dj < m_nnu) {
                            m_Tdu.push_back(Eigen::Triplet<double>(
                                di - m_nni, dj, elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else if (dj < m_nni) {
                            m_Tdp.push_back(Eigen::Triplet<double>(
                                di - m_nni,
                                dj - m_nnu,
                                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                        else {
                            m_Tdd.push_back(Eigen::Triplet<double>(
                                di - m_nni,
                                dj - m_nni,
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
    m_Aud.setFromTriplets(m_Tud.begin(), m_Tud.end());
    m_Apd.setFromTriplets(m_Tpd.begin(), m_Tpd.end());
    m_Adu.setFromTriplets(m_Tdu.begin(), m_Tdu.end());
    m_Adp.setFromTriplets(m_Tdp.begin(), m_Tdp.end());
    m_Add.setFromTriplets(m_Tdd.begin(), m_Tdd.end());
    m_changed = true;
}

inline Eigen::VectorXd
MatrixPartitionedTyings::AsDofs_u(const array_type::tensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    Eigen::VectorXd dofval_u(m_nnu, 1);

#pragma omp parallel for
    for (size_t d = 0; d < m_nnu; ++d) {
        dofval_u(d) = dofval(m_iiu(d));
    }

    return dofval_u;
}

inline Eigen::VectorXd
MatrixPartitionedTyings::AsDofs_u(const array_type::tensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

    Eigen::VectorXd dofval_u = Eigen::VectorXd::Zero(m_nnu, 1);

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_dofs(m, i) < m_nnu) {
                dofval_u(m_dofs(m, i)) = nodevec(m, i);
            }
        }
    }

    return dofval_u;
}

inline Eigen::VectorXd
MatrixPartitionedTyings::AsDofs_p(const array_type::tensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    Eigen::VectorXd dofval_p(m_nnp, 1);

#pragma omp parallel for
    for (size_t d = 0; d < m_nnp; ++d) {
        dofval_p(d) = dofval(m_iip(d));
    }

    return dofval_p;
}

inline Eigen::VectorXd
MatrixPartitionedTyings::AsDofs_p(const array_type::tensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

    Eigen::VectorXd dofval_p = Eigen::VectorXd::Zero(m_nnp, 1);

#pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            if (m_dofs(m, i) >= m_nnu && m_dofs(m, i) < m_nni) {
                dofval_p(m_dofs(m, i) - m_nnu) = nodevec(m, i);
            }
        }
    }

    return dofval_p;
}

inline Eigen::VectorXd
MatrixPartitionedTyings::AsDofs_d(const array_type::tensor<double, 1>& dofval) const
{
    GOOSEFEM_ASSERT(dofval.size() == m_ndof);

    Eigen::VectorXd dofval_d(m_nnd, 1);

#pragma omp parallel for
    for (size_t d = 0; d < m_nnd; ++d) {
        dofval_d(d) = dofval(m_iip(d));
    }

    return dofval_d;
}

inline Eigen::VectorXd
MatrixPartitionedTyings::AsDofs_d(const array_type::tensor<double, 2>& nodevec) const
{
    GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

    Eigen::VectorXd dofval_d = Eigen::VectorXd::Zero(m_nnd, 1);

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

template <class Solver>
inline void MatrixPartitionedTyingsSolver<Solver>::factorize(MatrixPartitionedTyings& matrix)
{
    if (!matrix.m_changed && !m_factor) {
        return;
    }

    matrix.m_ACuu = matrix.m_Auu + matrix.m_Aud * matrix.m_Cdu + matrix.m_Cud * matrix.m_Adu +
                    matrix.m_Cud * matrix.m_Add * matrix.m_Cdu;

    matrix.m_ACup = matrix.m_Aup + matrix.m_Aud * matrix.m_Cdp + matrix.m_Cud * matrix.m_Adp +
                    matrix.m_Cud * matrix.m_Add * matrix.m_Cdp;

    // matrix.m_ACpu = matrix.m_Apu + matrix.m_Apd * matrix.m_Cdu + matrix.m_Cpd * matrix.m_Adu
    //     + matrix.m_Cpd * matrix.m_Add * matrix.m_Cdu;

    // matrix.m_ACpp = matrix.m_App + matrix.m_Apd * matrix.m_Cdp + matrix.m_Cpd * matrix.m_Adp
    //     + matrix.m_Cpd * matrix.m_Add * matrix.m_Cdp;

    m_solver.compute(matrix.m_ACuu);
    m_factor = false;
    matrix.m_changed = false;
}

template <class Solver>
inline void MatrixPartitionedTyingsSolver<Solver>::solve(
    MatrixPartitionedTyings& matrix,
    const array_type::tensor<double, 2>& b,
    array_type::tensor<double, 2>& x)
{
    GOOSEFEM_ASSERT(xt::has_shape(b, {matrix.m_nnode, matrix.m_ndim}));
    GOOSEFEM_ASSERT(xt::has_shape(x, {matrix.m_nnode, matrix.m_ndim}));

    this->factorize(matrix);

    Eigen::VectorXd B_u = matrix.AsDofs_u(b);
    Eigen::VectorXd B_d = matrix.AsDofs_d(b);
    Eigen::VectorXd X_p = matrix.AsDofs_p(x);

    B_u += matrix.m_Cud * B_d;

    Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - matrix.m_ACup * X_p));
    Eigen::VectorXd X_d = matrix.m_Cdu * X_u + matrix.m_Cdp * X_p;

#pragma omp parallel for
    for (size_t m = 0; m < matrix.m_nnode; ++m) {
        for (size_t i = 0; i < matrix.m_ndim; ++i) {
            if (matrix.m_dofs(m, i) < matrix.m_nnu) {
                x(m, i) = X_u(matrix.m_dofs(m, i));
            }
            else if (matrix.m_dofs(m, i) >= matrix.m_nni) {
                x(m, i) = X_d(matrix.m_dofs(m, i) - matrix.m_nni);
            }
        }
    }
}

template <class Solver>
inline void MatrixPartitionedTyingsSolver<Solver>::solve(
    MatrixPartitionedTyings& matrix,
    const array_type::tensor<double, 1>& b,
    array_type::tensor<double, 1>& x)
{
    GOOSEFEM_ASSERT(b.size() == matrix.m_ndof);
    GOOSEFEM_ASSERT(x.size() == matrix.m_ndof);

    this->factorize(matrix);

    Eigen::VectorXd B_u = matrix.AsDofs_u(b);
    Eigen::VectorXd B_d = matrix.AsDofs_d(b);
    Eigen::VectorXd X_p = matrix.AsDofs_p(x);

    Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - matrix.m_ACup * X_p));
    Eigen::VectorXd X_d = matrix.m_Cdu * X_u + matrix.m_Cdp * X_p;

#pragma omp parallel for
    for (size_t d = 0; d < matrix.m_nnu; ++d) {
        x(matrix.m_iiu(d)) = X_u(d);
    }

#pragma omp parallel for
    for (size_t d = 0; d < matrix.m_nnd; ++d) {
        x(matrix.m_iid(d)) = X_d(d);
    }
}

template <class Solver>
inline void MatrixPartitionedTyingsSolver<Solver>::solve_u(
    MatrixPartitionedTyings& matrix,
    const array_type::tensor<double, 1>& b_u,
    const array_type::tensor<double, 1>& b_d,
    const array_type::tensor<double, 1>& x_p,
    array_type::tensor<double, 1>& x_u)
{
    UNUSED(b_d);
    GOOSEFEM_ASSERT(b_u.size() == matrix.m_nnu);
    GOOSEFEM_ASSERT(b_d.size() == matrix.m_nnd);
    GOOSEFEM_ASSERT(x_p.size() == matrix.m_nnp);
    GOOSEFEM_ASSERT(x_u.size() == matrix.m_nnu);

    this->factorize(matrix);

    Eigen::Map<Eigen::VectorXd>(x_u.data(), x_u.size()).noalias() = m_solver.solve(Eigen::VectorXd(
        Eigen::Map<const Eigen::VectorXd>(b_u.data(), b_u.size()) -
        matrix.m_ACup * Eigen::Map<const Eigen::VectorXd>(x_p.data(), x_p.size())));
}

template <class Solver>
inline array_type::tensor<double, 2> MatrixPartitionedTyingsSolver<Solver>::Solve(
    MatrixPartitionedTyings& matrix,
    const array_type::tensor<double, 2>& b,
    const array_type::tensor<double, 2>& x)
{
    array_type::tensor<double, 2> ret = x;
    this->solve(matrix, b, ret);
    return ret;
}

template <class Solver>
inline array_type::tensor<double, 1> MatrixPartitionedTyingsSolver<Solver>::Solve(
    MatrixPartitionedTyings& matrix,
    const array_type::tensor<double, 1>& b,
    const array_type::tensor<double, 1>& x)
{
    array_type::tensor<double, 1> ret = x;
    this->solve(matrix, b, ret);
    return ret;
}

template <class Solver>
inline array_type::tensor<double, 1> MatrixPartitionedTyingsSolver<Solver>::Solve_u(
    MatrixPartitionedTyings& matrix,
    const array_type::tensor<double, 1>& b_u,
    const array_type::tensor<double, 1>& b_d,
    const array_type::tensor<double, 1>& x_p)
{
    array_type::tensor<double, 1> x_u = xt::empty<double>({matrix.m_nnu});
    this->solve_u(matrix, b_u, b_d, x_p, x_u);
    return x_u;
}

} // namespace GooseFEM

#endif

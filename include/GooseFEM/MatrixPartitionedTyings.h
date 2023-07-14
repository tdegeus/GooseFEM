/**
 * Sparse matrix that is partitioned in:
 * -   unknown DOFs
 * -   prescribed DOFs
 * -   tied DOFs
 *
 * @file MatrixPartitionedTyings.h
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#ifndef GOOSEFEM_MATRIXPARTITIONEDTYINGS_H
#define GOOSEFEM_MATRIXPARTITIONEDTYINGS_H

#include "Matrix.h"
#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace GooseFEM {

// forward declaration
template <class>
class MatrixPartitionedTyingsSolver;

/**
 * Sparse matrix from with dependent DOFs are eliminated,
 * and the remaining (small) independent system is partitioned in an unknown and a prescribed part.
 * In particular:
 *
 * \f$ A_{ii} = \begin{bmatrix} A_{uu} & A_{up} \\ A_{pu} & A_{pp} \end{bmatrix} \f$
 *
 * See VectorPartitionedTyings() for bookkeeping definitions.
 */
class MatrixPartitionedTyings : public MatrixPartitionedTyingsBase<MatrixPartitionedTyings> {
private:
    friend MatrixBase<MatrixPartitionedTyings>;
    friend MatrixPartitionedBase<MatrixPartitionedTyings>;
    friend MatrixPartitionedTyingsBase<MatrixPartitionedTyings>;

protected:
    Eigen::SparseMatrix<double> m_Auu; ///< The matrix.
    Eigen::SparseMatrix<double> m_Aup; ///< The matrix.
    Eigen::SparseMatrix<double> m_Apu; ///< The matrix.
    Eigen::SparseMatrix<double> m_App; ///< The matrix.
    Eigen::SparseMatrix<double> m_Aud; ///< The matrix.
    Eigen::SparseMatrix<double> m_Apd; ///< The matrix.
    Eigen::SparseMatrix<double> m_Adu; ///< The matrix.
    Eigen::SparseMatrix<double> m_Adp; ///< The matrix.
    Eigen::SparseMatrix<double> m_Add; ///< The matrix.
    Eigen::SparseMatrix<double> m_ACuu; ///< // The matrix for which the tyings have been applied.
    Eigen::SparseMatrix<double> m_ACup; ///< // The matrix for which the tyings have been applied.
    Eigen::SparseMatrix<double> m_ACpu; ///< // The matrix for which the tyings have been applied.
    Eigen::SparseMatrix<double> m_ACpp; ///< // The matrix for which the tyings have been applied.
    std::vector<Eigen::Triplet<double>> m_Tuu; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tup; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tpu; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tpp; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tud; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tpd; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tdu; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tdp; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tdd; ///< Matrix entries.
    Eigen::SparseMatrix<double> m_Cdu; ///< Tying matrix, see Tyings::Periodic::Cdu().
    Eigen::SparseMatrix<double> m_Cdp; ///< Tying matrix, see Tyings::Periodic::Cdp().
    Eigen::SparseMatrix<double> m_Cud; ///< Transpose of "m_Cdu".
    Eigen::SparseMatrix<double> m_Cpd; ///< Transpose of "m_Cdp".

    // grant access to solver class
    template <class>
    friend class MatrixPartitionedTyingsSolver;

public:
    MatrixPartitionedTyings() = default;

    /**
     * Constructor.
     *
     * @param conn connectivity [#nelem, #nne].
     * @param dofs DOFs per node [#nnode, #ndim].
     * @param Cdu See Tyings::Periodic::Cdu().
     * @param Cdp See Tyings::Periodic::Cdp().
     */
    MatrixPartitionedTyings(
        const array_type::tensor<size_t, 2>& conn,
        const array_type::tensor<size_t, 2>& dofs,
        const Eigen::SparseMatrix<double>& Cdu,
        const Eigen::SparseMatrix<double>& Cdp
    )
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
        m_iii = xt::arange<size_t>(m_nni);
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

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_uu() const
    {
        return m_Auu;
    }

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_up() const
    {
        return m_Aup;
    }

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_pu() const
    {
        return m_Apu;
    }

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_pp() const
    {
        return m_App;
    }

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_ud() const
    {
        return m_Aud;
    }

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_pd() const
    {
        return m_Apd;
    }

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_du() const
    {
        return m_Adu;
    }

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_dp() const
    {
        return m_Adp;
    }

    /**
     * Pointer to data.
     */
    const Eigen::SparseMatrix<double>& data_dd() const
    {
        return m_Add;
    }

private:
    template <class T>
    void assemble_impl(const T& elemmat)
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
                                    di, dj, elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
                            }
                            else if (di < m_nnu && dj < m_nni) {
                                m_Tup.push_back(Eigen::Triplet<double>(
                                    di, dj - m_nnu, elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
                            }
                            else if (di < m_nnu) {
                                m_Tud.push_back(Eigen::Triplet<double>(
                                    di, dj - m_nni, elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
                            }
                            else if (di < m_nni && dj < m_nnu) {
                                m_Tpu.push_back(Eigen::Triplet<double>(
                                    di - m_nnu, dj, elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
                            }
                            else if (di < m_nni && dj < m_nni) {
                                m_Tpp.push_back(Eigen::Triplet<double>(
                                    di - m_nnu,
                                    dj - m_nnu,
                                    elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
                            }
                            else if (di < m_nni) {
                                m_Tpd.push_back(Eigen::Triplet<double>(
                                    di - m_nnu,
                                    dj - m_nni,
                                    elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
                            }
                            else if (dj < m_nnu) {
                                m_Tdu.push_back(Eigen::Triplet<double>(
                                    di - m_nni, dj, elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
                            }
                            else if (dj < m_nni) {
                                m_Tdp.push_back(Eigen::Triplet<double>(
                                    di - m_nni,
                                    dj - m_nnu,
                                    elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
                            }
                            else {
                                m_Tdd.push_back(Eigen::Triplet<double>(
                                    di - m_nni,
                                    dj - m_nni,
                                    elemmat(e, m * m_ndim + i, n * m_ndim + j)
                                ));
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

    // todo: test
    template <class T>
    void todense_impl(T& ret) const
    {
        ret.fill(0.0);

        for (int k = 0; k < m_Auu.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_Auu, k); it; ++it) {
                ret(m_iiu(it.row()), m_iiu(it.col())) = it.value();
            }
        }

        for (int k = 0; k < m_Aup.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_Aup, k); it; ++it) {
                ret(m_iiu(it.row()), m_iip(it.col())) = it.value();
            }
        }

        for (int k = 0; k < m_Apu.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_Apu, k); it; ++it) {
                ret(m_iip(it.row()), m_iiu(it.col())) = it.value();
            }
        }

        for (int k = 0; k < m_App.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_App, k); it; ++it) {
                ret(m_iip(it.row()), m_iip(it.col())) = it.value();
            }
        }

        for (int k = 0; k < m_Adu.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_Adu, k); it; ++it) {
                ret(m_iid(it.row()), m_iiu(it.col())) = it.value();
            }
        }

        for (int k = 0; k < m_Adp.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_Adp, k); it; ++it) {
                ret(m_iid(it.row()), m_iip(it.col())) = it.value();
            }
        }

        for (int k = 0; k < m_Aud.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_Aud, k); it; ++it) {
                ret(m_iiu(it.row()), m_iid(it.col())) = it.value();
            }
        }

        for (int k = 0; k < m_Apd.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_Apd, k); it; ++it) {
                ret(m_iip(it.row()), m_iid(it.col())) = it.value();
            }
        }

        for (int k = 0; k < m_Add.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_Add, k); it; ++it) {
                ret(m_iid(it.row()), m_iid(it.col())) = it.value();
            }
        }
    }

    template <class T>
    void dot_nodevec_impl(const T& x, T& b) const
    {
        UNUSED(x);
        UNUSED(b);
        throw std::runtime_error("Not yet implemented");
    }

    template <class T>
    void dot_dofval_impl(const T& x, T& b) const
    {
        UNUSED(x);
        UNUSED(b);
        throw std::runtime_error("Not yet implemented");
    }

    template <class T>
    void reaction_nodevec_impl(const T& x, T& b) const
    {
        UNUSED(x);
        UNUSED(b);
        throw std::runtime_error("Not yet implemented");
    }

    template <class T>
    void reaction_dofval_impl(const T& x, T& b) const
    {
        UNUSED(x);
        UNUSED(b);
        throw std::runtime_error("Not yet implemented");
    }

    void reaction_p_impl(
        const array_type::tensor<double, 1>& x_u,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& b_p
    ) const
    {
        UNUSED(x_u);
        UNUSED(x_p);
        UNUSED(b_p);
        throw std::runtime_error("Not yet implemented");
    }

private:
    // Convert arrays (Eigen version of VectorPartitioned, which contains public functions)
    Eigen::VectorXd AsDofs_u(const array_type::tensor<double, 1>& dofval) const
    {
        GOOSEFEM_ASSERT(dofval.size() == m_ndof);

        Eigen::VectorXd dofval_u(m_nnu, 1);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            dofval_u(d) = dofval(m_iiu(d));
        }

        return dofval_u;
    }

    Eigen::VectorXd AsDofs_u(const array_type::tensor<double, 2>& nodevec) const
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

    Eigen::VectorXd AsDofs_p(const array_type::tensor<double, 1>& dofval) const
    {
        GOOSEFEM_ASSERT(dofval.size() == m_ndof);

        Eigen::VectorXd dofval_p(m_nnp, 1);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            dofval_p(d) = dofval(m_iip(d));
        }

        return dofval_p;
    }

    Eigen::VectorXd AsDofs_p(const array_type::tensor<double, 2>& nodevec) const
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

    Eigen::VectorXd AsDofs_d(const array_type::tensor<double, 1>& dofval) const
    {
        GOOSEFEM_ASSERT(dofval.size() == m_ndof);

        Eigen::VectorXd dofval_d(m_nnd, 1);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnd; ++d) {
            dofval_d(d) = dofval(m_iip(d));
        }

        return dofval_d;
    }

    Eigen::VectorXd AsDofs_d(const array_type::tensor<double, 2>& nodevec) const
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
};

/**
 * Solver for MatrixPartitionedTyings().
 * This solver class can be used to solve for multiple right-hand-sides using one factorisation.
 *
 * Solving proceeds as follows:
 *
 * \f$ A' = A_{ii} + A_{id} * C_{di} + C_{di}^T * A_{di} + C_{di}^T * A_{dd} * C_{di} \f$
 *
 * \f$ b' = b_i + C_{di}^T * b_d \f$
 *
 * \f$ x_u = A'_{uu} \ ( b'_u - A'_{up} * x_p ) \f$
 *
 * \f$ x_i = \begin{bmatrix} x_u \\ x_p \end{bmatrix} \f$
 *
 * \f$ x_d = C_{di} * x_i \f$
 */
template <class Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
class MatrixPartitionedTyingsSolver
    : public MatrixSolverBase<MatrixPartitionedTyingsSolver<Solver>>,
      public MatrixSolverPartitionedBase<MatrixPartitionedTyingsSolver<Solver>> {
private:
    friend MatrixSolverBase<MatrixPartitionedTyingsSolver<Solver>>;
    friend MatrixSolverPartitionedBase<MatrixPartitionedTyingsSolver<Solver>>;

public:
    MatrixPartitionedTyingsSolver() = default;

private:
    template <class T>
    void solve_nodevec_impl(MatrixPartitionedTyings& A, const T& b, T& x)
    {
        this->factorize(A);

        Eigen::VectorXd B_u = A.AsDofs_u(b);
        Eigen::VectorXd B_d = A.AsDofs_d(b);
        Eigen::VectorXd X_p = A.AsDofs_p(x);

        B_u += A.m_Cud * B_d;

        Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - A.m_ACup * X_p));
        Eigen::VectorXd X_d = A.m_Cdu * X_u + A.m_Cdp * X_p;

#pragma omp parallel for
        for (size_t m = 0; m < A.m_nnode; ++m) {
            for (size_t i = 0; i < A.m_ndim; ++i) {
                if (A.m_dofs(m, i) < A.m_nnu) {
                    x(m, i) = X_u(A.m_dofs(m, i));
                }
                else if (A.m_dofs(m, i) >= A.m_nni) {
                    x(m, i) = X_d(A.m_dofs(m, i) - A.m_nni);
                }
            }
        }
    }

    template <class T>
    void solve_dofval_impl(MatrixPartitionedTyings& A, const T& b, T& x)
    {
        this->factorize(A);

        Eigen::VectorXd B_u = A.AsDofs_u(b);
        Eigen::VectorXd B_d = A.AsDofs_d(b);
        Eigen::VectorXd X_p = A.AsDofs_p(x);

        Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - A.m_ACup * X_p));
        Eigen::VectorXd X_d = A.m_Cdu * X_u + A.m_Cdp * X_p;

#pragma omp parallel for
        for (size_t d = 0; d < A.m_nnu; ++d) {
            x(A.m_iiu(d)) = X_u(d);
        }

#pragma omp parallel for
        for (size_t d = 0; d < A.m_nnd; ++d) {
            x(A.m_iid(d)) = X_d(d);
        }
    }

public:
    /**
     * Same as
     * Solve(MatrixPartitionedTyings&, const array_type::tensor<double, 2>&, const
     * array_type::tensor<double, 2>&), but with partitioned input and output.
     *
     * @param A sparse matrix, see MatrixPartitionedTyings().
     * @param b_u unknown dofval [nnu].
     * @param b_d dependent dofval [nnd].
     * @param x_p prescribed dofval [nnp]
     * @return x_u unknown dofval [nnu].
     */
    array_type::tensor<double, 1> Solve_u(
        MatrixPartitionedTyings& A,
        const array_type::tensor<double, 1>& b_u,
        const array_type::tensor<double, 1>& b_d,
        const array_type::tensor<double, 1>& x_p
    )
    {
        array_type::tensor<double, 1> x_u = xt::empty<double>({A.m_nnu});
        this->solve_u(A, b_u, b_d, x_p, x_u);
        return x_u;
    }

    /**
     * Same as
     * Solve_u(MatrixPartitionedTyings&, const array_type::tensor<double, 1>&, const
     * array_type::tensor<double, 1>&, const array_type::tensor<double, 1>&), but writing to
     * pre-allocated output.
     *
     * @param A sparse matrix, see MatrixPartitionedTyings().
     * @param b_u unknown dofval [nnu].
     * @param b_d dependent dofval [nnd].
     * @param x_p prescribed dofval [nnp]
     * @param x_u (overwritten) unknown dofval [nnu].
     */
    void solve_u(
        MatrixPartitionedTyings& A,
        const array_type::tensor<double, 1>& b_u,
        const array_type::tensor<double, 1>& b_d,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& x_u
    )
    {
        UNUSED(b_d);
        GOOSEFEM_ASSERT(b_u.size() == A.m_nnu);
        GOOSEFEM_ASSERT(b_d.size() == A.m_nnd);
        GOOSEFEM_ASSERT(x_p.size() == A.m_nnp);
        GOOSEFEM_ASSERT(x_u.size() == A.m_nnu);

        this->factorize(A);

        Eigen::Map<Eigen::VectorXd>(x_u.data(), x_u.size()).noalias() =
            m_solver.solve(Eigen::VectorXd(
                Eigen::Map<const Eigen::VectorXd>(b_u.data(), b_u.size()) -
                A.m_ACup * Eigen::Map<const Eigen::VectorXd>(x_p.data(), x_p.size())
            ));
    }

private:
    Solver m_solver; ///< solver
    bool m_factor = true; ///< signal to force factorization

    /**
     * compute inverse (evaluated by "solve")
     */
    void factorize(MatrixPartitionedTyings& A)
    {
        if (!A.m_changed && !m_factor) {
            return;
        }

        A.m_ACuu = A.m_Auu + A.m_Aud * A.m_Cdu + A.m_Cud * A.m_Adu + A.m_Cud * A.m_Add * A.m_Cdu;

        A.m_ACup = A.m_Aup + A.m_Aud * A.m_Cdp + A.m_Cud * A.m_Adp + A.m_Cud * A.m_Add * A.m_Cdp;

        // A.m_ACpu = A.m_Apu + A.m_Apd * A.m_Cdu + A.m_Cpd * A.m_Adu
        //     + A.m_Cpd * A.m_Add * A.m_Cdu;

        // A.m_ACpp = A.m_App + A.m_Apd * A.m_Cdp + A.m_Cpd * A.m_Adp
        //     + A.m_Cpd * A.m_Add * A.m_Cdp;

        m_solver.compute(A.m_ACuu);
        m_factor = false;
        A.m_changed = false;
    }
};

} // namespace GooseFEM

#endif

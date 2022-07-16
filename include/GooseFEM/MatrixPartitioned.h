/**
Sparse matrix that is partitioned in:
-   unknown DOFs
-   prescribed DOFs

\file MatrixPartitioned.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXPARTITIONED_H
#define GOOSEFEM_MATRIXPARTITIONED_H

#include "Matrix.h"
#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace GooseFEM {

// forward declaration
template <class>
class MatrixPartitionedSolver;

/**
Sparse matrix partitioned in an unknown and a prescribed part.
In particular:
\f$ \begin{bmatrix} A_{uu} & A_{up} \\ A_{pu} & A_{pp} \end{bmatrix} \f$

See VectorPartitioned() for bookkeeping definitions.
*/
class MatrixPartitioned : public MatrixPartitionedBase<MatrixPartitioned> {
private:
    friend MatrixBase<MatrixPartitioned>;
    friend MatrixPartitionedBase<MatrixPartitioned>;

protected:
    Eigen::SparseMatrix<double> m_Auu; ///< The matrix.
    Eigen::SparseMatrix<double> m_Aup; ///< The matrix.
    Eigen::SparseMatrix<double> m_Apu; ///< The matrix.
    Eigen::SparseMatrix<double> m_App; ///< The matrix.

    std::vector<Eigen::Triplet<double>> m_Tuu; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tup; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tpu; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tpp; ///< Matrix entries.

    /**
    Renumbered DOFs per node, such that:

        iiu = arange(nnu)
        iip = nnu + arange(nnp)

    making is much simpler to slice.
    */
    array_type::tensor<size_t, 2> m_part;

    /**
    Map real DOF to DOF in partitioned system.
    The partitioned system is defined as:

        iiu = arange(nnu)
        iip = nnu + arange(nnp)

    Similar to `m_part` but for a 1d sequential list of DOFs.
    */
    array_type::tensor<size_t, 1> m_part1d;

    /**
    Class to solve the system (allowing single factorisation for multiple right-hand-sides).
    */
    template <class>
    friend class MatrixPartitionedSolver;

public:
    MatrixPartitioned() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    \param iip prescribed DOFs [#nnp].
    */
    MatrixPartitioned(
        const array_type::tensor<size_t, 2>& conn,
        const array_type::tensor<size_t, 2>& dofs,
        const array_type::tensor<size_t, 1>& iip)
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
        m_part1d = Mesh::Reorder({m_iiu, m_iip}).apply(xt::eval(xt::arange<size_t>(m_ndof)));
        m_Tuu.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
        m_Tup.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
        m_Tpu.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
        m_Tpp.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
        m_Auu.resize(m_nnu, m_nnu);
        m_Aup.resize(m_nnu, m_nnp);
        m_Apu.resize(m_nnp, m_nnu);
        m_App.resize(m_nnp, m_nnp);
    }

    /**
    Pointer to data.
    */
    const Eigen::SparseMatrix<double>& data_uu() const
    {
        return m_Auu;
    }

    /**
    Pointer to data.
    */
    const Eigen::SparseMatrix<double>& data_up() const
    {
        return m_Aup;
    }

    /**
    Pointer to data.
    */
    const Eigen::SparseMatrix<double>& data_pu() const
    {
        return m_Apu;
    }

    /**
    Pointer to data.
    */
    const Eigen::SparseMatrix<double>& data_pp() const
    {
        return m_App;
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

public:
    /**
    Overwrite matrix.

    \param rows Row numbers [m].
    \param cols Column numbers [n].
    \param matrix Data entries `matrix(i, j)` for `rows(i), cols(j)` [m, n].
    */
    void
    set(const array_type::tensor<size_t, 1>& rows,
        const array_type::tensor<size_t, 1>& cols,
        const array_type::tensor<double, 2>& matrix)
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
                size_t di = m_part1d(rows(i));
                size_t dj = m_part1d(cols(j));
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

    /**
    Add matrix.

    \param rows Row numbers [m].
    \param cols Column numbers [n].
    \param matrix Data entries `matrix(i, j)` for `rows(i), cols(j)` [m, n].
    */
    void
    add(const array_type::tensor<size_t, 1>& rows,
        const array_type::tensor<size_t, 1>& cols,
        const array_type::tensor<double, 2>& matrix)
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
                size_t di = m_part1d(rows(i));
                size_t dj = m_part1d(cols(j));
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

private:
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
    }

    template <class T>
    void dot_nodevec_impl(const T& x, T& b) const
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
                else {
                    b(m, i) = B_p(m_part(m, i) - m_nnu);
                }
            }
        }
    }

    template <class T>
    void dot_dofval_impl(const T& x, T& b) const
    {
        GOOSEFEM_ASSERT(b.size() == m_ndof);
        GOOSEFEM_ASSERT(x.size() == m_ndof);

        Eigen::VectorXd X_u = this->AsDofs_u(x);
        Eigen::VectorXd X_p = this->AsDofs_p(x);

        Eigen::VectorXd B_u = m_Auu * X_u + m_Aup * X_p;
        Eigen::VectorXd B_p = m_Apu * X_u + m_App * X_p;

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            b(m_iiu(d)) = B_u(d);
        }

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            b(m_iip(d)) = B_p(d);
        }
    }

    template <class T>
    void reaction_nodevec_impl(const T& x, T& b) const
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

    template <class T>
    void reaction_dofval_impl(const T& x, T& b) const
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

    void reaction_p_impl(
        const array_type::tensor<double, 1>& x_u,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& b_p) const
    {
        GOOSEFEM_ASSERT(x_u.size() == m_nnu);
        GOOSEFEM_ASSERT(x_p.size() == m_nnp);
        GOOSEFEM_ASSERT(b_p.size() == m_nnp);

        Eigen::Map<Eigen::VectorXd>(b_p.data(), b_p.size()).noalias() =
            m_Apu * Eigen::Map<const Eigen::VectorXd>(x_u.data(), x_u.size()) +
            m_App * Eigen::Map<const Eigen::VectorXd>(x_p.data(), x_p.size());
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
                if (m_part(m, i) < m_nnu) {
                    dofval_u(m_part(m, i)) = nodevec(m, i);
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
                if (m_part(m, i) >= m_nnu) {
                    dofval_p(m_part(m, i) - m_nnu) = nodevec(m, i);
                }
            }
        }

        return dofval_p;
    }
};

/**
Solve \f$ x_u = A_{uu}^{-1} (b_u - A_{up} * x_p) \f$ for `A` of the MatrixPartitioned() class.
You can solve for multiple right-hand-sides using one factorisation.

For "nodevec" input `x` is used to read \f$ x_p \f$, while \f$ x_u \f$ is written.
See MatrixPartitioned::Reaction() to get \f$ b_p \f$.
*/
template <class Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
class MatrixPartitionedSolver
    : public MatrixSolverBase<MatrixPartitionedSolver<Solver>>,
      public MatrixSolverPartitionedBase<MatrixPartitionedSolver<Solver>> {
private:
    friend MatrixSolverBase<MatrixPartitionedSolver<Solver>>;
    friend MatrixSolverPartitionedBase<MatrixPartitionedSolver<Solver>>;

public:
    MatrixPartitionedSolver() = default;

    /**
    Solve \f$ x = A^{-1} b \f$.

    \param A GooseFEM (sparse) matrix, see e.g. GooseFEM::MatrixPartitioned().
    \param b_u unknown dofval [nnu].
    \param x_p prescribed dofval [nnp]
    \return x_u unknown dofval [nnu].
    */
    template <class M>
    array_type::tensor<double, 1> Solve_u(
        M& A,
        const array_type::tensor<double, 1>& b_u,
        const array_type::tensor<double, 1>& x_p)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b_u, {A.nnu()}));
        GOOSEFEM_ASSERT(xt::has_shape(x_p, {A.nnp()}));
        array_type::tensor<double, 1> x_u = xt::empty_like(b_u);
        this->solve_u_impl(A, b_u, x_p, x_u);
        return x_u;
    }

    /**
    Same as
    Solve \f$ x = A^{-1} b \f$.

    \param A GooseFEM (sparse) matrix, see e.g. GooseFEM::MatrixPartitioned().
    \param b_u unknown dofval [nnu].
    \param x_p prescribed dofval [nnp]
    \param x_u (overwritten) unknown dofval [nnu].
    */
    template <class M>
    void solve_u(
        M& A,
        const array_type::tensor<double, 1>& b_u,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& x_u)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b_u, {A.nnu()}));
        GOOSEFEM_ASSERT(xt::has_shape(x_p, {A.nnp()}));
        GOOSEFEM_ASSERT(xt::has_shape(x_u, {A.nnu()}));
        this->solve_u_impl(A, b_u, x_p, x_u);
    }

private:
    template <class T>
    void solve_nodevec_impl(MatrixPartitioned& A, const T& b, T& x)
    {
        this->factorize(A);
        Eigen::VectorXd B_u = A.AsDofs_u(b);
        Eigen::VectorXd X_p = A.AsDofs_p(x);
        Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - A.m_Aup * X_p));

#pragma omp parallel for
        for (size_t m = 0; m < A.m_nnode; ++m) {
            for (size_t i = 0; i < A.m_ndim; ++i) {
                if (A.m_part(m, i) < A.m_nnu) {
                    x(m, i) = X_u(A.m_part(m, i));
                }
            }
        }
    }

    template <class T>
    void solve_dofval_impl(MatrixPartitioned& A, const T& b, T& x)
    {
        this->factorize(A);
        Eigen::VectorXd B_u = A.AsDofs_u(b);
        Eigen::VectorXd X_p = A.AsDofs_p(x);
        Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - A.m_Aup * X_p));

#pragma omp parallel for
        for (size_t d = 0; d < A.m_nnu; ++d) {
            x(A.m_iiu(d)) = X_u(d);
        }
    }

    template <class T>
    void solve_u_impl(MatrixPartitioned& A, const T& b_u, const T& x_p, T& x_u)
    {
        this->factorize(A);

        Eigen::Map<Eigen::VectorXd>(x_u.data(), x_u.size()).noalias() =
            m_solver.solve(Eigen::VectorXd(
                Eigen::Map<const Eigen::VectorXd>(b_u.data(), b_u.size()) -
                A.m_Aup * Eigen::Map<const Eigen::VectorXd>(x_p.data(), x_p.size())));
    }

private:
    Solver m_solver; ///< solver
    bool m_factor = true; ///< signal to force factorization

    /**
    compute inverse (evaluated by "solve")
    */
    void factorize(MatrixPartitioned& A)
    {
        if (!A.m_changed && !m_factor) {
            return;
        }
        m_solver.compute(A.m_Auu);
        m_factor = false;
        A.m_changed = false;
    }
};

} // namespace GooseFEM

#endif

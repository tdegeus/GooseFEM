/**
Sparse matrix that is partitioned in:
-   unknown DOFs
-   prescribed DOFs
-   tied DOFs

\file MatrixPartitionedTyings.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXPARTITIONEDTYINGS_H
#define GOOSEFEM_MATRIXPARTITIONEDTYINGS_H

#include "config.h"
#include "Matrix.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace GooseFEM {

// forward declaration
template <class> class MatrixPartitionedTyingsSolver;

/**
Sparse matrix from with dependent DOFs are eliminated,
and the remaining (small) independent system is partitioned in an unknown and a prescribed part.
In particular:

\f$ A_{ii} = \begin{bmatrix} A_{uu} & A_{up} \\ A_{pu} & A_{pp} \end{bmatrix} \f$

See VectorPartitionedTyings() for bookkeeping definitions.
*/
class MatrixPartitionedTyings : public Matrix {
public:

    MatrixPartitionedTyings() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    \param Cdu See Tyings::Periodic::Cdu().
    \param Cdp See Tyings::Periodic::Cdp().
    */
    MatrixPartitionedTyings(
        const xt::xtensor<size_t, 2>& conn,
        const xt::xtensor<size_t, 2>& dofs,
        const Eigen::SparseMatrix<double>& Cdu,
        const Eigen::SparseMatrix<double>& Cdp);

    /**
    \return Number of dependent DOFs.
    */
    size_t nnd() const;

    /**
    \return Number of independent DOFs.
    */
    size_t nni() const;

    /**
    \return Number of independent unknown DOFs.
    */
    size_t nnu() const;

    /**
    \return Number of independent prescribed DOFs.
    */
    size_t nnp() const;

    /**
    Dependent DOFs.

    \return List of DOF numbers.
    */
    xt::xtensor<size_t, 1> iid() const;

    /**
    Independent DOFs.

    \return List of DOF numbers.
    */
    xt::xtensor<size_t, 1> iii() const;

    /**
    Independent unknown DOFs.

    \return List of DOF numbers.
    */
    xt::xtensor<size_t, 1> iiu() const;

    /**
    Independent prescribed DOFs.

    \return List of DOF numbers.
    */
    xt::xtensor<size_t, 1> iip() const;

    void assemble(const xt::xtensor<double, 3>& elemmat) override;

private:

    using Matrix::set;
    using Matrix::add;
    using Matrix::Todense;
    using Matrix::todense;
    using Matrix::Dot;
    using Matrix::dot;

private:

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
    xt::xtensor<size_t, 1> m_iiu; ///< See iiu()
    xt::xtensor<size_t, 1> m_iip; ///< See iip()
    xt::xtensor<size_t, 1> m_iid; ///< See iid()
    size_t m_nnu; ///< See #nnu
    size_t m_nnp; ///< See #nnp
    size_t m_nni; ///< See #nni
    size_t m_nnd; ///< See #nnd
    Eigen::SparseMatrix<double> m_Cdu; ///< Tying matrix, see Tyings::Periodic::Cdu().
    Eigen::SparseMatrix<double> m_Cdp; ///< Tying matrix, see Tyings::Periodic::Cdp().
    Eigen::SparseMatrix<double> m_Cud; ///< Transpose of "m_Cdu".
    Eigen::SparseMatrix<double> m_Cpd; ///< Transpose of "m_Cdp".

    // grant access to solver class
    template <class> friend class MatrixPartitionedTyingsSolver;

private:

    // Convert arrays (Eigen version of VectorPartitioned, which contains public functions)
    Eigen::VectorXd AsDofs_u(const xt::xtensor<double, 1>& dofval) const;
    Eigen::VectorXd AsDofs_u(const xt::xtensor<double, 2>& nodevec) const;
    Eigen::VectorXd AsDofs_p(const xt::xtensor<double, 1>& dofval) const;
    Eigen::VectorXd AsDofs_p(const xt::xtensor<double, 2>& nodevec) const;
    Eigen::VectorXd AsDofs_d(const xt::xtensor<double, 1>& dofval) const;
    Eigen::VectorXd AsDofs_d(const xt::xtensor<double, 2>& nodevec) const;
};

/**
Solver for MatrixPartitionedTyings().
The idea is that this solver class can be used to solve for multiple right-hand-sides
using one factorisation.
*/
template <class Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
class MatrixPartitionedTyingsSolver {
public:

    MatrixPartitionedTyingsSolver() = default;

    /**
    Solve as follows.

    \f$ A' = A_{ii} + A_{id} * C_{di} + C_{di}^T * A_{di} + C_{di}^T * A_{dd} * C_{di} \f$

    \f$ b' = b_i + C_{di}^T * b_d \f$

    \f$ x_u = A'_{uu} \ ( b'_u - A'_{up} * x_p ) \f$

    \f$ x_i = \begin{bmatrix} x_u \\ x_p \end{bmatrix} \f$

    \f$ x_d = C_{di} * x_i \f$

    \param A sparse matrix, see MatrixPartitionedTyings().
    \param b nodevec [nelem, ndim].
    \param x nodevec [nelem, ndim], used to read \f$ x_p \f$.
    \return x nodevec [nelem, ndim], \f$ x_u \f$ filled, \f$ x_p \f$ copied.
    */
    xt::xtensor<double, 2> Solve(
        MatrixPartitionedTyings& A,
        const xt::xtensor<double, 2>& b,
        const xt::xtensor<double, 2>& x);

    /**
    Same as
    Solve(MatrixPartitionedTyings&, const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but filling \f$ x_u \f$ and \f$ x_d \f$ in place.

    \param A sparse matrix, see MatrixPartitionedTyings().
    \param b nodevec [nelem, ndim].
    \param x nodevec [nelem, ndim], \f$ x_p \f$ read, \f$ x_u \f$ and \f$ x_d \f$ filled.
    */
    void solve(
        MatrixPartitionedTyings& A,
        const xt::xtensor<double, 2>& b,
        xt::xtensor<double, 2>& x);

    /**
    Same as
    Solve(MatrixPartitionedTyings&, const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but for "dofval" input and output.

    \param A sparse matrix, see MatrixPartitionedTyings().
    \param b dofval [ndof].
    \param x dofval [ndof], used to read \f$ x_p \f$.
    \return x dofval [ndof], \f$ x_u \f$ and \f$ x_d \f$ filled, \f$ x_p \f$ copied.
    */
    xt::xtensor<double, 1> Solve(
        MatrixPartitionedTyings& A,
        const xt::xtensor<double, 1>& b,
        const xt::xtensor<double, 1>& x);

    /**
    Same as
    Solve(MatrixPartitionedTyings&, const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but filling \f$ x_u \f$ and \f$ x_d \f$ in place.

    \param A sparse matrix, see MatrixPartitionedTyings().
    \param b dofval [ndof].
    \param x dofval [ndof], \f$ x_p \f$ read, \f$ x_u \f$ and \f$ x_d \f$ filled.
    */
    void solve(
        MatrixPartitionedTyings& A,
        const xt::xtensor<double, 1>& b,
        xt::xtensor<double, 1>& x);

    /**
    Same as
    Solve(MatrixPartitionedTyings&, const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but with partitioned input and output.

    \param A sparse matrix, see MatrixPartitionedTyings().
    \param b_u unknown dofval [nnu].
    \param b_d dependent dofval [nnd].
    \param x_p prescribed dofval [nnp]
    \return x_u unknown dofval [nnu].
    */
    xt::xtensor<double, 1> Solve_u(
        MatrixPartitionedTyings& A,
        const xt::xtensor<double, 1>& b_u,
        const xt::xtensor<double, 1>& b_d,
        const xt::xtensor<double, 1>& x_p);

    /**
    Same as
    Solve_u(MatrixPartitionedTyings&, const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but writing to pre-allocated output.

    \param A sparse matrix, see MatrixPartitionedTyings().
    \param b_u unknown dofval [nnu].
    \param b_d dependent dofval [nnd].
    \param x_p prescribed dofval [nnp]
    \param x_u (overwritten) unknown dofval [nnu].
    */
    void solve_u(
        MatrixPartitionedTyings& A,
        const xt::xtensor<double, 1>& b_u,
        const xt::xtensor<double, 1>& b_d,
        const xt::xtensor<double, 1>& x_p,
        xt::xtensor<double, 1>& x_u);

private:
    Solver m_solver; ///< solver
    bool m_factor = true; ///< signal to force factorization
    void factorize(MatrixPartitionedTyings& matrix); ///< compute inverse (evaluated by "solve")
};

} // namespace GooseFEM

#include "MatrixPartitionedTyings.hpp"

#endif

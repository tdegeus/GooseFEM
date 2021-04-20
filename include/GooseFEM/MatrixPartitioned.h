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

#include "config.h"
#include "Matrix.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace GooseFEM {

// forward declaration
template <class> class MatrixPartitionedSolver;

/**
Sparse matrix partitioned in an unknown and a prescribed part. In particular:

\f$ \begin{bmatrix} A_{uu} & A_{up} \\ A_{pu} & A_{pp} \end{bmatrix} \f$

See VectorPartitioned() for bookkeeping definitions.
*/
class MatrixPartitioned : public Matrix {
public:

    MatrixPartitioned() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    \param iip prescribed DOFs [#nnp].
    */
    MatrixPartitioned(
        const xt::xtensor<size_t, 2>& conn,
        const xt::xtensor<size_t, 2>& dofs,
        const xt::xtensor<size_t, 1>& iip);

    /**
    Number of unknown DOFs.
    */
    size_t nnu() const;

    /**
    Number of prescribed DOFs.
    */
    size_t nnp() const;

    /**
    Unknown DOFs [#nnu].
    */
    xt::xtensor<size_t, 1> iiu() const;

    /**
    Prescribed DOFs [#nnp].
    */
    xt::xtensor<size_t, 1> iip() const;

    void assemble(const xt::xtensor<double, 3>& elemmat) override;

    void set(
        const xt::xtensor<size_t, 1>& rows,
        const xt::xtensor<size_t, 1>& cols,
        const xt::xtensor<double, 2>& matrix) override;

    void add(
        const xt::xtensor<size_t, 1>& rows,
        const xt::xtensor<size_t, 1>& cols,
        const xt::xtensor<double, 2>& matrix) override;

    void todense(xt::xtensor<double, 2>& ret) const override;
    void dot(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const override;
    void dot(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const override;

    /**
    Get right-hand-size for corresponding to the prescribed DOFs.

    \f$ b_p = A_{pu} * x_u + A_{pp} * x_p \f$

    and assemble them to the appropriate places in "nodevec".

    \param x "nodevec" [#nnode, #ndim].
    \param b "nodevec" [#nnode, #ndim].
    \return Copy of `b` with \f$ b_p \f$ overwritten.
    */
    xt::xtensor<double, 2> Reaction(
        const xt::xtensor<double, 2>& x,
        const xt::xtensor<double, 2>& b) const;

    /**
    Same as Reaction(const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but inserting in-place.

    \param x "nodevec" [#nnode, #ndim].
    \param b "nodevec" [#nnode, #ndim], \f$ b_p \f$ overwritten.
    */
    void reaction(
        const xt::xtensor<double, 2>& x,
              xt::xtensor<double, 2>& b) const;

    /**
    Same as Reaction(const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but of "dofval" input and output.

    \param x "dofval" [#ndof].
    \param b "dofval" [#ndof].
    \return Copy of `b` with \f$ b_p \f$ overwritten.
    */
    xt::xtensor<double, 1> Reaction(
        const xt::xtensor<double, 1>& x,
        const xt::xtensor<double, 1>& b) const;

    /**
    Same as Reaction(const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but inserting in-place.

    \param x "dofval" [#ndof].
    \param b "dofval" [#ndof], \f$ b_p \f$ overwritten.
    */
    void reaction(
        const xt::xtensor<double, 1>& x,
              xt::xtensor<double, 1>& b) const;

    /**
    Same as Reaction(const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but with partitioned input and output.

    \param x_u unknown "dofval" [#nnu].
    \param x_p prescribed "dofval" [#nnp].
    \return b_p prescribed "dofval" [#nnp].
    */
    xt::xtensor<double, 1> Reaction_p(
        const xt::xtensor<double, 1>& x_u,
        const xt::xtensor<double, 1>& x_p) const;

    /**
    Same as Reaction_p(const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but writing to preallocated output.

    \param x_u unknown "dofval" [#nnu].
    \param x_p prescribed "dofval" [#nnp].
    \param b_p (overwritten) prescribed "dofval" [#nnp].
    */
    void reaction_p(
        const xt::xtensor<double, 1>& x_u,
        const xt::xtensor<double, 1>& x_p,
              xt::xtensor<double, 1>& b_p) const;

private:

    Eigen::SparseMatrix<double> m_Auu; ///< The matrix.
    Eigen::SparseMatrix<double> m_Aup; ///< The matrix.
    Eigen::SparseMatrix<double> m_Apu; ///< The matrix.
    Eigen::SparseMatrix<double> m_App; ///< The matrix.
    std::vector<Eigen::Triplet<double>> m_Tuu; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tup; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tpu; ///< Matrix entries.
    std::vector<Eigen::Triplet<double>> m_Tpp; ///< Matrix entries.
    xt::xtensor<size_t, 1> m_iiu; ///< See iiu()
    xt::xtensor<size_t, 1> m_iip; ///< See iip()
    size_t m_nnu; ///< See #nnu
    size_t m_nnp; ///< See #nnp

    /**
    Renumbered DOFs per node, such that

        iiu = arange(nnu)
        iip = nnu + arange(nnp)

    making is much simpler to slice.
    */
    xt::xtensor<size_t, 2> m_part;

    // grant access to solver class
    template <class> friend class MatrixPartitionedSolver;

private:

    // Convert arrays (Eigen version of VectorPartitioned, which contains public functions)
    Eigen::VectorXd AsDofs_u(const xt::xtensor<double, 1>& dofval) const;
    Eigen::VectorXd AsDofs_u(const xt::xtensor<double, 2>& nodevec) const;
    Eigen::VectorXd AsDofs_p(const xt::xtensor<double, 1>& dofval) const;
    Eigen::VectorXd AsDofs_p(const xt::xtensor<double, 2>& nodevec) const;
};

/**
Solver for MatrixPartitioned().
The idea is that this solver class can be used to solve for multiple right-hand-sides
using one factorisation.
*/
template <class Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
class MatrixPartitionedSolver {
public:

    MatrixPartitionedSolver() = default;

    /**
    Solve \f$ x_u = A_{uu}^{-1} (b_u - A_{up} * x_p) \f$.

    \param A sparse matrix, see MatrixPartitioned().
    \param b nodevec [nelem, ndim].
    \param x nodevec [nelem, ndim], used to read \f$ x_p \f$.
    \return x nodevec [nelem, ndim], \f$ x_u \f$ filled, \f$ x_p \f$ copied.
    */
    xt::xtensor<double, 2> Solve(
        MatrixPartitioned& A,
        const xt::xtensor<double, 2>& b,
        const xt::xtensor<double, 2>& x);

    /**
    Same as Solve(MatrixPartitioned&, const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but filling \f$ x_u \f$ in place.

    \param A sparse matrix, see MatrixPartitioned().
    \param b nodevec [nelem, ndim].
    \param x nodevec [nelem, ndim], \f$ x_p \f$ read, \f$ x_u \f$ filled.
    */
    void solve(MatrixPartitioned& A, const xt::xtensor<double, 2>& b, xt::xtensor<double, 2>& x);

    /**
    Same as Solve(MatrixPartitioned&, const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but for "dofval" input and output.

    \param A sparse matrix, see MatrixPartitioned().
    \param b dofval [ndof].
    \param x dofval [ndof], used to read \f$ x_p \f$.
    \return x dofval [ndof], \f$ x_u \f$ filled, \f$ x_p \f$ copied.
    */
    xt::xtensor<double, 1> Solve(
        MatrixPartitioned& A,
        const xt::xtensor<double, 1>& b,
        const xt::xtensor<double, 1>& x);

    /**
    Same as Solve(MatrixPartitioned&, const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but filling \f$ x_u \f$ in place.

    \param A sparse matrix, see MatrixPartitioned().
    \param b dofval [ndof].
    \param x dofval [ndof], \f$ x_p \f$ read, \f$ x_u \f$ filled.
    */
    void solve(MatrixPartitioned& A, const xt::xtensor<double, 1>& b, xt::xtensor<double, 1>& x);

    /**
    Same as Solve(MatrixPartitioned&, const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but with partitioned input and output.

    \param A sparse matrix, see MatrixPartitioned().
    \param b_u unknown dofval [nnu].
    \param x_p prescribed dofval [nnp]
    \return x_u unknown dofval [nnu].
    */
    xt::xtensor<double, 1> Solve_u(
        MatrixPartitioned& A,
        const xt::xtensor<double, 1>& b_u,
        const xt::xtensor<double, 1>& x_p);

    /**
    Same as
    Solve_u(MatrixPartitioned&, const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but writing to pre-allocated output.

    \param A sparse matrix, see MatrixPartitioned().
    \param b_u unknown dofval [nnu].
    \param x_p prescribed dofval [nnp]
    \param x_u (overwritten) unknown dofval [nnu].
    */
    void solve_u(
        MatrixPartitioned& A,
        const xt::xtensor<double, 1>& b_u,
        const xt::xtensor<double, 1>& x_p,
        xt::xtensor<double, 1>& x_u);

private:
    Solver m_solver; ///< solver
    bool m_factor = true; ///< signal to force factorization
    void factorize(MatrixPartitioned& matrix); ///< compute inverse (evaluated by "solve")
};

} // namespace GooseFEM

#include "MatrixPartitioned.hpp"

#endif

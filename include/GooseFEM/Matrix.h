/**
Sparse matrix.

\file Matrix.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIX_H
#define GOOSEFEM_MATRIX_H

#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace GooseFEM {

// forward declaration
template <class>
class MatrixSolver;

/**
Sparse matrix.

See Vector() for bookkeeping definitions.
*/
class Matrix {
public:
    Matrix() = default;

    virtual ~Matrix() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    */
    Matrix(const array_type::tensor<size_t, 2>& conn, const array_type::tensor<size_t, 2>& dofs);

    /**
    \return Number of elements.
    */
    size_t nelem() const;

    /**
    \return Number of nodes per element.
    */
    size_t nne() const;

    /**
    \return Number of nodes.
    */
    size_t nnode() const;

    /**
    \return Number of dimensions.
    */
    size_t ndim() const;

    /**
    \return Number of DOFs.
    */
    size_t ndof() const;

    /**
    \return DOFs per node [#nnode, #ndim]
    */
    array_type::tensor<size_t, 2> dofs() const;

    /**
    Assemble from matrices stored per element.

    \param elemmat [#nelem, #nne * #ndim, #nne * #ndim].
    */
    virtual void assemble(const array_type::tensor<double, 3>& elemmat);

    /**
    Overwrite matrix with (sparse) matrix.

    \param rows Row numbers in the matrix [n].
    \param cols Column numbers in the matrix [n].
    \param matrix Data entries on (rows, cols) [n].
    */
    virtual void
    set(const array_type::tensor<size_t, 1>& rows,
        const array_type::tensor<size_t, 1>& cols,
        const array_type::tensor<double, 2>& matrix);

    /**
    Add a (sparse) matrix to the matrix.

    \param rows Row numbers in the matrix [n].
    \param cols Column numbers in the matrix [n].
    \param matrix Data entries on (rows, cols) [n].
    */
    virtual void
    add(const array_type::tensor<size_t, 1>& rows,
        const array_type::tensor<size_t, 1>& cols,
        const array_type::tensor<double, 2>& matrix);

    /**
    \return Matrix as dense matrix [#ndof, #ndof].
    */
    array_type::tensor<double, 2> Todense() const;

    /**
    Covert matrix to dense matrix.

    \param ret overwritten [#ndof, #ndof].
    */
    virtual void todense(array_type::tensor<double, 2>& ret) const;

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x nodevec [#nelem, #ndim].
    \return b nodevec overwritten [#nelem, #ndim].
    */
    array_type::tensor<double, 2> Dot(const array_type::tensor<double, 2>& x) const;

    /**
    Same as Dot(const array_type::tensor<double, 2>&, array_type::tensor<double, 2>& b)
    but writing to preallocated data.

    \param x nodevec [#nelem, #ndim].
    \param b nodevec overwritten [#nelem, #ndim].
    */
    virtual void
    dot(const array_type::tensor<double, 2>& x, array_type::tensor<double, 2>& b) const;

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x dofval [#ndof].
    \return b dofval overwritten [#ndof].
    */
    array_type::tensor<double, 1> Dot(const array_type::tensor<double, 1>& x) const;

    /**
    Same as Dot(const array_type::tensor<double, 1>&, array_type::tensor<double, 1>& b)
    but writing to preallocated data.

    \param x dofval [#ndof].
    \param b dofval overwritten [#ndof].
    */
    virtual void
    dot(const array_type::tensor<double, 1>& x, array_type::tensor<double, 1>& b) const;

protected:
    array_type::tensor<size_t, 2> m_conn; ///< Connectivity [#nelem, #nne].
    array_type::tensor<size_t, 2> m_dofs; ///< DOF-numbers per node [#nnode, #ndim].
    size_t m_nelem; ///< See nelem().
    size_t m_nne; ///< See nne().
    size_t m_nnode; ///< See nnode().
    size_t m_ndim; ///< See ndim().
    size_t m_ndof; ///< See ndof().
    bool m_changed = true; ///< Signal changes to data.

private:
    Eigen::SparseMatrix<double> m_A; ///< The matrix.
    std::vector<Eigen::Triplet<double>> m_T; ///< Matrix entries.
    template <class>
    friend class MatrixSolver; ///< Grant access to solver class

    /**
    Convert "nodevec" to "dofval" (overwrite entries that occur more than once).

    \param nodevec input [#nnode, #ndim]
    \return dofval output [#ndof]
    */
    Eigen::VectorXd AsDofs(const array_type::tensor<double, 2>& nodevec) const;

    /**
    Convert "dofval" to "nodevec" (overwrite entries that occur more than once).

    \param dofval input [#ndof]
    \param nodevec output [#nnode, #ndim]
    */
    void asNode(const Eigen::VectorXd& dofval, array_type::tensor<double, 2>& nodevec) const;
};

/**
Solver for Matrix().
The idea is that this solver class can be used to solve for multiple right-hand-sides
using one factorisation.
*/
template <class Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
class MatrixSolver {
public:
    MatrixSolver() = default;

    /**
    Solve \f$ x = A^{-1} b \f$.

    \param A sparse matrix, see Matrix().
    \param b nodevec [nelem, ndim].
    \return x nodevec [nelem, ndim].
    */
    array_type::tensor<double, 2> Solve(Matrix& A, const array_type::tensor<double, 2>& b);

    /**
    Same as Solve(Matrix&, const array_type::tensor<double, 2>&)
    but writing to preallocated data.

    \param A sparse matrix, see Matrix().
    \param b nodevec [nelem, ndim].
    \param x nodevec overwritten [nelem, ndim].
    */
    void solve(Matrix& A, const array_type::tensor<double, 2>& b, array_type::tensor<double, 2>& x);

    /**
    Same as Solve(Matrix&, const array_type::tensor<double, 2>&)
    but for "dofval" input and output.

    \param A sparse matrix, see Matrix().
    \param b dofval [ndof].
    \return x dofval [ndof].
    */
    array_type::tensor<double, 1> Solve(Matrix& A, const array_type::tensor<double, 1>& b);

    /**
    Same as Solve(Matrix&, const array_type::tensor<double, 1>&)
    but writing to preallocated data.

    \param A sparse matrix, see Matrix().
    \param b dofval [ndof].
    \param x dofval overwritten [ndof].
    */
    void solve(Matrix& A, const array_type::tensor<double, 1>& b, array_type::tensor<double, 1>& x);

private:
    Solver m_solver; ///< Solver.
    bool m_factor = true; ///< Signal to force factorization.
    void factorize(Matrix& A); ///< Compute inverse (evaluated by "solve").
};

} // namespace GooseFEM

#include "Matrix.hpp"

#endif

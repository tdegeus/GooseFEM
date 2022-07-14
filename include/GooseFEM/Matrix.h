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
    Matrix(const array_type::tensor<size_t, 2>& conn, const array_type::tensor<size_t, 2>& dofs)
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

    /**
    \return Number of elements.
    */
    size_t nelem() const
    {
        return m_nelem;
    }

    /**
    \return Number of nodes per element.
    */
    size_t nne() const
    {
        return m_nne;
    }

    /**
    \return Number of nodes.
    */
    size_t nnode() const
    {
        return m_nnode;
    }

    /**
    \return Number of dimensions.
    */
    size_t ndim() const
    {
        return m_ndim;
    }

    /**
    \return Number of DOFs.
    */
    size_t ndof() const
    {
        return m_ndof;
    }

    /**
    \return DOFs per node [#nnode, #ndim]
    */
    array_type::tensor<size_t, 2> dofs() const
    {
        return m_dofs;
    }

    /**
    Assemble from matrices stored per element.

    \param elemmat [#nelem, #nne * #ndim, #nne * #ndim].
    */
    virtual void assemble(const array_type::tensor<double, 3>& elemmat)
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

    /**
    Overwrite matrix with (sparse) matrix.

    \param rows Row numbers in the matrix [n].
    \param cols Column numbers in the matrix [n].
    \param matrix Data entries on (rows, cols) [n].
    */
    virtual void
    set(const array_type::tensor<size_t, 1>& rows,
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

    /**
    Add a (sparse) matrix to the matrix.

    \param rows Row numbers in the matrix [n].
    \param cols Column numbers in the matrix [n].
    \param matrix Data entries on (rows, cols) [n].
    */
    virtual void
    add(const array_type::tensor<size_t, 1>& rows,
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

    /**
    \return Matrix as dense matrix [#ndof, #ndof].
    */
    array_type::tensor<double, 2> Todense() const
    {
        array_type::tensor<double, 2> ret = xt::empty<double>({m_ndof, m_ndof});
        this->todense(ret);
        return ret;
    }

    /**
    Covert matrix to dense matrix.

    \param ret overwritten [#ndof, #ndof].
    */
    virtual void todense(array_type::tensor<double, 2>& ret) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(ret, {m_ndof, m_ndof}));

        ret.fill(0.0);

        for (int k = 0; k < m_A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_A, k); it; ++it) {
                ret(it.row(), it.col()) = it.value();
            }
        }
    }

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x nodevec [#nelem, #ndim].
    \return b nodevec overwritten [#nelem, #ndim].
    */
    array_type::tensor<double, 2> Dot(const array_type::tensor<double, 2>& x) const
    {
        array_type::tensor<double, 2> b = xt::empty<double>({m_nnode, m_ndim});
        this->dot(x, b);
        return b;
    }

    /**
    Same as Dot(const array_type::tensor<double, 2>&, array_type::tensor<double, 2>& b)
    but writing to preallocated data.

    \param x nodevec [#nelem, #ndim].
    \param b nodevec overwritten [#nelem, #ndim].
    */
    virtual void dot(const array_type::tensor<double, 2>& x, array_type::tensor<double, 2>& b) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
        this->asNode(m_A * this->AsDofs(x), b);
    }

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x dofval [#ndof].
    \return b dofval overwritten [#ndof].
    */
    array_type::tensor<double, 1> Dot(const array_type::tensor<double, 1>& x) const
    {
        array_type::tensor<double, 1> b = xt::empty<double>({m_ndof});
        this->dot(x, b);
        return b;
    }

    /**
    Same as Dot(const array_type::tensor<double, 1>&, array_type::tensor<double, 1>& b)
    but writing to preallocated data.

    \param x dofval [#ndof].
    \param b dofval overwritten [#ndof].
    */
    virtual void dot(const array_type::tensor<double, 1>& x, array_type::tensor<double, 1>& b) const
    {
        GOOSEFEM_ASSERT(b.size() == m_ndof);
        GOOSEFEM_ASSERT(x.size() == m_ndof);

        Eigen::Map<Eigen::VectorXd>(b.data(), b.size()).noalias() =
            m_A * Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
    }

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
    Eigen::VectorXd AsDofs(const array_type::tensor<double, 2>& nodevec) const
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

    /**
    Convert "dofval" to "nodevec" (overwrite entries that occur more than once).

    \param dofval input [#ndof]
    \param nodevec output [#nnode, #ndim]
    */
    void asNode(const Eigen::VectorXd& dofval, array_type::tensor<double, 2>& nodevec) const
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
    array_type::tensor<double, 2> Solve(Matrix& A, const array_type::tensor<double, 2>& b)
    {
        array_type::tensor<double, 2> x = xt::empty<double>({A.m_nnode, A.m_ndim});
        this->solve(A, b, x);
        return x;
    }

    /**
    Same as Solve(Matrix&, const array_type::tensor<double, 2>&)
    but writing to preallocated data.

    \param A sparse matrix, see Matrix().
    \param b nodevec [nelem, ndim].
    \param x nodevec overwritten [nelem, ndim].
    */
    void solve(Matrix& A, const array_type::tensor<double, 2>& b, array_type::tensor<double, 2>& x)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b, {A.m_nnode, A.m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(x, {A.m_nnode, A.m_ndim}));
        this->factorize(A);
        Eigen::VectorXd X = m_solver.solve(A.AsDofs(b));
        A.asNode(X, x);
    }

    /**
    Same as Solve(Matrix&, const array_type::tensor<double, 2>&)
    but for "dofval" input and output.

    \param A sparse matrix, see Matrix().
    \param b dofval [ndof].
    \return x dofval [ndof].
    */
    array_type::tensor<double, 1> Solve(Matrix& A, const array_type::tensor<double, 1>& b)
    {
        array_type::tensor<double, 1> x = xt::empty<double>({A.m_ndof});
        this->solve(A, b, x);
        return x;
    }

    /**
    Same as Solve(Matrix&, const array_type::tensor<double, 1>&)
    but writing to preallocated data.

    \param A sparse matrix, see Matrix().
    \param b dofval [ndof].
    \param x dofval overwritten [ndof].
    */
    void solve(Matrix& A, const array_type::tensor<double, 1>& b, array_type::tensor<double, 1>& x)
    {
        GOOSEFEM_ASSERT(b.size() == A.m_ndof);
        GOOSEFEM_ASSERT(x.size() == A.m_ndof);
        this->factorize(A);
        Eigen::Map<Eigen::VectorXd>(x.data(), x.size()).noalias() =
            m_solver.solve(Eigen::Map<const Eigen::VectorXd>(b.data(), A.m_ndof));
    }

private:
    Solver m_solver; ///< Solver.
    bool m_factor = true; ///< Signal to force factorization.

    /**
    Compute inverse (evaluated by "solve").
    */
    void factorize(Matrix& A)
    {
        if (!A.m_changed && !m_factor) {
            return;
        }
        m_solver.compute(A.m_A);
        m_factor = false;
        A.m_changed = false;
    }
};

} // namespace GooseFEM

#endif

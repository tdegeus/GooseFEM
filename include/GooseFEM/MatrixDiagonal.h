/**
Diagonal matrix.

\file MatrixDiagonal.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXDIAGONAL_H
#define GOOSEFEM_MATRIXDIAGONAL_H

#include "config.h"
#include "Element.h"

namespace GooseFEM {

/**
Diagonal matrix.

See Vector() for bookkeeping definitions.
*/
class MatrixDiagonal {
public:

    MatrixDiagonal() = default;

    virtual ~MatrixDiagonal() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    */
    MatrixDiagonal(const xt::xtensor<size_t, 2>& conn, const xt::xtensor<size_t, 2>& dofs);

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
    xt::xtensor<size_t, 2> dofs() const;

    /**
    Assemble from matrices stored per element.
    \warning Ignores any off-diagonal terms.

    \param elemmat [#nelem, #nne * #ndim, #nne * #ndim].
    */
    virtual void assemble(const xt::xtensor<double, 3>& elemmat);

    /**
    Set all (diagonal) matrix components.

    \param A The matrix [#ndof].
    */
    virtual void set(const xt::xtensor<double, 1>& A);

    /**
    Return matrix as diagonal matrix.

    \param [#ndof].
    */
    virtual xt::xtensor<double, 1> Todiagonal() const;

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x nodevec [#nelem, #ndim].
    \return b nodevec overwritten [#nelem, #ndim].
    */
    xt::xtensor<double, 2> Dot(const xt::xtensor<double, 2>& x) const;

    /**
    Same as Dot(const xt::xtensor<double, 2>&, xt::xtensor<double, 2>& b)
    but writing to preallocated data.

    \param x nodevec [#nelem, #ndim].
    \param b nodevec overwritten [#nelem, #ndim].
    */
    virtual void dot(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const;

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x dofval [#ndof].
    \return b dofval overwritten [#ndof].
    */
    xt::xtensor<double, 1> Dot(const xt::xtensor<double, 1>& x) const;

    /**
    Same as Dot(const xt::xtensor<double, 1>&, xt::xtensor<double, 1>& b)
    but writing to preallocated data.

    \param x dofval [#ndof].
    \param b dofval overwritten [#ndof].
    */
    virtual void dot(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const;

    /**
    Solve \f$ x = A^{-1} b \f$.

    \param b nodevec [nelem, ndim].
    \return x nodevec [nelem, ndim].
    */
    xt::xtensor<double, 2> Solve(const xt::xtensor<double, 2>& b);

    /**
    Same as Solve(const xt::xtensor<double, 2>&)
    but writing to preallocated data.

    \param b nodevec [nelem, ndim].
    \param x nodevec overwritten [nelem, ndim].
    */
    virtual void solve(const xt::xtensor<double, 2>& b, xt::xtensor<double, 2>& x);

    /**
    Same as Solve(const xt::xtensor<double, 2>&)
    but for "dofval" input and output.

    \param b dofval [ndof].
    \return x dofval [ndof].
    */
    xt::xtensor<double, 1> Solve(const xt::xtensor<double, 1>& b);

    /**
    Same as Solve(const xt::xtensor<double, 1>&)
    but writing to preallocated data.

    \param b dofval [ndof].
    \param x dofval overwritten [ndof].
    */
    virtual void solve(const xt::xtensor<double, 1>& b, xt::xtensor<double, 1>& x);

protected:

    xt::xtensor<size_t, 2> m_conn; ///< Connectivity [#nelem, #nne].
    xt::xtensor<size_t, 2> m_dofs; ///< DOF-numbers per node [#nnode, #ndim].
    size_t m_nelem; ///< See nelem().
    size_t m_nne; ///< See nne().
    size_t m_nnode; ///< See nnode().
    size_t m_ndim; ///< See ndim().
    size_t m_ndof; ///< See ndof().
    bool m_changed = true; ///< Signal changes to data.

private:
    xt::xtensor<double, 1> m_A; ///< The matrix.
    xt::xtensor<double, 1> m_inv; /// Inverse of the matrix.
    void factorize(); ///< Compute inverse (automatically evaluated by "solve").
};

} // namespace GooseFEM

#include "MatrixDiagonal.hpp"

#endif

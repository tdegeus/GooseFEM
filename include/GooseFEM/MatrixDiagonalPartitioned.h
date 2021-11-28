/**
Diagonal matrix that is partitioned in:
-   unknown DOFs
-   prescribed DOFs

\file MatrixDiagonalPartitioned.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXDIAGONALPARTITIONED_H
#define GOOSEFEM_MATRIXDIAGONALPARTITIONED_H

#include "MatrixDiagonal.h"
#include "config.h"

namespace GooseFEM {

/**
Diagonal and partitioned matrix.

See Vector() for bookkeeping definitions.
*/
class MatrixDiagonalPartitioned : public MatrixDiagonal {
public:
    MatrixDiagonalPartitioned() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    \param iip prescribed DOFs [#nnp].
    */
    MatrixDiagonalPartitioned(
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

    void set(const xt::xtensor<double, 1>& A) override;

    xt::xtensor<double, 1> Todiagonal() const override;

    void dot(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const override;
    void dot(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const override;

    /**
    \param x_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \return b_u dofval [#nnu].
    */
    xt::xtensor<double, 1>
    Dot_u(const xt::xtensor<double, 1>& x_u, const xt::xtensor<double, 1>& x_p) const;

    /**
    \param x_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \param b_u (overwritten) dofval [#nnu].
    */
    void dot_u(
        const xt::xtensor<double, 1>& x_u,
        const xt::xtensor<double, 1>& x_p,
        xt::xtensor<double, 1>& b_u) const;

    /**
    \param x_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \return b_p dofval [#nnp].
    */
    xt::xtensor<double, 1>
    Dot_p(const xt::xtensor<double, 1>& x_u, const xt::xtensor<double, 1>& x_p) const;

    /**
    \param x_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \param b_p (overwritten) dofval [#nnp].
    */
    void dot_p(
        const xt::xtensor<double, 1>& x_u,
        const xt::xtensor<double, 1>& x_p,
        xt::xtensor<double, 1>& b_p) const;

    /**
    Solve \f$ x_u = A_{uu}^{-1} (b_u - A_{up} * x_p) \equiv A_{uu}^{-1} b_u \f$.

    \param b nodevec [#nelem, #ndim].
    \param x nodevec, modified with `x_u` [#nelem, #ndim].
    */
    void solve(const xt::xtensor<double, 2>& b, xt::xtensor<double, 2>& x) override;

    /**
    Solve \f$ x_u = A_{uu}^{-1} (b_u - A_{up} * x_p) \equiv A_{uu}^{-1} b_u \f$.

    \param b dofval [#ndof].
    \param x dofval, modified with `x_u` [#ndof].
    */
    void solve(const xt::xtensor<double, 1>& b, xt::xtensor<double, 1>& x) override;

    /**
    \param b_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \return x_u dofval [#nnu].
    */
    xt::xtensor<double, 1>
    Solve_u(const xt::xtensor<double, 1>& b_u, const xt::xtensor<double, 1>& x_p);

    /**
    \param b_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \param x_u (overwritten) dofval [#nnu].
    */
    void solve_u(
        const xt::xtensor<double, 1>& b_u,
        const xt::xtensor<double, 1>& x_p,
        xt::xtensor<double, 1>& x_u);

    /**
    Get right-hand-size for corresponding to the prescribed DOFs.

    \f$ b_p = A_{pu} * x_u + A_{pp} * x_p = A_{pp} * x_p \equiv A_{pp} * x_p \f$

    and assemble them to the appropriate places in "nodevec".

    \param x "nodevec" [#nnode, #ndim].
    \param b "nodevec" [#nnode, #ndim].
    \return Copy of `b` with \f$ b_p \f$ overwritten.
    */
    xt::xtensor<double, 2>
    Reaction(const xt::xtensor<double, 2>& x, const xt::xtensor<double, 2>& b) const;

    /**
    Same as Reaction(const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but inserting in-place.

    \param x "nodevec" [#nnode, #ndim].
    \param b "nodevec" [#nnode, #ndim], \f$ b_p \f$ overwritten.
    */
    void reaction(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const;

    /**
    Same as Reaction(const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&),
    but of "dofval" input and output.

    \param x "dofval" [#ndof].
    \param b "dofval" [#ndof].
    \return Copy of `b` with \f$ b_p \f$ overwritten.
    */
    xt::xtensor<double, 1>
    Reaction(const xt::xtensor<double, 1>& x, const xt::xtensor<double, 1>& b) const;

    /**
    Same as Reaction(const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but inserting in-place.

    \param x "dofval" [#ndof].
    \param b "dofval" [#ndof], \f$ b_p \f$ overwritten.
    */
    void reaction(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b)
        const; // modified with "b_p"

    /**
    Same as Reaction(const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&),
    but with partitioned input and output.

    \param x_u unknown "dofval" [#nnu].
    \param x_p prescribed "dofval" [#nnp].
    \return b_p prescribed "dofval" [#nnp].
    */
    xt::xtensor<double, 1>
    Reaction_p(const xt::xtensor<double, 1>& x_u, const xt::xtensor<double, 1>& x_p) const;

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
    // The diagonal matrix, and its inverse (re-used to solve different RHS)
    xt::xtensor<double, 1> m_Auu;
    xt::xtensor<double, 1> m_App;
    xt::xtensor<double, 1> m_inv_uu;

    // Bookkeeping
    xt::xtensor<size_t, 2> m_part; // DOF-numbers per node, renumbered  [nnode, ndim]
    xt::xtensor<size_t, 1> m_iiu; // DOF-numbers that are unknown      [nnu]
    xt::xtensor<size_t, 1> m_iip; // DOF-numbers that are prescribed   [nnp]

    // Dimensions
    size_t m_nnu; // number of unknown DOFs
    size_t m_nnp; // number of prescribed DOFs

    // Compute inverse (automatically evaluated by "solve")
    void factorize();
};

} // namespace GooseFEM

#include "MatrixDiagonalPartitioned.hpp"

#endif

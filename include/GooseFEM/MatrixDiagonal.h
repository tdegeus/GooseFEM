/**
Diagonal matrix.

\file Matrix.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXDIAGONAL_H
#define GOOSEFEM_MATRIXDIAGONAL_H

#include "config.h"

namespace GooseFEM {

class MatrixDiagonal {
public:
    // Constructors
    MatrixDiagonal() = default;
    MatrixDiagonal(const xt::xtensor<size_t, 2>& conn, const xt::xtensor<size_t, 2>& dofs);

    // Dimensions
    size_t nelem() const; // number of elements
    size_t nne() const;   // number of nodes per element
    size_t nnode() const; // number of nodes
    size_t ndim() const;  // number of dimensions
    size_t ndof() const;  // number of DOFs

    // DOF lists
    xt::xtensor<size_t, 2> dofs() const; // DOFs

    // Set matrix components
    void set(const xt::xtensor<double, 1>& A);

    // assemble from matrices stored per element [nelem, nne*ndim, nne*ndim]
    // WARNING: ignores any off-diagonal terms
    void assemble(const xt::xtensor<double, 3>& elemmat);

    // Dot-product:
    // b_i = A_ij * x_j
    void dot(const xt::xtensor<double, 2>& x, xt::xtensor<double, 2>& b) const;
    void dot(const xt::xtensor<double, 1>& x, xt::xtensor<double, 1>& b) const;

    // Solve:
    // x = A \ b
    void solve(const xt::xtensor<double, 2>& b, xt::xtensor<double, 2>& x);
    void solve(const xt::xtensor<double, 1>& b, xt::xtensor<double, 1>& x);

    // Return matrix as diagonal matrix (column)
    xt::xtensor<double, 1> Todiagonal() const;

    // Auto-allocation of the functions above
    xt::xtensor<double, 2> Dot(const xt::xtensor<double, 2>& x) const;
    xt::xtensor<double, 1> Dot(const xt::xtensor<double, 1>& x) const;
    xt::xtensor<double, 2> Solve(const xt::xtensor<double, 2>& b);
    xt::xtensor<double, 1> Solve(const xt::xtensor<double, 1>& b);

private:
    // The diagonal matrix, and its inverse (re-used to solve different RHS)
    xt::xtensor<double, 1> m_A;
    xt::xtensor<double, 1> m_inv;

    // Signal changes to data compare to the last inverse
    bool m_factor = true;

    // Bookkeeping
    xt::xtensor<size_t, 2> m_conn; // connectivity [nelem, nne]
    xt::xtensor<size_t, 2> m_dofs; // DOF-numbers per node [nnode, ndim]

    // Dimensions
    size_t m_nelem; // number of elements
    size_t m_nne;   // number of nodes per element
    size_t m_nnode; // number of nodes
    size_t m_ndim;  // number of dimensions
    size_t m_ndof;  // number of DOFs

    // Compute inverse (automatically evaluated by "solve")
    void factorize();
};

} // namespace GooseFEM

#include "MatrixDiagonal.hpp"

#endif

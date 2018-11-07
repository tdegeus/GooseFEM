/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_MATRIXDIAGONALPARTITIONED_H
#define XGOOSEFEM_MATRIXDIAGONALPARTITIONED_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// =========================================== GooseFEM ============================================

namespace xGooseFEM {

// -------------------------------------------------------------------------------------------------

class MatrixDiagonalPartitioned
{
public:

  // default constructor
  MatrixDiagonalPartitioned() = default;

  // constructor
  MatrixDiagonalPartitioned(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs,
    const xt::xtensor<size_t,1> &iip);

  // dimensions
  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t nnode() const; // number of nodes
  size_t ndim()  const; // number of dimensions
  size_t ndof()  const; // number of DOFs
  size_t nnu()   const; // number of unknown DOFs
  size_t nnp()   const; // number of prescribed DOFs

  // DOF lists
  xt::xtensor<size_t,2> dofs() const; // DOFs
  xt::xtensor<size_t,1> iiu()  const; // unknown    DOFs
  xt::xtensor<size_t,1> iip()  const; // prescribed DOFs

  // product: b_i = A_ij * x_j
  //   b   = A    * x
  xt::xtensor<double,1> dot  (const xt::xtensor<double,1> &x                                    ) const;
  //   b_u = A_uu * x_u + A_up * x_p
  xt::xtensor<double,1> dot_u(const xt::xtensor<double,1> &x_u, const xt::xtensor<double,1> &x_p) const;
  //   b_p = A_pu * x_u + A_pp * x_p
  xt::xtensor<double,1> dot_p(const xt::xtensor<double,1> &x_u, const xt::xtensor<double,1> &x_p) const;

  // assemble from matrices stored per element [nelem, nne*ndim, nne*ndim]
  // WARNING: ignores any off-diagonal terms
  void assemble(const xt::xtensor<double,3> &elemmat);

  // solve: x_u = A_uu \ ( b_u - A_up * x_p )
  xt::xtensor<double,1> solve(const xt::xtensor<double,1> &b_u, const xt::xtensor<double,1> &x_p);

  // return matrix as diagonal matrix
  xt::xtensor<double,1> asDiagonal() const;

private:

  // the diagonal matrix, and its inverse (re-used to solve different RHS)
  xt::xtensor<double,1> m_data_uu;
  xt::xtensor<double,1> m_data_pp;
  xt::xtensor<double,1> m_inv_uu;

  // signal changes to data compare to the last inverse
  bool m_change=false;

  // bookkeeping
  xt::xtensor<size_t,2> m_conn; // connectivity                      [nelem, nne ]
  xt::xtensor<size_t,2> m_dofs; // DOF-numbers per node              [nnode, ndim]
  xt::xtensor<size_t,2> m_part; // DOF-numbers per node, renumbered  [nnode, ndim]
  xt::xtensor<size_t,1> m_iiu;  // DOF-numbers that are unknown      [nnu]
  xt::xtensor<size_t,1> m_iip;  // DOF-numbers that are prescribed   [nnp]

  // dimensions
  size_t m_nelem; // number of elements
  size_t m_nne;   // number of nodes per element
  size_t m_nnode; // number of nodes
  size_t m_ndim;  // number of dimensions
  size_t m_ndof;  // number of DOFs
  size_t m_nnu;   // number of unknown DOFs
  size_t m_nnp;   // number of prescribed DOFs

  // compute inverse (automatically evaluated by "solve")
  void factorize();
};

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

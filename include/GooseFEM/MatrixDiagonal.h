/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIXDIAGONAL_H
#define GOOSEFEM_MATRIXDIAGONAL_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// =========================================== GooseFEM ============================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

class MatrixDiagonal
{
private:

  // data
  xt::xtensor<double,1> m_data;          // the diagonal matrix (not-partitioned)
  xt::xtensor<double,1> m_inv;           // inverse of "m_data", can be re-used to solve different right-hand-sides
  bool                  m_change=false;  // signal changes to data compare to the last inverse

  // information
  xt::xtensor<size_t,2> m_conn; // connectivity                               [nelem, nne ]
  xt::xtensor<size_t,2> m_dofs; // DOF-numbers per node                       [nnode, ndim]
  xt::xtensor<size_t,2> m_part; // DOF-numbers per node, after partitioning   [nnode, ndim]
  xt::xtensor<size_t,1> m_iiu;  // DOF-numbers that are unknown               [nnu]
  xt::xtensor<size_t,1> m_iip;  // DOF-numbers that are prescribed            [nnp]

  // dimensions
  size_t m_nelem; // number of elements
  size_t m_nne;   // number of nodes per element
  size_t m_nnode; // number of nodes
  size_t m_ndim;  // number of dimensions
  size_t m_ndof;  // number of DOFs
  size_t m_nnu;   // number of unknown DOFs
  size_t m_nnp;   // number of prescribed DOFs

public:

  // constructor
  MatrixDiagonal() = default;
  MatrixDiagonal(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs);
  MatrixDiagonal(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs, const xt::xtensor<size_t,1> &iip);

  // index operators: access plain storage
  double&       operator[](size_t i);
  const double& operator[](size_t i) const;

  // index operators: access using matrix indices
  double&       operator()(size_t a);
  const double& operator()(size_t a) const;
  double&       operator()(size_t a, size_t b);
  const double& operator()(size_t a, size_t b) const;

  // dimensions
  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t nnode() const; // number of nodes
  size_t ndim()  const; // number of dimensions
  size_t ndof()  const; // number of DOFs
  size_t nnu()   const; // number of unknown DOFs
  size_t nnp()   const; // number of prescribed DOFs

  // DOF lists
  xt::xtensor<size_t,1> iiu() const; // unknown    DOFs
  xt::xtensor<size_t,1> iip() const; // prescribed DOFs

  // product: c_i = A_ij * b_j
  xt::xtensor<double,1> dot  (const xt::xtensor<double,1> &b                                    ) const; // c   = A    * b
  xt::xtensor<double,1> dot_u(const xt::xtensor<double,1> &b                                    ) const; // c_u = A_uu * b_u + A_up * b_p
  xt::xtensor<double,1> dot_u(const xt::xtensor<double,1> &b_u, const xt::xtensor<double,1> &b_p) const; // c_u = A_uu * b_u + A_up * b_p
  xt::xtensor<double,1> dot_p(const xt::xtensor<double,1> &b                                    ) const; // c_p = A_pu * b_u + A_pp * b_p
  xt::xtensor<double,1> dot_p(const xt::xtensor<double,1> &b_u, const xt::xtensor<double,1> &b_p) const; // c_p = A_pu * b_u + A_pp * b_p

  // check structure of the matrices stored per element [nelem, nne*ndim, nne*ndim]
  void check_diagonal(const xt::xtensor<double,3> &elemmat) const;

  // assemble from matrices stored per element [nelem, nne*ndim, nne*ndim]
  // WARNING: ignores any off-diagonal terms
  void assemble(const xt::xtensor<double,3> &elemmat);

  // set matrix components from externally assembled object
  void set   (const xt::xtensor<double,1> &matrix   ) const; // diagonal [ndof]
  void set_uu(const xt::xtensor<double,1> &matrix_uu) const; // diagonal [nnu]
  void set_pp(const xt::xtensor<double,1> &matrix_pp) const; // diagonal [nnp]

  // solve
  // ("u_p" is useless as all cross-terms are zero, it is here for aesthetics)
  xt::xtensor<double,1> solve  (const xt::xtensor<double,1> &rhs                                    ); // [ndof] -> [ndof]
  xt::xtensor<double,1> solve  (const xt::xtensor<double,1> &rhs  , const xt::xtensor<double,1> &u_p); // [ndof] -> [ndof]
  xt::xtensor<double,1> solve_u(const xt::xtensor<double,1> &rhs_u                                  ); // [nnu]  -> [nnu]
  xt::xtensor<double,1> solve_u(const xt::xtensor<double,1> &rhs_u, const xt::xtensor<double,1> &u_p); // [nnu]  -> [nnu]

  // get the right-hand-side of prescribed DOFs: [nnu], [nnp] -> [nnp]
  // ("u_u" is useless as all cross-terms are zero, it is here for aesthetics)
  xt::xtensor<double,1> rhs_p(const xt::xtensor<double,1> &u_u, const xt::xtensor<double,1> &u_p) const;

  // return (sub-)matrix as diagonal matrix (column)
  xt::xtensor<double,1> asDiagonal   () const; // [ndpf]
  xt::xtensor<double,1> asDiagonal_uu() const; // [nnu]
  xt::xtensor<double,1> asDiagonal_pp() const; // [nnp]

  // return (sub-)matrix as sparse matrix
  SpMatD asSparse   () const; // [ndof,ndof]
  SpMatD asSparse_uu() const; // [nnu ,nnu ]
  SpMatD asSparse_up() const; // [nnu ,nnp ]
  SpMatD asSparse_pu() const; // [nnp ,nnu ]
  SpMatD asSparse_pp() const; // [nnp ,nnp ]

  // return (sub-)matrix as dense matrix
  xt::xtensor<double,2> asDense   () const; // [ndof,ndof]
  xt::xtensor<double,2> asDense_uu() const; // [nnu ,nnu ]
  xt::xtensor<double,2> asDense_up() const; // [nnu ,nnp ]
  xt::xtensor<double,2> asDense_pu() const; // [nnp ,nnu ]
  xt::xtensor<double,2> asDense_pp() const; // [nnp ,nnp ]

};

// -------------------------------------------------------------------------------------------------

// matrix/column product: '[ndof,ndof]' * [ndof] -> [ndof]
inline xt::xtensor<double,1> operator* (const MatrixDiagonal &A, const xt::xtensor<double,1> &b);

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

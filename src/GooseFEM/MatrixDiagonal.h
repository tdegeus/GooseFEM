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

class DiagonalMatrix
{
private:

  // data
  ColD m_mat; // the diagonal matrix (not-partitioned)
  ColD m_inv; // inverse of "m_mat", can be re-used to solve different right-hand-sides

  // information
  MatS m_conn; // connectivity                               [nelem, nne ]
  MatS m_dofs; // DOF-numbers per node                       [nnode, ndim]
  MatS m_part; // DOF-numbers per node, after partitioning   [nnode, ndim]
  ColS m_iiu;  // DOF-numbers that are unknown               [nnu]
  ColS m_iip;  // DOF-numbers that are prescribed            [nnp]

  // dimensions
  size_t m_nelem; // number of elements
  size_t m_nne;   // number of nodes per element
  size_t m_nnode; // number of nodes
  size_t m_ndim;  // number of dimensions
  size_t m_ndof;  // number of DOFs
  size_t m_nnu;   // number of unknown DOFs
  size_t m_nnp;   // number of prescribed DOFs

  // checks
  bool m_init=false; // "true" only if the matrix has been assembled at least once

public:

  // constructor
  DiagonalMatrix(const MatS &conn, const MatS &dofs, const ColS &iip=ColS());

  // dimensions
  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t nnode() const; // number of nodes
  size_t ndim()  const; // number of dimensions
  size_t ndof()  const; // number of DOFs
  size_t nnu()   const; // number of unknown DOFs
  size_t nnp()   const; // number of prescribed DOFs

  // DOF lists
  ColS iiu() const; // unknown    DOFs
  ColS iip() const; // prescribed DOFs

  // product: c_i = A_ij * b_j
  ColD dot  (const ColD &b                   ) const; // c   = A    * b
  ColD dot_u(const ColD &b                   ) const; // c_u = A_uu * b_u + A_up * b_p
  ColD dot_u(const ColD &b_u, const ColD &b_p) const; // c_u = A_uu * b_u + A_up * b_p
  ColD dot_p(const ColD &b                   ) const; // c_p = A_pu * b_u + A_pp * b_p
  ColD dot_p(const ColD &b_u, const ColD &b_p) const; // c_p = A_pu * b_u + A_pp * b_p

  // check structure of the matrices stored per element [nelem, nne*ndim, nne*ndim]
  void check_diagonal(const ArrD &elemmat) const;

  // assemble from matrices stored per element [nelem, nne*ndim, nne*ndim]
  // WARNING: ignores any off-diagonal terms
  void assemble(const ArrD &elemmat);

  // set matrix components from externally assembled object
  void set   (const ColD &matrix   ) const; // diagonal [ndof]
  void set_uu(const ColD &matrix_uu) const; // diagonal [nnu]
  void set_pp(const ColD &matrix_pp) const; // diagonal [nnp]

  // solve
  // ("u_p" is useless as all cross-terms are zero, it is here for aesthetics)
  ColD solve  (const ColD &rhs, const ColD &u_p=ColS()) const; // [ndof] -> [ndof]
  ColD solve_u(const ColD &rhs, const ColD &u_p=ColS()) const; // [nnu]  -> [nnu]

  // get the right-hand-side of prescribed DOFs: [nnu], [nnp] -> [nnp]
  // ("u_u" is useless as all cross-terms are zero, it is here for aesthetics)
  // TODO write function
  ColD rhs_p(const ColD &u_u, const ColD &u_p) const;

  // return (sub-)matrix as diagonal matrix (column)
  // TODO write function
  ColD asDiagonal   () const; // [ndpf]
  ColD asDiagonal_uu() const; // [nnu]
  ColD asDiagonal_pp() const; // [nnp]

  // return (sub-)matrix as sparse matrix
  // TODO write function
  SpMatD asSparse   () const; // [ndof,ndof]
  SpMatD asSparse_uu() const; // [nnu ,nnu ]
  SpMatD asSparse_up() const; // [nnu ,nnp ]
  SpMatD asSparse_pu() const; // [nnp ,nnu ]
  SpMatD asSparse_pp() const; // [nnp ,nnp ]

  // return (sub-)matrix as dense matrix
  // TODO write function
  MatD asDense   () const; // [ndof,ndof]
  MatD asDense_uu() const; // [nnu ,nnu ]
  MatD asDense_up() const; // [nnu ,nnp ]
  MatD asDense_pu() const; // [nnp ,nnu ]
  MatD asDense_pp() const; // [nnp ,nnp ]

};

// -------------------------------------------------------------------------------------------------

// matrix/column product: '[ndof,ndof]' * [ndof] -> [ndof]
inline ColD operator* (const DiagonalMatrix &A, const ColD &b);

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

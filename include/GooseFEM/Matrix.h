/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIX_H
#define GOOSEFEM_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// =========================================== GooseFEM ============================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

class Matrix
{
public:

  // constructors

  Matrix() = default;
  Matrix(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs);

  // dimensions

  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t nnode() const; // number of nodes
  size_t ndim()  const; // number of dimensions
  size_t ndof()  const; // number of DOFs

  // DOF lists

  xt::xtensor<size_t,2> dofs() const; // DOFs

  // assemble from matrices stored per element [nelem, nne*ndim, nne*ndim]

  void assemble(const xt::xtensor<double,3> &elemmat);

  // solve: x = A \ b
  //   x_u = A_uu \ ( b_u - A_up * x_p )
  //   b_p = A_pu * x_u + A_pp * x_p

  void solve(const xt::xtensor<double,2> &b,
    xt::xtensor<double,2> &x);

  void solve(const xt::xtensor<double,1> &b,
    xt::xtensor<double,1> &x);

  // auto allocation of the functions above

  xt::xtensor<double,2> solve(const xt::xtensor<double,2> &b);

  xt::xtensor<double,1> solve(const xt::xtensor<double,1> &b);

private:

  // the matrix
  Eigen::SparseMatrix<double> m_data;

  // the matrix to assemble
  std::vector<TripD> m_trip;

  // solver (re-used to solve different RHS)
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_solver;

  // signal changes to data compare to the last inverse
  bool m_change=false;

  // bookkeeping
  xt::xtensor<size_t,2> m_conn; // connectivity         [nelem, nne ]
  xt::xtensor<size_t,2> m_dofs; // DOF-numbers per node [nnode, ndim]

  // dimensions
  size_t m_nelem; // number of elements
  size_t m_nne;   // number of nodes per element
  size_t m_nnode; // number of nodes
  size_t m_ndim;  // number of dimensions
  size_t m_ndof;  // number of DOFs

  // compute inverse (automatically evaluated by "solve")

  void factorize();

  // convert arrays (see VectorPartitioned, which contains public functions)

  Eigen::VectorXd asDofs(const xt::xtensor<double,2> &nodevec) const;

  void asNode(const Eigen::VectorXd &dofval, xt::xtensor<double,2> &nodevec) const;

};

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIXPARTITIONED_H
#define GOOSEFEM_MATRIXPARTITIONED_H

// -------------------------------------------------------------------------------------------------

#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

// =================================================================================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

template <class Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
class MatrixPartitioned
{
public:

  // Constructors

  MatrixPartitioned() = default;

  MatrixPartitioned(
    const xt::xtensor<size_t,2> &conn,
    const xt::xtensor<size_t,2> &dofs,
    const xt::xtensor<size_t,1> &iip);

  // Dimensions

  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t nnode() const; // number of nodes
  size_t ndim()  const; // number of dimensions
  size_t ndof()  const; // number of DOFs
  size_t nnu()   const; // number of unknown DOFs
  size_t nnp()   const; // number of prescribed DOFs

  // DOF lists

  xt::xtensor<size_t,2> dofs() const; // DOFs
  xt::xtensor<size_t,1> iiu()  const; // unknown DOFs
  xt::xtensor<size_t,1> iip()  const; // prescribed DOFs

  // Assemble from matrices stored per element [nelem, nne*ndim, nne*ndim]

  void assemble(const xt::xtensor<double,3> &elemmat);

  // Solve:
  // x_u = A_uu \ ( b_u - A_up * x_p )

  void solve(
    const xt::xtensor<double,2> &b,
          xt::xtensor<double,2> &x); // modified with "x_u"

  void solve(
    const xt::xtensor<double,1> &b,
          xt::xtensor<double,1> &x); // modified with "x_u"

  void solve_u(
    const xt::xtensor<double,1> &b_u,
    const xt::xtensor<double,1> &x_p,
          xt::xtensor<double,1> &x_u); // overwritten

  // Get right-hand-size for corresponding to the prescribed DOFs:
  // b_p = A_pu * x_u + A_pp * x_p = A_pp * x_p

  void reaction(
    const xt::xtensor<double,2> &x,
          xt::xtensor<double,2> &b) const; // modified with "b_p"

  void reaction(
    const xt::xtensor<double,1> &x,
          xt::xtensor<double,1> &b) const; // modified with "b_p"

  void reaction_p(
    const xt::xtensor<double,1> &x_u,
    const xt::xtensor<double,1> &x_p,
          xt::xtensor<double,1> &b_p) const; // overwritten

  // Auto-allocation of the functions above

  xt::xtensor<double,2> Solve(
    const xt::xtensor<double,2> &b,
    const xt::xtensor<double,2> &x);

  xt::xtensor<double,1> Solve(
    const xt::xtensor<double,1> &b,
    const xt::xtensor<double,1> &x);

  xt::xtensor<double,1> Solve_u(
    const xt::xtensor<double,1> &b_u,
    const xt::xtensor<double,1> &x_p);

  xt::xtensor<double,2> Reaction(
    const xt::xtensor<double,2> &x,
    const xt::xtensor<double,2> &b) const;

  xt::xtensor<double,1> Reaction(
    const xt::xtensor<double,1> &x,
    const xt::xtensor<double,1> &b) const;

  xt::xtensor<double,1> Reaction_p(
    const xt::xtensor<double,1> &x_u,
    const xt::xtensor<double,1> &x_p) const;

private:

  // The matrix
  Eigen::SparseMatrix<double> m_Auu;
  Eigen::SparseMatrix<double> m_Aup;
  Eigen::SparseMatrix<double> m_Apu;
  Eigen::SparseMatrix<double> m_App;

  // Matrix entries
  std::vector<Eigen::Triplet<double>> m_Tuu;
  std::vector<Eigen::Triplet<double>> m_Tup;
  std::vector<Eigen::Triplet<double>> m_Tpu;
  std::vector<Eigen::Triplet<double>> m_Tpp;

  // Solver (re-used to solve different RHS)
  Solver m_solver;

  // Signal changes to data compare to the last inverse
  bool m_factor=false;

  // Bookkeeping
  xt::xtensor<size_t,2> m_conn; // connectivity                      [nelem, nne ]
  xt::xtensor<size_t,2> m_dofs; // DOF-numbers per node              [nnode, ndim]
  xt::xtensor<size_t,2> m_part; // DOF-numbers per node, renumbered  [nnode, ndim]
  xt::xtensor<size_t,1> m_iiu;  // unknown    DOFs                   [nnu]
  xt::xtensor<size_t,1> m_iip;  // prescribed DOFs                   [nnp]

  // Dimensions
  size_t m_nelem; // number of elements
  size_t m_nne;   // number of nodes per element
  size_t m_nnode; // number of nodes
  size_t m_ndim;  // number of dimensions
  size_t m_ndof;  // number of DOFs
  size_t m_nnu;   // number of unknown DOFs
  size_t m_nnp;   // number of prescribed DOFs

  // Compute inverse (automatically evaluated by "solve")
  void factorize();

  // Convert arrays (Eigen version of VectorPartitioned, which contains public functions)
  Eigen::VectorXd asDofs_u(const xt::xtensor<double,1> &dofval ) const;
  Eigen::VectorXd asDofs_u(const xt::xtensor<double,2> &nodevec) const;
  Eigen::VectorXd asDofs_p(const xt::xtensor<double,1> &dofval ) const;
  Eigen::VectorXd asDofs_p(const xt::xtensor<double,2> &nodevec) const;

};

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#include "MatrixPartitioned.hpp"

// =================================================================================================

#endif

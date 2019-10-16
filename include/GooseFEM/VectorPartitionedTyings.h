/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_VECTORPARTITIONEDTYINGS_H
#define GOOSEFEM_VECTORPARTITIONEDTYINGS_H

// -------------------------------------------------------------------------------------------------

#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>

// =================================================================================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

/*
  "nodevec"   -  nodal vectors                            -  [nnode, ndim]
  "elemvec"   -  nodal vectors stored per element         -  [nelem, nne, ndim]
  "dofval"    -  DOF values                               -  [ndof]
  "dofval_u"  -  DOF values (Unknown)    "== dofval[iiu]" -  [nnu]
  "dofval_p"  -  DOF values (Prescribed) "== dofval[iiu]" -  [nnp]
*/

class VectorPartitionedTyings
{
public:

  // Constructor

  VectorPartitionedTyings() = default;

  VectorPartitionedTyings(
    const xt::xtensor<size_t,2>& conn,
    const xt::xtensor<size_t,2>& dofs,
    const Eigen::SparseMatrix<double>& Cdu,
    const Eigen::SparseMatrix<double>& Cdp,
    const Eigen::SparseMatrix<double>& Cdi);

  // Dimensions

  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t nnode() const; // number of nodes
  size_t ndim()  const; // number of dimensions
  size_t ndof()  const; // number of DOFs
  size_t nnu()   const; // number of independent, unknown DOFs
  size_t nnp()   const; // number of independent, prescribed DOFs
  size_t nni()   const; // number of independent DOFs
  size_t nnd()   const; // number of dependent DOFs

  // DOF lists

  xt::xtensor<size_t,2> dofs() const; // DOFs
  xt::xtensor<size_t,1> iiu()  const; // independent, unknown DOFs
  xt::xtensor<size_t,1> iip()  const; // independent, prescribed DOFs
  xt::xtensor<size_t,1> iii()  const; // independent DOFs
  xt::xtensor<size_t,1> iid()  const; // dependent DOFs

  // Copy (part of) nodevec/dofval to another nodevec/dofval

  void copy_p(
    const xt::xtensor<double,1>& dofval_src,
          xt::xtensor<double,1>& dofval_dest) const; // "iip"  updated

  // Convert to "dofval" (overwrite entries that occur more than once)

  void asDofs_i(
    const xt::xtensor<double,2>& nodevec,
          xt::xtensor<double,1>& dofval_i,
          bool apply_tyings=true) const; // overwritten

  // Convert to "nodevec" (overwrite entries that occur more than once) -- (auto allocation below)

  void asNode(
    const xt::xtensor<double,1>& dofval,
          xt::xtensor<double,2>& nodevec) const;

  // Convert to "elemvec" (overwrite entries that occur more than once) -- (auto allocation below)

  void asElement(
    const xt::xtensor<double,2>& nodevec,
          xt::xtensor<double,3>& elemvec) const;

  // Assemble "dofval" (adds entries that occur more that once) -- (auto allocation below)

  void assembleDofs(
    const xt::xtensor<double,3>& elemvec,
          xt::xtensor<double,1>& dofval) const;

  // Assemble "nodevec" (adds entries that occur more that once) -- (auto allocation below)

  void assembleNode(
    const xt::xtensor<double,3>& elemvec,
          xt::xtensor<double,2>& nodevec) const;

  // Auto-allocation of the functions above

  xt::xtensor<double,1> AsDofs_i(
    const xt::xtensor<double,2>& nodevec) const;

  xt::xtensor<double,2> AsNode(
    const xt::xtensor<double,1>& dofval) const;

  xt::xtensor<double,3> AsElement(
    const xt::xtensor<double,2>& nodevec) const;

  xt::xtensor<double,1> AssembleDofs(
    const xt::xtensor<double,3>& elemvec) const;

  xt::xtensor<double,2> AssembleNode(
    const xt::xtensor<double,3>& elemvec) const;

private:

  // Bookkeeping
  xt::xtensor<size_t,2> m_conn; // connectivity          [nelem, nne ]
  xt::xtensor<size_t,2> m_dofs; // DOF-numbers per node  [nnode, ndim]
  xt::xtensor<size_t,1> m_iiu;  // unknown     DOFs      [nnu]
  xt::xtensor<size_t,1> m_iip;  // prescribed  DOFs      [nnp]
  xt::xtensor<size_t,1> m_iid;  // dependent   DOFs      [nnd]

  // Dimensions
  size_t m_nelem; // number of elements
  size_t m_nne;   // number of nodes per element
  size_t m_nnode; // number of nodes
  size_t m_ndim;  // number of dimensions
  size_t m_ndof;  // number of DOFs
  size_t m_nnu;   // number of independent, unknown DOFs
  size_t m_nnp;   // number of independent, prescribed DOFs
  size_t m_nni;   // number of independent DOFs
  size_t m_nnd;   // number of dependent DOFs

  // Tyings
  Eigen::SparseMatrix<double> m_Cdu;
  Eigen::SparseMatrix<double> m_Cdp;
  Eigen::SparseMatrix<double> m_Cdi;
  Eigen::SparseMatrix<double> m_Cud;
  Eigen::SparseMatrix<double> m_Cpd;
  Eigen::SparseMatrix<double> m_Cid;

  // equivalent Eigen functions

  Eigen::VectorXd Eigen_asDofs_d(
    const xt::xtensor<double,2>& nodevec) const;

};

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#include "VectorPartitionedTyings.hpp"

// =================================================================================================

#endif

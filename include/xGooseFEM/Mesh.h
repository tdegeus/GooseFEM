/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_MESH_H
#define XGOOSEFEM_MESH_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ======================================== GooseFEM::Mesh =========================================

namespace xGooseFEM {
namespace Mesh {

// -------------------------------------------------------------------------------------------------

// list with DOF-numbers in sequential order
inline xt::xtensor<size_t,2> dofs(size_t nnode, size_t ndim);

// renumber to lowest possible numbers (e.g. [0,3,4,2] -> [0,2,3,1])
inline xt::xtensor<size_t,2> renumber(const xt::xtensor<size_t,2> &dofs);

// get the list needed to renumber: dofs_renumbered(i,j) = index(dofs(i,j))
inline xt::xtensor<size_t,1> renumber_index(const xt::xtensor<size_t,2> &dofs);

// renumber such that certain indices "iip" are are moved to the beginning or the end
// (get the lowest or the highest indices); if "iiu" are the remaining indices, after renumbering:
// iiu = arange(nnu), iip = nnu + arange(nnp)
inline xt::xtensor<size_t,2> reorder(const xt::xtensor<size_t,2> &dofs,
  const xt::xtensor<size_t,1> &iip, const std::string &location="end");

// get the list needed to reorder: dofs_reordered(i,j) = index(dofs(i,j))
inline xt::xtensor<size_t,1> reorder_index(const xt::xtensor<size_t,2> &dofs,
  const xt::xtensor<size_t,1> &iip, const std::string &location="end");

// number of elements connected to each node
inline xt::xtensor<size_t,1> coordination(const xt::xtensor<size_t,2> &conn);

// elements connected to each node
inline SpMatS elem2node(const xt::xtensor<size_t,2> &conn);

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

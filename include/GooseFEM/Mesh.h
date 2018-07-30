/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESH_H
#define GOOSEFEM_MESH_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ======================================== GooseFEM::Mesh =========================================

namespace GooseFEM {
namespace Mesh {

// -------------------------------------------------------------------------------------------------

// list with DOF-numbers in sequential order
inline xt::xtensor<size_t,2> dofs(size_t nnode, size_t ndim);

// -------------------------------------------------------------------------------------------------

// renumber to lowest possible numbers (e.g. [0,3,4,2] -> [0,2,3,1])
// - core
template<class InputIterator, class OutputIterator>
inline void renumber(
  const InputIterator first, const InputIterator last, const OutputIterator result
);
// - interface
inline xt::xtensor<size_t,2> renumber(const xt::xtensor<size_t,2> &dofs);

// -------------------------------------------------------------------------------------------------

// reorder (and renumber) such that certain indices "idx" are moved to the beginning or the end
// - core
template<class InputIterator, class OutputIterator, class IndexIterator>
inline void reorder(
  const InputIterator first, const InputIterator last, const OutputIterator result,
  const IndexIterator idx_first, const IndexIterator idx_last, std::string location
);
// - interface
inline xt::xtensor<size_t,2> reorder(const xt::xtensor<size_t,2> &dofs,
  const xt::xtensor<size_t,1> &idx, std::string location="end");

// -------------------------------------------------------------------------------------------------

// elements connected to each node:
// out[: ,0  ] = number of elements connected to each node
// out[j ,i+1] = "i"th element connected to node "j"
inline SpMatS elem2node(const xt::xtensor<size_t,2> &conn);

// -------------------------------------------------------------------------------------------------

// number of elements connected to each node
inline xt::xtensor<size_t,1> coordination(const xt::xtensor<size_t,2> &conn);

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESH_H
#define GOOSEFEM_MESH_H

// -------------------------------------------------------------------------------------------------

#include "config.h"

// =================================================================================================

namespace GooseFEM {
namespace Mesh {

// -------------------------------------------------------------------------------------------------

// Renumber to lowest possible index. For example [0,3,4,2] -> [0,2,3,1]

class Renumber
{
public:

  // constructors

  Renumber() = default;

  Renumber(const xt::xarray<size_t>& dofs);

  // get renumbered DOFs (same as "Renumber::apply(dofs)")

  xt::xtensor<size_t,2> get(const xt::xtensor<size_t,2>& dofs) const;

  // apply renumbering to other set

  template <class T> T apply(const T& list) const;

  // get the list needed to renumber, e.g.:
  //   dofs_renumbered(i,j) = index(dofs(i,j))

  xt::xtensor<size_t,1> index() const;

private:

  xt::xtensor<size_t,1> m_renum;

};

// -------------------------------------------------------------------------------------------------

// Reorder to lowest possible index, in specific order.
//
// For example for "Reorder({iiu,iip})" after reordering:
//
//   iiu = xt::range<size_t>(nnu);
//   iip = xt::range<size_t>(nnp) + nnu;

class Reorder
{
public:

  // constructors

  Reorder() = default;

  Reorder(const std::initializer_list<xt::xtensor<size_t,1>> args);

  // get reordered DOFs (same as "Reorder::apply(dofs)")

  xt::xtensor<size_t,2> get(const xt::xtensor<size_t,2>& dofs) const;

  // apply renumbering to other set

  template <class T> T apply(const T& list) const;

  // apply renumbering

  // get the list needed to reorder, e.g.:
  //   dofs_reordered(i,j) = index(dofs(i,j))

  xt::xtensor<size_t,1> index() const;

private:

  xt::xtensor<size_t,1> m_renum;

};

// -------------------------------------------------------------------------------------------------

// list with DOF-numbers in sequential order
inline xt::xtensor<size_t,2> dofs(size_t nnode, size_t ndim);

// number of elements connected to each node
inline xt::xtensor<size_t,1> coordination(const xt::xtensor<size_t,2> &conn);

// elements connected to each node
inline std::vector<std::vector<size_t>> elem2node(const xt::xtensor<size_t,2> &conn);

// -------------------------------------------------------------------------------------------------

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

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#include "Mesh.hpp"

// =================================================================================================

#endif

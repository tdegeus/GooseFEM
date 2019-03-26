/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESH_HPP
#define GOOSEFEM_MESH_HPP

// -------------------------------------------------------------------------------------------------

#include "Mesh.h"

// =================================================================================================

namespace GooseFEM {
namespace Mesh {

// -------------------------------------------------------------------------------------------------

inline Renumber::Renumber(const xt::xarray<size_t>& dofs)
{
  size_t n = xt::amax(dofs)[0]+1;
  size_t i = 0;

  xt::xtensor<size_t,1> unique = xt::unique(dofs);

  m_renum = xt::empty<size_t>({n});

  for (auto& j : unique) {
    m_renum(j) = i;
    ++i;
  }
}

// -------------------------------------------------------------------------------------------------

// apply renumbering, e.g. for a matrix:
//
//   out(i,j) = renum(list(i,j))

template <class T>
T Renumber::apply(const T& list) const
{
  T out = T::from_shape(list.shape());

  auto jt = out.begin();

  for ( auto it = list.begin() ; it != list.end() ; ++it, ++jt )
    *jt = m_renum(*it);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> Renumber::get(const xt::xtensor<size_t,2>& dofs) const
{
  return this->apply(dofs);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Renumber::index() const
{
  return m_renum;
}

// -------------------------------------------------------------------------------------------------

inline Reorder::Reorder(const std::initializer_list<xt::xtensor<size_t,1>> args)
{
  size_t n = 0;
  size_t i = 0;

  for (auto& arg : args) {
    if (arg.size() == 0) continue;
    n = std::max(n, xt::amax(arg)[0]+1);
  }

  #ifdef GOOSEFEM_ENABLE_ASSERT
    for (auto& arg : args)
      GOOSEFEM_ASSERT(xt::unique(arg) == xt::sort(arg));
  #endif

  m_renum = xt::empty<size_t>({n});

  for (auto& arg : args)
  {
    for (auto& j : arg)
    {
      m_renum(j) = i;
      ++i;
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> Reorder::get(const xt::xtensor<size_t,2>& dofs) const
{
  return this->apply(dofs);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Reorder::index() const
{
  return m_renum;
}

// -------------------------------------------------------------------------------------------------

// apply renumbering, e.g. for a matrix:
//
//   out(i,j) = renum(list(i,j))

template <class T>
T Reorder::apply(const T& list) const
{
  T out = T::from_shape(list.shape());

  auto jt = out.begin();

  for ( auto it = list.begin() ; it != list.end() ; ++it, ++jt )
    *jt = m_renum(*it);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> renumber(const xt::xtensor<size_t,2> &dofs)
{
  return Renumber(dofs).get(dofs);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> dofs(size_t nnode, size_t ndim)
{
  return xt::reshape_view(xt::arange<size_t>(nnode*ndim),{nnode,ndim});
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> coordination(const xt::xtensor<size_t,2> &conn)
{
  // get number of nodes
  size_t nnode = xt::amax(conn)[0] + 1;

  // number of elements connected to each node
  // - allocate
  xt::xtensor<size_t,1> N = xt::zeros<size_t>({nnode});
  // - fill from connectivity
  for ( auto it = conn.begin(); it != conn.end(); ++it ) N(*it) += 1;

  return N;
}

// -------------------------------------------------------------------------------------------------

inline std::vector<std::vector<size_t>> elem2node(const xt::xtensor<size_t,2> &conn)
{
  // get coordination per node
  auto N = coordination(conn);

  // get number of nodes
  auto nnode = N.size();

  // allocate
  std::vector<std::vector<size_t>> out;
  // reserve outer size
  out.resize(nnode);
  // reserve inner sizes
  for ( size_t i = 0 ; i < nnode ; ++i )
    out[i].reserve(N(i));

  // fill
  for ( size_t e = 0 ; e < conn.shape()[0] ; ++e )
    for ( size_t m = 0 ; m < conn.shape()[1] ; ++m )
      out[conn(e,m)].push_back(e);

  return out;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

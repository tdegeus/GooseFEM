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

inline xt::xtensor<size_t,2> dofs(size_t nnode, size_t ndim)
{
  return xt::reshape_view(xt::arange<size_t>(nnode*ndim),{nnode,ndim});
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> renumber_index(const xt::xtensor<size_t,2> &dofs)
{
  // get unique list of DOFs
  xt::xtensor<size_t,1> unique = xt::unique(dofs);

  // allocate list to renumber "dofs"
  xt::xtensor<size_t,1> renum = xt::empty<size_t>({xt::amax(dofs)[0]+1});

  // define renumbering
  for ( size_t i = 0 ; i < unique.size() ; ++i ) renum[unique[i]] = i;

  return renum;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> renumber(const xt::xtensor<size_t,2> &dofs)
{
  // list to renumber "dofs"
  auto renum = renumber_index(dofs);

  // allocate reordered DOFs
  xt::xtensor<size_t,2> dofs_renumbered(dofs.shape());

  // iterator for loop below
  auto jt = dofs_renumbered.begin();

  // loop to renumber: dofs_renumbered(i,j) = renum(dofs(i,j))
  for ( auto it = dofs.begin() ; it != dofs.end() ; ++it, ++jt )
    (*jt) = renum((*it));

  return dofs_renumbered;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> reorder_index(const xt::xtensor<size_t,2> &dofs,
  const xt::xtensor<size_t,1> &iip, const std::string &location)
{
  // check "iip" to be a unique set
  assert( xt::unique(iip).size() == iip.size() );

  // get remaining DOFs
  auto iiu = xt::setdiff1d(dofs, iip);

  // get sizes
  auto nnu = iiu.size();
  auto nnp = iip.size();

  // original set of DOFs
  auto old = xt::unique(dofs);

  // sanity check
  assert( iiu.size() + iip.size() == dofs.size() );

  // list to renumber "dofs"
  // - allocate
  xt::xtensor<size_t,1> renum = xt::empty<size_t>({xt::amax(dofs)[0]+1});
  // - fill
  if ( location == "end" ) {
    for ( size_t i = 0 ; i < iiu.size() ; ++i ) renum(iiu(i)) = i    ;
    for ( size_t i = 0 ; i < iip.size() ; ++i ) renum(iip(i)) = i+nnu;
  }
  else if ( location == "begin" or location == "beginning" ) {
    for ( size_t i = 0 ; i < iip.size() ; ++i ) renum(iip(i)) = i    ;
    for ( size_t i = 0 ; i < iiu.size() ; ++i ) renum(iiu(i)) = i+nnp;
  }
  else {
    throw std::runtime_error("Unknown reorder location '" + location + "'");
  }

  return renum;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> reorder(const xt::xtensor<size_t,2> &dofs,
  const xt::xtensor<size_t,1> &iip, const std::string &location)
{
  // list to renumber "dofs"
  auto renum = reorder_index(dofs, iip, location);

  // allocate reordered DOFs
  xt::xtensor<size_t,2> dofs_reordered(dofs.shape());

  // iterator for loop below
  auto jt = dofs_reordered.begin();

  // loop to renumber: dofs_reordered(i,j) = renum(dofs(i,j))
  for ( auto it = dofs.begin() ; it != dofs.end() ; ++it, ++jt )
    *jt = renum(*it);

  return dofs_reordered;
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

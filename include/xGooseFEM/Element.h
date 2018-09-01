/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_ELEMENT_H
#define XGOOSEFEM_ELEMENT_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ======================================= GooseFEM::Element =======================================

namespace xGooseFEM {
namespace Element {

// -------------------------------------------------------------------------------------------------

// convert nodal vector [nnode, ndim] to nodal vector stored per element [nelem, nne, ndim]
inline xt::xtensor<double,3> asElementVector(
  const xt::xtensor<size_t,2> &conn, const xt::xtensor<double,2> &nodevec);

// assemble nodal vector stored per element [nelem, nne, ndim] to nodal vector [nnode, ndim]
inline xt::xtensor<double,2> assembleNodeVector(
  const xt::xtensor<size_t,2> &conn, const xt::xtensor<double,3> &elemvec);

// check DOFs to leave no holes
template<class E> inline bool isSequential(const E &dofs);

// check structure of the matrices stored per element [nelem, nne*ndim, nne*ndim]
bool isDiagonal(const xt::xtensor<double,3> &elemmat);

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

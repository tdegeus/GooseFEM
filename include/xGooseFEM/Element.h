/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_ELEMENT_H
#define XGOOSEFEM_ELEMENT_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ======================================= xGooseFEM::Element =======================================

namespace xGooseFEM {
namespace Element {

// -------------------------------------------------------------------------------------------------

// convert nodal vector [nnode, ndim] to nodal vector stored per element [nelem, nne, ndim]
inline xt::xtensor<double,3> asElementVector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<double,2> &nodevec);

// assemble nodal vector stored per element [nelem, nne, ndim] to nodal vector [nnode, ndim]
inline xt::xtensor<double,2> assembleNodeVector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<double,3> &elemvec);

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

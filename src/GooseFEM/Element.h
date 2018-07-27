/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENT_H
#define GOOSEFEM_ELEMENT_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ======================================= GooseFEM::Element =======================================

namespace GooseFEM {
namespace Element {

// -------------------------------------------------------------------------------------------------

// convert nodal vector [nnode, ndim] to nodal vector stored per element [nelem, nne, ndim]
inline ArrD asElementVector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<double,2> &nodevec);

// assemble nodal vector stored per element [nelem, nne, ndim] to nodal vector [nnode, ndim]
inline xt::xtensor<double,2> assembleNodeVector(const xt::xtensor<size_t,2> &conn, const ArrD &elemvec);

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

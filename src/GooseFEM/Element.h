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

// global vector (e.g. nodal positions) to element vector [nelem, nne, ndim]
ArrD elementVector(const MatS &conn, const MatD &global);

MatD assembleElementVector(const MatS &conn, const ArrD &local);

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

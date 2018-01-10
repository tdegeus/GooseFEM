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

// =================================================================================================

// elements connected to each node:
// out[ : , 0   ] = number of elements connected to each node
// out[ j , i+1 ] = "i"th element connected to node "j"
inline MatS elem2node ( const MatS &conn );

// list with DOF-numbers in sequential order
inline MatS dofs ( size_t nnode , size_t ndim );

// renumber list of indices
// - renumber to lowest possible numbers (e.g. [0,3,4,2] -> [0,2,3,1])
inline MatS renumber ( const MatS &in );
// - renumber to begin [ idx , ... ] or end [ ... , idx ] (e.g. [0,1,3,2] -> [3,0,2,1]; with idx=0)
inline MatS renumber ( const MatS &in , const ColS &idx , std::string location="end" );

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

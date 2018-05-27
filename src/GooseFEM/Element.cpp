/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENT_CPP
#define GOOSEFEM_ELEMENT_CPP

// -------------------------------------------------------------------------------------------------

#include "Element.h"

// ======================================= GooseFEM::Element =======================================

namespace GooseFEM {
namespace Element {

// -------------------------------------------------------------------------------------------------

inline ArrD asElementVector(const MatS &conn, const MatD &nodevec)
{
  // extract dimensions
  size_t nelem = static_cast<size_t>(conn   .rows());
  size_t nne   = static_cast<size_t>(conn   .cols());
  size_t ndim  = static_cast<size_t>(nodevec.cols());

  // allocate output: nodal vectors stored per element
  ArrD elemvec = ArrD::Zero({nelem, nne, ndim});

  // read from nodal vectors
  #pragma omp parallel for
  for ( size_t e = 0 ; e < nelem ; ++e )
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        elemvec(e,m,i) = nodevec(conn(e,m),i);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline MatD assembleNodeVector(const MatS &conn, const ArrD &elemvec)
{
  // check input
  assert( elemvec.rank() == 3 ); // nodal vector stored per element [nelem, nne, ndim]

  // extract dimensions
  size_t nelem = static_cast<size_t>(conn.rows()); // number of elements
  size_t nne   = static_cast<size_t>(conn.cols()); // number of nodes per element
  size_t ndim  = elemvec.shape(2);                 // number of dimensions
  size_t nnode = conn.maxCoeff()+1;                // number of nodes

  // check input
  assert( elemvec.shape(0) == nelem );
  assert( elemvec.shape(1) == nne   );

  // zero-initialize output: nodal vectors
  MatD nodevec = MatD::Zero(nnode, ndim);

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output: nodal vectors
    MatD t_nodevec = MatD::Zero(nnode, ndim);

    // assemble from nodal vectors stored per element
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
      for ( size_t m = 0 ; m < nne ; ++m )
        for ( size_t i = 0 ; i < ndim ; ++i )
          t_nodevec(conn(e,m),i) += elemvec(e,m,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      nodevec += t_nodevec;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

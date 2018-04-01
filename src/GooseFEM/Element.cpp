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

ArrD asElementVector(const MatS &conn, const MatD &nodevec)
{
  // extract dimensions
  size_t nelem = static_cast<size_t>(conn   .rows());
  size_t nne   = static_cast<size_t>(conn   .cols());
  size_t ndim  = static_cast<size_t>(nodevec.cols());

  // allocate output: nodal vectors stored per element
  ArrD elemvec({nelem, nne, ndim});

  // read from nodal vectors
  #pragma omp parallel for
  for ( size_t e = 0 ; e < nelem ; ++e )
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        elemvec(e,m,i) = nodevec(conn(e,m),i);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

MatD assembleElementVector(const MatS &conn, const ArrD &elemvec)
{
  // check input
  assert( elemvec.ndim() == 3 ); // nodal vector stored per element [nelem, nne, ndim]

  // extract dimensions
  size_t nelem = static_cast<size_t>(conn.rows()); // number of elements
  size_t nne   = static_cast<size_t>(conn.cols()); // number of nodes per element
  size_t ndim  = elemvec.shape(2);                 // number of dimensions
  size_t nnode = conn.maxCoeff()+1;                // number of nodes

  // check input
  assert( elemvec.shape(0) == nelem );
  assert( elemvec.shape(1) == nne   );

  // allocate output: nodal vectors
  MatD nodevec(nnode, ndim);

  // zero-initialize output
  nodevec.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; allocate output: nodal vectors
    MatD t_nodevec(nnode, ndim);

    // - per thread; zero-initialize output
    t_nodevec.setZero();

    // - per thread; assemble from nodal vectors stored per element
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
      for ( size_t m = 0 ; m < nne ; ++m )
        for ( size_t i = 0 ; i < ndim ; ++i )
          t_nodevec(conn(e,m),i) += elemvec(e,m,i);

    // - reduce: combine result obtained on the different threads
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

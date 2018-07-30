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

inline xt::xtensor<double,3> asElementVector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<double,2> &nodevec)
{
  // extract dimensions
  size_t nelem = conn   .shape()[0];
  size_t nne   = conn   .shape()[1];
  size_t ndim  = nodevec.shape()[1];

  // allocate output: nodal vectors stored per element
  xt::xtensor<double,3> elemvec = xt::zeros<double>({nelem, nne, ndim});

  // read from nodal vectors
  #pragma omp parallel for
  for ( size_t e = 0 ; e < nelem ; ++e )
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        elemvec(e,m,i) = nodevec(conn(e,m),i);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> assembleNodeVector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<double,3> &elemvec)
{
  // extract dimensions
  size_t nelem = conn.shape()[0];     // number of elements
  size_t nne   = conn.shape()[1];     // number of nodes per element
  size_t ndim  = elemvec.shape()[2];  // number of dimensions
  size_t nnode = xt::amax(conn)[0]+1; // number of nodes

  // check input
  assert( elemvec.shape()[0] == nelem );
  assert( elemvec.shape()[1] == nne   );

  // zero-initialize output: nodal vectors
  xt::xtensor<double,2> nodevec = xt::zeros<double>({nnode, ndim});

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output: nodal vectors
    xt::xtensor<double,2> t_nodevec = xt::zeros<double>({nnode, ndim});

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

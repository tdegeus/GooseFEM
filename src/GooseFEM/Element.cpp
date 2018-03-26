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

ArrD elementVector(const MatS &conn, const MatD &global)
{
  // extract dimensions
  size_t nelem = static_cast<size_t>(conn  .rows());
  size_t nne   = static_cast<size_t>(conn  .cols());
  size_t ndim  = static_cast<size_t>(global.cols());

  // allocate output element vector
  ArrD local({nelem, nne, ndim});

  // extract element vector for global vector
  #pragma omp parallel for
  for ( auto e = 0 ; e < nelem ; ++e )
    for ( auto m = 0 ; m < nne ; ++m )
      for ( auto i = 0 ; i < ndim ; ++i )
        local(e,m,i) = global(conn(e,m),i);

  return local;
}

// -------------------------------------------------------------------------------------------------

MatD assembleElementVector(const MatS &conn, const ArrD &local)
{
  // check input
  assert( local.ndim() == 3 ); // shape of the matrix [nelem, nne, ndim]

  // extract dimensions
  size_t nelem = static_cast<size_t>(conn.rows()); // number of elements
  size_t nne   = static_cast<size_t>(conn.cols()); // number of nodes per element
  size_t ndim  = local.shape(2);                   // number of dimensions
  size_t nnode = conn.maxCoeff()+1;                // number of nodes

  // check input
  assert( local.shape(0) == nelem );
  assert( local.shape(1) == nne   );

  // allocate/initialize global vector
  MatD out(nnode, ndim);
  out.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // assemble
  #pragma omp parallel
  {
    // - allocate/initialize global vector, per thread
    MatD global(nnode, ndim);
    global.setZero();

    // - assemble, per thread
    #pragma omp for
    for ( auto e = 0 ; e < nelem ; ++e )
      for ( auto m = 0 ; m < nne ; ++m )
        for ( auto i = 0 ; i < ndim ; ++i )
          global(conn(e,m), i) += local(e,m,i);

    // - reduce "global" per thread to total "out"
    #pragma omp critical
      out += global;
  }

  // automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return out;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

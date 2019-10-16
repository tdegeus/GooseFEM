/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENT_HPP
#define GOOSEFEM_ELEMENT_HPP

// -------------------------------------------------------------------------------------------------

#include "Element.h"

// =================================================================================================

namespace GooseFEM {
namespace Element {

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> asElementVector(
  const xt::xtensor<size_t,2>& conn,
  const xt::xtensor<double,2>& nodevec)
{
  size_t nelem = conn.shape(0);
  size_t nne = conn.shape(1);
  size_t ndim = nodevec.shape(1);

  xt::xtensor<double,3> elemvec = xt::empty<double>({nelem, nne, ndim});

  #pragma omp parallel for
  for (size_t e = 0 ; e < nelem ; ++e)
    for (size_t m = 0 ; m < nne ; ++m)
      for (size_t i = 0 ; i < ndim ; ++i)
        elemvec(e,m,i) = nodevec(conn(e,m),i);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> assembleNodeVector(
  const xt::xtensor<size_t,2>& conn,
  const xt::xtensor<double,3>& elemvec)
{
  size_t nelem = conn.shape(0);
  size_t nne = conn.shape(1);
  size_t ndim = elemvec.shape(2);
  size_t nnode = xt::amax(conn)[0]+1;

  GOOSEFEM_ASSERT(elemvec.shape(0) == nelem);
  GOOSEFEM_ASSERT(elemvec.shape(1) == nne);

  xt::xtensor<double,2> nodevec = xt::zeros<double>({nnode, ndim});

  for (size_t e = 0 ; e < nelem ; ++e)
    for (size_t m = 0 ; m < nne ; ++m)
      for (size_t i = 0 ; i < ndim ; ++i)
        nodevec(conn(e,m),i) += elemvec(e,m,i);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

template<class E>
inline bool isSequential(const E& dofs)
{
  size_t ndof = xt::amax(dofs)[0] + 1;

  xt::xtensor<int,1> exists = xt::zeros<int>({ndof});

  for (auto& i: dofs)
    exists[i]++;

  for (auto& i: dofs)
    if (exists[i] == 0)
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

inline bool isDiagonal(const xt::xtensor<double,3>& elemmat)
{
  GOOSEFEM_ASSERT(elemmat.shape(1) == elemmat.shape(2));

  size_t nelem = elemmat.shape(0);
  size_t N = elemmat.shape(1);

  double eps = std::numeric_limits<double>::epsilon();

  #pragma omp parallel for
  for (size_t e = 0 ; e < nelem ; ++e)
    for (size_t i = 0 ; i < N ; ++i)
      for (size_t j = 0 ; j < N ; ++j)
        if (i != j)
          if (std::abs(elemmat(e,i,j)) > eps)
            return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

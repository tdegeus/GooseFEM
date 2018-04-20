/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ITERATE_H
#define GOOSEFEM_ITERATE_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ======================================= GooseFEM::Iterate =======================================

namespace GooseFEM {
namespace Iterate {

// -------------------------------------------------------------------------------------------------

class StopList
{
private:

  // list with residuals
  std::vector<double> m_res;

  // stopping criterion
  double m_norm;

public:

  // constructors
  StopList(){};
  StopList(double norm, size_t n=1);

  // update list of norms, return true if all norms are below the tolerance
  bool stop(double res);

};

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

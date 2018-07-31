/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_DYNAMICS_H
#define XGOOSEFEM_DYNAMICS_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ======================================= xGooseFEM::Dynamics =======================================

namespace xGooseFEM {
namespace Dynamics {

// ------------------------------------------ dummy class ------------------------------------------

class Geometry
{
public:

  // solve for DOF-accelerations [ndof]
  virtual xt::xtensor<double,1> solve_A() { return xt::empty<double>({0}); };
  virtual xt::xtensor<double,1> solve_V() { return xt::empty<double>({0}); };

  // return nodal vectors [nnode, ndim]
  virtual xt::xtensor<double,2> u() const { return xt::empty<double>({0,0}); };
  virtual xt::xtensor<double,2> v() const { return xt::empty<double>({0,0}); };
  virtual xt::xtensor<double,2> a() const { return xt::empty<double>({0,0}); };

  // return DOF values [ndof]
  virtual xt::xtensor<double,1> dofs_u() const { return xt::empty<double>({0}); };
  virtual xt::xtensor<double,1> dofs_v() const { return xt::empty<double>({0}); };
  virtual xt::xtensor<double,1> dofs_a() const { return xt::empty<double>({0}); };

  // overwrite nodal vectors [nnode, ndim]
  virtual void set_u(const xt::xtensor<double,2> &nodevec) { UNUSED(nodevec); return; };

  // overwrite nodal vectors, reconstructed from DOF values [ndof]
  virtual void set_u(const xt::xtensor<double,1> &dofval) { UNUSED(dofval); return; };
  virtual void set_v(const xt::xtensor<double,1> &dofval) { UNUSED(dofval); return; };
  virtual void set_a(const xt::xtensor<double,1> &dofval) { UNUSED(dofval); return; };

};

// ------------------------------------ evaluate one time step -------------------------------------

inline void Verlet        (Geometry &geometry, double dt, size_t nstep=1);
inline void velocityVerlet(Geometry &geometry, double dt, size_t nstep=1);

namespace Overdamped {
inline void forwardEuler  (Geometry &geometry, double dt, size_t nstep=1);
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

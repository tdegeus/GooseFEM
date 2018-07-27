/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_DYNAMICS_H
#define GOOSEFEM_DYNAMICS_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ======================================= GooseFEM::Dynamics =======================================

namespace GooseFEM {
namespace Dynamics {

// ------------------------------------------ dummy class ------------------------------------------

class Geometry
{
public:

  // solve for DOF-accelerations [ndof]
  virtual xt::xtensor<double,1> solve_A() { xt::xtensor<double,1> out = xt::empty<double>({0}); return out; };
  virtual xt::xtensor<double,1> solve_V() { xt::xtensor<double,1> out = xt::empty<double>({0}); return out; };

  // return nodal vectors [nnode, ndim]
  virtual xt::xtensor<double,2> u() const { xt::xtensor<double,2> out = xt::empty<double>({0,0}); return out; };
  virtual xt::xtensor<double,2> v() const { xt::xtensor<double,2> out = xt::empty<double>({0,0}); return out; };
  virtual xt::xtensor<double,2> a() const { xt::xtensor<double,2> out = xt::empty<double>({0,0}); return out; };

  // return DOF values [ndof]
  virtual xt::xtensor<double,1> dofs_u() const { xt::xtensor<double,1> out = xt::empty<double>({0}); return out; };
  virtual xt::xtensor<double,1> dofs_v() const { xt::xtensor<double,1> out = xt::empty<double>({0}); return out; };
  virtual xt::xtensor<double,1> dofs_a() const { xt::xtensor<double,1> out = xt::empty<double>({0}); return out; };

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

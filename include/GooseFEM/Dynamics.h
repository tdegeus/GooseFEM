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
  virtual ColD solve_A() { return ColD(); };
  virtual ColD solve_V() { return ColD(); };

  // return nodal vectors [nnode, ndim]
  virtual MatD u() const { return MatD(); };
  virtual MatD v() const { return MatD(); };
  virtual MatD a() const { return MatD(); };

  // return DOF values [ndof]
  virtual ColD dofs_u() const { return ColD(); };
  virtual ColD dofs_v() const { return ColD(); };
  virtual ColD dofs_a() const { return ColD(); };

  // overwrite nodal vectors [nnode, ndim]
  virtual void set_u(const MatD &nodevec) { UNUSED(nodevec); return; };

  // overwrite nodal vectors, reconstructed from DOF values [ndof]
  virtual void set_u(const ColD &dofval) { UNUSED(dofval); return; };
  virtual void set_v(const ColD &dofval) { UNUSED(dofval); return; };
  virtual void set_a(const ColD &dofval) { UNUSED(dofval); return; };

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

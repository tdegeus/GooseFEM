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

// -------------------------------------------------------------------------------------------------

class Geometry
{
public:

  // solve for DOF-accelerations [ndof]
  virtual ColD solve() const { return ColD(); };

  // reset residuals, check for convergence
  virtual void reset()          { return;       }
  virtual bool stop(double tol) { UNUSED(tol); return false; };

  // process time-step
  virtual void timestep(double dt) { UNUSED(dt); return; };

  // return nodal vectors [nnode, ndim]
  virtual MatD x() const { return MatD(); };
  virtual MatD v() const { return MatD(); };
  virtual MatD a() const { return MatD(); };

  // return DOF values [ndof]
  virtual ColD dofs_v() const { return ColD(); };
  virtual ColD dofs_a() const { return ColD(); };

  // overwrite nodal vectors [nnode, ndim]
  virtual void set_x(const MatD &pvector) { UNUSED(pvector); return; };

  // overwrite nodal vectors, reconstructed from DOF values [ndof]
  virtual void set_v(const ColD &dofval) { UNUSED(dofval); return; };
  virtual void set_a(const ColD &dofval) { UNUSED(dofval); return; };

};

// -------------------------------------------------------------------------------------------------

// evaluate one time step
inline void velocityVerlet(Geometry &geometry, double dt);

// iterate until all nodes have come to a rest
inline size_t quasiStaticVelocityVerlet(Geometry &geometry, double dt, double tol);

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_DYNAMCS_CPP
#define GOOSEFEM_DYNAMCS_CPP

// -------------------------------------------------------------------------------------------------

#include "Dynamics.h"

// ======================================= GooseFEM::Dynamics =======================================

namespace GooseFEM {
namespace Dynamics {

// -------------------------------------------------------------------------------------------------

inline void Verlet(Geometry &g, double dt, size_t nstep)
{
  for ( size_t istep = 0 ; istep < nstep ; ++istep )
  {
    // history

    ColD V;
    ColD A;
    ColD V_n = g.dofs_v();
    ColD A_n = g.dofs_a();

    // new displacement

    g.set_u(MatD( g.u() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a() ));

    // new acceleration

    A = g.solve_A();

    g.set_a(ColD( A ));

    // new velocity

    g.set_v(ColD( V_n + .5 * dt * ( A_n + A ) ));
  }
}

// -------------------------------------------------------------------------------------------------

inline void velocityVerlet(Geometry &g, double dt, size_t nstep)
{
  for ( size_t istep = 0 ; istep < nstep ; ++istep )
  {
    // history

    ColD V;
    ColD A;
    ColD V_n = g.dofs_v();
    ColD A_n = g.dofs_a();

    // new displacement

    g.set_u(MatD( g.u() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a() ));

    // estimate new velocity

    g.set_v(ColD( V_n + dt * A_n ));

    A = g.solve_A();

    g.set_v(ColD( V_n + .5 * dt * ( A_n + A ) ));

    // new velocity

    A = g.solve_A();

    g.set_v(ColD( V_n + .5 * dt * ( A_n + A ) ));

    // new acceleration

    A = g.solve_A();

    g.set_a(ColD( A ));
  }
}

// -------------------------------------------------------------------------------------------------

namespace Overdamped
{
inline void forwardEuler(Geometry &g, double dt, size_t nstep)
{
  for ( size_t istep = 0 ; istep < nstep ; ++istep )
  {
    g.set_u(ColD( g.dofs_u() + dt * g.solve_V() ));
  }
}
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

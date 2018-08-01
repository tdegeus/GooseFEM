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
    // variables & history

    MatD u;
    ColD V;
    ColD A;
    ColD V_n = g.dofs_v();
    ColD A_n = g.dofs_a();

    // new displacement

    u = g.u() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a();

    g.set_u(u);

    // new acceleration

    A = g.solve_A();

    g.set_a(A);

    // new velocity

    V = V_n + .5 * dt * ( A_n + A );

    g.set_v(V);
  }
}

// -------------------------------------------------------------------------------------------------

inline void velocityVerlet(Geometry &g, double dt, size_t nstep)
{
  for ( size_t istep = 0 ; istep < nstep ; ++istep )
  {
    // variables & history

    MatD u;
    ColD V;
    ColD A;
    ColD V_n = g.dofs_v();
    ColD A_n = g.dofs_a();

    // new displacement

    u = g.u() + dt * g.v() + 0.5 * std::pow(dt,2.) * g.a();

    g.set_u(u);

    // estimate new velocity

    V = V_n + dt * A_n;

    g.set_v(V);

    A = g.solve_A();

    V = V_n + .5 * dt * ( A_n + A );

    g.set_v(V);

    // new velocity

    A = g.solve_A();

    V = V_n + .5 * dt * ( A_n + A );

    g.set_v(V);

    // new acceleration

    A = g.solve_A();

    g.set_a(A);
  }
}

// -------------------------------------------------------------------------------------------------

namespace Overdamped
{
inline void forwardEuler(Geometry &g, double dt, size_t nstep)
{
  for ( size_t istep = 0 ; istep < nstep ; ++istep )
  {
    ColD U = g.dofs_u() + dt * g.solve_V();

    g.set_u(U);
  }
}
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

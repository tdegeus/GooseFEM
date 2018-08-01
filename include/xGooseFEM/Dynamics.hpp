/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_DYNAMCS_CPP
#define XGOOSEFEM_DYNAMCS_CPP

// -------------------------------------------------------------------------------------------------

#include "Dynamics.h"

// ======================================= GooseFEM::Dynamics =======================================

namespace xGooseFEM {
namespace Dynamics {

// -------------------------------------------------------------------------------------------------

inline void Verlet(Geometry &g, double dt, size_t nstep)
{
  for ( size_t istep = 0 ; istep < nstep ; ++istep )
  {
    // variables & history

    xt::xtensor<double,2> u;
    xt::xtensor<double,1> V;
    xt::xtensor<double,1> A;
    xt::xtensor<double,1> V_n = g.dofs_v();
    xt::xtensor<double,1> A_n = g.dofs_a();

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

    xt::xtensor<double,2> u;
    xt::xtensor<double,1> V;
    xt::xtensor<double,1> A;
    xt::xtensor<double,1> V_n = g.dofs_v();
    xt::xtensor<double,1> A_n = g.dofs_a();

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
    xt::xtensor<double,1> U = g.dofs_u() + dt * g.solve_V();

    g.set_u(U);
  }
}
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// =================================================================================================

#endif

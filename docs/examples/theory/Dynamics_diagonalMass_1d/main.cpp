
#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <iostream>
#include <math.h>
#include <Eigen/Eigen>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1,              Eigen::ColMajor> ColD;
typedef Eigen::Matrix<size_t, Eigen::Dynamic, 1,              Eigen::ColMajor> ColS;

// -------------------------------------------------------------------------------------------------

ColD computeM(const ColD &x, double rho)
{
  ColD dNdxi(2), dNdx(2), xe(2), M(x.size());

  // integration variables
  double w = 1; // weight (N.B. such that total volume is correct)
  double J, V;

  // shape function gradients : local coordinates, constant over the element
  dNdxi(0) = -0.5;
  dNdxi(1) = +0.5;

  // extract number of elements : number of nodes - 1
  size_t nelem = x.size() - 1;

  // zero initialize force
  M.setZero();

  // loop over elements (and integration points, here equal to one and omitted)
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // element coordinates
    xe(0) = x(e  );
    xe(1) = x(e+1);

    // Jacobian
    J = dNdxi(0) * xe(0) + dNdxi(1) * xe(1);

    // element size (actually the integration point size)
    V = J * w;

    // add to force
    M(e  ) += rho * V;
    M(e+1) += rho * V;
  }

  return M;
}

// -------------------------------------------------------------------------------------------------

ColD computeFu(const ColD &x, const ColD &u, double G)
{
  ColD dNdxi(2), dNdx(2), xe(2), ue(2), F(x.size());

  // integration variables
  double w = 2; // weight (N.B. xi == 0, but irrelevant because dNdxi is constant)
  double J, V;

  // local constitutive response
  double sig, eps;

  // shape function gradients : local coordinates, constant over the element
  dNdxi(0) = -0.5;
  dNdxi(1) = +0.5;

  // extract number of elements : number of nodes - 1
  size_t nelem = x.size() - 1;

  // zero initialize force
  F.setZero();

  // loop over elements (and integration points, here equal to one and omitted)
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // element coordinates
    xe(0) = x(e  );
    xe(1) = x(e+1);

    // element displacements
    ue(0) = u(e  );
    ue(1) = u(e+1);

    // Jacobian
    J = dNdxi(0) * xe(0) + dNdxi(1) * xe(1);

    // element size (actually the integration point size)
    V = J * w;

    // shape function gradients : global coordinates (constant over the element)
    dNdx = ( 1./J ) * dNdxi;

    // strain (constant over the element)
    eps = dNdx(0) * ue(0) + dNdx(1) * ue(1);

    // stress (constant over the element)
    sig = G * eps;

    // add to force
    F(e  ) += dNdx(0) * sig * V;
    F(e+1) += dNdx(1) * sig * V;
  }

  return F;
}

// -------------------------------------------------------------------------------------------------

ColD computeFv(const ColD &x, const ColD &v, double eta)
{
  ColD dNdxi(2), dNdx(2), xe(2), ve(2), F(x.size());

  // integration variables
  double w = 2; // weight (N.B. xi == 0, but irrelevant because dNdxi is constant)
  double J, V;

  // local constitutive response
  double sig, epsdot;

  // shape function gradients : local coordinates, constant over the element
  dNdxi(0) = -0.5;
  dNdxi(1) = +0.5;

  // extract number of elements : number of nodes - 1
  size_t nelem = x.size() - 1;

  // zero initialize force
  F.setZero();

  // loop over elements (and integration points, here equal to one and omitted)
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // element coordinates
    xe(0) = x(e  );
    xe(1) = x(e+1);

    // element displacements
    ve(0) = v(e  );
    ve(1) = v(e+1);

    // Jacobian
    J = dNdxi(0) * xe(0) + dNdxi(1) * xe(1);

    // element size (actually the integration point size)
    V = J * w;

    // shape function gradients : global coordinates (constant over the element)
    dNdx = ( 1./J ) * dNdxi;

    // strain (constant over the element)
    epsdot = dNdx(0) * ve(0) + dNdx(1) * ve(1);

    // stress (constant over the element)
    sig = eta * epsdot;

    // add to force
    F(e  ) += dNdx(0) * sig * V;
    F(e+1) += dNdx(1) * sig * V;
  }

  return F;
}

// -------------------------------------------------------------------------------------------------

MatD velocityVerlet(
  double rho, double G, double eta, double h, const ColD &Fext,
  double dt, size_t ninc, size_t save_every, const ColS &save_nodes
)
{
  // get dimensions from input
  size_t nel   = static_cast<size_t>(Fext.size())-1;
  double L     = h * static_cast<double>(nel);
  double t     = 0.;
  size_t isave = 0;
  size_t jsave = static_cast<size_t>(std::floor(static_cast<double>(ninc)/static_cast<double>(save_every)));
  size_t nsave = static_cast<size_t>(save_nodes.size());

  // allocate output (displacement)
  MatD out(jsave,nsave);

  // linear system
  ColD x   (nel+1);
  ColD u   (nel+1);
  ColD v   (nel+1);
  ColD v_n (nel+1);
  ColD a   (nel+1);
  ColD a_n (nel+1);
  ColD M   (nel+1);
  ColD Minv(nel+1);
  ColD Fu  (nel+1);
  ColD Fv  (nel+1);

  // nodal coordinates
  x = ColD::LinSpaced( nel+1, 0.0, L );

  // zero initialize
  u  .setZero();
  v  .setZero();
  a  .setZero();
  a_n.setZero();
  v_n.setZero();

  // mass matrix and its inverse (constant in time)
  M    = computeM(x,rho);
  Minv = M.cwiseInverse();

  for ( size_t inc = 0 ; inc < ninc ; ++inc )
  {
    // (1) new nodal positions (displacements)
    // - apply update (nodal) : x_{n+1} = x_n + dt * v_n + .5 * dt^2 * a_n"
    u += dt * v + ( .5 * std::pow(dt,2.) ) * a;
    // - compute forces that are dependent on the displacement, but not on the velocity
    Fu = computeFu(x,u,G);

    // (2a) estimate nodal velocities
    // - update velocities (DOFs)
    v.noalias() = v_n + dt * a;
    // - compute forces that are dependent on the velocity
    Fv = computeFv(x,v,eta);
    // - solve for accelerations (DOFs)
    a.noalias() = Minv.cwiseProduct( Fext - Fu - Fv );
    // - apply boundary conditions
    a(0  ) = 0.;
    a(nel) = 0.;
    // - update velocities (DOFs)
    v.noalias() = v_n + ( .5 * dt ) * ( a_n + a );
    // - compute forces that are dependent on the velocity
    Fv = computeFv(x,v,eta);

    // (2b) new nodal velocities
    // - solve for accelerations (DOFs)
    a.noalias() = Minv.cwiseProduct( Fext - Fu - Fv );
    // - apply boundary conditions
    a(0  ) = 0.;
    a(nel) = 0.;
    // - update velocities (DOFs)
    v.noalias() = v_n + ( .5 * dt ) * ( a_n + a );
    // - compute forces that are dependent on the velocity
    Fv = computeFv(x,v,eta);

    // (3) new nodal accelerations
    // - solve for accelerations (DOFs)
    a.noalias() = Minv.cwiseProduct( Fext - Fu - Fv );
    // - apply boundary conditions
    a(0  ) = 0.;
    a(nel) = 0.;

    // store history
    a_n = a;  // accelerations (DOFs)
    v_n = v;  // nodal velocities
    t  += dt; // current time instance

    // store output
    if ( inc % save_every == 0 )
    {
      for ( size_t i=0; i<nsave; ++i ) out(isave,i) = u(save_nodes(i));
      ++isave;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

namespace py = pybind11;

PYBIND11_MODULE(main,m)
{
  m.def(
    "velocityVerlet", &velocityVerlet,
    py::arg("rho"), py::arg("G"), py::arg("eta"), py::arg("h"), py::arg("Fext"), py::arg("dt"),
    py::arg("ninc"), py::arg("save_every"), py::arg("save_nodes")
  );
}


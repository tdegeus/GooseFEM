
#include <GooseFEM/GooseFEM.h>
#include <GMatElastoPlasticQPot/Cartesian2d.h>
#include <LowFive.h>

// -------------------------------------------------------------------------------------------------

namespace GF = GooseFEM;
namespace QD = GooseFEM::Element::Quad4;
namespace GM = GMatElastoPlasticQPot::Cartesian2d;

// -------------------------------------------------------------------------------------------------

inline double sqdot(const xt::xtensor<double,1> &M, const xt::xtensor<double,1> &V)
{
  double out = 0.;

  for ( size_t i = 0 ; i < M.size() ; ++i )
    out += M(i) * V(i) * V(i);

  return out;
}

// -------------------------------------------------------------------------------------------------

int main()
{
  // simulation parameters

  double T     = 60.  ; // total time
  double dt    = 1.e-2; // time increment
  size_t nx    = 60   ; // number of elements in both directions
  double gamma = .05  ; // displacement step

  // get mesh & quadrature

  GF::Mesh::Quad4::Regular mesh(nx,nx,1.);

  xt::xtensor<double,2> coor = mesh.coor();
  xt::xtensor<size_t,2> conn = mesh.conn();
  xt::xtensor<size_t,2> dofs = mesh.dofsPeriodic();

  GF::Vector vector(conn, dofs);

  QD::Quadrature quad(vector.asElement(coor));

  size_t nnode = mesh.nnode();
  size_t ndim  = mesh.ndim();
  size_t nelem = mesh.nelem();
  size_t nip   = quad.nip();

  xt::xtensor<double,2> u    = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> v    = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> a    = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> v_n  = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> a_n  = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fint = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fext = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fres = xt::zeros<double>(coor.shape());

  xt::xtensor<double,4> Eps  = xt::zeros<double>({nelem, nip, ndim, ndim});
  xt::xtensor<double,4> Sig  = xt::zeros<double>({nelem, nip, ndim, ndim});

  // material definition

  GM::Matrix material({nelem, nip});

  xt::xtensor<size_t,2> Ihard = xt::zeros<size_t>({nelem, nip});
  xt::xtensor<size_t,2> Isoft = xt::ones <size_t>({nelem, nip});

  for ( size_t e = 0 ; e < nx*nx/4 ; ++e )
    for ( size_t q = 0 ; q < nip ; ++q )
      Ihard(e,q) = 1;

  Isoft -= Ihard;

  material.setElastic(Ihard, 100., 10.);
  material.setElastic(Isoft, 100.,  1.);

  material.check();

  // mass matrix

  QD::Quadrature nodalQuad(vector.asElement(coor), QD::Nodal::xi(), QD::Nodal::w());

  xt::xtensor<double,2> rho = 1.0 * xt::ones<double>({nelem, nip});

  GF::MatrixDiagonal M(conn, dofs);

  M.assemble( nodalQuad.int_N_scalar_NT_dV(rho) );

  xt::xtensor<double,1> mass = M.asDiagonal();

  // update in macroscopic deformation gradient

  xt::xtensor<double,2> dFbar = xt::zeros<double>({2,2});

  dFbar(0,1) = gamma;

  for ( size_t n = 0 ; n < nnode ; ++n )
    for ( size_t j = 0 ; j < ndim ; ++j )
      for ( size_t k = 0 ; k < ndim ; ++k )
        u(n,j) += dFbar(j,k) * ( coor(n,k) - coor(0,k) );

  // output variables

  xt::xtensor<double,1> Epot = xt::zeros<double>({static_cast<size_t>(T/dt)});
  xt::xtensor<double,1> Ekin = xt::zeros<double>({static_cast<size_t>(T/dt)});
  xt::xtensor<double,1> t    = xt::zeros<double>({static_cast<size_t>(T/dt)});

  // loop over increments

  for ( size_t inc = 0 ; inc < static_cast<size_t>(Epot.size()) ; ++inc )
  {
    // store history

    xt::noalias(v_n) = v;
    xt::noalias(a_n) = a;

    // new displacement

    xt::noalias(u) = u + dt * v + 0.5 * std::pow(dt,2.) * a;

    // compute strain/strain, and corresponding internal

    quad.symGradN_vector(vector.asElement(u), Eps);
    material.Sig(Eps, Sig);
    vector.assembleNode(quad.int_gradN_dot_tensor2_dV(Sig), fint);

    // estimate new velocity

    xt::noalias(v) = v_n + dt * a_n;

    // compute residual force & solve

    xt::noalias(fres) = fext - fint;

    M.solve(fres, a);

    // re-estimate new velocity

    xt::noalias(v) = v_n + .5 * dt * ( a_n + a );

    // compute residual force & solve

    xt::noalias(fres) = fext - fint;

    M.solve(fres, a);

    // new velocity

    xt::noalias(v) = v_n + .5 * dt * ( a_n + a );

    // compute residual force & solve

    xt::noalias(fres) = fext - fint;

    M.solve(fres, a);

    // store output variables

    xt::xtensor<double,2> E  = material.energy(Eps);
    xt::xtensor<double,2> dV = quad.dV();
    xt::xtensor<double,1> V  = vector.asDofs(v);

    t   (inc) = static_cast<double>(inc) * dt;
    Ekin(inc) = 0.5 * sqdot(mass,V);
    Epot(inc) = xt::sum(E*dV)[0];
  }

  // write output variables to file

  HighFive::File file("example.hdf5", HighFive::File::Overwrite);

  LowFive::xtensor::dump(file, "/global/Epot",Epot );
  LowFive::xtensor::dump(file, "/global/Ekin",Ekin );
  LowFive::xtensor::dump(file, "/global/t"   ,t    );
  LowFive::xtensor::dump(file, "/mesh/conn"  ,conn );
  LowFive::xtensor::dump(file, "/mesh/coor"  ,coor );
  LowFive::xtensor::dump(file, "/mesh/disp"  ,u    );

  return 0;
}

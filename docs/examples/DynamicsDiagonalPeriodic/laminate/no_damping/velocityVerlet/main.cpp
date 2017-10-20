
#include <HDF5pp.h>
#include <Eigen/Eigen>
#include <cppmat/cppmat.h>
#include <GooseFEM/GooseFEM.h>
#include <GooseMaterial/AmorphousSolid/LinearStrain/Elastic/Cartesian2d.h>

// -------------------------------------------------------------------------------------------------

using     MatS = GooseFEM::MatS;
using     MatD = GooseFEM::MatD;
using     ColD = GooseFEM::ColD;

using     vec  = cppmat::cartesian2d::vector  <double>;
using     T2   = cppmat::cartesian2d::tensor2 <double>;
using     T2s  = cppmat::cartesian2d::tensor2s<double>;
using     T2d  = cppmat::cartesian2d::tensor2d<double>;

namespace GM   = GooseMaterial::AmorphousSolid::LinearStrain::Elastic::Cartesian2d;

// =================================================================================================

class Quadrature
{
public:
  T2s eps, epsdot, sig;

  size_t nhard;
  GM::Material hard, soft;

  double Ebar, Vbar;

  Quadrature(size_t nhard);

  double density             (size_t elem, size_t k, double V);
  void   stressStrain        (size_t elem, size_t k, double V);
  void   stressStrainRate    (size_t elem, size_t k, double V);
  void   stressStrainPost    (size_t elem, size_t k, double V);
  void   stressStrainRatePost(size_t elem, size_t k, double V);
};

// -------------------------------------------------------------------------------------------------

Quadrature::Quadrature(size_t _nhard)
{
  nhard    = _nhard;
  hard     = GM::Material(100.,10.);
  soft     = GM::Material(100., 1.);
}

// -------------------------------------------------------------------------------------------------

double Quadrature::density(size_t elem, size_t k, double V)
{
  return 1.0;
}
// -------------------------------------------------------------------------------------------------

void Quadrature::stressStrain(size_t elem, size_t k, double V)
{
  if ( elem < nhard ) sig = hard.stress(eps);
  else                sig = soft.stress(eps);
}
// -------------------------------------------------------------------------------------------------

void Quadrature::stressStrainRate(size_t elem, size_t k, double V)
{
  sig.zeros();
}
// -------------------------------------------------------------------------------------------------

void Quadrature::stressStrainPost(size_t elem, size_t k, double V)
{
  Vbar += V;

  if ( elem < nhard ) Ebar += hard.energy(eps) * V;
  else                Ebar += soft.energy(eps) * V;
}

// -------------------------------------------------------------------------------------------------

void Quadrature::stressStrainRatePost(size_t elem, size_t k, double V)
{
}

// =================================================================================================

int main()
{
  // set simulation parameters
  double T     = 20.  ; // total time
  double dt    = 1.e-2; // time increment
  size_t nx    = 40   ; // number of elements in both directions
  double gamma = .05  ; // displacement step

  // class which provides the mesh
  GooseFEM::Mesh::Quad4::Regular mesh(nx,nx,static_cast<double>(nx),static_cast<double>(nx));
  // reference node
  size_t nodeRef = mesh.nodesRef();

  // class which provides the constitutive response at each quadrature point
  auto  quadrature = std::make_shared<Quadrature>(nx*nx/4);

  // class which provides the response of each element
  using Elem = GooseFEM::Dynamics::Diagonal::LinearStrain::Quad4<Quadrature>;
  auto  elem = std::make_shared<Elem>(quadrature);

  // class which provides the system and an increment
  GooseFEM::Dynamics::Diagonal::Periodic<Elem> sim(
    elem,
    mesh.coor(),
    mesh.conn(),
    mesh.dofsPeriodic(),
    dt,
    0.0
  );

  // define update in macroscopic deformation gradient
  T2 dFbar(0.);
  dFbar(0,1) = gamma;

  // update the displacement according to the macroscopic deformation gradient update
  for ( size_t i = 0 ; i < sim.nnode ; ++i )
    for ( size_t j = 0 ; j < sim.ndim ; ++j )
      for ( size_t k = 0 ; k < sim.ndim ; ++k )
        sim.u(i,j) += dFbar(j,k) * ( sim.x0(i,k) - sim.x0(nodeRef,k) );

  // output variables
  ColD Epot(static_cast<int>(T/dt)); Epot.setZero(); // potential energy
  ColD Ekin(static_cast<int>(T/dt)); Ekin.setZero();
  ColD t   (static_cast<int>(T/dt)); t   .setZero();

  // loop over increments
  for ( size_t inc = 0 ; inc < static_cast<size_t>(Epot.size()) ; ++inc )
  {
    // - compute increment
    sim.velocityVerlet();

    // - post: energy based on nodes
    // -- store time
    t(inc) = sim.t;
    // -- kinetic energy
    for ( size_t i = 0 ; i < sim.ndof ; ++i )
      Ekin(inc) += .5 / sim.Minv(i) * std::pow( sim.V(i) , 2. );
    // -- potential energy
    quadrature->Ebar = 0.0;
    quadrature->Vbar = 0.0;
    sim.post();
    Epot(inc) = quadrature->Ebar;
  }

  // write to output file
  H5p::File f = H5p::File("example.hdf5");
  f.write("/global/Epot",Epot       );
  f.write("/global/Ekin",Ekin       );
  f.write("/global/t"   ,t          );
  f.write("/mesh/coor"  ,mesh.coor());
  f.write("/mesh/conn"  ,mesh.conn());
  f.write("/mesh/disp"  ,sim.u      );

  return 0;
}

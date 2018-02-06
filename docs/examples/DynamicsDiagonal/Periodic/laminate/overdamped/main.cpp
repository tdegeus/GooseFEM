
#include <Eigen/Eigen>
#include <cppmat/cppmat.h>
#include <GooseFEM/GooseFEM.h>
#include <GooseMaterial/GooseMaterial.h>
#include <HDF5pp.h>

// -------------------------------------------------------------------------------------------------

// introduce aliases
using ColD = GooseFEM::ColD;
using T2   = cppmat::cartesian2d::tensor2<double>;
using Mat  = GooseMaterial::AmorphousSolid::LinearStrain::Elastic::Cartesian2d::Material;

// =================================================================================================

class Material
{
public:

  // class variables
  // ---------------

  // strain(-rate), stress, background damping [mandatory]
  cppmat::matrix<double> eps, sig, alpha;

  // mesh dimensions
  size_t nelem, nne=4, ndim=2, nip=4;

  // constitutive response [customize]
  std::vector<Mat> material;

  // class functions
  // ---------------

  // constructor [customize]
  Material(size_t nelem, size_t nhard);

  // compute stress for a certain element "e" and integration point "k" [mandatory]
  void updated_eps(size_t e, size_t k);

  // post process variables / functions [customize]
  // - potential energy
  cppmat::matrix<double> Epot;
  // - compute potential energy
  void post();
};

// -------------------------------------------------------------------------------------------------

Material::Material(size_t _nelem, size_t _nhard)
{
  // copy from input [customize]
  nelem = _nelem;

  // allocate symmetric tensors and scalars per element, per integration point [mandatory]
  eps   .resize({nelem,nip,3});
  sig   .resize({nelem,nip,3});
  alpha .resize({nelem,nne  });

  // set mass-density and background damping [customize]
  alpha.setConstant(0.2);

  // post variables [customize]
  Epot.resize({nelem,nip});

  // constitutive response per element (per integration post) [customize]
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    if ( e < _nhard )
    {
      Mat mat = Mat(100.,10.);
      material.push_back(mat);
    }
    else
    {
      Mat mat = Mat(100.,1.);
      material.push_back(mat);
    }
  }
}

// -------------------------------------------------------------------------------------------------

void Material::updated_eps(size_t e, size_t k)
{
  // local views of the global arrays (speeds up indexing, and increases readability)
  cppmat::cartesian2d::tensor2s<double> Eps, Sig;

  // point to correct position in matrix of tensors
  Eps.map(&eps(e,k));

  // compute stress
  Sig = material[e].stress(Eps);

  // copy to matrix of tensors
  std::copy(Sig.begin(), Sig.end(), &sig(e,k));
}

// -------------------------------------------------------------------------------------------------

void Material::post()
{
#pragma omp parallel
{
  // local views of the global arrays (speeds up indexing, and increases readability)
  cppmat::cartesian2d::tensor2s<double> Eps;

  // loop over all elements / integration points
  #pragma omp for
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    for ( size_t k = 0 ; k < nip ; ++k )
    {
      // point to correct position in matrix of tensors
      Eps.map(&eps(e,k));

      // compute energy
      Epot(e,k) = material[e].energy(Eps);
    }
  }
}
}

// =================================================================================================

// introduce aliases
using Mesh       = GooseFEM::Mesh::Quad4::Regular;
using Element    = GooseFEM::OverdampedDynamics::Diagonal::SmallStrain::Quad4<Material>;
using Simulation = GooseFEM::OverdampedDynamics::Diagonal::Periodic<Element>;

// =================================================================================================

int main()
{
  // set simulation parameters
  double T     = 60.  ; // total time
  double dt    = 1.e-3; // time increment
  size_t nx    = 40   ; // number of elements in both directions
  double gamma = .05  ; // displacement step

  // damping and mass-density (to allow conversion)
  double alpha = 0.2;
  double rho   = 1.0;

  // class which provides the mesh
  Mesh mesh(nx,nx,1.);
  // extract information
  size_t nodeOrigin = mesh.nodesOrigin();
  size_t nelem      = mesh.nelem();
  size_t nhard      = nx*nx/4;

  // simulation class
  Simulation sim = Simulation(
    std::make_unique<Element>(std::make_unique<Material>(nelem,nhard),nelem),
    mesh.coor(),
    mesh.conn(),
    mesh.dofsPeriodic(),
    dt
  );

  // define update in macroscopic deformation gradient
  T2 dFbar(0.);
  dFbar(0,1) = gamma;

  // update the displacement according to the macroscopic deformation gradient update
  for ( size_t i = 0 ; i < sim.nnode ; ++i )
    for ( size_t j = 0 ; j < sim.ndim ; ++j )
      for ( size_t k = 0 ; k < sim.ndim ; ++k )
        sim.u(i,j) += dFbar(j,k) * ( sim.x(i,k) - sim.x(nodeOrigin,k) );

  // process updates for displacement dependent variables
  sim.updated_u();

  // output variables
  ColD Epot(static_cast<int>(T/dt)); Epot.setZero();
  ColD Ekin(static_cast<int>(T/dt)); Ekin.setZero();
  ColD t   (static_cast<int>(T/dt)); t   .setZero();

  // mass matrix
  ColD M = rho/alpha * sim.D;

  // loop over increments
  for ( size_t inc = 0 ; inc < static_cast<size_t>(Epot.size()) ; ++inc )
  {
    // - compute increment
    sim.forwardEuler();

    // - store time
    t(inc) = sim.t;

    // - store total kinetic energy
    for ( size_t i = 0 ; i < sim.ndof ; ++i )
      Ekin(inc) += .5 * M(i) * std::pow(sim.V(i),2.);

    // - store total potential energy
    sim.elem->mat->post();
    Epot(inc) = sim.elem->mat->Epot.average(sim.elem->vol) * sim.elem->vol.sum();
  }

  // write to output file
  H5p::File f = H5p::File("example.hdf5","w");
  f.write("/global/Epot",Epot       );
  f.write("/global/Ekin",Ekin       );
  f.write("/global/t"   ,t          );
  f.write("/mesh/coor"  ,mesh.coor());
  f.write("/mesh/conn"  ,mesh.conn());
  f.write("/mesh/disp"  ,sim.u      );

  return 0;
}

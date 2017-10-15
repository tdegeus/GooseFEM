
#include "GooseFEM.h"

#include <GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian2d.h>
#include <cppmat/cppmat.h>

#include <memory>
#include <iostream>

// -------------------------------------------------------------------------------------------------

namespace GM = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d;

using MatD = GooseFEM::MatD;
using MatS = GooseFEM::MatS;
using ColD = GooseFEM::ColD;
using ColS = GooseFEM::ColS;

using T2  = cppmat::cartesian2d::tensor2 <double>;
using T2s = cppmat::cartesian2d::tensor2s<double>;

// =================================================================================================

class Material
{
private:
  GM::Material m_hard;
  GM::Material m_soft;
  GM::Material m_damp;

public:
  double Ebar;
  double Vbar;
  bool   store=false;

  Material();

  double density(size_t elem, size_t quad){return 1.;};

  T2s stressStrain    (size_t elem, size_t quad, const T2s &eps   , double V);
  T2s stressStrainRate(size_t elem, size_t quad, const T2s &epsdot, double V);
};

// =================================================================================================

Material::Material()
{
  m_hard = GM::Material(100.,10.);
  m_soft = GM::Material(100., 1.);
  m_damp = GM::Material(  1., 1.);
}

// =================================================================================================

T2s Material::stressStrain(size_t elem, size_t quad, const T2s &eps, double V)
{
  //50
  // store integration point quantities
  if ( store )
  {
    if ( elem < 2 ) Ebar += m_hard.energy(eps);
    else            Ebar += m_soft.energy(eps);

    Vbar += V;
  }

  // constitutive response
  if ( elem < 2 ) return m_hard.stress(eps);
  else            return m_soft.stress(eps);
}

// =================================================================================================

T2s Material::stressStrainRate(size_t elem, size_t quad, const T2s &epsdot, double V)
{
  return m_damp.stress(epsdot);
}

// =================================================================================================

int main()
{
  // -----------------------------------------------------------------------------------------------
  // define problem
  // -----------------------------------------------------------------------------------------------

  // define a mesh
  GooseFEM::Mesh::Quad4::Regular mesh(2,2);
  // extract a reference node (whose position will stay constant)
  size_t nodeRef = mesh.nodesRef();

  // initialize the material definition (see above)
  auto mat = std::make_shared<Material>();

  // initialize the simulation
  GooseFEM::Dynamics::Periodic::LinearStrain::DiagonalMass::Simulation\
  <Material,GooseFEM::Quad4,cppmat::cartesian2d::tensor2s<double>>\
  sim(
    mat,
    std::make_shared<GooseFEM::Quad4>(),
    mesh.coor(),
    mesh.conn(),
    mesh.dofsPeriodic(),
    1.e-3
  );

  // -----------------------------------------------------------------------------------------------
  // apply an affine deformation
  // -----------------------------------------------------------------------------------------------

  // update in macroscopic deformation gradient
  // - allocate null tensor
  T2 dFbar(0.);
  // - set non-zero components
  dFbar(0,1) = .01;

  // compute displacement increment that results from the macroscopic deformation gradient
  for ( size_t i = 0 ; i < sim.nnode ; ++i )
    for ( size_t j = 0 ; j < sim.ndim ; ++j )
      for ( size_t k = 0 ; k < sim.ndim ; ++k )
        sim.u(i,j) += dFbar(j,k) * ( sim.x0(i,k) - sim.x0(nodeRef,k) );




  // std::cout << sim.u << std::endl;



  mat->store = true;
  mat->Vbar  = 0.;
  mat->Ebar  = 0.;

  sim.computeInternalForce();

  std::cout << sim.Fint << std::endl;

  mat->store = false;


  return 0;
}

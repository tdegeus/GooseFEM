/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_OVERDAMPEDDYNAMICS_DIAGONAL_H
#define GOOSEFEM_OVERDAMPEDDYNAMICS_DIAGONAL_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ============================ GooseFEM::OverdampedDynamics::Diagonal =============================

namespace GooseFEM {
namespace OverdampedDynamics {
namespace Diagonal {

// ===================================== Simulation - Periodic =====================================

template<class Element>
class Periodic
{
public:

  // variables
  // ---------

  // element/quadrature/material definition
  std::unique_ptr<Element> elem;

  // mesh: nodal position/displacement/velocity, DOF-numbers, connectivity, dimensions
  size_t nnode, nelem, nne, ndim, ndof;
  MatS   conn, dofs;
  MatD   x, u, u_n, v;

  // linear system: columns (also "D" which is diagonal)
  ColD   D, Dinv, F, V;

  // time integration
  double dt, t=0.0;

  // constructor
  // -----------

  Periodic(){};

  Periodic(
    std::unique_ptr<Element> elem, const MatD &x, const MatS &conn, const MatS &dofs, double dt
  );

  // functions
  // ---------

  void forwardEuler();             // one time step of the time integrator
  void midpoint();                 // one time step of the time integrator
  void updated_x();                // process update in "x"
  void updated_u(bool init=false); // process update in "u", if init all possible updates are made
  void updated_v(bool init=false); // process update in "v", if init all possible updates are made
  void assemble_D();               // assemble the damping matrix from the element damping matrices
  void assemble_F();               // assemble the force from the element forces
};

// =================================== Simulation - SemiPeriodic ===================================

template<class Element>
class SemiPeriodic
{
public:

  // variables
  // ---------

  // element/quadrature/material definition
  std::unique_ptr<Element> elem;

  // mesh: nodal position/displacement/velocity, DOF-numbers, connectivity, dimensions
  size_t nnode, nelem, nne, ndim, ndof;
  MatS   conn, dofs;
  MatD   x, u, u_n, v;

  // fixed DOFs: prescribed velocity, DOF-numbers, dimensions
  size_t nfixed;
  ColS   fixedDofs;
  ColD   fixedV;

  // linear system: columns (also "D" which is diagonal)
  ColD   D, Dinv, F, V;

  // time integration
  double dt, t=0.0;

  // constructor
  // -----------

  SemiPeriodic(){};

  SemiPeriodic(
    std::unique_ptr<Element> elem, const MatD &x, const MatS &conn, const MatS &dofs,
    const ColS &fixedDofs, double dt
  );

  // functions
  // ---------

  void forwardEuler();             // one time step of the time integrator
  void midpoint();                 // one time step of the time integrator
  void updated_x();                // process update in "x"
  void updated_u(bool init=false); // process update in "u", if init all possible updates are made
  void updated_v(bool init=false); // process update in "v", if init all possible updates are made
  void assemble_D();               // assemble the damping matrix from the element damping matrices
  void assemble_F();               // assemble the force from the element forces
};

// ======================================== Element - Quad4 ========================================

namespace SmallStrain {

template<class Material>
class Quad4
{
public:

  // class variables
  // ---------------

  // constitutive response per integration point, per element
  std::unique_ptr<Material> mat;

  // matrices to store the element data
  cppmat::matrix<double> x, u, v, f, D, dNx, dNxi, w, vol, dNxi_n, w_n, vol_n;

  // dimensions
  size_t nelem, nip=4, nne=4, ndim=2;

  // signals for solvers
  bool changed_f=false, changed_D=true;

  // class functions
  // ---------------

  // constructor
  Quad4(std::unique_ptr<Material> mat, size_t nelem);

  // recompute relevant quantities after "x", "u", or "v" have been externally updated
  void updated_x();
  void updated_u();
  void updated_v();

  // return volume average stress and strain
  cppmat::cartesian2d::tensor2s<double> mean_eps(size_t e);   // of element "e"
  cppmat::cartesian2d::tensor2s<double> mean_sig(size_t e);   // of element "e"
  cppmat::cartesian2d::tensor2s<double> mean_eps();           // of all elements
  cppmat::cartesian2d::tensor2s<double> mean_sig();           // of all elements
};

} // namespace ...

// =================================================================================================

}}} // namespace ...

// =================================================================================================

#endif

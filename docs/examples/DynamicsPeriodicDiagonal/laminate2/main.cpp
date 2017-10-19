
#include <Eigen/Eigen>
#include <HDF5pp.h>
#include <cppmat/cppmat.h>
#include <GooseFEM/GooseFEM.h>
#include <GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian2d.h>

using MatD = GooseFEM::MatD;
using ColD = GooseFEM::ColD;



namespace GM = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d;

using vec = cppmat::cartesian2d::vector  <double>;
using T2  = cppmat::cartesian2d::tensor2 <double>;
using T2s = cppmat::cartesian2d::tensor2s<double>;
using T2d = cppmat::cartesian2d::tensor2d<double>;


// =================================================================================================

class Element
{
public:
  // arrays / matrices
  cppmat::tiny::matrix2<double,4,2> xe, ue, ve, xi, xi_n, dNdxi, dNdx;
  cppmat::tiny::vector <double,4>   N, w, w_n;
  cppmat::tiny::matrix2<double,8,8> M;
  cppmat::tiny::vector <double,8>   fu, fv;
  // tensors
  cppmat::cartesian2d::tensor2 <double> J, Jinv, gradu, gradv;
  cppmat::cartesian2d::tensor2s<double> eps, epsdot, sig;
  cppmat::cartesian2d::vector  <double> v;
  // scalars
  double Jdet, V;
  // sizes (nodes per element, dimensions, quadrature points)
  size_t nne=4, ndim=2, nk=4;

  // problem definition
  // - density
  double rho;
  // - constitutive model / Rayleigh damping
  GM::Material hard, soft, rayleigh;
  // - number of elements
  size_t nelem;

  // output quantities
  double Epotbar; // potential energy
  double Ekinbar; // kinetic   energy
  double Vbar;    // volume

  Element(double eta, size_t nelem);

  void computeMassMatrix   (size_t elem);
  void computeInternalForce(size_t elem);
  void computeDampingForce (size_t elem);
  void postProcess         (size_t elem);

};

// -------------------------------------------------------------------------------------------------

Element::Element(double eta, size_t _nelem)
{
  // constitutive definition
  hard     = GM::Material(100.,10.);
  soft     = GM::Material(100., 1.);
  rayleigh = GM::Material(eta ,eta);
  rho      = 1.;
  nelem    = _nelem;

  // integration point coordinates/weights: normal Gauss integration for Quad4
  xi(0,0) = -1./std::sqrt(3.); xi(0,1) = -1./std::sqrt(3.); w(0) = 1.;
  xi(1,0) = +1./std::sqrt(3.); xi(1,1) = -1./std::sqrt(3.); w(1) = 1.;
  xi(2,0) = +1./std::sqrt(3.); xi(2,1) = +1./std::sqrt(3.); w(2) = 1.;
  xi(3,0) = -1./std::sqrt(3.); xi(3,1) = +1./std::sqrt(3.); w(3) = 1.;

  // integration point coordinates/weights: integration at the nodes for Quad4
  xi_n(0,0) = -1.; xi_n(0,1) = -1.; w_n(0) = 1.;
  xi_n(1,0) = +1.; xi_n(1,1) = -1.; w_n(1) = 1.;
  xi_n(2,0) = +1.; xi_n(2,1) = +1.; w_n(2) = 1.;
  xi_n(3,0) = -1.; xi_n(3,1) = +1.; w_n(3) = 1.;

}

// -------------------------------------------------------------------------------------------------

void Element::computeMassMatrix(size_t elem)
{
  // zero-initialize element mass matrix
  M.zeros();

  // loop over integration points (coincide with the nodes to get a diagonal mass matrix)
  for ( size_t k = 0 ; k < nne ; ++k )
  {
    // - shape function gradients (local coordinates)
    dNdxi(0,0) = -.25*(1.-xi_n(k,1)); dNdxi(0,1) = -.25*(1.-xi_n(k,0));
    dNdxi(1,0) = +.25*(1.-xi_n(k,1)); dNdxi(1,1) = -.25*(1.+xi_n(k,0));
    dNdxi(2,0) = +.25*(1.+xi_n(k,1)); dNdxi(2,1) = +.25*(1.+xi_n(k,0));
    dNdxi(3,0) = -.25*(1.+xi_n(k,1)); dNdxi(3,1) = +.25*(1.-xi_n(k,0));

    // - Jacobian
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume (== volume associated with the node, in this element)
    V = w_n(k) * Jdet;

    // - assemble to element mass matrix (use the delta properties of the shape functions)
    for ( size_t i = 0 ; i < ndim ; ++i )
      M(i,i) += rho * V;
  }
}

// -------------------------------------------------------------------------------------------------

void Element::computeInternalForce(size_t elem)
{
  // zero-initialize element force
  fu.zeros();

  // loop over integration points
  for ( size_t k = 0 ; k < nk ; ++k )
  {
    // - shape function gradients (local coordinates)
    dNdxi(0,0) = -.25*(1.-xi(k,1)); dNdxi(0,1) = -.25*(1.-xi(k,0));
    dNdxi(1,0) = +.25*(1.-xi(k,1)); dNdxi(1,1) = -.25*(1.+xi(k,0));
    dNdxi(2,0) = +.25*(1.+xi(k,1)); dNdxi(2,1) = +.25*(1.+xi(k,0));
    dNdxi(3,0) = -.25*(1.+xi(k,1)); dNdxi(3,1) = +.25*(1.-xi(k,0));

    // - Jacobian
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;

    // - shape function gradients (global coordinates)
    dNdx.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          dNdx(m,i) += Jinv(i,j) * dNdxi(m,j);

    // - displacement gradient
    gradu.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradu(i,j) += dNdx(m,i) * ue(m,j);

    // - strain tensor (symmetric part of "gradu")
    eps = gradu.astensor2s();

    // - constitutive response
    if ( elem < (nelem/4)*nelem ) sig = hard.stress(eps);
    else                          sig = soft.stress(eps);

    // - assemble to element force
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          fu(m*ndim+j) += dNdx(m,i) * sig(i,j) * V;
  }
}

// -------------------------------------------------------------------------------------------------

void Element::computeDampingForce(size_t elem)
{
  // zero-initialize element force
  fv.zeros();

  // loop over integration points
  for ( size_t k = 0 ; k < nk ; ++k )
  {
    // - shape function gradients (local coordinates)
    dNdxi(0,0) = -.25*(1.-xi(k,1)); dNdxi(0,1) = -.25*(1.-xi(k,0));
    dNdxi(1,0) = +.25*(1.-xi(k,1)); dNdxi(1,1) = -.25*(1.+xi(k,0));
    dNdxi(2,0) = +.25*(1.+xi(k,1)); dNdxi(2,1) = +.25*(1.+xi(k,0));
    dNdxi(3,0) = -.25*(1.+xi(k,1)); dNdxi(3,1) = +.25*(1.-xi(k,0));

    // - Jacobian
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;

    // - shape function gradients (global coordinates)
    dNdx.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          dNdx(m,i) += Jinv(i,j) * dNdxi(m,j);

    // - velocity gradient
    gradv.zeros();

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradv(i,j) += dNdx(m,i) * ve(m,j);

    // - strain rate tensor (symmetric part of "gradv")
    epsdot = gradv.astensor2s();

    // - constitutive response
    sig = rayleigh.stress(epsdot);

    // - assemble to element force
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          fv(m*ndim+j) += dNdx(m,i) * sig(i,j) * V;
  }
}

// -------------------------------------------------------------------------------------------------

void Element::postProcess(size_t elem)
{
  // loop over integration points
  for ( size_t k = 0 ; k < nk ; ++k )
  {
    // - shape functions
    N(0) = .25*(1.-xi(k,0))*(1.-xi(k,1));
    N(1) = .25*(1.+xi(k,0))*(1.-xi(k,1));
    N(2) = .25*(1.+xi(k,0))*(1.+xi(k,1));
    N(3) = .25*(1.-xi(k,0))*(1.+xi(k,1));

    // - shape function gradients (local coordinates)
    dNdxi(0,0) = -.25*(1.-xi(k,1)); dNdxi(0,1) = -.25*(1.-xi(k,0));
    dNdxi(1,0) = +.25*(1.-xi(k,1)); dNdxi(1,1) = -.25*(1.+xi(k,0));
    dNdxi(2,0) = +.25*(1.+xi(k,1)); dNdxi(2,1) = +.25*(1.+xi(k,0));
    dNdxi(3,0) = -.25*(1.+xi(k,1)); dNdxi(3,1) = +.25*(1.-xi(k,0));

    // - Jacobian
    J.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    // - determinant and inverse of the Jacobian
    Jdet = J.det();
    Jinv = J.inv();

    // - integration point volume
    V = w(k) * Jdet;
    // - add to total volume
    Vbar += V;

    // - shape function gradients (global coordinates)
    dNdx.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          dNdx(m,i) += Jinv(i,j) * dNdxi(m,j);

    // - displacement gradient
    gradu.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradu(i,j) += dNdx(m,i) * ue(m,j);

    // - strain tensor (symmetric part of "gradu")
    eps = gradu.astensor2s();

    // - add to total strain energy
    if ( elem < (nelem/4)*nelem ) Epotbar += hard.energy(eps) * V;
    else                          Epotbar += soft.energy(eps) * V;

    // - velocity
    v.zeros();
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        v(i) += N(m) * ve(m,i);

    // - add to total kinetic energy
    Ekinbar += .5 * rho * v.dot(v) * V;

  }
}

// =================================================================================================

int main(int argc, const char** argv)
{
  // get "eta" and "name" of the output file from the input
  std::string num1 = argv[1];  size_t N     = std::stod(num1);
  std::string num2 = argv[2];  double eta   = std::stod(num2);
  std::string num3 = argv[3];  double alpha = std::stod(num3);
  std::string num4 = argv[4];  double dt    = std::stod(num4);
  std::string name = argv[5];

  // define a mesh
  GooseFEM::Mesh::Quad4::Regular mesh(N,N,static_cast<double>(N),static_cast<double>(N));
  // extract a reference node (whose position will stay constant)
  size_t nodeRef = mesh.nodesRef();

  // define the problem specific "Element" routines (see above)
  auto elem = std::make_shared<Element>(eta,N);

  // define the entire problem
  GooseFEM::Dynamics::Periodic::DiagonalMass<Element> sim(
    elem,
    mesh.coor(),
    mesh.conn(),
    mesh.dofsPeriodic(),
    dt,
    alpha
  );

  // define update in macroscopic deformation gradient
  T2 dFbar(0.);
  dFbar(0,1) = .01;

  // update the displacement according to the macroscopic deformation gradient update
  for ( size_t i = 0 ; i < sim.nnode ; ++i )
    for ( size_t j = 0 ; j < sim.ndim ; ++j )
      for ( size_t k = 0 ; k < sim.ndim ; ++k )
        sim.u(i,j) += dFbar(j,k) * ( sim.x0(i,k) - sim.x0(nodeRef,k) );

  // output : total energy and potential energy
  ColD Epotbar(static_cast<int>(15./dt));
  ColD Ekinbar(static_cast<int>(15./dt));
  ColD t      (static_cast<int>(15./dt));

  // loop over increments
  for ( size_t i = 0 ; i < static_cast<size_t>(Epotbar.size()) ; ++i )
  {
    // - compute increment
    sim.Verlet();

    // - post-process
    // -- zero-initialize averages
    elem->Vbar    = 0.;
    elem->Epotbar = 0.;
    elem->Ekinbar = 0.;
    // -- compute sums
    sim.postProcess();
    // -- normalize and store averages
    Epotbar(i) = elem->Epotbar;
    Ekinbar(i) = elem->Ekinbar;
    t      (i) = sim.t;
  }

  // write to output file
  H5p::File f = H5p::File(name);
  f.write("/global/Epot",Epotbar    );
  f.write("/global/Ekin",Ekinbar    );
  f.write("/global/t"   ,t          );
  f.write("/mesh/coor"  ,mesh.coor());
  f.write("/mesh/conn"  ,mesh.conn());
  f.write("/mesh/disp"  ,sim.u      );

  return 0;
}

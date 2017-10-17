
#include <Eigen/Eigen>
#include <HDF5pp.h>
#include <cppmat/cppmat.h>
#include <GooseFEM/GooseFEM.h>
#include <GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian2d.h>

using MatD = GooseFEM::MatD;
using ColD = GooseFEM::ColD;



namespace GM = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d;

using T2  = cppmat::cartesian2d::tensor2 <double>;
using T2s = cppmat::cartesian2d::tensor2s<double>;
using T2d = cppmat::cartesian2d::tensor2d<double>;


// =================================================================================================

class Element
{
public:
  cppmat::tiny::matrix2<double,4,2> xe, ue, ve, xi, xi_n, dNdxi, dNdx;
  cppmat::tiny::vector <double,4>   w, w_n;
  cppmat::tiny::matrix2<double,8,8> M;
  cppmat::tiny::vector <double,8>   fu, fv;

  cppmat::cartesian2d::tensor2<double> J, Jinv, gradu, gradv;
  cppmat::cartesian2d::tensor2s<double> eps, epsdot, sig;

  double Jdet, V;
  size_t nne=4, ndim=2, nk=4;

  double rho;
  GM::Material hard, soft, visco;

  double Ebar;
  double Vbar;

  Element(double eta);

  void computeMassMatrix   (size_t elem);
  void computeInternalForce(size_t elem);
  void computeDampingForce (size_t elem);

};

// -------------------------------------------------------------------------------------------------

Element::Element(double eta)
{
  hard  = GM::Material(100.,10.);
  soft  = GM::Material(100., 1.);
  visco = GM::Material(eta ,eta);
  rho   = 1.;

  xi(0,0) = -1./std::sqrt(3.); xi(0,1) = -1./std::sqrt(3.); w(0) = 1.;
  xi(1,0) = +1./std::sqrt(3.); xi(1,1) = -1./std::sqrt(3.); w(1) = 1.;
  xi(2,0) = +1./std::sqrt(3.); xi(2,1) = +1./std::sqrt(3.); w(2) = 1.;
  xi(3,0) = -1./std::sqrt(3.); xi(3,1) = +1./std::sqrt(3.); w(3) = 1.;

  xi_n(0,0) = -1.; xi_n(0,1) = -1.; w_n(0) = 1.;
  xi_n(1,0) = +1.; xi_n(1,1) = -1.; w_n(1) = 1.;
  xi_n(2,0) = +1.; xi_n(2,1) = +1.; w_n(2) = 1.;
  xi_n(3,0) = -1.; xi_n(3,1) = +1.; w_n(3) = 1.;

}

// -------------------------------------------------------------------------------------------------

void Element::computeMassMatrix(size_t elem)
{
  M.zeros();

  for ( size_t k = 0 ; k < nne ; ++k )
  {
    dNdxi(0,0) = -.25*(1.-xi_n(k,1)); dNdxi(0,1) = -.25*(1.-xi_n(k,0));
    dNdxi(1,0) = +.25*(1.-xi_n(k,1)); dNdxi(1,1) = -.25*(1.+xi_n(k,0));
    dNdxi(2,0) = +.25*(1.+xi_n(k,1)); dNdxi(2,1) = +.25*(1.+xi_n(k,0));
    dNdxi(3,0) = -.25*(1.+xi_n(k,1)); dNdxi(3,1) = +.25*(1.-xi_n(k,0));

    J.zeros();

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    Jdet = J.det();
    Jinv = J.inv();

    V = w_n(k) * Jdet;

    for ( size_t i = 0 ; i < ndim ; ++i )
      M(i,i) += rho * V;
  }

}

// -------------------------------------------------------------------------------------------------

void Element::computeInternalForce(size_t elem)
{
  fu.zeros();

  for ( size_t k = 0 ; k < nk ; ++k )
  {
    dNdxi(0,0) = -.25*(1.-xi(k,1)); dNdxi(0,1) = -.25*(1.-xi(k,0));
    dNdxi(1,0) = +.25*(1.-xi(k,1)); dNdxi(1,1) = -.25*(1.+xi(k,0));
    dNdxi(2,0) = +.25*(1.+xi(k,1)); dNdxi(2,1) = +.25*(1.+xi(k,0));
    dNdxi(3,0) = -.25*(1.+xi(k,1)); dNdxi(3,1) = +.25*(1.-xi(k,0));

    J.zeros();

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    Jdet = J.det();
    Jinv = J.inv();

    V = w(k) * Jdet;

    dNdx.zeros();

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          dNdx(m,i) += Jinv(i,j) * dNdxi(m,j);

    gradu.zeros();

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradu(i,j) += dNdx(m,i) * ue(m,j);

    eps = gradu.astensor2s();

    if ( elem < 10 ) sig   = hard.stress(eps);
    else             sig   = soft.stress(eps);

    if ( elem < 10 ) Ebar += hard.energy(eps);
    else             Ebar += soft.energy(eps);

    Vbar += V;

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          fu(m*ndim+j) += dNdx(m,i) * sig(i,j) * V;
  }
}

// -------------------------------------------------------------------------------------------------

void Element::computeDampingForce(size_t elem)
{
  fv.zeros();

  for ( size_t k = 0 ; k < nk ; ++k )
  {
    dNdxi(0,0) = -.25*(1.-xi(k,1)); dNdxi(0,1) = -.25*(1.-xi(k,0));
    dNdxi(1,0) = +.25*(1.-xi(k,1)); dNdxi(1,1) = -.25*(1.+xi(k,0));
    dNdxi(2,0) = +.25*(1.+xi(k,1)); dNdxi(2,1) = +.25*(1.+xi(k,0));
    dNdxi(3,0) = -.25*(1.+xi(k,1)); dNdxi(3,1) = +.25*(1.-xi(k,0));

    J.zeros();

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          J(i,j) += dNdxi(m,i) * xe(m,j);

    Jdet = J.det();
    Jinv = J.inv();

    V = w(k) * Jdet;

    dNdx.zeros();

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          dNdx(m,i) += Jinv(i,j) * dNdxi(m,j);

    gradv.zeros();

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          gradv(i,j) += dNdx(m,i) * ve(m,j);

    epsdot = gradv.astensor2s();

    sig = visco.stress(epsdot);

    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        for ( size_t j = 0 ; j < ndim ; ++j )
          fv(m*ndim+j) += dNdx(m,i) * sig(i,j) * V;
  }
}

// =================================================================================================

int main(int argc, const char** argv)
{
  std::string num    = argv[1];
  double      eta    = std::stod(num);
  std::string name   = argv[2];

  // define a mesh
  GooseFEM::Mesh::Quad4::Regular mesh(5,5,5.,5.);
  // extract a reference node (whose position will stay constant)
  size_t nodeRef = mesh.nodesRef();

  auto elem = std::make_shared<Element>(eta);

  GooseFEM::Dynamics::Periodic::DiagonalMass<Element> sim(
    elem,
    mesh.coor(),
    mesh.conn(),
    mesh.dofsPeriodic(),
    nodeRef,
    5.e-3
  );

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

  // sim.computeInternalForce();

  H5p::File f = H5p::File(name);




  ColD Ebar(30000);

  for ( size_t i = 0 ; i < static_cast<size_t>(Ebar.size()) ; ++i )
  {
    elem->Vbar = 0.;
    elem->Ebar = 0.;

    sim.increment();

    Ebar(i) = elem->Ebar / elem->Vbar;
  }

  f.write("/global/E",Ebar);

  f.write("/mesh/coor",mesh.coor());
  f.write("/mesh/conn",mesh.conn());
  f.write("/mesh/disp",sim.u      );

  // f.write("/global/E",Ebar);


  // // MatD xi = GooseFEM::Quadrature::Gauss::Positions::Quad4();

  // // std::cout << xi << std::endl;

  // GooseFEM::Element::Quad4 elem = GooseFEM::Element::Quad4(
  //   GooseFEM::Quadrature::Nodal::Positions::Quad4(),
  //   GooseFEM::Quadrature::Nodal::Weights::Quad4()
  // );

  // std::cout << GooseFEM::Quadrature::Nodal::Positions::Quad4() << std::endl;

  // MatD xe(4,2);

  // xe(0,0) = 0.; xe(0,1) = 0.;
  // xe(1,0) = 1.; xe(1,1) = 0.;
  // xe(2,0) = 1.; xe(2,1) = 1.;
  // xe(3,0) = 0.; xe(3,1) = 1.;

  // elem.setNodalPositions(xe);

  // MatD ue(4,2);

  // ue(0,0) = 0.; ue(0,1) = 0.;
  // ue(1,0) = 0.; ue(1,1) = 0.;
  // ue(2,0) = .1; ue(2,1) = 0.;
  // ue(3,0) = .1; ue(3,1) = 0.;

  // MatD gradu = elem.gradNT_nodalVectors(0,ue);

  // std::cout << gradu << std::endl;
  // std::cout << elem.V(0) << std::endl;




  return 0;
}

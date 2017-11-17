
#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

// -------------------------------------------------------------------------------------------------

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatS;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;

// -------------------------------------------------------------------------------------------------

using T2  = cppmat::cartesian2d::tensor2 <double>;
using T2s = cppmat::cartesian2d::tensor2s<double>;
using T2d = cppmat::cartesian2d::tensor2d<double>;

// =================================================================================================
// simulation: global arrays
// =================================================================================================

template<class Element>
class Simulation
{
public:

  // variables
  // ---------

  // element/quadrature/material definition
  std::unique_ptr<Element> elem;

  // mesh : dimensions
  size_t nnode, nelem, nne, ndim, ndof;

  // mesh : nodal quantities and connectivity
  MatS dofs;  // DOF-numbers of each node     [nnode,ndim]
  MatS conn;  // node numbers of each element [nelem,nne ]
  MatD x;     // positions of each node       [nnode,ndim]
  MatD u;     // displacements of each node   [nnode,ndim]

  // linear system
  ColD F;     // internal force               [ndof]

  // functions
  // ---------

  // constructor
  Simulation(std::unique_ptr<Element> elem, const MatD &x, const MatS &conn, const MatS &dofs);

  // process update in "x" or "u"
  void updated_x();
  void updated_u();

  // assemble the force from the element forces
  void assemble_F();
};

// -------------------------------------------------------------------------------------------------

template<class Element>
Simulation<Element>::Simulation(
  std::unique_ptr<Element> _elem, const MatD &_x, const MatS &_conn, const MatS &_dofs
)
{
  // copy input
  elem  = std::move(_elem);
  x     = _x;
  conn  = _conn;
  dofs  = _dofs;

  // compute sizes
  nnode = static_cast<size_t>(x.rows());
  ndim  = static_cast<size_t>(x.cols());
  nelem = static_cast<size_t>(conn.rows());
  nne   = static_cast<size_t>(conn.cols());
  ndof  = dofs.maxCoeff()+1;

  // allocate and zero-initialize nodal quantities
  u.conservativeResize(nnode,ndim);
  u.setZero();

  // allocate and zero-initialize linear system
  F.conservativeResize(ndof);
  F.setZero();

  // initialize all fields
  updated_x();
  updated_u();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
void Simulation<Element>::updated_x()
{
  // set the nodal positions of all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < nelem ; ++e )
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        elem->x(e,m,i) = x(conn(e,m),i);

  // signal update
  elem->updated_x();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
void Simulation<Element>::updated_u()
{
  // set the nodal displacements of all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < nelem ; ++e )
    for ( size_t m = 0 ; m < nne ; ++m )
      for ( size_t i = 0 ; i < ndim ; ++i )
        elem->u(e,m,i) = u(conn(e,m),i);

  // signal update
  elem->updated_u();

  // update
  assemble_F();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
void Simulation<Element>::assemble_F()
{
  // zero-initialize
  F.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);
  // assemble
  #pragma omp parallel
  {
    // - force, per thread
    ColD F_(ndof);
    F_.setZero();

    // - assemble force, per thread
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
      for ( size_t m = 0 ; m < nne ; ++m )
        for ( size_t i = 0 ; i < ndim ; ++i )
          F_(dofs(conn(e,m),i)) += elem->f(e,m,i);

    // - reduce "F_" per thread to total "F"
    #pragma omp critical
      F += F_;
  }
  // automatic parallelization by Eigen
  Eigen::setNbThreads(0);
}

// =================================================================================================
// all the element operations (fully decoupled), depends on the constitutive response in "Material"
// =================================================================================================

template<class Material>
class Element
{
public:

  // class variables
  // ---------------

  // constitutive response per integration point, per element
  std::unique_ptr<Material> mat;

  // matrices to store the element data
  cppmat::matrix<double> x, u, f, dNx, dNxi, w, V;

  // dimensions
  size_t nelem, nip=4, nne=4, ndim=2;

  // class functions
  // ---------------

  // constructor
  Element(std::unique_ptr<Material> mat, size_t nelem);

  // recompute relevant quantities after "x", "u" have been externally updated
  void updated_x();
  void updated_u();
};

// -------------------------------------------------------------------------------------------------

template<class Material>
Element<Material>::Element(std::unique_ptr<Material> _mat, size_t _nelem)
{
  // copy from input
  nelem = _nelem;
  mat   = std::move(_mat);

  // allocate matrices
  // -
  x   .resize({nelem,nne,ndim});
  u   .resize({nelem,nne,ndim});
  f   .resize({nelem,nne,ndim});
  // -
  dNx .resize({nelem,nip,nne,ndim});
  // -
  V   .resize({nelem,nip});
  // -
  dNxi.resize({nip,nne,ndim});
  // -
  w   .resize({nip});

  // shape function gradient at all Gauss points, in local coordinates
  // - k == 0
  dNxi(0,0,0) = -.25*(1.+1./std::sqrt(3.));     dNxi(0,0,1) = -.25*(1.+1./std::sqrt(3.));
  dNxi(0,1,0) = +.25*(1.+1./std::sqrt(3.));     dNxi(0,1,1) = -.25*(1.-1./std::sqrt(3.));
  dNxi(0,2,0) = +.25*(1.-1./std::sqrt(3.));     dNxi(0,2,1) = +.25*(1.-1./std::sqrt(3.));
  dNxi(0,3,0) = -.25*(1.-1./std::sqrt(3.));     dNxi(0,3,1) = +.25*(1.+1./std::sqrt(3.));
  // - k == 1
  dNxi(1,0,0) = -.25*(1.+1./std::sqrt(3.));     dNxi(1,0,1) = -.25*(1.-1./std::sqrt(3.));
  dNxi(1,1,0) = +.25*(1.+1./std::sqrt(3.));     dNxi(1,1,1) = -.25*(1.+1./std::sqrt(3.));
  dNxi(1,2,0) = +.25*(1.-1./std::sqrt(3.));     dNxi(1,2,1) = +.25*(1.+1./std::sqrt(3.));
  dNxi(1,3,0) = -.25*(1.-1./std::sqrt(3.));     dNxi(1,3,1) = +.25*(1.-1./std::sqrt(3.));
  // - k == 2
  dNxi(2,0,0) = -.25*(1.-1./std::sqrt(3.));     dNxi(2,0,1) = -.25*(1.-1./std::sqrt(3.));
  dNxi(2,1,0) = +.25*(1.-1./std::sqrt(3.));     dNxi(2,1,1) = -.25*(1.+1./std::sqrt(3.));
  dNxi(2,2,0) = +.25*(1.+1./std::sqrt(3.));     dNxi(2,2,1) = +.25*(1.+1./std::sqrt(3.));
  dNxi(2,3,0) = -.25*(1.+1./std::sqrt(3.));     dNxi(2,3,1) = +.25*(1.-1./std::sqrt(3.));
  // - k == 3
  dNxi(3,0,0) = -.25*(1.-1./std::sqrt(3.));     dNxi(3,0,1) = -.25*(1.+1./std::sqrt(3.));
  dNxi(3,1,0) = +.25*(1.-1./std::sqrt(3.));     dNxi(3,1,1) = -.25*(1.-1./std::sqrt(3.));
  dNxi(3,2,0) = +.25*(1.+1./std::sqrt(3.));     dNxi(3,2,1) = +.25*(1.-1./std::sqrt(3.));
  dNxi(3,3,0) = -.25*(1.+1./std::sqrt(3.));     dNxi(3,3,1) = +.25*(1.+1./std::sqrt(3.));

  // integration point weight at all Gauss points
  w(0) = 1.;
  w(1) = 1.;
  w(2) = 1.;
  w(3) = 1.;

  // Note, the above is a specialization of the following:
  //
  // - Shape function gradients
  //
  //    dNxi(0,0) = -.25*(1.-xi(k,1)); dNxi(0,1) = -.25*(1.-xi(k,0));
  //    dNxi(1,0) = +.25*(1.-xi(k,1)); dNxi(1,1) = -.25*(1.+xi(k,0));
  //    dNxi(2,0) = +.25*(1.+xi(k,1)); dNxi(2,1) = +.25*(1.+xi(k,0));
  //    dNxi(3,0) = -.25*(1.+xi(k,1)); dNxi(3,1) = +.25*(1.-xi(k,0));
  //
  // - Gauss point coordinates and weights
  //
  //    xi(0,0) = -1./std::sqrt(3.); xi(0,1) = -1./std::sqrt(3.); w(0) = 1.;
  //    xi(1,0) = +1./std::sqrt(3.); xi(1,1) = -1./std::sqrt(3.); w(1) = 1.;
  //    xi(2,0) = +1./std::sqrt(3.); xi(2,1) = +1./std::sqrt(3.); w(2) = 1.;
  //    xi(3,0) = -1./std::sqrt(3.); xi(3,1) = +1./std::sqrt(3.); w(3) = 1.;
}

// -------------------------------------------------------------------------------------------------

template<class Material>
void Element<Material>::updated_x()
{
#pragma omp parallel
{
  // intermediate quantities
  T2     J_, Jinv_;
  double Jdet_;
  // local views of the global arrays (speeds up indexing, and increases readability)
  cppmat::tiny::matrix2<double,8,8> M_, D_;
  cppmat::tiny::matrix2<double,4,2> dNxi_, dNx_, x_;
  cppmat::tiny::vector <double,4>   w_, V_;

  // loop over all elements (in parallel)
  #pragma omp for
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // pointer to element positions, element integration volume and weight
    x_.map(&x(e));
    V_.map(&V(e));
    w_.map(&w(0));

    // loop over Gauss points
    for ( size_t k = 0 ; k < nip ; ++k )
    {
      // - pointer to the shape function gradients (local/global coordinates)
      dNxi_.map(&dNxi(  k));
      dNx_ .map(&dNx (e,k));

      // - Jacobian
      //   J(i,j) += dNxi(m,i) * xe(m,j)
      J_(0,0) = dNxi_(0,0)*x_(0,0) + dNxi_(1,0)*x_(1,0) + dNxi_(2,0)*x_(2,0) + dNxi_(3,0)*x_(3,0);
      J_(0,1) = dNxi_(0,0)*x_(0,1) + dNxi_(1,0)*x_(1,1) + dNxi_(2,0)*x_(2,1) + dNxi_(3,0)*x_(3,1);
      J_(1,0) = dNxi_(0,1)*x_(0,0) + dNxi_(1,1)*x_(1,0) + dNxi_(2,1)*x_(2,0) + dNxi_(3,1)*x_(3,0);
      J_(1,1) = dNxi_(0,1)*x_(0,1) + dNxi_(1,1)*x_(1,1) + dNxi_(2,1)*x_(2,1) + dNxi_(3,1)*x_(3,1);

      // - determinant and inverse of the Jacobian
      Jdet_ = J_.det();
      Jinv_ = J_.inv();

      // - integration point volume
      V_(k) = w_(k) * Jdet_;

      // - shape function gradients (global coordinates)
      //   dNx(m,i) += Jinv(i,j) * dNxi(m,j)
      for ( size_t m = 0 ; m < nne ; ++m )
      {
        dNx_(m,0) = Jinv_(0,0) * dNxi_(m,0) + Jinv_(0,1) * dNxi_(m,1);
        dNx_(m,1) = Jinv_(1,0) * dNxi_(m,0) + Jinv_(1,1) * dNxi_(m,1);
      }
    }
  }
} // #pragma omp parallel
}

// =================================================================================================

template<class Material>
void Element<Material>::updated_u()
{
#pragma omp parallel
{
  // intermediate quantities
  T2  gradu_;
  T2s eps_, sig_;
  // local views of the global arrays (speeds up indexing, and increases readability)
  cppmat::tiny::matrix2<double,4,2> dNx_, u_, f_;
  cppmat::tiny::vector <double,4>   V_;

  // loop over all elements (in parallel)
  #pragma omp for
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    // pointer to element forces, displacements, and integration volume
    f_.map(&f(e));
    u_.map(&u(e));
    V_.map(&V(e));

    // zero initialize forces
    f_.setZero();

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < nip ; ++k )
    {
      // - pointer to the shape function gradients, strain and stress tensor (stored symmetric)
      dNx_.map(&dNx     (e,k));
      eps_.map(&mat->eps(e,k));
      sig_.map(&mat->sig(e,k));

      // - displacement gradient
      //   gradu_(i,j) += dNx(m,i) * ue(m,j)
      gradu_(0,0) = dNx_(0,0)*u_(0,0) + dNx_(1,0)*u_(1,0) + dNx_(2,0)*u_(2,0) + dNx_(3,0)*u_(3,0);
      gradu_(0,1) = dNx_(0,0)*u_(0,1) + dNx_(1,0)*u_(1,1) + dNx_(2,0)*u_(2,1) + dNx_(3,0)*u_(3,1);
      gradu_(1,0) = dNx_(0,1)*u_(0,0) + dNx_(1,1)*u_(1,0) + dNx_(2,1)*u_(2,0) + dNx_(3,1)*u_(3,0);
      gradu_(1,1) = dNx_(0,1)*u_(0,1) + dNx_(1,1)*u_(1,1) + dNx_(2,1)*u_(2,1) + dNx_(3,1)*u_(3,1);

      // - strain (stored symmetric)
      //   eps(i,j) = .5 * ( gradu_(i,j) + gradu_(j,i) )
      eps_(0,0) =        gradu_(0,0);
      eps_(0,1) = .5 * ( gradu_(0,1) + gradu_(1,0) );
      eps_(1,1) =        gradu_(1,1);

      // - constitutive response
      mat->updated_eps(e,k);

      // - assemble to element force
      //   f(m,j) += dNx(m,i) * sig(i,j) * V;
      for ( size_t m = 0 ; m < nne ; ++m )
      {
        f_(m,0) += dNx_(m,0) * sig_(0,0) * V_(k) + dNx_(m,1) * sig_(1,0) * V_(k);
        f_(m,1) += dNx_(m,0) * sig_(0,1) * V_(k) + dNx_(m,1) * sig_(1,1) * V_(k);
      }
    }
  }
}
}

// =================================================================================================
// material model
// =================================================================================================

class Elasticity
{
private:
  double m_K; // bulk  modulus
  double m_G; // shear modulus

public:
  Elasticity(){};
  Elasticity(double K, double G);

  // compute stress or the energy at "Eps"
  T2s stress(const T2s &Eps);
};

// ===================================== IMPLEMENTATION : CORE =====================================

Elasticity::Elasticity(double K, double G)
{
  // copy input - elastic moduli
  m_K = K;
  m_G = G;
}

// -------------------------------------------------------------------------------------------------

T2s Elasticity::stress(const T2s &Eps)
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cppmat::cartesian2d::identity2();
  double epsm = Eps.trace()/2.;
  T2s    Epsd = Eps - epsm*I;

  // return stress tensor
  return ( m_K * epsm ) * I + m_G * Epsd;
}

// =================================================================================================
// constitutive response of each integration point (decoupled)
// =================================================================================================

class Material
{
public:

  // class variables
  // ---------------

  // strain and stress
  cppmat::matrix<double> eps, sig;

  // dimensions
  size_t nelem, nne=4, ndim=2, nip=4;

  // constitutive response
  std::vector<Elasticity> material;

  // class functions
  // ---------------

  // constructor
  Material(size_t nelem);

  // compute stress for one integration point
  void updated_eps(size_t e, size_t k);
};

// -------------------------------------------------------------------------------------------------

Material::Material(size_t _nelem)
{
  // copy from input
  nelem = _nelem;

  // allocate symmetric tensors of each integration point
  eps.resize({nelem,nip,3});
  sig.resize({nelem,nip,3});

  // constitutive response per element
  for ( size_t e = 0 ; e < nelem ; ++e )
  {
    Elasticity mat = Elasticity(1.,1.);
    material.push_back(mat);
  }
}

// -------------------------------------------------------------------------------------------------

void Material::updated_eps(size_t e, size_t k)
{
  // local views of the global arrays (speeds up indexing, and increases readability)
  T2s Eps, Sig;

  // pointer to stress/strain
  Eps.map(&eps(e,k));
  Sig.map(&sig(e,k));

  // compute stress
  Sig.copy( material[e].stress(Eps) );
}

// =================================================================================================
// some example, with geometry
//
//      3       4       5
//   <==+-------+-------+==>
//      |       |       |
//      |   0   |   1   |
//      |       |       |
//   <==+-------+-------+==>
//      0       1       2
//
// =================================================================================================

int main()
{
  // dimensions
  size_t ndim  = 2;
  size_t nne   = 4;
  size_t nelem = 2;
  size_t nnode = 6;

  // allocate mesh
  MatS conn(nelem,nne );
  MatD coor(nnode,ndim);
  MatS dofs(nnode,ndim);

  // specify mesh
  // - nodal coordinates
  coor(0,0) = 0.; coor(0,1) = 0.;
  coor(1,0) = 1.; coor(1,1) = 0.;
  coor(2,0) = 2.; coor(2,1) = 0.;
  coor(3,0) = 0.; coor(3,1) = 1.;
  coor(4,0) = 1.; coor(4,1) = 1.;
  coor(5,0) = 2.; coor(5,1) = 1.;
  // - DOF numbers per node
  dofs(0,0) =  0; dofs(0,1) =  1;
  dofs(1,0) =  2; dofs(1,1) =  3;
  dofs(2,0) =  4; dofs(2,1) =  5;
  dofs(3,0) =  6; dofs(3,1) =  7;
  dofs(4,0) =  8; dofs(4,1) =  9;
  dofs(5,0) = 10; dofs(5,1) = 11;
  // - connectivity
  conn(0,0) = 0; conn(0,1) = 1; conn(0,2) = 4; conn(0,3) = 3;
  conn(1,0) = 1; conn(1,1) = 2; conn(1,2) = 5; conn(1,3) = 4;

  // define "simulation"
  Simulation<Element<Material>> sim(
    std::make_unique<Element<Material>>(
      std::make_unique<Material>(nelem),
      nelem
    ),
    coor,
    conn,
    dofs
  );

  // set hypothetical displacement
  sim.u(0,0) = -0.1; sim.u(3,0) = -0.1;
  sim.u(1,0) =  0.0; sim.u(4,0) =  0.0;
  sim.u(2,0) = +0.1; sim.u(5,0) = +0.1;

  // compute the internal force
  sim.updated_u();

  // print result
  std::cout << sim.F << std::endl;

}

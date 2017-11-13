
#include <GooseFEM/GooseFEM.h>
#include <cppmat/cppmat.h>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatS;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;
typedef Eigen::Matrix<size_t, Eigen::Dynamic,              1, Eigen::ColMajor> ColS;

// =================================================================================================
// simulation: global arrays
// =================================================================================================

template<class Element>
class Simulation
{
public:

  // class variables
  // ---------------

  // decoupled element response
  std::unique_ptr<Element> elem;

  // mesh : nodal quantities and connectivity
  MatS   dofs;  // DOF-numbers of each node     [nnode,ndim]
  MatS   conn;  // node numbers of each element [nelem,nne ]
  MatD   x;     // positions of each node       [nnode,ndim]
  MatD   u;     // displacements of each node   [nnode,ndim]

  // mesh : size
  size_t nnode, ndim, nelem, nne, ndof;

  // linear system
  ColD F;       // internal force [ndof]

  // class functions
  // ---------------

  // constructor
  Simulation(std::unique_ptr<Element> elem, const MatD &x, const MatS &conn, const MatS &dofs);

  // convert global arrays to arrays per element (making them fully decoupled)
  void system2elem_x();
  void system2elem_u();

  // compute global internal force
  void compute_F();
};

// -------------------------------------------------------------------------------------------------

template<class Element>
Simulation<Element>::Simulation(
  std::unique_ptr<Element> _elem, const MatD &_x, const MatS &_conn, const MatS &_dofs
)
{
  // copy from input
  x     = _x;
  conn  = _conn;
  dofs  = _dofs;
  elem  = std::move(_elem);

  // extract mesh size
  nnode = static_cast<size_t>(x  .rows());
  ndim  = static_cast<size_t>(x  .cols());
  nelem = static_cast<size_t>(conn.rows());
  nne   = static_cast<size_t>(conn.cols());
  ndof  = dofs.maxCoeff()+1;

  // allocate and zero-initialize nodal quantities
  u.conservativeResize(nnode,ndim);
  u.setZero();

  // allocate and zero-initialize linear system
  F.conservativeResize(ndof);
  F.setZero();
}

// -------------------------------------------------------------------------------------------------

template<class Element>
void Simulation<Element>::system2elem_x()
{
  #pragma omp parallel
  {
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
    {
      for ( size_t m = 0 ; m < nne ; ++m )
      {
        for ( size_t i = 0 ; i < ndim ; ++i )
        {
          elem->x(e,m,i) = x(conn(e,m),i);
        }
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class Element>
void Simulation<Element>::system2elem_u()
{
  #pragma omp parallel
  {
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
    {
      for ( size_t m = 0 ; m < nne ; ++m )
      {
        for ( size_t i = 0 ; i < ndim ; ++i )
        {
          elem->u(e,m,i) = u(conn(e,m),i);
        }
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class Element>
void Simulation<Element>::compute_F()
{
  // convert global arrays to arrays per element (making them fully decoupled)
  system2elem_x();
  system2elem_u();

  // compute element response (gradient of the shape functions -> strain -> stress -> force)
  elem->compute_dNdx();
  elem->compute_eps();
  elem->compute_sig();
  elem->compute_f();

  // global internal force
  // - zero-initialize
  F.setZero();
  // - assemble
  #pragma omp parallel
  {
    // -- total force, per thread
    ColD Ft(ndof);
    Ft.setZero();

    // -- compute total force by looping over elements
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
    {
      for ( size_t m = 0 ; m < nne ; ++m )
      {
        for ( size_t i = 0 ; i < ndim ; ++i )
        {
          Ft(dofs(conn(e,m),i)) += elem->f(e,m,i);
        }
      }
    }

    // -- reduce "Ft" per thread to one column "F"
    #pragma omp critical
    {
      F += Ft;
    }
  }
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

  std::unique_ptr<Material> mat;
  cppmat::matrix<double> x, u, f, dNdx, dNdxi, w, V;
  size_t nelem, nip=4, nne=4, ndim=2;

  // class functions
  // ---------------

  Element(std::unique_ptr<Material> mat, size_t nelem);

  void compute_dNdx();
  void compute_eps();
  void compute_sig();
  void compute_f();
};

// -------------------------------------------------------------------------------------------------

template<class Material>
Element<Material>::Element(std::unique_ptr<Material> _mat, size_t _nelem)
{
  // copy from input
  nelem = _nelem;
  mat   = std::move(_mat);

  // allocate matrices
  // - element vectors
  x    .resize({nelem,nne,ndim});
  u    .resize({nelem,nne,ndim});
  f    .resize({nelem,nne,ndim});
  // - element vectors at the integration point
  dNdx .resize({nelem,nip,nne,ndim});
  // - element scalars at the integration point
  V    .resize({nelem,nip});
  // - constant integration point vectors
  dNdxi.resize({nip,nne,ndim});
  // - constant integration point scalars
  w    .resize({nip});

  // set shape function gradients in local coordinates
  // - k == 0
  dNdxi(0,0,0) = -.25*(1.+1./std::sqrt(3.)); dNdxi(0,0,1) = -.25*(1.+1./std::sqrt(3.));
  dNdxi(0,1,0) = +.25*(1.+1./std::sqrt(3.)); dNdxi(0,1,1) = -.25*(1.-1./std::sqrt(3.));
  dNdxi(0,2,0) = +.25*(1.-1./std::sqrt(3.)); dNdxi(0,2,1) = +.25*(1.-1./std::sqrt(3.));
  dNdxi(0,3,0) = -.25*(1.-1./std::sqrt(3.)); dNdxi(0,3,1) = +.25*(1.+1./std::sqrt(3.));
  // - k == 1
  dNdxi(1,0,0) = -.25*(1.+1./std::sqrt(3.)); dNdxi(1,0,1) = -.25*(1.-1./std::sqrt(3.));
  dNdxi(1,1,0) = +.25*(1.+1./std::sqrt(3.)); dNdxi(1,1,1) = -.25*(1.+1./std::sqrt(3.));
  dNdxi(1,2,0) = +.25*(1.-1./std::sqrt(3.)); dNdxi(1,2,1) = +.25*(1.+1./std::sqrt(3.));
  dNdxi(1,3,0) = -.25*(1.-1./std::sqrt(3.)); dNdxi(1,3,1) = +.25*(1.-1./std::sqrt(3.));
  // - k == 2
  dNdxi(2,0,0) = -.25*(1.-1./std::sqrt(3.)); dNdxi(2,0,1) = -.25*(1.-1./std::sqrt(3.));
  dNdxi(2,1,0) = +.25*(1.-1./std::sqrt(3.)); dNdxi(2,1,1) = -.25*(1.+1./std::sqrt(3.));
  dNdxi(2,2,0) = +.25*(1.+1./std::sqrt(3.)); dNdxi(2,2,1) = +.25*(1.+1./std::sqrt(3.));
  dNdxi(2,3,0) = -.25*(1.+1./std::sqrt(3.)); dNdxi(2,3,1) = +.25*(1.-1./std::sqrt(3.));
  // - k == 3
  dNdxi(3,0,0) = -.25*(1.-1./std::sqrt(3.)); dNdxi(3,0,1) = -.25*(1.+1./std::sqrt(3.));
  dNdxi(3,1,0) = +.25*(1.-1./std::sqrt(3.)); dNdxi(3,1,1) = -.25*(1.-1./std::sqrt(3.));
  dNdxi(3,2,0) = +.25*(1.+1./std::sqrt(3.)); dNdxi(3,2,1) = +.25*(1.-1./std::sqrt(3.));
  dNdxi(3,3,0) = -.25*(1.+1./std::sqrt(3.)); dNdxi(3,3,1) = +.25*(1.+1./std::sqrt(3.));

  // set integration point weights
  w(0) = 1.;
  w(1) = 1.;
  w(2) = 1.;
  w(3) = 1.;
}

// -------------------------------------------------------------------------------------------------

template<class Material>
void Element<Material>::compute_dNdx()
{
  #pragma omp parallel
  {
    cppmat::cartesian2d::tensor2<double> J, Jinv;
    double Jdet;

    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
    {
      for ( size_t k = 0 ; k < nip ; ++k )
      {
        // - Jacobian
        //   J(i,j) += dNdxi(m,i) * xe(m,j)
        J.zeros();
        for ( size_t m = 0 ; m < nne ; ++m )
        {
          J(0,0) += dNdxi(k,m,0) * x(e,m,0);
          J(0,1) += dNdxi(k,m,0) * x(e,m,1);
          J(1,0) += dNdxi(k,m,1) * x(e,m,0);
          J(1,1) += dNdxi(k,m,1) * x(e,m,1);
        }

        // - determinant and inverse of the Jacobian
        Jdet = J.det();
        Jinv = J.inv();

        // - integration point volume
        V(e,k) = w(k) * Jdet;

        // - shape function gradients (global coordinates)
        //   dNdx(m,i) += Jinv(i,j) * dNdxi(m,j)
        for ( size_t m = 0 ; m < nne ; ++m )
        {
          dNdx(e,k,m,0) = Jinv(0,0) * dNdxi(k,m,0) + Jinv(0,1) * dNdxi(k,m,1);
          dNdx(e,k,m,1) = Jinv(1,0) * dNdxi(k,m,0) + Jinv(1,1) * dNdxi(k,m,1);
        }
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class Material>
void Element<Material>::compute_eps()
{
  #pragma omp parallel
  {
    cppmat::cartesian2d::tensor2<double> gradu;

    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
    {
      for ( size_t k = 0 ; k < nip ; ++k )
      {
        // - displacement gradient
        //   gradu(i,j) += dNdx(m,i) * ue(m,j)
        gradu.zeros();
        for ( size_t m = 0 ; m < nne ; ++m )
        {
          gradu(0,0) += dNdx(e,k,m,0) * u(e,m,0);
          gradu(0,1) += dNdx(e,k,m,0) * u(e,m,1);
          gradu(1,0) += dNdx(e,k,m,1) * u(e,m,0);
          gradu(1,1) += dNdx(e,k,m,1) * u(e,m,1);
        }

        // - strain
        //   eps(i,j) = .5 * ( gradu(i,j) + gradu(j,i) )
        mat->eps(e,k,0,0) =        gradu(0,0);
        mat->eps(e,k,0,1) = .5 * ( gradu(0,1) + gradu(1,0) );
        mat->eps(e,k,1,0) = .5 * ( gradu(0,1) + gradu(1,0) );
        mat->eps(e,k,1,1) =        gradu(1,1);
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class Material>
void Element<Material>::compute_sig()
{
  mat->compute_sig();
}

// -------------------------------------------------------------------------------------------------

template<class Material>
void Element<Material>::compute_f()
{
  f.setZero();

  #pragma omp parallel
  {
    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
    {
      for ( size_t k = 0 ; k < nip ; ++k )
      {
        // - assemble to element force
        //   fe(m,j) += dNdx(m,i) * sig(i,j) * V;
        for ( size_t m = 0 ; m < nne ; ++m )
        {
          f(e,m,0) += dNdx(e,k,m,0)*mat->sig(e,k,0,0)*V(e,k) + dNdx(e,k,m,1)*mat->sig(e,k,1,0)*V(e,k);
          f(e,m,1) += dNdx(e,k,m,0)*mat->sig(e,k,0,1)*V(e,k) + dNdx(e,k,m,1)*mat->sig(e,k,1,1)*V(e,k);
        }
      }
    }
  }
}

// =================================================================================================
// constitutive response of each integration point (decoupled)
// =================================================================================================

class Material
{
public:

  // class variables
  // ---------------

  cppmat::matrix<double> eps, sig;
  size_t nip, nelem, ndim=2;

  // class functions
  // ---------------

  Material(size_t nelem, size_t nip);

  void compute_sig();
};

// -------------------------------------------------------------------------------------------------

Material::Material(size_t _nelem, size_t _nip)
{
  // copy from input
  nip   = _nip;
  nelem = _nelem;

  // allocate data
  eps.resize({nelem,nip,ndim,ndim});
  sig.resize({nelem,nip,ndim,ndim});
}

// -------------------------------------------------------------------------------------------------

void Material::compute_sig()
{
  #pragma omp parallel
  {
    cppmat::cartesian2d::tensor2d<double> I = cppmat::cartesian2d::identity2();
    cppmat::cartesian2d::tensor2s<double> Eps, Epsd, Sig;
    double Epsm;
    double K=1., G=1.;

    #pragma omp for
    for ( size_t e = 0 ; e < nelem ; ++e )
    {
      for ( size_t k = 0 ; k < nip ; ++k )
      {
        // - copy to local tensor
        Eps(0,0) = eps(e,k,0,0);
        Eps(0,1) = eps(e,k,0,1);
        Eps(1,1) = eps(e,k,1,1);

        // - decompose strain
        Epsm = Eps.trace() / 2.;
        Epsd = Eps - Epsm * I;

        // - compute stress
        Sig = ( K * Epsm ) * I + G * Epsd;

        // - store in global array
        sig(e,k,0,0) = Sig(0,0);
        sig(e,k,0,1) = Sig(0,1);
        sig(e,k,1,0) = Sig(1,0);
        sig(e,k,1,1) = Sig(1,1);
      }
    }
  }
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
  size_t ndim  = 2;
  size_t nne   = 4;
  size_t nip   = 4;
  size_t nelem = 2;
  size_t nnode = 6;

  MatS conn(nelem,nne );
  MatD coor(nnode,ndim);
  MatS dofs(nnode,ndim);

  coor(0,0) = 0.; coor(0,1) = 0.;
  coor(1,0) = 1.; coor(1,1) = 0.;
  coor(2,0) = 2.; coor(2,1) = 0.;
  coor(3,0) = 0.; coor(3,1) = 1.;
  coor(4,0) = 1.; coor(4,1) = 1.;
  coor(5,0) = 2.; coor(5,1) = 1.;

  dofs(0,0) =  0; dofs(0,1) =  1;
  dofs(1,0) =  2; dofs(1,1) =  3;
  dofs(2,0) =  4; dofs(2,1) =  5;
  dofs(3,0) =  6; dofs(3,1) =  7;
  dofs(4,0) =  8; dofs(4,1) =  9;
  dofs(5,0) = 10; dofs(5,1) = 11;

  conn(0,0) = 0; conn(0,1) = 1; conn(0,2) = 4; conn(0,3) = 3;
  conn(1,0) = 1; conn(1,1) = 2; conn(1,2) = 5; conn(1,3) = 4;

  using Elem = Element<Material>;
  using Sim  = Simulation<Elem>;

  Simulation<Elem> sim(
    std::move(std::make_unique<Elem>(
      std::move(std::make_unique<Material>(nelem,nip)),
      nelem
    )),
    coor,
    conn,
    dofs
  );

  sim.u(0,0) = -0.1; sim.u(3,0) = -0.1;
  sim.u(1,0) =  0.0; sim.u(4,0) =  0.0;
  sim.u(2,0) = +0.1; sim.u(5,0) = +0.1;

  sim.compute_F();

  std::cout << sim.F << std::endl;

}

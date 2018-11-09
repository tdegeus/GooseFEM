
#include <GooseFEM/GooseFEM.h>
#include <xElastoPlasticQPot/ElastoPlasticQPot.h>
#include <HDF5pp.h>

// -------------------------------------------------------------------------------------------------

namespace GF = GooseFEM;
namespace QD = GooseFEM::Element::Quad4;
namespace GM = xElastoPlasticQPot::Cartesian2d;

using T2 = xt::xtensor_fixed<double, xt::xshape<2,2>>;

// -------------------------------------------------------------------------------------------------

inline double sqdot(const xt::xtensor<double,1> &M, const xt::xtensor<double,1> &V)
{
  double out = 0.;

  for ( size_t i = 0 ; i < M.size() ; ++i )
    out += M(i) * V(i) * V(i);

  return out;
}

// -------------------------------------------------------------------------------------------------

class Geometry : public GooseFEM::Dynamics::Geometry
{
private:

  // dimensions
  size_t m_nnode;
  size_t m_ndim;
  size_t m_nelem;
  size_t m_nne;
  size_t m_nip;

  // mesh
  xt::xtensor<size_t,2> m_conn;
  xt::xtensor<double,2> m_coor;
  xt::xtensor<double,2> m_u;
  xt::xtensor<double,2> m_v;
  xt::xtensor<double,2> m_a;

  // vector-definition: transform nodal vectors <-> DOF values
  GF::Vector m_vec;

  // mass matrix (diagonal)
  GF::MatrixDiagonal m_M;

  // numerical quadrature
  QD::Quadrature m_quad;

  // material definition
  GM::Matrix m_mat;

public:

  // constructor
  Geometry(size_t nx);

  // apply update in macroscopic deformation gradient
  void add_dFbar(const T2 &dFbar);

  // compute total kinetic and potential energy
  double Ekin() const;
  double Epot() const;

  // return mesh
  xt::xtensor<size_t,2> conn() const { return m_conn; }
  xt::xtensor<double,2> coor() const { return m_coor; }

  // return nodal vectors
  xt::xtensor<double,2> u() const { return m_u; }
  xt::xtensor<double,2> v() const { return m_v; }
  xt::xtensor<double,2> a() const { return m_a; }

  // return DOF values
  xt::xtensor<double,1> dofs_v() const { return m_vec.asDofs(m_v); }
  xt::xtensor<double,1> dofs_a() const { return m_vec.asDofs(m_a); }

  // overwrite nodal vectors
  void set_u(const xt::xtensor<double,2> &nodevec) { m_u = nodevec; };

  // overwrite nodal vectors, reconstructed from DOF values
  void set_v(const xt::xtensor<double,1> &dofval) { m_v = m_vec.asNode(dofval); };
  void set_a(const xt::xtensor<double,1> &dofval) { m_a = m_vec.asNode(dofval); };

  // solve for DOF-accelerations
  xt::xtensor<double,1> solve_A();
};

// -------------------------------------------------------------------------------------------------

inline Geometry::Geometry(size_t nx)
{
  // get mesh
  GF::Mesh::Quad4::Regular mesh(nx,nx,1.);

  // vector-definition
  m_vec   = GF::Vector(mesh.conn(), mesh.dofsPeriodic());

  // quadrature
  m_quad  = QD::Quadrature(m_vec.asElement(mesh.coor()));

  // dimensions, connectivity, and coordinates
  m_nnode = mesh.nnode();
  m_ndim  = mesh.ndim();
  m_nelem = mesh.nelem();
  m_nne   = mesh.nne();
  m_nip   = m_quad.nip();
  m_conn  = mesh.conn();
  m_coor  = mesh.coor();

  // zero-initialize displacement, velocity, acceleration
  m_u = xt::zeros<double>({m_nnode, m_ndim});
  m_v = xt::zeros<double>({m_nnode, m_ndim});
  m_a = xt::zeros<double>({m_nnode, m_ndim});

  // material definition
  // - allocate
  m_mat = GM::Matrix({m_nelem, m_nip});
  // - phase indicators
  xt::xtensor<size_t,2> Ihard = xt::zeros<size_t>({m_nelem, m_nip});
  xt::xtensor<size_t,2> Isoft = xt::ones <size_t>({m_nelem, m_nip});
  // - set hard indicator
  for ( size_t e = 0 ; e < nx*nx/4 ; ++e )
    for ( size_t k = 0 ; k < m_nip ; ++k )
      Ihard(e,k) = 1;
  // - set soft indicator
  Isoft -= Ihard;
  // - set material definition
  m_mat.setElastic(Ihard, 100., 10.);
  m_mat.setElastic(Isoft, 100.,  1.);
  // - check that all points have been set (not strictly needed)
  m_mat.check();

  // mass matrix
  // - nodal quadrature
  QD::Quadrature q(m_vec.asElement(m_coor), QD::Nodal::xi(), QD::Nodal::w());
  // - set density
  xt::xtensor<double,2> rho = 1.0 * xt::ones<double>({m_nelem, m_nip});
  // - allocate mass matrix
  m_M = GF::MatrixDiagonal(m_conn, mesh.dofsPeriodic());
  // - compute (constant hereafter)
  m_M.assemble( q.int_N_scalar_NT_dV(rho) );
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Geometry::solve_A()
{
  xt::xtensor<double,4> Eps = m_quad.symGradN_vector( m_vec.asElement(m_u) );
  xt::xtensor<double,4> Sig = m_mat.Sig(Eps);
  xt::xtensor<double,1> F   = m_vec.assembleDofs( m_quad.int_gradN_dot_tensor2_dV(Sig) );

  return m_M.solve( -F );
}

// -------------------------------------------------------------------------------------------------

inline void Geometry::add_dFbar(const T2 &dFbar)
{
  for ( size_t i = 0 ; i < m_nnode ; ++i )
    for ( size_t j = 0 ; j < m_ndim ; ++j )
      for ( size_t k = 0 ; k < m_ndim ; ++k )
        m_u(i,j) += dFbar(j,k) * ( m_coor(i,k) - m_coor(0,k) );
}

// -------------------------------------------------------------------------------------------------

inline double Geometry::Ekin() const
{
  xt::xtensor<double,1> V = m_vec.asDofs(m_v);
  xt::xtensor<double,1> M = m_M.asDiagonal();

  return 0.5 * sqdot(M,V);
}

// -------------------------------------------------------------------------------------------------

inline double Geometry::Epot() const
{
  xt::xtensor<double,4> Eps = m_quad.symGradN_vector( m_vec.asElement(m_u) );
  xt::xtensor<double,2> E   = m_mat.energy(Eps);
  xt::xtensor<double,2> dV  = m_quad.dV();

  return ( xt::sum( E*dV ) )[0];
}

// -------------------------------------------------------------------------------------------------

int main()
{
  // set simulation parameters
  double T     = 60.  ; // total time
  double dt    = 1.e-2; // time increment
  size_t nx    = 60   ; // number of elements in both directions
  double gamma = .05  ; // displacement step

  // define geometry
  Geometry geometry(nx);

  // define update in macroscopic deformation gradient
  T2 dFbar = xt::zeros<double>({2,2});
  dFbar(0,1) = gamma;

  // update displacement
  geometry.add_dFbar(dFbar);

  // output variables
  xt::xtensor<double,1> Epot = xt::zeros<double>({static_cast<size_t>(T/dt)});
  xt::xtensor<double,1> Ekin = xt::zeros<double>({static_cast<size_t>(T/dt)});
  xt::xtensor<double,1> t    = xt::zeros<double>({static_cast<size_t>(T/dt)});

  // loop over increments
  for ( size_t inc = 0 ; inc < static_cast<size_t>(Epot.size()) ; ++inc )
  {
    // - compute increment
    GF::Dynamics::velocityVerlet(geometry, dt);

    // - store output
    t   (inc) = static_cast<double>(inc) * dt;
    Ekin(inc) = geometry.Ekin();
    Epot(inc) = geometry.Epot();
  }

  // write to output file
  H5p::File f = H5p::File("example.hdf5","w");
  f.write("/global/Epot",Epot           );
  f.write("/global/Ekin",Ekin           );
  f.write("/global/t"   ,t              );
  f.write("/mesh/conn"  ,geometry.conn());
  f.write("/mesh/coor"  ,geometry.coor());
  f.write("/mesh/disp"  ,geometry.u()   );

  return 0;
}

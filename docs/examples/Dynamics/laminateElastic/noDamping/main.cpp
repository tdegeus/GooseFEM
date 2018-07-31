
#include <GooseFEM/GooseFEM.h>
#include <ElastoPlasticQPot/ElastoPlasticQPot.h>
// #include <HDF5pp.h>

// -------------------------------------------------------------------------------------------------

namespace GF = GooseFEM;
namespace QD = GooseFEM::Element::Quad4;
namespace GM = ElastoPlasticQPot::Cartesian2d;

using T2 = cppmat::tiny::cartesian::tensor2<double,2>;

// -------------------------------------------------------------------------------------------------

inline double sqdot(const GF::ColD &M, const GF::ColD &V)
{
  double out = 0.;

  for ( auto i = 0 ; i < M.size() ; ++i )
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
  GF::MatS m_conn;
  GF::MatD m_coor;
  GF::MatD m_u;
  GF::MatD m_v;
  GF::MatD m_a;

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
  GF::MatS conn() const { return m_conn; }
  GF::MatD coor() const { return m_coor; }

  // return nodal vectors
  GF::MatD u() const { return m_u; }
  GF::MatD v() const { return m_v; }
  GF::MatD a() const { return m_a; }

  // return DOF values
  GF::ColD dofs_v() const { return m_vec.asDofs(m_v); }
  GF::ColD dofs_a() const { return m_vec.asDofs(m_a); }

  // overwrite nodal vectors
  void set_u(const GF::MatD &nodevec) { m_u = nodevec; };

  // overwrite nodal vectors, reconstructed from DOF values
  void set_v(const GF::ColD &dofval) { m_v = m_vec.asNode(dofval); };
  void set_a(const GF::ColD &dofval) { m_a = m_vec.asNode(dofval); };

  // solve for DOF-accelerations
  GF::ColD solve_A();
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
  m_u = GF::MatD::Zero(m_nnode, m_ndim);
  m_v = GF::MatD::Zero(m_nnode, m_ndim);
  m_a = GF::MatD::Zero(m_nnode, m_ndim);

  // material definition
  // - allocate
  m_mat = GM::Matrix({m_nelem, m_nip});
  // - phase indicators
  GM::ArrS Ihard = GM::ArrS::Zero({m_nelem, m_nip});
  GM::ArrS Isoft = GM::ArrS::Ones({m_nelem, m_nip});
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
  GF::ArrD rho = 1.0 * GF::ArrD::Ones({m_nelem, m_nip});
  // - allocate mass matrix
  m_M = GF::MatrixDiagonal(m_conn, mesh.dofsPeriodic());
  // - compute (constant hereafter)
  m_M.assemble( q.int_N_scalar_NT_dV(rho) );
}

// -------------------------------------------------------------------------------------------------

inline GF::ColD Geometry::solve_A()
{
  GF::ArrD Eps = m_quad.symGradN_vector( m_vec.asElement(m_u) );
  GF::ArrD Sig = m_mat.Sig(Eps);
  GF::ColD F = m_vec.assembleDofs( m_quad.int_gradN_dot_tensor2_dV(Sig) );

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
  GF::ColD V = m_vec.asDofs(m_v);
  GF::ColD M = m_M.asDiagonal();

  return 0.5 * sqdot(M,V);
}

// -------------------------------------------------------------------------------------------------

inline double Geometry::Epot() const
{
  GF::ArrD Eps = m_quad.symGradN_vector( m_vec.asElement(m_u) );
  GF::ArrD E   = m_mat.energy(Eps);
  GF::ArrD dV  = m_quad.dV();

  return E.average(dV,false);
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
  T2 dFbar = T2::Zero();
  dFbar(0,1) = gamma;

  // update displacement
  geometry.add_dFbar(dFbar);

  // output variables
  GF::ColD Epot = GF::ColD::Zero(static_cast<int>(T/dt));
  GF::ColD Ekin = GF::ColD::Zero(static_cast<int>(T/dt));
  GF::ColD t    = GF::ColD::Zero(static_cast<int>(T/dt));

  // loop over increments
  for ( size_t inc = 0 ; inc < static_cast<size_t>(Epot.size()) ; ++inc )
  // for ( size_t inc = 0 ; inc < 2000 ; ++inc )
  {
    // - compute increment
    GF::Dynamics::velocityVerlet(geometry, dt);

    // - store output
    t   (inc) = static_cast<double>(inc) * dt;
    Ekin(inc) = geometry.Ekin();
    Epot(inc) = geometry.Epot();
  }

  // // write to output file
  // H5p::File f = H5p::File("example.hdf5","w");
  // f.write("/global/Epot",Epot           );
  // f.write("/global/Ekin",Ekin           );
  // f.write("/global/t"   ,t              );
  // f.write("/mesh/conn"  ,geometry.conn());
  // f.write("/mesh/coor"  ,geometry.coor());
  // f.write("/mesh/disp"  ,geometry.u()   );

  return 0;
}

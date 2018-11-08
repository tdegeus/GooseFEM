#include <xGooseFEM/GooseFEM.h>
#include <GMatLinearElastic/GMatLinearElastic.h>

int main()
{
  // mesh
  // ----

  xGooseFEM::Mesh::Quad4::Regular mesh(5,5);

  xt::xtensor<double,2> coor = mesh.coor();
  xt::xtensor<size_t,2> conn = mesh.conn();
  xt::xtensor<size_t,2> dofs = mesh.dofs();
  xt::xtensor<double,2> disp = xt::zeros<double>(coor.shape());

  xt::xtensor<size_t,1> nodesLeft   = mesh.nodesLeftEdge();
  xt::xtensor<size_t,1> nodesRight  = mesh.nodesRightEdge();
  xt::xtensor<size_t,1> nodesTop    = mesh.nodesTopEdge();
  xt::xtensor<size_t,1> nodesBottom = mesh.nodesBottomEdge();

  // boundary displacements
  // ----------------------

  xt::xtensor<size_t,1> iip = xt::empty<size_t>({2*nodesLeft.size()+2*nodesBottom.size()});
  xt::xtensor<double,1> u_p = xt::zeros<double>(iip.shape());

  size_t i = 0;
  for ( auto &n : nodesRight  ) { iip(i) = dofs(n,0); u_p(i) =  0.1; ++i; }
  for ( auto &n : nodesTop    ) { iip(i) = dofs(n,1); u_p(i) = -0.1; ++i; }
  for ( auto &n : nodesLeft   ) { iip(i) = dofs(n,0); u_p(i) =  0.0; ++i; }
  for ( auto &n : nodesBottom ) { iip(i) = dofs(n,1); u_p(i) =  0.0; ++i; }

  // element definition
  // ------------------

  xGooseFEM::VectorPartitioned vector(conn, dofs, iip);
  xGooseFEM::MatrixPartitioned K     (conn, dofs, iip);

  xGooseFEM::Element::Quad4::QuadraturePlanar elem(vector.asElement(coor));

  GMatLinearElastic::Cartesian3d::Material mat(mesh.nelem(), elem.nip(), 1., 1.);

  // solve
  // -----

  // strain
  xt::xtensor<double,4> Eps = elem.symGradN_vector(vector.asElement(disp));

  // stress
  xt::xtensor<double,4> Sig = mat.Sig(Eps);

  // tangent
  xt::xtensor<double,6> C = mat.Tangent();

  // internal force
  xt::xtensor<double,2> fint = vector.asNode(elem.int_gradN_dot_tensor2_dV(Sig));

  // stiffness matrix
  K.assemble(elem.int_gradN_dot_tensor4_dot_gradNT_dV(C));

  // solve
  xt::xtensor<double,1> u_u = K.solve(vector.asDofs_u(fint), u_p);

  // assemble to nodal vector
  disp = vector.asNode(u_u, u_p);

  // print result
  std::cout << disp << std::endl;

  return 0;
}

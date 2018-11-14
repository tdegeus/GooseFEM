#include <GooseFEM/GooseFEM.h>
#include <GMatLinearElastic/Cartesian3d.h>

int main()
{
  // mesh
  // ----

  GooseFEM::Mesh::Quad4::Regular mesh(5,5);

  xt::xtensor<double,2> coor = mesh.coor();
  xt::xtensor<size_t,2> conn = mesh.conn();
  xt::xtensor<size_t,2> dofs = mesh.dofs();
  xt::xtensor<double,2> disp = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fint = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fext = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fres = xt::zeros<double>(coor.shape());

  xt::xtensor<size_t,1> nodesLeft   = mesh.nodesLeftEdge();
  xt::xtensor<size_t,1> nodesRight  = mesh.nodesRightEdge();
  xt::xtensor<size_t,1> nodesTop    = mesh.nodesTopEdge();
  xt::xtensor<size_t,1> nodesBottom = mesh.nodesBottomEdge();

  // fixed displacements DOFs & displacements
  // ----------------------------------------

  xt::xtensor<size_t,1> iip = xt::empty<size_t>({2*nodesLeft.size()+2*nodesBottom.size()});
  xt::xtensor<double,1> u_p = xt::zeros<double>(iip.shape());

  {
    size_t i = 0;
    for ( auto &n : nodesRight  ) { iip(i) = dofs(n,0); u_p(i) =  0.1; ++i; }
    for ( auto &n : nodesTop    ) { iip(i) = dofs(n,1); u_p(i) = -0.1; ++i; }
    for ( auto &n : nodesLeft   ) { iip(i) = dofs(n,0); u_p(i) =  0.0; ++i; }
    for ( auto &n : nodesBottom ) { iip(i) = dofs(n,1); u_p(i) =  0.0; ++i; }
  }

  // element definition
  // ------------------

  xt::xtensor<double,4> Eps, Sig;
  xt::xtensor<double,6> C;

  GooseFEM::VectorPartitioned vector(conn, dofs, iip);
  GooseFEM::MatrixPartitioned K     (conn, dofs, iip);

  GooseFEM::Element::Quad4::QuadraturePlanar elem(vector.asElement(coor));

  GMatLinearElastic::Cartesian3d::Matrix mat(mesh.nelem(), elem.nip(), 1., 1.);

  // solve
  // -----

  // strain
  Eps = elem.symGradN_vector(vector.asElement(disp));

  // stress & tangent
  std::tie(Sig,C) = mat.Tangent(Eps);

  // internal force
  fint = vector.assembleNode(elem.int_gradN_dot_tensor2_dV(Sig));

  // stiffness matrix
  K.assemble(elem.int_gradN_dot_tensor4_dot_gradNT_dV(C));

  // compute residual
  xt::noalias(fres) = fext - fint;

  // solve
  xt::xtensor<double,1> u_u = K.solve_u(vector.asDofs_u(fres), u_p);

  // assemble to nodal vector
  disp = vector.asNode(u_u, u_p);

  // print result
  std::cout << disp << std::endl;

  return 0;
}

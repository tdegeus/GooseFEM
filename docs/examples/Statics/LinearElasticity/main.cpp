#include <GooseFEM/GooseFEM.h>
#include <GMatLinearElastic/Cartesian3d.h>
#include <LowFive.h>

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

  // fixed displacements DOFs
  // ------------------------

  xt::xtensor<size_t,1> iip = xt::empty<size_t>({2*nodesLeft.size()+2*nodesBottom.size()});

  {
    size_t i = 0;
    for ( auto &n : nodesRight  ) { iip(i) = dofs(n,0); ++i; }
    for ( auto &n : nodesTop    ) { iip(i) = dofs(n,1); ++i; }
    for ( auto &n : nodesLeft   ) { iip(i) = dofs(n,0); ++i; }
    for ( auto &n : nodesBottom ) { iip(i) = dofs(n,1); ++i; }
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

  // set fixed displacements
  for ( auto &n : nodesRight  ) disp(n,0) =  0.1;
  for ( auto &n : nodesTop    ) disp(n,1) = -0.1;
  for ( auto &n : nodesLeft   ) disp(n,0) =  0.0;
  for ( auto &n : nodesBottom ) disp(n,1) =  0.0;

  // compute residual
  xt::noalias(fres) = fext - fint;

  // solve
  K.solve(fres, disp);

  // print result
  std::cout << disp << std::endl;

  // store result

  Eps = elem.symGradN_vector(vector.asElement(disp));
  Sig = mat.Sig(Eps);

  xt::xtensor<double,4> dV = elem.dVtensor();
  xt::xtensor<double,3> SigBar = xt::average(Sig, dV, {1});

  HighFive::File file("main.h5", HighFive::File::Overwrite);

  LowFive::xtensor::dump(file, "/coor", coor);
  LowFive::xtensor::dump(file, "/conn", conn);
  LowFive::xtensor::dump(file, "/disp", disp);
  LowFive::xtensor::dump(file, "/Sig" , SigBar);

  return 0;
}

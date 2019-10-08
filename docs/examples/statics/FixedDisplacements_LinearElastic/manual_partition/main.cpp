#include <GooseFEM/GooseFEM.h>
#include <GooseFEM/MatrixPartitioned.h>
#include <GMatElastic/Cartesian3d.h>
#include <highfive/H5Easy.hpp>

int main()
{
  // mesh
  // ----

  // define mesh
  GooseFEM::Mesh::Quad4::Regular mesh(5,5);

  // mesh dimensions
  size_t nelem = mesh.nelem();
  size_t nne   = mesh.nne();
  size_t ndim  = mesh.ndim();

  // mesh definitions
  xt::xtensor<double,2> coor = mesh.coor();
  xt::xtensor<size_t,2> conn = mesh.conn();
  xt::xtensor<size_t,2> dofs = mesh.dofs();

  // node sets
  xt::xtensor<size_t,1> nodesLeft   = mesh.nodesLeftEdge();
  xt::xtensor<size_t,1> nodesRight  = mesh.nodesRightEdge();
  xt::xtensor<size_t,1> nodesTop    = mesh.nodesTopEdge();
  xt::xtensor<size_t,1> nodesBottom = mesh.nodesBottomEdge();

  // fixed displacements DOFs
  // ------------------------

  xt::xtensor<size_t,1> iip = xt::concatenate(xt::xtuple(
    xt::view(dofs, xt::keep(nodesRight ), 0),
    xt::view(dofs, xt::keep(nodesTop   ), 1),
    xt::view(dofs, xt::keep(nodesLeft  ), 0),
    xt::view(dofs, xt::keep(nodesBottom), 1)
  ));

  // simulation variables
  // --------------------

  // vector definition
  GooseFEM::VectorPartitioned vector(conn, dofs, iip);

  // allocate system matrix
  GooseFEM::MatrixPartitioned<> K(conn, dofs, iip);

  // nodal quantities
  xt::xtensor<double,2> disp = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fint = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fext = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fres = xt::zeros<double>(coor.shape());

  // DOF values
  xt::xtensor<double,1> u_u    = xt::zeros<double>({vector.nnu()});
  xt::xtensor<double,1> fres_u = xt::zeros<double>({vector.nnu()});
  xt::xtensor<double,1> fext_p = xt::zeros<double>({vector.nnp()});

  // element vectors
  xt::xtensor<double,3> ue = xt::empty<double>({nelem, nne, ndim});
  xt::xtensor<double,3> fe = xt::empty<double>({nelem, nne, ndim});
  xt::xtensor<double,3> Ke = xt::empty<double>({nelem, nne*ndim, nne*ndim});

  // element/material definition
  // ---------------------------

  // element definition
  GooseFEM::Element::Quad4::QuadraturePlanar elem(vector.AsElement(coor));
  size_t nip = elem.nip();

  // material definition
  GMatElastic::Cartesian3d::Matrix mat(nelem, nip, 1., 1.);

  // integration point tensors
  size_t d = 3;
  xt::xtensor<double,4> Eps = xt::empty<double>({nelem, nip, d, d      });
  xt::xtensor<double,4> Sig = xt::empty<double>({nelem, nip, d, d      });
  xt::xtensor<double,6> C   = xt::empty<double>({nelem, nip, d, d, d, d});

  // solve
  // -----

  // strain
  vector.asElement(disp, ue);
  elem.symGradN_vector(ue, Eps);

  // stress & tangent
  mat.tangent(Eps, Sig, C);

  // internal force
  elem.int_gradN_dot_tensor2_dV(Sig, fe);
  vector.assembleNode(fe, fint);

  // stiffness matrix
  elem.int_gradN_dot_tensor4_dot_gradNT_dV(C, Ke);
  K.assemble(Ke);

  // set fixed displacements
  xt::xtensor<double,1> u_p = xt::concatenate(xt::xtuple(
    +0.1 * xt::ones<double>({nodesRight .size()}),
    -0.1 * xt::ones<double>({nodesTop   .size()}),
     0.0 * xt::ones<double>({nodesLeft  .size()}),
     0.0 * xt::ones<double>({nodesBottom.size()})
  ));

  // residual
  xt::noalias(fres) = fext - fint;

  // partition
  vector.asDofs_u(fres, fres_u);

  // solve
  K.solve_u(fres_u, u_p, u_u);

  // assemble to nodal vector
  vector.asNode(u_u, u_p, disp);

  // post-process
  // ------------

  // compute strain and stress
  vector.asElement(disp, ue);
  elem.symGradN_vector(ue, Eps);
  mat.stress(Eps, Sig);

  // internal force
  elem.int_gradN_dot_tensor2_dV(Sig, fe);
  vector.assembleNode(fe, fint);

  // apply reaction force
  vector.asDofs_p(fint, fext_p);

  // residual
  xt::noalias(fres) = fext - fint;

  // partition
  vector.asDofs_u(fres, fres_u);

  // print residual
  std::cout << xt::sum(xt::abs(fres_u))[0] / xt::sum(xt::abs(fext_p))[0] << std::endl;

  // average stress per node
  xt::xtensor<double,4> dV = elem.DV(2);
  xt::xtensor<double,3> SigAv = xt::average(Sig, dV, {1});

  // write output
  H5Easy::File file("main.h5", H5Easy::File::Overwrite);
  H5Easy::dump(file, "/coor", coor);
  H5Easy::dump(file, "/conn", conn);
  H5Easy::dump(file, "/disp", disp);
  H5Easy::dump(file, "/Sig" , SigAv);

  return 0;
}

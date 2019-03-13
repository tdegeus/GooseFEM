#include <GooseFEM/GooseFEM.h>
#include <GooseFEM/TyingsPeriodic.h>
#include <GooseFEM/MatrixPartitionedTyings.h>
#include <GooseFEM/VectorPartitionedTyings.h>
#include <GMatLinearElastic/Cartesian3d.h>
#include <xtensor-io/xhighfive.hpp>

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

  // periodicity and fixed displacements DOFs
  // ----------------------------------------

  // add control nodes/DOFs

  xt::xtensor<size_t,2> control_dofs = xt::arange<size_t>(ndim*ndim).reshape({ndim,ndim});

  control_dofs += xt::amax(dofs)[0] + 1;

  xt::xtensor<size_t,1> control_nodes = mesh.nnode() + xt::arange(ndim);

  coor = xt::concatenate(xt::xtuple(coor, xt::zeros<double>({ndim,ndim})));

  dofs = xt::concatenate(xt::xtuple(dofs, control_dofs));

  // extract fixed DOFs

  xt::xtensor<size_t,1> iip = xt::concatenate(xt::xtuple(
    xt::reshape_view(control_dofs, {ndim*ndim}),
    xt::reshape_view(xt::view(dofs, xt::keep(mesh.nodesOrigin()), xt::all()), {ndim})
  ));

  // get DOF-tyings, reorganise system

  GooseFEM::Tyings::Periodic tyings(coor, dofs, control_dofs, mesh.nodesPeriodic(), iip);

  dofs = tyings.dofs();

  // simulation variables
  // --------------------

  // vector definition

  GooseFEM::VectorPartitionedTyings vector(conn, dofs, tyings.Cdu(), tyings.Cdp(), tyings.Cdi());

  // nodal quantities

  xt::xtensor<double,2> disp = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fint = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fext = xt::zeros<double>(coor.shape());
  xt::xtensor<double,2> fres = xt::zeros<double>(coor.shape());

  // element vectors

  xt::xtensor<double,3> ue = xt::empty<double>({nelem, nne, ndim});
  xt::xtensor<double,3> fe = xt::empty<double>({nelem, nne, ndim});
  xt::xtensor<double,3> Ke = xt::empty<double>({nelem, nne*ndim, nne*ndim});

  // element/material definition
  // ---------------------------

  GooseFEM::Element::Quad4::QuadraturePlanar elem(vector.AsElement(coor));

  size_t nip = elem.nip();

  GMatLinearElastic::Cartesian3d::Matrix mat(nelem, nip);

  xt::xtensor<size_t,2> Ihard = xt::zeros<size_t>({nelem, nip});

  xt::view(Ihard, xt::keep(0,1,5,6), xt::all()) = 1;

  xt::xtensor<size_t,2> Isoft = xt::ones<size_t>({nelem, nip}) - Ihard;

  mat.set(Isoft,10.,1. );
  mat.set(Ihard,10.,10.);

  // solve
  // -----

  // allocate tensors
  size_t d = 3;
  xt::xtensor<double,4> Eps = xt::empty<double>({nelem, nip, d, d      });
  xt::xtensor<double,4> Sig = xt::empty<double>({nelem, nip, d, d      });
  xt::xtensor<double,6> C   = xt::empty<double>({nelem, nip, d, d, d, d});

  // allocate system matrix
  GooseFEM::MatrixPartitionedTyings K(conn, dofs, tyings.Cdu(), tyings.Cdp());

  // strain
  vector.asElement(disp, ue);
  elem.symGradN_vector(ue, Eps);

  // stress & tangent
  mat.Tangent(Eps, Sig, C);

  // internal force
  elem.int_gradN_dot_tensor2_dV(Sig, fe);
  vector.assembleNode(fe, fint);

  // stiffness matrix
  elem.int_gradN_dot_tensor4_dot_gradNT_dV(C, Ke);
  K.assemble(Ke);

  // set fixed displacements
  disp(control_nodes(0),1) = 0.1;

  // residual
  xt::noalias(fres) = fext - fint;

  // solve
  K.solve(fres, disp);

  // post-process
  // ------------

  // compute strain and stress
  vector.asElement(disp, ue);
  elem.symGradN_vector(ue, Eps);
  mat.Sig(Eps, Sig);

  // internal force
  elem.int_gradN_dot_tensor2_dV(Sig, fe);
  vector.assembleNode(fe, fint);

  // allocate DOF-list
  xt::xtensor<double,1> Fext = xt::zeros<double>({tyings.nni()});
  xt::xtensor<double,1> Fint = xt::zeros<double>({tyings.nni()});

  // internal/external force on DOFs
  vector.asDofs_i(fext, Fext);
  vector.asDofs_i(fint, Fint);

  // apply reaction force
  vector.copy_p(Fint, Fext);

  // print residual
  std::cout << xt::sum(xt::abs(Fext-Fint))[0] / xt::sum(xt::abs(Fext))[0] << std::endl;

  // average stress per node
  xt::xtensor<double,4> dV = elem.DV(2);
  xt::xtensor<double,3> SigAv = xt::average(Sig, dV, {1});

  // write output

  HighFive::File file("main.h5", HighFive::File::Overwrite);

  xt::dump(file, "/coor", coor);
  xt::dump(file, "/conn", conn);
  xt::dump(file, "/disp", disp);
  xt::dump(file, "/Sig" , SigAv);

  return 0;
}

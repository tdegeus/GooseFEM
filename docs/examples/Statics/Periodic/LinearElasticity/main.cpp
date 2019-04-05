#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>
#include <GMatElastic/Cartesian3d.h>
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

  // add control nodes
  GooseFEM::Tyings::Control control(coor, dofs);
  coor = control.coor();
  dofs = control.dofs();
  xt::xtensor<size_t,2> control_dofs = control.controlDofs();
  xt::xtensor<size_t,1> control_nodes = control.controlNodes();

  // extract fixed DOFs:
  // - all control nodes: to prescribe the deformation gradient
  // - one node of the mesh: to remove rigid body modes
  xt::xtensor<size_t,1> iip = xt::concatenate(xt::xtuple(
    xt::reshape_view(control_dofs, {ndim*ndim}),
    xt::reshape_view(xt::view(dofs, xt::keep(mesh.nodesOrigin()), xt::all()), {ndim})
  ));

  // get DOF-tyings, reorganise system
  GooseFEM::Tyings::Periodic tyings(coor, dofs, control_dofs, mesh.nodesPeriodic(), iip);
  dofs = tyings.dofs();

  // simulation variables
  // --------------------

  // vector definition:
  // provides methods to switch between dofval/nodeval/elemvec, or to manipulate a part of them
  GooseFEM::VectorPartitionedTyings vector(conn, dofs, tyings.Cdu(), tyings.Cdp(), tyings.Cdi());

  // nodal quantities
  xt::xtensor<double,2> disp = xt::zeros<double>(coor.shape()); // nodal displacement
  xt::xtensor<double,2> fint = xt::zeros<double>(coor.shape()); // internal force
  xt::xtensor<double,2> fext = xt::zeros<double>(coor.shape()); // external force
  xt::xtensor<double,2> fres = xt::zeros<double>(coor.shape()); // residual force

  // element vectors / matrix
  xt::xtensor<double,3> ue = xt::empty<double>({nelem, nne, ndim});
  xt::xtensor<double,3> fe = xt::empty<double>({nelem, nne, ndim});
  xt::xtensor<double,3> Ke = xt::empty<double>({nelem, nne*ndim, nne*ndim});

  // element/material definition
  // ---------------------------

  // FEM quadrature
  GooseFEM::Element::Quad4::QuadraturePlanar elem(vector.AsElement(coor));
  size_t nip = elem.nip();

  // material model
  // even though the problem is 2-d, the material model is 3-d, plane strain is implicitly assumed
  GMatElastic::Cartesian3d::Matrix mat(nelem, nip);
  size_t tdim = mat.ndim();

  // some artificial material definition
  xt::xtensor<size_t,2> Ihard = xt::zeros<size_t>({nelem, nip});
  xt::view(Ihard, xt::keep(0,1,5,6), xt::all()) = 1;
  xt::xtensor<size_t,2> Isoft = xt::ones<size_t>({nelem, nip}) - Ihard;

  mat.setElastic(Isoft, 10.0,  1.0);
  mat.setElastic(Ihard, 10.0, 10.0);

  // solve
  // -----

  // allocate tensors
  xt::xtensor<double,4> Eps = xt::empty<double>({nelem, nip, tdim, tdim});
  xt::xtensor<double,4> Sig = xt::empty<double>({nelem, nip, tdim, tdim});
  xt::xtensor<double,6> C   = xt::empty<double>({nelem, nip, tdim, tdim, tdim, tdim});

  // allocate system matrix
  GooseFEM::MatrixPartitionedTyings K(conn, dofs, tyings.Cdu(), tyings.Cdp());

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


  // residual
  xt::noalias(fres) = fext - fint;

  // set fixed displacements
  disp(control_nodes(0),1) = 0.1;
  
  // solve
  K.solve(fres, disp);

  // post-process
  // ------------

  // compute strain and stress
  vector.asElement(disp, ue);
  elem.symGradN_vector(ue, Eps);
  mat.stress(Eps, Sig);

  // average stress per node
  xt::xtensor<double,4> dV = elem.DV(2);
  xt::xtensor<double,3> SigAv = xt::average(Sig, dV, {1});

  // write output
  HighFive::File file("main.h5", HighFive::File::Overwrite);
  xt::dump(file, "/coor", coor);
  xt::dump(file, "/conn", conn);
  xt::dump(file, "/disp", disp);
  xt::dump(file, "/Sig", SigAv);

  return 0;
}

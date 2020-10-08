#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>
#include <GooseFEM/ParaView.h>
#include <GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h>
#include <highfive/H5Easy.hpp>

namespace GM = GMatElastoPlasticFiniteStrainSimo::Cartesian3d;
namespace GF = GooseFEM;
namespace PV = GooseFEM::ParaView::HDF5;
namespace H5 = H5Easy;

int main()
{
  // mesh
  // ----

  // define mesh
  GF::Mesh::Quad4::Regular mesh(5, 5);

  // mesh dimensions
  size_t nelem = mesh.nelem();
  size_t nne   = mesh.nne();
  size_t ndim  = mesh.ndim();

  // mesh definitions
  xt::xtensor<double,2> coor  = mesh.coor();
  xt::xtensor<size_t,2> conn  = mesh.conn();
  xt::xtensor<size_t,2> dofs  = mesh.dofs();
  xt::xtensor<size_t,2> elmat = mesh.elementMatrix();

  // periodicity and fixed displacements DOFs
  // ----------------------------------------

  // add control nodes
  GF::Tyings::Control control(coor, dofs);
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
  GF::Tyings::Periodic tyings(coor, dofs, control_dofs, mesh.nodesPeriodic(), iip);
  dofs = tyings.dofs();

  // simulation variables
  // --------------------

  // vector definition:
  // provides methods to switch between dofval/nodeval/elemvec, or to manipulate a part of them
  GF::VectorPartitionedTyings vector(conn, dofs, tyings.Cdu(), tyings.Cdp(), tyings.Cdi());

  // nodal quantities
  xt::xtensor<double,2> disp = xt::zeros<double>(coor.shape()); // nodal displacement
  xt::xtensor<double,2> du   = xt::zeros<double>(coor.shape()); // iterative displacement update
  xt::xtensor<double,2> fint = xt::zeros<double>(coor.shape()); // internal force
  xt::xtensor<double,2> fext = xt::zeros<double>(coor.shape()); // external force
  xt::xtensor<double,2> fres = xt::zeros<double>(coor.shape()); // residual force

  // element vectors / matrix
  xt::xtensor<double,3> ue = xt::empty<double>({nelem, nne, ndim});
  xt::xtensor<double,3> fe = xt::empty<double>({nelem, nne, ndim});
  xt::xtensor<double,3> Ke = xt::empty<double>({nelem, nne*ndim, nne*ndim});

  // DOF values
  xt::xtensor<double,1> Fext = xt::zeros<double>({tyings.nni()});
  xt::xtensor<double,1> Fint = xt::zeros<double>({tyings.nni()});

  // element/material definition
  // ---------------------------

  // FEM quadrature
  GF::Element::Quad4::QuadraturePlanar elem0(vector.AsElement(coor));
  GF::Element::Quad4::QuadraturePlanar elem (vector.AsElement(coor));
  size_t nip = elem0.nip();

  // material model
  // even though the problem is 2-d, the material model is 3-d, plane strain is implicitly assumed
  GM::Matrix mat(nelem, nip);
  size_t tdim = mat.ndim();

  // some artificial material definition
  xt::xtensor<size_t,1> ehard = xt::ravel(xt::view(elmat, xt::range(0,2), xt::range(0,2)));
  xt::xtensor<size_t,2> Ihard = xt::zeros<size_t>({nelem, nip});
  xt::view(Ihard, xt::keep(ehard), xt::all()) = 1ul;
  xt::xtensor<size_t,2> Isoft = xt::ones<size_t>({nelem, nip}) - Ihard;

  mat.setLinearHardening(Isoft, 1.0, 1.0, 0.05, 0.05);
  mat.setElastic(Ihard, 1.0, 1.0);

  // solve
  // -----

  // allocate tensors
  xt::xtensor<double,4> I2  = mat.I2();
  xt::xtensor<double,4> F   = xt::empty<double>({nelem, nip, tdim, tdim});
  xt::xtensor<double,4> Eps = xt::empty<double>({nelem, nip, tdim, tdim});
  xt::xtensor<double,4> Sig = xt::empty<double>({nelem, nip, tdim, tdim});
  xt::xtensor<double,6> C   = xt::empty<double>({nelem, nip, tdim, tdim, tdim, tdim});

  // allocate system matrix
  GF::MatrixPartitionedTyings K(conn, dofs, tyings.Cdu(), tyings.Cdp());
  GF::MatrixPartitionedTyingsSolver<> Solver;

  // allocate internal variables
  double res;

  // some shear strain history
  xt::xtensor<double,1> dgamma = 0.001 * xt::ones<double>({101});
  dgamma(0) = 0.0;

  // initialise output
  // - output-file containing data
  HighFive::File file("main.h5", HighFive::File::Overwrite);
  // - ParaView meta-data
  PV::TimeSeries xdmf;
  // - write mesh
  H5::dump(file, "/conn", conn);
  H5::dump(file, "/coor", coor);

  // loop over increments
  for (size_t inc = 0; inc < dgamma.size(); ++inc)
  {
    // update history
    mat.increment();

    // iterate
    for (size_t iter = 0; ; ++iter)
    {
      // deformation gradient tensor
      vector.asElement(disp, ue);
      elem0.gradN_vector_T(ue, F);
      F += I2;

      // stress & tangent
      mat.tangent(F, Sig, C);

      // internal force
      elem.int_gradN_dot_tensor2_dV(Sig, fe);
      vector.assembleNode(fe, fint);

      // stiffness matrix
      elem.int_gradN_dot_tensor4_dot_gradNT_dV(C, Ke);
      K.assemble(Ke);

      // residual
      xt::noalias(fres) = fext - fint;

      // check for convergence (skip the zeroth iteration, as the residual still vanishes)
      if (iter > 0)
      {
        // - internal/external force as DOFs (account for periodicity)
        vector.asDofs_i(fext, Fext);
        vector.asDofs_i(fint, Fint);
        // - extract reaction force
        vector.copy_p(Fint, Fext);
        // - norm of the residual and the reaction force
        double nfres = xt::sum(xt::abs(Fext-Fint))[0];
        double nfext = xt::sum(xt::abs(Fext))[0];
        // - relative residual, for convergence check
        if (nfext)
          res = nfres / nfext;
        else
          res = nfres;
        // - print progress to screen
        std::cout << "inc = " << inc << ", iter = " << iter << ", res = " << res << std::endl;
        // - check for convergence
        if (res < 1.0e-5)
          break;
        // - safe-guard from infinite loop
        if (iter > 20)
          throw std::runtime_error("Maximal number of iterations exceeded");
      }

      // initialise displacement update
      du.fill(0.0);

      // set fixed displacements
      if (iter == 0)
        du(control_nodes(0),1) = dgamma(inc);

      // solve
      Solver.solve(K, fres, du);

      // add displacement update
      disp += du;

      // update shape functions
      elem.update_x(vector.AsElement(coor+disp));
    }

    // post-process
    // - compute strain and stress
    vector.asElement(disp, ue);
    elem0.gradN_vector_T(ue, F);
    F += I2;
    GM::strain(F, Eps);
    mat.stress(F, Sig);
    // - integration point volume
    xt::xtensor<double,4> dV = elem.AsTensor<2>(elem.dV());
    // - element average stress
    xt::xtensor<double,3> Sigelem = xt::average(Sig, dV, {1});
    xt::xtensor<double,3> Epselem = xt::average(Eps, dV, {1});
    // - macroscopic strain and stress
    xt::xtensor_fixed<double, xt::xshape<3,3>> Sigbar = xt::average(Sig, dV, {0,1});
    xt::xtensor_fixed<double, xt::xshape<3,3>> Epsbar = xt::average(Eps, dV, {0,1});
    // - write to output-file: increment numbers
    H5::dump(file, "/stored", inc, {inc});
    // - write to output-file: macroscopic response
    H5::dump(file, "/macroscopic/sigeq", GM::Sigeq(Sigbar), {inc});
    H5::dump(file, "/macroscopic/epseq", GM::Epseq(Epsbar), {inc});
    // - write to output-file: element quantities
    H5::dump(file, "/sigeq/" + std::to_string(inc), GM::Sigeq(Sigelem));
    H5::dump(file, "/epseq/" + std::to_string(inc), GM::Epseq(Epselem));
    H5::dump(file, "/disp/"  + std::to_string(inc), PV::as3d(disp));
    // - update ParaView meta-data
    xdmf.push_back(PV::Increment(
      PV::Connectivity(file, "/conn", mesh.getElementType()),
      PV::Coordinates (file, "/coor"),
      {
        PV::Attribute(file, "/disp/"  + std::to_string(inc), "Displacement", PV::AttributeType::Node),
        PV::Attribute(file, "/sigeq/" + std::to_string(inc), "Eq. stress"  , PV::AttributeType::Cell),
        PV::Attribute(file, "/epseq/" + std::to_string(inc), "Eq. strain"  , PV::AttributeType::Cell),
      }
    ));
  }

  // write ParaView meta-data
  xdmf.write("main.xdmf");

  return 0;
}

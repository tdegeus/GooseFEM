#include <GMatElastic/Cartesian3d.h>
#include <GooseFEM/GooseFEM.h>
#include <GooseFEM/MatrixPartitioned.h>
#include <highfive/H5Easy.hpp>

int main()
{
    // mesh
    // ----

    // define mesh
    GooseFEM::Mesh::Quad4::Regular mesh(5, 5);

    // mesh dimensions
    size_t nelem = mesh.nelem();
    size_t nne = mesh.nne();
    size_t ndim = mesh.ndim();

    // mesh definitions
    xt::xtensor<double, 2> coor = mesh.coor();
    xt::xtensor<size_t, 2> conn = mesh.conn();
    xt::xtensor<size_t, 2> dofs = mesh.dofs();

    // node sets
    xt::xtensor<size_t, 1> nodesLft = mesh.nodesLeftOpenEdge();
    xt::xtensor<size_t, 1> nodesRgt = mesh.nodesRightOpenEdge();
    xt::xtensor<size_t, 1> nodesTop = mesh.nodesTopEdge();
    xt::xtensor<size_t, 1> nodesBot = mesh.nodesBottomEdge();

    // periodicity and fixed displacements DOFs
    // ----------------------------------------

    xt::view(dofs, xt::keep(nodesRgt)) = xt::view(dofs, xt::keep(nodesLft));

    dofs = GooseFEM::Mesh::renumber(dofs);

    xt::xtensor<size_t, 1> iip = xt::concatenate(xt::xtuple(
        xt::view(dofs, xt::keep(nodesBot), 0),
        xt::view(dofs, xt::keep(nodesBot), 1),
        xt::view(dofs, xt::keep(nodesTop), 0),
        xt::view(dofs, xt::keep(nodesTop), 1)));

    // simulation variables
    // --------------------

    // vector definition
    GooseFEM::VectorPartitioned vector(conn, dofs, iip);

    // allocate system matrix
    GooseFEM::MatrixPartitioned K(conn, dofs, iip);
    GooseFEM::MatrixPartitionedSolver<> Solver;

    // nodal quantities
    xt::xtensor<double, 2> disp = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> fint = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> fext = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> fres = xt::zeros<double>(coor.shape());

    // element vectors
    xt::xtensor<double, 3> ue = xt::empty<double>({nelem, nne, ndim});
    xt::xtensor<double, 3> fe = xt::empty<double>({nelem, nne, ndim});
    xt::xtensor<double, 3> Ke = xt::empty<double>({nelem, nne * ndim, nne * ndim});

    // element/material definition
    // ---------------------------

    // element definition
    GooseFEM::Element::Quad4::QuadraturePlanar elem(vector.AsElement(coor));
    size_t nip = elem.nip();

    // material definition
    GMatElastic::Cartesian3d::Matrix mat(nelem, nip);
    xt::xtensor<size_t, 2> Ihard = xt::zeros<size_t>({nelem, nip});
    xt::view(Ihard, xt::keep(0, 1, 5, 6), xt::all()) = 1;
    xt::xtensor<size_t, 2> Isoft = xt::ones<size_t>({nelem, nip}) - Ihard;
    mat.setElastic(Isoft, 10.0, 1.0);
    mat.setElastic(Ihard, 10.0, 10.0);

    // integration point tensors
    xt::xtensor<double, 4> Eps = xt::empty<double>({nelem, nip, 3ul, 3ul});
    xt::xtensor<double, 4> Sig = xt::empty<double>({nelem, nip, 3ul, 3ul});
    xt::xtensor<double, 6> C = xt::empty<double>({nelem, nip, 3ul, 3ul, 3ul, 3ul});

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
    xt::view(disp, xt::keep(nodesTop), 0) = +0.1;

    // residual
    xt::noalias(fres) = fext - fint;

    // solve
    Solver.solve(K, fres, disp);

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
    vector.copy_p(fint, fext);

    // residual
    xt::noalias(fres) = fext - fint;

    // print residual
    std::cout << xt::sum(xt::abs(fres))[0] / xt::sum(xt::abs(fext))[0] << std::endl;

    // average stress per node
    xt::xtensor<double, 4> dV = elem.AsTensor<2>(elem.dV());
    xt::xtensor<double, 3> SigAv = xt::average(Sig, dV, {1});

    // write output
    H5Easy::File file("output.h5", H5Easy::File::Overwrite);
    H5Easy::dump(file, "/coor", coor);
    H5Easy::dump(file, "/conn", conn);
    H5Easy::dump(file, "/disp", disp);
    H5Easy::dump(file, "/Sig", SigAv);

    return 0;
}

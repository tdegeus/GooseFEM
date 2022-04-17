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

    // mesh definition, displacement, external forces
    xt::xtensor<double, 2> coor = mesh.coor();
    xt::xtensor<size_t, 2> conn = mesh.conn();
    xt::xtensor<size_t, 2> dofs = mesh.dofs();
    auto disp = xt::zeros_like(coor);
    auto fext = xt::zeros_like(coor);

    // node sets
    xt::xtensor<size_t, 1> nodesLft = mesh.nodesLeftEdge();
    xt::xtensor<size_t, 1> nodesRgt = mesh.nodesRightEdge();
    xt::xtensor<size_t, 1> nodesTop = mesh.nodesTopEdge();
    xt::xtensor<size_t, 1> nodesBot = mesh.nodesBottomEdge();

    // fixed displacements DOFs
    // ------------------------

    xt::xtensor<size_t, 1> iip = xt::concatenate(xt::xtuple(
        xt::view(dofs, xt::keep(nodesRgt), 0),
        xt::view(dofs, xt::keep(nodesTop), 1),
        xt::view(dofs, xt::keep(nodesLft), 0),
        xt::view(dofs, xt::keep(nodesBot), 1)));

    // simulation variables
    // --------------------

    // vector definition
    GooseFEM::VectorPartitioned vector(conn, dofs, iip);

    // allocate system matrix
    GooseFEM::MatrixPartitioned K(conn, dofs, iip);
    GooseFEM::MatrixPartitionedSolver<> Solver;

    // element/material definition
    // ---------------------------

    // element definition
    GooseFEM::Element::Quad4::QuadraturePlanar elem(vector.AsElement(coor));
    size_t nip = elem.nip();

    // material definition
    GMatElastic::Cartesian3d::Array<2> mat({nelem, nip}, 1.0, 1.0);

    // solve
    // -----

    // strain
    auto Eps = elem.SymGradN_vector(vector.AsElement(disp));

    // stress & tangent
    mat.setStrain(Eps);
    auto Sig = mat.Stress();
    auto C = mat.Tangent();

    // internal force
    auto fint = vector.AssembleNode(elem.Int_gradN_dot_tensor2_dV(Sig));

    // stiffness matrix
    K.assemble(elem.Int_gradN_dot_tensor4_dot_gradNT_dV(C));

    // set fixed displacements
    xt::view(disp, xt::keep(nodesRgt), 0) = +0.1;
    xt::view(disp, xt::keep(nodesTop), 1) = -0.1;
    xt::view(disp, xt::keep(nodesLft), 0) = 0.0;
    xt::view(disp, xt::keep(nodesBot), 1) = 0.0;

    // residual
    auto fres = fext - fint;

    // solve
    Solver.solve(K, fres, disp);

    // post-process
    // ------------

    // compute strain and stress
    Eps = elem.SymGradN_vector(vector.AsElement(disp));
    mat.setStrain(Eps);
    Sig = mat.Stress();

    // internal force
    fint = vector.AssembleNode(elem.Int_gradN_dot_tensor2_dV(Sig));

    // apply reaction force
    vector.copy_p(fint, fext);

    // residual
    xt::noalias(fres) = fext - fint;

    // print residual
    std::cout << xt::sum(xt::abs(fres))[0] / xt::sum(xt::abs(fext))[0] << std::endl;

    // average stress per element
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
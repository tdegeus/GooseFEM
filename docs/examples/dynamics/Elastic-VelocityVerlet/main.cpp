#include <GMatElastoPlasticQPot/Cartesian2d.h>
#include <GooseFEM/GooseFEM.h>
#include <highfive/H5Easy.hpp>

namespace GM = GMatElastoPlasticQPot::Cartesian2d;
namespace GF = GooseFEM;
namespace QD = GooseFEM::Element::Quad4;
namespace H5 = H5Easy;

inline double sqdot(const xt::xtensor<double, 1>& M, const xt::xtensor<double, 1>& V)
{
    double ret = 0.;

    for (size_t i = 0; i < M.size(); ++i) {
        ret += M(i) * V(i) * V(i);
    }

    return ret;
}

int main()
{
    // simulation parameters

    double T = 60.0;     // total time
    double dt = 1.0e-2;  // time increment
    size_t nx = 60;      // number of elements in both directions
    double gamma = 0.05; // displacement step

    // get mesh & quadrature

    GF::Mesh::Quad4::Regular mesh(nx, nx, 1.0);

    xt::xtensor<double, 2> coor = mesh.coor();
    xt::xtensor<size_t, 2> conn = mesh.conn();
    xt::xtensor<size_t, 2> dofs = mesh.dofsPeriodic();

    GF::Vector vector(conn, dofs);

    QD::Quadrature quad(vector.AsElement(coor));

    size_t ndim = mesh.ndim();
    size_t nne = mesh.nne();
    size_t nelem = mesh.nelem();
    size_t nip = quad.nip();

    xt::xtensor<double, 2> u = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> v = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> a = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> v_n = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> a_n = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> fint = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> fext = xt::zeros<double>(coor.shape());
    xt::xtensor<double, 2> fres = xt::zeros<double>(coor.shape());

    xt::xtensor<double, 3> ue = xt::zeros<double>({nelem, nne, ndim});
    xt::xtensor<double, 3> fe = xt::zeros<double>({nelem, nne, ndim});

    xt::xtensor<double, 4> Eps = xt::zeros<double>({nelem, nip, ndim, ndim});
    xt::xtensor<double, 4> Sig = xt::zeros<double>({nelem, nip, ndim, ndim});

    // material definition

    GM::Matrix material({nelem, nip});

    xt::xtensor<size_t, 2> Ihard = xt::zeros<size_t>({nelem, nip});
    xt::xtensor<size_t, 2> Isoft = xt::ones<size_t>({nelem, nip});
    xt::view(Ihard, xt::range(0, nx * nx / 4), xt::all()) = 1ul;
    Isoft -= Ihard;

    material.setElastic(Ihard, 100.0, 10.0);
    material.setElastic(Isoft, 100.0, 1.0);

    material.check();

    // mass matrix

    QD::Quadrature nodalQuad(vector.AsElement(coor), QD::Nodal::xi(), QD::Nodal::w());

    xt::xtensor<double, 2> rho = 1.0 * xt::ones<double>({nelem, nip});

    GF::MatrixDiagonal M(conn, dofs);

    M.assemble(nodalQuad.Int_N_scalar_NT_dV(rho));

    xt::xtensor<double, 1> mass = M.AsDiagonal();

    // update in macroscopic deformation gradient

    xt::xtensor<double, 2> dFbar = xt::zeros<double>({2, 2});

    dFbar(0, 1) = gamma;

    for (size_t j = 0; j < ndim; ++j) {
        for (size_t k = 0; k < ndim; ++k) {
            xt::view(u, xt::all(), j) += dFbar(j, k) * (xt::view(coor, xt::all(), k) - coor(0, k));
        }
    }

    // output variables

    xt::xtensor<double, 1> Epot = xt::zeros<double>({static_cast<size_t>(T / dt)});
    xt::xtensor<double, 1> Ekin = xt::zeros<double>({static_cast<size_t>(T / dt)});
    xt::xtensor<double, 1> t = xt::zeros<double>({static_cast<size_t>(T / dt)});

    xt::xtensor<double, 2> dV = quad.DV();

    // loop over increments

    for (size_t inc = 0; inc < static_cast<size_t>(Epot.size()); ++inc) {

        // store history

        xt::noalias(v_n) = v;
        xt::noalias(a_n) = a;

        // new displacement

        xt::noalias(u) = u + dt * v + 0.5 * std::pow(dt, 2.0) * a;

        // compute strain/strain, and corresponding internal force

        vector.asElement(u, ue);
        quad.symGradN_vector(ue, Eps);
        material.stress(Eps, Sig);
        quad.int_gradN_dot_tensor2_dV(Sig, fe);
        vector.assembleNode(fe, fint);

        // estimate new velocity

        xt::noalias(v) = v_n + dt * a_n;

        // compute residual force & solve

        xt::noalias(fres) = fext - fint;

        M.solve(fres, a);

        // re-estimate new velocity

        xt::noalias(v) = v_n + 0.5 * dt * (a_n + a);

        // compute residual force & solve

        xt::noalias(fres) = fext - fint;

        M.solve(fres, a);

        // new velocity

        xt::noalias(v) = v_n + 0.5 * dt * (a_n + a);

        // compute residual force & solve

        xt::noalias(fres) = fext - fint;

        M.solve(fres, a);

        // store output variables

        xt::xtensor<double, 2> E = material.Energy(Eps);
        xt::xtensor<double, 1> V = vector.AsDofs(v);

        t(inc) = static_cast<double>(inc) * dt;
        Ekin(inc) = 0.5 * sqdot(mass, V);
        Epot(inc) = xt::sum(E * dV)[0];
    }

    // write output variables to file
    H5::File file("example.hdf5", H5::File::Overwrite);
    H5::dump(file, "/global/Epot", Epot);
    H5::dump(file, "/global/Ekin", Ekin);
    H5::dump(file, "/global/t", t);
    H5::dump(file, "/mesh/conn", conn);
    H5::dump(file, "/mesh/coor", coor);
    H5::dump(file, "/mesh/disp", u);

    return 0;
}

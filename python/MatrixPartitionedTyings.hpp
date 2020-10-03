/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pyxtensor/pyxtensor.hpp>
#include "../include/GooseFEM/GooseFEM.h"

// =================================================================================================

namespace py = pybind11;

// =================================================================================================

void init_MatrixPartitionedTyings(py::module& m)
{

py::class_<GooseFEM::MatrixPartitionedTyings>(m, "MatrixPartitionedTyings")

    .def(py::init<
            const xt::xtensor<size_t,2>&,
            const xt::xtensor<size_t,2>&,
            const Eigen::SparseMatrix<double>&,
            const Eigen::SparseMatrix<double>&>(),
        "Sparse, partitioned, matrix",
        py::arg("conn"),
        py::arg("dofs"),
        py::arg("Cdu"),
        py::arg("Cdp"))

    .def("nelem",
        &GooseFEM::MatrixPartitionedTyings::nelem,
        "Number of element")

    .def("nne",
        &GooseFEM::MatrixPartitionedTyings::nne,
        "Number of nodes per element")

    .def("nnode",
        &GooseFEM::MatrixPartitionedTyings::nnode,
        "Number of nodes")

    .def("ndim",
        &GooseFEM::MatrixPartitionedTyings::ndim,
        "Number of dimensions")

    .def("ndof",
        &GooseFEM::MatrixPartitionedTyings::ndof,
        "Number of degrees-of-freedom")

    .def("nnu",
        &GooseFEM::MatrixPartitionedTyings::nnu,
        "Number of unknown degrees-of-freedom")

    .def("nnp",
        &GooseFEM::MatrixPartitionedTyings::nnp,
        "Number of prescribed degrees-of-freedom")

    .def("nni",
        &GooseFEM::MatrixPartitionedTyings::nni,
        "Number of independent degrees-of-freedom")

    .def("nnd",
        &GooseFEM::MatrixPartitionedTyings::nnd,
        "Number of dependent degrees-of-freedom")

    .def("assemble",
        &GooseFEM::MatrixPartitionedTyings::assemble,
        "Assemble matrix from 'elemmat",
        py::arg("elemmat"))

    .def("dofs",
        &GooseFEM::MatrixPartitionedTyings::dofs,
        "Return degrees-of-freedom")

    .def("iiu",
        &GooseFEM::MatrixPartitionedTyings::iiu,
        "Return unknown degrees-of-freedom")

    .def("iip",
        &GooseFEM::MatrixPartitionedTyings::iip,
        "Return prescribed degrees-of-freedom")

    .def("iii",
        &GooseFEM::MatrixPartitionedTyings::iii,
        "Return independent degrees-of-freedom")

    .def("iid",
        &GooseFEM::MatrixPartitionedTyings::iid,
        "Return dependent degrees-of-freedom")

    .def("__repr__",
        [](const GooseFEM::MatrixPartitionedTyings&){
            return "<GooseFEM.MatrixPartitionedTyings>"; });


py::class_<GooseFEM::MatrixPartitionedTyingsSolver<>>(m, "MatrixPartitionedTyingsSolver")

    .def(py::init<>(),
        "Sparse, partitioned, matrix solver")

    .def("Solve",
        py::overload_cast<
            GooseFEM::MatrixPartitionedTyings&,
            const xt::xtensor<double,1>&,
            const xt::xtensor<double,1>&>(
                &GooseFEM::MatrixPartitionedTyingsSolver<>::Solve),
        "Solve",
        py::arg("matrix"),
        py::arg("b"),
        py::arg("x"))

    .def("Solve",
        py::overload_cast<
            GooseFEM::MatrixPartitionedTyings&,
            const xt::xtensor<double,2>&,
            const xt::xtensor<double,2>&>(
                &GooseFEM::MatrixPartitionedTyingsSolver<>::Solve),
        "Solve",
        py::arg("matrix"),
        py::arg("b"),
        py::arg("x"))

    .def("Solve_u",
        py::overload_cast<
            GooseFEM::MatrixPartitionedTyings&,
            const xt::xtensor<double,1>&,
            const xt::xtensor<double,1>&,
            const xt::xtensor<double,1>&>(
                &GooseFEM::MatrixPartitionedTyingsSolver<>::Solve_u),
        "Solve_u",
        py::arg("matrix"),
        py::arg("b_u"),
        py::arg("b_d"),
        py::arg("x_p"))

    .def("__repr__",
        [](const GooseFEM::MatrixPartitionedTyingsSolver<>&){
            return "<GooseFEM.MatrixPartitionedTyingsSolver>"; });

}

// =================================================================================================


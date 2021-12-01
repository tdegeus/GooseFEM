/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_MatrixPartitionedTyings(py::module& m)
{

    py::class_<GooseFEM::MatrixPartitionedTyings>(m, "MatrixPartitionedTyings")

        .def(
            py::init<
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 2>&,
                const Eigen::SparseMatrix<double>&,
                const Eigen::SparseMatrix<double>&>(),
            "See :cpp:class:`GooseFEM::MatrixPartitionedTyings`.",
            py::arg("conn"),
            py::arg("dofs"),
            py::arg("Cdu"),
            py::arg("Cdp"))

        .def("nelem", &GooseFEM::MatrixPartitionedTyings::nelem)
        .def("nne", &GooseFEM::MatrixPartitionedTyings::nne)
        .def("nnode", &GooseFEM::MatrixPartitionedTyings::nnode)
        .def("ndim", &GooseFEM::MatrixPartitionedTyings::ndim)
        .def("ndof", &GooseFEM::MatrixPartitionedTyings::ndof)
        .def("nnu", &GooseFEM::MatrixPartitionedTyings::nnu)
        .def("nnp", &GooseFEM::MatrixPartitionedTyings::nnp)
        .def("nni", &GooseFEM::MatrixPartitionedTyings::nni)
        .def("nnd", &GooseFEM::MatrixPartitionedTyings::nnd)
        .def("assemble", &GooseFEM::MatrixPartitionedTyings::assemble, py::arg("elemmat"))
        .def("dofs", &GooseFEM::MatrixPartitionedTyings::dofs)
        .def("iiu", &GooseFEM::MatrixPartitionedTyings::iiu)
        .def("iip", &GooseFEM::MatrixPartitionedTyings::iip)
        .def("iii", &GooseFEM::MatrixPartitionedTyings::iii)
        .def("iid", &GooseFEM::MatrixPartitionedTyings::iid)

        .def("__repr__", [](const GooseFEM::MatrixPartitionedTyings&) {
            return "<GooseFEM.MatrixPartitionedTyings>";
        });

    py::class_<GooseFEM::MatrixPartitionedTyingsSolver<>>(m, "MatrixPartitionedTyingsSolver")

        .def(py::init<>(), "See :cpp:class:`GooseFEM::MatrixPartitionedTyingsSolver`.")

        .def(
            "Solve",
            py::overload_cast<
                GooseFEM::MatrixPartitionedTyings&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&>(&GooseFEM::MatrixPartitionedTyingsSolver<>::Solve),
            py::arg("matrix"),
            py::arg("b"),
            py::arg("x"))

        .def(
            "Solve",
            py::overload_cast<
                GooseFEM::MatrixPartitionedTyings&,
                const xt::xtensor<double, 2>&,
                const xt::xtensor<double, 2>&>(&GooseFEM::MatrixPartitionedTyingsSolver<>::Solve),
            py::arg("matrix"),
            py::arg("b"),
            py::arg("x"))

        .def(
            "Solve_u",
            py::overload_cast<
                GooseFEM::MatrixPartitionedTyings&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&>(&GooseFEM::MatrixPartitionedTyingsSolver<>::Solve_u),
            py::arg("matrix"),
            py::arg("b_u"),
            py::arg("b_d"),
            py::arg("x_p"))

        .def("__repr__", [](const GooseFEM::MatrixPartitionedTyingsSolver<>&) {
            return "<GooseFEM.MatrixPartitionedTyingsSolver>";
        });
}

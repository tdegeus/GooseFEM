/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#include <Eigen/Eigen>
#include <GooseFEM/Matrix.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_Matrix(py::module& m)
{

    py::class_<GooseFEM::Matrix>(m, "Matrix")

        .def(
            py::init<const xt::xtensor<size_t, 2>&, const xt::xtensor<size_t, 2>&>(),
            "See :cpp:class:`GooseFEM::Matrix`.",
            py::arg("conn"),
            py::arg("dofs"))

        .def("nelem", &GooseFEM::Matrix::nelem)
        .def("nne", &GooseFEM::Matrix::nne)
        .def("nnode", &GooseFEM::Matrix::nnode)
        .def("ndim", &GooseFEM::Matrix::ndim)
        .def("ndof", &GooseFEM::Matrix::ndof)
        .def("dofs", &GooseFEM::Matrix::dofs)
        .def("assemble", &GooseFEM::Matrix::assemble, py::arg("elemmat"))
        .def("set", &GooseFEM::Matrix::set, py::arg("rows"), py::arg("cols"), py::arg("matrix"))
        .def("add", &GooseFEM::Matrix::add, py::arg("rows"), py::arg("cols"), py::arg("matrix"))
        .def("Todense", &GooseFEM::Matrix::Todense)

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 1>&>(&GooseFEM::Matrix::Dot, py::const_),
            py::arg("x"))

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 2>&>(&GooseFEM::Matrix::Dot, py::const_),
            py::arg("x"))

        .def("__repr__", [](const GooseFEM::Matrix&) { return "<GooseFEM.Matrix>"; });

    py::class_<GooseFEM::MatrixSolver<>>(m, "MatrixSolver")

        .def(py::init<>(), "See :cpp:class:`GooseFEM::MatrixSolver`.")

        .def(
            "Solve",
            py::overload_cast<GooseFEM::Matrix&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixSolver<>::Solve),
            py::arg("A"),
            py::arg("b"))

        .def(
            "Solve",
            py::overload_cast<GooseFEM::Matrix&, const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixSolver<>::Solve),
            py::arg("A"),
            py::arg("b"))

        .def("__repr__", [](const GooseFEM::MatrixSolver<>&) { return "<GooseFEM.MatrixSolver>"; });
}

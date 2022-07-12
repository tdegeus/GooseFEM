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
            py::init<const xt::pytensor<size_t, 2>&, const xt::pytensor<size_t, 2>&>(),
            "See :cpp:class:`GooseFEM::Matrix`.",
            py::arg("conn"),
            py::arg("dofs"))

        .def_property_readonly("nelem", &GooseFEM::Matrix::nelem)
        .def_property_readonly("nne", &GooseFEM::Matrix::nne)
        .def_property_readonly("nnode", &GooseFEM::Matrix::nnode)
        .def_property_readonly("ndim", &GooseFEM::Matrix::ndim)
        .def_property_readonly("ndof", &GooseFEM::Matrix::ndof)
        .def_property_readonly("dofs", &GooseFEM::Matrix::dofs)
        .def("assemble", &GooseFEM::Matrix::assemble, py::arg("elemmat"))
        .def("set", &GooseFEM::Matrix::set, py::arg("rows"), py::arg("cols"), py::arg("matrix"))
        .def("add", &GooseFEM::Matrix::add, py::arg("rows"), py::arg("cols"), py::arg("matrix"))
        .def("Todense", &GooseFEM::Matrix::Todense)

        .def(
            "Dot",
            py::overload_cast<const xt::pytensor<double, 1>&>(&GooseFEM::Matrix::Dot, py::const_),
            py::arg("x"))

        .def(
            "Dot",
            py::overload_cast<const xt::pytensor<double, 2>&>(&GooseFEM::Matrix::Dot, py::const_),
            py::arg("x"))

        .def("__repr__", [](const GooseFEM::Matrix&) { return "<GooseFEM.Matrix>"; });

    py::class_<GooseFEM::MatrixSolver<>>(m, "MatrixSolver")

        .def(py::init<>(), "See :cpp:class:`GooseFEM::MatrixSolver`.")

        .def(
            "Solve",
            py::overload_cast<GooseFEM::Matrix&, const xt::pytensor<double, 1>&>(
                &GooseFEM::MatrixSolver<>::Solve),
            py::arg("A"),
            py::arg("b"))

        .def(
            "Solve",
            py::overload_cast<GooseFEM::Matrix&, const xt::pytensor<double, 2>&>(
                &GooseFEM::MatrixSolver<>::Solve),
            py::arg("A"),
            py::arg("b"))

        .def("__repr__", [](const GooseFEM::MatrixSolver<>&) { return "<GooseFEM.MatrixSolver>"; });
}

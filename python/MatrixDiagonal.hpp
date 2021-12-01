/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_MatrixDiagonal(py::module& m)
{

    py::class_<GooseFEM::MatrixDiagonal>(m, "MatrixDiagonal")

        .def(
            py::init<const xt::xtensor<size_t, 2>&, const xt::xtensor<size_t, 2>&>(),
            "See :cpp:class:`GooseFEM::MatrixDiagonal`.",
            py::arg("conn"),
            py::arg("dofs"))

        .def("nelem", &GooseFEM::MatrixDiagonal::nelem)
        .def("nne", &GooseFEM::MatrixDiagonal::nne)
        .def("nnode", &GooseFEM::MatrixDiagonal::nnode)
        .def("ndim", &GooseFEM::MatrixDiagonal::ndim)
        .def("ndof", &GooseFEM::MatrixDiagonal::ndof)
        .def("dofs", &GooseFEM::MatrixDiagonal::dofs)
        .def("assemble", &GooseFEM::MatrixDiagonal::assemble, py::arg("elemmat"))
        .def("set", &GooseFEM::MatrixDiagonal::set, py::arg("A"))
        .def("Todiagonal", &GooseFEM::MatrixDiagonal::Todiagonal)

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonal::Dot, py::const_),
            py::arg("x"))

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixDiagonal::Dot, py::const_),
            py::arg("x"))

        .def(
            "Solve",
            py::overload_cast<const xt::xtensor<double, 1>&>(&GooseFEM::MatrixDiagonal::Solve),
            py::arg("b"))

        .def(
            "Solve",
            py::overload_cast<const xt::xtensor<double, 2>&>(&GooseFEM::MatrixDiagonal::Solve),
            py::arg("b"))

        .def("__repr__", [](const GooseFEM::MatrixDiagonal&) {
            return "<GooseFEM.MatrixDiagonal>";
        });
}

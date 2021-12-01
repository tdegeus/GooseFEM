/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#include <GooseFEM/MatrixDiagonalPartitioned.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_MatrixDiagonalPartitioned(py::module& m)
{

    py::class_<GooseFEM::MatrixDiagonalPartitioned, GooseFEM::MatrixDiagonal>(
        m, "MatrixDiagonalPartitioned")

        .def(
            py::init<
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 1>&>(),
            "See :cpp:class:`GooseFEM::MatrixDiagonalPartitioned`.",
            py::arg("conn"),
            py::arg("dofs"),
            py::arg("iip"))

        .def("nnu", &GooseFEM::MatrixDiagonalPartitioned::nnu)
        .def("nnp", &GooseFEM::MatrixDiagonalPartitioned::nnp)
        .def("iiu", &GooseFEM::MatrixDiagonalPartitioned::iiu)
        .def("iip", &GooseFEM::MatrixDiagonalPartitioned::iip)

        .def(
            "Dot_u",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Dot_u, py::const_),
            py::arg("x_u"),
            py::arg("x_p"))

        .def(
            "Dot_p",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Dot_p, py::const_),
            py::arg("x_u"),
            py::arg("x_p"))

        .def(
            "Solve_u",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Solve_u),
            py::arg("b_u"),
            py::arg("x_p"))

        .def(
            "Reaction",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Reaction, py::const_),
            py::arg("x"),
            py::arg("b"))

        .def(
            "Reaction",
            py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Reaction, py::const_),
            py::arg("x"),
            py::arg("b"))

        .def(
            "Reaction_p",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Reaction_p, py::const_),
            py::arg("x_u"),
            py::arg("x_p"))

        .def("__repr__", [](const GooseFEM::MatrixDiagonalPartitioned&) {
            return "<GooseFEM.MatrixDiagonalPartitioned>";
        });
}

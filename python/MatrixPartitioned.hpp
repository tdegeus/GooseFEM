/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#include <Eigen/Eigen>
#include <GooseFEM/MatrixPartitioned.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include "Matrix.hpp"

namespace py = pybind11;

void init_MatrixPartitioned(py::module& m)
{
    py::class_<GooseFEM::MatrixPartitioned> cls(m, "MatrixPartitioned");
    register_Matrix_MatrixBase<GooseFEM::MatrixPartitioned>(cls);
    register_Matrix_MatrixPartitionedBase<GooseFEM::MatrixPartitioned>(cls);

    cls.def(
        py::init<
            const xt::pytensor<size_t, 2>&,
            const xt::pytensor<size_t, 2>&,
            const xt::pytensor<size_t, 1>&>(),
        "See :cpp:class:`GooseFEM::MatrixPartitioned`.",
        py::arg("conn"),
        py::arg("dofs"),
        py::arg("iip"));

    cls.def("__repr__", [](const GooseFEM::MatrixPartitioned&) {
        return "<GooseFEM.MatrixPartitioned>";
    });

    py::class_<GooseFEM::MatrixPartitionedSolver<>>(m, "MatrixPartitionedSolver")

        .def(py::init<>(), "See :cpp:class:`GooseFEM::MatrixPartitionedSolver`.")

        .def(
            "Solve",
            py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::pytensor<double, 1>&,
                const xt::pytensor<double, 1>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve),
            py::arg("matrix"),
            py::arg("b"),
            py::arg("x"))

        .def(
            "Solve",
            py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::pytensor<double, 2>&,
                const xt::pytensor<double, 2>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve),
            py::arg("matrix"),
            py::arg("b"),
            py::arg("x"))

        .def(
            "Solve_u",
            py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::pytensor<double, 1>&,
                const xt::pytensor<double, 1>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve_u),
            py::arg("matrix"),
            py::arg("b_u"),
            py::arg("x_p"))

        .def("__repr__", [](const GooseFEM::MatrixPartitionedSolver<>&) {
            return "<GooseFEM.MatrixPartitionedSolver>";
        });
}

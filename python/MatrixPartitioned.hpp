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

    py::class_<GooseFEM::MatrixPartitionedSolver<>> slv(m, "MatrixPartitionedSolver");

    register_MatrixSolver_MatrixSolverBase<
        GooseFEM::MatrixPartitionedSolver<>,
        GooseFEM::MatrixPartitioned>(slv);

    register_MatrixSolver_MatrixSolverPartitionedBase<
        GooseFEM::MatrixPartitionedSolver<>,
        GooseFEM::MatrixPartitioned>(slv);

    slv.def(py::init<>(), "See :cpp:class:`GooseFEM::MatrixPartitionedSolver`.");

    slv.def("__repr__", [](const GooseFEM::MatrixPartitionedSolver<>&) {
        return "<GooseFEM.MatrixPartitionedSolver>";
    });
}

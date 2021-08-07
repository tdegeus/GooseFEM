/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <GooseFEM/MatrixPartitioned.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_MatrixPartitioned(py::module& m)
{

    py::class_<GooseFEM::MatrixPartitioned, GooseFEM::Matrix>(m, "MatrixPartitioned")

        .def(py::init<
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 1>&>(),
             "Sparse, partitioned, matrix",
             py::arg("conn"),
             py::arg("dofs"),
             py::arg("iip"))

        .def("nnu",
             &GooseFEM::MatrixPartitioned::nnu,
             "See :cpp:func:`GooseFEM::MatrixPartitioned::nnu`.")

        .def("nnp",
             &GooseFEM::MatrixPartitioned::nnp,
             "See :cpp:func:`GooseFEM::MatrixPartitioned::nnp`.")

        .def("iiu",
             &GooseFEM::MatrixPartitioned::iiu,
             "See :cpp:func:`GooseFEM::MatrixPartitioned::iiu`.")

        .def("iip",
             &GooseFEM::MatrixPartitioned::iip,
             "See :cpp:func:`GooseFEM::MatrixPartitioned::iip`.")

        .def("Reaction",
             py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixPartitioned::Reaction, py::const_),
             "See :cpp:func:`GooseFEM::MatrixPartitioned::Reaction`.",
             py::arg("x"),
             py::arg("b"))

        .def("Reaction",
             py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixPartitioned::Reaction, py::const_),
             "See :cpp:func:`GooseFEM::MatrixPartitioned::Reaction`.",
             py::arg("x"),
             py::arg("b"))

        .def("Reaction_p",
             py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixPartitioned::Reaction_p, py::const_),
             "See :cpp:func:`GooseFEM::MatrixPartitioned::Reaction_p`.",
             py::arg("x_u"),
             py::arg("x_p"))

        .def("__repr__", [](const GooseFEM::MatrixPartitioned&) {
            return "<GooseFEM.MatrixPartitioned>";
        });

    py::class_<GooseFEM::MatrixPartitionedSolver<>>(m, "MatrixPartitionedSolver")

        .def(py::init<>(), "Sparse, partitioned, matrix solver")

        .def("Solve",
             py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve),
             "See :cpp:func:`GooseFEM::MatrixPartitionedSolver::Solve`.",
             py::arg("matrix"),
             py::arg("b"),
             py::arg("x"))

        .def("Solve",
             py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::xtensor<double, 2>&,
                const xt::xtensor<double, 2>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve),
             "See :cpp:func:`GooseFEM::MatrixPartitionedSolver::Solve`.",
             py::arg("matrix"),
             py::arg("b"),
             py::arg("x"))

        .def("Solve_u",
             py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve_u),
             "See :cpp:func:`GooseFEM::MatrixPartitionedSolver::Solve_u`.",
             py::arg("matrix"),
             py::arg("b_u"),
             py::arg("x_p"))

        .def("__repr__", [](const GooseFEM::MatrixPartitionedSolver<>&) {
            return "<GooseFEM.MatrixPartitionedSolver>";
        });
}

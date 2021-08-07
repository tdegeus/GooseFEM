/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

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
            "Sparse matrix",
            py::arg("conn"),
            py::arg("dofs"))

        .def("nelem",
             &GooseFEM::Matrix::nelem,
             "See :cpp:func:`GooseFEM::Matrix::nelem`.")

        .def("nne",
             &GooseFEM::Matrix::nne,
             "See :cpp:func:`GooseFEM::Matrix::nne`.")

        .def("nnode",
             &GooseFEM::Matrix::nnode,
             "See :cpp:func:`GooseFEM::Matrix::nnode`.")

        .def("ndim",
             &GooseFEM::Matrix::ndim,
             "See :cpp:func:`GooseFEM::Matrix::ndim`.")

        .def("ndof",
             &GooseFEM::Matrix::ndof,
             "See :cpp:func:`GooseFEM::Matrix::ndof`.")

        .def("dofs",
             &GooseFEM::Matrix::dofs,
             "See :cpp:func:`GooseFEM::Matrix::dofs`.")

        .def("assemble",
             &GooseFEM::Matrix::assemble,
             "See :cpp:func:`GooseFEM::Matrix::assemble`.",
             py::arg("elemmat"))

        .def("set",
             &GooseFEM::Matrix::set,
             "See :cpp:func:`GooseFEM::Matrix::set`.",
             py::arg("rows"),
             py::arg("cols"),
             py::arg("matrix"))

        .def("add",
             &GooseFEM::Matrix::add,
             "See :cpp:func:`GooseFEM::Matrix::add`.",
             py::arg("rows"),
             py::arg("cols"),
             py::arg("matrix"))

        .def("Todense",
             &GooseFEM::Matrix::Todense,
             "See :cpp:func:`GooseFEM::Matrix::Todense`.")

        .def("Dot",
             py::overload_cast<const xt::xtensor<double, 1>&>(&GooseFEM::Matrix::Dot, py::const_),
             "See :cpp:func:`GooseFEM::Matrix::Dot`.",
             py::arg("x"))

        .def("Dot",
             py::overload_cast<const xt::xtensor<double, 2>&>(&GooseFEM::Matrix::Dot, py::const_),
             "See :cpp:func:`GooseFEM::Matrix::Dot`.",
             py::arg("x"))

        .def("__repr__", [](const GooseFEM::Matrix&) { return "<GooseFEM.Matrix>"; });

    py::class_<GooseFEM::MatrixSolver<>>(m, "MatrixSolver")

        .def(py::init<>(), "Sparse matrix solver")

        .def("Solve",
             py::overload_cast<GooseFEM::Matrix&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixSolver<>::Solve),
             "See :cpp:func:`GooseFEM::MatrixSolver::Solve`.",
             py::arg("A"),
             py::arg("b"))

        .def("Solve",
             py::overload_cast<GooseFEM::Matrix&, const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixSolver<>::Solve),
             "See :cpp:func:`GooseFEM::MatrixSolver::Solve`.",
             py::arg("A"),
             py::arg("b"))

        .def("__repr__", [](const GooseFEM::MatrixSolver<>&) { return "<GooseFEM.MatrixSolver>"; });
}

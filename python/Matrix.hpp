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

void init_Matrix(py::module& m)
{

py::class_<GooseFEM::Matrix<>>(m, "Matrix")

    .def(py::init<
            const xt::xtensor<size_t,2>&,
            const xt::xtensor<size_t,2>&>(),
        "Sparse matrix",
        py::arg("conn"),
        py::arg("dofs"))

    .def("nelem",
        &GooseFEM::Matrix<>::nelem,
        "Number of element")

    .def("nne",
        &GooseFEM::Matrix<>::nne,
        "Number of nodes per element")

    .def("nnode",
        &GooseFEM::Matrix<>::nnode,
        "Number of nodes")

    .def("ndim",
        &GooseFEM::Matrix<>::ndim,
        "Number of dimensions")

    .def("ndof",
        &GooseFEM::Matrix<>::ndof,
        "Number of degrees-of-freedom")

    .def("assemble",
        &GooseFEM::Matrix<>::assemble,
        "Assemble matrix from 'elemmat",
        py::arg("elemmat"))

    .def("dofs",
        &GooseFEM::Matrix<>::dofs,
        "Return degrees-of-freedom")

    .def("Dot",
        py::overload_cast<const xt::xtensor<double,1>&>(&GooseFEM::Matrix<>::Dot, py::const_),
        "Dot",
        py::arg("x"))

    .def("Dot",
        py::overload_cast<const xt::xtensor<double,2>&>(&GooseFEM::Matrix<>::Dot, py::const_),
        "Dot",
        py::arg("x"))

    .def("Solve",
        py::overload_cast<const xt::xtensor<double,1>&>(&GooseFEM::Matrix<>::Solve),
        "Solve",
        py::arg("b"))

    .def("Solve",
        py::overload_cast<const xt::xtensor<double,2>&>(&GooseFEM::Matrix<>::Solve),
        "Solve",
        py::arg("b"))

    .def("__repr__",
        [](const GooseFEM::Matrix<>&){
            return "<GooseFEM.Matrix>"; });

}

// =================================================================================================


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

void init_MatrixPartitioned(py::module &m)
{

py::class_<GooseFEM::MatrixPartitioned<>>(m, "MatrixPartitioned")

  .def(py::init<const xt::xtensor<size_t,2>&,
                const xt::xtensor<size_t,2>&,
                const xt::xtensor<size_t,1>&>(),
    "Sparse, partitioned, matrix",
    py::arg("conn"),
    py::arg("dofs"),
    py::arg("iip"))

  .def("nelem",
    &GooseFEM::MatrixPartitioned<>::nelem,
    "Number of element")

  .def("nne",
    &GooseFEM::MatrixPartitioned<>::nne,
    "Number of nodes per element")

  .def("nnode",
    &GooseFEM::MatrixPartitioned<>::nnode,
    "Number of nodes")

  .def("ndim",
    &GooseFEM::MatrixPartitioned<>::ndim,
    "Number of dimensions")

  .def("ndof",
    &GooseFEM::MatrixPartitioned<>::ndof,
    "Number of degrees-of-freedom")

  .def("nnu",
    &GooseFEM::MatrixPartitioned<>::nnu,
    "Number of unknown degrees-of-freedom")

  .def("nnp",
    &GooseFEM::MatrixPartitioned<>::nnp,
    "Number of prescribed degrees-of-freedom")

  .def("assemble",
    &GooseFEM::MatrixPartitioned<>::assemble,
    "Assemble matrix from 'elemmat",
    py::arg("elemmat"))

  .def("dofs",
    &GooseFEM::MatrixPartitioned<>::dofs,
    "Return degrees-of-freedom")

  .def("iiu",
    &GooseFEM::MatrixPartitioned<>::iiu,
    "Return unknown degrees-of-freedom")

  .def("iip",
    &GooseFEM::MatrixPartitioned<>::iip,
    "Return prescribed degrees-of-freedom")

  .def("Solve",
    py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(
      &GooseFEM::MatrixPartitioned<>::Solve),
    "Solve",
    py::arg("b"),
    py::arg("x"))

  .def("Solve",
    py::overload_cast<const xt::xtensor<double,2>&, const xt::xtensor<double,2>&>(
      &GooseFEM::MatrixPartitioned<>::Solve),
    "Solve",
    py::arg("b"),
    py::arg("x"))

  .def("Solve_u",
    py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(
      &GooseFEM::MatrixPartitioned<>::Solve_u),
    "Solve_u",
    py::arg("b_u"),
    py::arg("x_p"))

  .def("Reaction",
    py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(
      &GooseFEM::MatrixPartitioned<>::Reaction, py::const_),
    "Reaction",
    py::arg("x"),
    py::arg("b"))

  .def("Reaction",
    py::overload_cast<const xt::xtensor<double,2>&, const xt::xtensor<double,2>&>(
      &GooseFEM::MatrixPartitioned<>::Reaction, py::const_),
    "Reaction",
    py::arg("x"),
    py::arg("b"))

  .def("Reaction_p",
    py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(
      &GooseFEM::MatrixPartitioned<>::Reaction_p, py::const_),
    "Reaction_p",
    py::arg("x_u"),
    py::arg("x_p"))

  .def("__repr__", [](const GooseFEM::MatrixPartitioned<>&){
    return "<GooseFEM.MatrixPartitioned>";});

}

// =================================================================================================


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

  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(), " matrix", py::arg("conn"), py::arg("dofs"), py::arg("iip"))

  .def("nelem", &GooseFEM::MatrixPartitioned<>::nelem, "Return number of element")
  .def("nne"  , &GooseFEM::MatrixPartitioned<>::nne  , "Return number of nodes per element")
  .def("nnode", &GooseFEM::MatrixPartitioned<>::nnode, "Return number of nodes")
  .def("ndim" , &GooseFEM::MatrixPartitioned<>::ndim , "Return number of dimensions")
  .def("ndof" , &GooseFEM::MatrixPartitioned<>::ndof , "Return number of degrees-of-freedom")
  .def("nnu"  , &GooseFEM::MatrixPartitioned<>::nnu  , "Return number of unknown degrees-of-freedom")
  .def("nnp"  , &GooseFEM::MatrixPartitioned<>::nnp  , "Return number of prescribed degrees-of-freedom")

  .def("assemble", &GooseFEM::MatrixPartitioned<>::assemble, "Assemble matrix from 'elemmat", py::arg("elemmat"))

  .def("dofs" , &GooseFEM::MatrixPartitioned<>::dofs , "Return degrees-of-freedom")
  .def("iiu"  , &GooseFEM::MatrixPartitioned<>::iiu  , "Return unknown degrees-of-freedom")
  .def("iip"  , &GooseFEM::MatrixPartitioned<>::iip  , "Return prescribed degrees-of-freedom")

  .def("Solve", py::overload_cast<const xt::xtensor<double,1> &, const xt::xtensor<double,1> &>(&GooseFEM::MatrixPartitioned<>::Solve), "Solve")
  .def("Solve", py::overload_cast<const xt::xtensor<double,2> &, const xt::xtensor<double,2> &>(&GooseFEM::MatrixPartitioned<>::Solve), "Solve")

  .def("__repr__", [](const GooseFEM::MatrixPartitioned<> &){ return "<GooseFEM.MatrixPartitioned>"; });

}

// =================================================================================================


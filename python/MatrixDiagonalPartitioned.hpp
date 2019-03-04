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

void init_MatrixDiagonalPartitioned(py::module &m)
{

py::class_<GooseFEM::MatrixDiagonalPartitioned>(m, "MatrixDiagonalPartitioned")

  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(), "Diagonal matrix", py::arg("conn"), py::arg("dofs"), py::arg("iip"))

  .def("nelem", &GooseFEM::MatrixDiagonalPartitioned::nelem, "Return number of element")
  .def("nne"  , &GooseFEM::MatrixDiagonalPartitioned::nne  , "Return number of nodes per element")
  .def("nnode", &GooseFEM::MatrixDiagonalPartitioned::nnode, "Return number of nodes")
  .def("ndim" , &GooseFEM::MatrixDiagonalPartitioned::ndim , "Return number of dimensions")
  .def("ndof" , &GooseFEM::MatrixDiagonalPartitioned::ndof , "Return number of degrees-of-freedom")
  .def("nnu"  , &GooseFEM::MatrixDiagonalPartitioned::nnu  , "Return number of unknown degrees-of-freedom")
  .def("nnp"  , &GooseFEM::MatrixDiagonalPartitioned::nnp  , "Return number of prescribed degrees-of-freedom")

  .def("assemble", &GooseFEM::MatrixDiagonalPartitioned::assemble, "Assemble matrix from 'elemmat", py::arg("elemmat"))

  .def("dofs" , &GooseFEM::MatrixDiagonalPartitioned::dofs , "Return degrees-of-freedom")
  .def("iiu"  , &GooseFEM::MatrixDiagonalPartitioned::iiu  , "Return unknown degrees-of-freedom")
  .def("iip"  , &GooseFEM::MatrixDiagonalPartitioned::iip  , "Return prescribed degrees-of-freedom")

  .def("dot"  , py::overload_cast<const xt::xtensor<double,1>&                              >(&GooseFEM::MatrixDiagonalPartitioned::dot  , py::const_), "Dot product 'b_i = A_ij * x_j"                                              , py::arg("x"))
  .def("dot_u", py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(&GooseFEM::MatrixDiagonalPartitioned::dot_u, py::const_), "Dot product 'b_i = A_ij * x_j (b_u = A_uu * x_u + A_up * x_p == A_uu * x_u)", py::arg("x_u"), py::arg("x_p"))
  .def("dot_p", py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(&GooseFEM::MatrixDiagonalPartitioned::dot_p, py::const_), "Dot product 'b_i = A_ij * x_j (b_p = A_pu * x_u + A_pp * x_p == A_pp * x_p)", py::arg("x_u"), py::arg("x_p"))

  .def("asDiagonal", &GooseFEM::MatrixDiagonalPartitioned::asDiagonal, "Return as diagonal matrix (column)")

  .def("__repr__", [](const GooseFEM::MatrixDiagonalPartitioned &){ return "<GooseFEM.MatrixDiagonalPartitioned>"; });

}

// =================================================================================================


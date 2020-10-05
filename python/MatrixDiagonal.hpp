/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_MatrixDiagonal(py::module& m)
{

    py::class_<GooseFEM::MatrixDiagonal>(m, "MatrixDiagonal")

        .def(
            py::init<const xt::xtensor<size_t, 2>&, const xt::xtensor<size_t, 2>&>(),
            "Diagonal matrix",
            py::arg("conn"),
            py::arg("dofs"))

        .def("nelem", &GooseFEM::MatrixDiagonal::nelem, "Number of element")

        .def("nne", &GooseFEM::MatrixDiagonal::nne, "Number of nodes per element")

        .def("nnode", &GooseFEM::MatrixDiagonal::nnode, "Number of nodes")

        .def("ndim", &GooseFEM::MatrixDiagonal::ndim, "Number of dimensions")

        .def("ndof", &GooseFEM::MatrixDiagonal::ndof, "Number of degrees-of-freedom")

        .def("set", &GooseFEM::MatrixDiagonal::set, "Set matrix components", py::arg("A"))

        .def(
            "assemble",
            &GooseFEM::MatrixDiagonal::assemble,
            "Assemble matrix from 'elemmat",
            py::arg("elemmat"))

        .def("dofs", &GooseFEM::MatrixDiagonal::dofs, "Return degrees-of-freedom")

        .def(
            "AsDiagonal",
            &GooseFEM::MatrixDiagonal::AsDiagonal,
            "Return as diagonal matrix (column)")

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonal::Dot, py::const_),
            "Dot product 'b_i = A_ij * x_j",
            py::arg("x"))

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixDiagonal::Dot, py::const_),
            "Dot product 'b_i = A_ij * x_j",
            py::arg("x"))

        .def(
            "Solve",
            py::overload_cast<const xt::xtensor<double, 1>&>(&GooseFEM::MatrixDiagonal::Solve),
            "Solve",
            py::arg("b"))

        .def(
            "Solve",
            py::overload_cast<const xt::xtensor<double, 2>&>(&GooseFEM::MatrixDiagonal::Solve),
            "Solve",
            py::arg("b"))

        .def("__repr__", [](const GooseFEM::MatrixDiagonal&) {
            return "<GooseFEM.MatrixDiagonal>";
        });
}

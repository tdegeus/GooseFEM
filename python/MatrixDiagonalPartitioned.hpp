/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_MatrixDiagonalPartitioned(py::module& m)
{

    py::class_<GooseFEM::MatrixDiagonalPartitioned>(m, "MatrixDiagonalPartitioned")

        .def(
            py::init<
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 1>&>(),
            "Diagonal, partitioned, matrix",
            py::arg("conn"),
            py::arg("dofs"),
            py::arg("iip"))

        .def("nelem", &GooseFEM::MatrixDiagonalPartitioned::nelem, "Number of element")

        .def("nne", &GooseFEM::MatrixDiagonalPartitioned::nne, "Number of nodes per element")

        .def("nnode", &GooseFEM::MatrixDiagonalPartitioned::nnode, "Number of nodes")

        .def("ndim", &GooseFEM::MatrixDiagonalPartitioned::ndim, "Number of dimensions")

        .def("ndof", &GooseFEM::MatrixDiagonalPartitioned::ndof, "Number of degrees-of-freedom")

        .def(
            "nnu",
            &GooseFEM::MatrixDiagonalPartitioned::nnu,
            "Number of unknown degrees-of-freedom")

        .def(
            "nnp",
            &GooseFEM::MatrixDiagonalPartitioned::nnp,
            "Number of prescribed degrees-of-freedom")

        .def(
            "assemble",
            &GooseFEM::MatrixDiagonalPartitioned::assemble,
            "Assemble matrix from 'elemmat",
            py::arg("elemmat"))

        .def("dofs", &GooseFEM::MatrixDiagonalPartitioned::dofs, "Return degrees-of-freedom")

        .def("iiu", &GooseFEM::MatrixDiagonalPartitioned::iiu, "Return unknown degrees-of-freedom")

        .def(
            "iip",
            &GooseFEM::MatrixDiagonalPartitioned::iip,
            "Return prescribed degrees-of-freedom")

        .def(
            "AsDiagonal",
            &GooseFEM::MatrixDiagonalPartitioned::AsDiagonal,
            "Return as diagonal matrix (column)")

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Dot, py::const_),
            "Dot product 'b_i = A_ij * x_j",
            py::arg("x"))

        .def(
            "Dot_u",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Dot_u, py::const_),
            "Dot product 'b_i = A_ij * x_j (b_u = A_uu * x_u + A_up * x_p == A_uu * x_u)",
            py::arg("x_u"),
            py::arg("x_p"))

        .def(
            "Dot_p",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Dot_p, py::const_),
            "Dot product 'b_i = A_ij * x_j (b_p = A_pu * x_u + A_pp * x_p == A_pp * x_p)",
            py::arg("x_u"),
            py::arg("x_p"))

        .def(
            "Solve",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Solve),
            "Solve",
            py::arg("b"),
            py::arg("x"))

        .def(
            "Solve",
            py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Solve),
            "Solve",
            py::arg("b"),
            py::arg("x"))

        .def(
            "Solve_u",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Solve_u),
            "Solve_u",
            py::arg("b_u"),
            py::arg("x_p"))

        .def(
            "Reaction",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Reaction, py::const_),
            "Reaction",
            py::arg("x"),
            py::arg("b"))

        .def(
            "Reaction",
            py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Reaction, py::const_),
            "Reaction",
            py::arg("x"),
            py::arg("b"))

        .def(
            "Reaction_p",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonalPartitioned::Reaction_p, py::const_),
            "Reaction_p",
            py::arg("x_u"),
            py::arg("x_p"))

        .def("__repr__", [](const GooseFEM::MatrixDiagonalPartitioned&) {
            return "<GooseFEM.MatrixDiagonalPartitioned>";
        });
}

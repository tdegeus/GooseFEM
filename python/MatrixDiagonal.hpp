/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>
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

        .def("nelem",
             &GooseFEM::MatrixDiagonal::nelem,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::nelem`.")

        .def("nne",
             &GooseFEM::MatrixDiagonal::nne,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::nne`.")

        .def("nnode",
             &GooseFEM::MatrixDiagonal::nnode,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::nnode`.")

        .def("ndim",
             &GooseFEM::MatrixDiagonal::ndim,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::ndim`.")

        .def("ndof",
             &GooseFEM::MatrixDiagonal::ndof,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::ndof`.")

        .def("dofs",
             &GooseFEM::MatrixDiagonal::dofs,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::dofs`.")

        .def("assemble",
             &GooseFEM::MatrixDiagonal::assemble,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::assemble`.",
             py::arg("elemmat"))

        .def("set",
             &GooseFEM::MatrixDiagonal::set,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::set`.",
             py::arg("A"))

       .def("Todiagonal",
             &GooseFEM::MatrixDiagonal::Todiagonal,
             "See :cpp:func:`GooseFEM::MatrixDiagonal::Todiagonal`.")

        .def("Dot",
             py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixDiagonal::Dot, py::const_),
             "See :cpp:func:`GooseFEM::MatrixDiagonal::Dot`.",
             py::arg("x"))

        .def("Dot",
             py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixDiagonal::Dot, py::const_),
             "See :cpp:func:`GooseFEM::MatrixDiagonal::Dot`.",
             py::arg("x"))

        .def("Solve",
             py::overload_cast<const xt::xtensor<double, 1>&>(&GooseFEM::MatrixDiagonal::Solve),
             "See :cpp:func:`GooseFEM::MatrixDiagonal::Solve`.",
             py::arg("b"))

        .def("Solve",
             py::overload_cast<const xt::xtensor<double, 2>&>(&GooseFEM::MatrixDiagonal::Solve),
             "See :cpp:func:`GooseFEM::MatrixDiagonal::Solve`.",
             py::arg("b"))

        .def("__repr__", [](const GooseFEM::MatrixDiagonal&) {
            return "<GooseFEM.MatrixDiagonal>";
        });
}

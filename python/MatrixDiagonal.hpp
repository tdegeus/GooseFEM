/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_MATRIXDIAGONAL_H
#define PYGOOSEFEM_MATRIXDIAGONAL_H

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include "Matrix.hpp"

namespace py = pybind11;

template <class C, class P>
void register_Matrix_MatrixDiagonalBase(P& cls)
{
    cls.def(
        "Solve",
        py::overload_cast<const xt::pytensor<double, 1>&>(&C::Solve),
        "Solve",
        py::arg("b"));

    cls.def(
        "Solve",
        py::overload_cast<const xt::pytensor<double, 2>&>(&C::Solve),
        "Solve",
        py::arg("b"));

    cls.def(
        "solve",
        py::overload_cast<const xt::pytensor<double, 1>&, xt::pytensor<double, 1>&>(&C::solve),
        "Solve (write to x)",
        py::arg("b"),
        py::arg("x"));

    cls.def(
        "solve",
        py::overload_cast<const xt::pytensor<double, 2>&, xt::pytensor<double, 2>&>(&C::solve),
        "Solve (write to x)",
        py::arg("b"),
        py::arg("x"));
}

void init_MatrixDiagonal(py::module& m)
{
    py::class_<GooseFEM::MatrixDiagonal> cls(m, "MatrixDiagonal");
    register_Matrix_MatrixBase<GooseFEM::MatrixDiagonal>(cls);
    register_Matrix_MatrixDiagonalBase<GooseFEM::MatrixDiagonal>(cls);

    cls.def(
        py::init<const xt::pytensor<size_t, 2>&, const xt::pytensor<size_t, 2>&>(),
        "See :cpp:class:`GooseFEM::MatrixDiagonal`.",
        py::arg("conn"),
        py::arg("dofs"));

    cls.def("set", &GooseFEM::MatrixDiagonal::set, py::arg("A"));
    cls.def_property_readonly("data", &GooseFEM::MatrixDiagonal::data);

    cls.def(
        "__repr__", [](const GooseFEM::MatrixDiagonal&) { return "<GooseFEM.MatrixDiagonal>"; });
}

#endif

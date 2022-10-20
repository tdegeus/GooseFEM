/**
 * @file
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * \license This project is released under the GNU Public License (GPLv3).
 */

#include <GooseFEM/MatrixDiagonalPartitioned.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_MatrixDiagonalPartitioned(py::module& m)
{
    py::class_<GooseFEM::MatrixDiagonalPartitioned> cls(m, "MatrixDiagonalPartitioned");
    register_Matrix_MatrixBase<GooseFEM::MatrixDiagonalPartitioned>(cls);
    register_Matrix_MatrixPartitionedBase<GooseFEM::MatrixDiagonalPartitioned>(cls);

    cls.def(
        py::init<
            const xt::pytensor<size_t, 2>&,
            const xt::pytensor<size_t, 2>&,
            const xt::pytensor<size_t, 1>&>(),
        "See :cpp:class:`GooseFEM::MatrixDiagonalPartitioned`.",
        py::arg("conn"),
        py::arg("dofs"),
        py::arg("iip"));

    cls.def("data", &GooseFEM::MatrixDiagonalPartitioned::data, "Copy to assemble diagonal matrix");
    cls.def_property_readonly("data_uu", &GooseFEM::MatrixDiagonalPartitioned::data_uu);
    cls.def_property_readonly("data_pp", &GooseFEM::MatrixDiagonalPartitioned::data_pp);

    cls.def(
        "Dot_u",
        py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
            &GooseFEM::MatrixDiagonalPartitioned::Dot_u, py::const_),
        py::arg("x_u"),
        py::arg("x_p"));

    cls.def(
        "Dot_p",
        py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
            &GooseFEM::MatrixDiagonalPartitioned::Dot_p, py::const_),
        py::arg("x_u"),
        py::arg("x_p"));

    cls.def(
        "Solve_u",
        py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
            &GooseFEM::MatrixDiagonalPartitioned::Solve_u),
        py::arg("b_u"),
        py::arg("x_p"));

    cls.def("__repr__", [](const GooseFEM::MatrixDiagonalPartitioned&) {
        return "<GooseFEM.MatrixDiagonalPartitioned>";
    });
}

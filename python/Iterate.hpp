/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_ITERATE_H
#define PYGOOSEFEM_ITERATE_H

#include <GooseFEM/Iterate.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_Iterate(py::module& mod)
{
    py::class_<GooseFEM::Iterate::StopList> cls(mod, "Iterate");

    cls.def(
        py::init<size_t>(),
        "Class to perform a residual check based on the last 'n' iterations."
        "See :cpp:class:`GooseFEM::Iterate::StopList`.",
        py::arg("n") = 1);

    cls.def(
        "reset",
        py::overload_cast<>(&GooseFEM::Iterate::StopList::reset),
        "Reset."
        "See :cpp:func:`GooseFEM::Iterate::StopList::reset`.");

    cls.def(
        "roll_insert",
        &GooseFEM::Iterate::StopList::roll_insert,
        "Roll the list with the residuals, and add a new residual to the end."
        "See :cpp:func:`GooseFEM::Iterate::StopList::roll_insert`.",
        py::arg("res"));

    cls.def(
        "descending",
        &GooseFEM::Iterate::StopList::descending,
        "Check of the sequence of `n` residuals is in descending order."
        "See :cpp:func:`GooseFEM::Iterate::StopList::descending`.");

    cls.def(
        "all_less",
        &GooseFEM::Iterate::StopList::all_less,
        "Check of the sequence of `n` residuals are all below a tolerance."
        "See :cpp:func:`GooseFEM::Iterate::StopList::all_less`.",
        py::arg("tol"));

    cls.def(
        "stop_simple",
        &GooseFEM::Iterate::StopList::stop_simple,
        "Update list of residuals, return `true` if all residuals are below the tolerance."
        "See :cpp:func:`GooseFEM::Iterate::StopList::stop_simple`.",
        py::arg("res"),
        py::arg("tol"));

    cls.def(
        "stop",
        &GooseFEM::Iterate::StopList::stop,
        "Update list of residuals, return `true` if all residuals are sorted and below the "
        "tolerance."
        "See :cpp:func:`GooseFEM::Iterate::StopList::stop`.",
        py::arg("res"),
        py::arg("tol"));

    cls.def(
        "get",
        &GooseFEM::Iterate::StopList::get,
        "Get the list of residuals."
        "See :cpp:func:`GooseFEM::Iterate::StopList::get`.");

    cls.def("__repr__", [](const GooseFEM::Iterate::StopList&) {
        return "<GooseFEM.Iterate.StopList>";
    });
}

#endif

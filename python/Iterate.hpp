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
    py::class_<GooseFEM::Iterate::StopList> cls(mod, "StopList");

    cls.def(py::init<size_t>(), "See :cpp:class:`GooseFEM::Iterate::StopList`.", py::arg("n") = 1);

    cls.def("reset", py::overload_cast<>(&GooseFEM::Iterate::StopList::reset));
    cls.def("roll_insert", &GooseFEM::Iterate::StopList::roll_insert, py::arg("res"));
    cls.def("descending", &GooseFEM::Iterate::StopList::descending);
    cls.def("all_less", &GooseFEM::Iterate::StopList::all_less, py::arg("tol"));
    cls.def_property_readonly("data", &GooseFEM::Iterate::StopList::data);

    cls.def("__repr__", [](const GooseFEM::Iterate::StopList&) {
        return "<GooseFEM.Iterate.StopList>";
    });
}

#endif

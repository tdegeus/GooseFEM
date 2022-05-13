/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_TYINGSPERIODIC_H
#define PYGOOSEFEM_TYINGSPERIODIC_H

#include <GooseFEM/TyingsPeriodic.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_TyingsPeriodic(py::module& mod)
{
    {
        py::class_<GooseFEM::Tyings::Periodic> cls(mod, "Periodic");

        cls.def(
            py::init<
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 2>&>(),
            "See :cpp:class:`GooseFEM::Tyings::Periodic`.",
            py::arg("coor"),
            py::arg("dofs"),
            py::arg("control_dofs"),
            py::arg("nodal_tyings"));

        cls.def(
            py::init<
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 1>&>(),
            "See :cpp:class:`GooseFEM::Tyings::Periodic`.",
            py::arg("coor"),
            py::arg("dofs"),
            py::arg("control_dofs"),
            py::arg("nodal_tyings"),
            py::arg("iip"));

        cls.def("nnd", &GooseFEM::Tyings::Periodic::nnd);
        cls.def("nni", &GooseFEM::Tyings::Periodic::nni);
        cls.def("nnu", &GooseFEM::Tyings::Periodic::nnu);
        cls.def("nnp", &GooseFEM::Tyings::Periodic::nnp);
        cls.def("dofs", &GooseFEM::Tyings::Periodic::dofs);
        cls.def("control", &GooseFEM::Tyings::Periodic::control);
        cls.def("nodal_tyings", &GooseFEM::Tyings::Periodic::nodal_tyings);
        cls.def("iid", &GooseFEM::Tyings::Periodic::iid);
        cls.def("iii", &GooseFEM::Tyings::Periodic::iii);
        cls.def("iiu", &GooseFEM::Tyings::Periodic::iiu);
        cls.def("iip", &GooseFEM::Tyings::Periodic::iip);
        cls.def("Cdi", &GooseFEM::Tyings::Periodic::Cdi);
        cls.def("Cdu", &GooseFEM::Tyings::Periodic::Cdu);
        cls.def("Cdp", &GooseFEM::Tyings::Periodic::Cdp);

        cls.def("__repr__", [](const GooseFEM::Tyings::Periodic&) {
            return "<GooseFEM.Tyings.Periodic>";
        });
    }

    {
        py::class_<GooseFEM::Tyings::Control> cls(mod, "Control");

        cls.def(
            py::init<const xt::pytensor<double, 2>&, const xt::pytensor<size_t, 2>&>(),
            "See :cpp:class:`GooseFEM::Tyings::Control`.",
            py::arg("coor"),
            py::arg("dofs"));

        cls.def("coor", &GooseFEM::Tyings::Control::coor);
        cls.def("dofs", &GooseFEM::Tyings::Control::dofs);
        cls.def("controlDofs", &GooseFEM::Tyings::Control::controlDofs);
        cls.def("controlNodes", &GooseFEM::Tyings::Control::controlNodes);

        cls.def("__repr__", [](const GooseFEM::Tyings::Control&) {
            return "<GooseFEM.Tyings.Control>";
        });
    }
}

#endif

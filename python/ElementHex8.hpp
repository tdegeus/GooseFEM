/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_ELEMENTHEX8_H
#define PYGOOSEFEM_ELEMENTHEX8_H

#include <GooseFEM/ElementHex8.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pytensor.hpp>

#include "Element.hpp"

namespace py = pybind11;

void init_ElementHex8(py::module& m)
{
    py::class_<GooseFEM::Element::Hex8::Quadrature> cls(m, "Quadrature");

    cls.def(
        py::init<const xt::pytensor<double, 3>&>(),
        "See :cpp:class:`GooseFEM::Element::Hex8::Quadrature`.",
        py::arg("x"));

    cls.def(
        py::init<
            const xt::pytensor<double, 3>&,
            const xt::pytensor<double, 2>&,
            const xt::pytensor<double, 1>&>(),
        "See :cpp:class:`GooseFEM::Element::Hex8::Quadrature`.",
        py::arg("x"),
        py::arg("xi"),
        py::arg("w"));

    register_Mesh_QuadratureBase<GooseFEM::Element::Hex8::Quadrature>(cls);
    register_Mesh_QuadratureBaseCartesian<GooseFEM::Element::Hex8::Quadrature>(cls);

    cls.def_property_readonly(
        "GradN",
        &GooseFEM::Element::Hex8::Quadrature::GradN,
        "Shape function gradients  [nelem, nip, nne, ndim]");

    cls.def("__repr__", [](const GooseFEM::Element::Hex8::Quadrature&) {
        return "<GooseFEM.Element.Hex8.Quadrature>";
    });
}

void init_ElementHex8Gauss(py::module& m)
{
    m.def("nip", &GooseFEM::Element::Hex8::Gauss::nip);
    m.def("xi", &GooseFEM::Element::Hex8::Gauss::xi);
    m.def("w", &GooseFEM::Element::Hex8::Gauss::w);
}

void init_ElementHex8Nodal(py::module& m)
{
    m.def("nip", &GooseFEM::Element::Hex8::Nodal::nip);
    m.def("xi", &GooseFEM::Element::Hex8::Nodal::xi);
    m.def("w", &GooseFEM::Element::Hex8::Nodal::w);
}

#endif

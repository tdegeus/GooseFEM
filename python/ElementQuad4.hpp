/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_ELEMENTQUAD4_H
#define PYGOOSEFEM_ELEMENTQUAD4_H

#include <GooseFEM/ElementQuad4.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pytensor.hpp>

#include "Element.hpp"

namespace py = pybind11;

void init_ElementQuad4(py::module& m)
{
    py::class_<GooseFEM::Element::Quad4::Quadrature> cls(m, "Quadrature");

    cls.def(py::init<const xt::pytensor<double, 3>&>(), "Quadrature", py::arg("x"));

    cls.def(
        py::init<
            const xt::pytensor<double, 3>&,
            const xt::pytensor<double, 2>&,
            const xt::pytensor<double, 1>&>(),
        "Quadrature",
        py::arg("x"),
        py::arg("xi"),
        py::arg("w"));

    register_Element_QuadratureBase<GooseFEM::Element::Quad4::Quadrature>(cls);
    register_Element_QuadratureBaseCartesian<GooseFEM::Element::Quad4::Quadrature>(cls);

    cls.def("GradN", &GooseFEM::Element::Quad4::Quadrature::GradN, "Shape function gradients");

    cls.def("__repr__", [](const GooseFEM::Element::Quad4::Quadrature&) {
        return "<GooseFEM.Element.Quad4.Quadrature>";
    });
}

void init_ElementQuad4Gauss(py::module& m)
{
    m.def("nip", &GooseFEM::Element::Quad4::Gauss::nip, "Return number of integration point");
    m.def("xi", &GooseFEM::Element::Quad4::Gauss::xi, "Return integration point coordinates");
    m.def("w", &GooseFEM::Element::Quad4::Gauss::w, "Return integration point weights");
}

void init_ElementQuad4Nodal(py::module& m)
{
    m.def("nip", &GooseFEM::Element::Quad4::Nodal::nip, "Return number of integration point");
    m.def("xi", &GooseFEM::Element::Quad4::Nodal::xi, "Return integration point coordinates");
    m.def("w", &GooseFEM::Element::Quad4::Nodal::w, "Return integration point weights");
}

void init_ElementQuad4MidPoint(py::module& m)
{
    m.def("nip", &GooseFEM::Element::Quad4::MidPoint::nip, "Return number of integration point");
    m.def("xi", &GooseFEM::Element::Quad4::MidPoint::xi, "Return integration point coordinates");
    m.def("w", &GooseFEM::Element::Quad4::MidPoint::w, "Return integration point weights");
}

#endif

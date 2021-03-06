/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_ELEMENTQUAD4PLANAR_H
#define PYGOOSEFEM_ELEMENTQUAD4PLANAR_H

#include <GooseFEM/ElementQuad4Planar.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pytensor.hpp>

#include "Element.hpp"

namespace py = pybind11;

void init_ElementQuad4Planar(py::module& m)
{
    py::class_<GooseFEM::Element::Quad4::QuadraturePlanar> cls(m, "QuadraturePlanar");

    cls.def(py::init<const xt::pytensor<double, 3>&, double>(), "QuadraturePlanar", py::arg("x"), py::arg("thick") = 1.0);

    cls.def(py::init<const xt::pytensor<double, 3>&,
                     const xt::pytensor<double, 2>&,
                     const xt::pytensor<double, 1>&,
                     double>(),
            "QuadraturePlanar",
            py::arg("x"),
            py::arg("xi"),
            py::arg("w"),
            py::arg("thick") = 1.0);

    register_Element_QuadratureBase<GooseFEM::Element::Quad4::QuadraturePlanar>(cls);
    register_Element_QuadratureBaseCartesian<GooseFEM::Element::Quad4::QuadraturePlanar>(cls);

    cls.def("GradN", &GooseFEM::Element::Quad4::QuadraturePlanar::GradN, "Shape function gradients");

    cls.def("__repr__", [](const GooseFEM::Element::Quad4::QuadraturePlanar&) {
        return "<GooseFEM.Element.Quad4.QuadraturePlanar>";
    });
}

#endif

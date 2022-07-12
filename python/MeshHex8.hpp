/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_MESHHEX8_H
#define PYGOOSEFEM_MESHHEX8_H

#include <GooseFEM/MeshHex8.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include "Mesh.hpp"

namespace py = pybind11;

void init_MeshHex8(py::module& mod)
{
    {
        py::class_<GooseFEM::Mesh::Hex8::Regular> cls(mod, "Regular");

        cls.def(
            py::init<size_t, size_t, size_t, double>(),
            "See :cpp:class:`GooseFEM::Mesh::Hex8::Regular`.",
            py::arg("nx"),
            py::arg("ny"),
            py::arg("nz"),
            py::arg("h") = 1.);

        register_Element_RegularBase<GooseFEM::Mesh::Hex8::Regular>(cls);
        register_Element_RegularBase3d<GooseFEM::Mesh::Hex8::Regular>(cls);

        cls.def("__repr__", [](const GooseFEM::Mesh::Hex8::Regular&) {
            return "<GooseFEM.Mesh.Hex8.Regular>";
        });
    }

    {
        py::class_<GooseFEM::Mesh::Hex8::FineLayer> cls(mod, "FineLayer");

        cls.def(
            py::init<size_t, size_t, size_t, double, size_t>(),
            "See :cpp:class:`GooseFEM::Mesh::Hex8::FineLayer`.",
            py::arg("nx"),
            py::arg("ny"),
            py::arg("nz"),
            py::arg("h") = 1.0,
            py::arg("nfine") = 1);

        cls.def_property_readonly(
            "elementsMiddleLayer", &GooseFEM::Mesh::Hex8::FineLayer::elementsMiddleLayer);

        cls.def("__repr__", [](const GooseFEM::Mesh::Hex8::FineLayer&) {
            return "<GooseFEM.Mesh.Hex8.FineLayer>";
        });
    }
}

#endif

/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_MESHTRI3_H
#define PYGOOSEFEM_MESHTRI3_H

#include <GooseFEM/MeshTri3.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>
#include <pyxtensor/pyxtensor.hpp>

#include "Mesh.hpp"

namespace py = pybind11;

void init_MeshTri3(py::module& mod)
{
    py::class_<GooseFEM::Mesh::Tri3::Regular> cls(mod, "Regular");

    cls.def(py::init<size_t, size_t, double>(),
            "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
            py::arg("nx"),
            py::arg("ny"),
            py::arg("h") = 1.);

    register_Element_RegularBase<GooseFEM::Mesh::Tri3::Regular>(cls);
    register_Element_RegularBase2d<GooseFEM::Mesh::Tri3::Regular>(cls);

    cls.def("__repr__", [](const GooseFEM::Mesh::Tri3::Regular&) {
        return "<GooseFEM.Mesh.Tri3.Regular>";
    });

    mod.def("getOrientation",
            &GooseFEM::Mesh::Tri3::getOrientation,
            "Get the orientation of each element",
            py::arg("coor"),
            py::arg("conn"));

    mod.def("setOrientation",
            py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&, int>(
                &GooseFEM::Mesh::Tri3::setOrientation),
            "Set the orientation of each element",
            py::arg("coor"),
            py::arg("conn"),
            py::arg("orientation"));

    mod.def("setOrientation",
            py::overload_cast<
                const xt::xtensor<double, 2>&,
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<int, 1>&,
                int>(&GooseFEM::Mesh::Tri3::setOrientation),
            "Set the orientation of each element",
            py::arg("coor"),
            py::arg("conn"),
            py::arg("val"),
            py::arg("orientation"));
}

#endif

/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef PYGOOSEFEM_MESHQUAD4_H
#define PYGOOSEFEM_MESHQUAD4_H

#include <GooseFEM/MeshQuad4.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include "Mesh.hpp"

namespace py = pybind11;

void init_MeshQuad4(py::module& m)
{
    {
        py::class_<GooseFEM::Mesh::Quad4::Regular> cls(m, "Regular");

        cls.def(
            py::init<size_t, size_t, double>(),
            "See :cpp:class:`GooseFEM::Mesh::Quad4::Regular`.",
            py::arg("nx"),
            py::arg("ny"),
            py::arg("h") = 1.0);

        register_Element_RegularBase<GooseFEM::Mesh::Quad4::Regular>(cls);
        register_Element_RegularBase2d<GooseFEM::Mesh::Quad4::Regular>(cls);

        cls.def("elementgrid", &GooseFEM::Mesh::Quad4::Regular::elementgrid);

        cls.def("__repr__", [](const GooseFEM::Mesh::Quad4::Regular&) {
            return "<GooseFEM.Mesh.Quad4.Regular>";
        });
    }

    {
        py::class_<GooseFEM::Mesh::Quad4::FineLayer> cls(m, "FineLayer");

        cls.def(
            py::init<size_t, size_t, double, size_t>(),
            "See :cpp:class:`GooseFEM::Mesh::Quad4::FineLayer`.",
            py::arg("nx"),
            py::arg("ny"),
            py::arg("h") = 1.,
            py::arg("nfine") = 1);

        cls.def(
            py::init<const xt::pytensor<double, 2>&, const xt::pytensor<size_t, 2>&>(),
            "See :cpp:class:`GooseFEM::Mesh::Quad4::FineLayer`.",
            py::arg("coor"),
            py::arg("conn"));

        register_Element_RegularBase<GooseFEM::Mesh::Quad4::FineLayer>(cls);
        register_Element_RegularBase2d<GooseFEM::Mesh::Quad4::FineLayer>(cls);

        cls.def("elemrow_nhx", &GooseFEM::Mesh::Quad4::FineLayer::elemrow_nhx);
        cls.def("elemrow_nhy", &GooseFEM::Mesh::Quad4::FineLayer::elemrow_nhy);
        cls.def("elemrow_type", &GooseFEM::Mesh::Quad4::FineLayer::elemrow_type);
        cls.def("elemrow_nelem", &GooseFEM::Mesh::Quad4::FineLayer::elemrow_nelem);
        cls.def("elementsMiddleLayer", &GooseFEM::Mesh::Quad4::FineLayer::elementsMiddleLayer);
        cls.def("elementsLayer", &GooseFEM::Mesh::Quad4::FineLayer::elementsLayer);

        cls.def(
            "elementgrid_ravel",
            &GooseFEM::Mesh::Quad4::FineLayer::elementgrid_ravel,
            py::arg("rows_range"),
            py::arg("cols_range"));

        cls.def(
            "elementgrid_around_ravel",
            &GooseFEM::Mesh::Quad4::FineLayer::elementgrid_around_ravel,
            py::arg("element"),
            py::arg("size"),
            py::arg("periodic") = true);

        cls.def(
            "elementgrid_leftright",
            &GooseFEM::Mesh::Quad4::FineLayer::elementgrid_leftright,
            py::arg("element"),
            py::arg("left"),
            py::arg("right"),
            py::arg("periodic") = true);

        cls.def("roll", &GooseFEM::Mesh::Quad4::FineLayer::roll);

        cls.def("__repr__", [](const GooseFEM::Mesh::Quad4::FineLayer&) {
            return "<GooseFEM.Mesh.Quad4.FineLayer>";
        });
    }
}

void init_MeshQuad4Map(py::module& m)
{

    py::class_<GooseFEM::Mesh::Quad4::Map::RefineRegular>(m, "RefineRegular")

        .def(
            py::init<const GooseFEM::Mesh::Quad4::Regular&, size_t, size_t>(),
            "See :cpp:class:`GooseFEM::Mesh::Quad4::Map::RefineRegular`.",
            py::arg("mesh"),
            py::arg("nx"),
            py::arg("ny"))

        .def("getCoarseMesh", &GooseFEM::Mesh::Quad4::Map::RefineRegular::getCoarseMesh)
        .def("getFineMesh", &GooseFEM::Mesh::Quad4::Map::RefineRegular::getFineMesh)
        .def("getMap", &GooseFEM::Mesh::Quad4::Map::RefineRegular::getMap)

        .def(
            "meanToCoarse",
            py::overload_cast<const xt::pytensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::meanToCoarse<double, 1>, py::const_))

        .def(
            "meanToCoarse",
            py::overload_cast<const xt::pytensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::meanToCoarse<double, 2>, py::const_))

        .def(
            "meanToCoarse",
            py::overload_cast<const xt::pytensor<double, 3>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::meanToCoarse<double, 3>, py::const_))

        .def(
            "meanToCoarse",
            py::overload_cast<const xt::pytensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::meanToCoarse<double, 4>, py::const_))

        .def(
            "averageToCoarse",
            py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::averageToCoarse<double, 1, double>,
                py::const_))

        .def(
            "averageToCoarse",
            py::overload_cast<const xt::pytensor<double, 2>&, const xt::pytensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::averageToCoarse<double, 2, double>,
                py::const_))

        .def(
            "averageToCoarse",
            py::overload_cast<const xt::pytensor<double, 3>&, const xt::pytensor<double, 3>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::averageToCoarse<double, 3, double>,
                py::const_))

        .def(
            "averageToCoarse",
            py::overload_cast<const xt::pytensor<double, 4>&, const xt::pytensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::averageToCoarse<double, 4, double>,
                py::const_))

        .def(
            "mapToFine",
            py::overload_cast<const xt::pytensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToFine<double, 1>, py::const_))

        .def(
            "mapToFine",
            py::overload_cast<const xt::pytensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToFine<double, 2>, py::const_))

        .def(
            "mapToFine",
            py::overload_cast<const xt::pytensor<double, 3>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToFine<double, 3>, py::const_))

        .def(
            "mapToFine",
            py::overload_cast<const xt::pytensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToFine<double, 4>, py::const_))

        .def("__repr__", [](const GooseFEM::Mesh::Quad4::Map::RefineRegular&) {
            return "<GooseFEM.Mesh.Quad4.Map.RefineRegular>";
        });

    py::class_<GooseFEM::Mesh::Quad4::Map::FineLayer2Regular>(m, "FineLayer2Regular")

        .def(
            py::init<const GooseFEM::Mesh::Quad4::FineLayer&>(),
            "See :cpp:class:`GooseFEM::Mesh::Quad4::Map::FineLayer2Regular`.",
            py::arg("mesh"))

        .def("getRegularMesh", &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::getRegularMesh)
        .def("getFineLayerMesh", &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::getFineLayerMesh)
        .def("getMap", &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::getMap)
        .def("getMapFraction", &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::getMapFraction)

        .def(
            "mapToRegular",
            py::overload_cast<const xt::pytensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular<double, 1>,
                py::const_))

        .def(
            "mapToRegular",
            py::overload_cast<const xt::pytensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular<double, 2>,
                py::const_))

        .def(
            "mapToRegular",
            py::overload_cast<const xt::pytensor<double, 3>&>(
                &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular<double, 3>,
                py::const_))

        .def(
            "mapToRegular",
            py::overload_cast<const xt::pytensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular<double, 4>,
                py::const_))

        .def("__repr__", [](const GooseFEM::Mesh::Quad4::Map::FineLayer2Regular&) {
            return "<GooseFEM.Mesh.Quad4.Map.FineLayer2Regular>";
        });
}

#endif

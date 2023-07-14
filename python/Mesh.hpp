/**
 * @file
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#ifndef PYGOOSEFEM_MESH_H
#define PYGOOSEFEM_MESH_H

#include <GooseFEM/Mesh.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

template <class C, class P>
void register_Mesh_RegularBase(P& cls)
{
    cls.def_property_readonly("nelem", &C::nelem);
    cls.def_property_readonly("nnode", &C::nnode);
    cls.def_property_readonly("nne", &C::nne);
    cls.def_property_readonly("ndim", &C::ndim);
    cls.def_property_readonly("nelx", &C::nelx);
    cls.def_property_readonly("nely", &C::nely);
    cls.def_property_readonly("h", &C::h);
    cls.def_property_readonly("elementType", &C::getElementType);
    cls.def_property_readonly("coor", &C::coor);
    cls.def_property_readonly("conn", &C::conn);
    cls.def_property_readonly("dofs", &C::dofs);
    cls.def_property_readonly("dofsPeriodic", &C::dofsPeriodic);
    cls.def_property_readonly("nodesPeriodic", &C::nodesPeriodic);
    cls.def_property_readonly("nodesOrigin", &C::nodesOrigin);
}

template <class C, class P>
void register_Mesh_RegularBase2d(P& cls)
{
    cls.def_property_readonly("nodesBottomEdge", &C::nodesBottomEdge);
    cls.def_property_readonly("nodesTopEdge", &C::nodesTopEdge);
    cls.def_property_readonly("nodesLeftEdge", &C::nodesLeftEdge);
    cls.def_property_readonly("nodesRightEdge", &C::nodesRightEdge);
    cls.def_property_readonly("nodesBottomOpenEdge", &C::nodesBottomOpenEdge);
    cls.def_property_readonly("nodesTopOpenEdge", &C::nodesTopOpenEdge);
    cls.def_property_readonly("nodesLeftOpenEdge", &C::nodesLeftOpenEdge);
    cls.def_property_readonly("nodesRightOpenEdge", &C::nodesRightOpenEdge);
    cls.def_property_readonly("nodesBottomLeftCorner", &C::nodesBottomLeftCorner);
    cls.def_property_readonly("nodesBottomRightCorner", &C::nodesBottomRightCorner);
    cls.def_property_readonly("nodesTopLeftCorner", &C::nodesTopLeftCorner);
    cls.def_property_readonly("nodesTopRightCorner", &C::nodesTopRightCorner);
    cls.def_property_readonly("nodesLeftBottomCorner", &C::nodesLeftBottomCorner);
    cls.def_property_readonly("nodesLeftTopCorner", &C::nodesLeftTopCorner);
    cls.def_property_readonly("nodesRightBottomCorner", &C::nodesRightBottomCorner);
    cls.def_property_readonly("nodesRightTopCorner", &C::nodesRightTopCorner);
}

template <class C, class P>
void register_Mesh_RegularBase3d(P& cls)
{
    cls.def_property_readonly("nodesFront", &C::nodesFront);
    cls.def_property_readonly("nodesBack", &C::nodesBack);
    cls.def_property_readonly("nodesLeft", &C::nodesLeft);
    cls.def_property_readonly("nodesRight", &C::nodesRight);
    cls.def_property_readonly("nodesBottom", &C::nodesBottom);
    cls.def_property_readonly("nodesTop", &C::nodesTop);
    cls.def_property_readonly("nodesFrontFace", &C::nodesFrontFace);
    cls.def_property_readonly("nodesBackFace", &C::nodesBackFace);
    cls.def_property_readonly("nodesLeftFace", &C::nodesLeftFace);
    cls.def_property_readonly("nodesRightFace", &C::nodesRightFace);
    cls.def_property_readonly("nodesBottomFace", &C::nodesBottomFace);
    cls.def_property_readonly("nodesTopFace", &C::nodesTopFace);
    cls.def_property_readonly("nodesFrontBottomEdge", &C::nodesFrontBottomEdge);
    cls.def_property_readonly("nodesFrontTopEdge", &C::nodesFrontTopEdge);
    cls.def_property_readonly("nodesFrontLeftEdge", &C::nodesFrontLeftEdge);
    cls.def_property_readonly("nodesFrontRightEdge", &C::nodesFrontRightEdge);
    cls.def_property_readonly("nodesBackBottomEdge", &C::nodesBackBottomEdge);
    cls.def_property_readonly("nodesBackTopEdge", &C::nodesBackTopEdge);
    cls.def_property_readonly("nodesBackLeftEdge", &C::nodesBackLeftEdge);
    cls.def_property_readonly("nodesBackRightEdge", &C::nodesBackRightEdge);
    cls.def_property_readonly("nodesBottomLeftEdge", &C::nodesBottomLeftEdge);
    cls.def_property_readonly("nodesBottomRightEdge", &C::nodesBottomRightEdge);
    cls.def_property_readonly("nodesTopLeftEdge", &C::nodesTopLeftEdge);
    cls.def_property_readonly("nodesTopRightEdge", &C::nodesTopRightEdge);
    cls.def_property_readonly("nodesBottomFrontEdge", &C::nodesBottomFrontEdge);
    cls.def_property_readonly("nodesBottomBackEdge", &C::nodesBottomBackEdge);
    cls.def_property_readonly("nodesTopFrontEdge", &C::nodesTopFrontEdge);
    cls.def_property_readonly("nodesTopBackEdge", &C::nodesTopBackEdge);
    cls.def_property_readonly("nodesLeftBottomEdge", &C::nodesLeftBottomEdge);
    cls.def_property_readonly("nodesLeftFrontEdge", &C::nodesLeftFrontEdge);
    cls.def_property_readonly("nodesLeftBackEdge", &C::nodesLeftBackEdge);
    cls.def_property_readonly("nodesLeftTopEdge", &C::nodesLeftTopEdge);
    cls.def_property_readonly("nodesRightBottomEdge", &C::nodesRightBottomEdge);
    cls.def_property_readonly("nodesRightTopEdge", &C::nodesRightTopEdge);
    cls.def_property_readonly("nodesRightFrontEdge", &C::nodesRightFrontEdge);
    cls.def_property_readonly("nodesRightBackEdge", &C::nodesRightBackEdge);
    cls.def_property_readonly("nodesFrontBottomOpenEdge", &C::nodesFrontBottomOpenEdge);
    cls.def_property_readonly("nodesFrontTopOpenEdge", &C::nodesFrontTopOpenEdge);
    cls.def_property_readonly("nodesFrontLeftOpenEdge", &C::nodesFrontLeftOpenEdge);
    cls.def_property_readonly("nodesFrontRightOpenEdge", &C::nodesFrontRightOpenEdge);
    cls.def_property_readonly("nodesBackBottomOpenEdge", &C::nodesBackBottomOpenEdge);
    cls.def_property_readonly("nodesBackTopOpenEdge", &C::nodesBackTopOpenEdge);
    cls.def_property_readonly("nodesBackLeftOpenEdge", &C::nodesBackLeftOpenEdge);
    cls.def_property_readonly("nodesBackRightOpenEdge", &C::nodesBackRightOpenEdge);
    cls.def_property_readonly("nodesBottomLeftOpenEdge", &C::nodesBottomLeftOpenEdge);
    cls.def_property_readonly("nodesBottomRightOpenEdge", &C::nodesBottomRightOpenEdge);
    cls.def_property_readonly("nodesTopLeftOpenEdge", &C::nodesTopLeftOpenEdge);
    cls.def_property_readonly("nodesTopRightOpenEdge", &C::nodesTopRightOpenEdge);
    cls.def_property_readonly("nodesBottomFrontOpenEdge", &C::nodesBottomFrontOpenEdge);
    cls.def_property_readonly("nodesBottomBackOpenEdge", &C::nodesBottomBackOpenEdge);
    cls.def_property_readonly("nodesTopFrontOpenEdge", &C::nodesTopFrontOpenEdge);
    cls.def_property_readonly("nodesTopBackOpenEdge", &C::nodesTopBackOpenEdge);
    cls.def_property_readonly("nodesLeftBottomOpenEdge", &C::nodesLeftBottomOpenEdge);
    cls.def_property_readonly("nodesLeftFrontOpenEdge", &C::nodesLeftFrontOpenEdge);
    cls.def_property_readonly("nodesLeftBackOpenEdge", &C::nodesLeftBackOpenEdge);
    cls.def_property_readonly("nodesLeftTopOpenEdge", &C::nodesLeftTopOpenEdge);
    cls.def_property_readonly("nodesRightBottomOpenEdge", &C::nodesRightBottomOpenEdge);
    cls.def_property_readonly("nodesRightTopOpenEdge", &C::nodesRightTopOpenEdge);
    cls.def_property_readonly("nodesRightFrontOpenEdge", &C::nodesRightFrontOpenEdge);
    cls.def_property_readonly("nodesRightBackOpenEdge", &C::nodesRightBackOpenEdge);
    cls.def_property_readonly("nodesFrontBottomLeftCorner", &C::nodesFrontBottomLeftCorner);
    cls.def_property_readonly("nodesFrontBottomRightCorner", &C::nodesFrontBottomRightCorner);
    cls.def_property_readonly("nodesFrontTopLeftCorner", &C::nodesFrontTopLeftCorner);
    cls.def_property_readonly("nodesFrontTopRightCorner", &C::nodesFrontTopRightCorner);
    cls.def_property_readonly("nodesBackBottomLeftCorner", &C::nodesBackBottomLeftCorner);
    cls.def_property_readonly("nodesBackBottomRightCorner", &C::nodesBackBottomRightCorner);
    cls.def_property_readonly("nodesBackTopLeftCorner", &C::nodesBackTopLeftCorner);
    cls.def_property_readonly("nodesBackTopRightCorner", &C::nodesBackTopRightCorner);
    cls.def_property_readonly("nodesFrontLeftBottomCorner", &C::nodesFrontLeftBottomCorner);
    cls.def_property_readonly("nodesBottomFrontLeftCorner", &C::nodesBottomFrontLeftCorner);
    cls.def_property_readonly("nodesBottomLeftFrontCorner", &C::nodesBottomLeftFrontCorner);
    cls.def_property_readonly("nodesLeftFrontBottomCorner", &C::nodesLeftFrontBottomCorner);
    cls.def_property_readonly("nodesLeftBottomFrontCorner", &C::nodesLeftBottomFrontCorner);
    cls.def_property_readonly("nodesFrontRightBottomCorner", &C::nodesFrontRightBottomCorner);
    cls.def_property_readonly("nodesBottomFrontRightCorner", &C::nodesBottomFrontRightCorner);
    cls.def_property_readonly("nodesBottomRightFrontCorner", &C::nodesBottomRightFrontCorner);
    cls.def_property_readonly("nodesRightFrontBottomCorner", &C::nodesRightFrontBottomCorner);
    cls.def_property_readonly("nodesRightBottomFrontCorner", &C::nodesRightBottomFrontCorner);
    cls.def_property_readonly("nodesFrontLeftTopCorner", &C::nodesFrontLeftTopCorner);
    cls.def_property_readonly("nodesTopFrontLeftCorner", &C::nodesTopFrontLeftCorner);
    cls.def_property_readonly("nodesTopLeftFrontCorner", &C::nodesTopLeftFrontCorner);
    cls.def_property_readonly("nodesLeftFrontTopCorner", &C::nodesLeftFrontTopCorner);
    cls.def_property_readonly("nodesLeftTopFrontCorner", &C::nodesLeftTopFrontCorner);
    cls.def_property_readonly("nodesFrontRightTopCorner", &C::nodesFrontRightTopCorner);
    cls.def_property_readonly("nodesTopFrontRightCorner", &C::nodesTopFrontRightCorner);
    cls.def_property_readonly("nodesTopRightFrontCorner", &C::nodesTopRightFrontCorner);
    cls.def_property_readonly("nodesRightFrontTopCorner", &C::nodesRightFrontTopCorner);
    cls.def_property_readonly("nodesRightTopFrontCorner", &C::nodesRightTopFrontCorner);
    cls.def_property_readonly("nodesBackLeftBottomCorner", &C::nodesBackLeftBottomCorner);
    cls.def_property_readonly("nodesBottomBackLeftCorner", &C::nodesBottomBackLeftCorner);
    cls.def_property_readonly("nodesBottomLeftBackCorner", &C::nodesBottomLeftBackCorner);
    cls.def_property_readonly("nodesLeftBackBottomCorner", &C::nodesLeftBackBottomCorner);
    cls.def_property_readonly("nodesLeftBottomBackCorner", &C::nodesLeftBottomBackCorner);
    cls.def_property_readonly("nodesBackRightBottomCorner", &C::nodesBackRightBottomCorner);
    cls.def_property_readonly("nodesBottomBackRightCorner", &C::nodesBottomBackRightCorner);
    cls.def_property_readonly("nodesBottomRightBackCorner", &C::nodesBottomRightBackCorner);
    cls.def_property_readonly("nodesRightBackBottomCorner", &C::nodesRightBackBottomCorner);
    cls.def_property_readonly("nodesRightBottomBackCorner", &C::nodesRightBottomBackCorner);
    cls.def_property_readonly("nodesBackLeftTopCorner", &C::nodesBackLeftTopCorner);
    cls.def_property_readonly("nodesTopBackLeftCorner", &C::nodesTopBackLeftCorner);
    cls.def_property_readonly("nodesTopLeftBackCorner", &C::nodesTopLeftBackCorner);
    cls.def_property_readonly("nodesLeftBackTopCorner", &C::nodesLeftBackTopCorner);
    cls.def_property_readonly("nodesLeftTopBackCorner", &C::nodesLeftTopBackCorner);
    cls.def_property_readonly("nodesBackRightTopCorner", &C::nodesBackRightTopCorner);
    cls.def_property_readonly("nodesTopBackRightCorner", &C::nodesTopBackRightCorner);
    cls.def_property_readonly("nodesTopRightBackCorner", &C::nodesTopRightBackCorner);
    cls.def_property_readonly("nodesRightBackTopCorner", &C::nodesRightBackTopCorner);
    cls.def_property_readonly("nodesRightTopBackCorner", &C::nodesRightTopBackCorner);
}

void init_Mesh(py::module& mod)
{
    py::enum_<GooseFEM::Mesh::ElementType>(
        mod, "ElementType", "See :cpp:enum:`GooseFEM::Mesh::ElementType`.")
        .value("Unknown", GooseFEM::Mesh::ElementType::Unknown)
        .value("Tri3", GooseFEM::Mesh::ElementType::Tri3)
        .value("Quad4", GooseFEM::Mesh::ElementType::Quad4)
        .value("Hex8", GooseFEM::Mesh::ElementType::Hex8)
        .export_values();

    mod.def(
        "overlapping",
        &GooseFEM::Mesh::overlapping<xt::pytensor<double, 2>, xt::pytensor<double, 2>>,
        "See :cpp:func:`GooseFEM::Mesh::overlapping`.",
        py::arg("coor_a"),
        py::arg("coor_b"),
        py::arg("rtol") = 1e-5,
        py::arg("atol") = 1e-8);

    py::class_<GooseFEM::Mesh::ManualStitch>(mod, "ManualStitch")

        .def(
            py::init<
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 1>&,
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 1>&,
                bool,
                double,
                double>(),
            "See :cpp:class:`GooseFEM::Mesh::ManualStitch`.",
            py::arg("coor_a"),
            py::arg("conn_a"),
            py::arg("overlapping_nodes_a"),
            py::arg("coor_b"),
            py::arg("conn_b"),
            py::arg("overlapping_nodes_b"),
            py::arg("check_position") = true,
            py::arg("rtol") = 1e-5,
            py::arg("atol") = 1e-8)

        .def_property_readonly("nmesh", &GooseFEM::Mesh::ManualStitch::nmesh)
        .def_property_readonly("nnode", &GooseFEM::Mesh::ManualStitch::nnode)
        .def_property_readonly("nelem", &GooseFEM::Mesh::ManualStitch::nelem)
        .def_property_readonly("ndim", &GooseFEM::Mesh::ManualStitch::ndim)
        .def_property_readonly("nne", &GooseFEM::Mesh::ManualStitch::nne)
        .def_property_readonly("coor", &GooseFEM::Mesh::ManualStitch::coor)
        .def_property_readonly("conn", &GooseFEM::Mesh::ManualStitch::conn)
        .def_property_readonly("dofs", &GooseFEM::Mesh::ManualStitch::dofs)
        .def("nodemap", py::overload_cast<>(&GooseFEM::Mesh::ManualStitch::nodemap, py::const_))
        .def("elemmap", py::overload_cast<>(&GooseFEM::Mesh::ManualStitch::elemmap, py::const_))

        .def(
            "nodemap",
            py::overload_cast<size_t>(&GooseFEM::Mesh::ManualStitch::nodemap, py::const_),
            py::arg("mesh_index"))

        .def(
            "elemmap",
            py::overload_cast<size_t>(&GooseFEM::Mesh::ManualStitch::elemmap, py::const_),
            py::arg("mesh_index"))

        .def(
            "nodeset",
            &GooseFEM::Mesh::ManualStitch::nodeset<xt::pytensor<size_t, 1>>,
            py::arg("set"),
            py::arg("mesh_index"))

        .def(
            "elemset",
            &GooseFEM::Mesh::ManualStitch::elemset<xt::pytensor<size_t, 1>>,
            py::arg("set"),
            py::arg("mesh_index"))

        .def("__repr__", [](const GooseFEM::Mesh::ManualStitch&) {
            return "<GooseFEM.Mesh.ManualStitch>";
        });

    py::class_<GooseFEM::Mesh::Stitch>(mod, "Stitch")

        .def(
            py::init<double, double>(),
            "See :cpp:class:`GooseFEM::Mesh::Stitch`.",
            py::arg("rtol") = 1e-5,
            py::arg("atol") = 1e-8)

        .def(
            "push_back",
            &GooseFEM::Mesh::Stitch::push_back<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>,
            py::arg("coor"),
            py::arg("conn"))

        .def_property_readonly("nmesh", &GooseFEM::Mesh::Stitch::nmesh)
        .def_property_readonly("nnode", &GooseFEM::Mesh::Stitch::nnode)
        .def_property_readonly("nelem", &GooseFEM::Mesh::Stitch::nelem)
        .def_property_readonly("ndim", &GooseFEM::Mesh::Stitch::ndim)
        .def_property_readonly("nne", &GooseFEM::Mesh::Stitch::nne)
        .def_property_readonly("coor", &GooseFEM::Mesh::Stitch::coor)
        .def_property_readonly("conn", &GooseFEM::Mesh::Stitch::conn)
        .def_property_readonly("dofs", &GooseFEM::Mesh::Stitch::dofs)
        .def("nodemap", py::overload_cast<>(&GooseFEM::Mesh::Stitch::nodemap, py::const_))
        .def("elemmap", py::overload_cast<>(&GooseFEM::Mesh::Stitch::elemmap, py::const_))

        .def(
            "nodemap",
            py::overload_cast<size_t>(&GooseFEM::Mesh::Stitch::nodemap, py::const_),
            py::arg("mesh_index"))

        .def(
            "elemmap",
            py::overload_cast<size_t>(&GooseFEM::Mesh::Stitch::elemmap, py::const_),
            py::arg("mesh_index"))

        .def(
            "nodeset",
            static_cast<xt::pytensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const xt::pytensor<size_t, 1>&, size_t) const>(&GooseFEM::Mesh::Stitch::nodeset),
            py::arg("set"),
            py::arg("mesh_index"))

        .def(
            "elemset",
            static_cast<xt::pytensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const xt::pytensor<size_t, 1>&, size_t) const>(&GooseFEM::Mesh::Stitch::elemset),
            py::arg("set"),
            py::arg("mesh_index"))

        .def(
            "nodeset",
            static_cast<xt::pytensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const std::vector<xt::pytensor<size_t, 1>>&) const>(
                &GooseFEM::Mesh::Stitch::nodeset),
            py::arg("set"))

        .def(
            "elemset",
            static_cast<xt::pytensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const std::vector<xt::pytensor<size_t, 1>>&) const>(
                &GooseFEM::Mesh::Stitch::elemset),
            py::arg("set"))

        .def("__repr__", [](const GooseFEM::Mesh::Stitch&) { return "<GooseFEM.Mesh.Stitch>"; });

    py::class_<GooseFEM::Mesh::Vstack, GooseFEM::Mesh::Stitch>(mod, "Vstack")

        .def(
            py::init<bool, double, double>(),
            "See :cpp:class:`GooseFEM::Mesh::Vstack`.",
            py::arg("check_overlap") = true,
            py::arg("rtol") = 1e-5,
            py::arg("atol") = 1e-8)

        .def(
            "push_back",
            &GooseFEM::Mesh::Vstack::push_back<
                xt::pytensor<double, 2>,
                xt::pytensor<size_t, 2>,
                xt::pytensor<size_t, 1>>,
            py::arg("coor"),
            py::arg("conn"),
            py::arg("nodes_bottom"),
            py::arg("nodes_top"))

        .def("__repr__", [](const GooseFEM::Mesh::Vstack&) { return "<GooseFEM.Mesh.Vstack>"; });

    py::class_<GooseFEM::Mesh::Renumber>(mod, "Renumber")

        .def(
            py::init<const xt::pyarray<size_t>&>(),
            "See :cpp:class:`GooseFEM::Mesh::Renumber`.",
            py::arg("dofs"))

        .def("apply", &GooseFEM::Mesh::Renumber::apply<xt::pyarray<size_t>>)
        .def_property_readonly("index", &GooseFEM::Mesh::Renumber::index)

        .def(
            "__repr__", [](const GooseFEM::Mesh::Renumber&) { return "<GooseFEM.Mesh.Renumber>"; });

    py::class_<GooseFEM::Mesh::Reorder>(mod, "Reorder")

        .def(
            py::init([](xt::pytensor<size_t, 1>& a) { return new GooseFEM::Mesh::Reorder({a}); }),
            "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(
            py::init([](xt::pytensor<size_t, 1>& a, xt::pytensor<size_t, 1>& b) {
                return new GooseFEM::Mesh::Reorder({a, b});
            }),
            "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(
            py::init([](xt::pytensor<size_t, 1>& a,
                        xt::pytensor<size_t, 1>& b,
                        xt::pytensor<size_t, 1>& c) {
                return new GooseFEM::Mesh::Reorder({a, b, c});
            }),
            "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(
            py::init([](xt::pytensor<size_t, 1>& a,
                        xt::pytensor<size_t, 1>& b,
                        xt::pytensor<size_t, 1>& c,
                        xt::pytensor<size_t, 1>& d) {
                return new GooseFEM::Mesh::Reorder({a, b, c, d});
            }),
            "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def("apply", &GooseFEM::Mesh::Reorder::apply<xt::pyarray<size_t>>)
        .def_property_readonly("index", &GooseFEM::Mesh::Reorder::index)

        .def("__repr__", [](const GooseFEM::Mesh::Reorder&) { return "<GooseFEM.Mesh.Reorder>"; });

    mod.def(
        "dofs",
        &GooseFEM::Mesh::dofs,
        "See :cpp:func:`GooseFEM::Mesh::dofs`.",
        py::arg("nnode"),
        py::arg("ndim"));

    mod.def(
        "nodaltyings",
        &GooseFEM::Mesh::nodaltyings<xt::pytensor<size_t, 2>>,
        "See :cpp:func:`GooseFEM::Mesh::nodaltyings`.",
        py::arg("dofs"));

    mod.def(
        "renumber",
        &GooseFEM::Mesh::renumber<xt::pytensor<size_t, 2>>,
        "See :cpp:func:`GooseFEM::Mesh::renumber`.",
        py::arg("dofs"));

    mod.def(
        "coordination",
        &GooseFEM::Mesh::coordination<xt::pytensor<size_t, 2>>,
        "See :cpp:func:`GooseFEM::Mesh::coordination`.",
        py::arg("conn"));

    mod.def(
        "node2dof",
        &GooseFEM::Mesh::node2dof<xt::pytensor<size_t, 2>>,
        "See :cpp:func:`GooseFEM::Mesh::node2dof`.",
        py::arg("dofs"),
        py::arg("sorted") = true);

    mod.def(
        "elem2node",
        py::overload_cast<const xt::pytensor<size_t, 2>&, bool>(
            &GooseFEM::Mesh::elem2node<xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::elem2node`.",
        py::arg("conn"),
        py::arg("sorted") = true);

    mod.def(
        "elem2node",
        py::overload_cast<const xt::pytensor<size_t, 2>&, const xt::pytensor<size_t, 2>&, bool>(
            &GooseFEM::Mesh::elem2node<xt::pytensor<size_t, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::elem2node`.",
        py::arg("conn"),
        py::arg("dofs"),
        py::arg("sorted") = true);

    mod.def(
        "edgesize",
        py::overload_cast<const xt::pytensor<double, 2>&, const xt::pytensor<size_t, 2>&>(
            &GooseFEM::Mesh::edgesize<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::edgesize`.",
        py::arg("coor"),
        py::arg("conn"));

    mod.def(
        "edgesize",
        py::overload_cast<
            const xt::pytensor<double, 2>&,
            const xt::pytensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(
            &GooseFEM::Mesh::edgesize<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::edgesize`.",
        py::arg("coor"),
        py::arg("conn"),
        py::arg("type"));

    mod.def(
        "centers",
        py::overload_cast<const xt::pytensor<double, 2>&, const xt::pytensor<size_t, 2>&>(
            &GooseFEM::Mesh::centers<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::centers`.",
        py::arg("coor"),
        py::arg("conn"));

    mod.def(
        "centers",
        py::overload_cast<
            const xt::pytensor<double, 2>&,
            const xt::pytensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(
            &GooseFEM::Mesh::centers<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::centers`.",
        py::arg("coor"),
        py::arg("conn"),
        py::arg("type"));

    mod.def(
        "elemmap2nodemap",
        py::overload_cast<
            const xt::pytensor<size_t, 1>&,
            const xt::pytensor<double, 2>&,
            const xt::pytensor<size_t, 2>&>(&GooseFEM::Mesh::elemmap2nodemap<
                                            xt::pytensor<size_t, 1>,
                                            xt::pytensor<double, 2>,
                                            xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::elemmap2nodemap`.",
        py::arg("elem_map"),
        py::arg("coor"),
        py::arg("conn"));

    mod.def(
        "elemmap2nodemap",
        py::overload_cast<
            const xt::pytensor<size_t, 1>&,
            const xt::pytensor<double, 2>&,
            const xt::pytensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::elemmap2nodemap<
                                         xt::pytensor<size_t, 1>,
                                         xt::pytensor<double, 2>,
                                         xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::elemmap2nodemap`.",
        py::arg("elem_map"),
        py::arg("coor"),
        py::arg("conn"),
        py::arg("type"));

    mod.def(
        "nodal_mass",
        py::overload_cast<const xt::pytensor<double, 2>&, const xt::pytensor<size_t, 2>&>(
            &GooseFEM::Mesh::nodal_mass<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::nodal_mass`.",
        py::arg("coor"),
        py::arg("conn"));

    mod.def(
        "nodal_mass",
        py::overload_cast<
            const xt::pytensor<double, 2>&,
            const xt::pytensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(
            &GooseFEM::Mesh::nodal_mass<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::nodal_mass`.",
        py::arg("coor"),
        py::arg("conn"),
        py::arg("type"));

    mod.def(
        "center_of_gravity",
        py::overload_cast<const xt::pytensor<double, 2>&, const xt::pytensor<size_t, 2>&>(
            &GooseFEM::Mesh::center_of_gravity<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::center_of_gravity`.",
        py::arg("coor"),
        py::arg("conn"));

    mod.def(
        "center_of_gravity",
        py::overload_cast<
            const xt::pytensor<double, 2>&,
            const xt::pytensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(
            &GooseFEM::Mesh::center_of_gravity<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
        "See :cpp:func:`GooseFEM::Mesh::center_of_gravity`.",
        py::arg("coor"),
        py::arg("conn"),
        py::arg("type"));
}

#endif

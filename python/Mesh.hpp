/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
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
void register_Element_RegularBase(P& cls)
{
    cls.def("nelem", &C::nelem);
    cls.def("nnode", &C::nnode);
    cls.def("nne", &C::nne);
    cls.def("ndim", &C::ndim);
    cls.def("nelx", &C::nelx);
    cls.def("nely", &C::nely);
    cls.def("h", &C::h);
    cls.def("getElementType", &C::getElementType);
    cls.def("coor", &C::coor);
    cls.def("conn", &C::conn);
    cls.def("dofs", &C::dofs);
    cls.def("dofsPeriodic", &C::dofsPeriodic);
    cls.def("nodesPeriodic", &C::nodesPeriodic);
    cls.def("nodesOrigin", &C::nodesOrigin);
}

template <class C, class P>
void register_Element_RegularBase2d(P& cls)
{
     cls.def("nodesBottomEdge", &C::nodesBottomEdge);
     cls.def("nodesTopEdge", &C::nodesTopEdge);
     cls.def("nodesLeftEdge", &C::nodesLeftEdge);
     cls.def("nodesRightEdge", &C::nodesRightEdge);
     cls.def("nodesBottomOpenEdge", &C::nodesBottomOpenEdge);
     cls.def("nodesTopOpenEdge", &C::nodesTopOpenEdge);
     cls.def("nodesLeftOpenEdge", &C::nodesLeftOpenEdge);
     cls.def("nodesRightOpenEdge", &C::nodesRightOpenEdge);
     cls.def("nodesBottomLeftCorner", &C::nodesBottomLeftCorner);
     cls.def("nodesBottomRightCorner", &C::nodesBottomRightCorner);
     cls.def("nodesTopLeftCorner", &C::nodesTopLeftCorner);
     cls.def("nodesTopRightCorner", &C::nodesTopRightCorner);
     cls.def("nodesLeftBottomCorner", &C::nodesLeftBottomCorner);
     cls.def("nodesLeftTopCorner", &C::nodesLeftTopCorner);
     cls.def("nodesRightBottomCorner", &C::nodesRightBottomCorner);
     cls.def("nodesRightTopCorner", &C::nodesRightTopCorner);
}

template <class C, class P>
void register_Element_RegularBase3d(P& cls)
{
     cls.def("nodesFront", &C::nodesFront);
     cls.def("nodesBack", &C::nodesBack);
     cls.def("nodesLeft", &C::nodesLeft);
     cls.def("nodesRight", &C::nodesRight);
     cls.def("nodesBottom", &C::nodesBottom);
     cls.def("nodesTop", &C::nodesTop);
     cls.def("nodesFrontFace", &C::nodesFrontFace);
     cls.def("nodesBackFace", &C::nodesBackFace);
     cls.def("nodesLeftFace", &C::nodesLeftFace);
     cls.def("nodesRightFace", &C::nodesRightFace);
     cls.def("nodesBottomFace", &C::nodesBottomFace);
     cls.def("nodesTopFace", &C::nodesTopFace);
     cls.def("nodesFrontBottomEdge", &C::nodesFrontBottomEdge);
     cls.def("nodesFrontTopEdge", &C::nodesFrontTopEdge);
     cls.def("nodesFrontLeftEdge", &C::nodesFrontLeftEdge);
     cls.def("nodesFrontRightEdge", &C::nodesFrontRightEdge);
     cls.def("nodesBackBottomEdge", &C::nodesBackBottomEdge);
     cls.def("nodesBackTopEdge", &C::nodesBackTopEdge);
     cls.def("nodesBackLeftEdge", &C::nodesBackLeftEdge);
     cls.def("nodesBackRightEdge", &C::nodesBackRightEdge);
     cls.def("nodesBottomLeftEdge", &C::nodesBottomLeftEdge);
     cls.def("nodesBottomRightEdge", &C::nodesBottomRightEdge);
     cls.def("nodesTopLeftEdge", &C::nodesTopLeftEdge);
     cls.def("nodesTopRightEdge", &C::nodesTopRightEdge);
     cls.def("nodesBottomFrontEdge", &C::nodesBottomFrontEdge);
     cls.def("nodesBottomBackEdge", &C::nodesBottomBackEdge);
     cls.def("nodesTopFrontEdge", &C::nodesTopFrontEdge);
     cls.def("nodesTopBackEdge", &C::nodesTopBackEdge);
     cls.def("nodesLeftBottomEdge", &C::nodesLeftBottomEdge);
     cls.def("nodesLeftFrontEdge", &C::nodesLeftFrontEdge);
     cls.def("nodesLeftBackEdge", &C::nodesLeftBackEdge);
     cls.def("nodesLeftTopEdge", &C::nodesLeftTopEdge);
     cls.def("nodesRightBottomEdge", &C::nodesRightBottomEdge);
     cls.def("nodesRightTopEdge", &C::nodesRightTopEdge);
     cls.def("nodesRightFrontEdge", &C::nodesRightFrontEdge);
     cls.def("nodesRightBackEdge", &C::nodesRightBackEdge);
     cls.def("nodesFrontBottomOpenEdge", &C::nodesFrontBottomOpenEdge);
     cls.def("nodesFrontTopOpenEdge", &C::nodesFrontTopOpenEdge);
     cls.def("nodesFrontLeftOpenEdge", &C::nodesFrontLeftOpenEdge);
     cls.def("nodesFrontRightOpenEdge", &C::nodesFrontRightOpenEdge);
     cls.def("nodesBackBottomOpenEdge", &C::nodesBackBottomOpenEdge);
     cls.def("nodesBackTopOpenEdge", &C::nodesBackTopOpenEdge);
     cls.def("nodesBackLeftOpenEdge", &C::nodesBackLeftOpenEdge);
     cls.def("nodesBackRightOpenEdge", &C::nodesBackRightOpenEdge);
     cls.def("nodesBottomLeftOpenEdge", &C::nodesBottomLeftOpenEdge);
     cls.def("nodesBottomRightOpenEdge", &C::nodesBottomRightOpenEdge);
     cls.def("nodesTopLeftOpenEdge", &C::nodesTopLeftOpenEdge);
     cls.def("nodesTopRightOpenEdge", &C::nodesTopRightOpenEdge);
     cls.def("nodesBottomFrontOpenEdge", &C::nodesBottomFrontOpenEdge);
     cls.def("nodesBottomBackOpenEdge", &C::nodesBottomBackOpenEdge);
     cls.def("nodesTopFrontOpenEdge", &C::nodesTopFrontOpenEdge);
     cls.def("nodesTopBackOpenEdge", &C::nodesTopBackOpenEdge);
     cls.def("nodesLeftBottomOpenEdge", &C::nodesLeftBottomOpenEdge);
     cls.def("nodesLeftFrontOpenEdge", &C::nodesLeftFrontOpenEdge);
     cls.def("nodesLeftBackOpenEdge", &C::nodesLeftBackOpenEdge);
     cls.def("nodesLeftTopOpenEdge", &C::nodesLeftTopOpenEdge);
     cls.def("nodesRightBottomOpenEdge", &C::nodesRightBottomOpenEdge);
     cls.def("nodesRightTopOpenEdge", &C::nodesRightTopOpenEdge);
     cls.def("nodesRightFrontOpenEdge", &C::nodesRightFrontOpenEdge);
     cls.def("nodesRightBackOpenEdge", &C::nodesRightBackOpenEdge);
     cls.def("nodesFrontBottomLeftCorner", &C::nodesFrontBottomLeftCorner);
     cls.def("nodesFrontBottomRightCorner", &C::nodesFrontBottomRightCorner);
     cls.def("nodesFrontTopLeftCorner", &C::nodesFrontTopLeftCorner);
     cls.def("nodesFrontTopRightCorner", &C::nodesFrontTopRightCorner);
     cls.def("nodesBackBottomLeftCorner", &C::nodesBackBottomLeftCorner);
     cls.def("nodesBackBottomRightCorner", &C::nodesBackBottomRightCorner);
     cls.def("nodesBackTopLeftCorner", &C::nodesBackTopLeftCorner);
     cls.def("nodesBackTopRightCorner", &C::nodesBackTopRightCorner);
     cls.def("nodesFrontLeftBottomCorner", &C::nodesFrontLeftBottomCorner);
     cls.def("nodesBottomFrontLeftCorner", &C::nodesBottomFrontLeftCorner);
     cls.def("nodesBottomLeftFrontCorner", &C::nodesBottomLeftFrontCorner);
     cls.def("nodesLeftFrontBottomCorner", &C::nodesLeftFrontBottomCorner);
     cls.def("nodesLeftBottomFrontCorner", &C::nodesLeftBottomFrontCorner);
     cls.def("nodesFrontRightBottomCorner", &C::nodesFrontRightBottomCorner);
     cls.def("nodesBottomFrontRightCorner", &C::nodesBottomFrontRightCorner);
     cls.def("nodesBottomRightFrontCorner", &C::nodesBottomRightFrontCorner);
     cls.def("nodesRightFrontBottomCorner", &C::nodesRightFrontBottomCorner);
     cls.def("nodesRightBottomFrontCorner", &C::nodesRightBottomFrontCorner);
     cls.def("nodesFrontLeftTopCorner", &C::nodesFrontLeftTopCorner);
     cls.def("nodesTopFrontLeftCorner", &C::nodesTopFrontLeftCorner);
     cls.def("nodesTopLeftFrontCorner", &C::nodesTopLeftFrontCorner);
     cls.def("nodesLeftFrontTopCorner", &C::nodesLeftFrontTopCorner);
     cls.def("nodesLeftTopFrontCorner", &C::nodesLeftTopFrontCorner);
     cls.def("nodesFrontRightTopCorner", &C::nodesFrontRightTopCorner);
     cls.def("nodesTopFrontRightCorner", &C::nodesTopFrontRightCorner);
     cls.def("nodesTopRightFrontCorner", &C::nodesTopRightFrontCorner);
     cls.def("nodesRightFrontTopCorner", &C::nodesRightFrontTopCorner);
     cls.def("nodesRightTopFrontCorner", &C::nodesRightTopFrontCorner);
     cls.def("nodesBackLeftBottomCorner", &C::nodesBackLeftBottomCorner);
     cls.def("nodesBottomBackLeftCorner", &C::nodesBottomBackLeftCorner);
     cls.def("nodesBottomLeftBackCorner", &C::nodesBottomLeftBackCorner);
     cls.def("nodesLeftBackBottomCorner", &C::nodesLeftBackBottomCorner);
     cls.def("nodesLeftBottomBackCorner", &C::nodesLeftBottomBackCorner);
     cls.def("nodesBackRightBottomCorner", &C::nodesBackRightBottomCorner);
     cls.def("nodesBottomBackRightCorner", &C::nodesBottomBackRightCorner);
     cls.def("nodesBottomRightBackCorner", &C::nodesBottomRightBackCorner);
     cls.def("nodesRightBackBottomCorner", &C::nodesRightBackBottomCorner);
     cls.def("nodesRightBottomBackCorner", &C::nodesRightBottomBackCorner);
     cls.def("nodesBackLeftTopCorner", &C::nodesBackLeftTopCorner);
     cls.def("nodesTopBackLeftCorner", &C::nodesTopBackLeftCorner);
     cls.def("nodesTopLeftBackCorner", &C::nodesTopLeftBackCorner);
     cls.def("nodesLeftBackTopCorner", &C::nodesLeftBackTopCorner);
     cls.def("nodesLeftTopBackCorner", &C::nodesLeftTopBackCorner);
     cls.def("nodesBackRightTopCorner", &C::nodesBackRightTopCorner);
     cls.def("nodesTopBackRightCorner", &C::nodesTopBackRightCorner);
     cls.def("nodesTopRightBackCorner", &C::nodesTopRightBackCorner);
     cls.def("nodesRightBackTopCorner", &C::nodesRightBackTopCorner);
     cls.def("nodesRightTopBackCorner", &C::nodesRightTopBackCorner);
}

void init_Mesh(py::module& mod)
{
    py::enum_<GooseFEM::Mesh::ElementType>(mod,
            "ElementType",
            "ElementType."
            "See :cpp:enum:`GooseFEM::Mesh::ElementType`.")
        .value("Unknown", GooseFEM::Mesh::ElementType::Unknown)
        .value("Tri3", GooseFEM::Mesh::ElementType::Tri3)
        .value("Quad4", GooseFEM::Mesh::ElementType::Quad4)
        .value("Hex8", GooseFEM::Mesh::ElementType::Hex8)
        .export_values();

    mod.def("overlapping",
            &GooseFEM::Mesh::overlapping<xt::pytensor<double, 2>, xt::pytensor<double, 2>>,
            "Find overlapping nodes."
            "See :cpp:func:`GooseFEM::Mesh::overlapping`.",
            py::arg("coor_a"),
            py::arg("coor_b"),
            py::arg("rtol") = 1e-5,
            py::arg("atol") = 1e-8);

    py::class_<GooseFEM::Mesh::ManualStitch>(mod, "ManualStitch")

        .def(py::init<
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 1>&,
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 1>&,
                bool,
                double,
                double>(),
             "Manually stitch meshes."
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

        .def("nmesh",
             &GooseFEM::Mesh::ManualStitch::nmesh,
             "Number of sub meshes."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::nmesh`.")

        .def("nnode",
             &GooseFEM::Mesh::ManualStitch::nnode,
             "Number of nodes of stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::nnode`.")

        .def("nelem",
             &GooseFEM::Mesh::ManualStitch::nelem,
             "Number of elements of stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::nelem`.")

        .def("ndim",
             &GooseFEM::Mesh::ManualStitch::ndim,
             "Number of dimensions (of stitched mesh)."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::ndim`.")

        .def("nne",
             &GooseFEM::Mesh::ManualStitch::nne,
             "Number of nodes per element (of stitched mesh)."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::nne`.")

        .def("coor",
             &GooseFEM::Mesh::ManualStitch::coor,
             "Coordinates of stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::coor`.")

        .def("conn",
             &GooseFEM::Mesh::ManualStitch::conn,
             "Connectivity of stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::conn`.")

        .def("dofs",
             &GooseFEM::Mesh::ManualStitch::dofs,
             "DOF numbers per node."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::dofs`.")

        .def("nodemap",
             py::overload_cast<>(
                &GooseFEM::Mesh::ManualStitch::nodemap, py::const_),
             "Node-map for givall sub-meshes."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::nodemap`.")

        .def("elemmap",
             py::overload_cast<>(
                &GooseFEM::Mesh::ManualStitch::elemmap, py::const_),
             "Element-map for all sub-meshes."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::elemmap`.")

        .def("nodemap",
             py::overload_cast<size_t>(
                &GooseFEM::Mesh::ManualStitch::nodemap, py::const_),
             "Node-map for given sub-mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::nodemap`.",
             py::arg("mesh_index"))

        .def("elemmap",
             py::overload_cast<size_t>(
                &GooseFEM::Mesh::ManualStitch::elemmap, py::const_),
             "Element-map for given sub-mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::elemmap`.",
             py::arg("mesh_index"))

        .def("nodeset",
             &GooseFEM::Mesh::ManualStitch::nodeset<xt::pytensor<size_t, 1>>,
             "Convert node-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::nodeset`.",
             py::arg("set"),
             py::arg("mesh_index"))

        .def("elemset",
             &GooseFEM::Mesh::ManualStitch::elemset<xt::pytensor<size_t, 1>>,
             "Convert element-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::elemset`.",
             py::arg("set"),
             py::arg("mesh_index"))

        .def("__repr__",
            [](const GooseFEM::Mesh::ManualStitch&) { return "<GooseFEM.Mesh.ManualStitch>"; });

    py::class_<GooseFEM::Mesh::Stitch>(mod, "Stitch")

        .def(py::init<double, double>(),
             "Stitch meshes."
             "See :cpp:class:`GooseFEM::Mesh::Stitch`.",
             py::arg("rtol") = 1e-5,
             py::arg("atol") = 1e-8)

        .def("push_back",
             &GooseFEM::Mesh::Stitch::push_back<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>,
             "Add mesh to be stitched."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::push_back`.",
             py::arg("coor"),
             py::arg("conn"))

        .def("nmesh",
             &GooseFEM::Mesh::Stitch::nmesh,
             "Number of sub meshes."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nmesh`.")

        .def("nnode",
             &GooseFEM::Mesh::Stitch::nnode,
             "Number of nodes of stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nnode`.")

        .def("nelem",
             &GooseFEM::Mesh::Stitch::nelem,
             "Number of elements of stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nelem`.")

        .def("ndim",
             &GooseFEM::Mesh::Stitch::ndim,
             "Number of dimensions (of stitched mesh)."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::ndim`.")

        .def("nne",
             &GooseFEM::Mesh::Stitch::nne,
             "Number of nodes per element (of stitched mesh)."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nne`.")

        .def("coor",
             &GooseFEM::Mesh::Stitch::coor,
             "Coordinates of stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::coor`.")

        .def("conn",
             &GooseFEM::Mesh::Stitch::conn,
             "Connectivity of stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::conn`.")

        .def("dofs",
             &GooseFEM::Mesh::Stitch::dofs,
             "DOF numbers per node."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::dofs`.")

        .def("nodemap",
             py::overload_cast<>(
                &GooseFEM::Mesh::Stitch::nodemap, py::const_),
             "Node-map for givall sub-meshes."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nodemap`.")

        .def("elemmap",
             py::overload_cast<>(
                &GooseFEM::Mesh::Stitch::elemmap, py::const_),
             "Element-map for all sub-meshes."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::elemmap`.")

        .def("nodemap",
             py::overload_cast<size_t>(
                &GooseFEM::Mesh::Stitch::nodemap, py::const_),
             "Node-map for given sub-mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nodemap`.",
             py::arg("mesh_index"))

        .def("elemmap",
             py::overload_cast<size_t>(
                &GooseFEM::Mesh::Stitch::elemmap, py::const_),
             "Element-map for given sub-mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::elemmap`.",
             py::arg("mesh_index"))

        .def("nodeset",
             static_cast<xt::pytensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const xt::pytensor<size_t, 1>&, size_t) const>(&GooseFEM::Mesh::Stitch::nodeset),
             "Convert node-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nodeset`.",
             py::arg("set"),
             py::arg("mesh_index"))

        .def("elemset",
             static_cast<xt::pytensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const xt::pytensor<size_t, 1>&, size_t) const>(&GooseFEM::Mesh::Stitch::elemset),
             "Convert element-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::elemset`.",
             py::arg("set"),
             py::arg("mesh_index"))

        .def("nodeset",
             static_cast<xt::pytensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const std::vector<xt::pytensor<size_t, 1>>&) const>(&GooseFEM::Mesh::Stitch::nodeset),
             "Convert node-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nodeset`.",
             py::arg("set"))

        .def("elemset",
             static_cast<xt::pytensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const std::vector<xt::pytensor<size_t, 1>>&) const>(&GooseFEM::Mesh::Stitch::elemset),
             "Convert element-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::elemset`.",
             py::arg("set"))

        .def("__repr__", [](const GooseFEM::Mesh::Stitch&) { return "<GooseFEM.Mesh.Stitch>"; });

    py::class_<GooseFEM::Mesh::Vstack, GooseFEM::Mesh::Stitch>(mod, "Vstack")

        .def(py::init<bool, double, double>(),
             "Vstack meshes."
             "See :cpp:class:`GooseFEM::Mesh::Vstack`.",
             py::arg("check_overlap") = true,
             py::arg("rtol") = 1e-5,
             py::arg("atol") = 1e-8)

        .def("push_back",
             &GooseFEM::Mesh::Vstack::push_back<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>, xt::pytensor<size_t, 1>>,
             "Add mesh to be stitched."
             "See :cpp:func:`GooseFEM::Mesh::Vstack::push_back`.",
             py::arg("coor"),
             py::arg("conn"),
             py::arg("nodes_bottom"),
             py::arg("nodes_top"))

        .def("__repr__", [](const GooseFEM::Mesh::Vstack&) { return "<GooseFEM.Mesh.Vstack>"; });

    py::class_<GooseFEM::Mesh::Renumber>(mod, "Renumber")

        .def(py::init<const xt::pyarray<size_t>&>(),
             "Renumber to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Renumber`.",
             py::arg("dofs"))

        .def("get",
             &GooseFEM::Mesh::Renumber::apply<xt::pyarray<size_t>>,
             "Get renumbered DOFs."
             "See :cpp:func:`GooseFEM::Mesh::Renumber::get`.")

        .def("apply",
             &GooseFEM::Mesh::Renumber::apply<xt::pyarray<size_t>>,
             "Get renumbered list."
             "See :cpp:func:`GooseFEM::Mesh::Renumber::apply`.")

        .def("index",
             &GooseFEM::Mesh::Renumber::index,
             "Get index list to apply renumbering. Apply renumbering using ``index[dofs]``."
             "See :cpp:func:`GooseFEM::Mesh::Renumber::index`.")

        .def("__repr__",
             [](const GooseFEM::Mesh::Renumber&) { return "<GooseFEM.Mesh.Renumber>"; });

    py::class_<GooseFEM::Mesh::Reorder>(mod, "Reorder")

        .def(py::init([](xt::pytensor<size_t, 1>& a) { return new GooseFEM::Mesh::Reorder({a}); }),
             "Reorder to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(py::init([](xt::pytensor<size_t, 1>& a, xt::pytensor<size_t, 1>& b) {
                return new GooseFEM::Mesh::Reorder({a, b});
             }),
             "Reorder to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(py::init(
                [](xt::pytensor<size_t, 1>& a, xt::pytensor<size_t, 1>& b, xt::pytensor<size_t, 1>& c) {
                    return new GooseFEM::Mesh::Reorder({a, b, c});
                }),
             "Reorder to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(py::init([](xt::pytensor<size_t, 1>& a,
                         xt::pytensor<size_t, 1>& b,
                         xt::pytensor<size_t, 1>& c,
                         xt::pytensor<size_t, 1>& d) {
                return new GooseFEM::Mesh::Reorder({a, b, c, d});
                }),
             "Reorder to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def("get",
             &GooseFEM::Mesh::Reorder::apply<xt::pyarray<size_t>>,
             "Reorder matrix (e.g. ``dofs``)."
             "See :cpp:func:`GooseFEM::Mesh::Reorder::get`.",
             py::arg("dofs"))

        .def("apply",
             &GooseFEM::Mesh::Reorder::apply<xt::pyarray<size_t>>,
             "Get reordered list."
             "See :cpp:func:`GooseFEM::Mesh::Reorder::apply`.")

        .def("index",
             &GooseFEM::Mesh::Reorder::index,
             "Get index list to apply renumbering. Apply renumbering using ``index[dofs]``."
             "See :cpp:func:`GooseFEM::Mesh::Reorder::index`.")

        .def("__repr__", [](const GooseFEM::Mesh::Reorder&) { return "<GooseFEM.Mesh.Reorder>"; });

    mod.def("dofs",
            &GooseFEM::Mesh::dofs,
            "List with DOF-numbers in sequential order."
            "See :cpp:func:`GooseFEM::Mesh::dofs`.",
            py::arg("nnode"),
            py::arg("ndim"));

    mod.def("renumber",
            &GooseFEM::Mesh::renumber<xt::pytensor<size_t, 2>>,
            "Renumber to lowest possible indices."
            "See :cpp:func:`GooseFEM::Mesh::renumber`.",
            py::arg("dofs"));

    mod.def("coordination",
            &GooseFEM::Mesh::coordination<xt::pytensor<size_t, 2>>,
            "Coordination number of each node."
            "See :cpp:func:`GooseFEM::Mesh::coordination`.",
            py::arg("conn"));

    mod.def("elem2node",
            &GooseFEM::Mesh::elem2node<xt::pytensor<size_t, 2>>,
            "Element-numbers connected to each node."
            "See :cpp:func:`GooseFEM::Mesh::elem2node`.",
            py::arg("conn"),
            py::arg("sorted") = true);

    mod.def("edgesize",
            py::overload_cast<const xt::pytensor<double, 2>&, const xt::pytensor<size_t, 2>&>(
                &GooseFEM::Mesh::edgesize<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
            "Get the edge size of all elements."
            "See :cpp:func:`GooseFEM::Mesh::edgesize`.",
            py::arg("coor"),
            py::arg("conn"));

    mod.def("edgesize",
            py::overload_cast<
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::edgesize<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
            "Get the edge size of all elements."
            "See :cpp:func:`GooseFEM::Mesh::edgesize`.",
            py::arg("coor"),
            py::arg("conn"),
            py::arg("type"));

    mod.def("centers",
            py::overload_cast<const xt::pytensor<double, 2>&, const xt::pytensor<size_t, 2>&>(
                &GooseFEM::Mesh::centers<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
            "Coordinates of the center of each element."
            "See :cpp:func:`GooseFEM::Mesh::centers`.",
            py::arg("coor"),
            py::arg("conn"));

    mod.def("centers",
            py::overload_cast<
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::centers<xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
            "Coordinates of the center of each element."
            "See :cpp:func:`GooseFEM::Mesh::centers`.",
            py::arg("coor"),
            py::arg("conn"),
            py::arg("type"));

    mod.def("elemmap2nodemap",
            py::overload_cast<
                const xt::pytensor<size_t, 1>&,
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&>(&GooseFEM::Mesh::elemmap2nodemap<
                    xt::pytensor<size_t, 1>, xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
            "Convert an element-map to a node-map."
            "See :cpp:func:`GooseFEM::Mesh::elemmap2nodemap`.",
            py::arg("elem_map"),
            py::arg("coor"),
            py::arg("conn"));

    mod.def("elemmap2nodemap",
            py::overload_cast<
                const xt::pytensor<size_t, 1>&,
                const xt::pytensor<double, 2>&,
                const xt::pytensor<size_t, 2>&,
                GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::elemmap2nodemap<
                    xt::pytensor<size_t, 1>, xt::pytensor<double, 2>, xt::pytensor<size_t, 2>>),
            "Convert an element-map to a node-map."
            "See :cpp:func:`GooseFEM::Mesh::elemmap2nodemap`.",
            py::arg("elem_map"),
            py::arg("coor"),
            py::arg("conn"),
            py::arg("type"));
}

#endif

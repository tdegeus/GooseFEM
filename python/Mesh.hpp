/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_Mesh(py::module& m)
{

    py::enum_<GooseFEM::Mesh::ElementType>(m,
            "ElementType",
            "ElementType."
            "See :cpp:enum:`GooseFEM::Mesh::ElementType`.")
        .value("Unknown", GooseFEM::Mesh::ElementType::Unknown)
        .value("Tri3", GooseFEM::Mesh::ElementType::Tri3)
        .value("Quad4", GooseFEM::Mesh::ElementType::Quad4)
        .value("Hex8", GooseFEM::Mesh::ElementType::Hex8)
        .export_values();

    m.def("overlapping",
          &GooseFEM::Mesh::overlapping<xt::xtensor<double, 2>, xt::xtensor<double, 2>>,
          "Find overlapping nodes."
          "See :cpp:func:`GooseFEM::Mesh::overlapping`.",
          py::arg("coor_a"),
          py::arg("coor_b"),
          py::arg("rtol") = 1e-5,
          py::arg("atol") = 1e-8);

    py::class_<GooseFEM::Mesh::ManualStitch>(m, "ManualStitch")

        .def(py::init<
                const xt::xtensor<double, 2>&,
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 1>&,
                const xt::xtensor<double, 2>&,
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 1>&,
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
             &GooseFEM::Mesh::ManualStitch::nodeset<xt::xtensor<size_t, 1>>,
             "Convert node-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::nodeset`.",
             py::arg("set"),
             py::arg("mesh_index"))

        .def("elemset",
             &GooseFEM::Mesh::ManualStitch::elemset<xt::xtensor<size_t, 1>>,
             "Convert element-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::ManualStitch::elemset`.",
             py::arg("set"),
             py::arg("mesh_index"))

        .def("__repr__",
            [](const GooseFEM::Mesh::ManualStitch&) { return "<GooseFEM.Mesh.ManualStitch>"; });

    py::class_<GooseFEM::Mesh::Stitch>(m, "Stitch")

        .def(py::init<double, double>(),
             "Stitch meshes."
             "See :cpp:class:`GooseFEM::Mesh::Stitch`.",
             py::arg("rtol") = 1e-5,
             py::arg("atol") = 1e-8)

        .def("push_back",
             &GooseFEM::Mesh::Stitch::push_back<xt::xtensor<double, 2>, xt::xtensor<size_t, 2>>,
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
             static_cast<xt::xtensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const xt::xtensor<size_t, 1>&, size_t) const>(&GooseFEM::Mesh::Stitch::nodeset),
             "Convert node-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nodeset`.",
             py::arg("set"),
             py::arg("mesh_index"))

        .def("elemset",
             static_cast<xt::xtensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const xt::xtensor<size_t, 1>&, size_t) const>(&GooseFEM::Mesh::Stitch::elemset),
             "Convert element-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::elemset`.",
             py::arg("set"),
             py::arg("mesh_index"))

        .def("nodeset",
             static_cast<xt::xtensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const std::vector<xt::xtensor<size_t, 1>>&) const>(&GooseFEM::Mesh::Stitch::nodeset),
             "Convert node-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::nodeset`.",
             py::arg("set"))

        .def("elemset",
             static_cast<xt::xtensor<size_t, 1> (GooseFEM::Mesh::Stitch::*)(
                const std::vector<xt::xtensor<size_t, 1>>&) const>(&GooseFEM::Mesh::Stitch::elemset),
             "Convert element-set to the stitched mesh."
             "See :cpp:func:`GooseFEM::Mesh::Stitch::elemset`.",
             py::arg("set"))

        .def("__repr__", [](const GooseFEM::Mesh::Stitch&) { return "<GooseFEM.Mesh.Stitch>"; });

    py::class_<GooseFEM::Mesh::Vstack, GooseFEM::Mesh::Stitch>(m, "Vstack")

        .def(py::init<bool, double, double>(),
             "Vstack meshes."
             "See :cpp:class:`GooseFEM::Mesh::Vstack`.",
             py::arg("check_overlap") = true,
             py::arg("rtol") = 1e-5,
             py::arg("atol") = 1e-8)

        .def("push_back",
             &GooseFEM::Mesh::Vstack::push_back<xt::xtensor<double, 2>, xt::xtensor<size_t, 2>, xt::xtensor<size_t, 1>>,
             "Add mesh to be stitched."
             "See :cpp:func:`GooseFEM::Mesh::Vstack::push_back`.",
             py::arg("coor"),
             py::arg("conn"),
             py::arg("nodes_bottom"),
             py::arg("nodes_top"))

        .def("__repr__", [](const GooseFEM::Mesh::Vstack&) { return "<GooseFEM.Mesh.Vstack>"; });

    py::class_<GooseFEM::Mesh::Renumber>(m, "Renumber")

        .def(py::init<const xt::xarray<size_t>&>(),
             "Renumber to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Renumber`.",
             py::arg("dofs"))

        .def("get",
             &GooseFEM::Mesh::Renumber::apply<xt::xarray<size_t>>,
             "Get renumbered DOFs."
             "See :cpp:func:`GooseFEM::Mesh::Renumber::get`.")

        .def("apply",
             &GooseFEM::Mesh::Renumber::apply<xt::xarray<size_t>>,
             "Get renumbered list."
             "See :cpp:func:`GooseFEM::Mesh::Renumber::apply`.")

        .def("index",
             &GooseFEM::Mesh::Renumber::index,
             "Get index list to apply renumbering. Apply renumbering using ``index[dofs]``."
             "See :cpp:func:`GooseFEM::Mesh::Renumber::index`.")

        .def("__repr__",
             [](const GooseFEM::Mesh::Renumber&) { return "<GooseFEM.Mesh.Renumber>"; });

    py::class_<GooseFEM::Mesh::Reorder>(m, "Reorder")

        .def(py::init([](xt::xtensor<size_t, 1>& a) { return new GooseFEM::Mesh::Reorder({a}); }),
             "Reorder to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(py::init([](xt::xtensor<size_t, 1>& a, xt::xtensor<size_t, 1>& b) {
                return new GooseFEM::Mesh::Reorder({a, b});
             }),
             "Reorder to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(py::init(
                [](xt::xtensor<size_t, 1>& a, xt::xtensor<size_t, 1>& b, xt::xtensor<size_t, 1>& c) {
                    return new GooseFEM::Mesh::Reorder({a, b, c});
                }),
             "Reorder to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def(py::init([](xt::xtensor<size_t, 1>& a,
                         xt::xtensor<size_t, 1>& b,
                         xt::xtensor<size_t, 1>& c,
                         xt::xtensor<size_t, 1>& d) {
                return new GooseFEM::Mesh::Reorder({a, b, c, d});
                }),
             "Reorder to lowest possible index."
             "See :cpp:class:`GooseFEM::Mesh::Reorder`.")

        .def("get",
             &GooseFEM::Mesh::Reorder::apply<xt::xarray<size_t>>,
             "Reorder matrix (e.g. ``dofs``)."
             "See :cpp:func:`GooseFEM::Mesh::Reorder::get`.",
             py::arg("dofs"))

        .def("apply",
             &GooseFEM::Mesh::Reorder::apply<xt::xarray<size_t>>,
             "Get reordered list."
             "See :cpp:func:`GooseFEM::Mesh::Reorder::apply`.")

        .def("index",
             &GooseFEM::Mesh::Reorder::index,
             "Get index list to apply renumbering. Apply renumbering using ``index[dofs]``."
             "See :cpp:func:`GooseFEM::Mesh::Reorder::index`.")

        .def("__repr__", [](const GooseFEM::Mesh::Reorder&) { return "<GooseFEM.Mesh.Reorder>"; });

    m.def("dofs",
          &GooseFEM::Mesh::dofs,
          "List with DOF-numbers in sequential order."
          "See :cpp:func:`GooseFEM::Mesh::dofs`.",
          py::arg("nnode"),
          py::arg("ndim"));

    m.def("renumber",
          &GooseFEM::Mesh::renumber<xt::xtensor<size_t, 2>>,
          "Renumber to lowest possible indices."
          "See :cpp:func:`GooseFEM::Mesh::renumber`.",
          py::arg("dofs"));

    m.def("coordination",
          &GooseFEM::Mesh::coordination<xt::xtensor<size_t, 2>>,
          "Coordination number of each node."
          "See :cpp:func:`GooseFEM::Mesh::coordination`.",
          py::arg("conn"));

    m.def("elem2node",
          &GooseFEM::Mesh::elem2node<xt::xtensor<size_t, 2>>,
          "Element-numbers connected to each node."
          "See :cpp:func:`GooseFEM::Mesh::elem2node`.",
          py::arg("conn"),
          py::arg("sorted") = true);

    m.def("edgesize",
          py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&>(
            &GooseFEM::Mesh::edgesize<xt::xtensor<double, 2>, xt::xtensor<size_t, 2>>),
          "Get the edge size of all elements."
          "See :cpp:func:`GooseFEM::Mesh::edgesize`.",
          py::arg("coor"),
          py::arg("conn"));

    m.def("edgesize",
          py::overload_cast<
            const xt::xtensor<double, 2>&,
            const xt::xtensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::edgesize<xt::xtensor<double, 2>, xt::xtensor<size_t, 2>>),
          "Get the edge size of all elements."
          "See :cpp:func:`GooseFEM::Mesh::edgesize`.",
          py::arg("coor"),
          py::arg("conn"),
          py::arg("type"));

    m.def("centers",
          py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&>(
            &GooseFEM::Mesh::centers<xt::xtensor<double, 2>, xt::xtensor<size_t, 2>>),
          "Coordinates of the center of each element."
          "See :cpp:func:`GooseFEM::Mesh::centers`.",
          py::arg("coor"),
          py::arg("conn"));

    m.def("centers",
          py::overload_cast<
            const xt::xtensor<double, 2>&,
            const xt::xtensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::centers<xt::xtensor<double, 2>, xt::xtensor<size_t, 2>>),
          "Coordinates of the center of each element."
          "See :cpp:func:`GooseFEM::Mesh::centers`.",
          py::arg("coor"),
          py::arg("conn"),
          py::arg("type"));

    m.def("elemmap2nodemap",
          py::overload_cast<
            const xt::xtensor<size_t, 1>&,
            const xt::xtensor<double, 2>&,
            const xt::xtensor<size_t, 2>&>(&GooseFEM::Mesh::elemmap2nodemap<
                xt::xtensor<size_t, 1>, xt::xtensor<double, 2>, xt::xtensor<size_t, 2>>),
          "Convert an element-map to a node-map."
          "See :cpp:func:`GooseFEM::Mesh::elemmap2nodemap`.",
          py::arg("elem_map"),
          py::arg("coor"),
          py::arg("conn"));

    m.def("elemmap2nodemap",
          py::overload_cast<
            const xt::xtensor<size_t, 1>&,
            const xt::xtensor<double, 2>&,
            const xt::xtensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::elemmap2nodemap<
                xt::xtensor<size_t, 1>, xt::xtensor<double, 2>, xt::xtensor<size_t, 2>>),
          "Convert an element-map to a node-map."
          "See :cpp:func:`GooseFEM::Mesh::elemmap2nodemap`.",
          py::arg("elem_map"),
          py::arg("coor"),
          py::arg("conn"),
          py::arg("type"));
}

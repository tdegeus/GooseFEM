/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_Mesh(py::module& m)
{

    py::enum_<GooseFEM::Mesh::ElementType>(m, "ElementType", "ElementType")
        .value("Tri3", GooseFEM::Mesh::ElementType::Tri3)
        .value("Quad4", GooseFEM::Mesh::ElementType::Quad4)
        .value("Hex8", GooseFEM::Mesh::ElementType::Hex8)
        .export_values();

    py::class_<GooseFEM::Mesh::Renumber>(m, "Renumber")

        .def(
            py::init<const xt::xarray<size_t>&>(),
            "Class to renumber to the lowest possible indices. Use ``Renumber(...).get()`` to get "
            "the renumbered result",
            py::arg("dofs"))

        .def("get", &GooseFEM::Mesh::Renumber::get, "Get result of renumbering")

        .def(
            "index",
            &GooseFEM::Mesh::Renumber::index,
            "Get index list to apply renumbering. Apply renumbering using ``index[dofs]``")

        .def(
            "__repr__", [](const GooseFEM::Mesh::Renumber&) { return "<GooseFEM.Mesh.Renumber>"; });

    py::class_<GooseFEM::Mesh::Reorder>(m, "Reorder")

        .def(py::init([](xt::xtensor<size_t, 1>& a) { return new GooseFEM::Mesh::Reorder({a}); }))

        .def(py::init([](xt::xtensor<size_t, 1>& a, xt::xtensor<size_t, 1>& b) {
            return new GooseFEM::Mesh::Reorder({a, b});
        }))

        .def(py::init(
            [](xt::xtensor<size_t, 1>& a, xt::xtensor<size_t, 1>& b, xt::xtensor<size_t, 1>& c) {
                return new GooseFEM::Mesh::Reorder({a, b, c});
            }))

        .def(py::init([](xt::xtensor<size_t, 1>& a,
                         xt::xtensor<size_t, 1>& b,
                         xt::xtensor<size_t, 1>& c,
                         xt::xtensor<size_t, 1>& d) {
            return new GooseFEM::Mesh::Reorder({a, b, c, d});
        }))

        .def("get", &GooseFEM::Mesh::Reorder::get, "Reorder matrix (e.g. ``dofs``)")

        .def(
            "index",
            &GooseFEM::Mesh::Reorder::index,
            "Get index list to apply renumbering. Apply renumbering using ``index[dofs]``")

        .def("__repr__", [](const GooseFEM::Mesh::Reorder&) { return "<GooseFEM.Mesh.Reorder>"; });

    m.def(
        "dofs",
        &GooseFEM::Mesh::dofs,
        "List with DOF-numbers (in sequential order)",
        py::arg("nnode"),
        py::arg("ndim"));

    m.def(
        "renumber",
        &GooseFEM::Mesh::renumber,
        "Renumber to lowest possible indices. Use ``GooseFEM.Mesh.Renumber`` for advanced "
        "functionality",
        py::arg("dofs"));

    m.def(
        "coordination",
        &GooseFEM::Mesh::coordination,
        "Coordination number of each node (number of elements connected to each node)",
        py::arg("conn"));

    m.def(
        "elem2node",
        &GooseFEM::Mesh::elem2node,
        "Element-numbers connected to each node",
        py::arg("conn"),
        py::arg("sorted") = true);

    m.def(
        "edgesize",
        py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&>(
            &GooseFEM::Mesh::edgesize),
        "Get the edge size of all elements",
        py::arg("coor"),
        py::arg("conn"));

    m.def(
        "edgesize",
        py::overload_cast<
            const xt::xtensor<double, 2>&,
            const xt::xtensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::edgesize),
        "Get the edge size of all elements",
        py::arg("coor"),
        py::arg("conn"),
        py::arg("type"));

    m.def(
        "centers",
        py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&>(
            &GooseFEM::Mesh::centers),
        "Coordinates of the center of each element",
        py::arg("coor"),
        py::arg("conn"));

    m.def(
        "centers",
        py::overload_cast<
            const xt::xtensor<double, 2>&,
            const xt::xtensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::centers),
        "Coordinates of the center of each element",
        py::arg("coor"),
        py::arg("conn"),
        py::arg("type"));

    m.def(
        "elemmap2nodemap",
        py::overload_cast<
            const xt::xtensor<size_t, 1>&,
            const xt::xtensor<double, 2>&,
            const xt::xtensor<size_t, 2>&>(&GooseFEM::Mesh::elemmap2nodemap),
        "Convert an element-map to a node-map",
        py::arg("elem_map"),
        py::arg("coor"),
        py::arg("conn"));

    m.def(
        "elemmap2nodemap",
        py::overload_cast<
            const xt::xtensor<size_t, 1>&,
            const xt::xtensor<double, 2>&,
            const xt::xtensor<size_t, 2>&,
            GooseFEM::Mesh::ElementType>(&GooseFEM::Mesh::elemmap2nodemap),
        "Convert an element-map to a node-map",
        py::arg("elem_map"),
        py::arg("coor"),
        py::arg("conn"),
        py::arg("type"));
}

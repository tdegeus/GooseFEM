/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#include <GooseFEM/Vector.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_Vector(py::module& m)
{

    py::class_<GooseFEM::Vector>(m, "Vector")

        .def(
            py::init<const xt::pytensor<size_t, 2>&, const xt::pytensor<size_t, 2>&>(),
            "See :cpp:class:`GooseFEM::Vector`.",
            py::arg("conn"),
            py::arg("dofs"))

        .def_property_readonly("nelem", &GooseFEM::Vector::nelem)
        .def_property_readonly("nne", &GooseFEM::Vector::nne)
        .def_property_readonly("nnode", &GooseFEM::Vector::nnode)
        .def_property_readonly("ndim", &GooseFEM::Vector::ndim)
        .def_property_readonly("ndof", &GooseFEM::Vector::ndof)
        .def_property_readonly("conn", &GooseFEM::Vector::conn)
        .def_property_readonly("dofs", &GooseFEM::Vector::dofs)

        .def(
            "copy",
            &GooseFEM::Vector::copy<xt::pyarray<double>>,
            py::arg("nodevec_src"),
            py::arg("nodevec_dest"))

        .def(
            "Copy",
            &GooseFEM::Vector::Copy<xt::pyarray<double>>,
            py::arg("nodevec_src"),
            py::arg("nodevec_dest"))

        .def("AsDofs", &GooseFEM::Vector::AsDofs<xt::pyarray<double>>, py::arg("arg"))

        .def(
            "asDofs",
            &GooseFEM::Vector::asDofs<xt::pyarray<double>, xt::pytensor<double, 1>>,
            py::arg("arg"),
            py::arg("ret"))

        .def("AsNode", &GooseFEM::Vector::AsNode<xt::pyarray<double>>, py::arg("arg"))

        .def(
            "asNode",
            &GooseFEM::Vector::asNode<xt::pyarray<double>, xt::pytensor<double, 2>>,
            py::arg("arg"),
            py::arg("ret"))

        .def("AsElement", &GooseFEM::Vector::AsElement<xt::pyarray<double>>, py::arg("arg"))

        .def(
            "asElement",
            &GooseFEM::Vector::asElement<xt::pyarray<double>, xt::pytensor<double, 3>>,
            py::arg("arg"),
            py::arg("ret"))

        .def("AssembleDofs", &GooseFEM::Vector::AssembleDofs<xt::pyarray<double>>, py::arg("arg"))

        .def(
            "assembleDofs",
            &GooseFEM::Vector::assembleDofs<xt::pyarray<double>, xt::pytensor<double, 1>>,
            py::arg("arg"),
            py::arg("ret"))

        .def("AssembleNode", &GooseFEM::Vector::AssembleNode<xt::pyarray<double>>, py::arg("arg"))

        .def(
            "assembleNode",
            &GooseFEM::Vector::assembleNode<xt::pyarray<double>, xt::pytensor<double, 2>>,
            py::arg("arg"),
            py::arg("ret"))

        .def_property_readonly("shape_dofval", &GooseFEM::Vector::shape_dofval)
        .def_property_readonly("shape_nodevec", &GooseFEM::Vector::shape_nodevec)
        .def_property_readonly("shape_elemvec", &GooseFEM::Vector::shape_elemvec)
        .def_property_readonly("shape_elemmat", &GooseFEM::Vector::shape_elemmat)

        .def("allocate_dofval", py::overload_cast<>(&GooseFEM::Vector::allocate_dofval, py::const_))

        .def(
            "allocate_dofval",
            py::overload_cast<double>(&GooseFEM::Vector::allocate_dofval, py::const_))

        .def(
            "allocate_nodevec",
            py::overload_cast<>(&GooseFEM::Vector::allocate_nodevec, py::const_))

        .def(
            "allocate_nodevec",
            py::overload_cast<double>(&GooseFEM::Vector::allocate_nodevec, py::const_))

        .def(
            "allocate_elemvec",
            py::overload_cast<>(&GooseFEM::Vector::allocate_elemvec, py::const_))

        .def(
            "allocate_elemvec",
            py::overload_cast<double>(&GooseFEM::Vector::allocate_elemvec, py::const_))

        .def(
            "allocate_elemmat",
            py::overload_cast<>(&GooseFEM::Vector::allocate_elemmat, py::const_))

        .def(
            "allocate_elemmat",
            py::overload_cast<double>(&GooseFEM::Vector::allocate_elemmat, py::const_))

        .def("__repr__", [](const GooseFEM::Vector&) { return "<GooseFEM.Vector>"; });
}

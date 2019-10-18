/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>
#define GOOSEFEM_NO_HIGHFIVE
#include "../include/GooseFEM/GooseFEM.h"
#include "../include/GooseFEM/ParaView.h"

// =================================================================================================

namespace py = pybind11;

// =================================================================================================

void init_ParaViewHDF5(py::module& m)
{

m.def("as3d", &GooseFEM::ParaView::HDF5::as3d);

// -------------------------------------------------------------------------------------------------

py::enum_<GooseFEM::ParaView::HDF5::ElementType>(m, "ElementType", "ElementType")
    .value("Triangle", GooseFEM::ParaView::HDF5::ElementType::Triangle)
    .value("Quadrilateral", GooseFEM::ParaView::HDF5::ElementType::Quadrilateral)
    .value("Hexahedron", GooseFEM::ParaView::HDF5::ElementType::Hexahedron)
    .export_values();

// -------------------------------------------------------------------------------------------------

py::enum_<GooseFEM::ParaView::HDF5::AttributeType>(m, "AttributeType", "AttributeType")
    .value("Cell", GooseFEM::ParaView::HDF5::AttributeType::Cell)
    .value("Node", GooseFEM::ParaView::HDF5::AttributeType::Node)
    .export_values();

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::ParaView::HDF5::Connectivity>(m, "Connectivity")

    .def(py::init<
            const std::string&,
            const std::string&,
            GooseFEM::ParaView::HDF5::ElementType,
            const std::vector<size_t>&>(),
        "Connectivity",
        py::arg("fname"),
        py::arg("dataset"),
        py::arg("ElementType"),
        py::arg("shape"))

    .def(py::init<
            const std::string&,
            const std::string&,
            GooseFEM::Mesh::ElementType,
            const std::vector<size_t>&>(),
        "Connectivity",
        py::arg("fname"),
        py::arg("dataset"),
        py::arg("ElementType"),
        py::arg("shape"))

    .def("nelem",
        &GooseFEM::ParaView::HDF5::Connectivity::nelem)

    .def("nne",
        &GooseFEM::ParaView::HDF5::Connectivity::nne)

    .def("shape",
        &GooseFEM::ParaView::HDF5::Connectivity::shape)

    .def("fname",
        &GooseFEM::ParaView::HDF5::Connectivity::fname)

    .def("xdmf",
        &GooseFEM::ParaView::HDF5::Connectivity::xdmf,
        py::arg("indent")=4)

    .def("__repr__",
        [](const GooseFEM::ParaView::HDF5::Connectivity &){
            return "<GooseFEM.ParaView.HDF5.Connectivity>"; });

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::ParaView::HDF5::Coordinates>(m, "Coordinates")

    .def(py::init<
            const std::string&,
            const std::string&,
            const std::vector<size_t>&>(),
        "Coordinates",
        py::arg("fname"),
        py::arg("dataset"),
        py::arg("shape"))

    .def("nnode", &GooseFEM::ParaView::HDF5::Coordinates::nnode)
    .def("ndim", &GooseFEM::ParaView::HDF5::Coordinates::ndim)
    .def("shape", &GooseFEM::ParaView::HDF5::Coordinates::shape)
    .def("fname", &GooseFEM::ParaView::HDF5::Coordinates::fname)
    .def("xdmf", &GooseFEM::ParaView::HDF5::Coordinates::xdmf, py::arg("indent")=4)

    .def("__repr__",
        [](const GooseFEM::ParaView::HDF5::Coordinates &){
            return "<GooseFEM.ParaView.HDF5.Coordinates>"; });

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::ParaView::HDF5::Attribute>(m, "Attribute")

    .def(py::init<
            const std::string&,
            const std::string&,
            const std::string&, GooseFEM::ParaView::HDF5::AttributeType,
            const std::vector<size_t>&>(),
        "Attribute",
        py::arg("fname"),
        py::arg("dataset"),
        py::arg("name"),
        py::arg("AttributeType"),
        py::arg("shape"))

    .def("shape",
        &GooseFEM::ParaView::HDF5::Attribute::shape)

    .def("fname",
        &GooseFEM::ParaView::HDF5::Attribute::fname)

    .def("xdmf",
        &GooseFEM::ParaView::HDF5::Attribute::xdmf,
        py::arg("indent")=4)

    .def("__repr__",
        [](const GooseFEM::ParaView::HDF5::Attribute &){
            return "<GooseFEM.ParaView.HDF5.Attribute>"; });

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::ParaView::HDF5::Mesh>(m, "Mesh")

    .def(py::init<
            const GooseFEM::ParaView::HDF5::Connectivity&,
            const GooseFEM::ParaView::HDF5::Coordinates&>(),
        "Mesh",
        py::arg("conn"),
        py::arg("coor"))

    .def("push_back",
        py::overload_cast<const GooseFEM::ParaView::HDF5::Attribute&>(
            &GooseFEM::ParaView::HDF5::Mesh::push_back))

    .def("xdmf",
        &GooseFEM::ParaView::HDF5::Mesh::xdmf,
        py::arg("indent")=4)

    .def("write",
        &GooseFEM::ParaView::HDF5::Mesh::write,
        py::arg("fname"),
        py::arg("indent")=4)

    .def("__repr__",
        [](const GooseFEM::ParaView::HDF5::Mesh &){
            return "<GooseFEM.ParaView.HDF5.Mesh>"; });

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::ParaView::HDF5::Increment>(m, "Increment")

    .def(py::init<
            const GooseFEM::ParaView::HDF5::Connectivity&,
            const GooseFEM::ParaView::HDF5::Coordinates&>(),
        "Increment",
        py::arg("conn"),
        py::arg("coor"))

    .def("push_back",
        py::overload_cast<const GooseFEM::ParaView::HDF5::Connectivity&>(
            &GooseFEM::ParaView::HDF5::Increment::push_back))

    .def("push_back",
        py::overload_cast<const GooseFEM::ParaView::HDF5::Coordinates&>(
            &GooseFEM::ParaView::HDF5::Increment::push_back))

    .def("push_back",
        py::overload_cast<const GooseFEM::ParaView::HDF5::Attribute&>(
            &GooseFEM::ParaView::HDF5::Increment::push_back))

    .def("xdmf",
        &GooseFEM::ParaView::HDF5::Increment::xdmf,
        py::arg("indent")=4)

    .def("__repr__",
        [](const GooseFEM::ParaView::HDF5::Increment &){
            return "<GooseFEM.ParaView.HDF5.Increment>"; });

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::ParaView::HDF5::TimeSeries>(m, "TimeSeries")

    .def(py::init<>(),
        "TimeSeries")

    .def("push_back",
        &GooseFEM::ParaView::HDF5::TimeSeries::push_back)

    .def("xdmf",
        &GooseFEM::ParaView::HDF5::TimeSeries::xdmf,
        py::arg("indent")=4)

    .def("write",
        &GooseFEM::ParaView::HDF5::TimeSeries::write,
        py::arg("fname"),
        py::arg("indent")=4)

    .def("__repr__",
        [](const GooseFEM::ParaView::HDF5::TimeSeries &){
            return "<GooseFEM.ParaView.HDF5.TimeSeries>"; });

}

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <pyxtensor/pyxtensor.hpp>

#include "../include/GooseFEM/GooseFEM.h"

// =================================================================================================

namespace py = pybind11;

// =================================================================================================

void init_MeshTri3(py::module& m)
{

py::class_<GooseFEM::Mesh::Tri3::Regular>(m, "Regular")

    .def(
        py::init<size_t,size_t,double>(),
        "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
        py::arg("nx"),
        py::arg("ny"),
        py::arg("h")=1.)

    .def("coor",
        &GooseFEM::Mesh::Tri3::Regular::coor)

    .def("conn",
        &GooseFEM::Mesh::Tri3::Regular::conn)

    .def("nelem",
        &GooseFEM::Mesh::Tri3::Regular::nelem)

    .def("nnode",
        &GooseFEM::Mesh::Tri3::Regular::nnode)

    .def("nne",
        &GooseFEM::Mesh::Tri3::Regular::nne)

    .def("ndim",
        &GooseFEM::Mesh::Tri3::Regular::ndim)

    .def("nodesBottomEdge",
        &GooseFEM::Mesh::Tri3::Regular::nodesBottomEdge)

    .def("nodesTopEdge",
        &GooseFEM::Mesh::Tri3::Regular::nodesTopEdge)

    .def("nodesLeftEdge",
        &GooseFEM::Mesh::Tri3::Regular::nodesLeftEdge)

    .def("nodesRightEdge",
        &GooseFEM::Mesh::Tri3::Regular::nodesRightEdge)

    .def("nodesBottomOpenEdge",
        &GooseFEM::Mesh::Tri3::Regular::nodesBottomOpenEdge)

    .def("nodesTopOpenEdge",
        &GooseFEM::Mesh::Tri3::Regular::nodesTopOpenEdge)

    .def("nodesLeftOpenEdge",
        &GooseFEM::Mesh::Tri3::Regular::nodesLeftOpenEdge)

    .def("nodesRightOpenEdge",
        &GooseFEM::Mesh::Tri3::Regular::nodesRightOpenEdge)

    .def("nodesBottomLeftCorner",
        &GooseFEM::Mesh::Tri3::Regular::nodesBottomLeftCorner)

    .def("nodesBottomRightCorner",
        &GooseFEM::Mesh::Tri3::Regular::nodesBottomRightCorner)

    .def("nodesTopLeftCorner",
        &GooseFEM::Mesh::Tri3::Regular::nodesTopLeftCorner)

    .def("nodesTopRightCorner",
        &GooseFEM::Mesh::Tri3::Regular::nodesTopRightCorner)

    .def("nodesLeftBottomCorner",
        &GooseFEM::Mesh::Tri3::Regular::nodesLeftBottomCorner)

    .def("nodesLeftTopCorner",
        &GooseFEM::Mesh::Tri3::Regular::nodesLeftTopCorner)

    .def("nodesRightBottomCorner",
        &GooseFEM::Mesh::Tri3::Regular::nodesRightBottomCorner)

    .def("nodesRightTopCorner",
        &GooseFEM::Mesh::Tri3::Regular::nodesRightTopCorner)

    .def("nodesPeriodic",
        &GooseFEM::Mesh::Tri3::Regular::nodesPeriodic)

    .def("nodesOrigin",
        &GooseFEM::Mesh::Tri3::Regular::nodesOrigin)

    .def("dofs",
        &GooseFEM::Mesh::Tri3::Regular::dofs)

    .def("dofsPeriodic",
        &GooseFEM::Mesh::Tri3::Regular::dofsPeriodic)

    .def("__repr__",
        [](const GooseFEM::Mesh::Tri3::Regular&){ return "<GooseFEM.Mesh.Tri3.Regular>"; });

// -------------------------------------------------------------------------------------------------

m.def("getOrientation",
    &GooseFEM::Mesh::Tri3::getOrientation,
    "Get the orientation of each element",
    py::arg("coor"),
    py::arg("conn"));

m.def("retriangulate",
    &GooseFEM::Mesh::Tri3::retriangulate,
    "Re-triangulate existing mesh",
    py::arg("coor"),
    py::arg("conn"),
    py::arg("orientation")=-1);

}

// =================================================================================================


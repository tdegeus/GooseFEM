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

void init_MeshHex8(py::module& m)
{

py::class_<GooseFEM::Mesh::Hex8::Regular>(m, "Regular")

    .def(py::init<size_t,size_t,size_t,double>(),
        "Mesh with nx*ny*nz 'pixels' and edge size h",
        py::arg("nx"),
        py::arg("ny"),
        py::arg("nz"),
        py::arg("h")=1.)

    .def("nelem",
        &GooseFEM::Mesh::Hex8::Regular::nelem)

    .def("nnode",
        &GooseFEM::Mesh::Hex8::Regular::nnode)

    .def("nne",
        &GooseFEM::Mesh::Hex8::Regular::nne)

    .def("ndim",
        &GooseFEM::Mesh::Hex8::Regular::ndim)

    .def("coor",
        &GooseFEM::Mesh::Hex8::Regular::coor)

    .def("conn",
        &GooseFEM::Mesh::Hex8::Regular::conn)

    .def("getElementType",
        &GooseFEM::Mesh::Hex8::Regular::getElementType)

    .def("nodesFront",
        &GooseFEM::Mesh::Hex8::Regular::nodesFront)

    .def("nodesBack",
        &GooseFEM::Mesh::Hex8::Regular::nodesBack)

    .def("nodesLeft",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeft)

    .def("nodesRight",
        &GooseFEM::Mesh::Hex8::Regular::nodesRight)

    .def("nodesBottom",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottom)

    .def("nodesTop",
        &GooseFEM::Mesh::Hex8::Regular::nodesTop)

    .def("nodesFrontFace",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontFace)

    .def("nodesBackFace",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackFace)

    .def("nodesLeftFace",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftFace)

    .def("nodesRightFace",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightFace)

    .def("nodesBottomFace",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomFace)

    .def("nodesTopFace",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopFace)

    .def("nodesFrontBottomEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomEdge)

    .def("nodesFrontTopEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopEdge)

    .def("nodesFrontLeftEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftEdge)

    .def("nodesFrontRightEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightEdge)

    .def("nodesBackBottomEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomEdge)

    .def("nodesBackTopEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackTopEdge)

    .def("nodesBackLeftEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftEdge)

    .def("nodesBackRightEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackRightEdge)

    .def("nodesBottomLeftEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftEdge)

    .def("nodesBottomRightEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightEdge)

    .def("nodesTopLeftEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftEdge)

    .def("nodesTopRightEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopRightEdge)

    .def("nodesBottomFrontEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontEdge)

    .def("nodesBottomBackEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackEdge)

    .def("nodesTopFrontEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontEdge)

    .def("nodesTopBackEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopBackEdge)

    .def("nodesLeftBottomEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomEdge)

    .def("nodesLeftFrontEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontEdge)

    .def("nodesLeftBackEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackEdge)

    .def("nodesLeftTopEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopEdge)

    .def("nodesRightBottomEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomEdge)

    .def("nodesRightTopEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightTopEdge)

    .def("nodesRightFrontEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontEdge)

    .def("nodesRightBackEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightBackEdge)

    .def("nodesFrontBottomOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomOpenEdge)

    .def("nodesFrontTopOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopOpenEdge)

    .def("nodesFrontLeftOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftOpenEdge)

    .def("nodesFrontRightOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightOpenEdge)

    .def("nodesBackBottomOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomOpenEdge)

    .def("nodesBackTopOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackTopOpenEdge)

    .def("nodesBackLeftOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftOpenEdge)

    .def("nodesBackRightOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackRightOpenEdge)

    .def("nodesBottomLeftOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftOpenEdge)

    .def("nodesBottomRightOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightOpenEdge)

    .def("nodesTopLeftOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftOpenEdge)

    .def("nodesTopRightOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopRightOpenEdge)

    .def("nodesBottomFrontOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontOpenEdge)

    .def("nodesBottomBackOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackOpenEdge)

    .def("nodesTopFrontOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontOpenEdge)

    .def("nodesTopBackOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopBackOpenEdge)

    .def("nodesLeftBottomOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomOpenEdge)

    .def("nodesLeftFrontOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontOpenEdge)

    .def("nodesLeftBackOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackOpenEdge)

    .def("nodesLeftTopOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopOpenEdge)

    .def("nodesRightBottomOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomOpenEdge)

    .def("nodesRightTopOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightTopOpenEdge)

    .def("nodesRightFrontOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontOpenEdge)

    .def("nodesRightBackOpenEdge",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightBackOpenEdge)

    .def("nodesFrontBottomLeftCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomLeftCorner)

    .def("nodesFrontBottomRightCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomRightCorner)

    .def("nodesFrontTopLeftCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopLeftCorner)

    .def("nodesFrontTopRightCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopRightCorner)

    .def("nodesBackBottomLeftCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomLeftCorner)

    .def("nodesBackBottomRightCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomRightCorner)

    .def("nodesBackTopLeftCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackTopLeftCorner)

    .def("nodesBackTopRightCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackTopRightCorner)

    .def("nodesFrontLeftBottomCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftBottomCorner)

    .def("nodesBottomFrontLeftCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontLeftCorner)

    .def("nodesBottomLeftFrontCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftFrontCorner)

    .def("nodesLeftFrontBottomCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontBottomCorner)

    .def("nodesLeftBottomFrontCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomFrontCorner)

    .def("nodesFrontRightBottomCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightBottomCorner)

    .def("nodesBottomFrontRightCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontRightCorner)

    .def("nodesBottomRightFrontCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightFrontCorner)

    .def("nodesRightFrontBottomCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontBottomCorner)

    .def("nodesRightBottomFrontCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomFrontCorner)

    .def("nodesFrontLeftTopCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftTopCorner)

    .def("nodesTopFrontLeftCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontLeftCorner)

    .def("nodesTopLeftFrontCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftFrontCorner)

    .def("nodesLeftFrontTopCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontTopCorner)

    .def("nodesLeftTopFrontCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopFrontCorner)

    .def("nodesFrontRightTopCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightTopCorner)

    .def("nodesTopFrontRightCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontRightCorner)

    .def("nodesTopRightFrontCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopRightFrontCorner)

    .def("nodesRightFrontTopCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontTopCorner)

    .def("nodesRightTopFrontCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightTopFrontCorner)

    .def("nodesBackLeftBottomCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftBottomCorner)

    .def("nodesBottomBackLeftCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackLeftCorner)

    .def("nodesBottomLeftBackCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftBackCorner)

    .def("nodesLeftBackBottomCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackBottomCorner)

    .def("nodesLeftBottomBackCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomBackCorner)

    .def("nodesBackRightBottomCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackRightBottomCorner)

    .def("nodesBottomBackRightCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackRightCorner)

    .def("nodesBottomRightBackCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightBackCorner)

    .def("nodesRightBackBottomCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightBackBottomCorner)

    .def("nodesRightBottomBackCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomBackCorner)

    .def("nodesBackLeftTopCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftTopCorner)

    .def("nodesTopBackLeftCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopBackLeftCorner)

    .def("nodesTopLeftBackCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftBackCorner)

    .def("nodesLeftBackTopCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackTopCorner)

    .def("nodesLeftTopBackCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopBackCorner)

    .def("nodesBackRightTopCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesBackRightTopCorner)

    .def("nodesTopBackRightCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopBackRightCorner)

    .def("nodesTopRightBackCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesTopRightBackCorner)

    .def("nodesRightBackTopCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightBackTopCorner)

    .def("nodesRightTopBackCorner",
        &GooseFEM::Mesh::Hex8::Regular::nodesRightTopBackCorner)

    .def("nodesPeriodic",
        &GooseFEM::Mesh::Hex8::Regular::nodesPeriodic)

    .def("nodesOrigin",
        &GooseFEM::Mesh::Hex8::Regular::nodesOrigin)

    .def("dofs",
        &GooseFEM::Mesh::Hex8::Regular::dofs)

    .def("dofsPeriodic",
        &GooseFEM::Mesh::Hex8::Regular::dofsPeriodic)

    .def("__repr__",
        [](const GooseFEM::Mesh::Hex8::Regular &){ return "<GooseFEM.Mesh.Hex8.Regular>"; });

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Hex8::FineLayer>(m, "FineLayer")

    .def(py::init<size_t,size_t,size_t,double,size_t>(),
        "Mesh with nx*ny*nz 'pixels' and edge size h",
        py::arg("nx"),
        py::arg("ny"),
        py::arg("nz"),
        py::arg("h")=1.,
        py::arg("nfine")=1)

    .def("nelem",
        &GooseFEM::Mesh::Hex8::FineLayer::nelem)

    .def("nnode",
        &GooseFEM::Mesh::Hex8::FineLayer::nnode)

    .def("nne",
        &GooseFEM::Mesh::Hex8::FineLayer::nne)

    .def("ndim",
        &GooseFEM::Mesh::Hex8::FineLayer::ndim)

    .def("shape",
        &GooseFEM::Mesh::Hex8::FineLayer::shape)

    .def("coor",
        &GooseFEM::Mesh::Hex8::FineLayer::coor)

    .def("conn",
        &GooseFEM::Mesh::Hex8::FineLayer::conn)

    .def("getElementType",
        &GooseFEM::Mesh::Hex8::FineLayer::getElementType)

    .def("elementsMiddleLayer",
        &GooseFEM::Mesh::Hex8::FineLayer::elementsMiddleLayer)

    .def("nodesFront",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFront)

    .def("nodesBack",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBack)

    .def("nodesLeft",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeft)

    .def("nodesRight",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRight)

    .def("nodesBottom",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottom)

    .def("nodesTop",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTop)

    .def("nodesFrontFace",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontFace)

    .def("nodesBackFace",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackFace)

    .def("nodesLeftFace",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFace)

    .def("nodesRightFace",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFace)

    .def("nodesBottomFace",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFace)

    .def("nodesTopFace",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFace)

    .def("nodesFrontBottomEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomEdge)

    .def("nodesFrontTopEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopEdge)

    .def("nodesFrontLeftEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftEdge)

    .def("nodesFrontRightEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightEdge)

    .def("nodesBackBottomEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomEdge)

    .def("nodesBackTopEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopEdge)

    .def("nodesBackLeftEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftEdge)

    .def("nodesBackRightEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightEdge)

    .def("nodesBottomLeftEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftEdge)

    .def("nodesBottomRightEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightEdge)

    .def("nodesTopLeftEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftEdge)

    .def("nodesTopRightEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightEdge)

    .def("nodesBottomFrontEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontEdge)

    .def("nodesBottomBackEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackEdge)

    .def("nodesTopFrontEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontEdge)

    .def("nodesTopBackEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackEdge)

    .def("nodesLeftBottomEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomEdge)

    .def("nodesLeftFrontEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontEdge)

    .def("nodesLeftBackEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackEdge)

    .def("nodesLeftTopEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopEdge)

    .def("nodesRightBottomEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomEdge)

    .def("nodesRightTopEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopEdge)

    .def("nodesRightFrontEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontEdge)

    .def("nodesRightBackEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackEdge)

    .def("nodesFrontBottomOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomOpenEdge)

    .def("nodesFrontTopOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopOpenEdge)

    .def("nodesFrontLeftOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftOpenEdge)

    .def("nodesFrontRightOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightOpenEdge)

    .def("nodesBackBottomOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomOpenEdge)

    .def("nodesBackTopOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopOpenEdge)

    .def("nodesBackLeftOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftOpenEdge)

    .def("nodesBackRightOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightOpenEdge)

    .def("nodesBottomLeftOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftOpenEdge)

    .def("nodesBottomRightOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightOpenEdge)

    .def("nodesTopLeftOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftOpenEdge)

    .def("nodesTopRightOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightOpenEdge)

    .def("nodesBottomFrontOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontOpenEdge)

    .def("nodesBottomBackOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackOpenEdge)

    .def("nodesTopFrontOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontOpenEdge)

    .def("nodesTopBackOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackOpenEdge)

    .def("nodesLeftBottomOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomOpenEdge)

    .def("nodesLeftFrontOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontOpenEdge)

    .def("nodesLeftBackOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackOpenEdge)

    .def("nodesLeftTopOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopOpenEdge)

    .def("nodesRightBottomOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomOpenEdge)

    .def("nodesRightTopOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopOpenEdge)

    .def("nodesRightFrontOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontOpenEdge)

    .def("nodesRightBackOpenEdge",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackOpenEdge)

    .def("nodesFrontBottomLeftCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomLeftCorner)

    .def("nodesFrontBottomRightCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomRightCorner)

    .def("nodesFrontTopLeftCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopLeftCorner)

    .def("nodesFrontTopRightCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopRightCorner)

    .def("nodesBackBottomLeftCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomLeftCorner)

    .def("nodesBackBottomRightCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomRightCorner)

    .def("nodesBackTopLeftCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopLeftCorner)

    .def("nodesBackTopRightCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopRightCorner)

    .def("nodesFrontLeftBottomCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftBottomCorner)

    .def("nodesBottomFrontLeftCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontLeftCorner)

    .def("nodesBottomLeftFrontCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftFrontCorner)

    .def("nodesLeftFrontBottomCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontBottomCorner)

    .def("nodesLeftBottomFrontCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomFrontCorner)

    .def("nodesFrontRightBottomCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightBottomCorner)

    .def("nodesBottomFrontRightCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontRightCorner)

    .def("nodesBottomRightFrontCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightFrontCorner)

    .def("nodesRightFrontBottomCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontBottomCorner)

    .def("nodesRightBottomFrontCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomFrontCorner)

    .def("nodesFrontLeftTopCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftTopCorner)

    .def("nodesTopFrontLeftCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontLeftCorner)

    .def("nodesTopLeftFrontCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftFrontCorner)

    .def("nodesLeftFrontTopCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontTopCorner)

    .def("nodesLeftTopFrontCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopFrontCorner)

    .def("nodesFrontRightTopCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightTopCorner)

    .def("nodesTopFrontRightCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontRightCorner)

    .def("nodesTopRightFrontCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightFrontCorner)

    .def("nodesRightFrontTopCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontTopCorner)

    .def("nodesRightTopFrontCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopFrontCorner)

    .def("nodesBackLeftBottomCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftBottomCorner)

    .def("nodesBottomBackLeftCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackLeftCorner)

    .def("nodesBottomLeftBackCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftBackCorner)

    .def("nodesLeftBackBottomCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackBottomCorner)

    .def("nodesLeftBottomBackCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomBackCorner)

    .def("nodesBackRightBottomCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightBottomCorner)

    .def("nodesBottomBackRightCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackRightCorner)

    .def("nodesBottomRightBackCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightBackCorner)

    .def("nodesRightBackBottomCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackBottomCorner)

    .def("nodesRightBottomBackCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomBackCorner)

    .def("nodesBackLeftTopCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftTopCorner)

    .def("nodesTopBackLeftCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackLeftCorner)

    .def("nodesTopLeftBackCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftBackCorner)

    .def("nodesLeftBackTopCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackTopCorner)

    .def("nodesLeftTopBackCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopBackCorner)

    .def("nodesBackRightTopCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightTopCorner)

    .def("nodesTopBackRightCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackRightCorner)

    .def("nodesTopRightBackCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightBackCorner)

    .def("nodesRightBackTopCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackTopCorner)

    .def("nodesRightTopBackCorner",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopBackCorner)

    .def("nodesPeriodic",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesPeriodic)

    .def("nodesOrigin",
        &GooseFEM::Mesh::Hex8::FineLayer::nodesOrigin)

    .def("dofs",
        &GooseFEM::Mesh::Hex8::FineLayer::dofs)

    .def("dofsPeriodic",
        &GooseFEM::Mesh::Hex8::FineLayer::dofsPeriodic)

    .def("__repr__",
        [](const GooseFEM::Mesh::Hex8::FineLayer &){ return "<GooseFEM.Mesh.Hex8.FineLayer>"; });

}

// =================================================================================================


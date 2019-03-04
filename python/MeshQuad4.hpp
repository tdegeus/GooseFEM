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

void init_MeshQuad4(py::module &m)
{

py::class_<GooseFEM::Mesh::Quad4::Regular>(m, "Regular")

  .def(py::init<size_t,size_t,double>(), "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'", py::arg("nx"), py::arg("ny"), py::arg("h")=1.)

  .def("coor"                  , &GooseFEM::Mesh::Quad4::Regular::coor                  )
  .def("conn"                  , &GooseFEM::Mesh::Quad4::Regular::conn                  )
  .def("nelem"                 , &GooseFEM::Mesh::Quad4::Regular::nelem                 )
  .def("nnode"                 , &GooseFEM::Mesh::Quad4::Regular::nnode                 )
  .def("nne"                   , &GooseFEM::Mesh::Quad4::Regular::nne                   )
  .def("ndim"                  , &GooseFEM::Mesh::Quad4::Regular::ndim                  )

  .def("nodesBottomEdge"       , &GooseFEM::Mesh::Quad4::Regular::nodesBottomEdge       )
  .def("nodesTopEdge"          , &GooseFEM::Mesh::Quad4::Regular::nodesTopEdge          )
  .def("nodesLeftEdge"         , &GooseFEM::Mesh::Quad4::Regular::nodesLeftEdge         )
  .def("nodesRightEdge"        , &GooseFEM::Mesh::Quad4::Regular::nodesRightEdge        )
  .def("nodesBottomOpenEdge"   , &GooseFEM::Mesh::Quad4::Regular::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &GooseFEM::Mesh::Quad4::Regular::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &GooseFEM::Mesh::Quad4::Regular::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &GooseFEM::Mesh::Quad4::Regular::nodesRightOpenEdge    )

  .def("nodesBottomLeftCorner" , &GooseFEM::Mesh::Quad4::Regular::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &GooseFEM::Mesh::Quad4::Regular::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &GooseFEM::Mesh::Quad4::Regular::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &GooseFEM::Mesh::Quad4::Regular::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &GooseFEM::Mesh::Quad4::Regular::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &GooseFEM::Mesh::Quad4::Regular::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &GooseFEM::Mesh::Quad4::Regular::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &GooseFEM::Mesh::Quad4::Regular::nodesRightTopCorner   )

  .def("nodesPeriodic"         , &GooseFEM::Mesh::Quad4::Regular::nodesPeriodic         )
  .def("nodesOrigin"           , &GooseFEM::Mesh::Quad4::Regular::nodesOrigin           )

  .def("dofs"                  , &GooseFEM::Mesh::Quad4::Regular::dofs                  )
  .def("dofsPeriodic"          , &GooseFEM::Mesh::Quad4::Regular::dofsPeriodic          )

  .def("__repr__", [](const GooseFEM::Mesh::Quad4::Regular &){ return "<GooseFEM.Mesh.Quad4.Regular>"; });

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Quad4::FineLayer>(m, "FineLayer")

  .def(
    py::init<size_t,size_t,double,size_t>(),
    "FineLayer mesh: 'nx' pixels in horizontal direction (length 'Lx'), idem in vertical direction",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.,
    py::arg("nfine")=1
  )

  .def("shape"                 , &GooseFEM::Mesh::Quad4::FineLayer::shape                 )
  .def("coor"                  , &GooseFEM::Mesh::Quad4::FineLayer::coor                  )
  .def("conn"                  , &GooseFEM::Mesh::Quad4::FineLayer::conn                  )
  .def("nelem"                 , &GooseFEM::Mesh::Quad4::FineLayer::nelem                 )
  .def("nnode"                 , &GooseFEM::Mesh::Quad4::FineLayer::nnode                 )
  .def("nne"                   , &GooseFEM::Mesh::Quad4::FineLayer::nne                   )
  .def("ndim"                  , &GooseFEM::Mesh::Quad4::FineLayer::ndim                  )
  .def("elementsMiddleLayer"   , &GooseFEM::Mesh::Quad4::FineLayer::elementsMiddleLayer   )
  .def("nodesBottomEdge"       , &GooseFEM::Mesh::Quad4::FineLayer::nodesBottomEdge       )
  .def("nodesTopEdge"          , &GooseFEM::Mesh::Quad4::FineLayer::nodesTopEdge          )
  .def("nodesLeftEdge"         , &GooseFEM::Mesh::Quad4::FineLayer::nodesLeftEdge         )
  .def("nodesRightEdge"        , &GooseFEM::Mesh::Quad4::FineLayer::nodesRightEdge        )
  .def("nodesBottomOpenEdge"   , &GooseFEM::Mesh::Quad4::FineLayer::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &GooseFEM::Mesh::Quad4::FineLayer::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &GooseFEM::Mesh::Quad4::FineLayer::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &GooseFEM::Mesh::Quad4::FineLayer::nodesRightOpenEdge    )
  .def("nodesBottomLeftCorner" , &GooseFEM::Mesh::Quad4::FineLayer::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &GooseFEM::Mesh::Quad4::FineLayer::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &GooseFEM::Mesh::Quad4::FineLayer::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &GooseFEM::Mesh::Quad4::FineLayer::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &GooseFEM::Mesh::Quad4::FineLayer::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &GooseFEM::Mesh::Quad4::FineLayer::nodesRightTopCorner   )
  .def("nodesPeriodic"         , &GooseFEM::Mesh::Quad4::FineLayer::nodesPeriodic         )
  .def("nodesOrigin"           , &GooseFEM::Mesh::Quad4::FineLayer::nodesOrigin           )
  .def("dofs"                  , &GooseFEM::Mesh::Quad4::FineLayer::dofs                  )
  .def("dofsPeriodic"          , &GooseFEM::Mesh::Quad4::FineLayer::dofsPeriodic          )

  .def("__repr__",
    [](const GooseFEM::Mesh::Quad4::FineLayer &){ return "<GooseFEM.Mesh.Quad4.FineLayer>"; }
  );

}

// =================================================================================================


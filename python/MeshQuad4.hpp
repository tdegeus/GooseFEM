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
namespace MD = GooseFEM::Mesh::Quad4;

// =================================================================================================

void init_MeshQuad4(py::module &m)
{

py::class_<MD::Regular>(m, "Regular")

  .def(py::init<size_t,size_t,double>(),
    "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.
  )

  .def("coor" , &MD::Regular::coor )
  .def("conn" , &MD::Regular::conn )
  .def("nelem", &MD::Regular::nelem)
  .def("nnode", &MD::Regular::nnode)
  .def("nne"  , &MD::Regular::nne  )
  .def("ndim" , &MD::Regular::ndim )
  .def("nelx" , &MD::Regular::nelx )
  .def("nely" , &MD::Regular::nely )
  .def("h"    , &MD::Regular::h    )

  .def("nodesBottomEdge"    , &MD::Regular::nodesBottomEdge    )
  .def("nodesTopEdge"       , &MD::Regular::nodesTopEdge       )
  .def("nodesLeftEdge"      , &MD::Regular::nodesLeftEdge      )
  .def("nodesRightEdge"     , &MD::Regular::nodesRightEdge     )
  .def("nodesBottomOpenEdge", &MD::Regular::nodesBottomOpenEdge)
  .def("nodesTopOpenEdge"   , &MD::Regular::nodesTopOpenEdge   )
  .def("nodesLeftOpenEdge"  , &MD::Regular::nodesLeftOpenEdge  )
  .def("nodesRightOpenEdge" , &MD::Regular::nodesRightOpenEdge )

  .def("nodesBottomLeftCorner" , &MD::Regular::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &MD::Regular::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &MD::Regular::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &MD::Regular::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &MD::Regular::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &MD::Regular::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &MD::Regular::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &MD::Regular::nodesRightTopCorner   )

  .def("dofs", &MD::Regular::dofs)

  .def("nodesPeriodic", &MD::Regular::nodesPeriodic)
  .def("nodesOrigin"  , &MD::Regular::nodesOrigin  )
  .def("dofsPeriodic" , &MD::Regular::dofsPeriodic )

  .def("elementMatrix", &MD::Regular::elementMatrix)

  .def("__repr__",
    [](const MD::Regular &){ return "<GooseFEM.Mesh.Quad4.Regular>"; });

// -------------------------------------------------------------------------------------------------

py::class_<MD::FineLayer>(m, "FineLayer")

  .def(
    py::init<size_t,size_t,double,size_t>(),
    "FineLayer mesh: 'nx' pixels in horizontal direction (length 'Lx'), idem in vertical direction",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.,
    py::arg("nfine")=1
  )

  .def("coor" , &MD::FineLayer::coor )
  .def("conn" , &MD::FineLayer::conn )
  .def("nelem", &MD::FineLayer::nelem)
  .def("nnode", &MD::FineLayer::nnode)
  .def("nne"  , &MD::FineLayer::nne  )
  .def("ndim" , &MD::FineLayer::ndim )
  .def("nelx" , &MD::FineLayer::nelx )
  .def("nely" , &MD::FineLayer::nely )
  .def("h"    , &MD::FineLayer::h    )

  .def("elementsMiddleLayer", &MD::FineLayer::elementsMiddleLayer)

  .def("nodesBottomEdge"    , &MD::FineLayer::nodesBottomEdge    )
  .def("nodesTopEdge"       , &MD::FineLayer::nodesTopEdge       )
  .def("nodesLeftEdge"      , &MD::FineLayer::nodesLeftEdge      )
  .def("nodesRightEdge"     , &MD::FineLayer::nodesRightEdge     )
  .def("nodesBottomOpenEdge", &MD::FineLayer::nodesBottomOpenEdge)
  .def("nodesTopOpenEdge"   , &MD::FineLayer::nodesTopOpenEdge   )
  .def("nodesLeftOpenEdge"  , &MD::FineLayer::nodesLeftOpenEdge  )
  .def("nodesRightOpenEdge" , &MD::FineLayer::nodesRightOpenEdge )

  .def("nodesBottomLeftCorner" , &MD::FineLayer::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &MD::FineLayer::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &MD::FineLayer::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &MD::FineLayer::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &MD::FineLayer::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &MD::FineLayer::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &MD::FineLayer::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &MD::FineLayer::nodesRightTopCorner   )

  .def("dofs", &MD::FineLayer::dofs)

  .def("nodesPeriodic", &MD::FineLayer::nodesPeriodic)
  .def("nodesOrigin"  , &MD::FineLayer::nodesOrigin  )
  .def("dofsPeriodic" , &MD::FineLayer::dofsPeriodic )

  .def("__repr__",
    [](const MD::FineLayer &){ return "<GooseFEM.Mesh.Quad4.FineLayer>"; }
  );

}

// =================================================================================================

void init_MeshQuad4Map(py::module &m)
{

// -------------------------------------------------------------------------------------------------

py::class_<MD::Map::RefineRegular>(m, "RefineRegular")

  .def(py::init<const MD::Regular &, size_t, size_t>(),
    "Refine a regular mesh",
    py::arg("mesh"),
    py::arg("nx"),
    py::arg("ny")
  )

  .def("getCoarseMesh", &MD::Map::RefineRegular::getCoarseMesh)

  .def("getFineMesh", &MD::Map::RefineRegular::getFineMesh)

  .def("getMap", &MD::Map::RefineRegular::getMap)

  .def("mapToCoarse",
    py::overload_cast<const xt::xtensor<double,1>&>(&MD::Map::RefineRegular::mapToCoarse, py::const_))

  .def("mapToCoarse",
    py::overload_cast<const xt::xtensor<double,2>&>(&MD::Map::RefineRegular::mapToCoarse, py::const_))

  .def("mapToCoarse",
    py::overload_cast<const xt::xtensor<double,4>&>(&MD::Map::RefineRegular::mapToCoarse, py::const_))

  .def("mapToFine",
    py::overload_cast<const xt::xtensor<double,1>&>(&MD::Map::RefineRegular::mapToFine, py::const_))

  .def("mapToFine",
    py::overload_cast<const xt::xtensor<double,2>&>(&MD::Map::RefineRegular::mapToFine, py::const_))

  .def("mapToFine",
    py::overload_cast<const xt::xtensor<double,4>&>(&MD::Map::RefineRegular::mapToFine, py::const_))

  .def("__repr__",
    [](const MD::Map::RefineRegular &){ return "<GooseFEM.Mesh.Quad4.Map.RefineRegular>"; });

// -------------------------------------------------------------------------------------------------

py::class_<MD::Map::FineLayer2Regular>(m, "FineLayer2Regular")

  .def(py::init<const MD::FineLayer &>(),
    "Map a FineLayer mesh to a Regular mesh",
    py::arg("mesh")
  )

  .def("getRegularMesh", &MD::Map::FineLayer2Regular::getRegularMesh)

  .def("getFineLayerMesh", &MD::Map::FineLayer2Regular::getFineLayerMesh)

  .def("getMap", &MD::Map::FineLayer2Regular::getMap)

  .def("getMapFraction", &MD::Map::FineLayer2Regular::getMapFraction)

  .def("mapToRegular",
    py::overload_cast<const xt::xtensor<double,1>&>(&MD::Map::FineLayer2Regular::mapToRegular, py::const_))

  .def("mapToRegular",
    py::overload_cast<const xt::xtensor<double,2>&>(&MD::Map::FineLayer2Regular::mapToRegular, py::const_))

  .def("mapToRegular",
    py::overload_cast<const xt::xtensor<double,4>&>(&MD::Map::FineLayer2Regular::mapToRegular, py::const_))

  .def("__repr__",
    [](const MD::Map::FineLayer2Regular &){ return "<GooseFEM.Mesh.Quad4.Map.FineLayer2Regular>"; });

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================


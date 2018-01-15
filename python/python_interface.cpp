/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

#include "../src/GooseFEM/GooseFEM.h"

// alias for short-hand notation below
namespace py = pybind11;

PYBIND11_MODULE(GooseFEM, m) {

// =================================================================================================

m.doc() = "Some simple finite element meshes and operations";

// define submodules "mXXX"
py::module mMesh       = m    .def_submodule("Mesh"  , "Generic mesh routines"                  );
py::module mMeshTri3   = mMesh.def_submodule("Tri3"  , "Linear triangular elements (2D)"        );
py::module mMeshQuad4  = mMesh.def_submodule("Quad4" , "Linear quadrilateral elements (2D)"     );
py::module mMeshHex8   = mMesh.def_submodule("Hex8"  , "Linear hexahedron (brick) elements (3D)");

// ======================================= GooseFEM/Mesh.h ========================================

mMesh.def("elem2node",
  &GooseFEM::Mesh::elem2node,
  "Elements connect to each node: [ number of elements , element numbers ]",
  py::arg("conn")
);

mMesh.def("dofs",
  &GooseFEM::Mesh::dofs,
  "List with DOF-numbers (in sequential order)",
  py::arg("nnode"),
  py::arg("ndim")
);

mMesh.def("renumber",
  py::overload_cast<const GooseFEM::MatS&>(&GooseFEM::Mesh::renumber),
  "Renumber DOF-list to use the lowest possible index",
  py::arg("dofs")
);

mMesh.def("renumber",
  py::overload_cast<const GooseFEM::MatS&,const GooseFEM::ColS&,std::string>(&GooseFEM::Mesh::renumber),
  "Renumber DOF-list to begin or end with 'idx'",
  py::arg("dofs"),
  py::arg("idx"),
  py::arg("location")="end"
);

// ===================================== GooseFEM/MeshHex8.h =====================================

py::class_<GooseFEM::Mesh::Hex8::Regular>(mMeshHex8, "Regular")

  .def(
    py::init<size_t,size_t,size_t,double>(),
    "mesh with nx*ny*nz 'pixels' and edge size h",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("nz"),
    py::arg("h")=1.
  )

  .def("nelem"                       , &GooseFEM::Mesh::Hex8::Regular::nelem                      )
  .def("nnode"                       , &GooseFEM::Mesh::Hex8::Regular::nnode                      )
  .def("nne"                         , &GooseFEM::Mesh::Hex8::Regular::nne                        )
  .def("ndim"                        , &GooseFEM::Mesh::Hex8::Regular::ndim                       )
  .def("coor"                        , &GooseFEM::Mesh::Hex8::Regular::coor                       )
  .def("conn"                        , &GooseFEM::Mesh::Hex8::Regular::conn                       )
  .def("nodesBottom"                 , &GooseFEM::Mesh::Hex8::Regular::nodesBottom                )
  .def("nodesTop"                    , &GooseFEM::Mesh::Hex8::Regular::nodesTop                   )
  .def("nodesLeft"                   , &GooseFEM::Mesh::Hex8::Regular::nodesLeft                  )
  .def("nodesRight"                  , &GooseFEM::Mesh::Hex8::Regular::nodesRight                 )
  .def("nodesFront"                  , &GooseFEM::Mesh::Hex8::Regular::nodesFront                 )
  .def("nodesBack"                   , &GooseFEM::Mesh::Hex8::Regular::nodesBack                  )
  .def("nodesBottomFace"             , &GooseFEM::Mesh::Hex8::Regular::nodesBottomFace            )
  .def("nodesTopFace"                , &GooseFEM::Mesh::Hex8::Regular::nodesTopFace               )
  .def("nodesLeftFace"               , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFace              )
  .def("nodesRightFace"              , &GooseFEM::Mesh::Hex8::Regular::nodesRightFace             )
  .def("nodesFrontFace"              , &GooseFEM::Mesh::Hex8::Regular::nodesFrontFace             )
  .def("nodesBackFace"               , &GooseFEM::Mesh::Hex8::Regular::nodesBackFace              )
  .def("nodesBottomFrontEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontEdge       )
  .def("nodesBottomBackEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackEdge        )
  .def("nodesBottomLeftEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftEdge        )
  .def("nodesBottomRightEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightEdge       )
  .def("nodesTopFrontEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontEdge          )
  .def("nodesTopBackEdge"            , &GooseFEM::Mesh::Hex8::Regular::nodesTopBackEdge           )
  .def("nodesTopLeftEdge"            , &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftEdge           )
  .def("nodesTopRightEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesTopRightEdge          )
  .def("nodesFrontLeftEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftEdge         )
  .def("nodesFrontRightEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightEdge        )
  .def("nodesBackLeftEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftEdge          )
  .def("nodesBackRightEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesBackRightEdge         )
  .def("nodesFrontBottomEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomEdge       )
  .def("nodesFrontTopEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopEdge          )
  .def("nodesBackBottomEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomEdge        )
  .def("nodesBackTopEdge"            , &GooseFEM::Mesh::Hex8::Regular::nodesBackTopEdge           )
  .def("nodesLeftFrontEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontEdge         )
  .def("nodesLeftBottomEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomEdge        )
  .def("nodesLeftTopEdge"            , &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopEdge           )
  .def("nodesLeftBackEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackEdge          )
  .def("nodesRightFrontEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontEdge        )
  .def("nodesRightBackEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesRightBackEdge         )
  .def("nodesRightBottomEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomEdge       )
  .def("nodesRightTopEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesRightTopEdge          )
  .def("nodesBottomFrontLeftCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontLeftCorner )
  .def("nodesBottomFrontRightCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontRightCorner)
  .def("nodesBottomBackLeftCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackLeftCorner  )
  .def("nodesBottomBackRightCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackRightCorner )
  .def("nodesTopFrontLeftCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontLeftCorner    )
  .def("nodesTopFrontRightCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontRightCorner   )
  .def("nodesTopBackLeftCorner"      , &GooseFEM::Mesh::Hex8::Regular::nodesTopBackLeftCorner     )
  .def("nodesTopBackRightCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesTopBackRightCorner    )
  .def("nodesBottomLeftFrontCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftFrontCorner )
  .def("nodesFrontBottomLeftCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomLeftCorner )
  .def("nodesFrontLeftBottomCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftBottomCorner )
  .def("nodesLeftBottomFrontCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomFrontCorner )
  .def("nodesLeftFrontBottomCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontBottomCorner )
  .def("nodesBottomRightFrontCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightFrontCorner)
  .def("nodesFrontBottomRightCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomRightCorner)
  .def("nodesFrontRightBottomCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightBottomCorner)
  .def("nodesRightBottomFrontCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomFrontCorner)
  .def("nodesRightFrontBottomCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontBottomCorner)
  .def("nodesBottomLeftBackCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftBackCorner  )
  .def("nodesBackBottomLeftCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomLeftCorner  )
  .def("nodesBackLeftBottomCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftBottomCorner  )
  .def("nodesLeftBottomBackCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomBackCorner  )
  .def("nodesLeftBackBottomCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackBottomCorner  )
  .def("nodesBottomRightBackCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightBackCorner )
  .def("nodesBackBottomRightCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomRightCorner )
  .def("nodesBackRightBottomCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBackRightBottomCorner )
  .def("nodesRightBottomBackCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomBackCorner )
  .def("nodesRightBackBottomCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesRightBackBottomCorner )
  .def("nodesTopLeftFrontCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftFrontCorner    )
  .def("nodesFrontTopLeftCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopLeftCorner    )
  .def("nodesFrontLeftTopCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftTopCorner    )
  .def("nodesLeftTopFrontCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopFrontCorner    )
  .def("nodesLeftFrontTopCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontTopCorner    )
  .def("nodesTopRightFrontCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesTopRightFrontCorner   )
  .def("nodesFrontTopRightCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopRightCorner   )
  .def("nodesFrontRightTopCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightTopCorner   )
  .def("nodesRightTopFrontCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesRightTopFrontCorner   )
  .def("nodesRightFrontTopCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontTopCorner   )
  .def("nodesTopLeftBackCorner"      , &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftBackCorner     )
  .def("nodesBackTopLeftCorner"      , &GooseFEM::Mesh::Hex8::Regular::nodesBackTopLeftCorner     )
  .def("nodesBackLeftTopCorner"      , &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftTopCorner     )
  .def("nodesLeftTopBackCorner"      , &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopBackCorner     )
  .def("nodesLeftBackTopCorner"      , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackTopCorner     )
  .def("nodesTopRightBackCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesTopRightBackCorner    )
  .def("nodesBackTopRightCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesBackTopRightCorner    )
  .def("nodesBackRightTopCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesBackRightTopCorner    )
  .def("nodesRightTopBackCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesRightTopBackCorner    )
  .def("nodesRightBackTopCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesRightBackTopCorner    )
  .def("nodesPeriodic"               , &GooseFEM::Mesh::Hex8::Regular::nodesPeriodic              )
  .def("nodesOrigin"                 , &GooseFEM::Mesh::Hex8::Regular::nodesOrigin                )
  .def("dofs"                        , &GooseFEM::Mesh::Hex8::Regular::dofs                       )
  .def("dofsPeriodic"                , &GooseFEM::Mesh::Hex8::Regular::dofsPeriodic               )

  .def("__repr__",
    [](const GooseFEM::Mesh::Hex8::Regular &a){ return "<GooseFEM.Mesh.Hex8.Regular>"; }
  );

// ===================================== GooseFEM/MeshQuad4.h =====================================

py::class_<GooseFEM::Mesh::Quad4::Regular>(mMeshQuad4, "Regular")

  .def(
    py::init<size_t,size_t,double>(),
    "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.
  )

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

  .def("__repr__",
    [](const GooseFEM::Mesh::Quad4::Regular &a){ return "<GooseFEM.Mesh.Quad4.Regular>"; }
  );

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Quad4::FineLayer>(mMeshQuad4, "FineLayer")

  .def(
    py::init<size_t,size_t,double,size_t,size_t>(),
    "FineLayer mesh: 'nx' pixels in horizontal direction (length 'Lx'), idem in vertical direction",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.,
    py::arg("nfine")=1,
    py::arg("nskip")=0
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
    [](const GooseFEM::Mesh::Quad4::FineLayer &a){ return "<GooseFEM.Mesh.Quad4.FineLayer>"; }
  );

// ===================================== GooseFEM/MeshTri3.h ======================================

py::class_<GooseFEM::Mesh::Tri3::Regular>(mMeshTri3, "Regular")

  .def(
    py::init<size_t,size_t,double>(),
    "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.
  )

  .def("coor"                  , &GooseFEM::Mesh::Tri3::Regular::coor                  )
  .def("conn"                  , &GooseFEM::Mesh::Tri3::Regular::conn                  )
  .def("nelem"                 , &GooseFEM::Mesh::Tri3::Regular::nelem                 )
  .def("nnode"                 , &GooseFEM::Mesh::Tri3::Regular::nnode                 )
  .def("nne"                   , &GooseFEM::Mesh::Tri3::Regular::nne                   )
  .def("ndim"                  , &GooseFEM::Mesh::Tri3::Regular::ndim                  )
  .def("nodesBottomEdge"       , &GooseFEM::Mesh::Tri3::Regular::nodesBottomEdge       )
  .def("nodesTopEdge"          , &GooseFEM::Mesh::Tri3::Regular::nodesTopEdge          )
  .def("nodesLeftEdge"         , &GooseFEM::Mesh::Tri3::Regular::nodesLeftEdge         )
  .def("nodesRightEdge"        , &GooseFEM::Mesh::Tri3::Regular::nodesRightEdge        )
  .def("nodesBottomLeftCorner" , &GooseFEM::Mesh::Tri3::Regular::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &GooseFEM::Mesh::Tri3::Regular::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &GooseFEM::Mesh::Tri3::Regular::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &GooseFEM::Mesh::Tri3::Regular::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &GooseFEM::Mesh::Tri3::Regular::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &GooseFEM::Mesh::Tri3::Regular::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &GooseFEM::Mesh::Tri3::Regular::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &GooseFEM::Mesh::Tri3::Regular::nodesRightTopCorner   )
  .def("nodesPeriodic"         , &GooseFEM::Mesh::Tri3::Regular::nodesPeriodic         )
  .def("nodesOrigin"           , &GooseFEM::Mesh::Tri3::Regular::nodesOrigin           )
  .def("dofs"                  , &GooseFEM::Mesh::Tri3::Regular::dofs                  )
  .def("dofsPeriodic"          , &GooseFEM::Mesh::Tri3::Regular::dofsPeriodic          )

  .def("__repr__",
    [](const GooseFEM::Mesh::Tri3::Regular &a){ return "<GooseFEM.Mesh.Tri3.Regular>"; }
  );

// -------------------------------------------------------------------------------------------------

mMeshTri3.def("getOrientation",&GooseFEM::Mesh::Tri3::getOrientation,
  "Get the orientation of each element",
  py::arg("coor"),
  py::arg("conn")
);

// -------------------------------------------------------------------------------------------------

mMeshTri3.def("retriangulate",&GooseFEM::Mesh::Tri3::retriangulate,
  "Re-triangulate existing mesh",
  py::arg("coor"),
  py::arg("conn"),
  py::arg("orientation")=-1
);

// =================================================================================================

}


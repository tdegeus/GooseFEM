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

  .def("nelem"                      , &GooseFEM::Mesh::Hex8::Regular::nelem                      )
  .def("nnode"                      , &GooseFEM::Mesh::Hex8::Regular::nnode                      )
  .def("nne"                        , &GooseFEM::Mesh::Hex8::Regular::nne                        )
  .def("ndim"                       , &GooseFEM::Mesh::Hex8::Regular::ndim                       )
  .def("coor"                       , &GooseFEM::Mesh::Hex8::Regular::coor                       )
  .def("conn"                       , &GooseFEM::Mesh::Hex8::Regular::conn                       )
  .def("nodesFront"                 , &GooseFEM::Mesh::Hex8::Regular::nodesFront                 )
  .def("nodesBack"                  , &GooseFEM::Mesh::Hex8::Regular::nodesBack                  )
  .def("nodesLeft"                  , &GooseFEM::Mesh::Hex8::Regular::nodesLeft                  )
  .def("nodesRight"                 , &GooseFEM::Mesh::Hex8::Regular::nodesRight                 )
  .def("nodesBottom"                , &GooseFEM::Mesh::Hex8::Regular::nodesBottom                )
  .def("nodesTop"                   , &GooseFEM::Mesh::Hex8::Regular::nodesTop                   )
  .def("nodesFrontFace"             , &GooseFEM::Mesh::Hex8::Regular::nodesFrontFace             )
  .def("nodesBackFace"              , &GooseFEM::Mesh::Hex8::Regular::nodesBackFace              )
  .def("nodesLeftFace"              , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFace              )
  .def("nodesRightFace"             , &GooseFEM::Mesh::Hex8::Regular::nodesRightFace             )
  .def("nodesBottomFace"            , &GooseFEM::Mesh::Hex8::Regular::nodesBottomFace            )
  .def("nodesTopFace"               , &GooseFEM::Mesh::Hex8::Regular::nodesTopFace               )
  .def("nodesFrontBottomEdge"       , &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomEdge       )
  .def("nodesFrontTopEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopEdge          )
  .def("nodesFrontLeftEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftEdge         )
  .def("nodesFrontRightEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightEdge        )
  .def("nodesBackBottomEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomEdge        )
  .def("nodesBackTopEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesBackTopEdge           )
  .def("nodesBackLeftEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftEdge          )
  .def("nodesBackRightEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesBackRightEdge         )
  .def("nodesBottomLeftEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftEdge        )
  .def("nodesBottomRightEdge"       , &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightEdge       )
  .def("nodesTopLeftEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftEdge           )
  .def("nodesTopRightEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesTopRightEdge          )
  .def("nodesBottomFrontEdge"       , &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontEdge       )
  .def("nodesBottomBackEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackEdge        )
  .def("nodesTopFrontEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontEdge          )
  .def("nodesTopBackEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesTopBackEdge           )
  .def("nodesLeftBottomEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomEdge        )
  .def("nodesLeftFrontEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontEdge         )
  .def("nodesLeftBackEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackEdge          )
  .def("nodesLeftTopEdge"           , &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopEdge           )
  .def("nodesRightBottomEdge"       , &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomEdge       )
  .def("nodesRightTopEdge"          , &GooseFEM::Mesh::Hex8::Regular::nodesRightTopEdge          )
  .def("nodesRightFrontEdge"        , &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontEdge        )
  .def("nodesRightBackEdge"         , &GooseFEM::Mesh::Hex8::Regular::nodesRightBackEdge         )
  .def("nodesFrontBottomLeftCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomLeftCorner )
  .def("nodesFrontBottomRightCorner", &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomRightCorner)
  .def("nodesFrontTopLeftCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopLeftCorner    )
  .def("nodesFrontTopRightCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopRightCorner   )
  .def("nodesBackBottomLeftCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomLeftCorner  )
  .def("nodesBackBottomRightCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomRightCorner )
  .def("nodesBackTopLeftCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesBackTopLeftCorner     )
  .def("nodesBackTopRightCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesBackTopRightCorner    )
  .def("nodesFrontLeftBottomCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftBottomCorner )
  .def("nodesBottomFrontLeftCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontLeftCorner )
  .def("nodesBottomLeftFrontCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftFrontCorner )
  .def("nodesLeftFrontBottomCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontBottomCorner )
  .def("nodesLeftBottomFrontCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomFrontCorner )
  .def("nodesFrontRightBottomCorner", &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightBottomCorner)
  .def("nodesBottomFrontRightCorner", &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontRightCorner)
  .def("nodesBottomRightFrontCorner", &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightFrontCorner)
  .def("nodesRightFrontBottomCorner", &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontBottomCorner)
  .def("nodesRightBottomFrontCorner", &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomFrontCorner)
  .def("nodesFrontLeftTopCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftTopCorner    )
  .def("nodesTopFrontLeftCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontLeftCorner    )
  .def("nodesTopLeftFrontCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftFrontCorner    )
  .def("nodesLeftFrontTopCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontTopCorner    )
  .def("nodesLeftTopFrontCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopFrontCorner    )
  .def("nodesFrontRightTopCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightTopCorner   )
  .def("nodesTopFrontRightCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontRightCorner   )
  .def("nodesTopRightFrontCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesTopRightFrontCorner   )
  .def("nodesRightFrontTopCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontTopCorner   )
  .def("nodesRightTopFrontCorner"   , &GooseFEM::Mesh::Hex8::Regular::nodesRightTopFrontCorner   )
  .def("nodesBackLeftBottomCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftBottomCorner  )
  .def("nodesBottomBackLeftCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackLeftCorner  )
  .def("nodesBottomLeftBackCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftBackCorner  )
  .def("nodesLeftBackBottomCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackBottomCorner  )
  .def("nodesLeftBottomBackCorner"  , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomBackCorner  )
  .def("nodesBackRightBottomCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesBackRightBottomCorner )
  .def("nodesBottomBackRightCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackRightCorner )
  .def("nodesBottomRightBackCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightBackCorner )
  .def("nodesRightBackBottomCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesRightBackBottomCorner )
  .def("nodesRightBottomBackCorner" , &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomBackCorner )
  .def("nodesBackLeftTopCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftTopCorner     )
  .def("nodesTopBackLeftCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesTopBackLeftCorner     )
  .def("nodesTopLeftBackCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftBackCorner     )
  .def("nodesLeftBackTopCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackTopCorner     )
  .def("nodesLeftTopBackCorner"     , &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopBackCorner     )
  .def("nodesBackRightTopCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesBackRightTopCorner    )
  .def("nodesTopBackRightCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesTopBackRightCorner    )
  .def("nodesTopRightBackCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesTopRightBackCorner    )
  .def("nodesRightBackTopCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesRightBackTopCorner    )
  .def("nodesRightTopBackCorner"    , &GooseFEM::Mesh::Hex8::Regular::nodesRightTopBackCorner    )
  .def("nodesPeriodic"              , &GooseFEM::Mesh::Hex8::Regular::nodesPeriodic              )
  .def("nodesOrigin"                , &GooseFEM::Mesh::Hex8::Regular::nodesOrigin                )
  .def("dofs"                       , &GooseFEM::Mesh::Hex8::Regular::dofs                       )
  .def("dofsPeriodic"               , &GooseFEM::Mesh::Hex8::Regular::dofsPeriodic               )

  .def("__repr__",
    [](const GooseFEM::Mesh::Hex8::Regular &a){ return "<GooseFEM.Mesh.Hex8.Regular>"; }
  );

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Hex8::FineLayer>(mMeshHex8, "FineLayer")

  .def(
    py::init<size_t,size_t,size_t,double,size_t>(),
    "mesh with nx*ny*nz 'pixels' and edge size h",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("nz"),
    py::arg("h")=1.,
    py::arg("nfine")=1
  )

  .def("nelem"                      , &GooseFEM::Mesh::Hex8::FineLayer::nelem                      )
  .def("nnode"                      , &GooseFEM::Mesh::Hex8::FineLayer::nnode                      )
  .def("nne"                        , &GooseFEM::Mesh::Hex8::FineLayer::nne                        )
  .def("ndim"                       , &GooseFEM::Mesh::Hex8::FineLayer::ndim                       )
  .def("shape"                      , &GooseFEM::Mesh::Hex8::FineLayer::shape                      )
  .def("coor"                       , &GooseFEM::Mesh::Hex8::FineLayer::coor                       )
  .def("conn"                       , &GooseFEM::Mesh::Hex8::FineLayer::conn                       )
  // .def("nodesFront"                 , &GooseFEM::Mesh::Hex8::FineLayer::nodesFront                 )
  // .def("nodesBack"                  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBack                  )
  // .def("nodesLeft"                  , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeft                  )
  // .def("nodesRight"                 , &GooseFEM::Mesh::Hex8::FineLayer::nodesRight                 )
  // .def("nodesBottom"                , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottom                )
  // .def("nodesTop"                   , &GooseFEM::Mesh::Hex8::FineLayer::nodesTop                   )
  // .def("nodesFrontFace"             , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontFace             )
  // .def("nodesBackFace"              , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackFace              )
  // .def("nodesLeftFace"              , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFace              )
  // .def("nodesRightFace"             , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFace             )
  // .def("nodesBottomFace"            , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFace            )
  // .def("nodesTopFace"               , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFace               )
  // .def("nodesFrontBottomEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomEdge       )
  // .def("nodesFrontTopEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopEdge          )
  // .def("nodesFrontLeftEdge"         , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftEdge         )
  // .def("nodesFrontRightEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightEdge        )
  // .def("nodesBackBottomEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomEdge        )
  // .def("nodesBackTopEdge"           , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopEdge           )
  // .def("nodesBackLeftEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftEdge          )
  // .def("nodesBackRightEdge"         , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightEdge         )
  // .def("nodesBottomLeftEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftEdge        )
  // .def("nodesBottomRightEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightEdge       )
  // .def("nodesTopLeftEdge"           , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftEdge           )
  // .def("nodesTopRightEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightEdge          )
  // .def("nodesBottomFrontEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontEdge       )
  // .def("nodesBottomBackEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackEdge        )
  // .def("nodesTopFrontEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontEdge          )
  // .def("nodesTopBackEdge"           , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackEdge           )
  // .def("nodesLeftBottomEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomEdge        )
  // .def("nodesLeftFrontEdge"         , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontEdge         )
  // .def("nodesLeftBackEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackEdge          )
  // .def("nodesLeftTopEdge"           , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopEdge           )
  // .def("nodesRightBottomEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomEdge       )
  // .def("nodesRightTopEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopEdge          )
  // .def("nodesRightFrontEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontEdge        )
  // .def("nodesRightBackEdge"         , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackEdge         )
  // .def("nodesFrontBottomLeftCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomLeftCorner )
  // .def("nodesFrontBottomRightCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomRightCorner)
  // .def("nodesFrontTopLeftCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopLeftCorner    )
  // .def("nodesFrontTopRightCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopRightCorner   )
  // .def("nodesBackBottomLeftCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomLeftCorner  )
  // .def("nodesBackBottomRightCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomRightCorner )
  // .def("nodesBackTopLeftCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopLeftCorner     )
  // .def("nodesBackTopRightCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopRightCorner    )
  // .def("nodesFrontLeftBottomCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftBottomCorner )
  // .def("nodesBottomFrontLeftCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontLeftCorner )
  // .def("nodesBottomLeftFrontCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftFrontCorner )
  // .def("nodesLeftFrontBottomCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontBottomCorner )
  // .def("nodesLeftBottomFrontCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomFrontCorner )
  // .def("nodesFrontRightBottomCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightBottomCorner)
  // .def("nodesBottomFrontRightCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontRightCorner)
  // .def("nodesBottomRightFrontCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightFrontCorner)
  // .def("nodesRightFrontBottomCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontBottomCorner)
  // .def("nodesRightBottomFrontCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomFrontCorner)
  // .def("nodesFrontLeftTopCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftTopCorner    )
  // .def("nodesTopFrontLeftCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontLeftCorner    )
  // .def("nodesTopLeftFrontCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftFrontCorner    )
  // .def("nodesLeftFrontTopCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontTopCorner    )
  // .def("nodesLeftTopFrontCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopFrontCorner    )
  // .def("nodesFrontRightTopCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightTopCorner   )
  // .def("nodesTopFrontRightCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontRightCorner   )
  // .def("nodesTopRightFrontCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightFrontCorner   )
  // .def("nodesRightFrontTopCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontTopCorner   )
  // .def("nodesRightTopFrontCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopFrontCorner   )
  // .def("nodesBackLeftBottomCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftBottomCorner  )
  // .def("nodesBottomBackLeftCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackLeftCorner  )
  // .def("nodesBottomLeftBackCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftBackCorner  )
  // .def("nodesLeftBackBottomCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackBottomCorner  )
  // .def("nodesLeftBottomBackCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomBackCorner  )
  // .def("nodesBackRightBottomCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightBottomCorner )
  // .def("nodesBottomBackRightCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackRightCorner )
  // .def("nodesBottomRightBackCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightBackCorner )
  // .def("nodesRightBackBottomCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackBottomCorner )
  // .def("nodesRightBottomBackCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomBackCorner )
  // .def("nodesBackLeftTopCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftTopCorner     )
  // .def("nodesTopBackLeftCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackLeftCorner     )
  // .def("nodesTopLeftBackCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftBackCorner     )
  // .def("nodesLeftBackTopCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackTopCorner     )
  // .def("nodesLeftTopBackCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopBackCorner     )
  // .def("nodesBackRightTopCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightTopCorner    )
  // .def("nodesTopBackRightCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackRightCorner    )
  // .def("nodesTopRightBackCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightBackCorner    )
  // .def("nodesRightBackTopCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackTopCorner    )
  // .def("nodesRightTopBackCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopBackCorner    )
  // .def("nodesPeriodic"              , &GooseFEM::Mesh::Hex8::FineLayer::nodesPeriodic              )
  // .def("nodesOrigin"                , &GooseFEM::Mesh::Hex8::FineLayer::nodesOrigin                )
  // .def("dofs"                       , &GooseFEM::Mesh::Hex8::FineLayer::dofs                       )
  // .def("dofsPeriodic"               , &GooseFEM::Mesh::Hex8::FineLayer::dofsPeriodic               )

  .def("__repr__",
    [](const GooseFEM::Mesh::Hex8::FineLayer &a){ return "<GooseFEM.Mesh.Hex8.FineLayer>"; }
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


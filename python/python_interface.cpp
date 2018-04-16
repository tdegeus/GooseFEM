/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cppmat/pybind11.h>

#include "../src/GooseFEM/GooseFEM.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;
namespace M  = GooseFEM;

// abbreviate types(s)
typedef GooseFEM::ColD ColD;
typedef GooseFEM::ColS ColS;
typedef GooseFEM::MatD MatD;
typedef GooseFEM::MatS MatS;
typedef GooseFEM::ArrD ArrD;

typedef const GooseFEM::ColD cColD;
typedef const GooseFEM::ColS cColS;
typedef const GooseFEM::MatD cMatD;
typedef const GooseFEM::MatS cMatS;
typedef const GooseFEM::ArrD cArrD;

// =========================================== GooseFEM ============================================

PYBIND11_MODULE(GooseFEM, m) {

m.doc() = "Some simple finite element meshes and operations";

// ================================= GooseFEM - GooseFEM/Vector.h ==================================

py::class_<GooseFEM::Vector>(m, "Vector")
  // constructor
  .def(
    py::init<cMatS &, cMatS &, cColS &>(),
    "Class to switch between DOF/nodal/element views of vectors",
    py::arg("conn"),
    py::arg("dofs"),
    py::arg("iip")=GooseFEM::ColS()
  )
  // methods
  .def("nelem", &M::Vector::nelem)
  .def("nne"  , &M::Vector::nne  )
  .def("nnode", &M::Vector::nnode)
  .def("ndim" , &M::Vector::ndim )
  .def("ndof" , &M::Vector::ndof )
  .def("nnu"  , &M::Vector::nnu  )
  .def("nnp"  , &M::Vector::nnp  )
  // -
  .def("iiu"  , &M::Vector::iiu  )
  .def("iip"  , &M::Vector::iip  )
  // -
  .def("asDofs"    , py::overload_cast<cColD&,cColD&>(&M::Vector::asDofs    , py::const_))
  .def("asDofs"    , py::overload_cast<cMatD&       >(&M::Vector::asDofs    , py::const_))
  .def("asDofs"    , py::overload_cast<cArrD&       >(&M::Vector::asDofs    , py::const_))
  .def("asDofs_u"  , py::overload_cast<cMatD&       >(&M::Vector::asDofs_u  , py::const_))
  .def("asDofs_u"  , py::overload_cast<cArrD&       >(&M::Vector::asDofs_u  , py::const_))
  .def("asDofs_p"  , py::overload_cast<cMatD&       >(&M::Vector::asDofs_p  , py::const_))
  .def("asDofs_p"  , py::overload_cast<cArrD&       >(&M::Vector::asDofs_p  , py::const_))
  .def("asNode"    , py::overload_cast<cColD&       >(&M::Vector::asNode    , py::const_))
  .def("asNode"    , py::overload_cast<cColD&,cColD&>(&M::Vector::asNode    , py::const_))
  .def("asNode"    , py::overload_cast<cArrD&       >(&M::Vector::asNode    , py::const_))
  .def("asElement" , py::overload_cast<cColD&       >(&M::Vector::asElement , py::const_))
  .def("asElement" , py::overload_cast<cColD&,cColD&>(&M::Vector::asElement , py::const_))
  .def("asElement" , py::overload_cast<cMatD&       >(&M::Vector::asElement , py::const_))
  // -
  .def("assembleDofs"  , py::overload_cast<cMatD&>(&M::Vector::assembleDofs  , py::const_))
  .def("assembleDofs"  , py::overload_cast<cArrD&>(&M::Vector::assembleDofs  , py::const_))
  .def("assembleDofs_u", py::overload_cast<cMatD&>(&M::Vector::assembleDofs_u, py::const_))
  .def("assembleDofs_u", py::overload_cast<cArrD&>(&M::Vector::assembleDofs_u, py::const_))
  .def("assembleDofs_p", py::overload_cast<cMatD&>(&M::Vector::assembleDofs_p, py::const_))
  .def("assembleDofs_p", py::overload_cast<cArrD&>(&M::Vector::assembleDofs_p, py::const_))
  .def("assembleNode"  , py::overload_cast<cArrD&>(&M::Vector::assembleNode  , py::const_))
  // print to screen
  .def("__repr__",
    [](const GooseFEM::Vector &a){ return "<GooseFEM.Vector>"; }
  );

// ============================= GooseFEM - GooseFEM/MatrixDiagonal.h ==============================

py::class_<GooseFEM::DiagonalMatrix>(m, "DiagonalMatrix")
  // constructor
  .def(
    py::init<cMatS &, cMatS &, cColS &>(),
    "Class to switch between DOF/nodal/element views of vectors",
    py::arg("conn"),
    py::arg("dofs"),
    py::arg("iip")=GooseFEM::ColS()
  )
  // methods
  .def("nelem", &M::DiagonalMatrix::nelem)
  .def("nne"  , &M::DiagonalMatrix::nne  )
  .def("nnode", &M::DiagonalMatrix::nnode)
  .def("ndim" , &M::DiagonalMatrix::ndim )
  .def("ndof" , &M::DiagonalMatrix::ndof )
  .def("nnu"  , &M::DiagonalMatrix::nnu  )
  .def("nnp"  , &M::DiagonalMatrix::nnp  )
  // -
  .def("iiu"  , &M::DiagonalMatrix::iiu  )
  .def("iip"  , &M::DiagonalMatrix::iip  )
  // -
  .def("dot"  ,                                  &M::DiagonalMatrix::dot               )
  .def("dot_u", py::overload_cast<cColD&       >(&M::DiagonalMatrix::dot_u, py::const_))
  .def("dot_u", py::overload_cast<cColD&,cColD&>(&M::DiagonalMatrix::dot_u, py::const_))
  .def("dot_p", py::overload_cast<cColD&       >(&M::DiagonalMatrix::dot_p, py::const_))
  .def("dot_p", py::overload_cast<cColD&,cColD&>(&M::DiagonalMatrix::dot_p, py::const_))
  // -
  .def("check_diagonal", &M::DiagonalMatrix::check_diagonal)
  .def("assemble"      , &M::DiagonalMatrix::assemble      )
  // .def("set"           , &M::DiagonalMatrix::set           )
  // .def("set_uu"        , &M::DiagonalMatrix::set_uu        )
  // .def("set_pp"        , &M::DiagonalMatrix::set_pp        )
  .def("solve"         , &M::DiagonalMatrix::solve         )
  .def("solve_u"       , &M::DiagonalMatrix::solve_u       )
  // .def("rhs_p"         , &M::DiagonalMatrix::rhs_p         )
  .def("asDiagonal"    , &M::DiagonalMatrix::asDiagonal    )
  // .def("asDiagonal_uu" , &M::DiagonalMatrix::asDiagonal_uu )
  // .def("asDiagonal_pp" , &M::DiagonalMatrix::asDiagonal_pp )
  // .def("asSparse"      , &M::DiagonalMatrix::asSparse      )
  // .def("asSparse_uu"   , &M::DiagonalMatrix::asSparse_uu   )
  // .def("asSparse_up"   , &M::DiagonalMatrix::asSparse_up   )
  // .def("asSparse_pu"   , &M::DiagonalMatrix::asSparse_pu   )
  // .def("asSparse_pp"   , &M::DiagonalMatrix::asSparse_pp   )
  // .def("asDense"       , &M::DiagonalMatrix::asDense       )
  // .def("asDense_uu"    , &M::DiagonalMatrix::asDense_uu    )
  // .def("asDense_up"    , &M::DiagonalMatrix::asDense_up    )
  // .def("asDense_pu"    , &M::DiagonalMatrix::asDense_pu    )
  // .def("asDense_pp"    , &M::DiagonalMatrix::asDense_pp    )

  // print to screen
  .def("__repr__",
    [](const GooseFEM::DiagonalMatrix &a){ return "<GooseFEM.DiagonalMatrix>"; }
  );

// ============================ GooseFEM::Element - GooseFEM/Element.h =============================

py::module mElement = m.def_submodule("Element", "Generic element routines");

// -------------------------------------------------------------------------------------------------

mElement.def("asElementVector",
  &GooseFEM::Element::asElementVector,
  "convert nodal vector [nnode, ndim] to nodal vector stored per element [nelem, nne, ndim]",
  py::arg("conn"),
  py::arg("nodevec")
);

// -------------------------------------------------------------------------------------------------

mElement.def("assembleElementVector",
  &GooseFEM::Element::assembleNodeVector,
  "assemble nodal vector stored per element [nelem, nne, ndim] to nodal vector [nnode, ndim]",
  py::arg("conn"),
  py::arg("elemvec")
);

// ====================== GooseFEM::Element::Quad4 - GooseFEM/ElementQuad4.h =======================

{

// create submodule
py::module sm = mElement.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

// abbreviate name-space
namespace SM = GooseFEM::Element::Quad4;

// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------

}

// =============================== GooseFEM::Mesh - GooseFEM/Mesh.h ================================

py::module mMesh = m.def_submodule("Mesh", "Generic mesh routines");

// -------------------------------------------------------------------------------------------------

mMesh.def("elem2node",
  &GooseFEM::Mesh::elem2node,
  "Elements connect to each node: [ number of elements , element numbers ]",
  py::arg("conn")
);

// -------------------------------------------------------------------------------------------------

mMesh.def("dofs",
  &GooseFEM::Mesh::dofs,
  "List with DOF-numbers (in sequential order)",
  py::arg("nnode"),
  py::arg("ndim")
);

// -------------------------------------------------------------------------------------------------

using renumber = GooseFEM::MatS(cMatS &);

mMesh.def("renumber",
  py::overload_cast<cMatS &>((renumber*)&GooseFEM::Mesh::renumber),
  "Renumber DOF-list to use the lowest possible index",
  py::arg("dofs")
);

// -------------------------------------------------------------------------------------------------

using reorder = GooseFEM::MatS(cMatS &, cColS&, std::string);

mMesh.def("reorder",
  py::overload_cast<cMatS&,cColS&,std::string>(
    (reorder*)&GooseFEM::Mesh::reorder),
  "Renumber DOF-list to begin or end with 'idx'",
  py::arg("dofs"),
  py::arg("idx"),
  py::arg("location")="end"
);

// ========================== GooseFEM::Mesh::Hex8 - GooseFEM/MeshHex8.h ===========================

{

// create submodule
py::module sm = mMesh.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

// abbreviate name-space
namespace SM = GooseFEM::Mesh::Hex8;

// -------------------------------------------------------------------------------------------------

py::class_<SM::Regular>(sm, "Regular")
  // constructor
  .def(
    py::init<size_t,size_t,size_t,double>(),
    "mesh with nx*ny*nz 'pixels' and edge size h",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("nz"),
    py::arg("h")=1.
  )
  // sizes
  .def("nelem"                      , &SM::Regular::nelem                      )
  .def("nnode"                      , &SM::Regular::nnode                      )
  .def("nne"                        , &SM::Regular::nne                        )
  .def("ndim"                       , &SM::Regular::ndim                       )
  // mesh
  .def("coor"                       , &SM::Regular::coor                       )
  .def("conn"                       , &SM::Regular::conn                       )
  // boundary nodes: planes
  .def("nodesFront"                 , &SM::Regular::nodesFront                 )
  .def("nodesBack"                  , &SM::Regular::nodesBack                  )
  .def("nodesLeft"                  , &SM::Regular::nodesLeft                  )
  .def("nodesRight"                 , &SM::Regular::nodesRight                 )
  .def("nodesBottom"                , &SM::Regular::nodesBottom                )
  .def("nodesTop"                   , &SM::Regular::nodesTop                   )
  // boundary nodes: faces
  .def("nodesFrontFace"             , &SM::Regular::nodesFrontFace             )
  .def("nodesBackFace"              , &SM::Regular::nodesBackFace              )
  .def("nodesLeftFace"              , &SM::Regular::nodesLeftFace              )
  .def("nodesRightFace"             , &SM::Regular::nodesRightFace             )
  .def("nodesBottomFace"            , &SM::Regular::nodesBottomFace            )
  .def("nodesTopFace"               , &SM::Regular::nodesTopFace               )
  // boundary nodes: edges
  .def("nodesFrontBottomEdge"       , &SM::Regular::nodesFrontBottomEdge       )
  .def("nodesFrontTopEdge"          , &SM::Regular::nodesFrontTopEdge          )
  .def("nodesFrontLeftEdge"         , &SM::Regular::nodesFrontLeftEdge         )
  .def("nodesFrontRightEdge"        , &SM::Regular::nodesFrontRightEdge        )
  .def("nodesBackBottomEdge"        , &SM::Regular::nodesBackBottomEdge        )
  .def("nodesBackTopEdge"           , &SM::Regular::nodesBackTopEdge           )
  .def("nodesBackLeftEdge"          , &SM::Regular::nodesBackLeftEdge          )
  .def("nodesBackRightEdge"         , &SM::Regular::nodesBackRightEdge         )
  .def("nodesBottomLeftEdge"        , &SM::Regular::nodesBottomLeftEdge        )
  .def("nodesBottomRightEdge"       , &SM::Regular::nodesBottomRightEdge       )
  .def("nodesTopLeftEdge"           , &SM::Regular::nodesTopLeftEdge           )
  .def("nodesTopRightEdge"          , &SM::Regular::nodesTopRightEdge          )
  // boundary nodes: faces (aliases)
  .def("nodesBottomFrontEdge"       , &SM::Regular::nodesBottomFrontEdge       )
  .def("nodesBottomBackEdge"        , &SM::Regular::nodesBottomBackEdge        )
  .def("nodesTopFrontEdge"          , &SM::Regular::nodesTopFrontEdge          )
  .def("nodesTopBackEdge"           , &SM::Regular::nodesTopBackEdge           )
  .def("nodesLeftBottomEdge"        , &SM::Regular::nodesLeftBottomEdge        )
  .def("nodesLeftFrontEdge"         , &SM::Regular::nodesLeftFrontEdge         )
  .def("nodesLeftBackEdge"          , &SM::Regular::nodesLeftBackEdge          )
  .def("nodesLeftTopEdge"           , &SM::Regular::nodesLeftTopEdge           )
  .def("nodesRightBottomEdge"       , &SM::Regular::nodesRightBottomEdge       )
  .def("nodesRightTopEdge"          , &SM::Regular::nodesRightTopEdge          )
  .def("nodesRightFrontEdge"        , &SM::Regular::nodesRightFrontEdge        )
  .def("nodesRightBackEdge"         , &SM::Regular::nodesRightBackEdge         )
  // boundary nodes: edges, without corners
  .def("nodesFrontBottomOpenEdge"   , &SM::Regular::nodesFrontBottomOpenEdge   )
  .def("nodesFrontTopOpenEdge"      , &SM::Regular::nodesFrontTopOpenEdge      )
  .def("nodesFrontLeftOpenEdge"     , &SM::Regular::nodesFrontLeftOpenEdge     )
  .def("nodesFrontRightOpenEdge"    , &SM::Regular::nodesFrontRightOpenEdge    )
  .def("nodesBackBottomOpenEdge"    , &SM::Regular::nodesBackBottomOpenEdge    )
  .def("nodesBackTopOpenEdge"       , &SM::Regular::nodesBackTopOpenEdge       )
  .def("nodesBackLeftOpenEdge"      , &SM::Regular::nodesBackLeftOpenEdge      )
  .def("nodesBackRightOpenEdge"     , &SM::Regular::nodesBackRightOpenEdge     )
  .def("nodesBottomLeftOpenEdge"    , &SM::Regular::nodesBottomLeftOpenEdge    )
  .def("nodesBottomRightOpenEdge"   , &SM::Regular::nodesBottomRightOpenEdge   )
  .def("nodesTopLeftOpenEdge"       , &SM::Regular::nodesTopLeftOpenEdge       )
  .def("nodesTopRightOpenEdge"      , &SM::Regular::nodesTopRightOpenEdge      )
  // boundary nodes: edges, without corners (aliases)
  .def("nodesBottomFrontOpenEdge"   , &SM::Regular::nodesBottomFrontOpenEdge   )
  .def("nodesBottomBackOpenEdge"    , &SM::Regular::nodesBottomBackOpenEdge    )
  .def("nodesTopFrontOpenEdge"      , &SM::Regular::nodesTopFrontOpenEdge      )
  .def("nodesTopBackOpenEdge"       , &SM::Regular::nodesTopBackOpenEdge       )
  .def("nodesLeftBottomOpenEdge"    , &SM::Regular::nodesLeftBottomOpenEdge    )
  .def("nodesLeftFrontOpenEdge"     , &SM::Regular::nodesLeftFrontOpenEdge     )
  .def("nodesLeftBackOpenEdge"      , &SM::Regular::nodesLeftBackOpenEdge      )
  .def("nodesLeftTopOpenEdge"       , &SM::Regular::nodesLeftTopOpenEdge       )
  .def("nodesRightBottomOpenEdge"   , &SM::Regular::nodesRightBottomOpenEdge   )
  .def("nodesRightTopOpenEdge"      , &SM::Regular::nodesRightTopOpenEdge      )
  .def("nodesRightFrontOpenEdge"    , &SM::Regular::nodesRightFrontOpenEdge    )
  .def("nodesRightBackOpenEdge"     , &SM::Regular::nodesRightBackOpenEdge     )
  // boundary nodes: corners
  .def("nodesFrontBottomLeftCorner" , &SM::Regular::nodesFrontBottomLeftCorner )
  .def("nodesFrontBottomRightCorner", &SM::Regular::nodesFrontBottomRightCorner)
  .def("nodesFrontTopLeftCorner"    , &SM::Regular::nodesFrontTopLeftCorner    )
  .def("nodesFrontTopRightCorner"   , &SM::Regular::nodesFrontTopRightCorner   )
  .def("nodesBackBottomLeftCorner"  , &SM::Regular::nodesBackBottomLeftCorner  )
  .def("nodesBackBottomRightCorner" , &SM::Regular::nodesBackBottomRightCorner )
  .def("nodesBackTopLeftCorner"     , &SM::Regular::nodesBackTopLeftCorner     )
  .def("nodesBackTopRightCorner"    , &SM::Regular::nodesBackTopRightCorner    )
  // boundary nodes: corners (aliases)
  .def("nodesFrontLeftBottomCorner" , &SM::Regular::nodesFrontLeftBottomCorner )
  .def("nodesBottomFrontLeftCorner" , &SM::Regular::nodesBottomFrontLeftCorner )
  .def("nodesBottomLeftFrontCorner" , &SM::Regular::nodesBottomLeftFrontCorner )
  .def("nodesLeftFrontBottomCorner" , &SM::Regular::nodesLeftFrontBottomCorner )
  .def("nodesLeftBottomFrontCorner" , &SM::Regular::nodesLeftBottomFrontCorner )
  .def("nodesFrontRightBottomCorner", &SM::Regular::nodesFrontRightBottomCorner)
  .def("nodesBottomFrontRightCorner", &SM::Regular::nodesBottomFrontRightCorner)
  .def("nodesBottomRightFrontCorner", &SM::Regular::nodesBottomRightFrontCorner)
  .def("nodesRightFrontBottomCorner", &SM::Regular::nodesRightFrontBottomCorner)
  .def("nodesRightBottomFrontCorner", &SM::Regular::nodesRightBottomFrontCorner)
  .def("nodesFrontLeftTopCorner"    , &SM::Regular::nodesFrontLeftTopCorner    )
  .def("nodesTopFrontLeftCorner"    , &SM::Regular::nodesTopFrontLeftCorner    )
  .def("nodesTopLeftFrontCorner"    , &SM::Regular::nodesTopLeftFrontCorner    )
  .def("nodesLeftFrontTopCorner"    , &SM::Regular::nodesLeftFrontTopCorner    )
  .def("nodesLeftTopFrontCorner"    , &SM::Regular::nodesLeftTopFrontCorner    )
  .def("nodesFrontRightTopCorner"   , &SM::Regular::nodesFrontRightTopCorner   )
  .def("nodesTopFrontRightCorner"   , &SM::Regular::nodesTopFrontRightCorner   )
  .def("nodesTopRightFrontCorner"   , &SM::Regular::nodesTopRightFrontCorner   )
  .def("nodesRightFrontTopCorner"   , &SM::Regular::nodesRightFrontTopCorner   )
  .def("nodesRightTopFrontCorner"   , &SM::Regular::nodesRightTopFrontCorner   )
  .def("nodesBackLeftBottomCorner"  , &SM::Regular::nodesBackLeftBottomCorner  )
  .def("nodesBottomBackLeftCorner"  , &SM::Regular::nodesBottomBackLeftCorner  )
  .def("nodesBottomLeftBackCorner"  , &SM::Regular::nodesBottomLeftBackCorner  )
  .def("nodesLeftBackBottomCorner"  , &SM::Regular::nodesLeftBackBottomCorner  )
  .def("nodesLeftBottomBackCorner"  , &SM::Regular::nodesLeftBottomBackCorner  )
  .def("nodesBackRightBottomCorner" , &SM::Regular::nodesBackRightBottomCorner )
  .def("nodesBottomBackRightCorner" , &SM::Regular::nodesBottomBackRightCorner )
  .def("nodesBottomRightBackCorner" , &SM::Regular::nodesBottomRightBackCorner )
  .def("nodesRightBackBottomCorner" , &SM::Regular::nodesRightBackBottomCorner )
  .def("nodesRightBottomBackCorner" , &SM::Regular::nodesRightBottomBackCorner )
  .def("nodesBackLeftTopCorner"     , &SM::Regular::nodesBackLeftTopCorner     )
  .def("nodesTopBackLeftCorner"     , &SM::Regular::nodesTopBackLeftCorner     )
  .def("nodesTopLeftBackCorner"     , &SM::Regular::nodesTopLeftBackCorner     )
  .def("nodesLeftBackTopCorner"     , &SM::Regular::nodesLeftBackTopCorner     )
  .def("nodesLeftTopBackCorner"     , &SM::Regular::nodesLeftTopBackCorner     )
  .def("nodesBackRightTopCorner"    , &SM::Regular::nodesBackRightTopCorner    )
  .def("nodesTopBackRightCorner"    , &SM::Regular::nodesTopBackRightCorner    )
  .def("nodesTopRightBackCorner"    , &SM::Regular::nodesTopRightBackCorner    )
  .def("nodesRightBackTopCorner"    , &SM::Regular::nodesRightBackTopCorner    )
  .def("nodesRightTopBackCorner"    , &SM::Regular::nodesRightTopBackCorner    )
  // periodicity
  .def("nodesPeriodic"              , &SM::Regular::nodesPeriodic              )
  .def("nodesOrigin"                , &SM::Regular::nodesOrigin                )
  .def("dofs"                       , &SM::Regular::dofs                       )
  .def("dofsPeriodic"               , &SM::Regular::dofsPeriodic               )
  // print to screen
  .def("__repr__",
    [](const SM::Regular &a){ return "<GooseFEM.Mesh.Hex8.Regular>"; }
  );

// -------------------------------------------------------------------------------------------------

py::class_<SM::FineLayer>(sm, "FineLayer")
  // constructor
  .def(
    py::init<size_t,size_t,size_t,double,size_t>(),
    "mesh with nx*ny*nz 'pixels' and edge size h",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("nz"),
    py::arg("h")=1.,
    py::arg("nfine")=1
  )
  // sizes
  .def("nelem"                      , &SM::FineLayer::nelem                      )
  .def("nnode"                      , &SM::FineLayer::nnode                      )
  .def("nne"                        , &SM::FineLayer::nne                        )
  .def("ndim"                       , &SM::FineLayer::ndim                       )
  .def("shape"                      , &SM::FineLayer::shape                      )
  // mesh
  .def("coor"                       , &SM::FineLayer::coor                       )
  .def("conn"                       , &SM::FineLayer::conn                       )
  // element sets
  .def("elementsMiddleLayer"        , &SM::FineLayer::elementsMiddleLayer        )
  // boundary nodes: planes
  .def("nodesFront"                 , &SM::FineLayer::nodesFront                 )
  .def("nodesBack"                  , &SM::FineLayer::nodesBack                  )
  .def("nodesLeft"                  , &SM::FineLayer::nodesLeft                  )
  .def("nodesRight"                 , &SM::FineLayer::nodesRight                 )
  .def("nodesBottom"                , &SM::FineLayer::nodesBottom                )
  .def("nodesTop"                   , &SM::FineLayer::nodesTop                   )
  // boundary nodes: faces
  .def("nodesFrontFace"             , &SM::FineLayer::nodesFrontFace             )
  .def("nodesBackFace"              , &SM::FineLayer::nodesBackFace              )
  .def("nodesLeftFace"              , &SM::FineLayer::nodesLeftFace              )
  .def("nodesRightFace"             , &SM::FineLayer::nodesRightFace             )
  .def("nodesBottomFace"            , &SM::FineLayer::nodesBottomFace            )
  .def("nodesTopFace"               , &SM::FineLayer::nodesTopFace               )
  // boundary nodes: edges
  .def("nodesFrontBottomEdge"       , &SM::FineLayer::nodesFrontBottomEdge       )
  .def("nodesFrontTopEdge"          , &SM::FineLayer::nodesFrontTopEdge          )
  .def("nodesFrontLeftEdge"         , &SM::FineLayer::nodesFrontLeftEdge         )
  .def("nodesFrontRightEdge"        , &SM::FineLayer::nodesFrontRightEdge        )
  .def("nodesBackBottomEdge"        , &SM::FineLayer::nodesBackBottomEdge        )
  .def("nodesBackTopEdge"           , &SM::FineLayer::nodesBackTopEdge           )
  .def("nodesBackLeftEdge"          , &SM::FineLayer::nodesBackLeftEdge          )
  .def("nodesBackRightEdge"         , &SM::FineLayer::nodesBackRightEdge         )
  .def("nodesBottomLeftEdge"        , &SM::FineLayer::nodesBottomLeftEdge        )
  .def("nodesBottomRightEdge"       , &SM::FineLayer::nodesBottomRightEdge       )
  .def("nodesTopLeftEdge"           , &SM::FineLayer::nodesTopLeftEdge           )
  .def("nodesTopRightEdge"          , &SM::FineLayer::nodesTopRightEdge          )
  // boundary nodes: faces (aliases)
  .def("nodesBottomFrontEdge"       , &SM::FineLayer::nodesBottomFrontEdge       )
  .def("nodesBottomBackEdge"        , &SM::FineLayer::nodesBottomBackEdge        )
  .def("nodesTopFrontEdge"          , &SM::FineLayer::nodesTopFrontEdge          )
  .def("nodesTopBackEdge"           , &SM::FineLayer::nodesTopBackEdge           )
  .def("nodesLeftBottomEdge"        , &SM::FineLayer::nodesLeftBottomEdge        )
  .def("nodesLeftFrontEdge"         , &SM::FineLayer::nodesLeftFrontEdge         )
  .def("nodesLeftBackEdge"          , &SM::FineLayer::nodesLeftBackEdge          )
  .def("nodesLeftTopEdge"           , &SM::FineLayer::nodesLeftTopEdge           )
  .def("nodesRightBottomEdge"       , &SM::FineLayer::nodesRightBottomEdge       )
  .def("nodesRightTopEdge"          , &SM::FineLayer::nodesRightTopEdge          )
  .def("nodesRightFrontEdge"        , &SM::FineLayer::nodesRightFrontEdge        )
  .def("nodesRightBackEdge"         , &SM::FineLayer::nodesRightBackEdge         )
  // boundary nodes: edges, without corners
  .def("nodesFrontBottomOpenEdge"   , &SM::FineLayer::nodesFrontBottomOpenEdge   )
  .def("nodesFrontTopOpenEdge"      , &SM::FineLayer::nodesFrontTopOpenEdge      )
  .def("nodesFrontLeftOpenEdge"     , &SM::FineLayer::nodesFrontLeftOpenEdge     )
  .def("nodesFrontRightOpenEdge"    , &SM::FineLayer::nodesFrontRightOpenEdge    )
  .def("nodesBackBottomOpenEdge"    , &SM::FineLayer::nodesBackBottomOpenEdge    )
  .def("nodesBackTopOpenEdge"       , &SM::FineLayer::nodesBackTopOpenEdge       )
  .def("nodesBackLeftOpenEdge"      , &SM::FineLayer::nodesBackLeftOpenEdge      )
  .def("nodesBackRightOpenEdge"     , &SM::FineLayer::nodesBackRightOpenEdge     )
  .def("nodesBottomLeftOpenEdge"    , &SM::FineLayer::nodesBottomLeftOpenEdge    )
  .def("nodesBottomRightOpenEdge"   , &SM::FineLayer::nodesBottomRightOpenEdge   )
  .def("nodesTopLeftOpenEdge"       , &SM::FineLayer::nodesTopLeftOpenEdge       )
  .def("nodesTopRightOpenEdge"      , &SM::FineLayer::nodesTopRightOpenEdge      )
  // boundary nodes: edges, without corners (aliases)
  .def("nodesBottomFrontOpenEdge"   , &SM::FineLayer::nodesBottomFrontOpenEdge   )
  .def("nodesBottomBackOpenEdge"    , &SM::FineLayer::nodesBottomBackOpenEdge    )
  .def("nodesTopFrontOpenEdge"      , &SM::FineLayer::nodesTopFrontOpenEdge      )
  .def("nodesTopBackOpenEdge"       , &SM::FineLayer::nodesTopBackOpenEdge       )
  .def("nodesLeftBottomOpenEdge"    , &SM::FineLayer::nodesLeftBottomOpenEdge    )
  .def("nodesLeftFrontOpenEdge"     , &SM::FineLayer::nodesLeftFrontOpenEdge     )
  .def("nodesLeftBackOpenEdge"      , &SM::FineLayer::nodesLeftBackOpenEdge      )
  .def("nodesLeftTopOpenEdge"       , &SM::FineLayer::nodesLeftTopOpenEdge       )
  .def("nodesRightBottomOpenEdge"   , &SM::FineLayer::nodesRightBottomOpenEdge   )
  .def("nodesRightTopOpenEdge"      , &SM::FineLayer::nodesRightTopOpenEdge      )
  .def("nodesRightFrontOpenEdge"    , &SM::FineLayer::nodesRightFrontOpenEdge    )
  .def("nodesRightBackOpenEdge"     , &SM::FineLayer::nodesRightBackOpenEdge     )
  // boundary nodes: corners
  .def("nodesFrontBottomLeftCorner" , &SM::FineLayer::nodesFrontBottomLeftCorner )
  .def("nodesFrontBottomRightCorner", &SM::FineLayer::nodesFrontBottomRightCorner)
  .def("nodesFrontTopLeftCorner"    , &SM::FineLayer::nodesFrontTopLeftCorner    )
  .def("nodesFrontTopRightCorner"   , &SM::FineLayer::nodesFrontTopRightCorner   )
  .def("nodesBackBottomLeftCorner"  , &SM::FineLayer::nodesBackBottomLeftCorner  )
  .def("nodesBackBottomRightCorner" , &SM::FineLayer::nodesBackBottomRightCorner )
  .def("nodesBackTopLeftCorner"     , &SM::FineLayer::nodesBackTopLeftCorner     )
  .def("nodesBackTopRightCorner"    , &SM::FineLayer::nodesBackTopRightCorner    )
  // boundary nodes: corners (aliases)
  .def("nodesFrontLeftBottomCorner" , &SM::FineLayer::nodesFrontLeftBottomCorner )
  .def("nodesBottomFrontLeftCorner" , &SM::FineLayer::nodesBottomFrontLeftCorner )
  .def("nodesBottomLeftFrontCorner" , &SM::FineLayer::nodesBottomLeftFrontCorner )
  .def("nodesLeftFrontBottomCorner" , &SM::FineLayer::nodesLeftFrontBottomCorner )
  .def("nodesLeftBottomFrontCorner" , &SM::FineLayer::nodesLeftBottomFrontCorner )
  .def("nodesFrontRightBottomCorner", &SM::FineLayer::nodesFrontRightBottomCorner)
  .def("nodesBottomFrontRightCorner", &SM::FineLayer::nodesBottomFrontRightCorner)
  .def("nodesBottomRightFrontCorner", &SM::FineLayer::nodesBottomRightFrontCorner)
  .def("nodesRightFrontBottomCorner", &SM::FineLayer::nodesRightFrontBottomCorner)
  .def("nodesRightBottomFrontCorner", &SM::FineLayer::nodesRightBottomFrontCorner)
  .def("nodesFrontLeftTopCorner"    , &SM::FineLayer::nodesFrontLeftTopCorner    )
  .def("nodesTopFrontLeftCorner"    , &SM::FineLayer::nodesTopFrontLeftCorner    )
  .def("nodesTopLeftFrontCorner"    , &SM::FineLayer::nodesTopLeftFrontCorner    )
  .def("nodesLeftFrontTopCorner"    , &SM::FineLayer::nodesLeftFrontTopCorner    )
  .def("nodesLeftTopFrontCorner"    , &SM::FineLayer::nodesLeftTopFrontCorner    )
  .def("nodesFrontRightTopCorner"   , &SM::FineLayer::nodesFrontRightTopCorner   )
  .def("nodesTopFrontRightCorner"   , &SM::FineLayer::nodesTopFrontRightCorner   )
  .def("nodesTopRightFrontCorner"   , &SM::FineLayer::nodesTopRightFrontCorner   )
  .def("nodesRightFrontTopCorner"   , &SM::FineLayer::nodesRightFrontTopCorner   )
  .def("nodesRightTopFrontCorner"   , &SM::FineLayer::nodesRightTopFrontCorner   )
  .def("nodesBackLeftBottomCorner"  , &SM::FineLayer::nodesBackLeftBottomCorner  )
  .def("nodesBottomBackLeftCorner"  , &SM::FineLayer::nodesBottomBackLeftCorner  )
  .def("nodesBottomLeftBackCorner"  , &SM::FineLayer::nodesBottomLeftBackCorner  )
  .def("nodesLeftBackBottomCorner"  , &SM::FineLayer::nodesLeftBackBottomCorner  )
  .def("nodesLeftBottomBackCorner"  , &SM::FineLayer::nodesLeftBottomBackCorner  )
  .def("nodesBackRightBottomCorner" , &SM::FineLayer::nodesBackRightBottomCorner )
  .def("nodesBottomBackRightCorner" , &SM::FineLayer::nodesBottomBackRightCorner )
  .def("nodesBottomRightBackCorner" , &SM::FineLayer::nodesBottomRightBackCorner )
  .def("nodesRightBackBottomCorner" , &SM::FineLayer::nodesRightBackBottomCorner )
  .def("nodesRightBottomBackCorner" , &SM::FineLayer::nodesRightBottomBackCorner )
  .def("nodesBackLeftTopCorner"     , &SM::FineLayer::nodesBackLeftTopCorner     )
  .def("nodesTopBackLeftCorner"     , &SM::FineLayer::nodesTopBackLeftCorner     )
  .def("nodesTopLeftBackCorner"     , &SM::FineLayer::nodesTopLeftBackCorner     )
  .def("nodesLeftBackTopCorner"     , &SM::FineLayer::nodesLeftBackTopCorner     )
  .def("nodesLeftTopBackCorner"     , &SM::FineLayer::nodesLeftTopBackCorner     )
  .def("nodesBackRightTopCorner"    , &SM::FineLayer::nodesBackRightTopCorner    )
  .def("nodesTopBackRightCorner"    , &SM::FineLayer::nodesTopBackRightCorner    )
  .def("nodesTopRightBackCorner"    , &SM::FineLayer::nodesTopRightBackCorner    )
  .def("nodesRightBackTopCorner"    , &SM::FineLayer::nodesRightBackTopCorner    )
  .def("nodesRightTopBackCorner"    , &SM::FineLayer::nodesRightTopBackCorner    )
  // periodicity
  .def("nodesPeriodic"              , &SM::FineLayer::nodesPeriodic              )
  .def("nodesOrigin"                , &SM::FineLayer::nodesOrigin                )
  .def("dofs"                       , &SM::FineLayer::dofs                       )
  .def("dofsPeriodic"               , &SM::FineLayer::dofsPeriodic               )
  // print to screen
  .def("__repr__",
    [](const SM::FineLayer &a){ return "<GooseFEM.Mesh.Hex8.FineLayer>"; }
  );

// -------------------------------------------------------------------------------------------------

}

// ========================= GooseFEM::Mesh::Quad4 - GooseFEM/MeshQuad4.h ==========================

{

// create submodule
py::module sm = mMesh.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

// abbreviate name-space
namespace SM = GooseFEM::Mesh::Quad4;

// -------------------------------------------------------------------------------------------------

py::class_<SM::Regular>(sm, "Regular")

  .def(
    py::init<size_t,size_t,double>(),
    "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.
  )

  .def("coor"                  , &SM::Regular::coor                  )
  .def("conn"                  , &SM::Regular::conn                  )
  .def("nelem"                 , &SM::Regular::nelem                 )
  .def("nnode"                 , &SM::Regular::nnode                 )
  .def("nne"                   , &SM::Regular::nne                   )
  .def("ndim"                  , &SM::Regular::ndim                  )
  .def("nodesBottomEdge"       , &SM::Regular::nodesBottomEdge       )
  .def("nodesTopEdge"          , &SM::Regular::nodesTopEdge          )
  .def("nodesLeftEdge"         , &SM::Regular::nodesLeftEdge         )
  .def("nodesRightEdge"        , &SM::Regular::nodesRightEdge        )
  .def("nodesBottomOpenEdge"   , &SM::Regular::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &SM::Regular::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &SM::Regular::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &SM::Regular::nodesRightOpenEdge    )
  .def("nodesBottomLeftCorner" , &SM::Regular::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &SM::Regular::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &SM::Regular::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &SM::Regular::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &SM::Regular::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &SM::Regular::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &SM::Regular::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &SM::Regular::nodesRightTopCorner   )
  .def("nodesPeriodic"         , &SM::Regular::nodesPeriodic         )
  .def("nodesOrigin"           , &SM::Regular::nodesOrigin           )
  .def("dofs"                  , &SM::Regular::dofs                  )
  .def("dofsPeriodic"          , &SM::Regular::dofsPeriodic          )

  .def("__repr__",
    [](const SM::Regular &a){ return "<GooseFEM.Mesh.Quad4.Regular>"; }
  );

// -------------------------------------------------------------------------------------------------

py::class_<SM::FineLayer>(sm, "FineLayer")

  .def(
    py::init<size_t,size_t,double,size_t>(),
    "FineLayer mesh: 'nx' pixels in horizontal direction (length 'Lx'), idem in vertical direction",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.,
    py::arg("nfine")=1
  )

  .def("shape"                 , &SM::FineLayer::shape                 )
  .def("coor"                  , &SM::FineLayer::coor                  )
  .def("conn"                  , &SM::FineLayer::conn                  )
  .def("nelem"                 , &SM::FineLayer::nelem                 )
  .def("nnode"                 , &SM::FineLayer::nnode                 )
  .def("nne"                   , &SM::FineLayer::nne                   )
  .def("ndim"                  , &SM::FineLayer::ndim                  )
  .def("elementsMiddleLayer"   , &SM::FineLayer::elementsMiddleLayer   )
  .def("nodesBottomEdge"       , &SM::FineLayer::nodesBottomEdge       )
  .def("nodesTopEdge"          , &SM::FineLayer::nodesTopEdge          )
  .def("nodesLeftEdge"         , &SM::FineLayer::nodesLeftEdge         )
  .def("nodesRightEdge"        , &SM::FineLayer::nodesRightEdge        )
  .def("nodesBottomOpenEdge"   , &SM::FineLayer::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &SM::FineLayer::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &SM::FineLayer::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &SM::FineLayer::nodesRightOpenEdge    )
  .def("nodesBottomLeftCorner" , &SM::FineLayer::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &SM::FineLayer::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &SM::FineLayer::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &SM::FineLayer::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &SM::FineLayer::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &SM::FineLayer::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &SM::FineLayer::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &SM::FineLayer::nodesRightTopCorner   )
  .def("nodesPeriodic"         , &SM::FineLayer::nodesPeriodic         )
  .def("nodesOrigin"           , &SM::FineLayer::nodesOrigin           )
  .def("dofs"                  , &SM::FineLayer::dofs                  )
  .def("dofsPeriodic"          , &SM::FineLayer::dofsPeriodic          )

  .def("__repr__",
    [](const SM::FineLayer &a){ return "<GooseFEM.Mesh.Quad4.FineLayer>"; }
  );

// -------------------------------------------------------------------------------------------------

}

// ========================== GooseFEM::Mesh::Tri3 - GooseFEM/MeshTri3.h ===========================

{

// create submodule
py::module sm = mMesh.def_submodule("Tri3" , "Linear triangular elements (2D)");

// abbreviate name-space
namespace SM = GooseFEM::Mesh::Tri3;

// -------------------------------------------------------------------------------------------------

py::class_<SM::Regular>(sm, "Regular")

  .def(
    py::init<size_t,size_t,double>(),
    "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.
  )

  .def("coor"                  , &SM::Regular::coor                  )
  .def("conn"                  , &SM::Regular::conn                  )
  .def("nelem"                 , &SM::Regular::nelem                 )
  .def("nnode"                 , &SM::Regular::nnode                 )
  .def("nne"                   , &SM::Regular::nne                   )
  .def("ndim"                  , &SM::Regular::ndim                  )
  .def("nodesBottomEdge"       , &SM::Regular::nodesBottomEdge       )
  .def("nodesTopEdge"          , &SM::Regular::nodesTopEdge          )
  .def("nodesLeftEdge"         , &SM::Regular::nodesLeftEdge         )
  .def("nodesRightEdge"        , &SM::Regular::nodesRightEdge        )
  .def("nodesBottomOpenEdge"   , &SM::Regular::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &SM::Regular::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &SM::Regular::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &SM::Regular::nodesRightOpenEdge    )
  .def("nodesBottomLeftCorner" , &SM::Regular::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &SM::Regular::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &SM::Regular::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &SM::Regular::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &SM::Regular::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &SM::Regular::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &SM::Regular::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &SM::Regular::nodesRightTopCorner   )
  .def("nodesPeriodic"         , &SM::Regular::nodesPeriodic         )
  .def("nodesOrigin"           , &SM::Regular::nodesOrigin           )
  .def("dofs"                  , &SM::Regular::dofs                  )
  .def("dofsPeriodic"          , &SM::Regular::dofsPeriodic          )

  .def("__repr__",
    [](const SM::Regular &a){ return "<GooseFEM.Mesh.Tri3.Regular>"; }
  );

// -------------------------------------------------------------------------------------------------

sm.def("getOrientation",&SM::getOrientation,
  "Get the orientation of each element",
  py::arg("coor"),
  py::arg("conn")
);

// -------------------------------------------------------------------------------------------------

sm.def("retriangulate",&SM::retriangulate,
  "Re-triangulate existing mesh",
  py::arg("coor"),
  py::arg("conn"),
  py::arg("orientation")=-1
);

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}


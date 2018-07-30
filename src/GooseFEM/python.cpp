/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cppmat/pybind11.h>

#include <pyxtensor/pyxtensor.hpp>

#include "GooseFEM.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;
namespace M  = GooseFEM;

// =================================================================================================

class PyGeometry : public GooseFEM::Dynamics::Geometry
{
public:
  // inherit the constructors
  using GooseFEM::Dynamics::Geometry::Geometry;
  using Arr1 = xt::xtensor<double,1>;
  using Arr2 = xt::xtensor<double,2>;
  using Arr3 = xt::xtensor<double,3>;

  // trampoline
  xt::xtensor<double,1> solve_A()                  override { PYBIND11_OVERLOAD_PURE( Arr1, GooseFEM::Dynamics::Geometry, solve_A         ); }
  xt::xtensor<double,1> solve_V()                  override { PYBIND11_OVERLOAD_PURE( Arr1, GooseFEM::Dynamics::Geometry, solve_V         ); }
  xt::xtensor<double,2> u()      const             override { PYBIND11_OVERLOAD_PURE( Arr2, GooseFEM::Dynamics::Geometry, u               ); }
  xt::xtensor<double,2> v()      const             override { PYBIND11_OVERLOAD_PURE( Arr2, GooseFEM::Dynamics::Geometry, v               ); }
  xt::xtensor<double,2> a()      const             override { PYBIND11_OVERLOAD_PURE( Arr2, GooseFEM::Dynamics::Geometry, a               ); }
  xt::xtensor<double,1> dofs_u() const             override { PYBIND11_OVERLOAD_PURE( Arr1, GooseFEM::Dynamics::Geometry, dofs_u          ); }
  xt::xtensor<double,1> dofs_v() const             override { PYBIND11_OVERLOAD_PURE( Arr1, GooseFEM::Dynamics::Geometry, dofs_v          ); }
  xt::xtensor<double,1> dofs_a() const             override { PYBIND11_OVERLOAD_PURE( Arr1, GooseFEM::Dynamics::Geometry, dofs_a          ); }
  void set_u(const xt::xtensor<double,2> &nodevec) override { PYBIND11_OVERLOAD_PURE( void, GooseFEM::Dynamics::Geometry, set_u  , nodevec); }
  void set_u(const xt::xtensor<double,1> &dofval ) override { PYBIND11_OVERLOAD_PURE( void, GooseFEM::Dynamics::Geometry, set_u  , dofval ); }
  void set_v(const xt::xtensor<double,1> &dofval ) override { PYBIND11_OVERLOAD_PURE( void, GooseFEM::Dynamics::Geometry, set_v  , dofval ); }
  void set_a(const xt::xtensor<double,1> &dofval ) override { PYBIND11_OVERLOAD_PURE( void, GooseFEM::Dynamics::Geometry, set_a  , dofval ); }
};

// =========================================== GooseFEM ============================================

PYBIND11_MODULE(GooseFEM, m) {

m.doc() = "Some simple finite element meshes and operations";

// ================================= GooseFEM - GooseFEM/Vector.h ==================================

py::class_<GooseFEM::Vector>(m, "Vector")
  // constructor
  .def(
    py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &>(),
    "Class to switch between DOF/nodal/element views of vectors",
    py::arg("conn"),
    py::arg("dofs")
  )
  // constructor
  .def(
    py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(),
    "Class to switch between DOF/nodal/element views of vectors",
    py::arg("conn"),
    py::arg("dofs"),
    py::arg("iip")
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
  .def("asDofs"    , py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::Vector::asDofs    , py::const_))
  .def("asDofs"    , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::Vector::asDofs    , py::const_))
  .def("asDofs"    , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::Vector::asDofs    , py::const_))
  .def("asDofs_u"  , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::Vector::asDofs_u  , py::const_))
  .def("asDofs_u"  , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::Vector::asDofs_u  , py::const_))
  .def("asDofs_p"  , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::Vector::asDofs_p  , py::const_))
  .def("asDofs_p"  , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::Vector::asDofs_p  , py::const_))
  .def("asNode"    , py::overload_cast<const xt::xtensor<double,1>&                             >(&M::Vector::asNode    , py::const_))
  .def("asNode"    , py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::Vector::asNode    , py::const_))
  .def("asNode"    , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::Vector::asNode    , py::const_))
  .def("asElement" , py::overload_cast<const xt::xtensor<double,1>&                             >(&M::Vector::asElement , py::const_))
  .def("asElement" , py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::Vector::asElement , py::const_))
  .def("asElement" , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::Vector::asElement , py::const_))
  // -
  .def("assembleDofs"  , py::overload_cast<const xt::xtensor<double,2>&>(&M::Vector::assembleDofs  , py::const_))
  .def("assembleDofs"  , py::overload_cast<const xt::xtensor<double,3>&>(&M::Vector::assembleDofs  , py::const_))
  .def("assembleDofs_u", py::overload_cast<const xt::xtensor<double,2>&>(&M::Vector::assembleDofs_u, py::const_))
  .def("assembleDofs_u", py::overload_cast<const xt::xtensor<double,3>&>(&M::Vector::assembleDofs_u, py::const_))
  .def("assembleDofs_p", py::overload_cast<const xt::xtensor<double,2>&>(&M::Vector::assembleDofs_p, py::const_))
  .def("assembleDofs_p", py::overload_cast<const xt::xtensor<double,3>&>(&M::Vector::assembleDofs_p, py::const_))
  .def("assembleNode"  , py::overload_cast<const xt::xtensor<double,3>&>(&M::Vector::assembleNode  , py::const_))
  // print to screen
  .def("__repr__",
    [](const GooseFEM::Vector &){ return "<GooseFEM.Vector>"; }
  );

// ============================= GooseFEM - GooseFEM/MatrixDiagonal.h ==============================

py::class_<GooseFEM::MatrixDiagonal>(m, "MatrixDiagonal")
// constructor
  .def(
    py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &>(),
    "Class to switch between DOF/nodal/element views of vectors",
    py::arg("conn"),
    py::arg("dofs")
  )
  // constructor
  .def(
    py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(),
    "Class to switch between DOF/nodal/element views of vectors",
    py::arg("conn"),
    py::arg("dofs"),
    py::arg("iip")
  )
  // methods
  .def("nelem", &M::MatrixDiagonal::nelem)
  .def("nne"  , &M::MatrixDiagonal::nne  )
  .def("nnode", &M::MatrixDiagonal::nnode)
  .def("ndim" , &M::MatrixDiagonal::ndim )
  .def("ndof" , &M::MatrixDiagonal::ndof )
  .def("nnu"  , &M::MatrixDiagonal::nnu  )
  .def("nnp"  , &M::MatrixDiagonal::nnp  )
  // -
  .def("iiu"  , &M::MatrixDiagonal::iiu  )
  .def("iip"  , &M::MatrixDiagonal::iip  )
  // -
  .def("dot"  ,                                                                              &M::MatrixDiagonal::dot               )
  .def("dot_u", py::overload_cast<const xt::xtensor<double,1>&                             >(&M::MatrixDiagonal::dot_u, py::const_))
  .def("dot_u", py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::MatrixDiagonal::dot_u, py::const_))
  .def("dot_p", py::overload_cast<const xt::xtensor<double,1>&                             >(&M::MatrixDiagonal::dot_p, py::const_))
  .def("dot_p", py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::MatrixDiagonal::dot_p, py::const_))
  // -
  .def("check_diagonal", &M::MatrixDiagonal::check_diagonal)
  .def("assemble"      , &M::MatrixDiagonal::assemble      )
  // .def("set"           , &M::MatrixDiagonal::set           )
  // .def("set_uu"        , &M::MatrixDiagonal::set_uu        )
  // .def("set_pp"        , &M::MatrixDiagonal::set_pp        )
  .def("solve"         , py::overload_cast<const xt::xtensor<double,1>&                              >(&M::MatrixDiagonal::solve  ), "Solve", py::arg("rhs"  )                )
  .def("solve"         , py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(&M::MatrixDiagonal::solve  ), "Solve", py::arg("rhs"  ), py::arg("u_p"))
  .def("solve_u"       , py::overload_cast<const xt::xtensor<double,1>&                              >(&M::MatrixDiagonal::solve_u), "Solve", py::arg("rhs_u")                )
  .def("solve_u"       , py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(&M::MatrixDiagonal::solve_u), "Solve", py::arg("rhs_u"), py::arg("u_p"))
  // .def("rhs_p"         , &M::MatrixDiagonal::rhs_p         )
  .def("asDiagonal"    , &M::MatrixDiagonal::asDiagonal    )
  // .def("asDiagonal_uu" , &M::MatrixDiagonal::asDiagonal_uu )
  // .def("asDiagonal_pp" , &M::MatrixDiagonal::asDiagonal_pp )
  // .def("asSparse"      , &M::MatrixDiagonal::asSparse      )
  // .def("asSparse_uu"   , &M::MatrixDiagonal::asSparse_uu   )
  // .def("asSparse_up"   , &M::MatrixDiagonal::asSparse_up   )
  // .def("asSparse_pu"   , &M::MatrixDiagonal::asSparse_pu   )
  // .def("asSparse_pp"   , &M::MatrixDiagonal::asSparse_pp   )
  // .def("asDense"       , &M::MatrixDiagonal::asDense       )
  // .def("asDense_uu"    , &M::MatrixDiagonal::asDense_uu    )
  // .def("asDense_up"    , &M::MatrixDiagonal::asDense_up    )
  // .def("asDense_pu"    , &M::MatrixDiagonal::asDense_pu    )
  // .def("asDense_pp"    , &M::MatrixDiagonal::asDense_pp    )

  // print to screen
  .def("__repr__",
    [](const GooseFEM::MatrixDiagonal &){ return "<GooseFEM.MatrixDiagonal>"; }
  );

// =========================== GooseFEM::Dynamics - GooseFEM/Dynamics.h ============================

py::module mDynamics = m.def_submodule("Dynamics", "Solve routines for dynamic FEM");

// -------------------------------------------------------------------------------------------------

mDynamics.def("Verlet",
  &GooseFEM::Dynamics::Verlet,
  "Verlet time integration",
  py::arg("geometry"),
  py::arg("dt"),
  py::arg("nstep")=1
);

// -------------------------------------------------------------------------------------------------

mDynamics.def("velocityVerlet",
  &GooseFEM::Dynamics::velocityVerlet,
  "Velocity-Verlet time integration",
  py::arg("geometry"),
  py::arg("dt"),
  py::arg("nstep")=1
);

// // -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Dynamics::Geometry, PyGeometry>(mDynamics, "Geometry")
  // constructor
  .def(py::init<>())
  // methods
  .def("solve_A"   , &GooseFEM::Dynamics::Geometry::solve_A)
  .def("solve_V"   , &GooseFEM::Dynamics::Geometry::solve_V)
  .def("u"         , &GooseFEM::Dynamics::Geometry::u      )
  .def("v"         , &GooseFEM::Dynamics::Geometry::v      )
  .def("a"         , &GooseFEM::Dynamics::Geometry::a      )
  .def("dofs_u"    , &GooseFEM::Dynamics::Geometry::dofs_u )
  .def("dofs_v"    , &GooseFEM::Dynamics::Geometry::dofs_v )
  .def("dofs_a"    , &GooseFEM::Dynamics::Geometry::dofs_a )
  .def("set_v"     , &GooseFEM::Dynamics::Geometry::set_v  )
  .def("set_a"     , &GooseFEM::Dynamics::Geometry::set_a  )
  .def("set_u"     , py::overload_cast<const xt::xtensor<double,2> &>(&GooseFEM::Dynamics::Geometry::set_u))
  .def("set_u"     , py::overload_cast<const xt::xtensor<double,1> &>(&GooseFEM::Dynamics::Geometry::set_u))

  // print to screen
  .def("__repr__",
    [](const GooseFEM::Dynamics::Geometry &){ return "<GooseDEM.Dynamics.Geometry>"; }
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

// create sub-module
py::module sm = mElement.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

// abbreviate name-space
namespace SM = GooseFEM::Element::Quad4;

// -------------------------------------------------------------------------------------------------

py::class_<SM::Quadrature>(sm, "Quadrature")
  // constructor
  .def(
    py::init<const xt::xtensor<double,3> &>(),
    "Quadrature",
    py::arg("x")
  )
  // constructor
  .def(
    py::init<const xt::xtensor<double,3> &, const xt::xtensor<double,2> &, const xt::xtensor<double,1> &>(),
    "Quadrature",
    py::arg("x"),
    py::arg("xi"),
    py::arg("w")
  )
  // sizes
  .def("nelem"                    , &SM::Quadrature::nelem)
  .def("nne"                      , &SM::Quadrature::nne)
  .def("ndim"                     , &SM::Quadrature::ndim)
  .def("nip"                      , &SM::Quadrature::nip)
  .def("dV"                       , py::overload_cast<>(&SM::Quadrature::dV, py::const_))
  .def("gradN_vector"             , py::overload_cast<const xt::xtensor<double,3> &>(&SM::Quadrature::gradN_vector, py::const_))
  .def("gradN_vector_T"           , py::overload_cast<const xt::xtensor<double,3> &>(&SM::Quadrature::gradN_vector_T, py::const_))
  .def("symGradN_vector"          , py::overload_cast<const xt::xtensor<double,3> &>(&SM::Quadrature::symGradN_vector, py::const_))
  .def("int_N_scalar_NT_dV"       , py::overload_cast<const xt::xtensor<double,2> &>(&SM::Quadrature::int_N_scalar_NT_dV, py::const_))
  .def("int_gradN_dot_tensor2_dV" , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Quadrature::int_gradN_dot_tensor2_dV, py::const_))
  // print to screen
  .def("__repr__",
    [](const SM::Quadrature &){ return "<GooseFEM.Element.Quad4.Quadrature>"; }
  );

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Gauss", "Gauss quadrature");

namespace SSM = GooseFEM::Element::Quad4::Gauss;

ssm.def("nip", &SSM::nip);
ssm.def("xi" , &SSM::xi);
ssm.def("w"  , &SSM::w);

}

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Nodal", "Nodal quadrature");

namespace SSM = GooseFEM::Element::Quad4::Nodal;

ssm.def("nip", &SSM::nip);
ssm.def("xi" , &SSM::xi);
ssm.def("w"  , &SSM::w);

}

// -------------------------------------------------------------------------------------------------

}

// ====================== GooseFEM::Element::Hex8 - GooseFEM/ElementHex8.h =======================

{

// create sub-module
py::module sm = mElement.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

// abbreviate name-space
namespace SM = GooseFEM::Element::Hex8;

// -------------------------------------------------------------------------------------------------

py::class_<SM::Quadrature>(sm, "Quadrature")
  // constructor
  .def(
    py::init<const xt::xtensor<double,3> &>(),
    "Quadrature",
    py::arg("x")
  )
  // constructor
  .def(
    py::init<const xt::xtensor<double,3> &, const xt::xtensor<double,2> &, const xt::xtensor<double,1> &>(),
    "Quadrature",
    py::arg("x"),
    py::arg("xi"),
    py::arg("w")
  )
  // sizes
  .def("nelem"                    , &SM::Quadrature::nelem)
  .def("nne"                      , &SM::Quadrature::nne)
  .def("ndim"                     , &SM::Quadrature::ndim)
  .def("nip"                      , &SM::Quadrature::nip)
  .def("dV"                       , py::overload_cast<>(&SM::Quadrature::dV, py::const_))
  .def("gradN_vector"             , py::overload_cast<const xt::xtensor<double,3> &>(&SM::Quadrature::gradN_vector, py::const_))
  .def("gradN_vector_T"           , py::overload_cast<const xt::xtensor<double,3> &>(&SM::Quadrature::gradN_vector_T, py::const_))
  .def("symGradN_vector"          , py::overload_cast<const xt::xtensor<double,3> &>(&SM::Quadrature::symGradN_vector, py::const_))
  .def("int_N_scalar_NT_dV"       , py::overload_cast<const xt::xtensor<double,2> &>(&SM::Quadrature::int_N_scalar_NT_dV, py::const_))
  .def("int_gradN_dot_tensor2_dV" , py::overload_cast<const xt::xtensor<double,4> &>(&SM::Quadrature::int_gradN_dot_tensor2_dV, py::const_))
  // print to screen
  .def("__repr__",
    [](const SM::Quadrature &){ return "<GooseFEM.Element.Hex8.Quadrature>"; }
  );

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Gauss", "Gauss quadrature");

namespace SSM = GooseFEM::Element::Hex8::Gauss;

ssm.def("nip", &SSM::nip);
ssm.def("xi" , &SSM::xi);
ssm.def("w"  , &SSM::w);

}

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Nodal", "Nodal quadrature");

namespace SSM = GooseFEM::Element::Hex8::Nodal;

ssm.def("nip", &SSM::nip);
ssm.def("xi" , &SSM::xi);
ssm.def("w"  , &SSM::w);

}

// -------------------------------------------------------------------------------------------------

}

// =============================== GooseFEM::Mesh - GooseFEM/Mesh.h ================================

py::module mMesh = m.def_submodule("Mesh", "Generic mesh routines");

// -------------------------------------------------------------------------------------------------

mMesh.def("elem2node", &GooseFEM::Mesh::elem2node,
  "Elements connect to each node: [ number of elements , element numbers ]",
  py::arg("conn")
);

// -------------------------------------------------------------------------------------------------

mMesh.def("coordination", &GooseFEM::Mesh::coordination,
  "Get the coordination number of each node: the number of elements connected to it",
  py::arg("conn")
);

// -------------------------------------------------------------------------------------------------

mMesh.def("dofs", &GooseFEM::Mesh::dofs,
  "List with DOF-numbers (in sequential order)",
  py::arg("nnode"),
  py::arg("ndim")
);

// -------------------------------------------------------------------------------------------------

using renumber = xt::xtensor<size_t,2>(const xt::xtensor<size_t,2> &);

mMesh.def("renumber",
  py::overload_cast<const xt::xtensor<size_t,2> &>((renumber*)&GooseFEM::Mesh::renumber),
  "Renumber DOF-list to use the lowest possible index",
  py::arg("dofs")
);

// -------------------------------------------------------------------------------------------------

using reorder = xt::xtensor<size_t,2>(const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1>&, std::string);

mMesh.def("reorder",
  py::overload_cast<const xt::xtensor<size_t,2>&,const xt::xtensor<size_t,1>&,std::string>(
    (reorder*)&GooseFEM::Mesh::reorder),
  "Renumber DOF-list to begin or end with 'idx'",
  py::arg("dofs"),
  py::arg("idx"),
  py::arg("location")="end"
);

// ========================== GooseFEM::Mesh::Hex8 - GooseFEM/MeshHex8.h ===========================

{

// create sub-module
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
    [](const SM::Regular &){ return "<GooseFEM.Mesh.Hex8.Regular>"; }
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
    [](const SM::FineLayer &){ return "<GooseFEM.Mesh.Hex8.FineLayer>"; }
  );

// -------------------------------------------------------------------------------------------------

}

// ========================= GooseFEM::Mesh::Quad4 - GooseFEM/MeshQuad4.h ==========================

{

// create sub-module
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
    [](const SM::Regular &){ return "<GooseFEM.Mesh.Quad4.Regular>"; }
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
    [](const SM::FineLayer &){ return "<GooseFEM.Mesh.Quad4.FineLayer>"; }
  );

// -------------------------------------------------------------------------------------------------

}

// ========================== GooseFEM::Mesh::Tri3 - GooseFEM/MeshTri3.h ===========================

{

// create sub-module
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
    [](const SM::Regular &){ return "<GooseFEM.Mesh.Tri3.Regular>"; }
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


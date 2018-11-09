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

namespace py = pybind11;
namespace M  = GooseFEM;

// ======================================= trampoline class ========================================

class PyGeometry : public GooseFEM::Dynamics::Geometry
{
public:
  using GooseFEM::Dynamics::Geometry::Geometry;
  using Arr1 = xt::xtensor<double,1>;
  using Arr2 = xt::xtensor<double,2>;
  using Arr3 = xt::xtensor<double,3>;

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

// ======================================== GooseFEM.Vector ========================================

py::class_<GooseFEM::VectorPartitioned>(m, "VectorPartitioned")

  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(), "Switch between dofval/nodevec/elemvec", py::arg("conn"), py::arg("dofs"), py::arg("iip"))

  .def("nelem", &M::VectorPartitioned::nelem, "Return number of element")
  .def("nne"  , &M::VectorPartitioned::nne  , "Return number of nodes per element")
  .def("nnode", &M::VectorPartitioned::nnode, "Return number of nodes")
  .def("ndim" , &M::VectorPartitioned::ndim , "Return number of dimensions")
  .def("ndof" , &M::VectorPartitioned::ndof , "Return number of degrees-of-freedom")
  .def("nnu"  , &M::VectorPartitioned::nnu  , "Return number of unknown degrees-of-freedom")
  .def("nnp"  , &M::VectorPartitioned::nnp  , "Return number of prescribed degrees-of-freedom")

  .def("dofs" , &M::VectorPartitioned::dofs , "Return degrees-of-freedom")
  .def("iiu"  , &M::VectorPartitioned::iiu  , "Return unknown degrees-of-freedom")
  .def("iip"  , &M::VectorPartitioned::iip  , "Return prescribed degrees-of-freedom")

  .def("asDofs"   , py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::VectorPartitioned::asDofs   , py::const_), "Set 'dofval" , py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asDofs"   , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::VectorPartitioned::asDofs   , py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs"   , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::VectorPartitioned::asDofs   , py::const_), "Set 'dofval" , py::arg("elemvec"))
  .def("asDofs_u" , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::VectorPartitioned::asDofs_u , py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs_u" , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::VectorPartitioned::asDofs_u , py::const_), "Set 'dofval" , py::arg("elemvec"))
  .def("asDofs_p" , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::VectorPartitioned::asDofs_p , py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs_p" , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::VectorPartitioned::asDofs_p , py::const_), "Set 'dofval" , py::arg("elemvec"))

  .def("asNode"   , py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::VectorPartitioned::asNode   , py::const_), "Set 'nodevec", py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asNode"   , py::overload_cast<const xt::xtensor<double,1>&                             >(&M::VectorPartitioned::asNode   , py::const_), "Set 'nodevec", py::arg("dofval"))
  .def("asNode"   , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::VectorPartitioned::asNode   , py::const_), "Set 'nodevec", py::arg("elemvec"))

  .def("asElement", py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::VectorPartitioned::asElement, py::const_), "Set 'elemvec", py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asElement", py::overload_cast<const xt::xtensor<double,1>&                             >(&M::VectorPartitioned::asElement, py::const_), "Set 'elemvec", py::arg("dofval"))
  .def("asElement", py::overload_cast<const xt::xtensor<double,2>&                             >(&M::VectorPartitioned::asElement, py::const_), "Set 'elemvec", py::arg("nodevec"))

  .def("assembleDofs"  , py::overload_cast<const xt::xtensor<double,2>&>(&M::VectorPartitioned::assembleDofs  , py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs"  , py::overload_cast<const xt::xtensor<double,3>&>(&M::VectorPartitioned::assembleDofs  , py::const_), "Assemble 'dofval'" , py::arg("elemvec"))
  .def("assembleDofs_u", py::overload_cast<const xt::xtensor<double,2>&>(&M::VectorPartitioned::assembleDofs_u, py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs_u", py::overload_cast<const xt::xtensor<double,3>&>(&M::VectorPartitioned::assembleDofs_u, py::const_), "Assemble 'dofval'" , py::arg("elemvec"))
  .def("assembleDofs_p", py::overload_cast<const xt::xtensor<double,2>&>(&M::VectorPartitioned::assembleDofs_p, py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs_p", py::overload_cast<const xt::xtensor<double,3>&>(&M::VectorPartitioned::assembleDofs_p, py::const_), "Assemble 'dofval'" , py::arg("elemvec"))

  .def("assembleNode"  , py::overload_cast<const xt::xtensor<double,3>&>(&M::VectorPartitioned::assembleNode  , py::const_), "Assemble 'nodevec'", py::arg("elemvec"))

  .def("__repr__", [](const GooseFEM::VectorPartitioned &){ return "<GooseFEM.Vector>"; });

// ==================================== GooseFEM.MatrixDiagonalPartitioned ====================================

py::class_<GooseFEM::MatrixDiagonalPartitioned>(m, "MatrixDiagonalPartitioned")

  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(), "Diagonal matrix", py::arg("conn"), py::arg("dofs"), py::arg("iip"))

  .def("nelem", &M::MatrixDiagonalPartitioned::nelem, "Return number of element")
  .def("nne"  , &M::MatrixDiagonalPartitioned::nne  , "Return number of nodes per element")
  .def("nnode", &M::MatrixDiagonalPartitioned::nnode, "Return number of nodes")
  .def("ndim" , &M::MatrixDiagonalPartitioned::ndim , "Return number of dimensions")
  .def("ndof" , &M::MatrixDiagonalPartitioned::ndof , "Return number of degrees-of-freedom")
  .def("nnu"  , &M::MatrixDiagonalPartitioned::nnu  , "Return number of unknown degrees-of-freedom")
  .def("nnp"  , &M::MatrixDiagonalPartitioned::nnp  , "Return number of prescribed degrees-of-freedom")

  .def("dofs" , &M::MatrixDiagonalPartitioned::dofs , "Return degrees-of-freedom")
  .def("iiu"  , &M::MatrixDiagonalPartitioned::iiu  , "Return unknown degrees-of-freedom")
  .def("iip"  , &M::MatrixDiagonalPartitioned::iip  , "Return prescribed degrees-of-freedom")

  .def("dot"  , &M::MatrixDiagonalPartitioned::dot  , "Dot product 'b_i = A_ij * x_j", py::arg("x"))
  .def("dot_u", &M::MatrixDiagonalPartitioned::dot_u, "Dot product 'b_i = A_ij * x_j", py::arg("x_u"), py::arg("x_p"))
  .def("dot_p", &M::MatrixDiagonalPartitioned::dot_p, "Dot product 'b_i = A_ij * x_j", py::arg("x_u"), py::arg("x_p"))

  .def("assemble", &M::MatrixDiagonalPartitioned::assemble, "Assemble from 'elemmat", py::arg("elemmat"))

  .def("solve", &M::MatrixDiagonalPartitioned::solve, "Solve", py::arg("b_u"), py::arg("x_p"))

  .def("asDiagonal", &M::MatrixDiagonalPartitioned::asDiagonal, "Return as diagonal matrix (column)")

  .def("__repr__", [](const GooseFEM::MatrixDiagonalPartitioned &){ return "<GooseFEM.MatrixDiagonalPartitioned>"; });

// ======================================= GooseFEM.Dynamics =======================================

py::module mDynamics = m.def_submodule("Dynamics", "Solve routines for dynamic FEM");

// -------------------------------------------------------------------------------------------------

mDynamics.def("Verlet"        , &GooseFEM::Dynamics::Verlet        , "Verlet time integration"         , py::arg("geometry"), py::arg("dt"), py::arg("nstep")=1);
mDynamics.def("velocityVerlet", &GooseFEM::Dynamics::velocityVerlet, "Velocity-Verlet time integration", py::arg("geometry"), py::arg("dt"), py::arg("nstep")=1);

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Dynamics::Geometry, PyGeometry>(mDynamics, "Geometry")

  .def(py::init<>())

  .def("solve_A"   , &GooseFEM::Dynamics::Geometry::solve_A, "Solve for accelerations (dofval)" )
  .def("solve_V"   , &GooseFEM::Dynamics::Geometry::solve_V, "Solve for velocities    (dofval)" )

  .def("u"         , &GooseFEM::Dynamics::Geometry::u      , "Return displacements (nodevec)")
  .def("v"         , &GooseFEM::Dynamics::Geometry::v      , "Return velocities    (nodevec)")
  .def("a"         , &GooseFEM::Dynamics::Geometry::a      , "Return accelerations (nodevec)")

  .def("dofs_u"    , &GooseFEM::Dynamics::Geometry::dofs_u , "Return displacements (dofval)" )
  .def("dofs_v"    , &GooseFEM::Dynamics::Geometry::dofs_v , "Return velocities    (dofval)" )
  .def("dofs_a"    , &GooseFEM::Dynamics::Geometry::dofs_a , "Return accelerations (dofval)" )

  .def("set_u"     , py::overload_cast<const xt::xtensor<double,1> &>(&GooseFEM::Dynamics::Geometry::set_u), "Overwrite displacements", py::arg("dofval" ))
  .def("set_u"     , py::overload_cast<const xt::xtensor<double,2> &>(&GooseFEM::Dynamics::Geometry::set_u), "Overwrite displacements", py::arg("nodevec"))
  .def("set_v"     ,                                                  &GooseFEM::Dynamics::Geometry::set_v , "Overwrite velocities"   , py::arg("nodevec"))
  .def("set_a"     ,                                                  &GooseFEM::Dynamics::Geometry::set_a , "Overwrite accelerations", py::arg("nodevec"))

  .def("__repr__", [](const GooseFEM::Dynamics::Geometry &){ return "<GooseDEM.Dynamics.Geometry>"; });

// ======================================= GooseFEM.Element ========================================

py::module mElement = m.def_submodule("Element", "Generic element routines");

// -------------------------------------------------------------------------------------------------

mElement.def("asElementVector"      , &GooseFEM::Element::asElementVector   , "Covert to 'elemvec'", py::arg("conn"), py::arg("nodevec"));
mElement.def("assembleElementVector", &GooseFEM::Element::assembleNodeVector, "Assemble 'nodevec'" , py::arg("conn"), py::arg("elemvec"));

// ==================================== GooseFEM.Element.Quad4 =====================================

{

py::module sm = mElement.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Element::Quad4::Quadrature>(sm, "Quadrature")

  .def(py::init<const xt::xtensor<double,3> &>(), "Quadrature", py::arg("x"))

  .def(py::init<const xt::xtensor<double,3> &, const xt::xtensor<double,2> &, const xt::xtensor<double,1> &>(), "Quadrature", py::arg("x"), py::arg("xi"), py::arg("w"))

  .def("nelem", &GooseFEM::Element::Quad4::Quadrature::nelem, "Return number of elements"          )
  .def("nne"  , &GooseFEM::Element::Quad4::Quadrature::nne  , "Return number of nodes per element" )
  .def("ndim" , &GooseFEM::Element::Quad4::Quadrature::ndim , "Return number of dimensions"        )
  .def("nip"  , &GooseFEM::Element::Quad4::Quadrature::nip  , "Return number of integration points")

  .def("dV"      , py::overload_cast<>(&GooseFEM::Element::Quad4::Quadrature::dV      , py::const_), "Integration point volume (qscalar)")
  .def("dVtensor", py::overload_cast<>(&GooseFEM::Element::Quad4::Quadrature::dVtensor, py::const_), "Integration point volume (qtensor)")

  .def("gradN_vector"   , py::overload_cast<const xt::xtensor<double,3> &>(&GooseFEM::Element::Quad4::Quadrature::gradN_vector   , py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))
  .def("gradN_vector_T" , py::overload_cast<const xt::xtensor<double,3> &>(&GooseFEM::Element::Quad4::Quadrature::gradN_vector_T , py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))
  .def("symGradN_vector", py::overload_cast<const xt::xtensor<double,3> &>(&GooseFEM::Element::Quad4::Quadrature::symGradN_vector, py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))

  .def("int_N_scalar_NT_dV"      , py::overload_cast<const xt::xtensor<double,2> &>(&GooseFEM::Element::Quad4::Quadrature::int_N_scalar_NT_dV      , py::const_), "Integration, returns 'elemmat'", py::arg("qscalar"))
  .def("int_gradN_dot_tensor2_dV", py::overload_cast<const xt::xtensor<double,4> &>(&GooseFEM::Element::Quad4::Quadrature::int_gradN_dot_tensor2_dV, py::const_), "Integration, returns 'elemvec'", py::arg("qtensor"))

  .def("__repr__", [](const GooseFEM::Element::Quad4::Quadrature &){ return "<GooseFEM.Element.Quad4.Quadrature>"; });

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Gauss", "Gauss quadrature");

ssm.def("nip", &GooseFEM::Element::Quad4::Gauss::nip, "Return number of integration point"  );
ssm.def("xi" , &GooseFEM::Element::Quad4::Gauss::xi , "Return integration point coordinates");
ssm.def("w"  , &GooseFEM::Element::Quad4::Gauss::w  , "Return integration point weights"    );

}

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Nodal", "Nodal quadrature");

ssm.def("nip", &GooseFEM::Element::Quad4::Nodal::nip, "Return number of integration point"  );
ssm.def("xi" , &GooseFEM::Element::Quad4::Nodal::xi , "Return integration point coordinates");
ssm.def("w"  , &GooseFEM::Element::Quad4::Nodal::w  , "Return integration point weights"    );

}

// -------------------------------------------------------------------------------------------------

}

// ===================================== GooseFEM.Element.Hex8 =====================================

{

py::module sm = mElement.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Element::Hex8::Quadrature>(sm, "Quadrature")

  .def(py::init<const xt::xtensor<double,3> &>(), "Quadrature", py::arg("x"))

  .def(py::init<const xt::xtensor<double,3> &, const xt::xtensor<double,2> &, const xt::xtensor<double,1> &>(), "Quadrature", py::arg("x"), py::arg("xi"), py::arg("w"))

  .def("nelem", &GooseFEM::Element::Hex8::Quadrature::nelem, "Return number of elements"          )
  .def("nne"  , &GooseFEM::Element::Hex8::Quadrature::nne  , "Return number of nodes per element" )
  .def("ndim" , &GooseFEM::Element::Hex8::Quadrature::ndim , "Return number of dimensions"        )
  .def("nip"  , &GooseFEM::Element::Hex8::Quadrature::nip  , "Return number of integration points")

  .def("dV"      , py::overload_cast<>(&GooseFEM::Element::Hex8::Quadrature::dV      , py::const_), "Integration point volume (qscalar)")
  .def("dVtensor", py::overload_cast<>(&GooseFEM::Element::Hex8::Quadrature::dVtensor, py::const_), "Integration point volume (qtensor)")

  .def("gradN_vector"   , py::overload_cast<const xt::xtensor<double,3> &>(&GooseFEM::Element::Hex8::Quadrature::gradN_vector   , py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))
  .def("gradN_vector_T" , py::overload_cast<const xt::xtensor<double,3> &>(&GooseFEM::Element::Hex8::Quadrature::gradN_vector_T , py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))
  .def("symGradN_vector", py::overload_cast<const xt::xtensor<double,3> &>(&GooseFEM::Element::Hex8::Quadrature::symGradN_vector, py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))

  .def("int_N_scalar_NT_dV"      , py::overload_cast<const xt::xtensor<double,2> &>(&GooseFEM::Element::Hex8::Quadrature::int_N_scalar_NT_dV      , py::const_), "Integration, returns 'elemmat'", py::arg("qscalar"))
  .def("int_gradN_dot_tensor2_dV", py::overload_cast<const xt::xtensor<double,4> &>(&GooseFEM::Element::Hex8::Quadrature::int_gradN_dot_tensor2_dV, py::const_), "Integration, returns 'elemvec'", py::arg("qtensor"))

  .def("__repr__", [](const GooseFEM::Element::Hex8::Quadrature &){ return "<GooseFEM.Element.Hex8.Quadrature>"; });

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Gauss", "Gauss quadrature");

ssm.def("nip", &GooseFEM::Element::Hex8::Gauss::nip, "Return number of integration point"  );
ssm.def("xi" , &GooseFEM::Element::Hex8::Gauss::xi , "Return integration point coordinates");
ssm.def("w"  , &GooseFEM::Element::Hex8::Gauss::w  , "Return integration point weights"    );

}

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Nodal", "Nodal quadrature");

ssm.def("nip", &GooseFEM::Element::Hex8::Nodal::nip, "Return number of integration point"  );
ssm.def("xi" , &GooseFEM::Element::Hex8::Nodal::xi , "Return integration point coordinates");
ssm.def("w"  , &GooseFEM::Element::Hex8::Nodal::w  , "Return integration point weights"    );

}

// -------------------------------------------------------------------------------------------------

}

// ========================================= GooseFEM.Mesh =========================================

py::module mMesh = m.def_submodule("Mesh", "Generic mesh routines");

// -------------------------------------------------------------------------------------------------

mMesh.def("dofs", &GooseFEM::Mesh::dofs, "List with DOF-numbers (in sequential order)", py::arg("nnode"), py::arg("ndim"));

mMesh.def("renumber", &GooseFEM::Mesh::renumber, "Renumber DOF-list to use the lowest possible index", py::arg("dofs"));

mMesh.def("renumber_index", &GooseFEM::Mesh::renumber_index, "Index-list to renumber", py::arg("dofs"));

mMesh.def("reorder", &GooseFEM::Mesh::reorder, "Renumber DOF-list to begin or end with 'idx'", py::arg("dofs"), py::arg("idx"), py::arg("location")="end");

mMesh.def("reorder_index", &GooseFEM::Mesh::reorder_index, "Index-list to reorder", py::arg("dofs"), py::arg("idx"), py::arg("location")="end");

mMesh.def("coordination", &GooseFEM::Mesh::coordination, "Coordination number of each node", py::arg("conn"));

mMesh.def("elem2node", &GooseFEM::Mesh::elem2node, "Elements connect to each node", py::arg("conn"));

// ====================================== GooseFEM.Mesh.Hex8 =======================================

{

py::module sm = mMesh.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Hex8::Regular>(sm, "Regular")

  .def(py::init<size_t,size_t,size_t,double>(), "mesh with nx*ny*nz 'pixels' and edge size h", py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("h")=1.)

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

  .def("nodesFrontBottomOpenEdge"   , &GooseFEM::Mesh::Hex8::Regular::nodesFrontBottomOpenEdge   )
  .def("nodesFrontTopOpenEdge"      , &GooseFEM::Mesh::Hex8::Regular::nodesFrontTopOpenEdge      )
  .def("nodesFrontLeftOpenEdge"     , &GooseFEM::Mesh::Hex8::Regular::nodesFrontLeftOpenEdge     )
  .def("nodesFrontRightOpenEdge"    , &GooseFEM::Mesh::Hex8::Regular::nodesFrontRightOpenEdge    )
  .def("nodesBackBottomOpenEdge"    , &GooseFEM::Mesh::Hex8::Regular::nodesBackBottomOpenEdge    )
  .def("nodesBackTopOpenEdge"       , &GooseFEM::Mesh::Hex8::Regular::nodesBackTopOpenEdge       )
  .def("nodesBackLeftOpenEdge"      , &GooseFEM::Mesh::Hex8::Regular::nodesBackLeftOpenEdge      )
  .def("nodesBackRightOpenEdge"     , &GooseFEM::Mesh::Hex8::Regular::nodesBackRightOpenEdge     )
  .def("nodesBottomLeftOpenEdge"    , &GooseFEM::Mesh::Hex8::Regular::nodesBottomLeftOpenEdge    )
  .def("nodesBottomRightOpenEdge"   , &GooseFEM::Mesh::Hex8::Regular::nodesBottomRightOpenEdge   )
  .def("nodesTopLeftOpenEdge"       , &GooseFEM::Mesh::Hex8::Regular::nodesTopLeftOpenEdge       )
  .def("nodesTopRightOpenEdge"      , &GooseFEM::Mesh::Hex8::Regular::nodesTopRightOpenEdge      )

  .def("nodesBottomFrontOpenEdge"   , &GooseFEM::Mesh::Hex8::Regular::nodesBottomFrontOpenEdge   )
  .def("nodesBottomBackOpenEdge"    , &GooseFEM::Mesh::Hex8::Regular::nodesBottomBackOpenEdge    )
  .def("nodesTopFrontOpenEdge"      , &GooseFEM::Mesh::Hex8::Regular::nodesTopFrontOpenEdge      )
  .def("nodesTopBackOpenEdge"       , &GooseFEM::Mesh::Hex8::Regular::nodesTopBackOpenEdge       )
  .def("nodesLeftBottomOpenEdge"    , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBottomOpenEdge    )
  .def("nodesLeftFrontOpenEdge"     , &GooseFEM::Mesh::Hex8::Regular::nodesLeftFrontOpenEdge     )
  .def("nodesLeftBackOpenEdge"      , &GooseFEM::Mesh::Hex8::Regular::nodesLeftBackOpenEdge      )
  .def("nodesLeftTopOpenEdge"       , &GooseFEM::Mesh::Hex8::Regular::nodesLeftTopOpenEdge       )
  .def("nodesRightBottomOpenEdge"   , &GooseFEM::Mesh::Hex8::Regular::nodesRightBottomOpenEdge   )
  .def("nodesRightTopOpenEdge"      , &GooseFEM::Mesh::Hex8::Regular::nodesRightTopOpenEdge      )
  .def("nodesRightFrontOpenEdge"    , &GooseFEM::Mesh::Hex8::Regular::nodesRightFrontOpenEdge    )
  .def("nodesRightBackOpenEdge"     , &GooseFEM::Mesh::Hex8::Regular::nodesRightBackOpenEdge     )

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

  .def("__repr__", [](const GooseFEM::Mesh::Hex8::Regular &){ return "<GooseFEM.Mesh.Hex8.Regular>"; });

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Hex8::FineLayer>(sm, "FineLayer")

  .def(py::init<size_t,size_t,size_t,double,size_t>(), "mesh with nx*ny*nz 'pixels' and edge size h", py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("h")=1., py::arg("nfine")=1)

  .def("nelem"                      , &GooseFEM::Mesh::Hex8::FineLayer::nelem                      )
  .def("nnode"                      , &GooseFEM::Mesh::Hex8::FineLayer::nnode                      )
  .def("nne"                        , &GooseFEM::Mesh::Hex8::FineLayer::nne                        )
  .def("ndim"                       , &GooseFEM::Mesh::Hex8::FineLayer::ndim                       )
  .def("shape"                      , &GooseFEM::Mesh::Hex8::FineLayer::shape                      )

  .def("coor"                       , &GooseFEM::Mesh::Hex8::FineLayer::coor                       )
  .def("conn"                       , &GooseFEM::Mesh::Hex8::FineLayer::conn                       )

  .def("elementsMiddleLayer"        , &GooseFEM::Mesh::Hex8::FineLayer::elementsMiddleLayer        )

  .def("nodesFront"                 , &GooseFEM::Mesh::Hex8::FineLayer::nodesFront                 )
  .def("nodesBack"                  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBack                  )
  .def("nodesLeft"                  , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeft                  )
  .def("nodesRight"                 , &GooseFEM::Mesh::Hex8::FineLayer::nodesRight                 )
  .def("nodesBottom"                , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottom                )
  .def("nodesTop"                   , &GooseFEM::Mesh::Hex8::FineLayer::nodesTop                   )

  .def("nodesFrontFace"             , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontFace             )
  .def("nodesBackFace"              , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackFace              )
  .def("nodesLeftFace"              , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFace              )
  .def("nodesRightFace"             , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFace             )
  .def("nodesBottomFace"            , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFace            )
  .def("nodesTopFace"               , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFace               )

  .def("nodesFrontBottomEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomEdge       )
  .def("nodesFrontTopEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopEdge          )
  .def("nodesFrontLeftEdge"         , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftEdge         )
  .def("nodesFrontRightEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightEdge        )
  .def("nodesBackBottomEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomEdge        )
  .def("nodesBackTopEdge"           , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopEdge           )
  .def("nodesBackLeftEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftEdge          )
  .def("nodesBackRightEdge"         , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightEdge         )
  .def("nodesBottomLeftEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftEdge        )
  .def("nodesBottomRightEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightEdge       )
  .def("nodesTopLeftEdge"           , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftEdge           )
  .def("nodesTopRightEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightEdge          )

  .def("nodesBottomFrontEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontEdge       )
  .def("nodesBottomBackEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackEdge        )
  .def("nodesTopFrontEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontEdge          )
  .def("nodesTopBackEdge"           , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackEdge           )
  .def("nodesLeftBottomEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomEdge        )
  .def("nodesLeftFrontEdge"         , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontEdge         )
  .def("nodesLeftBackEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackEdge          )
  .def("nodesLeftTopEdge"           , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopEdge           )
  .def("nodesRightBottomEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomEdge       )
  .def("nodesRightTopEdge"          , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopEdge          )
  .def("nodesRightFrontEdge"        , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontEdge        )
  .def("nodesRightBackEdge"         , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackEdge         )

  .def("nodesFrontBottomOpenEdge"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomOpenEdge   )
  .def("nodesFrontTopOpenEdge"      , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopOpenEdge      )
  .def("nodesFrontLeftOpenEdge"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftOpenEdge     )
  .def("nodesFrontRightOpenEdge"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightOpenEdge    )
  .def("nodesBackBottomOpenEdge"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomOpenEdge    )
  .def("nodesBackTopOpenEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopOpenEdge       )
  .def("nodesBackLeftOpenEdge"      , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftOpenEdge      )
  .def("nodesBackRightOpenEdge"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightOpenEdge     )
  .def("nodesBottomLeftOpenEdge"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftOpenEdge    )
  .def("nodesBottomRightOpenEdge"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightOpenEdge   )
  .def("nodesTopLeftOpenEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftOpenEdge       )
  .def("nodesTopRightOpenEdge"      , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightOpenEdge      )

  .def("nodesBottomFrontOpenEdge"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontOpenEdge   )
  .def("nodesBottomBackOpenEdge"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackOpenEdge    )
  .def("nodesTopFrontOpenEdge"      , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontOpenEdge      )
  .def("nodesTopBackOpenEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackOpenEdge       )
  .def("nodesLeftBottomOpenEdge"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomOpenEdge    )
  .def("nodesLeftFrontOpenEdge"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontOpenEdge     )
  .def("nodesLeftBackOpenEdge"      , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackOpenEdge      )
  .def("nodesLeftTopOpenEdge"       , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopOpenEdge       )
  .def("nodesRightBottomOpenEdge"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomOpenEdge   )
  .def("nodesRightTopOpenEdge"      , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopOpenEdge      )
  .def("nodesRightFrontOpenEdge"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontOpenEdge    )
  .def("nodesRightBackOpenEdge"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackOpenEdge     )

  .def("nodesFrontBottomLeftCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomLeftCorner )
  .def("nodesFrontBottomRightCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomRightCorner)
  .def("nodesFrontTopLeftCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopLeftCorner    )
  .def("nodesFrontTopRightCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopRightCorner   )
  .def("nodesBackBottomLeftCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomLeftCorner  )
  .def("nodesBackBottomRightCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomRightCorner )
  .def("nodesBackTopLeftCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopLeftCorner     )
  .def("nodesBackTopRightCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackTopRightCorner    )

  .def("nodesFrontLeftBottomCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftBottomCorner )
  .def("nodesBottomFrontLeftCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontLeftCorner )
  .def("nodesBottomLeftFrontCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftFrontCorner )
  .def("nodesLeftFrontBottomCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontBottomCorner )
  .def("nodesLeftBottomFrontCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomFrontCorner )
  .def("nodesFrontRightBottomCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightBottomCorner)
  .def("nodesBottomFrontRightCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontRightCorner)
  .def("nodesBottomRightFrontCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightFrontCorner)
  .def("nodesRightFrontBottomCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontBottomCorner)
  .def("nodesRightBottomFrontCorner", &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomFrontCorner)
  .def("nodesFrontLeftTopCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftTopCorner    )
  .def("nodesTopFrontLeftCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontLeftCorner    )
  .def("nodesTopLeftFrontCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftFrontCorner    )
  .def("nodesLeftFrontTopCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontTopCorner    )
  .def("nodesLeftTopFrontCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopFrontCorner    )
  .def("nodesFrontRightTopCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightTopCorner   )
  .def("nodesTopFrontRightCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontRightCorner   )
  .def("nodesTopRightFrontCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightFrontCorner   )
  .def("nodesRightFrontTopCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontTopCorner   )
  .def("nodesRightTopFrontCorner"   , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopFrontCorner   )
  .def("nodesBackLeftBottomCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftBottomCorner  )
  .def("nodesBottomBackLeftCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackLeftCorner  )
  .def("nodesBottomLeftBackCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftBackCorner  )
  .def("nodesLeftBackBottomCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackBottomCorner  )
  .def("nodesLeftBottomBackCorner"  , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomBackCorner  )
  .def("nodesBackRightBottomCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightBottomCorner )
  .def("nodesBottomBackRightCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackRightCorner )
  .def("nodesBottomRightBackCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightBackCorner )
  .def("nodesRightBackBottomCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackBottomCorner )
  .def("nodesRightBottomBackCorner" , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomBackCorner )
  .def("nodesBackLeftTopCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftTopCorner     )
  .def("nodesTopBackLeftCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackLeftCorner     )
  .def("nodesTopLeftBackCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftBackCorner     )
  .def("nodesLeftBackTopCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackTopCorner     )
  .def("nodesLeftTopBackCorner"     , &GooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopBackCorner     )
  .def("nodesBackRightTopCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesBackRightTopCorner    )
  .def("nodesTopBackRightCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopBackRightCorner    )
  .def("nodesTopRightBackCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesTopRightBackCorner    )
  .def("nodesRightBackTopCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightBackTopCorner    )
  .def("nodesRightTopBackCorner"    , &GooseFEM::Mesh::Hex8::FineLayer::nodesRightTopBackCorner    )

  .def("nodesPeriodic"              , &GooseFEM::Mesh::Hex8::FineLayer::nodesPeriodic              )
  .def("nodesOrigin"                , &GooseFEM::Mesh::Hex8::FineLayer::nodesOrigin                )

  .def("dofs"                       , &GooseFEM::Mesh::Hex8::FineLayer::dofs                       )
  .def("dofsPeriodic"               , &GooseFEM::Mesh::Hex8::FineLayer::dofsPeriodic               )

  .def("__repr__", [](const GooseFEM::Mesh::Hex8::FineLayer &){ return "<GooseFEM.Mesh.Hex8.FineLayer>"; });

// -------------------------------------------------------------------------------------------------

}

// ====================================== GooseFEM.Mesh.Quad4 ======================================

{

py::module sm = mMesh.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Quad4::Regular>(sm, "Regular")

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

py::class_<GooseFEM::Mesh::Quad4::FineLayer>(sm, "FineLayer")

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

// -------------------------------------------------------------------------------------------------

}

// ====================================== GooseFEM.Mesh.Tri3 =======================================

{

py::module sm = mMesh.def_submodule("Tri3" , "Linear triangular elements (2D)");

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Tri3::Regular>(sm, "Regular")

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
  .def("nodesBottomOpenEdge"   , &GooseFEM::Mesh::Tri3::Regular::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &GooseFEM::Mesh::Tri3::Regular::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &GooseFEM::Mesh::Tri3::Regular::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &GooseFEM::Mesh::Tri3::Regular::nodesRightOpenEdge    )

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

  .def("__repr__", [](const GooseFEM::Mesh::Tri3::Regular &){ return "<GooseFEM.Mesh.Tri3.Regular>"; });

// -------------------------------------------------------------------------------------------------

sm.def("getOrientation", &GooseFEM::Mesh::Tri3::getOrientation, "Get the orientation of each element", py::arg("coor"), py::arg("conn"));

sm.def("retriangulate", &GooseFEM::Mesh::Tri3::retriangulate, "Re-triangulate existing mesh", py::arg("coor"), py::arg("conn"), py::arg("orientation")=-1);

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}

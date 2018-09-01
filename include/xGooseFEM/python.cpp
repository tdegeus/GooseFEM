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
namespace M  = xGooseFEM;

// ======================================= trampoline class ========================================

class PyGeometry : public xGooseFEM::Dynamics::Geometry
{
public:
  using xGooseFEM::Dynamics::Geometry::Geometry;
  using Arr1 = xt::xtensor<double,1>;
  using Arr2 = xt::xtensor<double,2>;
  using Arr3 = xt::xtensor<double,3>;

  xt::xtensor<double,1> solve_A()                  override { PYBIND11_OVERLOAD_PURE( Arr1, xGooseFEM::Dynamics::Geometry, solve_A         ); }
  xt::xtensor<double,1> solve_V()                  override { PYBIND11_OVERLOAD_PURE( Arr1, xGooseFEM::Dynamics::Geometry, solve_V         ); }
  xt::xtensor<double,2> u()      const             override { PYBIND11_OVERLOAD_PURE( Arr2, xGooseFEM::Dynamics::Geometry, u               ); }
  xt::xtensor<double,2> v()      const             override { PYBIND11_OVERLOAD_PURE( Arr2, xGooseFEM::Dynamics::Geometry, v               ); }
  xt::xtensor<double,2> a()      const             override { PYBIND11_OVERLOAD_PURE( Arr2, xGooseFEM::Dynamics::Geometry, a               ); }
  xt::xtensor<double,1> dofs_u() const             override { PYBIND11_OVERLOAD_PURE( Arr1, xGooseFEM::Dynamics::Geometry, dofs_u          ); }
  xt::xtensor<double,1> dofs_v() const             override { PYBIND11_OVERLOAD_PURE( Arr1, xGooseFEM::Dynamics::Geometry, dofs_v          ); }
  xt::xtensor<double,1> dofs_a() const             override { PYBIND11_OVERLOAD_PURE( Arr1, xGooseFEM::Dynamics::Geometry, dofs_a          ); }

  void set_u(const xt::xtensor<double,2> &nodevec) override { PYBIND11_OVERLOAD_PURE( void, xGooseFEM::Dynamics::Geometry, set_u  , nodevec); }
  void set_u(const xt::xtensor<double,1> &dofval ) override { PYBIND11_OVERLOAD_PURE( void, xGooseFEM::Dynamics::Geometry, set_u  , dofval ); }
  void set_v(const xt::xtensor<double,1> &dofval ) override { PYBIND11_OVERLOAD_PURE( void, xGooseFEM::Dynamics::Geometry, set_v  , dofval ); }
  void set_a(const xt::xtensor<double,1> &dofval ) override { PYBIND11_OVERLOAD_PURE( void, xGooseFEM::Dynamics::Geometry, set_a  , dofval ); }
};

// =========================================== GooseFEM ============================================

PYBIND11_MODULE(xGooseFEM, m) {

m.doc() = "Some simple finite element meshes and operations";

// ======================================== GooseFEM.Vector ========================================

py::class_<xGooseFEM::Vector>(m, "Vector")

  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &                               >(), "Switch between dofval/nodevec/elemvec", py::arg("conn"), py::arg("dofs"))
  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(), "Switch between dofval/nodevec/elemvec", py::arg("conn"), py::arg("dofs"), py::arg("iip"))

  .def("nelem", &M::Vector::nelem, "Return number of element")
  .def("nne"  , &M::Vector::nne  , "Return number of nodes per element")
  .def("nnode", &M::Vector::nnode, "Return number of nodes")
  .def("ndim" , &M::Vector::ndim , "Return number of dimensions")
  .def("ndof" , &M::Vector::ndof , "Return number of degrees-of-freedom")
  .def("nnu"  , &M::Vector::nnu  , "Return number of unknown degrees-of-freedom")
  .def("nnp"  , &M::Vector::nnp  , "Return number of prescribed degrees-of-freedom")

  .def("dofs" , &M::Vector::dofs , "Return degrees-of-freedom")
  .def("iiu"  , &M::Vector::iiu  , "Return unknown degrees-of-freedom")
  .def("iip"  , &M::Vector::iip  , "Return prescribed degrees-of-freedom")

  .def("asDofs"   , py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::Vector::asDofs   , py::const_), "Set 'dofval" , py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asDofs"   , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::Vector::asDofs   , py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs"   , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::Vector::asDofs   , py::const_), "Set 'dofval" , py::arg("elemvec"))
  .def("asDofs_u" , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::Vector::asDofs_u , py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs_u" , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::Vector::asDofs_u , py::const_), "Set 'dofval" , py::arg("elemvec"))
  .def("asDofs_p" , py::overload_cast<const xt::xtensor<double,2>&                             >(&M::Vector::asDofs_p , py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs_p" , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::Vector::asDofs_p , py::const_), "Set 'dofval" , py::arg("elemvec"))

  .def("asNode"   , py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::Vector::asNode   , py::const_), "Set 'nodevec", py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asNode"   , py::overload_cast<const xt::xtensor<double,1>&                             >(&M::Vector::asNode   , py::const_), "Set 'nodevec", py::arg("dofval"))
  .def("asNode"   , py::overload_cast<const xt::xtensor<double,3>&                             >(&M::Vector::asNode   , py::const_), "Set 'nodevec", py::arg("elemvec"))

  .def("asElement", py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&M::Vector::asElement, py::const_), "Set 'elemvec", py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asElement", py::overload_cast<const xt::xtensor<double,1>&                             >(&M::Vector::asElement, py::const_), "Set 'elemvec", py::arg("dofval"))
  .def("asElement", py::overload_cast<const xt::xtensor<double,2>&                             >(&M::Vector::asElement, py::const_), "Set 'elemvec", py::arg("nodevec"))

  .def("assembleDofs"  , py::overload_cast<const xt::xtensor<double,2>&>(&M::Vector::assembleDofs  , py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs"  , py::overload_cast<const xt::xtensor<double,3>&>(&M::Vector::assembleDofs  , py::const_), "Assemble 'dofval'" , py::arg("elemvec"))
  .def("assembleDofs_u", py::overload_cast<const xt::xtensor<double,2>&>(&M::Vector::assembleDofs_u, py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs_u", py::overload_cast<const xt::xtensor<double,3>&>(&M::Vector::assembleDofs_u, py::const_), "Assemble 'dofval'" , py::arg("elemvec"))
  .def("assembleDofs_p", py::overload_cast<const xt::xtensor<double,2>&>(&M::Vector::assembleDofs_p, py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs_p", py::overload_cast<const xt::xtensor<double,3>&>(&M::Vector::assembleDofs_p, py::const_), "Assemble 'dofval'" , py::arg("elemvec"))

  .def("assembleNode"  , py::overload_cast<const xt::xtensor<double,3>&>(&M::Vector::assembleNode  , py::const_), "Assemble 'nodevec'", py::arg("elemvec"))

  .def("__repr__", [](const xGooseFEM::Vector &){ return "<GooseFEM.Vector>"; });

// ==================================== GooseFEM.MatrixDiagonal ====================================

py::class_<xGooseFEM::MatrixDiagonal>(m, "MatrixDiagonal")

  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &                               >(), "Diagonal matrix", py::arg("conn"), py::arg("dofs"))
  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(), "Diagonal matrix", py::arg("conn"), py::arg("dofs"), py::arg("iip"))

  .def("nelem", &M::MatrixDiagonal::nelem, "Return number of element")
  .def("nne"  , &M::MatrixDiagonal::nne  , "Return number of nodes per element")
  .def("nnode", &M::MatrixDiagonal::nnode, "Return number of nodes")
  .def("ndim" , &M::MatrixDiagonal::ndim , "Return number of dimensions")
  .def("ndof" , &M::MatrixDiagonal::ndof , "Return number of degrees-of-freedom")
  .def("nnu"  , &M::MatrixDiagonal::nnu  , "Return number of unknown degrees-of-freedom")
  .def("nnp"  , &M::MatrixDiagonal::nnp  , "Return number of prescribed degrees-of-freedom")

  .def("dofs" , &M::MatrixDiagonal::dofs , "Return degrees-of-freedom")
  .def("iiu"  , &M::MatrixDiagonal::iiu  , "Return unknown degrees-of-freedom")
  .def("iip"  , &M::MatrixDiagonal::iip  , "Return prescribed degrees-of-freedom")

  .def("dot"  , &M::MatrixDiagonal::dot  , "Dot product 'b_i = A_ij * x_j", py::arg("x"))
  .def("dot_u", &M::MatrixDiagonal::dot_u, "Dot product 'b_i = A_ij * x_j", py::arg("x_u"), py::arg("x_p"))
  .def("dot_p", &M::MatrixDiagonal::dot_p, "Dot product 'b_i = A_ij * x_j", py::arg("x_u"), py::arg("x_p"))

  .def("assemble", &M::MatrixDiagonal::assemble, "Assemble from 'elemmat", py::arg("elemmat"))

  .def("solve"  , py::overload_cast<const xt::xtensor<double,1>&                              >(&M::MatrixDiagonal::solve  ), "Solve", py::arg("b"  )                )
  .def("solve"  , py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(&M::MatrixDiagonal::solve  ), "Solve", py::arg("b"  ), py::arg("x_p"))
  .def("solve_u", py::overload_cast<const xt::xtensor<double,1>&, const xt::xtensor<double,1>&>(&M::MatrixDiagonal::solve_u), "Solve", py::arg("b_u"), py::arg("x_p"))

  .def("asDiagonal", &M::MatrixDiagonal::asDiagonal, "Return as diagonal matrix (column)")

  .def("__repr__", [](const xGooseFEM::MatrixDiagonal &){ return "<GooseFEM.MatrixDiagonal>"; });

// ======================================= GooseFEM.Dynamics =======================================

py::module mDynamics = m.def_submodule("Dynamics", "Solve routines for dynamic FEM");

// -------------------------------------------------------------------------------------------------

mDynamics.def("Verlet"        , &xGooseFEM::Dynamics::Verlet        , "Verlet time integration"         , py::arg("geometry"), py::arg("dt"), py::arg("nstep")=1);
mDynamics.def("velocityVerlet", &xGooseFEM::Dynamics::velocityVerlet, "Velocity-Verlet time integration", py::arg("geometry"), py::arg("dt"), py::arg("nstep")=1);

// -------------------------------------------------------------------------------------------------

py::class_<xGooseFEM::Dynamics::Geometry, PyGeometry>(mDynamics, "Geometry")

  .def(py::init<>())

  .def("solve_A"   , &xGooseFEM::Dynamics::Geometry::solve_A, "Solve for accelerations (dofval)" )
  .def("solve_V"   , &xGooseFEM::Dynamics::Geometry::solve_V, "Solve for velocities    (dofval)" )

  .def("u"         , &xGooseFEM::Dynamics::Geometry::u      , "Return displacements (nodevec)")
  .def("v"         , &xGooseFEM::Dynamics::Geometry::v      , "Return velocities    (nodevec)")
  .def("a"         , &xGooseFEM::Dynamics::Geometry::a      , "Return accelerations (nodevec)")

  .def("dofs_u"    , &xGooseFEM::Dynamics::Geometry::dofs_u , "Return displacements (dofval)" )
  .def("dofs_v"    , &xGooseFEM::Dynamics::Geometry::dofs_v , "Return velocities    (dofval)" )
  .def("dofs_a"    , &xGooseFEM::Dynamics::Geometry::dofs_a , "Return accelerations (dofval)" )

  .def("set_u"     , py::overload_cast<const xt::xtensor<double,1> &>(&xGooseFEM::Dynamics::Geometry::set_u), "Overwrite displacements", py::arg("dofval" ))
  .def("set_u"     , py::overload_cast<const xt::xtensor<double,2> &>(&xGooseFEM::Dynamics::Geometry::set_u), "Overwrite displacements", py::arg("nodevec"))
  .def("set_v"     ,                                                  &xGooseFEM::Dynamics::Geometry::set_v , "Overwrite velocities"   , py::arg("nodevec"))
  .def("set_a"     ,                                                  &xGooseFEM::Dynamics::Geometry::set_a , "Overwrite accelerations", py::arg("nodevec"))

  .def("__repr__", [](const xGooseFEM::Dynamics::Geometry &){ return "<GooseDEM.Dynamics.Geometry>"; });

// ======================================= GooseFEM.Element ========================================

py::module mElement = m.def_submodule("Element", "Generic element routines");

// -------------------------------------------------------------------------------------------------

mElement.def("asElementVector"      , &xGooseFEM::Element::asElementVector   , "Covert to 'elemvec'", py::arg("conn"), py::arg("nodevec"));
mElement.def("assembleElementVector", &xGooseFEM::Element::assembleNodeVector, "Assemble 'nodevec'" , py::arg("conn"), py::arg("elemvec"));

// ==================================== GooseFEM.Element.Quad4 =====================================

{

py::module sm = mElement.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

// -------------------------------------------------------------------------------------------------

py::class_<xGooseFEM::Element::Quad4::Quadrature>(sm, "Quadrature")

  .def(py::init<const xt::xtensor<double,3> &>(), "Quadrature", py::arg("x"))

  .def(py::init<const xt::xtensor<double,3> &, const xt::xtensor<double,2> &, const xt::xtensor<double,1> &>(), "Quadrature", py::arg("x"), py::arg("xi"), py::arg("w"))

  .def("nelem", &xGooseFEM::Element::Quad4::Quadrature::nelem, "Return number of elements"          )
  .def("nne"  , &xGooseFEM::Element::Quad4::Quadrature::nne  , "Return number of nodes per element" )
  .def("ndim" , &xGooseFEM::Element::Quad4::Quadrature::ndim , "Return number of dimensions"        )
  .def("nip"  , &xGooseFEM::Element::Quad4::Quadrature::nip  , "Return number of integration points")

  .def("dV"      , py::overload_cast<>(&xGooseFEM::Element::Quad4::Quadrature::dV      , py::const_), "Integration point volume (qscalar)")
  .def("dVtensor", py::overload_cast<>(&xGooseFEM::Element::Quad4::Quadrature::dVtensor, py::const_), "Integration point volume (qtensor)")

  .def("gradN_vector"   , py::overload_cast<const xt::xtensor<double,3> &>(&xGooseFEM::Element::Quad4::Quadrature::gradN_vector   , py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))
  .def("gradN_vector_T" , py::overload_cast<const xt::xtensor<double,3> &>(&xGooseFEM::Element::Quad4::Quadrature::gradN_vector_T , py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))
  .def("symGradN_vector", py::overload_cast<const xt::xtensor<double,3> &>(&xGooseFEM::Element::Quad4::Quadrature::symGradN_vector, py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))

  .def("int_N_scalar_NT_dV"      , py::overload_cast<const xt::xtensor<double,2> &>(&xGooseFEM::Element::Quad4::Quadrature::int_N_scalar_NT_dV      , py::const_), "Integration, returns 'elemmat'", py::arg("qscalar"))
  .def("int_gradN_dot_tensor2_dV", py::overload_cast<const xt::xtensor<double,4> &>(&xGooseFEM::Element::Quad4::Quadrature::int_gradN_dot_tensor2_dV, py::const_), "Integration, returns 'elemvec'", py::arg("qtensor"))

  .def("__repr__", [](const xGooseFEM::Element::Quad4::Quadrature &){ return "<GooseFEM.Element.Quad4.Quadrature>"; });

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Gauss", "Gauss quadrature");

ssm.def("nip", &xGooseFEM::Element::Quad4::Gauss::nip, "Return number of integration point"  );
ssm.def("xi" , &xGooseFEM::Element::Quad4::Gauss::xi , "Return integration point coordinates");
ssm.def("w"  , &xGooseFEM::Element::Quad4::Gauss::w  , "Return integration point weights"    );

}

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Nodal", "Nodal quadrature");

ssm.def("nip", &xGooseFEM::Element::Quad4::Nodal::nip, "Return number of integration point"  );
ssm.def("xi" , &xGooseFEM::Element::Quad4::Nodal::xi , "Return integration point coordinates");
ssm.def("w"  , &xGooseFEM::Element::Quad4::Nodal::w  , "Return integration point weights"    );

}

// -------------------------------------------------------------------------------------------------

}

// ===================================== GooseFEM.Element.Hex8 =====================================

{

py::module sm = mElement.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

// -------------------------------------------------------------------------------------------------

py::class_<xGooseFEM::Element::Hex8::Quadrature>(sm, "Quadrature")

  .def(py::init<const xt::xtensor<double,3> &>(), "Quadrature", py::arg("x"))

  .def(py::init<const xt::xtensor<double,3> &, const xt::xtensor<double,2> &, const xt::xtensor<double,1> &>(), "Quadrature", py::arg("x"), py::arg("xi"), py::arg("w"))

  .def("nelem", &xGooseFEM::Element::Hex8::Quadrature::nelem, "Return number of elements"          )
  .def("nne"  , &xGooseFEM::Element::Hex8::Quadrature::nne  , "Return number of nodes per element" )
  .def("ndim" , &xGooseFEM::Element::Hex8::Quadrature::ndim , "Return number of dimensions"        )
  .def("nip"  , &xGooseFEM::Element::Hex8::Quadrature::nip  , "Return number of integration points")

  .def("dV"      , py::overload_cast<>(&xGooseFEM::Element::Hex8::Quadrature::dV      , py::const_), "Integration point volume (qscalar)")
  .def("dVtensor", py::overload_cast<>(&xGooseFEM::Element::Hex8::Quadrature::dVtensor, py::const_), "Integration point volume (qtensor)")

  .def("gradN_vector"   , py::overload_cast<const xt::xtensor<double,3> &>(&xGooseFEM::Element::Hex8::Quadrature::gradN_vector   , py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))
  .def("gradN_vector_T" , py::overload_cast<const xt::xtensor<double,3> &>(&xGooseFEM::Element::Hex8::Quadrature::gradN_vector_T , py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))
  .def("symGradN_vector", py::overload_cast<const xt::xtensor<double,3> &>(&xGooseFEM::Element::Hex8::Quadrature::symGradN_vector, py::const_), "Dyadic product, returns 'qtensor'", py::arg("elemvec"))

  .def("int_N_scalar_NT_dV"      , py::overload_cast<const xt::xtensor<double,2> &>(&xGooseFEM::Element::Hex8::Quadrature::int_N_scalar_NT_dV      , py::const_), "Integration, returns 'elemmat'", py::arg("qscalar"))
  .def("int_gradN_dot_tensor2_dV", py::overload_cast<const xt::xtensor<double,4> &>(&xGooseFEM::Element::Hex8::Quadrature::int_gradN_dot_tensor2_dV, py::const_), "Integration, returns 'elemvec'", py::arg("qtensor"))

  .def("__repr__", [](const xGooseFEM::Element::Hex8::Quadrature &){ return "<GooseFEM.Element.Hex8.Quadrature>"; });

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Gauss", "Gauss quadrature");

ssm.def("nip", &xGooseFEM::Element::Hex8::Gauss::nip, "Return number of integration point"  );
ssm.def("xi" , &xGooseFEM::Element::Hex8::Gauss::xi , "Return integration point coordinates");
ssm.def("w"  , &xGooseFEM::Element::Hex8::Gauss::w  , "Return integration point weights"    );

}

// -------------------------------------------------------------------------------------------------

{

py::module ssm = sm.def_submodule("Nodal", "Nodal quadrature");

ssm.def("nip", &xGooseFEM::Element::Hex8::Nodal::nip, "Return number of integration point"  );
ssm.def("xi" , &xGooseFEM::Element::Hex8::Nodal::xi , "Return integration point coordinates");
ssm.def("w"  , &xGooseFEM::Element::Hex8::Nodal::w  , "Return integration point weights"    );

}

// -------------------------------------------------------------------------------------------------

}

// ========================================= GooseFEM.Mesh =========================================

py::module mMesh = m.def_submodule("Mesh", "Generic mesh routines");

// -------------------------------------------------------------------------------------------------

mMesh.def("dofs", &xGooseFEM::Mesh::dofs, "List with DOF-numbers (in sequential order)", py::arg("nnode"), py::arg("ndim"));

mMesh.def("renumber", &xGooseFEM::Mesh::renumber, "Renumber DOF-list to use the lowest possible index", py::arg("dofs"));

mMesh.def("renumber_index", &xGooseFEM::Mesh::renumber_index, "Index-list to renumber", py::arg("dofs"));

mMesh.def("reorder", &xGooseFEM::Mesh::reorder, "Renumber DOF-list to begin or end with 'idx'", py::arg("dofs"), py::arg("idx"), py::arg("location")="end");

mMesh.def("reorder_index", &xGooseFEM::Mesh::reorder_index, "Index-list to reorder", py::arg("dofs"), py::arg("idx"), py::arg("location")="end");

mMesh.def("coordination", &xGooseFEM::Mesh::coordination, "Coordination number of each node", py::arg("conn"));

mMesh.def("elem2node", &xGooseFEM::Mesh::elem2node, "Elements connect to each node", py::arg("conn"));

// ====================================== GooseFEM.Mesh.Hex8 =======================================

{

py::module sm = mMesh.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

// -------------------------------------------------------------------------------------------------

py::class_<xGooseFEM::Mesh::Hex8::Regular>(sm, "Regular")

  .def(py::init<size_t,size_t,size_t,double>(), "mesh with nx*ny*nz 'pixels' and edge size h", py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("h")=1.)

  .def("nelem"                      , &xGooseFEM::Mesh::Hex8::Regular::nelem                      )
  .def("nnode"                      , &xGooseFEM::Mesh::Hex8::Regular::nnode                      )
  .def("nne"                        , &xGooseFEM::Mesh::Hex8::Regular::nne                        )
  .def("ndim"                       , &xGooseFEM::Mesh::Hex8::Regular::ndim                       )

  .def("coor"                       , &xGooseFEM::Mesh::Hex8::Regular::coor                       )
  .def("conn"                       , &xGooseFEM::Mesh::Hex8::Regular::conn                       )

  .def("nodesFront"                 , &xGooseFEM::Mesh::Hex8::Regular::nodesFront                 )
  .def("nodesBack"                  , &xGooseFEM::Mesh::Hex8::Regular::nodesBack                  )
  .def("nodesLeft"                  , &xGooseFEM::Mesh::Hex8::Regular::nodesLeft                  )
  .def("nodesRight"                 , &xGooseFEM::Mesh::Hex8::Regular::nodesRight                 )
  .def("nodesBottom"                , &xGooseFEM::Mesh::Hex8::Regular::nodesBottom                )
  .def("nodesTop"                   , &xGooseFEM::Mesh::Hex8::Regular::nodesTop                   )

  .def("nodesFrontFace"             , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontFace             )
  .def("nodesBackFace"              , &xGooseFEM::Mesh::Hex8::Regular::nodesBackFace              )
  .def("nodesLeftFace"              , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftFace              )
  .def("nodesRightFace"             , &xGooseFEM::Mesh::Hex8::Regular::nodesRightFace             )
  .def("nodesBottomFace"            , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomFace            )
  .def("nodesTopFace"               , &xGooseFEM::Mesh::Hex8::Regular::nodesTopFace               )

  .def("nodesFrontBottomEdge"       , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontBottomEdge       )
  .def("nodesFrontTopEdge"          , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontTopEdge          )
  .def("nodesFrontLeftEdge"         , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontLeftEdge         )
  .def("nodesFrontRightEdge"        , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontRightEdge        )
  .def("nodesBackBottomEdge"        , &xGooseFEM::Mesh::Hex8::Regular::nodesBackBottomEdge        )
  .def("nodesBackTopEdge"           , &xGooseFEM::Mesh::Hex8::Regular::nodesBackTopEdge           )
  .def("nodesBackLeftEdge"          , &xGooseFEM::Mesh::Hex8::Regular::nodesBackLeftEdge          )
  .def("nodesBackRightEdge"         , &xGooseFEM::Mesh::Hex8::Regular::nodesBackRightEdge         )
  .def("nodesBottomLeftEdge"        , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomLeftEdge        )
  .def("nodesBottomRightEdge"       , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomRightEdge       )
  .def("nodesTopLeftEdge"           , &xGooseFEM::Mesh::Hex8::Regular::nodesTopLeftEdge           )
  .def("nodesTopRightEdge"          , &xGooseFEM::Mesh::Hex8::Regular::nodesTopRightEdge          )

  .def("nodesBottomFrontEdge"       , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomFrontEdge       )
  .def("nodesBottomBackEdge"        , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomBackEdge        )
  .def("nodesTopFrontEdge"          , &xGooseFEM::Mesh::Hex8::Regular::nodesTopFrontEdge          )
  .def("nodesTopBackEdge"           , &xGooseFEM::Mesh::Hex8::Regular::nodesTopBackEdge           )
  .def("nodesLeftBottomEdge"        , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftBottomEdge        )
  .def("nodesLeftFrontEdge"         , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftFrontEdge         )
  .def("nodesLeftBackEdge"          , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftBackEdge          )
  .def("nodesLeftTopEdge"           , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftTopEdge           )
  .def("nodesRightBottomEdge"       , &xGooseFEM::Mesh::Hex8::Regular::nodesRightBottomEdge       )
  .def("nodesRightTopEdge"          , &xGooseFEM::Mesh::Hex8::Regular::nodesRightTopEdge          )
  .def("nodesRightFrontEdge"        , &xGooseFEM::Mesh::Hex8::Regular::nodesRightFrontEdge        )
  .def("nodesRightBackEdge"         , &xGooseFEM::Mesh::Hex8::Regular::nodesRightBackEdge         )

  .def("nodesFrontBottomOpenEdge"   , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontBottomOpenEdge   )
  .def("nodesFrontTopOpenEdge"      , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontTopOpenEdge      )
  .def("nodesFrontLeftOpenEdge"     , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontLeftOpenEdge     )
  .def("nodesFrontRightOpenEdge"    , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontRightOpenEdge    )
  .def("nodesBackBottomOpenEdge"    , &xGooseFEM::Mesh::Hex8::Regular::nodesBackBottomOpenEdge    )
  .def("nodesBackTopOpenEdge"       , &xGooseFEM::Mesh::Hex8::Regular::nodesBackTopOpenEdge       )
  .def("nodesBackLeftOpenEdge"      , &xGooseFEM::Mesh::Hex8::Regular::nodesBackLeftOpenEdge      )
  .def("nodesBackRightOpenEdge"     , &xGooseFEM::Mesh::Hex8::Regular::nodesBackRightOpenEdge     )
  .def("nodesBottomLeftOpenEdge"    , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomLeftOpenEdge    )
  .def("nodesBottomRightOpenEdge"   , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomRightOpenEdge   )
  .def("nodesTopLeftOpenEdge"       , &xGooseFEM::Mesh::Hex8::Regular::nodesTopLeftOpenEdge       )
  .def("nodesTopRightOpenEdge"      , &xGooseFEM::Mesh::Hex8::Regular::nodesTopRightOpenEdge      )

  .def("nodesBottomFrontOpenEdge"   , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomFrontOpenEdge   )
  .def("nodesBottomBackOpenEdge"    , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomBackOpenEdge    )
  .def("nodesTopFrontOpenEdge"      , &xGooseFEM::Mesh::Hex8::Regular::nodesTopFrontOpenEdge      )
  .def("nodesTopBackOpenEdge"       , &xGooseFEM::Mesh::Hex8::Regular::nodesTopBackOpenEdge       )
  .def("nodesLeftBottomOpenEdge"    , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftBottomOpenEdge    )
  .def("nodesLeftFrontOpenEdge"     , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftFrontOpenEdge     )
  .def("nodesLeftBackOpenEdge"      , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftBackOpenEdge      )
  .def("nodesLeftTopOpenEdge"       , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftTopOpenEdge       )
  .def("nodesRightBottomOpenEdge"   , &xGooseFEM::Mesh::Hex8::Regular::nodesRightBottomOpenEdge   )
  .def("nodesRightTopOpenEdge"      , &xGooseFEM::Mesh::Hex8::Regular::nodesRightTopOpenEdge      )
  .def("nodesRightFrontOpenEdge"    , &xGooseFEM::Mesh::Hex8::Regular::nodesRightFrontOpenEdge    )
  .def("nodesRightBackOpenEdge"     , &xGooseFEM::Mesh::Hex8::Regular::nodesRightBackOpenEdge     )

  .def("nodesFrontBottomLeftCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontBottomLeftCorner )
  .def("nodesFrontBottomRightCorner", &xGooseFEM::Mesh::Hex8::Regular::nodesFrontBottomRightCorner)
  .def("nodesFrontTopLeftCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontTopLeftCorner    )
  .def("nodesFrontTopRightCorner"   , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontTopRightCorner   )
  .def("nodesBackBottomLeftCorner"  , &xGooseFEM::Mesh::Hex8::Regular::nodesBackBottomLeftCorner  )
  .def("nodesBackBottomRightCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesBackBottomRightCorner )
  .def("nodesBackTopLeftCorner"     , &xGooseFEM::Mesh::Hex8::Regular::nodesBackTopLeftCorner     )
  .def("nodesBackTopRightCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesBackTopRightCorner    )

  .def("nodesFrontLeftBottomCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontLeftBottomCorner )
  .def("nodesBottomFrontLeftCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomFrontLeftCorner )
  .def("nodesBottomLeftFrontCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomLeftFrontCorner )
  .def("nodesLeftFrontBottomCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftFrontBottomCorner )
  .def("nodesLeftBottomFrontCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftBottomFrontCorner )
  .def("nodesFrontRightBottomCorner", &xGooseFEM::Mesh::Hex8::Regular::nodesFrontRightBottomCorner)
  .def("nodesBottomFrontRightCorner", &xGooseFEM::Mesh::Hex8::Regular::nodesBottomFrontRightCorner)
  .def("nodesBottomRightFrontCorner", &xGooseFEM::Mesh::Hex8::Regular::nodesBottomRightFrontCorner)
  .def("nodesRightFrontBottomCorner", &xGooseFEM::Mesh::Hex8::Regular::nodesRightFrontBottomCorner)
  .def("nodesRightBottomFrontCorner", &xGooseFEM::Mesh::Hex8::Regular::nodesRightBottomFrontCorner)
  .def("nodesFrontLeftTopCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontLeftTopCorner    )
  .def("nodesTopFrontLeftCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesTopFrontLeftCorner    )
  .def("nodesTopLeftFrontCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesTopLeftFrontCorner    )
  .def("nodesLeftFrontTopCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftFrontTopCorner    )
  .def("nodesLeftTopFrontCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftTopFrontCorner    )
  .def("nodesFrontRightTopCorner"   , &xGooseFEM::Mesh::Hex8::Regular::nodesFrontRightTopCorner   )
  .def("nodesTopFrontRightCorner"   , &xGooseFEM::Mesh::Hex8::Regular::nodesTopFrontRightCorner   )
  .def("nodesTopRightFrontCorner"   , &xGooseFEM::Mesh::Hex8::Regular::nodesTopRightFrontCorner   )
  .def("nodesRightFrontTopCorner"   , &xGooseFEM::Mesh::Hex8::Regular::nodesRightFrontTopCorner   )
  .def("nodesRightTopFrontCorner"   , &xGooseFEM::Mesh::Hex8::Regular::nodesRightTopFrontCorner   )
  .def("nodesBackLeftBottomCorner"  , &xGooseFEM::Mesh::Hex8::Regular::nodesBackLeftBottomCorner  )
  .def("nodesBottomBackLeftCorner"  , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomBackLeftCorner  )
  .def("nodesBottomLeftBackCorner"  , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomLeftBackCorner  )
  .def("nodesLeftBackBottomCorner"  , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftBackBottomCorner  )
  .def("nodesLeftBottomBackCorner"  , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftBottomBackCorner  )
  .def("nodesBackRightBottomCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesBackRightBottomCorner )
  .def("nodesBottomBackRightCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomBackRightCorner )
  .def("nodesBottomRightBackCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesBottomRightBackCorner )
  .def("nodesRightBackBottomCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesRightBackBottomCorner )
  .def("nodesRightBottomBackCorner" , &xGooseFEM::Mesh::Hex8::Regular::nodesRightBottomBackCorner )
  .def("nodesBackLeftTopCorner"     , &xGooseFEM::Mesh::Hex8::Regular::nodesBackLeftTopCorner     )
  .def("nodesTopBackLeftCorner"     , &xGooseFEM::Mesh::Hex8::Regular::nodesTopBackLeftCorner     )
  .def("nodesTopLeftBackCorner"     , &xGooseFEM::Mesh::Hex8::Regular::nodesTopLeftBackCorner     )
  .def("nodesLeftBackTopCorner"     , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftBackTopCorner     )
  .def("nodesLeftTopBackCorner"     , &xGooseFEM::Mesh::Hex8::Regular::nodesLeftTopBackCorner     )
  .def("nodesBackRightTopCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesBackRightTopCorner    )
  .def("nodesTopBackRightCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesTopBackRightCorner    )
  .def("nodesTopRightBackCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesTopRightBackCorner    )
  .def("nodesRightBackTopCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesRightBackTopCorner    )
  .def("nodesRightTopBackCorner"    , &xGooseFEM::Mesh::Hex8::Regular::nodesRightTopBackCorner    )

  .def("nodesPeriodic"              , &xGooseFEM::Mesh::Hex8::Regular::nodesPeriodic              )
  .def("nodesOrigin"                , &xGooseFEM::Mesh::Hex8::Regular::nodesOrigin                )
  .def("dofs"                       , &xGooseFEM::Mesh::Hex8::Regular::dofs                       )
  .def("dofsPeriodic"               , &xGooseFEM::Mesh::Hex8::Regular::dofsPeriodic               )

  .def("__repr__", [](const xGooseFEM::Mesh::Hex8::Regular &){ return "<GooseFEM.Mesh.Hex8.Regular>"; });

// -------------------------------------------------------------------------------------------------

py::class_<xGooseFEM::Mesh::Hex8::FineLayer>(sm, "FineLayer")

  .def(py::init<size_t,size_t,size_t,double,size_t>(), "mesh with nx*ny*nz 'pixels' and edge size h", py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("h")=1., py::arg("nfine")=1)

  .def("nelem"                      , &xGooseFEM::Mesh::Hex8::FineLayer::nelem                      )
  .def("nnode"                      , &xGooseFEM::Mesh::Hex8::FineLayer::nnode                      )
  .def("nne"                        , &xGooseFEM::Mesh::Hex8::FineLayer::nne                        )
  .def("ndim"                       , &xGooseFEM::Mesh::Hex8::FineLayer::ndim                       )
  .def("shape"                      , &xGooseFEM::Mesh::Hex8::FineLayer::shape                      )

  .def("coor"                       , &xGooseFEM::Mesh::Hex8::FineLayer::coor                       )
  .def("conn"                       , &xGooseFEM::Mesh::Hex8::FineLayer::conn                       )

  .def("elementsMiddleLayer"        , &xGooseFEM::Mesh::Hex8::FineLayer::elementsMiddleLayer        )

  .def("nodesFront"                 , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFront                 )
  .def("nodesBack"                  , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBack                  )
  .def("nodesLeft"                  , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeft                  )
  .def("nodesRight"                 , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRight                 )
  .def("nodesBottom"                , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottom                )
  .def("nodesTop"                   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTop                   )

  .def("nodesFrontFace"             , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontFace             )
  .def("nodesBackFace"              , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackFace              )
  .def("nodesLeftFace"              , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftFace              )
  .def("nodesRightFace"             , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightFace             )
  .def("nodesBottomFace"            , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomFace            )
  .def("nodesTopFace"               , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopFace               )

  .def("nodesFrontBottomEdge"       , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomEdge       )
  .def("nodesFrontTopEdge"          , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopEdge          )
  .def("nodesFrontLeftEdge"         , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftEdge         )
  .def("nodesFrontRightEdge"        , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightEdge        )
  .def("nodesBackBottomEdge"        , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomEdge        )
  .def("nodesBackTopEdge"           , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackTopEdge           )
  .def("nodesBackLeftEdge"          , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftEdge          )
  .def("nodesBackRightEdge"         , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackRightEdge         )
  .def("nodesBottomLeftEdge"        , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftEdge        )
  .def("nodesBottomRightEdge"       , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightEdge       )
  .def("nodesTopLeftEdge"           , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftEdge           )
  .def("nodesTopRightEdge"          , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopRightEdge          )

  .def("nodesBottomFrontEdge"       , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontEdge       )
  .def("nodesBottomBackEdge"        , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackEdge        )
  .def("nodesTopFrontEdge"          , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontEdge          )
  .def("nodesTopBackEdge"           , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopBackEdge           )
  .def("nodesLeftBottomEdge"        , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomEdge        )
  .def("nodesLeftFrontEdge"         , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontEdge         )
  .def("nodesLeftBackEdge"          , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackEdge          )
  .def("nodesLeftTopEdge"           , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopEdge           )
  .def("nodesRightBottomEdge"       , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomEdge       )
  .def("nodesRightTopEdge"          , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightTopEdge          )
  .def("nodesRightFrontEdge"        , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontEdge        )
  .def("nodesRightBackEdge"         , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightBackEdge         )

  .def("nodesFrontBottomOpenEdge"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomOpenEdge   )
  .def("nodesFrontTopOpenEdge"      , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopOpenEdge      )
  .def("nodesFrontLeftOpenEdge"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftOpenEdge     )
  .def("nodesFrontRightOpenEdge"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightOpenEdge    )
  .def("nodesBackBottomOpenEdge"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomOpenEdge    )
  .def("nodesBackTopOpenEdge"       , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackTopOpenEdge       )
  .def("nodesBackLeftOpenEdge"      , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftOpenEdge      )
  .def("nodesBackRightOpenEdge"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackRightOpenEdge     )
  .def("nodesBottomLeftOpenEdge"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftOpenEdge    )
  .def("nodesBottomRightOpenEdge"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightOpenEdge   )
  .def("nodesTopLeftOpenEdge"       , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftOpenEdge       )
  .def("nodesTopRightOpenEdge"      , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopRightOpenEdge      )

  .def("nodesBottomFrontOpenEdge"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontOpenEdge   )
  .def("nodesBottomBackOpenEdge"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackOpenEdge    )
  .def("nodesTopFrontOpenEdge"      , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontOpenEdge      )
  .def("nodesTopBackOpenEdge"       , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopBackOpenEdge       )
  .def("nodesLeftBottomOpenEdge"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomOpenEdge    )
  .def("nodesLeftFrontOpenEdge"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontOpenEdge     )
  .def("nodesLeftBackOpenEdge"      , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackOpenEdge      )
  .def("nodesLeftTopOpenEdge"       , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopOpenEdge       )
  .def("nodesRightBottomOpenEdge"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomOpenEdge   )
  .def("nodesRightTopOpenEdge"      , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightTopOpenEdge      )
  .def("nodesRightFrontOpenEdge"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontOpenEdge    )
  .def("nodesRightBackOpenEdge"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightBackOpenEdge     )

  .def("nodesFrontBottomLeftCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomLeftCorner )
  .def("nodesFrontBottomRightCorner", &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontBottomRightCorner)
  .def("nodesFrontTopLeftCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopLeftCorner    )
  .def("nodesFrontTopRightCorner"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontTopRightCorner   )
  .def("nodesBackBottomLeftCorner"  , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomLeftCorner  )
  .def("nodesBackBottomRightCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackBottomRightCorner )
  .def("nodesBackTopLeftCorner"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackTopLeftCorner     )
  .def("nodesBackTopRightCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackTopRightCorner    )

  .def("nodesFrontLeftBottomCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftBottomCorner )
  .def("nodesBottomFrontLeftCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontLeftCorner )
  .def("nodesBottomLeftFrontCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftFrontCorner )
  .def("nodesLeftFrontBottomCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontBottomCorner )
  .def("nodesLeftBottomFrontCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomFrontCorner )
  .def("nodesFrontRightBottomCorner", &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightBottomCorner)
  .def("nodesBottomFrontRightCorner", &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomFrontRightCorner)
  .def("nodesBottomRightFrontCorner", &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightFrontCorner)
  .def("nodesRightFrontBottomCorner", &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontBottomCorner)
  .def("nodesRightBottomFrontCorner", &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomFrontCorner)
  .def("nodesFrontLeftTopCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontLeftTopCorner    )
  .def("nodesTopFrontLeftCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontLeftCorner    )
  .def("nodesTopLeftFrontCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftFrontCorner    )
  .def("nodesLeftFrontTopCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftFrontTopCorner    )
  .def("nodesLeftTopFrontCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopFrontCorner    )
  .def("nodesFrontRightTopCorner"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesFrontRightTopCorner   )
  .def("nodesTopFrontRightCorner"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopFrontRightCorner   )
  .def("nodesTopRightFrontCorner"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopRightFrontCorner   )
  .def("nodesRightFrontTopCorner"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightFrontTopCorner   )
  .def("nodesRightTopFrontCorner"   , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightTopFrontCorner   )
  .def("nodesBackLeftBottomCorner"  , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftBottomCorner  )
  .def("nodesBottomBackLeftCorner"  , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackLeftCorner  )
  .def("nodesBottomLeftBackCorner"  , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomLeftBackCorner  )
  .def("nodesLeftBackBottomCorner"  , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackBottomCorner  )
  .def("nodesLeftBottomBackCorner"  , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftBottomBackCorner  )
  .def("nodesBackRightBottomCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackRightBottomCorner )
  .def("nodesBottomBackRightCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomBackRightCorner )
  .def("nodesBottomRightBackCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBottomRightBackCorner )
  .def("nodesRightBackBottomCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightBackBottomCorner )
  .def("nodesRightBottomBackCorner" , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightBottomBackCorner )
  .def("nodesBackLeftTopCorner"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackLeftTopCorner     )
  .def("nodesTopBackLeftCorner"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopBackLeftCorner     )
  .def("nodesTopLeftBackCorner"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopLeftBackCorner     )
  .def("nodesLeftBackTopCorner"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftBackTopCorner     )
  .def("nodesLeftTopBackCorner"     , &xGooseFEM::Mesh::Hex8::FineLayer::nodesLeftTopBackCorner     )
  .def("nodesBackRightTopCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesBackRightTopCorner    )
  .def("nodesTopBackRightCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopBackRightCorner    )
  .def("nodesTopRightBackCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesTopRightBackCorner    )
  .def("nodesRightBackTopCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightBackTopCorner    )
  .def("nodesRightTopBackCorner"    , &xGooseFEM::Mesh::Hex8::FineLayer::nodesRightTopBackCorner    )

  .def("nodesPeriodic"              , &xGooseFEM::Mesh::Hex8::FineLayer::nodesPeriodic              )
  .def("nodesOrigin"                , &xGooseFEM::Mesh::Hex8::FineLayer::nodesOrigin                )

  .def("dofs"                       , &xGooseFEM::Mesh::Hex8::FineLayer::dofs                       )
  .def("dofsPeriodic"               , &xGooseFEM::Mesh::Hex8::FineLayer::dofsPeriodic               )

  .def("__repr__", [](const xGooseFEM::Mesh::Hex8::FineLayer &){ return "<GooseFEM.Mesh.Hex8.FineLayer>"; });

// -------------------------------------------------------------------------------------------------

}

// ====================================== GooseFEM.Mesh.Quad4 ======================================

{

py::module sm = mMesh.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

// -------------------------------------------------------------------------------------------------

py::class_<xGooseFEM::Mesh::Quad4::Regular>(sm, "Regular")

  .def(py::init<size_t,size_t,double>(), "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'", py::arg("nx"), py::arg("ny"), py::arg("h")=1.)

  .def("coor"                  , &xGooseFEM::Mesh::Quad4::Regular::coor                  )
  .def("conn"                  , &xGooseFEM::Mesh::Quad4::Regular::conn                  )
  .def("nelem"                 , &xGooseFEM::Mesh::Quad4::Regular::nelem                 )
  .def("nnode"                 , &xGooseFEM::Mesh::Quad4::Regular::nnode                 )
  .def("nne"                   , &xGooseFEM::Mesh::Quad4::Regular::nne                   )
  .def("ndim"                  , &xGooseFEM::Mesh::Quad4::Regular::ndim                  )

  .def("nodesBottomEdge"       , &xGooseFEM::Mesh::Quad4::Regular::nodesBottomEdge       )
  .def("nodesTopEdge"          , &xGooseFEM::Mesh::Quad4::Regular::nodesTopEdge          )
  .def("nodesLeftEdge"         , &xGooseFEM::Mesh::Quad4::Regular::nodesLeftEdge         )
  .def("nodesRightEdge"        , &xGooseFEM::Mesh::Quad4::Regular::nodesRightEdge        )
  .def("nodesBottomOpenEdge"   , &xGooseFEM::Mesh::Quad4::Regular::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &xGooseFEM::Mesh::Quad4::Regular::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &xGooseFEM::Mesh::Quad4::Regular::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &xGooseFEM::Mesh::Quad4::Regular::nodesRightOpenEdge    )

  .def("nodesBottomLeftCorner" , &xGooseFEM::Mesh::Quad4::Regular::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &xGooseFEM::Mesh::Quad4::Regular::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &xGooseFEM::Mesh::Quad4::Regular::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &xGooseFEM::Mesh::Quad4::Regular::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &xGooseFEM::Mesh::Quad4::Regular::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &xGooseFEM::Mesh::Quad4::Regular::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &xGooseFEM::Mesh::Quad4::Regular::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &xGooseFEM::Mesh::Quad4::Regular::nodesRightTopCorner   )

  .def("nodesPeriodic"         , &xGooseFEM::Mesh::Quad4::Regular::nodesPeriodic         )
  .def("nodesOrigin"           , &xGooseFEM::Mesh::Quad4::Regular::nodesOrigin           )

  .def("dofs"                  , &xGooseFEM::Mesh::Quad4::Regular::dofs                  )
  .def("dofsPeriodic"          , &xGooseFEM::Mesh::Quad4::Regular::dofsPeriodic          )

  .def("__repr__", [](const xGooseFEM::Mesh::Quad4::Regular &){ return "<GooseFEM.Mesh.Quad4.Regular>"; });

// -------------------------------------------------------------------------------------------------

py::class_<xGooseFEM::Mesh::Quad4::FineLayer>(sm, "FineLayer")

  .def(
    py::init<size_t,size_t,double,size_t>(),
    "FineLayer mesh: 'nx' pixels in horizontal direction (length 'Lx'), idem in vertical direction",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.,
    py::arg("nfine")=1
  )

  .def("shape"                 , &xGooseFEM::Mesh::Quad4::FineLayer::shape                 )
  .def("coor"                  , &xGooseFEM::Mesh::Quad4::FineLayer::coor                  )
  .def("conn"                  , &xGooseFEM::Mesh::Quad4::FineLayer::conn                  )
  .def("nelem"                 , &xGooseFEM::Mesh::Quad4::FineLayer::nelem                 )
  .def("nnode"                 , &xGooseFEM::Mesh::Quad4::FineLayer::nnode                 )
  .def("nne"                   , &xGooseFEM::Mesh::Quad4::FineLayer::nne                   )
  .def("ndim"                  , &xGooseFEM::Mesh::Quad4::FineLayer::ndim                  )
  .def("elementsMiddleLayer"   , &xGooseFEM::Mesh::Quad4::FineLayer::elementsMiddleLayer   )
  .def("nodesBottomEdge"       , &xGooseFEM::Mesh::Quad4::FineLayer::nodesBottomEdge       )
  .def("nodesTopEdge"          , &xGooseFEM::Mesh::Quad4::FineLayer::nodesTopEdge          )
  .def("nodesLeftEdge"         , &xGooseFEM::Mesh::Quad4::FineLayer::nodesLeftEdge         )
  .def("nodesRightEdge"        , &xGooseFEM::Mesh::Quad4::FineLayer::nodesRightEdge        )
  .def("nodesBottomOpenEdge"   , &xGooseFEM::Mesh::Quad4::FineLayer::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &xGooseFEM::Mesh::Quad4::FineLayer::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &xGooseFEM::Mesh::Quad4::FineLayer::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &xGooseFEM::Mesh::Quad4::FineLayer::nodesRightOpenEdge    )
  .def("nodesBottomLeftCorner" , &xGooseFEM::Mesh::Quad4::FineLayer::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &xGooseFEM::Mesh::Quad4::FineLayer::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &xGooseFEM::Mesh::Quad4::FineLayer::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &xGooseFEM::Mesh::Quad4::FineLayer::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &xGooseFEM::Mesh::Quad4::FineLayer::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &xGooseFEM::Mesh::Quad4::FineLayer::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &xGooseFEM::Mesh::Quad4::FineLayer::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &xGooseFEM::Mesh::Quad4::FineLayer::nodesRightTopCorner   )
  .def("nodesPeriodic"         , &xGooseFEM::Mesh::Quad4::FineLayer::nodesPeriodic         )
  .def("nodesOrigin"           , &xGooseFEM::Mesh::Quad4::FineLayer::nodesOrigin           )
  .def("dofs"                  , &xGooseFEM::Mesh::Quad4::FineLayer::dofs                  )
  .def("dofsPeriodic"          , &xGooseFEM::Mesh::Quad4::FineLayer::dofsPeriodic          )

  .def("__repr__",
    [](const xGooseFEM::Mesh::Quad4::FineLayer &){ return "<GooseFEM.Mesh.Quad4.FineLayer>"; }
  );

// -------------------------------------------------------------------------------------------------

}

// ====================================== GooseFEM.Mesh.Tri3 =======================================

{

py::module sm = mMesh.def_submodule("Tri3" , "Linear triangular elements (2D)");

// -------------------------------------------------------------------------------------------------

py::class_<xGooseFEM::Mesh::Tri3::Regular>(sm, "Regular")

  .def(
    py::init<size_t,size_t,double>(),
    "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.
  )

  .def("coor"                  , &xGooseFEM::Mesh::Tri3::Regular::coor                  )
  .def("conn"                  , &xGooseFEM::Mesh::Tri3::Regular::conn                  )
  .def("nelem"                 , &xGooseFEM::Mesh::Tri3::Regular::nelem                 )
  .def("nnode"                 , &xGooseFEM::Mesh::Tri3::Regular::nnode                 )
  .def("nne"                   , &xGooseFEM::Mesh::Tri3::Regular::nne                   )
  .def("ndim"                  , &xGooseFEM::Mesh::Tri3::Regular::ndim                  )

  .def("nodesBottomEdge"       , &xGooseFEM::Mesh::Tri3::Regular::nodesBottomEdge       )
  .def("nodesTopEdge"          , &xGooseFEM::Mesh::Tri3::Regular::nodesTopEdge          )
  .def("nodesLeftEdge"         , &xGooseFEM::Mesh::Tri3::Regular::nodesLeftEdge         )
  .def("nodesRightEdge"        , &xGooseFEM::Mesh::Tri3::Regular::nodesRightEdge        )
  .def("nodesBottomOpenEdge"   , &xGooseFEM::Mesh::Tri3::Regular::nodesBottomOpenEdge   )
  .def("nodesTopOpenEdge"      , &xGooseFEM::Mesh::Tri3::Regular::nodesTopOpenEdge      )
  .def("nodesLeftOpenEdge"     , &xGooseFEM::Mesh::Tri3::Regular::nodesLeftOpenEdge     )
  .def("nodesRightOpenEdge"    , &xGooseFEM::Mesh::Tri3::Regular::nodesRightOpenEdge    )

  .def("nodesBottomLeftCorner" , &xGooseFEM::Mesh::Tri3::Regular::nodesBottomLeftCorner )
  .def("nodesBottomRightCorner", &xGooseFEM::Mesh::Tri3::Regular::nodesBottomRightCorner)
  .def("nodesTopLeftCorner"    , &xGooseFEM::Mesh::Tri3::Regular::nodesTopLeftCorner    )
  .def("nodesTopRightCorner"   , &xGooseFEM::Mesh::Tri3::Regular::nodesTopRightCorner   )
  .def("nodesLeftBottomCorner" , &xGooseFEM::Mesh::Tri3::Regular::nodesLeftBottomCorner )
  .def("nodesLeftTopCorner"    , &xGooseFEM::Mesh::Tri3::Regular::nodesLeftTopCorner    )
  .def("nodesRightBottomCorner", &xGooseFEM::Mesh::Tri3::Regular::nodesRightBottomCorner)
  .def("nodesRightTopCorner"   , &xGooseFEM::Mesh::Tri3::Regular::nodesRightTopCorner   )

  .def("nodesPeriodic"         , &xGooseFEM::Mesh::Tri3::Regular::nodesPeriodic         )
  .def("nodesOrigin"           , &xGooseFEM::Mesh::Tri3::Regular::nodesOrigin           )

  .def("dofs"                  , &xGooseFEM::Mesh::Tri3::Regular::dofs                  )
  .def("dofsPeriodic"          , &xGooseFEM::Mesh::Tri3::Regular::dofsPeriodic          )

  .def("__repr__", [](const xGooseFEM::Mesh::Tri3::Regular &){ return "<GooseFEM.Mesh.Tri3.Regular>"; });

// -------------------------------------------------------------------------------------------------

sm.def("getOrientation", &xGooseFEM::Mesh::Tri3::getOrientation, "Get the orientation of each element", py::arg("coor"), py::arg("conn"));

sm.def("retriangulate", &xGooseFEM::Mesh::Tri3::retriangulate, "Re-triangulate existing mesh", py::arg("coor"), py::arg("conn"), py::arg("orientation")=-1);

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

}


/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include <Eigen/Eigen>

#define GOOSEFEM_ENABLE_ASSERT
#define GOOSEFEM_ENABLE_WARNING_PYTHON

namespace py = pybind11;

#include "version.hpp"
#include "Allocate.hpp"
#include "Vector.hpp"
#include "VectorPartitioned.hpp"
#include "VectorPartitionedTyings.hpp"
#include "Matrix.hpp"
#include "MatrixPartitioned.hpp"
#include "MatrixPartitionedTyings.hpp"
#include "MatrixDiagonal.hpp"
#include "MatrixDiagonalPartitioned.hpp"
#include "Element.hpp"
#include "ElementQuad4.hpp"
#include "ElementQuad4Planar.hpp"
#include "ElementQuad4Axisymmetric.hpp"
#include "ElementHex8.hpp"
#include "Mesh.hpp"
#include "MeshTri3.hpp"
#include "MeshQuad4.hpp"
#include "MeshHex8.hpp"
#include "Iterate.hpp"

PYBIND11_MODULE(GooseFEM, m) {

xt::import_numpy();

// --------
// GooseFEM
// --------

m.doc() = "Some simple finite element meshes and operations";

init_version(m);
init_Allocate(m);
init_Vector(m);
init_VectorPartitioned(m);
init_VectorPartitionedTyings(m);
init_Matrix(m);
init_MatrixPartitioned(m);
init_MatrixPartitionedTyings(m);
init_MatrixDiagonal(m);
init_MatrixDiagonalPartitioned(m);

// ----------------
// GooseFEM.Iterate
// ----------------

py::module mIterate = m.def_submodule("Iterate", "Iteration support tools");

init_Iterate(mIterate);


// ----------------
// GooseFEM.Element
// ----------------

py::module mElement = m.def_submodule("Element", "Generic element routines");

init_Element(mElement);

// ----------------------
// GooseFEM.Element.Quad4
// ----------------------

py::module mElementQuad4 = mElement.def_submodule("Quad4", "Linear quadrilateral elements (2D)");
py::module mElementQuad4Gauss = mElementQuad4.def_submodule("Gauss", "Gauss quadrature");
py::module mElementQuad4Nodal = mElementQuad4.def_submodule("Nodal", "Nodal quadrature");
py::module mElementQuad4MidPoint = mElementQuad4.def_submodule("MidPoint", "MidPoint quadrature");

init_ElementQuad4(mElementQuad4);
init_ElementQuad4Planar(mElementQuad4);
init_ElementQuad4Axisymmetric(mElementQuad4);
init_ElementQuad4Gauss(mElementQuad4Gauss);
init_ElementQuad4Nodal(mElementQuad4Nodal);
init_ElementQuad4MidPoint(mElementQuad4MidPoint);

// ---------------------
// GooseFEM.Element.Hex8
// ---------------------

py::module mElementHex8 = mElement.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");
py::module mElementHex8Gauss = mElementHex8.def_submodule("Gauss", "Gauss quadrature");
py::module mElementHex8Nodal = mElementHex8.def_submodule("Nodal", "Nodal quadrature");

init_ElementHex8(mElementHex8);
init_ElementHex8Gauss(mElementHex8Gauss);
init_ElementHex8Nodal(mElementHex8Nodal);

// -------------
// GooseFEM.Mesh
// -------------

py::module mMesh = m.def_submodule("Mesh", "Generic mesh routines");

init_Mesh(mMesh);

// ------------------
// GooseFEM.Mesh.Tri3
// ------------------

py::module mMeshTri3 = mMesh.def_submodule("Tri3", "Linear triangular elements (2D)");

init_MeshTri3(mMeshTri3);

// -------------------
// GooseFEM.Mesh.Quad4
// -------------------

py::module mMeshQuad4 = mMesh.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

init_MeshQuad4(mMeshQuad4);

py::module mMeshQuad4Map = mMeshQuad4.def_submodule("Map", "Map mesh objects");

init_MeshQuad4Map(mMeshQuad4Map);

// ------------------
// GooseFEM.Mesh.Hex8
// ------------------

py::module mMeshHex8 = mMesh.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

init_MeshHex8(mMeshHex8);

}


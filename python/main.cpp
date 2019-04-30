/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <pyxtensor/pyxtensor.hpp>

#define GOOSEFEM_ENABLE_ASSERT
#include "../include/GooseFEM/GooseFEM.h"

// =================================================================================================

namespace py = pybind11;

// =================================================================================================

#include "Vector.hpp"
#include "VectorPartitioned.hpp"
#include "MatrixDiagonalPartitioned.hpp"
#include "Element.hpp"
#include "ElementQuad4.hpp"
#include "ElementHex8.hpp"
#include "Mesh.hpp"
#include "MeshTri3.hpp"
#include "MeshQuad4.hpp"
#include "MeshHex8.hpp"
#include "ParaView.hpp"

// =================================================================================================

PYBIND11_MODULE(GooseFEM, m) {

// -------------------------------------------------------------------------------------------------

m.doc() = "Some simple finite element meshes and operations";

// -------------------------------------------------------------------------------------------------

init_Vector(m);

// -------------------------------------------------------------------------------------------------

init_VectorPartitioned(m);

// -------------------------------------------------------------------------------------------------

init_MatrixDiagonalPartitioned(m);

// -------------------------------------------------------------------------------------------------

py::module mElement = m.def_submodule("Element", "Generic element routines");

init_Element(mElement);

// -------------------------------------------------------------------------------------------------

py::module mElementQuad4 = mElement.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

py::module mElementQuad4Gauss = mElementQuad4.def_submodule("Gauss", "Gauss quadrature");

py::module mElementQuad4Nodal = mElementQuad4.def_submodule("Nodal", "Nodal quadrature");

init_ElementQuad4(mElementQuad4);

init_ElementQuad4Gauss(mElementQuad4Gauss);

init_ElementQuad4Nodal(mElementQuad4Nodal);

// -------------------------------------------------------------------------------------------------

py::module mElementHex8 = mElement.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

py::module mElementHex8Gauss = mElementHex8.def_submodule("Gauss", "Gauss quadrature");

py::module mElementHex8Nodal = mElementHex8.def_submodule("Nodal", "Nodal quadrature");

init_ElementHex8(mElementHex8);

init_ElementHex8Gauss(mElementHex8Gauss);

init_ElementHex8Nodal(mElementHex8Nodal);

// -------------------------------------------------------------------------------------------------

py::module mMesh = m.def_submodule("Mesh", "Generic mesh routines");

init_Mesh(mMesh);

// -------------------------------------------------------------------------------------------------

py::module mMeshTri3 = mMesh.def_submodule("Tri3", "Linear triangular elements (2D)");

init_MeshTri3(mMeshTri3);

// -------------------------------------------------------------------------------------------------

py::module mMeshQuad4 = mMesh.def_submodule("Quad4", "Linear quadrilateral elements (2D)");

init_MeshQuad4(mMeshQuad4);

py::module mMeshQuad4Map = mMeshQuad4.def_submodule("Map", "Map mesh objects");

init_MeshQuad4Map(mMeshQuad4Map);

// -------------------------------------------------------------------------------------------------

py::module mMeshHex8 = mMesh.def_submodule("Hex8", "Linear hexahedron (brick) elements (3D)");

init_MeshHex8(mMeshHex8);

// -------------------------------------------------------------------------------------------------

py::module mParaView = m.def_submodule("ParaView", "ParaView output files");

// -------------------------------------------------------------------------------------------------

py::module mParaViewHDF5 = mParaView.def_submodule("HDF5", "ParaView/HDF5 support using XDMF files");

init_ParaViewHDF5(mParaViewHDF5);


// =================================================================================================

}


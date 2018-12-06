/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_H
#define GOOSEFEM_H

// =================================================================================================

#ifdef EIGEN_WORLD_VERSION
#define GOOSEFEM_EIGEN
#endif

// =================================================================================================

#include "Mesh.h"
#include "MeshTri3.h"
#include "MeshQuad4.h"
#include "MeshHex8.h"
#include "Element.h"
#include "ElementQuad4.h"
#include "ElementQuad4Planar.h"
#include "ElementQuad4Axisymmetric.h"
#include "ElementHex8.h"
#include "Vector.h"
#include "VectorPartitioned.h"
#include "MatrixDiagonal.h"
#include "MatrixDiagonalPartitioned.h"
#include "Iterate.h"

#ifdef GOOSEFEM_EIGEN
#include "Matrix.h"
#include "MatrixPartitioned.h"
#endif

// =================================================================================================

#endif

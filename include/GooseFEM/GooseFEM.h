/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_H
#define GOOSEFEM_H

#ifdef EIGEN_WORLD_VERSION
    #define GOOSEFEM_EIGEN
#endif

#include "Element.h"
#include "ElementHex8.h"
#include "ElementQuad4.h"
#include "ElementQuad4Axisymmetric.h"
#include "ElementQuad4Planar.h"
#include "Iterate.h"
#include "MatrixDiagonal.h"
#include "MatrixDiagonalPartitioned.h"
#include "Mesh.h"
#include "MeshHex8.h"
#include "MeshQuad4.h"
#include "MeshTri3.h"
#include "Vector.h"
#include "VectorPartitioned.h"

#ifdef GOOSEFEM_EIGEN
#include "Matrix.h"
#include "MatrixPartitioned.h"
#include "MatrixPartitionedTyings.h"
#include "TyingsPeriodic.h"
#include "VectorPartitionedTyings.h"
#endif

#endif

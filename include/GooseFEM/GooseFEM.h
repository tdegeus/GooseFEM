/**
Basic include of common methods.

\file GooseFEM.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_H
#define GOOSEFEM_H

#ifdef EIGEN_WORLD_VERSION
#define GOOSEFEM_EIGEN
#endif

#include "version.h"

#include "Allocate.h"
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

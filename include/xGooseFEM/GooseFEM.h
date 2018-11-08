/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_H
#define XGOOSEFEM_H

// =================================================================================================

#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <assert.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <limits>
#include <algorithm>
#include <math.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include <xtensor/xarray.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xstrided_view.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xutils.hpp>
#include <xtensor/xlayout.hpp>
#include <xtensor/xinfo.hpp>
#include <xtensor/xnoalias.hpp>

using namespace xt::placeholders;

// =================================================================================================

#define GOOSEFEM_WORLD_VERSION 0
#define GOOSEFEM_MAJOR_VERSION 1
#define GOOSEFEM_MINOR_VERSION 0

#define GOOSEFEM_VERSION_AT_LEAST(x,y,z) \
  (GOOSEFEM_WORLD_VERSION>x || (GOOSEFEM_WORLD_VERSION>=x && \
  (GOOSEFEM_MAJOR_VERSION>y || (GOOSEFEM_MAJOR_VERSION>=y && \
                                GOOSEFEM_MINOR_VERSION>=z))))

#define GOOSEFEM_VERSION(x,y,z) \
  (GOOSEFEM_WORLD_VERSION==x && \
   GOOSEFEM_MAJOR_VERSION==y && \
   GOOSEFEM_MINOR_VERSION==z)

// =================================================================================================

// dummy operation that can be use to suppress the "unused parameter" warnings
#define UNUSED(p) ( (void)(p) )

// =================================================================================================

// alias Eigen sparse matrices
namespace xGooseFEM
{
  typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMatD;
  typedef Eigen::SparseMatrix<size_t,Eigen::RowMajor> SpMatS;
  typedef Eigen::Triplet<double> TripD;
}

// =================================================================================================

#include "Mesh.h"
#include "Mesh.hpp"
#include "MeshTri3.h"
#include "MeshTri3.hpp"
#include "MeshQuad4.h"
#include "MeshQuad4.hpp"
#include "MeshHex8.h"
#include "MeshHex8.hpp"
#include "Element.h"
#include "Element.hpp"
#include "ElementQuad4.h"
#include "ElementQuad4.hpp"
#include "ElementQuad4Planar.h"
#include "ElementQuad4Planar.hpp"
#include "ElementHex8.h"
#include "ElementHex8.hpp"
#include "Vector.h"
#include "Vector.hpp"
#include "VectorPartitioned.h"
#include "VectorPartitioned.hpp"
#include "MatrixPartitioned.h"
#include "MatrixPartitioned.hpp"
#include "MatrixDiagonal.h"
#include "MatrixDiagonal.hpp"
#include "MatrixDiagonalPartitioned.h"
#include "MatrixDiagonalPartitioned.hpp"
#include "Iterate.h"
#include "Iterate.hpp"
#include "Dynamics.h"
#include "Dynamics.hpp"

// =================================================================================================

#endif

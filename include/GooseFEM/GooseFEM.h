/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_H
#define GOOSEFEM_H

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

#include <xtensor/xarray.hpp>
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
namespace GooseFEM
{
  typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMatD;
  typedef Eigen::SparseMatrix<size_t,Eigen::RowMajor> SpMatS;
}

// =================================================================================================

#include "Mesh.h"
#include "MeshTri3.h"
#include "MeshQuad4.h"
#include "MeshHex8.h"
#include "Element.h"
#include "ElementQuad4.h"
#include "ElementHex8.h"
#include "Vector.h"
#include "MatrixDiagonal.h"
#include "Iterate.h"
#include "Dynamics.h"

#include "Mesh.hpp"
#include "MeshTri3.hpp"
#include "MeshQuad4.hpp"
#include "MeshHex8.hpp"
#include "Element.hpp"
#include "ElementQuad4.hpp"
#include "ElementHex8.hpp"
#include "Vector.hpp"
#include "MatrixDiagonal.hpp"
#include "Iterate.hpp"
#include "Dynamics.hpp"

// =================================================================================================

#endif

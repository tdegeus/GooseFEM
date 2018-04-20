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
#include <cppmat/cppmat.h>

// =================================================================================================

#define GOOSEFEM_WORLD_VERSION 0
#define GOOSEFEM_MAJOR_VERSION 0
#define GOOSEFEM_MINOR_VERSION 8

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

// alias types
namespace GooseFEM
{
  // - alias Eigen dense matrices
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
  typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatS;
  typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;
  typedef Eigen::Matrix<size_t, Eigen::Dynamic,              1, Eigen::ColMajor> ColS;
  typedef Eigen::Matrix<int   , Eigen::Dynamic,              1, Eigen::ColMajor> ColI;
  // - alias Eigen sparse matrices
  typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMatD;
  typedef Eigen::SparseMatrix<size_t,Eigen::RowMajor> SpMatS;
  // - alias cppmat matrices
  typedef cppmat::matrix<double> ArrD;
  typedef cppmat::matrix<size_t> ArrS;
}

// =================================================================================================

#include "Mesh.h"
#include "MeshTri3.h"
#include "MeshQuad4.h"
#include "MeshHex8.h"
#include "Element.h"
#include "ElementQuad4.h"
#include "Vector.h"
#include "MatrixDiagonal.h"
#include "Iterate.h"

#include "Mesh.cpp"
#include "MeshTri3.cpp"
#include "MeshQuad4.cpp"
#include "MeshHex8.cpp"
#include "Element.cpp"
#include "ElementQuad4.cpp"
#include "Vector.cpp"
#include "MatrixDiagonal.cpp"
#include "Iterate.cpp"


// #include "DynamicsDiagonal.h"
// #include "OverdampedDynamicsDiagonal.h"
// #include "DynamicsDiagonal.cpp"
// #include "OverdampedDynamicsDiagonal.cpp"

// =================================================================================================

#endif

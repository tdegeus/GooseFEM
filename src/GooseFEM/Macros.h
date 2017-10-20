/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MACROS_H
#define GOOSEFEM_MACROS_H

#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <assert.h>
#include <cstdlib>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <Eigen/Dense>

#include <cppmat/tensor.h>
#include <cppmat/tensor2.h>
#include <cppmat/tensor3.h>

// =================================================================================================

#define GOOSEFEM_WORLD_VERSION 0
#define GOOSEFEM_MAJOR_VERSION 0
#define GOOSEFEM_MINOR_VERSION 2

#define GOOSEFEM_VERSION_AT_LEAST(x,y,z) \
  (GOOSEFEM_WORLD_VERSION>x || (GOOSEFEM_WORLD_VERSION>=x && \
  (GOOSEFEM_MAJOR_VERSION>y || (GOOSEFEM_MAJOR_VERSION>=y && \
                                GOOSEFEM_MINOR_VERSION>=z))))

#define GOOSEFEM_VERSION(x,y,z) \
  (GOOSEFEM_WORLD_VERSION==x && \
   GOOSEFEM_MAJOR_VERSION==y && \
   GOOSEFEM_MINOR_VERSION==z)

// =================================================================================================

namespace GooseFEM {

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatS;
typedef Eigen::Matrix<double, Eigen::Dynamic,              1, Eigen::ColMajor> ColD;
typedef Eigen::Matrix<size_t, Eigen::Dynamic,              1, Eigen::ColMajor> ColS;
typedef Eigen::Matrix<int   , Eigen::Dynamic,              1, Eigen::ColMajor> ColI;

// -------------------------------------------------------------------------------------------------

}

#endif

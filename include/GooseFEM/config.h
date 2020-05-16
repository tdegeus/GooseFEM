/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_CONFIG_H
#define GOOSEFEM_CONFIG_H

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

#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xinfo.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xlayout.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xstrided_view.hpp>
#include <xtensor/xutils.hpp>
#include <xtensor/xview.hpp>

using namespace xt::placeholders;

#define UNUSED(p) ((void)(p))

#ifdef GOOSEFEM_ENABLE_ASSERT
    #define GOOSEFEM_ASSERT(expr) GOOSEFEM_ASSERT_IMPL(expr, __FILE__, __LINE__)
    #define GOOSEFEM_ASSERT_IMPL(expr, file, line) \
        if (!(expr)) { \
            throw std::runtime_error( \
                std::string(file) + ':' + std::to_string(line) + \
                ": assertion failed (" #expr ") \n\t"); \
        }
#else
    #define GOOSEFEM_ASSERT(expr)
#endif

#define GOOSEFEM_CHECK(expr) GOOSEFEM_CHECK_IMPL(expr, __FILE__, __LINE__)
#define GOOSEFEM_CHECK_IMPL(expr, file, line) \
    if (!(expr)) { \
        throw std::runtime_error( \
            std::string(file) + ':' + std::to_string(line) + \
            ": assertion failed (" #expr ") \n\t"); \
    }

#define GOOSEFEM_VERSION_MAJOR 0
#define GOOSEFEM_VERSION_MINOR 3
#define GOOSEFEM_VERSION_PATCH 1

#define GOOSEFEM_VERSION_AT_LEAST(x, y, z) \
    (GOOSEFEM_VERSION_MAJOR > x || (GOOSEFEM_VERSION_MAJOR >= x && \
    (GOOSEFEM_VERSION_MINOR > y || (GOOSEFEM_VERSION_MINOR >= y && \
                                    GOOSEFEM_VERSION_PATCH >= z))))

#define GOOSEFEM_VERSION(x, y, z) \
    (GOOSEFEM_VERSION_MAJOR == x && \
     GOOSEFEM_VERSION_MINOR == y && \
     GOOSEFEM_VERSION_PATCH == z)

#endif

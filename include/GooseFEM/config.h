/**
Basic configuration:

-   Include general dependencies.
-   Define assertions.
-   Define version.
-   Define git commit hash/branch.

\file config.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_CONFIG_H
#define GOOSEFEM_CONFIG_H

#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <xtensor/xadapt.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xinfo.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xlayout.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xshape.hpp>
#include <xtensor/xstrided_view.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xutils.hpp>
#include <xtensor/xview.hpp>

using namespace xt::placeholders;

#define Q(x) #x
#define QUOTE(x) Q(x)

#define UNUSED(p) ((void)(p))

/**
All assertions are implementation as::

    GOOSEFEM_ASSERT(...)

They can be enabled by::

    #define GOOSEFEM_ENABLE_ASSERT

(before including GooseFEM).
The advantage is that:

-   File and line-number are displayed if the assertion fails.
-   GooseFEM's assertions can be enabled/disabled independently from those of other libraries.

\throw std::runtime_error
*/
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

/**
Assertion that cannot be switched of. Implement assertion by::

    GOOSEFEM_CHECK(...)

\throw std::runtime_error
*/
#define GOOSEFEM_CHECK(expr) GOOSEFEM_CHECK_IMPL(expr, __FILE__, __LINE__)
#define GOOSEFEM_CHECK_IMPL(expr, file, line) \
    if (!(expr)) { \
        throw std::runtime_error( \
            std::string(file) + ':' + std::to_string(line) + \
            ": assertion failed (" #expr ") \n\t"); \
    }

#define GOOSEFEM_WIP_ASSERT(expr) GOOSEFEM_CHECK_IMPL(expr, __FILE__, __LINE__)
#define GOOSEFEM_WIP_ASSERT_IMPL(expr, file, line) \
    if (!(expr)) { \
        throw std::runtime_error( \
            std::string(file) + ':' + std::to_string(line) + \
            ": WIP, please extend the code, assertion failed (" #expr ") \n\t"); \
    }

#endif

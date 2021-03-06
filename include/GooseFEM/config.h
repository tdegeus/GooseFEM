/**
Basic configuration:

-   Include general dependencies.
-   Define assertions.

\file config.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_CONFIG_H
#define GOOSEFEM_CONFIG_H

/**
\cond
*/
#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"
/**
\endcond
*/

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
#include <array>
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

/**
\cond
*/
using namespace xt::placeholders;

#define Q(x) #x
#define QUOTE(x) Q(x)

#define UNUSED(p) ((void)(p))

#define GOOSEFEM_WARNING_IMPL(message, file, line) \
    std::cout << \
        std::string(file) + ':' + std::to_string(line) + \
        ": " message ") \n\t"; \

#define GOOSEFEM_ASSERT_IMPL(expr, file, line) \
    if (!(expr)) { \
        throw std::runtime_error( \
            std::string(file) + ':' + std::to_string(line) + \
            ": assertion failed (" #expr ") \n\t"); \
    }

/**
\endcond
*/

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
#else
#define GOOSEFEM_ASSERT(expr)
#endif

/**
Assertion that cannot be switched off. Implement assertion by::

    GOOSEFEM_CHECK(...)

\throw std::runtime_error
*/
#define GOOSEFEM_CHECK(expr) GOOSEFEM_ASSERT_IMPL(expr, __FILE__, __LINE__)

/**
Assertion that concerns temporary implementation limitations.
Implement assertion by::

    GOOSEFEM_WIP_ASSERT(...)

\throw std::runtime_error
*/
#define GOOSEFEM_WIP_ASSERT(expr) GOOSEFEM_ASSERT_IMPL(expr, __FILE__, __LINE__)

/**
All warnings are implemented as::

    GOOSEFEM_WARNING(...)

They can be disabled by::

    #define GOOSEFEM_DISABLE_WARNING
*/
#ifdef GOOSEFEM_DISABLE_WARNING
#define GOOSEFEM_WARNING(message)
#else
#define GOOSEFEM_WARNING(message) GOOSEFEM_WARNING_IMPL(message, __FILE__, __LINE__)
#endif

/**
All warnings specific to the Python API are implemented as::

    GOOSEFEM_WARNING_PYTHON(...)

They can be enabled by::

    #define GOOSEFEM_ENABLE_WARNING_PYTHON
*/
#ifdef GOOSEFEM_ENABLE_WARNING_PYTHON
#define GOOSEFEM_WARNING_PYTHON(message) GOOSEFEM_WARNING_IMPL(message, __FILE__, __LINE__)
#else
#define GOOSEFEM_WARNING_PYTHON(message)
#endif

/**
Toolbox to perform finite element computations.
*/
namespace GooseFEM {}

#endif

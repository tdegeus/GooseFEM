/**
 *  Python API of `assertions.h`/
 *
 *  \file assertions.hpp
 *  \copyright Copyright 2017. Tom de Geus. All rights reserved.
 *  \license This project is released under the GNU Public License (GPLv3).
 */

#include <GooseFEM/assertions.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>

namespace py = pybind11;

void init_assertions(py::module& mod)
{
    mod.def(
        "is_unique",
        &GooseFEM::is_unique<xt::pyarray<long>>,
        "See :cpp:func:`GooseFEM::is_unique`.",
        py::arg("a"));
}

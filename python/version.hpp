/**
 * @file
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#include <GooseFEM/version.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_version(py::module& m)
{

    m.def("version", &GooseFEM::version, "See :cpp:func:`GooseFEM::version`.");

    m.def(
        "version_dependencies",
        &GooseFEM::version_dependencies,
        "See :cpp:func:`GooseFEM::version_dependencies`.");
}

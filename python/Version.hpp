/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_version(py::module& m)
{

    m.def("version",
          &GooseFEM::version,
          "Return version string."
          "See :cpp:class:`GooseFEM::version`.");

    m.def("version_dependencies",
          &GooseFEM::version_dependencies,
          "Return version information of library and its dependencies."
          "See :cpp:class:`GooseFEM::version_dependencies`.");

}

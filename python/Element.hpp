/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_Element(py::module& m)
{

    m.def(
        "asElementVector",
        &GooseFEM::Element::asElementVector,
        "Convert nodal vector [nnode, ndim] to nodal vector stored per element [nelem, nne, ndim]",
        py::arg("conn"),
        py::arg("nodevec"));

    m.def(
        "assembleElementVector",
        &GooseFEM::Element::assembleNodeVector,
        "Assemble nodal vector stored per element [nelem, nne, ndim] to nodal vector [nnode, ndim]",
        py::arg("conn"),
        py::arg("elemvec"));
}

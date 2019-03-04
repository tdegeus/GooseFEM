/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <pyxtensor/pyxtensor.hpp>

#include "../include/GooseFEM/GooseFEM.h"

// =================================================================================================

namespace py = pybind11;

// =================================================================================================

void init_Element(py::module &m)
{

m.def("asElementVector"      , &GooseFEM::Element::asElementVector   , "Covert to 'elemvec'", py::arg("conn"), py::arg("nodevec"));
m.def("assembleElementVector", &GooseFEM::Element::assembleNodeVector, "Assemble 'nodevec'" , py::arg("conn"), py::arg("elemvec"));

}

// =================================================================================================


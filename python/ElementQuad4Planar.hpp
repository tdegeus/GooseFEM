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

void init_ElementQuad4Planar(py::module &m)
{

py::class_<GooseFEM::Element::Quad4::QuadraturePlanar>(m, "QuadraturePlanar")

  .def(py::init<const xt::xtensor<double,3>&>(),
    "QuadraturePlanar",
    py::arg("x"))

  .def(py::init<const xt::xtensor<double,3>&,
                const xt::xtensor<double,2>&,
                const xt::xtensor<double,1>&>(),
    "QuadraturePlanar",
    py::arg("x"),
    py::arg("xi"),
    py::arg("w"))

  .def("update_x",
    &GooseFEM::Element::Quad4::QuadraturePlanar::update_x,
    "Update the nodal positions")

  .def("nelem",
    &GooseFEM::Element::Quad4::QuadraturePlanar::nelem,
    "Number of elements")

  .def("nne",
    &GooseFEM::Element::Quad4::QuadraturePlanar::nne,
    "Number of nodes per element")

  .def("ndim",
    &GooseFEM::Element::Quad4::QuadraturePlanar::ndim,
    "Number of dimensions")

  .def("nip",
    &GooseFEM::Element::Quad4::QuadraturePlanar::nip,
    "Number of integration points")

  .def("DV",
    py::overload_cast<size_t>(&GooseFEM::Element::Quad4::QuadraturePlanar::DV, py::const_),
    "Integration point volume (qtensor)")

  .def("DV",
    py::overload_cast<>(&GooseFEM::Element::Quad4::QuadraturePlanar::DV, py::const_),
    "Integration point volume (qscalar)")

  .def("GradN_vector",
    py::overload_cast<const xt::xtensor<double,3>&>(
      &GooseFEM::Element::Quad4::QuadraturePlanar::GradN_vector, py::const_),
    "Dyadic product, returns 'qtensor'",
    py::arg("elemvec"))

  .def("GradN_vector_T",
    py::overload_cast<const xt::xtensor<double,3>&>(
      &GooseFEM::Element::Quad4::QuadraturePlanar::GradN_vector_T, py::const_),
    "Dyadic product, returns 'qtensor'",
    py::arg("elemvec"))

  .def("SymGradN_vector",
    py::overload_cast<const xt::xtensor<double,3>&>(
      &GooseFEM::Element::Quad4::QuadraturePlanar::SymGradN_vector, py::const_),
    "Dyadic product, returns 'qtensor'",
    py::arg("elemvec"))

  .def("Int_N_scalar_NT_dV",
    py::overload_cast<const xt::xtensor<double,2>&>(
      &GooseFEM::Element::Quad4::QuadraturePlanar::Int_N_scalar_NT_dV, py::const_),
    "Integration, returns 'elemmat'",
    py::arg("qscalar"))

  .def("Int_gradN_dot_tensor2_dV",
    py::overload_cast<const xt::xtensor<double,4>&>(
      &GooseFEM::Element::Quad4::QuadraturePlanar::Int_gradN_dot_tensor2_dV, py::const_),
    "Integration, returns 'elemvec'",
    py::arg("qtensor"))

  .def("Int_gradN_dot_tensor4_dot_gradNT_dV",
    py::overload_cast<const xt::xtensor<double,6>&>(
      &GooseFEM::Element::Quad4::QuadraturePlanar::Int_gradN_dot_tensor4_dot_gradNT_dV, py::const_),
    "Integration, returns 'elemvec'",
    py::arg("qtensor"))

  .def("__repr__", [](const GooseFEM::Element::Quad4::QuadraturePlanar&){
    return "<GooseFEM.Element.Quad4.QuadraturePlanar>";});

}

// =================================================================================================


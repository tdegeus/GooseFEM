/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_ElementQuad4(py::module& m)
{

    py::class_<GooseFEM::Element::Quad4::Quadrature>(m, "Quadrature")

        .def(py::init<const xt::xtensor<double, 3>&>(), "Quadrature", py::arg("x"))

        .def(
            py::init<
                const xt::xtensor<double, 3>&,
                const xt::xtensor<double, 2>&,
                const xt::xtensor<double, 1>&>(),
            "Quadrature",
            py::arg("x"),
            py::arg("xi"),
            py::arg("w"))

        .def(
            "update_x",
            &GooseFEM::Element::Quad4::Quadrature::update_x,
            "Update the nodal positions")

        .def("nelem", &GooseFEM::Element::Quad4::Quadrature::nelem, "Number of elements")

        .def("nne", &GooseFEM::Element::Quad4::Quadrature::nne, "Number of nodes per element")

        .def("ndim", &GooseFEM::Element::Quad4::Quadrature::ndim, "Number of dimensions")

        .def("nip", &GooseFEM::Element::Quad4::Quadrature::nip, "Number of integration points")

        .def("dV", &GooseFEM::Element::Quad4::Quadrature::dV, "Integration point volume (qscalar)")

        .def(
            "Interp_N_vector",
            &GooseFEM::Element::Quad4::Quadrature::Interp_N_vector<double>,
            "See :cpp:class:`GooseFEM::Quad4::Quadrature::Interp_N_vector`.",
            py::arg("elemvec"))

        .def(
            "GradN_vector",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Element::Quad4::Quadrature::GradN_vector, py::const_),
            "Dyadic product, returns 'qtensor'",
            py::arg("elemvec"))

        .def(
            "GradN_vector_T",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Element::Quad4::Quadrature::GradN_vector_T, py::const_),
            "Dyadic product, returns 'qtensor'",
            py::arg("elemvec"))

        .def(
            "SymGradN_vector",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Element::Quad4::Quadrature::SymGradN_vector, py::const_),
            "Dyadic product, returns 'qtensor'",
            py::arg("elemvec"))

        .def(
            "Int_N_scalar_NT_dV",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::Element::Quad4::Quadrature::Int_N_scalar_NT_dV, py::const_),
            "Integration, returns 'elemmat'",
            py::arg("qscalar"))

        .def(
            "Int_gradN_dot_tensor2_dV",
            py::overload_cast<const xt::xtensor<double, 4>&>(
                &GooseFEM::Element::Quad4::Quadrature::Int_gradN_dot_tensor2_dV, py::const_),
            "Integration, returns 'elemvec'",
            py::arg("qtensor"))

        .def(
            "Int_gradN_dot_tensor4_dot_gradNT_dV",
            py::overload_cast<const xt::xtensor<double, 6>&>(
                &GooseFEM::Element::Quad4::Quadrature::Int_gradN_dot_tensor4_dot_gradNT_dV,
                py::const_),
            "Integration, returns 'elemvec'",
            py::arg("qtensor"))

        .def(
            "AsTensor",
            (xt::xarray<double>(GooseFEM::Element::Quad4::Quadrature::*)(
                size_t, const xt::xtensor<double, 2>&) const)
                &GooseFEM::Element::Quad4::Quadrature::AsTensor<double>,
            "Convert 'qscalar' to 'qtensor' of certain rank")

        .def(
            "AllocateQtensor",
            (xt::xarray<double>(GooseFEM::Element::Quad4::Quadrature::*)(
                size_t) const)
                &GooseFEM::Element::Quad4::Quadrature::AllocateQtensor<double>,
            "Allocate 'qtensor'",
            py::arg("rank"))

        .def(
            "AllocateQtensor",
            (xt::xarray<double>(GooseFEM::Element::Quad4::Quadrature::*)(
                size_t, double) const)
                &GooseFEM::Element::Quad4::Quadrature::AllocateQtensor<double>,
            "Allocate 'qtensor'",
            py::arg("rank"),
            py::arg("val"))

        .def(
            "AllocateQscalar",
            py::overload_cast<>(
                &GooseFEM::Element::Quad4::Quadrature::AllocateQscalar<double>, py::const_),
            "Allocate 'qscalar'")

        .def(
            "AllocateQscalar",
            py::overload_cast<double>(
                &GooseFEM::Element::Quad4::Quadrature::AllocateQscalar<double>, py::const_),
            "Allocate 'qscalar'",
            py::arg("val"))

        .def("__repr__", [](const GooseFEM::Element::Quad4::Quadrature&) {
            return "<GooseFEM.Element.Quad4.Quadrature>";
        });
}

void init_ElementQuad4Gauss(py::module& m)
{

    m.def("nip", &GooseFEM::Element::Quad4::Gauss::nip, "Return number of integration point");

    m.def("xi", &GooseFEM::Element::Quad4::Gauss::xi, "Return integration point coordinates");

    m.def("w", &GooseFEM::Element::Quad4::Gauss::w, "Return integration point weights");
}

void init_ElementQuad4Nodal(py::module& m)
{

    m.def("nip", &GooseFEM::Element::Quad4::Nodal::nip, "Return number of integration point");

    m.def("xi", &GooseFEM::Element::Quad4::Nodal::xi, "Return integration point coordinates");

    m.def("w", &GooseFEM::Element::Quad4::Nodal::w, "Return integration point weights");
}

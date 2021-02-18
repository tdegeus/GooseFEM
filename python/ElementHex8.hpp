/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_ElementHex8(py::module& m)
{

    py::class_<GooseFEM::Element::Hex8::Quadrature>(m, "Quadrature")

        .def(py::init<const xt::xtensor<double, 3>&>(), "Quadrature", py::arg("x"))

        .def(py::init<
                const xt::xtensor<double, 3>&,
                const xt::xtensor<double, 2>&,
                const xt::xtensor<double, 1>&>(),
            "Quadrature",
            py::arg("x"),
            py::arg("xi"),
            py::arg("w"))

        .def("update_x",
            &GooseFEM::Element::Hex8::Quadrature::update_x,
            "Update the nodal positions")

        .def("dV", &GooseFEM::Element::Hex8::Quadrature::dV, "Integration point volume (qscalar)")

        .def("GradN_vector",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Element::Hex8::Quadrature::GradN_vector, py::const_),
            "Dyadic product, returns 'qtensor'",
            py::arg("elemvec"))

        .def("GradN_vector_T",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Element::Hex8::Quadrature::GradN_vector_T, py::const_),
            "Dyadic product, returns 'qtensor'",
            py::arg("elemvec"))

        .def("SymGradN_vector",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Element::Hex8::Quadrature::SymGradN_vector, py::const_),
            "Dyadic product, returns 'qtensor'",
            py::arg("elemvec"))

        .def("Int_N_scalar_NT_dV",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::Element::Hex8::Quadrature::Int_N_scalar_NT_dV, py::const_),
            "Integration, returns 'elemmat'",
            py::arg("qscalar"))

        .def("Int_gradN_dot_tensor2_dV",
            py::overload_cast<const xt::xtensor<double, 4>&>(
                &GooseFEM::Element::Hex8::Quadrature::Int_gradN_dot_tensor2_dV, py::const_),
            "Integration, returns 'elemvec'",
            py::arg("qtensor"))

        .def("Int_gradN_dot_tensor4_dot_gradNT_dV",
            py::overload_cast<const xt::xtensor<double, 6>&>(
                &GooseFEM::Element::Hex8::Quadrature::Int_gradN_dot_tensor4_dot_gradNT_dV,
                py::const_),
            "Integration, returns 'elemvec'",
            py::arg("qtensor"))

        // Derived from QuadratureBase

        .def("nelem", &GooseFEM::Element::Hex8::Quadrature::nelem, "Number of elements")

        .def("nne", &GooseFEM::Element::Hex8::Quadrature::nne, "Number of nodes per element")

        .def("ndim", &GooseFEM::Element::Hex8::Quadrature::ndim, "Number of dimensions")

        .def("nip", &GooseFEM::Element::Hex8::Quadrature::nip, "Number of integration points")

        .def("AsTensor",
            (xt::xarray<double>(GooseFEM::Element::Hex8::Quadrature::*)(
                size_t, const xt::xtensor<double, 2>&) const)
                &GooseFEM::Element::Hex8::Quadrature::AsTensor<double>,
            "Convert 'qscalar' to 'qtensor' of certain rank")

        .def("AllocateQtensor",
            (xt::xarray<double>(GooseFEM::Element::Hex8::Quadrature::*)(size_t) const)
                &GooseFEM::Element::Hex8::Quadrature::AllocateQtensor<double>,
            "Allocate 'qtensor'",
            py::arg("rank"))

        .def("AllocateQtensor",
            (xt::xarray<double>(GooseFEM::Element::Hex8::Quadrature::*)(size_t, double) const)
                &GooseFEM::Element::Hex8::Quadrature::AllocateQtensor<double>,
            "Allocate 'qtensor'",
            py::arg("rank"),
            py::arg("val"))

        .def("AllocateQscalar",
            py::overload_cast<>(
                &GooseFEM::Element::Hex8::Quadrature::AllocateQscalar<double>, py::const_),
            "Allocate 'qscalar'")

        .def("AllocateQscalar",
            py::overload_cast<double>(
                &GooseFEM::Element::Hex8::Quadrature::AllocateQscalar<double>, py::const_),
            "Allocate 'qscalar'",
            py::arg("val"))

        .def("__repr__", [](const GooseFEM::Element::Hex8::Quadrature&) {
            return "<GooseFEM.Element.Hex8.Quadrature>";
        });
}

void init_ElementHex8Gauss(py::module& m)
{

    m.def("nip", &GooseFEM::Element::Hex8::Gauss::nip, "Return number of integration point");

    m.def("xi", &GooseFEM::Element::Hex8::Gauss::xi, "Return integration point coordinates");

    m.def("w", &GooseFEM::Element::Hex8::Gauss::w, "Return integration point weights");
}

void init_ElementHex8Nodal(py::module& m)
{

    m.def("nip", &GooseFEM::Element::Hex8::Nodal::nip, "Return number of integration point");

    m.def("xi", &GooseFEM::Element::Hex8::Nodal::xi, "Return integration point coordinates");

    m.def("w", &GooseFEM::Element::Hex8::Nodal::w, "Return integration point weights");
}

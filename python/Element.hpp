/**
 * @file
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * \license This project is released under the GNU Public License (GPLv3).
 */

#ifndef PYGOOSEFEM_ELEMENT_H
#define PYGOOSEFEM_ELEMENT_H

#include <GooseFEM/Element.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

template <class C, class P>
void register_Mesh_QuadratureBase(P& cls)
{
    cls.def_property_readonly("nelem", &C::nelem, "Number of elements");
    cls.def_property_readonly("nne", &C::nne, "Number of nodes per element");
    cls.def_property_readonly("ndim", &C::ndim, "Number of spatial dimensions");
    cls.def_property_readonly("tdim", &C::tdim, "Number of spatial dimensions for tensors");
    cls.def_property_readonly("nip", &C::nip, "Number of integration points per element");

    // todo: https://github.com/xtensor-stack/xtensor-python/issues/265
    // cls.def("asTensor",
    //         static_cast<void (C::*)(const xt::pyarray<double>&, xt::pyarray<double>&)
    //         const>(&C::asTensor), "Convert 'qscalar' to 'qtensor' of certain rank",
    //         py::arg("qscalar"),
    //         py::arg("qtensor"));

    cls.def(
        "AsTensor",
        static_cast<xt::xarray<double> (C::*)(size_t, const xt::pyarray<double>&) const>(
            &C::AsTensor),
        py::arg("rank"),
        py::arg("qscalar"));

    cls.def("shape_elemvec", static_cast<std::array<size_t, 3> (C::*)() const>(&C::shape_elemvec));

    cls.def(
        "shape_elemvec",
        static_cast<std::array<size_t, 3> (C::*)(size_t) const>(&C::shape_elemvec),
        py::arg("tdim"));

    cls.def("shape_elemmat", &C::shape_elemmat);

    cls.def(
        "shape_qtensor",
        static_cast<std::vector<size_t> (C::*)(size_t) const>(&C::shape_qtensor),
        py::arg("rank"));

    cls.def(
        "shape_qtensor",
        static_cast<std::vector<size_t> (C::*)(size_t, size_t) const>(&C::shape_qtensor),
        py::arg("rank"),
        py::arg("tdim"));

    cls.def("shape_qscalar", &C::shape_qscalar);

    cls.def("shape_qvector", static_cast<std::array<size_t, 3> (C::*)() const>(&C::shape_qvector));

    cls.def(
        "shape_qvector",
        static_cast<std::array<size_t, 3> (C::*)(size_t) const>(&C::shape_qvector),
        py::arg("tdim"));
}

template <class C, class P>
void register_Mesh_QuadratureBaseCartesian(P& cls)
{
    cls.def("update_x", &C::template update_x<xt::pytensor<double, 3>>, py::arg("x"));

    cls.def_property_readonly("dV", &C::dV, "Integration point volume (qscalar)");

    cls.def(
        "InterpQuad_vector",
        &C::template InterpQuad_vector<xt::pytensor<double, 3>>,
        py::arg("elemvec"));

    cls.def(
        "interpQuad_vector",
        &C::template interpQuad_vector<xt::pytensor<double, 3>, xt::pytensor<double, 3>>,
        py::arg("elemvec"),
        py::arg("qvector"));

    cls.def("GradN_vector", &C::template GradN_vector<xt::pytensor<double, 3>>, py::arg("elemvec"));

    cls.def(
        "gradN_vector",
        &C::template gradN_vector<xt::pytensor<double, 3>, xt::pytensor<double, 4>>,
        py::arg("elemvec"),
        py::arg("qtensor"));

    cls.def(
        "GradN_vector_T", &C::template GradN_vector_T<xt::pytensor<double, 3>>, py::arg("elemvec"));

    cls.def(
        "gradN_vector_T",
        &C::template gradN_vector_T<xt::pytensor<double, 3>, xt::pytensor<double, 4>>,
        py::arg("elemvec"),
        py::arg("qtensor"));

    cls.def(
        "SymGradN_vector",
        &C::template SymGradN_vector<xt::pytensor<double, 3>>,
        py::arg("elemvec"));

    cls.def(
        "symGradN_vector",
        &C::template symGradN_vector<xt::pytensor<double, 3>, xt::pytensor<double, 4>>,
        py::arg("elemvec"),
        py::arg("qtensor"));

    cls.def(
        "Int_N_vector_dV",
        &C::template Int_N_vector_dV<xt::pytensor<double, 3>>,
        py::arg("qvector"));

    cls.def(
        "int_N_vector_dV",
        &C::template int_N_vector_dV<xt::pytensor<double, 3>, xt::pytensor<double, 3>>,
        py::arg("qvector"),
        py::arg("elemvec"));

    cls.def(
        "Int_N_scalar_NT_dV",
        &C::template Int_N_scalar_NT_dV<xt::pytensor<double, 2>>,
        py::arg("qscalar"));

    cls.def(
        "int_N_scalar_NT_dV",
        &C::template int_N_scalar_NT_dV<xt::pytensor<double, 2>, xt::pytensor<double, 3>>,
        py::arg("qscalar"),
        py::arg("elemmat"));

    cls.def(
        "Int_gradN_dot_tensor2_dV",
        &C::template Int_gradN_dot_tensor2_dV<xt::pytensor<double, 4>>,
        py::arg("qtensor"));

    cls.def(
        "int_gradN_dot_tensor2_dV",
        &C::template int_gradN_dot_tensor2_dV<xt::pytensor<double, 4>, xt::pytensor<double, 3>>,
        py::arg("qtensor"),
        py::arg("elemvec"));

    cls.def(
        "Int_gradN_dot_tensor4_dot_gradNT_dV",
        &C::template Int_gradN_dot_tensor4_dot_gradNT_dV<xt::pytensor<double, 6>>,
        py::arg("qtensor"));

    cls.def(
        "int_gradN_dot_tensor4_dot_gradNT_dV",
        &C::template int_gradN_dot_tensor4_dot_gradNT_dV<
            xt::pytensor<double, 6>,
            xt::pytensor<double, 3>>,
        py::arg("qtensor"),
        py::arg("elemmat"));
}

void init_Element(py::module& m)
{

    m.def(
        "asElementVector",
        &GooseFEM::Element::asElementVector,
        "See :cpp:func:`GooseFEM::Element::asElementVector`.",
        py::arg("conn"),
        py::arg("nodevec"));

    m.def(
        "assembleElementVector",
        &GooseFEM::Element::assembleNodeVector,
        "See :cpp:func:`GooseFEM::Element::assembleNodeVector`.",
        py::arg("conn"),
        py::arg("elemvec"));
}

#endif

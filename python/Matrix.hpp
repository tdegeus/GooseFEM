/**
 * @file
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * \license This project is released under the GNU Public License (GPLv3).
 */

#ifndef PYGOOSEFEM_MATRIX_H
#define PYGOOSEFEM_MATRIX_H

#include <Eigen/Eigen>
#include <GooseFEM/Matrix.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

template <class C, class P>
void register_Matrix_MatrixBase(P& cls)
{
    cls.def_property_readonly("nelem", &C::nelem, "Number of elements");
    cls.def_property_readonly("nne", &C::nne, "Number of nodes per element");
    cls.def_property_readonly("nnode", &C::nnode, "Number of nodes");
    cls.def_property_readonly("ndim", &C::ndim, "Number of dimensions");
    cls.def_property_readonly("ndof", &C::ndof, "Number of DOFs");
    cls.def_property_readonly("dofs", &C::dofs, "DOFs [nnode, ndim]");
    cls.def_property_readonly("conn", &C::conn, "Connectivity [nelem, nne]");

    cls.def(
        "assemble",
        &C::template assemble<xt::pytensor<double, 3>>,
        "Assemble from elemmat",
        py::arg("elemmat"));

    cls.def("Todense", &C::Todense, "Return a dense matrix (copy)");

    cls.def(
        "todense",
        &C::template todense<xt::pytensor<double, 2>>,
        "Dense matrix (write to ret)",
        py::arg("ret"));

    cls.def(
        "Dot",
        py::overload_cast<const xt::pytensor<double, 1>&>(&C::Dot, py::const_),
        "Dot product.",
        py::arg("x"));

    cls.def(
        "Dot",
        py::overload_cast<const xt::pytensor<double, 2>&>(&C::Dot, py::const_),
        "Dot product.",
        py::arg("x"));

    cls.def(
        "dot",
        py::overload_cast<const xt::pytensor<double, 1>&, xt::pytensor<double, 1>&>(
            &C::dot, py::const_),
        "Dot product (write to b).",
        py::arg("x"),
        py::arg("b"));

    cls.def(
        "dot",
        py::overload_cast<const xt::pytensor<double, 2>&, xt::pytensor<double, 2>&>(
            &C::dot, py::const_),
        "Dot product (write to b).",
        py::arg("x"),
        py::arg("b"));
}

template <class C, class P>
void register_Matrix_MatrixPartitionedBase(P& cls)
{
    cls.def_property_readonly("nnu", &C::nnu, "Number of unknown DOFs");
    cls.def_property_readonly("nnp", &C::nnp, "Number of prescribed DOFs");
    cls.def_property_readonly("iiu", &C::iiu, "Unknown DOFs");
    cls.def_property_readonly("iip", &C::iip, "Prescribed DOFs");

    cls.def(
        "Reaction",
        py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
            &C::Reaction, py::const_),
        "Return ``b`` with correct right-hand-side",
        py::arg("x"),
        py::arg("b"));

    cls.def(
        "Reaction",
        py::overload_cast<const xt::pytensor<double, 2>&, const xt::pytensor<double, 2>&>(
            &C::Reaction, py::const_),
        "Return ``b`` with correct right-hand-side",
        py::arg("x"),
        py::arg("b"));

    cls.def(
        "reaction",
        py::overload_cast<const xt::pytensor<double, 1>&, xt::pytensor<double, 1>&>(
            &C::reaction, py::const_),
        "Update ``b`` with correct right-hand-side",
        py::arg("x"),
        py::arg("b"));

    cls.def(
        "reaction",
        py::overload_cast<const xt::pytensor<double, 2>&, xt::pytensor<double, 2>&>(
            &C::reaction, py::const_),
        "Update ``b`` with correct right-hand-side",
        py::arg("x"),
        py::arg("b"));

    cls.def(
        "Reaction_p",
        py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
            &C::Reaction_p, py::const_),
        "Return ``b_p``",
        py::arg("x_u"),
        py::arg("x_p"));
}

template <class C, class P>
void register_Matrix_MatrixPartitionedTyingsBase(P& cls)
{
    cls.def_property_readonly("nni", &C::nnu, "Number of independent DOFs");
    cls.def_property_readonly("nnd", &C::nnp, "Number of dependent DOFs");
    cls.def_property_readonly("iii", &C::iiu, "Independent DOFs");
    cls.def_property_readonly("iid", &C::iip, "Dependent DOFs");
}

template <class C, class M, class P>
void register_MatrixSolver_MatrixSolverBase(P& cls)
{
    cls.def(
        "solve",
        py::overload_cast<M&, const xt::pytensor<double, 1>&, xt::pytensor<double, 1>&>(
            &C::template solve<M>),
        "Solve system.",
        py::arg("A"),
        py::arg("b"),
        py::arg("x"));

    cls.def(
        "solve",
        py::overload_cast<M&, const xt::pytensor<double, 2>&, xt::pytensor<double, 2>&>(
            &C::template solve<M>),
        "Solve system.",
        py::arg("A"),
        py::arg("b"),
        py::arg("x"));
}

template <class C, class M, class P>
void register_MatrixSolver_MatrixSolverSingleBase(P& cls)
{
    cls.def(
        "Solve",
        py::overload_cast<M&, const xt::pytensor<double, 1>&>(&C::template Solve<M>),
        "Solve system.",
        py::arg("A"),
        py::arg("b"));

    cls.def(
        "Solve",
        py::overload_cast<M&, const xt::pytensor<double, 2>&>(&C::template Solve<M>),
        "Solve system.",
        py::arg("A"),
        py::arg("b"));
}

template <class C, class M, class P>
void register_MatrixSolver_MatrixSolverPartitionedBase(P& cls)
{
    cls.def(
        "Solve",
        py::overload_cast<M&, const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
            &C::template Solve<M>),
        "Solve system.",
        py::arg("A"),
        py::arg("b"),
        py::arg("x"));

    cls.def(
        "Solve",
        py::overload_cast<M&, const xt::pytensor<double, 2>&, const xt::pytensor<double, 2>&>(
            &C::template Solve<M>),
        "Solve system.",
        py::arg("A"),
        py::arg("b"),
        py::arg("x"));
}

void init_Matrix(py::module& m)
{
    // ---

    py::class_<GooseFEM::Matrix> cls(m, "Matrix");
    register_Matrix_MatrixBase<GooseFEM::Matrix>(cls);

    cls.def(
        py::init<const xt::pytensor<size_t, 2>&, const xt::pytensor<size_t, 2>&>(),
        "See :cpp:class:`GooseFEM::Matrix`.",
        py::arg("conn"),
        py::arg("dofs"));

    cls.def_property_readonly("data", &GooseFEM::Matrix::data);

    cls.def("set", &GooseFEM::Matrix::set, py::arg("rows"), py::arg("cols"), py::arg("matrix"));
    cls.def("add", &GooseFEM::Matrix::add, py::arg("rows"), py::arg("cols"), py::arg("matrix"));

    cls.def("__repr__", [](const GooseFEM::Matrix&) { return "<GooseFEM.Matrix>"; });

    // ---

    py::class_<GooseFEM::MatrixSolver<>> slv(m, "MatrixSolver");
    register_MatrixSolver_MatrixSolverBase<GooseFEM::MatrixSolver<>, GooseFEM::Matrix>(slv);
    register_MatrixSolver_MatrixSolverSingleBase<GooseFEM::MatrixSolver<>, GooseFEM::Matrix>(slv);

    slv.def(py::init<>(), "See :cpp:class:`GooseFEM::MatrixSolver`.");

    slv.def("__repr__", [](const GooseFEM::MatrixSolver<>&) { return "<GooseFEM.MatrixSolver>"; });
}

#endif

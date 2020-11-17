/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <GooseFEM/GooseFEM.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_MatrixPartitioned(py::module& m)
{

    py::class_<GooseFEM::MatrixPartitioned>(m, "MatrixPartitioned")

        .def(
            py::init<
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 1>&>(),
            "Sparse, partitioned, matrix",
            py::arg("conn"),
            py::arg("dofs"),
            py::arg("iip"))

        .def("nelem", &GooseFEM::MatrixPartitioned::nelem, "Number of element")

        .def("nne", &GooseFEM::MatrixPartitioned::nne, "Number of nodes per element")

        .def("nnode", &GooseFEM::MatrixPartitioned::nnode, "Number of nodes")

        .def("ndim", &GooseFEM::MatrixPartitioned::ndim, "Number of dimensions")

        .def("ndof", &GooseFEM::MatrixPartitioned::ndof, "Number of degrees-of-freedom")

        .def("nnu", &GooseFEM::MatrixPartitioned::nnu, "Number of unknown degrees-of-freedom")

        .def("nnp", &GooseFEM::MatrixPartitioned::nnp, "Number of prescribed degrees-of-freedom")

        .def(
            "assemble",
            &GooseFEM::MatrixPartitioned::assemble,
            "Assemble matrix from 'elemmat",
            py::arg("elemmat"))

        .def(
            "set",
            &GooseFEM::MatrixPartitioned::set,
            "Overwrite with a dense (sub-) matrix",
            py::arg("rows"),
            py::arg("cols"),
            py::arg("matrix"))

        .def(
            "add",
            &GooseFEM::MatrixPartitioned::add,
            "Add a dense (sub-) matrix to the current matrix",
            py::arg("rows"),
            py::arg("cols"),
            py::arg("matrix"))

        .def("dofs", &GooseFEM::MatrixPartitioned::dofs, "Return degrees-of-freedom")

        .def("iiu", &GooseFEM::MatrixPartitioned::iiu, "Return unknown degrees-of-freedom")

        .def("iip", &GooseFEM::MatrixPartitioned::iip, "Return prescribed degrees-of-freedom")

        .def(
            "Reaction",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixPartitioned::Reaction, py::const_),
            "Reaction",
            py::arg("x"),
            py::arg("b"))

        .def(
            "Reaction",
            py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&>(
                &GooseFEM::MatrixPartitioned::Reaction, py::const_),
            "Reaction",
            py::arg("x"),
            py::arg("b"))

        .def(
            "Reaction_p",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::MatrixPartitioned::Reaction_p, py::const_),
            "Reaction_p",
            py::arg("x_u"),
            py::arg("x_p"))

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 1>&>(&GooseFEM::MatrixPartitioned::Dot, py::const_),
            "Dot",
            py::arg("x"))

        .def(
            "Dot",
            py::overload_cast<const xt::xtensor<double, 2>&>(&GooseFEM::MatrixPartitioned::Dot, py::const_),
            "Dot",
            py::arg("x"))

        .def("__repr__", [](const GooseFEM::MatrixPartitioned&) {
            return "<GooseFEM.MatrixPartitioned>";
        });

    py::class_<GooseFEM::MatrixPartitionedSolver<>>(m, "MatrixPartitionedSolver")

        .def(py::init<>(), "Sparse, partitioned, matrix solver")

        .def(
            "Solve",
            py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve),
            "Solve",
            py::arg("matrix"),
            py::arg("b"),
            py::arg("x"))

        .def(
            "Solve",
            py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::xtensor<double, 2>&,
                const xt::xtensor<double, 2>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve),
            "Solve",
            py::arg("matrix"),
            py::arg("b"),
            py::arg("x"))

        .def(
            "Solve_u",
            py::overload_cast<
                GooseFEM::MatrixPartitioned&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&>(&GooseFEM::MatrixPartitionedSolver<>::Solve_u),
            "Solve_u",
            py::arg("matrix"),
            py::arg("b_u"),
            py::arg("x_p"))

        .def("__repr__", [](const GooseFEM::MatrixPartitionedSolver<>&) {
            return "<GooseFEM.MatrixPartitionedSolver>";
        });
}

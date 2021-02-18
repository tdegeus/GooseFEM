/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_VectorPartitionedTyings(py::module& m)
{

    py::class_<GooseFEM::VectorPartitionedTyings>(m, "VectorPartitionedTyings")

        .def(
            py::init<
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 2>&,
                const Eigen::SparseMatrix<double>&,
                const Eigen::SparseMatrix<double>&,
                const Eigen::SparseMatrix<double>&>(),
            "Switch between dofval/nodevec/elemvec",
            py::arg("conn"),
            py::arg("dofs"),
            py::arg("Cdu"),
            py::arg("Cdp"),
            py::arg("Cdi"))

        .def("nnu", &GooseFEM::VectorPartitionedTyings::nnu, "Number of unknown DOFs")

        .def("nnp", &GooseFEM::VectorPartitionedTyings::nnp, "Number of prescribed DOFs")

        .def("nni", &GooseFEM::VectorPartitionedTyings::nni, "Number of independent DOFs")

        .def("nnd", &GooseFEM::VectorPartitionedTyings::nnd, "Number of dependent DOFs")

        .def("iiu", &GooseFEM::VectorPartitionedTyings::iiu, "Unknown DOFs")

        .def("iip", &GooseFEM::VectorPartitionedTyings::iip, "Prescribed DOFs")

        .def("iii", &GooseFEM::VectorPartitionedTyings::iii, "Independent DOFs")

        .def("iid", &GooseFEM::VectorPartitionedTyings::iid, "Dependent DOFs")

        .def(
            "AsDofs_i",
            &GooseFEM::VectorPartitionedTyings::AsDofs_i,
            "Set 'dofval",
            py::arg("nodevec"))

        // Overloaded methods from Vector

        .def("nelem", &GooseFEM::VectorPartitionedTyings::nelem, "Number of element")

        .def("nne", &GooseFEM::VectorPartitionedTyings::nne, "Number of nodes per element")

        .def("nnode", &GooseFEM::VectorPartitionedTyings::nnode, "Number of nodes")

        .def("ndim", &GooseFEM::VectorPartitionedTyings::ndim, "Number of dimensions")

        .def("ndof", &GooseFEM::VectorPartitionedTyings::ndof, "Number of degrees-of-freedom")

        .def("dofs", &GooseFEM::VectorPartitionedTyings::dofs, "Return degrees-of-freedom")

        .def(
            "Copy",
            &GooseFEM::VectorPartitionedTyings::Copy,
            py::arg("nodevec_src"),
            py::arg("nodevec_dest"))

        .def(
            "AsDofs",
            py::overload_cast<const xt::xtensor<double, 2>&>(&GooseFEM::VectorPartitionedTyings::AsDofs, py::const_),
            "Set 'dofval",
            py::arg("nodevec"))

        .def(
            "AsDofs",
            py::overload_cast<const xt::xtensor<double, 3>&>(&GooseFEM::VectorPartitionedTyings::AsDofs, py::const_),
            "Set 'dofval",
            py::arg("elemvec"))

        .def(
            "AsNode",
            py::overload_cast<const xt::xtensor<double, 1>&>(&GooseFEM::VectorPartitionedTyings::AsNode, py::const_),
            "Set 'nodevec",
            py::arg("dofval"))

        .def(
            "AsNode",
            py::overload_cast<const xt::xtensor<double, 3>&>(&GooseFEM::VectorPartitionedTyings::AsNode, py::const_),
            "Set 'nodevec",
            py::arg("elemvec"))

        .def(
            "AsElement",
            py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::VectorPartitionedTyings::AsElement, py::const_),
            "Set 'elemvec",
            py::arg("dofval"))

        .def(
            "AsElement",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::VectorPartitionedTyings::AsElement, py::const_),
            "Set 'elemvec",
            py::arg("nodevec"))

        .def(
            "AssembleDofs",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::VectorPartitionedTyings::AssembleDofs, py::const_),
            "Assemble 'dofval'",
            py::arg("nodevec"))

        .def(
            "AssembleDofs",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::VectorPartitionedTyings::AssembleDofs, py::const_),
            "Assemble 'dofval'",
            py::arg("elemvec"))

        .def(
            "AssembleNode",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::VectorPartitionedTyings::AssembleNode, py::const_),
            "Assemble 'nodevec'",
            py::arg("elemvec"))

        .def(
            "ShapeDofval",
            &GooseFEM::VectorPartitionedTyings::ShapeDofval)

        .def(
            "ShapeNodevec",
            &GooseFEM::VectorPartitionedTyings::ShapeNodevec)

        .def(
            "ShapeElemvec",
            &GooseFEM::VectorPartitionedTyings::ShapeElemvec)

        .def(
            "ShapeElemmat",
            &GooseFEM::VectorPartitionedTyings::ShapeElemmat)

        .def(
            "AllocateDofval",
            py::overload_cast<>(&GooseFEM::VectorPartitionedTyings::AllocateDofval, py::const_))

        .def(
            "AllocateDofval",
            py::overload_cast<double>(&GooseFEM::VectorPartitionedTyings::AllocateDofval, py::const_))

        .def(
            "AllocateNodevec",
            py::overload_cast<>(&GooseFEM::VectorPartitionedTyings::AllocateNodevec, py::const_))

        .def(
            "AllocateNodevec",
            py::overload_cast<double>(&GooseFEM::VectorPartitionedTyings::AllocateNodevec, py::const_))

        .def(
            "AllocateElemvec",
            py::overload_cast<>(&GooseFEM::VectorPartitionedTyings::AllocateElemvec, py::const_))

        .def(
            "AllocateElemvec",
            py::overload_cast<double>(&GooseFEM::VectorPartitionedTyings::AllocateElemvec, py::const_))

        .def(
            "AllocateElemmat",
            py::overload_cast<>(&GooseFEM::VectorPartitionedTyings::AllocateElemmat, py::const_))

        .def(
            "AllocateElemmat",
            py::overload_cast<double>(&GooseFEM::VectorPartitionedTyings::AllocateElemmat, py::const_))

        .def("__repr__", [](const GooseFEM::VectorPartitionedTyings&) {
            return "<GooseFEM.VectorPartitionedTyings>";
        });
}

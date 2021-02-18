/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_VectorPartitioned(py::module& m)
{

    py::class_<GooseFEM::VectorPartitioned, GooseFEM::Vector>(m, "VectorPartitioned")

        .def(
            py::init<
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 2>&,
                const xt::xtensor<size_t, 1>&>(),
            "Switch between dofval/nodevec/elemvec",
            py::arg("conn"),
            py::arg("dofs"),
            py::arg("iip"))

        .def("nnu", &GooseFEM::VectorPartitioned::nnu, "Number of unknown degrees-of-freedom")

        .def("nnp", &GooseFEM::VectorPartitioned::nnp, "Number of prescribed degrees-of-freedom")

        .def("iiu", &GooseFEM::VectorPartitioned::iiu, "Return unknown degrees-of-freedom")

        .def("iip", &GooseFEM::VectorPartitioned::iip, "Return prescribed degrees-of-freedom")

        .def(
            "AsDofs",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::AsDofs, py::const_),
            "Set 'dofval",
            py::arg("dofval_u"),
            py::arg("dofval_p"))

        .def(
            "AsDofs_u",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::VectorPartitioned::AsDofs_u, py::const_),
            "Set 'dofval",
            py::arg("nodevec"))

        .def(
            "AsDofs_u",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::VectorPartitioned::AsDofs_u, py::const_),
            "Set 'dofval",
            py::arg("elemvec"))

        .def(
            "AsDofs_p",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::VectorPartitioned::AsDofs_p, py::const_),
            "Set 'dofval",
            py::arg("nodevec"))

        .def(
            "AsDofs_p",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::VectorPartitioned::AsDofs_p, py::const_),
            "Set 'dofval",
            py::arg("elemvec"))

        .def(
            "AsNode",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::AsNode, py::const_),
            "Set 'nodevec",
            py::arg("dofval_u"),
            py::arg("dofval_p"))

        .def(
            "AsElement",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::AsElement, py::const_),
            "Set 'elemvec",
            py::arg("dofval_u"),
            py::arg("dofval_p"))

        .def(
            "AssembleDofs_u",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::VectorPartitioned::AssembleDofs_u, py::const_),
            "Assemble 'dofval'",
            py::arg("nodevec"))

        .def(
            "AssembleDofs_u",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::VectorPartitioned::AssembleDofs_u, py::const_),
            "Assemble 'dofval'",
            py::arg("elemvec"))

        .def(
            "AssembleDofs_p",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::VectorPartitioned::AssembleDofs_p, py::const_),
            "Assemble 'dofval'",
            py::arg("nodevec"))

        .def(
            "AssembleDofs_p",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::VectorPartitioned::AssembleDofs_p, py::const_),
            "Assemble 'dofval'",
            py::arg("elemvec"))

        .def("Copy_u", &GooseFEM::VectorPartitioned::Copy_u, "Copy iiu")

        .def("Copy_p", &GooseFEM::VectorPartitioned::Copy_p, "Copy iip")

        // Overloaded method from Vector

        .def("nelem", &GooseFEM::VectorPartitioned::nelem, "Number of element")

        .def("nne", &GooseFEM::VectorPartitioned::nne, "Number of nodes per element")

        .def("nnode", &GooseFEM::VectorPartitioned::nnode, "Number of nodes")

        .def("ndim", &GooseFEM::VectorPartitioned::ndim, "Number of dimensions")

        .def("ndof", &GooseFEM::VectorPartitioned::ndof, "Number of degrees-of-freedom")

        .def("dofs", &GooseFEM::VectorPartitioned::dofs, "Return degrees-of-freedom")

        .def(
            "Copy",
            &GooseFEM::Vector::Copy,
            py::arg("nodevec_src"),
            py::arg("nodevec_dest"))

        .def(
            "AsDofs",
            py::overload_cast<const xt::xtensor<double, 2>&>(&GooseFEM::Vector::AsDofs, py::const_),
            "Set 'dofval",
            py::arg("nodevec"))

        .def(
            "AsDofs",
            py::overload_cast<const xt::xtensor<double, 3>&>(&GooseFEM::Vector::AsDofs, py::const_),
            "Set 'dofval",
            py::arg("elemvec"))

        .def(
            "AsNode",
            py::overload_cast<const xt::xtensor<double, 1>&>(&GooseFEM::Vector::AsNode, py::const_),
            "Set 'nodevec",
            py::arg("dofval"))

        .def(
            "AsNode",
            py::overload_cast<const xt::xtensor<double, 3>&>(&GooseFEM::Vector::AsNode, py::const_),
            "Set 'nodevec",
            py::arg("elemvec"))

        .def(
            "AsElement",
            py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::Vector::AsElement, py::const_),
            "Set 'elemvec",
            py::arg("dofval"))

        .def(
            "AsElement",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::Vector::AsElement, py::const_),
            "Set 'elemvec",
            py::arg("nodevec"))

        .def(
            "AssembleDofs",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::Vector::AssembleDofs, py::const_),
            "Assemble 'dofval'",
            py::arg("nodevec"))

        .def(
            "AssembleDofs",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Vector::AssembleDofs, py::const_),
            "Assemble 'dofval'",
            py::arg("elemvec"))

        .def(
            "AssembleNode",
            py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Vector::AssembleNode, py::const_),
            "Assemble 'nodevec'",
            py::arg("elemvec"))

        .def(
            "ShapeDofval",
            &GooseFEM::Vector::ShapeDofval)

        .def(
            "ShapeNodevec",
            &GooseFEM::Vector::ShapeNodevec)

        .def(
            "ShapeElemvec",
            &GooseFEM::Vector::ShapeElemvec)

        .def(
            "ShapeElemmat",
            &GooseFEM::Vector::ShapeElemmat)

        .def(
            "AllocateDofval",
            py::overload_cast<>(&GooseFEM::Vector::AllocateDofval, py::const_))

        .def(
            "AllocateDofval",
            py::overload_cast<double>(&GooseFEM::Vector::AllocateDofval, py::const_))

        .def(
            "AllocateNodevec",
            py::overload_cast<>(&GooseFEM::Vector::AllocateNodevec, py::const_))

        .def(
            "AllocateNodevec",
            py::overload_cast<double>(&GooseFEM::Vector::AllocateNodevec, py::const_))

        .def(
            "AllocateElemvec",
            py::overload_cast<>(&GooseFEM::Vector::AllocateElemvec, py::const_))

        .def(
            "AllocateElemvec",
            py::overload_cast<double>(&GooseFEM::Vector::AllocateElemvec, py::const_))

        .def(
            "AllocateElemmat",
            py::overload_cast<>(&GooseFEM::Vector::AllocateElemmat, py::const_))

        .def(
            "AllocateElemmat",
            py::overload_cast<double>(&GooseFEM::Vector::AllocateElemmat, py::const_))

        .def("__repr__", [](const GooseFEM::VectorPartitioned&) {
            return "<GooseFEM.VectorPartitioned>";
        });
}

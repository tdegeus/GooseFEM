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
            "DofsFromParitioned",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::DofsFromParitioned, py::const_),
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
            "NodeFromPartitioned",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::NodeFromPartitioned, py::const_),
            "Set 'nodevec",
            py::arg("dofval_u"),
            py::arg("dofval_p"))

        .def(
            "ElementFromPartitioned",
            py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::ElementFromPartitioned, py::const_),
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

        .def("__repr__", [](const GooseFEM::VectorPartitioned&) {
            return "<GooseFEM.VectorPartitioned>";
        });
}

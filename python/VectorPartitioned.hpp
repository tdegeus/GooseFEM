/**
 * @file
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_VectorPartitioned(py::module& m)
{

    py::class_<GooseFEM::VectorPartitioned, GooseFEM::Vector>(m, "VectorPartitioned")

        .def(
            py::init<
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 2>&,
                const xt::pytensor<size_t, 1>&>(),
            "See :cpp:class:`GooseFEM::VectorPartitioned`.",
            py::arg("conn"),
            py::arg("dofs"),
            py::arg("iip")
        )

        .def_property_readonly("nnu", &GooseFEM::VectorPartitioned::nnu)
        .def_property_readonly("nnp", &GooseFEM::VectorPartitioned::nnp)
        .def_property_readonly("iiu", &GooseFEM::VectorPartitioned::iiu)
        .def_property_readonly("iip", &GooseFEM::VectorPartitioned::iip)
        .def_property_readonly("dofs_is_u", &GooseFEM::VectorPartitioned::dofs_is_u)
        .def_property_readonly("dofs_is_p", &GooseFEM::VectorPartitioned::dofs_is_p)

        .def(
            "DofsFromParitioned",
            py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::DofsFromParitioned, py::const_
            ),
            py::arg("dofval_u"),
            py::arg("dofval_p")
        )

        .def(
            "AsDofs_u",
            py::overload_cast<const xt::pytensor<double, 2>&>(
                &GooseFEM::VectorPartitioned::AsDofs_u, py::const_
            ),
            py::arg("nodevec")
        )

        .def(
            "AsDofs_u",
            py::overload_cast<const xt::pytensor<double, 3>&>(
                &GooseFEM::VectorPartitioned::AsDofs_u, py::const_
            ),
            py::arg("elemvec")
        )

        .def(
            "AsDofs_p",
            py::overload_cast<const xt::pytensor<double, 2>&>(
                &GooseFEM::VectorPartitioned::AsDofs_p, py::const_
            ),
            py::arg("nodevec")
        )

        .def(
            "AsDofs_p",
            py::overload_cast<const xt::pytensor<double, 3>&>(
                &GooseFEM::VectorPartitioned::AsDofs_p, py::const_
            ),
            py::arg("elemvec")
        )

        .def(
            "NodeFromPartitioned",
            py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::NodeFromPartitioned, py::const_
            ),
            py::arg("dofval_u"),
            py::arg("dofval_p")
        )

        .def(
            "ElementFromPartitioned",
            py::overload_cast<const xt::pytensor<double, 1>&, const xt::pytensor<double, 1>&>(
                &GooseFEM::VectorPartitioned::ElementFromPartitioned, py::const_
            ),
            py::arg("dofval_u"),
            py::arg("dofval_p")
        )

        .def(
            "Copy_u",
            &GooseFEM::VectorPartitioned::Copy_u,
            py::arg("nodevec_src"),
            py::arg("nodevec_dest")
        )
        .def(
            "Copy_p",
            &GooseFEM::VectorPartitioned::Copy_p,
            py::arg("nodevec_src"),
            py::arg("nodevec_dest")
        )
        .def(
            "copy_u",
            &GooseFEM::VectorPartitioned::copy_u,
            py::arg("nodevec_src"),
            py::arg("nodevec_dest")
        )
        .def(
            "copy_p",
            &GooseFEM::VectorPartitioned::copy_p,
            py::arg("nodevec_src"),
            py::arg("nodevec_dest")
        )

        .def("__repr__", [](const GooseFEM::VectorPartitioned&) {
            return "<GooseFEM.VectorPartitioned>";
        });
}

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>
#include "../include/GooseFEM/GooseFEM.h"

// =================================================================================================

namespace py = pybind11;

// =================================================================================================

void init_VectorPartitionedTyings(py::module& m)
{

py::class_<GooseFEM::VectorPartitionedTyings>(m, "VectorPartitionedTyings")

    .def(py::init<
            const xt::xtensor<size_t,2>&,
            const xt::xtensor<size_t,2>&,
            const Eigen::SparseMatrix<double>&,
            const Eigen::SparseMatrix<double>&,
            const Eigen::SparseMatrix<double>&>(),
        "Switch between dofval/nodevec/elemvec",
        py::arg("conn"),
        py::arg("dofs"),
        py::arg("Cdu"),
        py::arg("Cdp"),
        py::arg("Cdi"))

    .def("nelem",
        &GooseFEM::VectorPartitionedTyings::nelem,
        "Number of element")

    .def("nne",
        &GooseFEM::VectorPartitionedTyings::nne,
        "Number of nodes per element")

    .def("nnode",
        &GooseFEM::VectorPartitionedTyings::nnode,
        "Number of nodes")

    .def("ndim",
        &GooseFEM::VectorPartitionedTyings::ndim,
        "Number of dimensions")

    .def("ndof",
        &GooseFEM::VectorPartitionedTyings::ndof,
        "Number of degrees-of-freedom")

    .def("nnu",
        &GooseFEM::VectorPartitionedTyings::nnu,
        "Number of unknown degrees-of-freedom")

    .def("nnp",
        &GooseFEM::VectorPartitionedTyings::nnp,
        "Number of prescribed degrees-of-freedom")

    .def("nni",
        &GooseFEM::VectorPartitionedTyings::nni,
        "Number of independent degrees-of-freedom")

    .def("nnd",
        &GooseFEM::VectorPartitionedTyings::nnd,
        "Number of dependent degrees-of-freedom")

    .def("dofs",
        &GooseFEM::VectorPartitionedTyings::dofs,
        "Return degrees-of-freedom")

    .def("iiu",
        &GooseFEM::VectorPartitionedTyings::iiu,
        "Return unknown degrees-of-freedom")

    .def("iip",
        &GooseFEM::VectorPartitionedTyings::iip,
        "Return prescribed degrees-of-freedom")

    .def("iii",
        &GooseFEM::VectorPartitionedTyings::iii,
        "Return independent degrees-of-freedom")

    .def("iid",
        &GooseFEM::VectorPartitionedTyings::iid,
        "Return dependent degrees-of-freedom")

    .def("AsDofs_i",
        &GooseFEM::VectorPartitionedTyings::AsDofs_i,
        "Set 'dofval",
        py::arg("nodevec"))

    .def("AsNode",
        &GooseFEM::VectorPartitionedTyings::AsNode,
        "Set 'nodevec",
        py::arg("dofval"))

    .def("AsElement",
        &GooseFEM::VectorPartitionedTyings::AsElement,
        "Set 'elemvec",
        py::arg("nodevec"))

    .def("AssembleDofs",
        &GooseFEM::VectorPartitionedTyings::AssembleDofs,
        "Assemble 'dofval'",
        py::arg("elemvec"))

    .def("AssembleNode",
        &GooseFEM::VectorPartitionedTyings::AssembleNode,
        "Assemble 'nodevec'",
        py::arg("elemvec"))

    .def("__repr__", [](
        const GooseFEM::VectorPartitionedTyings&){ return "<GooseFEM.VectorPartitionedTyings>"; });

}

// =================================================================================================


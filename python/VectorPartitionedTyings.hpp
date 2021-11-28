/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_VectorPartitionedTyings(py::module& m)
{

    py::class_<GooseFEM::VectorPartitionedTyings, GooseFEM::Vector>(m, "VectorPartitionedTyings")

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
            &GooseFEM::VectorPartitionedTyings::AsDofs_i<xt::xtensor<double, 2>>,
            "Set 'dofval",
            py::arg("nodevec"))

        .def("__repr__", [](const GooseFEM::VectorPartitionedTyings&) {
            return "<GooseFEM.VectorPartitionedTyings>";
        });
}

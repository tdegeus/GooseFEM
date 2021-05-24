/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/Allocate.h>
#include <xtensor/xarray.hpp>
#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_Allocate(py::module& m)
{
    m.def("AsTensor",
          py::overload_cast<
            const xt::xarray<double>&,
            const std::vector<size_t>&>(&GooseFEM::AsTensor<xt::xarray<double>, std::vector<size_t>>),
          "See :cpp:func:`GooseFEM::AsTensor`.",
          py::arg("arg"),
          py::arg("shape"));

    m.def("AsTensor",
          py::overload_cast<
            size_t,
            const xt::xarray<double>&,
            size_t>(&GooseFEM::AsTensor<xt::xarray<double>>),
          "See :cpp:func:`GooseFEM::AsTensor`.",
          py::arg("rank"),
          py::arg("arg"),
          py::arg("n"));

    m.def("as3d",
          &GooseFEM::as3d<xt::xarray<double>>,
          "See :cpp:func:`GooseFEM::as3d`.",
          py::arg("arg"));
}

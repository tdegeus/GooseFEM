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
          static_cast<xt::xarray<double> (*)(
          const xt::xarray<double>&, const std::vector<size_t>&)>(&GooseFEM::AsTensor),
          "See :cpp:func:`GooseFEM::AsTensor`.",
          py::arg("arg"),
          py::arg("shape"));

    m.def("AsTensor",
          static_cast<xt::xarray<double> (*)(
            size_t, const xt::xarray<double>&, size_t)>(&GooseFEM::AsTensor),
          "See :cpp:func:`GooseFEM::AsTensor`.",
          py::arg("rank"),
          py::arg("arg"),
          py::arg("n"));

    m.def("as3d",
          &GooseFEM::as3d<xt::xarray<double>>,
          "See :cpp:func:`GooseFEM::as3d`.",
          py::arg("arg"));
}

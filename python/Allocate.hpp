/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/Allocate.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>

namespace py = pybind11;

void init_Allocate(py::module& m)
{
    m.def(
        "AsTensor",
        static_cast<xt::xarray<double> (*)(const xt::pyarray<double>&, const std::vector<size_t>&)>(
            &GooseFEM::AsTensor),
        "See :cpp:func:`GooseFEM::AsTensor`.",
        py::arg("arg"),
        py::arg("shape"));

    m.def(
        "AsTensor",
        static_cast<xt::xarray<double> (*)(size_t, const xt::pyarray<double>&, size_t)>(
            &GooseFEM::AsTensor),
        "See :cpp:func:`GooseFEM::AsTensor`.",
        py::arg("rank"),
        py::arg("arg"),
        py::arg("n"));

    m.def(
        "as3d",
        &GooseFEM::as3d<xt::pyarray<double>>,
        "See :cpp:func:`GooseFEM::as3d`.",
        py::arg("arg"));
}

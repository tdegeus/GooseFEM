/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/Vector.h>
#include <pybind11/pybind11.h>
#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

namespace py = pybind11;

void init_Vector(py::module& m)
{

    py::class_<GooseFEM::Vector>(m, "Vector")

        .def(py::init<const xt::pytensor<size_t, 2>&, const xt::pytensor<size_t, 2>&>(),
             "Switch between dofval/nodevec/elemvec",
             py::arg("conn"),
             py::arg("dofs"))

        .def("nelem", &GooseFEM::Vector::nelem, "Number of element")
        .def("nne", &GooseFEM::Vector::nne, "Number of nodes per element")
        .def("nnode", &GooseFEM::Vector::nnode, "Number of nodes")
        .def("ndim", &GooseFEM::Vector::ndim, "Number of dimensions")
        .def("ndof", &GooseFEM::Vector::ndof, "Number of degrees-of-freedom")
        .def("conn", &GooseFEM::Vector::conn, "Return connectivity")
        .def("dofs", &GooseFEM::Vector::dofs, "Return degrees-of-freedom")

        .def("copy",
             &GooseFEM::Vector::copy<xt::pyarray<double>>,
             py::arg("nodevec_src"),
             py::arg("nodevec_dest"))

        .def("Copy",
             &GooseFEM::Vector::Copy<xt::pyarray<double>>,
             py::arg("nodevec_src"),
             py::arg("nodevec_dest"))

        .def("AsDofs",
             &GooseFEM::Vector::AsDofs<xt::pyarray<double>>,
             "Convert to 'dofval' (overwrite entries that occur more than once)",
             py::arg("arg"))

        .def("asDofs",
             &GooseFEM::Vector::asDofs<xt::pyarray<double>, xt::pytensor<double, 1>>,
             "Convert to 'dofval' (overwrite entries that occur more than once)",
             py::arg("arg"),
             py::arg("ret"))

        .def("AsNode",
             &GooseFEM::Vector::AsNode<xt::pyarray<double>>,
             "Convert to 'nodevec' (overwrite entries that occur more than once)",
             py::arg("arg"))

        .def("asNode",
             &GooseFEM::Vector::asNode<xt::pyarray<double>, xt::pytensor<double, 2>>,
             "Convert to 'nodevec' (overwrite entries that occur more than once)",
             py::arg("arg"),
             py::arg("ret"))

        .def("AsElement",
             &GooseFEM::Vector::AsElement<xt::pyarray<double>>,
             "Convert to 'elemvec' (overwrite entries that occur more than once)",
             py::arg("arg"))

        .def("asElement",
             &GooseFEM::Vector::asElement<xt::pyarray<double>, xt::pytensor<double, 3>>,
             "Convert to 'elemvec' (overwrite entries that occur more than once)",
             py::arg("arg"),
             py::arg("ret"))

        .def("AssembleDofs",
             &GooseFEM::Vector::AssembleDofs<xt::pyarray<double>>,
             "Assemble to 'dofval' (add entries that occur more than once)",
             py::arg("arg"))

        .def("assembleDofs",
             &GooseFEM::Vector::assembleDofs<xt::pyarray<double>, xt::pytensor<double, 1>>,
             "Assemble to 'dofval' (add entries that occur more than once)",
             py::arg("arg"),
             py::arg("ret"))

        .def("AssembleNode",
            &GooseFEM::Vector::AssembleNode<xt::pyarray<double>>,
             "Assemble to 'nodevec' (add entries that occur more than once)",
             py::arg("arg"))

        .def("assembleNode",
            &GooseFEM::Vector::assembleNode<xt::pyarray<double>, xt::pytensor<double, 2>>,
             "Assemble to 'nodevec' (add entries that occur more than once)",
             py::arg("arg"),
             py::arg("ret"))

        .def("shape_dofval", &GooseFEM::Vector::shape_dofval)
        .def("shape_nodevec", &GooseFEM::Vector::shape_nodevec)
        .def("shape_elemvec", &GooseFEM::Vector::shape_elemvec)
        .def("shape_elemmat", &GooseFEM::Vector::shape_elemmat)

        .def("allocate_dofval",
             py::overload_cast<>(&GooseFEM::Vector::allocate_dofval, py::const_))

        .def("allocate_dofval",
             py::overload_cast<double>(&GooseFEM::Vector::allocate_dofval, py::const_))

        .def("allocate_nodevec",
             py::overload_cast<>(&GooseFEM::Vector::allocate_nodevec, py::const_))

        .def("allocate_nodevec",
             py::overload_cast<double>(&GooseFEM::Vector::allocate_nodevec, py::const_))

        .def("allocate_elemvec",
             py::overload_cast<>(&GooseFEM::Vector::allocate_elemvec, py::const_))

        .def("allocate_elemvec",
             py::overload_cast<double>(&GooseFEM::Vector::allocate_elemvec, py::const_))

        .def("allocate_elemmat",
             py::overload_cast<>(&GooseFEM::Vector::allocate_elemmat, py::const_))

        .def("allocate_elemmat",
             py::overload_cast<double>(&GooseFEM::Vector::allocate_elemmat, py::const_))

        .def("__repr__", [](const GooseFEM::Vector&) { return "<GooseFEM.Vector>"; });
}

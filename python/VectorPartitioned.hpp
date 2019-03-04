/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <pyxtensor/pyxtensor.hpp>

#include "../include/GooseFEM/GooseFEM.h"

// =================================================================================================

namespace py = pybind11;

// =================================================================================================

void init_VectorPartitioned(py::module &m)
{

py::class_<GooseFEM::VectorPartitioned>(m, "VectorPartitioned")

  .def(py::init<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,1> &>(), "Switch between dofval/nodevec/elemvec", py::arg("conn"), py::arg("dofs"), py::arg("iip"))

  .def("nelem", &GooseFEM::VectorPartitioned::nelem, "Return number of element")
  .def("nne"  , &GooseFEM::VectorPartitioned::nne  , "Return number of nodes per element")
  .def("nnode", &GooseFEM::VectorPartitioned::nnode, "Return number of nodes")
  .def("ndim" , &GooseFEM::VectorPartitioned::ndim , "Return number of dimensions")
  .def("ndof" , &GooseFEM::VectorPartitioned::ndof , "Return number of degrees-of-freedom")
  .def("nnu"  , &GooseFEM::VectorPartitioned::nnu  , "Return number of unknown degrees-of-freedom")
  .def("nnp"  , &GooseFEM::VectorPartitioned::nnp  , "Return number of prescribed degrees-of-freedom")

  .def("dofs" , &GooseFEM::VectorPartitioned::dofs , "Return degrees-of-freedom")
  .def("iiu"  , &GooseFEM::VectorPartitioned::iiu  , "Return unknown degrees-of-freedom")
  .def("iip"  , &GooseFEM::VectorPartitioned::iip  , "Return prescribed degrees-of-freedom")

  .def("asDofs", py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&GooseFEM::VectorPartitioned::asDofs, py::const_), "Set 'dofval" , py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asDofs", py::overload_cast<const xt::xtensor<double,2>&                             >(&GooseFEM::VectorPartitioned::asDofs, py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs", py::overload_cast<const xt::xtensor<double,3>&                             >(&GooseFEM::VectorPartitioned::asDofs, py::const_), "Set 'dofval" , py::arg("elemvec"))

  .def("asDofs_u", py::overload_cast<const xt::xtensor<double,2>&>(&GooseFEM::VectorPartitioned::asDofs_u , py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs_u", py::overload_cast<const xt::xtensor<double,3>&>(&GooseFEM::VectorPartitioned::asDofs_u , py::const_), "Set 'dofval" , py::arg("elemvec"))
  .def("asDofs_p", py::overload_cast<const xt::xtensor<double,2>&>(&GooseFEM::VectorPartitioned::asDofs_p , py::const_), "Set 'dofval" , py::arg("nodevec"))
  .def("asDofs_p", py::overload_cast<const xt::xtensor<double,3>&>(&GooseFEM::VectorPartitioned::asDofs_p , py::const_), "Set 'dofval" , py::arg("elemvec"))

  .def("asNode", py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&GooseFEM::VectorPartitioned::asNode, py::const_), "Set 'nodevec", py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asNode", py::overload_cast<const xt::xtensor<double,1>&                             >(&GooseFEM::VectorPartitioned::asNode, py::const_), "Set 'nodevec", py::arg("dofval"))
  .def("asNode", py::overload_cast<const xt::xtensor<double,3>&                             >(&GooseFEM::VectorPartitioned::asNode, py::const_), "Set 'nodevec", py::arg("elemvec"))

  .def("asElement", py::overload_cast<const xt::xtensor<double,1>&,const xt::xtensor<double,1>&>(&GooseFEM::VectorPartitioned::asElement, py::const_), "Set 'elemvec", py::arg("dofval_u"), py::arg("dofval_p"))
  .def("asElement", py::overload_cast<const xt::xtensor<double,1>&                             >(&GooseFEM::VectorPartitioned::asElement, py::const_), "Set 'elemvec", py::arg("dofval"))
  .def("asElement", py::overload_cast<const xt::xtensor<double,2>&                             >(&GooseFEM::VectorPartitioned::asElement, py::const_), "Set 'elemvec", py::arg("nodevec"))

  .def("assembleDofs"  , py::overload_cast<const xt::xtensor<double,2>&>(&GooseFEM::VectorPartitioned::assembleDofs  , py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs"  , py::overload_cast<const xt::xtensor<double,3>&>(&GooseFEM::VectorPartitioned::assembleDofs  , py::const_), "Assemble 'dofval'" , py::arg("elemvec"))
  .def("assembleDofs_u", py::overload_cast<const xt::xtensor<double,2>&>(&GooseFEM::VectorPartitioned::assembleDofs_u, py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs_u", py::overload_cast<const xt::xtensor<double,3>&>(&GooseFEM::VectorPartitioned::assembleDofs_u, py::const_), "Assemble 'dofval'" , py::arg("elemvec"))
  .def("assembleDofs_p", py::overload_cast<const xt::xtensor<double,2>&>(&GooseFEM::VectorPartitioned::assembleDofs_p, py::const_), "Assemble 'dofval'" , py::arg("nodevec"))
  .def("assembleDofs_p", py::overload_cast<const xt::xtensor<double,3>&>(&GooseFEM::VectorPartitioned::assembleDofs_p, py::const_), "Assemble 'dofval'" , py::arg("elemvec"))

  .def("assembleNode"  , py::overload_cast<const xt::xtensor<double,3>&>(&GooseFEM::VectorPartitioned::assembleNode  , py::const_), "Assemble 'nodevec'", py::arg("elemvec"))

  .def("__repr__", [](const GooseFEM::VectorPartitioned &){ return "<GooseFEM.VectorPartitioned>"; });

}

// =================================================================================================


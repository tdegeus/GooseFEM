/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

#include "../src/GooseFEM/GooseFEM.h"

// alias for short-hand notation below
namespace py = pybind11;

PYBIND11_MODULE(GooseFEM, m) {

// =================================================================================================

m.doc() = "Some simple finite element meshes and operations";

// define submodules "mXXX"
py::module mMesh      = m    .def_submodule("Mesh" , "Generic mesh routines"        );
py::module mMeshTri3  = mMesh.def_submodule("Tri3" , "Linear triangular elements"   );
py::module mMeshQuad4 = mMesh.def_submodule("Quad4", "Linear quadrilateral elements");

// ======================================= GooseFEM/Mesh.h ========================================

mMesh.def("elem2node",
  &GooseFEM::Mesh::elem2node,
  "Elements connect to each node: [ number of elements , element numbers ]",
  py::arg("conn")
);

mMesh.def("dofs",
  &GooseFEM::Mesh::dofs,
  "List with DOF-numbers (in sequential order)",
  py::arg("nnode"),
  py::arg("ndim")
);

mMesh.def("renumber",
  py::overload_cast<const GooseFEM::MatS&>(&GooseFEM::Mesh::renumber),
  "Renumber DOF-list to use the lowest possible index",
  py::arg("dofs")
);

mMesh.def("renumber",
  py::overload_cast<const GooseFEM::MatS&,const GooseFEM::ColS&,std::string>(&GooseFEM::Mesh::renumber),
  "Renumber DOF-list to begin or end with 'idx'",
  py::arg("dofs"),
  py::arg("idx"),
  py::arg("location")="end"
);

// ===================================== GooseFEM/MeshQuad4.h =====================================

py::class_<GooseFEM::Mesh::Quad4::Regular>(mMeshQuad4, "Regular")

  .def(
    py::init<size_t,size_t,double>(),
    "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.
  )

  .def("coor"         ,&GooseFEM::Mesh::Quad4::Regular::coor         )
  .def("conn"         ,&GooseFEM::Mesh::Quad4::Regular::conn         )
  .def("nelem"        ,&GooseFEM::Mesh::Quad4::Regular::nelem        )
  .def("nnode"        ,&GooseFEM::Mesh::Quad4::Regular::nnode        )
  .def("nne"          ,&GooseFEM::Mesh::Quad4::Regular::nne          )
  .def("ndim"         ,&GooseFEM::Mesh::Quad4::Regular::ndim         )
  .def("nodesBottom"  ,&GooseFEM::Mesh::Quad4::Regular::nodesBottom  )
  .def("nodesTop"     ,&GooseFEM::Mesh::Quad4::Regular::nodesTop     )
  .def("nodesLeft"    ,&GooseFEM::Mesh::Quad4::Regular::nodesLeft    )
  .def("nodesRight"   ,&GooseFEM::Mesh::Quad4::Regular::nodesRight   )
  .def("nodesPeriodic",&GooseFEM::Mesh::Quad4::Regular::nodesPeriodic)
  .def("nodeOrigin"   ,&GooseFEM::Mesh::Quad4::Regular::nodeOrigin   )
  .def("dofs"         ,&GooseFEM::Mesh::Quad4::Regular::dofs         )
  .def("dofsPeriodic" ,&GooseFEM::Mesh::Quad4::Regular::dofsPeriodic )

  .def("__repr__",
    [](const GooseFEM::Mesh::Quad4::Regular &a){ return "<GooseFEM.Mesh.Quad4.Regular>"; }
  );

// -------------------------------------------------------------------------------------------------

py::class_<GooseFEM::Mesh::Quad4::FineLayer>(mMeshQuad4, "FineLayer")

  .def(
    py::init<size_t,size_t,double,size_t,size_t>(),
    "FineLayer mesh: 'nx' pixels in horizontal direction (length 'Lx'), idem in vertical direction",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.,
    py::arg("nfine")=1,
    py::arg("nskip")=0
  )

  .def("shape"         ,&GooseFEM::Mesh::Quad4::FineLayer::shape         )
  .def("coor"          ,&GooseFEM::Mesh::Quad4::FineLayer::coor          )
  .def("conn"          ,&GooseFEM::Mesh::Quad4::FineLayer::conn          )
  .def("nelem"         ,&GooseFEM::Mesh::Quad4::FineLayer::nelem         )
  .def("nnode"         ,&GooseFEM::Mesh::Quad4::FineLayer::nnode         )
  .def("nne"           ,&GooseFEM::Mesh::Quad4::FineLayer::nne           )
  .def("ndim"          ,&GooseFEM::Mesh::Quad4::FineLayer::ndim          )
  .def("elementsMiddle",&GooseFEM::Mesh::Quad4::FineLayer::elementsMiddle)
  .def("nodesBottom"   ,&GooseFEM::Mesh::Quad4::FineLayer::nodesBottom   )
  .def("nodesTop"      ,&GooseFEM::Mesh::Quad4::FineLayer::nodesTop      )
  .def("nodesLeft"     ,&GooseFEM::Mesh::Quad4::FineLayer::nodesLeft     )
  .def("nodesRight"    ,&GooseFEM::Mesh::Quad4::FineLayer::nodesRight    )
  .def("nodesPeriodic" ,&GooseFEM::Mesh::Quad4::FineLayer::nodesPeriodic )
  .def("nodeOrigin"    ,&GooseFEM::Mesh::Quad4::FineLayer::nodeOrigin    )
  .def("dofs"          ,&GooseFEM::Mesh::Quad4::FineLayer::dofs          )
  .def("dofsPeriodic"  ,&GooseFEM::Mesh::Quad4::FineLayer::dofsPeriodic  )

  .def("__repr__",
    [](const GooseFEM::Mesh::Quad4::FineLayer &a){ return "<GooseFEM.Mesh.Quad4.FineLayer>"; }
  );

// ===================================== GooseFEM/MeshTri3.h ======================================

py::class_<GooseFEM::Mesh::Tri3::Regular>(mMeshTri3, "Regular")

  .def(
    py::init<size_t,size_t,double>(),
    "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge size 'h'",
    py::arg("nx"),
    py::arg("ny"),
    py::arg("h")=1.
  )

  .def("coor"         ,&GooseFEM::Mesh::Tri3::Regular::coor         )
  .def("conn"         ,&GooseFEM::Mesh::Tri3::Regular::conn         )
  .def("nelem"        ,&GooseFEM::Mesh::Tri3::Regular::nelem        )
  .def("nnode"        ,&GooseFEM::Mesh::Tri3::Regular::nnode        )
  .def("nne"          ,&GooseFEM::Mesh::Tri3::Regular::nne          )
  .def("ndim"         ,&GooseFEM::Mesh::Tri3::Regular::ndim         )
  .def("nodesBottom"  ,&GooseFEM::Mesh::Tri3::Regular::nodesBottom  )
  .def("nodesTop"     ,&GooseFEM::Mesh::Tri3::Regular::nodesTop     )
  .def("nodesLeft"    ,&GooseFEM::Mesh::Tri3::Regular::nodesLeft    )
  .def("nodesRight"   ,&GooseFEM::Mesh::Tri3::Regular::nodesRight   )
  .def("nodesPeriodic",&GooseFEM::Mesh::Tri3::Regular::nodesPeriodic)
  .def("nodeOrigin"   ,&GooseFEM::Mesh::Tri3::Regular::nodeOrigin   )
  .def("dofs"         ,&GooseFEM::Mesh::Tri3::Regular::dofs         )
  .def("dofsPeriodic" ,&GooseFEM::Mesh::Tri3::Regular::dofsPeriodic )

  .def("__repr__",
    [](const GooseFEM::Mesh::Tri3::Regular &a){ return "<GooseFEM.Mesh.Tri3.Regular>"; }
  );

// -------------------------------------------------------------------------------------------------

mMeshTri3.def("getOrientation",&GooseFEM::Mesh::Tri3::getOrientation,
  "Get the orientation of each element",
  py::arg("coor"),
  py::arg("conn")
);

// -------------------------------------------------------------------------------------------------

mMeshTri3.def("retriangulate",&GooseFEM::Mesh::Tri3::retriangulate,
  "Re-triangulate existing mesh",
  py::arg("coor"),
  py::arg("conn"),
  py::arg("orientation")=-1
);

// =================================================================================================

}


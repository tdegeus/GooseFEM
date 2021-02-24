/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <GooseFEM/GooseFEM.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pyxtensor/pyxtensor.hpp>

namespace py = pybind11;

void init_MeshQuad4(py::module& m)
{

    py::class_<GooseFEM::Mesh::Quad4::Regular>(m, "Regular")

        .def(
            py::init<size_t, size_t, double>(),
            "Regular mesh: 'nx' pixels in horizontal direction, 'ny' in vertical direction, edge "
            "size 'h'",
            py::arg("nx"),
            py::arg("ny"),
            py::arg("h") = 1.)

        .def("coor", &GooseFEM::Mesh::Quad4::Regular::coor)
        .def("conn", &GooseFEM::Mesh::Quad4::Regular::conn)
        .def("nelem", &GooseFEM::Mesh::Quad4::Regular::nelem)
        .def("nnode", &GooseFEM::Mesh::Quad4::Regular::nnode)
        .def("nne", &GooseFEM::Mesh::Quad4::Regular::nne)
        .def("ndim", &GooseFEM::Mesh::Quad4::Regular::ndim)
        .def("nelx", &GooseFEM::Mesh::Quad4::Regular::nelx)
        .def("nely", &GooseFEM::Mesh::Quad4::Regular::nely)
        .def("h", &GooseFEM::Mesh::Quad4::Regular::h)
        .def("getElementType", &GooseFEM::Mesh::Quad4::Regular::getElementType)
        .def("nodesBottomEdge", &GooseFEM::Mesh::Quad4::Regular::nodesBottomEdge)
        .def("nodesTopEdge", &GooseFEM::Mesh::Quad4::Regular::nodesTopEdge)
        .def("nodesLeftEdge", &GooseFEM::Mesh::Quad4::Regular::nodesLeftEdge)
        .def("nodesRightEdge", &GooseFEM::Mesh::Quad4::Regular::nodesRightEdge)
        .def("nodesBottomOpenEdge", &GooseFEM::Mesh::Quad4::Regular::nodesBottomOpenEdge)
        .def("nodesTopOpenEdge", &GooseFEM::Mesh::Quad4::Regular::nodesTopOpenEdge)
        .def("nodesLeftOpenEdge", &GooseFEM::Mesh::Quad4::Regular::nodesLeftOpenEdge)
        .def("nodesRightOpenEdge", &GooseFEM::Mesh::Quad4::Regular::nodesRightOpenEdge)
        .def("nodesBottomLeftCorner", &GooseFEM::Mesh::Quad4::Regular::nodesBottomLeftCorner)
        .def("nodesBottomRightCorner", &GooseFEM::Mesh::Quad4::Regular::nodesBottomRightCorner)
        .def("nodesTopLeftCorner", &GooseFEM::Mesh::Quad4::Regular::nodesTopLeftCorner)
        .def("nodesTopRightCorner", &GooseFEM::Mesh::Quad4::Regular::nodesTopRightCorner)
        .def("nodesLeftBottomCorner", &GooseFEM::Mesh::Quad4::Regular::nodesLeftBottomCorner)
        .def("nodesLeftTopCorner", &GooseFEM::Mesh::Quad4::Regular::nodesLeftTopCorner)
        .def("nodesRightBottomCorner", &GooseFEM::Mesh::Quad4::Regular::nodesRightBottomCorner)
        .def("nodesRightTopCorner", &GooseFEM::Mesh::Quad4::Regular::nodesRightTopCorner)
        .def("dofs", &GooseFEM::Mesh::Quad4::Regular::dofs)
        .def("nodesPeriodic", &GooseFEM::Mesh::Quad4::Regular::nodesPeriodic)
        .def("nodesOrigin", &GooseFEM::Mesh::Quad4::Regular::nodesOrigin)
        .def("dofsPeriodic", &GooseFEM::Mesh::Quad4::Regular::dofsPeriodic)
        .def("elementgrid", &GooseFEM::Mesh::Quad4::Regular::elementgrid)

        .def("__repr__", [](const GooseFEM::Mesh::Quad4::Regular&) {
            return "<GooseFEM.Mesh.Quad4.Regular>";
        });

    py::class_<GooseFEM::Mesh::Quad4::FineLayer>(m, "FineLayer")

        .def(
            py::init<size_t, size_t, double, size_t>(),
            "FineLayer mesh: 'nx' pixels in horizontal direction (length 'Lx'), idem in vertical "
            "direction",
            py::arg("nx"),
            py::arg("ny"),
            py::arg("h") = 1.,
            py::arg("nfine") = 1)

        .def(
            py::init<const xt::xtensor<double, 2>&, const xt::xtensor<size_t, 2>&>(),
            "Map connectivity to generating FineLayer-object.",
            py::arg("coor"),
            py::arg("conn"))

        .def("coor", &GooseFEM::Mesh::Quad4::FineLayer::coor)
        .def("conn", &GooseFEM::Mesh::Quad4::FineLayer::conn)
        .def("nelem", &GooseFEM::Mesh::Quad4::FineLayer::nelem)
        .def("nnode", &GooseFEM::Mesh::Quad4::FineLayer::nnode)
        .def("nne", &GooseFEM::Mesh::Quad4::FineLayer::nne)
        .def("ndim", &GooseFEM::Mesh::Quad4::FineLayer::ndim)
        .def("nelx", &GooseFEM::Mesh::Quad4::FineLayer::nelx)
        .def("nely", &GooseFEM::Mesh::Quad4::FineLayer::nely)
        .def("h", &GooseFEM::Mesh::Quad4::FineLayer::h)
        .def("elemrow_nhx", &GooseFEM::Mesh::Quad4::FineLayer::elemrow_nhx)
        .def("elemrow_nhy", &GooseFEM::Mesh::Quad4::FineLayer::elemrow_nhy)
        .def("elemrow_nelem", &GooseFEM::Mesh::Quad4::FineLayer::elemrow_nelem)
        .def("getElementType", &GooseFEM::Mesh::Quad4::FineLayer::getElementType)
        .def("elementsMiddleLayer", &GooseFEM::Mesh::Quad4::FineLayer::elementsMiddleLayer)
        .def("elementsLayer", &GooseFEM::Mesh::Quad4::FineLayer::elementsLayer)

        .def(
            "elementgrid_ravel",
            &GooseFEM::Mesh::Quad4::FineLayer::elementgrid_ravel,
            py::arg("rows_range"),
            py::arg("cols_range"))

        .def(
            "elementgrid_around_ravel",
            &GooseFEM::Mesh::Quad4::FineLayer::elementgrid_around_ravel,
            py::arg("element"),
            py::arg("size"),
            py::arg("periodic") = true)

        .def(
            "elementgrid_leftright",
            &GooseFEM::Mesh::Quad4::FineLayer::elementgrid_leftright,
            py::arg("element"),
            py::arg("left"),
            py::arg("right"),
            py::arg("periodic") = true)

        .def("nodesBottomEdge", &GooseFEM::Mesh::Quad4::FineLayer::nodesBottomEdge)
        .def("nodesTopEdge", &GooseFEM::Mesh::Quad4::FineLayer::nodesTopEdge)
        .def("nodesLeftEdge", &GooseFEM::Mesh::Quad4::FineLayer::nodesLeftEdge)
        .def("nodesRightEdge", &GooseFEM::Mesh::Quad4::FineLayer::nodesRightEdge)
        .def("nodesBottomOpenEdge", &GooseFEM::Mesh::Quad4::FineLayer::nodesBottomOpenEdge)
        .def("nodesTopOpenEdge", &GooseFEM::Mesh::Quad4::FineLayer::nodesTopOpenEdge)
        .def("nodesLeftOpenEdge", &GooseFEM::Mesh::Quad4::FineLayer::nodesLeftOpenEdge)
        .def("nodesRightOpenEdge", &GooseFEM::Mesh::Quad4::FineLayer::nodesRightOpenEdge)
        .def("nodesBottomLeftCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesBottomLeftCorner)
        .def("nodesBottomRightCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesBottomRightCorner)
        .def("nodesTopLeftCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesTopLeftCorner)
        .def("nodesTopRightCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesTopRightCorner)
        .def("nodesLeftBottomCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesLeftBottomCorner)
        .def("nodesLeftTopCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesLeftTopCorner)
        .def("nodesRightBottomCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesRightBottomCorner)
        .def("nodesRightTopCorner", &GooseFEM::Mesh::Quad4::FineLayer::nodesRightTopCorner)
        .def("dofs", &GooseFEM::Mesh::Quad4::FineLayer::dofs)
        .def("nodesPeriodic", &GooseFEM::Mesh::Quad4::FineLayer::nodesPeriodic)
        .def("nodesOrigin", &GooseFEM::Mesh::Quad4::FineLayer::nodesOrigin)
        .def("dofsPeriodic", &GooseFEM::Mesh::Quad4::FineLayer::dofsPeriodic)
        .def("roll", &GooseFEM::Mesh::Quad4::FineLayer::roll)

        .def("__repr__", [](const GooseFEM::Mesh::Quad4::FineLayer&) {
            return "<GooseFEM.Mesh.Quad4.FineLayer>";
        });
}

void init_MeshQuad4Map(py::module& m)
{

    py::class_<GooseFEM::Mesh::Quad4::Map::RefineRegular>(m, "RefineRegular")

        .def(
            py::init<const GooseFEM::Mesh::Quad4::Regular&, size_t, size_t>(),
            "Refine a regular mesh",
            py::arg("mesh"),
            py::arg("nx"),
            py::arg("ny"))

        .def("getCoarseMesh", &GooseFEM::Mesh::Quad4::Map::RefineRegular::getCoarseMesh)

        .def("getFineMesh", &GooseFEM::Mesh::Quad4::Map::RefineRegular::getFineMesh)

        .def("getMap", &GooseFEM::Mesh::Quad4::Map::RefineRegular::getMap)

        .def(
            "mapToCoarse",
            py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToCoarse, py::const_))

        .def(
            "mapToCoarse",
            py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToCoarse, py::const_))

        .def(
            "mapToCoarse",
            py::overload_cast<const xt::xtensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToCoarse, py::const_))

        .def("meanToCoarse",
             py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::meanToCoarse<double, 1>, py::const_))

        .def("meanToCoarse",
             py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::meanToCoarse<double, 2>, py::const_))

        .def("meanToCoarse",
             py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::meanToCoarse<double, 3>, py::const_))

        .def("meanToCoarse",
             py::overload_cast<const xt::xtensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::meanToCoarse<double, 4>, py::const_))

        .def("averageToCoarse",
             py::overload_cast<const xt::xtensor<double, 1>&, const xt::xtensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::averageToCoarse<double, 1, double>, py::const_))

        .def("averageToCoarse",
             py::overload_cast<const xt::xtensor<double, 2>&, const xt::xtensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::averageToCoarse<double, 2, double>, py::const_))

        .def("averageToCoarse",
             py::overload_cast<const xt::xtensor<double, 3>&, const xt::xtensor<double, 3>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::averageToCoarse<double, 3, double>, py::const_))

        .def("averageToCoarse",
             py::overload_cast<const xt::xtensor<double, 4>&, const xt::xtensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::averageToCoarse<double, 4, double>, py::const_))

        .def("mapToFine",
             py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToFine<double, 1>, py::const_))

        .def("mapToFine",
             py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToFine<double, 2>, py::const_))

        .def("mapToFine",
             py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToFine<double, 3>, py::const_))

        .def("mapToFine",
             py::overload_cast<const xt::xtensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::RefineRegular::mapToFine<double, 4>, py::const_))

        .def("__repr__", [](const GooseFEM::Mesh::Quad4::Map::RefineRegular&) {
            return "<GooseFEM.Mesh.Quad4.Map.RefineRegular>";
        });

    py::class_<GooseFEM::Mesh::Quad4::Map::FineLayer2Regular>(m, "FineLayer2Regular")

        .def(
            py::init<const GooseFEM::Mesh::Quad4::FineLayer&>(),
            "Map a FineLayer mesh to a Regular mesh",
            py::arg("mesh"))

        .def("getRegularMesh", &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::getRegularMesh)

        .def("getFineLayerMesh", &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::getFineLayerMesh)

        .def("getMap", &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::getMap)

        .def("getMapFraction", &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::getMapFraction)

        .def("mapToRegular",
             py::overload_cast<const xt::xtensor<double, 1>&>(
                &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular<double, 1>, py::const_))

        .def("mapToRegular",
             py::overload_cast<const xt::xtensor<double, 2>&>(
                &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular<double, 2>, py::const_))

        .def("mapToRegular",
             py::overload_cast<const xt::xtensor<double, 3>&>(
                &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular<double, 3>, py::const_))

        .def("mapToRegular",
             py::overload_cast<const xt::xtensor<double, 4>&>(
                &GooseFEM::Mesh::Quad4::Map::FineLayer2Regular::mapToRegular<double, 4>, py::const_))

        .def("__repr__", [](const GooseFEM::Mesh::Quad4::Map::FineLayer2Regular&) {
            return "<GooseFEM.Mesh.Quad4.Map.FineLayer2Regular>";
        });
}

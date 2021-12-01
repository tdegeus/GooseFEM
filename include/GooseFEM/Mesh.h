/**
Generic mesh operations.

\file Mesh.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESH_H
#define GOOSEFEM_MESH_H

#include "config.h"
#include "Vector.h"
#include "MatrixDiagonal.h"

namespace GooseFEM {

/**
Generic mesh operations, and simple mesh definitions.
*/
namespace Mesh {

/**
Enumerator for element-types
*/
enum class ElementType {
    Unknown, ///< Unknown element-type
    Quad4, ///< Quadrilateral: 4-noded element in 2-d
    Hex8, ///< Hexahedron: 8-noded element in 3-d
    Tri3 ///< Triangle: 3-noded element in 2-d
};

/**
Extract the element type based on the connectivity.

\param coor Nodal coordinates [nnode, ndim].
\param conn Connectivity [nelem, nne].
\return ElementType().
*/
template <class S, class T>
inline ElementType defaultElementType(const S& coor, const T& conn);

/**
CRTP base class for regular meshes.
*/
template <class D>
class RegularBase {
public:
    /**
    Underlying type.
    */
    using derived_type = D;

    /**
    Number of elements.
    \return unsigned int
    */
    auto nelem() const;

    /**
    Number of nodes.
    \return unsigned int
    */
    auto nnode() const;

    /**
    Number of nodes-per-element == 4.
    \return unsigned int
    */
    auto nne() const;

    /**
    Number of dimensions == 2.
    \return unsigned int
    */
    auto ndim() const;

    /**
    Number of elements in x-direction == width of the mesh in units of #h.
    \return unsigned int
    */
    auto nelx() const;

    /**
    Number of elements in y-direction == height of the mesh, in units of #h,
    \return unsigned int
    */
    auto nely() const;

    /**
    Linear edge size of one 'block'.
    \return double
    */
    auto h() const;

    /**
    The ElementType().
    \return element type
    */
    auto getElementType() const;

    /**
    Nodal coordinates [#nnode, #ndim].
    \return coordinates per node
    */
    auto coor() const;

    /**
    Connectivity [#nelem, #nne].
    \return nodes per element
    */
    auto conn() const;

    /**
    DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
    \return DOFs per node
    */
    auto dofs() const;

    /**
    DOF-numbers for the case that the periodicity if fully eliminated.
    \return DOF numbers for each node [#nnode, #ndim].
    */
    auto dofsPeriodic() const;

    /**
    Periodic node pairs, in two columns: (independent, dependent).
    \return [ntyings, #ndim].
    */
    auto nodesPeriodic() const;

    /**
    Reference node to use for periodicity, because all corners are tied to it.
    \return Node number.
    */
    auto nodesOrigin() const;

private:
    auto derived_cast() -> derived_type&;
    auto derived_cast() const -> const derived_type&;
};

/**
CRTP base class for regular meshes in 2d.
*/
template <class D>
class RegularBase2d : public RegularBase<D> {
public:
    /**
    Underlying type.
    */
    using derived_type = D;

    /**
    Nodes along the bottom edge (y = 0), in order of increasing x.
    \return List of node numbers.
    */
    auto nodesBottomEdge() const;

    /**
    Nodes along the top edge (y = #nely * #h), in order of increasing x.
    \return List of node numbers.
    */
    auto nodesTopEdge() const;

    /**
    Nodes along the left edge (x = 0), in order of increasing y.
    \return List of node numbers.
    */
    auto nodesLeftEdge() const;

    /**
    Nodes along the right edge (x = #nelx * #h), in order of increasing y.
    \return List of node numbers.
    */
    auto nodesRightEdge() const;

    /**
    Nodes along the bottom edge (y = 0), without the corners (at x = 0 and x = #nelx * #h).
    Same as: nodesBottomEdge()[1: -1].
    \return List of node numbers.
    */
    auto nodesBottomOpenEdge() const;

    /**
    Nodes along the top edge (y = #nely * #h), without the corners (at x = 0 and x = #nelx * #h).
    Same as: nodesTopEdge()[1: -1].
    \return List of node numbers.
    */
    auto nodesTopOpenEdge() const;

    /**
    Nodes along the left edge (x = 0), without the corners (at y = 0 and y = #nely * #h).
    Same as: nodesLeftEdge()[1: -1].
    \return List of node numbers.
    */
    auto nodesLeftOpenEdge() const;

    /**
    Nodes along the right edge (x = #nelx * #h), without the corners (at y = 0 and y = #nely * #h).
    Same as: nodesRightEdge()[1: -1].
    \return List of node numbers.
    */
    auto nodesRightOpenEdge() const;

    /**
    The bottom-left corner node (at x = 0, y = 0).
    Same as nodesBottomEdge()[0] and nodesLeftEdge()[0].
    \return Node number.
    */
    auto nodesBottomLeftCorner() const;

    /**
    The bottom-right corner node (at x = #nelx * #h, y = 0).
    Same as nodesBottomEdge()[-1] and nodesRightEdge()[0].
    \return Node number.
    */
    auto nodesBottomRightCorner() const;

    /**
    The top-left corner node (at x = 0, y = #nely * #h).
    Same as nodesTopEdge()[0] and nodesRightEdge()[-1].
    \return Node number.
    */
    auto nodesTopLeftCorner() const;

    /**
    The top-right corner node (at x = #nelx * #h, y = #nely * #h).
    Same as nodesTopEdge()[-1] and nodesRightEdge()[-1].
    \return Node number.
    */
    auto nodesTopRightCorner() const;

    /**
    Alias of nodesBottomLeftCorner().
    \return Node number.
    */
    auto nodesLeftBottomCorner() const;

    /**
    Alias of nodesTopLeftCorner().
    \return Node number.
    */
    auto nodesLeftTopCorner() const;

    /**
    Alias of nodesBottomRightCorner().
    \return Node number.
    */
    auto nodesRightBottomCorner() const;

    /**
    Alias of nodesTopRightCorner().
    \return Node number.
    */
    auto nodesRightTopCorner() const;

private:
    auto derived_cast() -> derived_type&;
    auto derived_cast() const -> const derived_type&;

    friend class RegularBase<D>;

    xt::xtensor<size_t, 2> nodesPeriodic_impl() const;
    auto nodesOrigin_impl() const;
};

/**
CRTP base class for regular meshes in 3d.
*/
template <class D>
class RegularBase3d : public RegularBase<D> {
public:
    /**
    Underlying type.
    */
    using derived_type = D;

    /**
    Number of elements in y-direction == height of the mesh, in units of #h,
    \return unsigned int
    */
    auto nelz() const;

    /**
    Nodes along the bottom face (y = 0).
    \return List of node numbers.
    */
    auto nodesBottom() const;

    /**
    Nodes along the top face (y = #nely * #h).
    \return List of node numbers.
    */
    auto nodesTop() const;

    /**
    Nodes along the left face (x = 0).
    \return List of node numbers.
    */
    auto nodesLeft() const;

    /**
    Nodes along the right face (x = #nelx * #h).
    \return List of node numbers.
    */
    auto nodesRight() const;

    /**
    Nodes along the front face (z = 0).
    \return List of node numbers.
    */
    auto nodesFront() const;

    /**
    Nodes along the back face (z = #nelz * #h).
    \return List of node numbers.
    */
    auto nodesBack() const;

    /**
    Nodes along the edge at the intersection of the front and bottom faces
    (z = 0 and y = 0).
    \return List of node numbers.
    */
    auto nodesFrontBottomEdge() const;

    /**
    Nodes along the edge at the intersection of the front and top faces
    (z = 0 and y = #nely * #h).
    \return List of node numbers.
    */
    auto nodesFrontTopEdge() const;

    /**
    Nodes along the edge at the intersection of the front and left faces
    (z = 0 and x = 0).
    \return List of node numbers.
    */
    auto nodesFrontLeftEdge() const;

    /**
    Nodes along the edge at the intersection of the front and right faces
    (z = 0 and x = #nelx * #h).
    \return List of node numbers.
    */
    auto nodesFrontRightEdge() const;

    /**
    Nodes along the edge at the intersection of the back and bottom faces
    (z = #nelz * #h and y = #nely * #h).
    \return List of node numbers.
    */
    auto nodesBackBottomEdge() const;

    /**
    Nodes along the edge at the intersection of the back and top faces
    (z = #nelz * #h and x = 0).
    \return List of node numbers.
    */
    auto nodesBackTopEdge() const;

    /**
    Nodes along the edge at the intersection of the back and left faces
    (z = #nelz * #h and x = #nelx * #h).
    \return List of node numbers.
    */
    auto nodesBackLeftEdge() const;

    /**
    Nodes along the edge at the intersection of the back and right faces
    (? = #nelz * #h and ? = ?).
    \return List of node numbers.
    */
    auto nodesBackRightEdge() const;

    /**
    Nodes along the edge at the intersection of the bottom and left faces
    (y = 0 and x = 0).
    \return List of node numbers.
    */
    auto nodesBottomLeftEdge() const;

    /**
    Nodes along the edge at the intersection of the bottom and right faces
    (y = 0 and x = #nelx * #h).
    \return List of node numbers.
    */
    auto nodesBottomRightEdge() const;

    /**
    Nodes along the edge at the intersection of the top and left faces
    (y = 0 and x = #nelx * #h).
    \return List of node numbers.
    */
    auto nodesTopLeftEdge() const;

    /**
    Nodes along the edge at the intersection of the top and right faces
    (y = #nely * #h and x = #nelx * #h).
    \return List of node numbers.
    */
    auto nodesTopRightEdge() const;

    /**
    Alias of nodesFrontBottomEdge()
    \return List of node numbers.
    */
    auto nodesBottomFrontEdge() const;

    /**
    Alias of nodesBackBottomEdge()
    \return List of node numbers.
    */
    auto nodesBottomBackEdge() const;

    /**
    Alias of nodesFrontTopEdge()
    \return List of node numbers.
    */
    auto nodesTopFrontEdge() const;

    /**
    Alias of nodesBackTopEdge()
    \return List of node numbers.
    */
    auto nodesTopBackEdge() const;

    /**
    Alias of nodesBottomLeftEdge()
    \return List of node numbers.
    */
    auto nodesLeftBottomEdge() const;

    /**
    Alias of nodesFrontLeftEdge()
    \return List of node numbers.
    */
    auto nodesLeftFrontEdge() const;

    /**
    Alias of nodesBackLeftEdge()
    \return List of node numbers.
    */
    auto nodesLeftBackEdge() const;

    /**
    Alias of nodesTopLeftEdge()
    \return List of node numbers.
    */
    auto nodesLeftTopEdge() const;

    /**
    Alias of nodesBottomRightEdge()
    \return List of node numbers.
    */
    auto nodesRightBottomEdge() const;

    /**
    Alias of nodesTopRightEdge()
    \return List of node numbers.
    */
    auto nodesRightTopEdge() const;

    /**
    Alias of nodesFrontRightEdge()
    \return List of node numbers.
    */
    auto nodesRightFrontEdge() const;

    /**
    Alias of nodesBackRightEdge()
    \return List of node numbers.
    */
    auto nodesRightBackEdge() const;

    /**
    Nodes along the front face excluding edges.
    Same as different between nodesFront() and
    [nodesFrontBottomEdge(), nodesFrontTopEdge(), nodesFrontLeftEdge(), nodesFrontRightEdge()]
    \return list of node numbers.
    */
    auto nodesFrontFace() const;

    /**
    Nodes along the back face excluding edges.
    Same as different between nodesBack() and
    [nodesBackBottomEdge(), nodesBackTopEdge(), nodesBackLeftEdge(), nodesBackRightEdge()]
    \return list of node numbers.
    */
    auto nodesBackFace() const;

    /**
    Nodes along the left face excluding edges.
    Same as different between nodesLeft() and
    [nodesFrontLeftEdge(), nodesBackLeftEdge(), nodesBottomLeftEdge(), nodesTopLeftEdge()]
    \return list of node numbers.
    */
    auto nodesLeftFace() const;

    /**
    Nodes along the right face excluding edges.
    Same as different between nodesRight() and
    [nodesFrontRightEdge(), nodesBackRightEdge(), nodesBottomRightEdge(), nodesTopRightEdge()]
    \return list of node numbers.
    */
    auto nodesRightFace() const;

    /**
    Nodes along the bottom face excluding edges.
    Same as different between nodesBottom() and
    [nodesBackBottomEdge(), nodesBackTopEdge(), nodesBackLeftEdge(), nodesBackRightEdge()]
    \return list of node numbers.
    */
    auto nodesBottomFace() const;

    /**
    Nodes along the top face excluding edges.
    Same as different between nodesTop() and
    [nodesFrontBottomEdge(), nodesFrontTopEdge(), nodesFrontLeftEdge(), nodesFrontRightEdge()]
    \return list of node numbers.
    */
    auto nodesTopFace() const;

    /**
    Same as nodesFrontBottomEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesFrontBottomOpenEdge() const;

    /**
    Same as nodesFrontTopEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesFrontTopOpenEdge() const;

    /**
    Same as nodesFrontLeftEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesFrontLeftOpenEdge() const;

    /**
    Same as nodesFrontRightEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesFrontRightOpenEdge() const;

    /**
    Same as nodesBackBottomEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesBackBottomOpenEdge() const;

    /**
    Same as nodesBackTopEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesBackTopOpenEdge() const;

    /**
    Same as nodesBackLeftEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesBackLeftOpenEdge() const;

    /**
    Same as nodesBackRightEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesBackRightOpenEdge() const;

    /**
    Same as nodesBottomLeftEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesBottomLeftOpenEdge() const;

    /**
    Same as nodesBottomRightEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesBottomRightOpenEdge() const;

    /**
    Same as nodesTopLeftEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesTopLeftOpenEdge() const;

    /**
    Same as nodesTopRightEdge() but without corners.
    \return List of node numbers.
    */
    auto nodesTopRightOpenEdge() const;

    /**
    Alias of nodesFrontBottomOpenEdge().
    \return List of node numbers.
    */
    auto nodesBottomFrontOpenEdge() const;

    /**
    Alias of nodesBackBottomOpenEdge().
    \return List of node numbers.
    */
    auto nodesBottomBackOpenEdge() const;

    /**
    Alias of nodesFrontTopOpenEdge().
    \return List of node numbers.
    */
    auto nodesTopFrontOpenEdge() const;

    /**
    Alias of nodesBackTopOpenEdge().
    \return List of node numbers.
    */
    auto nodesTopBackOpenEdge() const;

    /**
    Alias of nodesBottomLeftOpenEdge().
    \return List of node numbers.
    */
    auto nodesLeftBottomOpenEdge() const;

    /**
    Alias of nodesFrontLeftOpenEdge().
    \return List of node numbers.
    */
    auto nodesLeftFrontOpenEdge() const;

    /**
    Alias of nodesBackLeftOpenEdge().
    \return List of node numbers.
    */
    auto nodesLeftBackOpenEdge() const;

    /**
    Alias of nodesTopLeftOpenEdge().
    \return List of node numbers.
    */
    auto nodesLeftTopOpenEdge() const;

    /**
    Alias of nodesBottomRightOpenEdge().
    \return List of node numbers.
    */
    auto nodesRightBottomOpenEdge() const;

    /**
    Alias of nodesTopRightOpenEdge().
    \return List of node numbers.
    */
    auto nodesRightTopOpenEdge() const;

    /**
    Alias of nodesFrontRightOpenEdge().
    \return List of node numbers.
    */
    auto nodesRightFrontOpenEdge() const;

    /**
    Alias of nodesBackRightOpenEdge().
    \return List of node numbers.
    */
    auto nodesRightBackOpenEdge() const;

    /**
    Front-Bottom-Left corner node.
    \return Node number.
    */
    auto nodesFrontBottomLeftCorner() const;

    /**
    Front-Bottom-Right corner node.
    \return Node number.
    */
    auto nodesFrontBottomRightCorner() const;

    /**
    Front-Top-Left corner node.
    \return Node number.
    */
    auto nodesFrontTopLeftCorner() const;

    /**
    Front-Top-Right corner node.
    \return Node number.
    */
    auto nodesFrontTopRightCorner() const;

    /**
    Back-Bottom-Left corner node.
    \return Node number.
    */
    auto nodesBackBottomLeftCorner() const;

    /**
    Back-Bottom-Right corner node.
    \return Node number.
    */
    auto nodesBackBottomRightCorner() const;

    /**
    Back-Top-Left corner node.
    \return Node number.
    */
    auto nodesBackTopLeftCorner() const;

    /**
    Back-Top-Right corner node.
    \return Node number.
    */
    auto nodesBackTopRightCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    \return Node number.
    */
    auto nodesFrontLeftBottomCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    \return Node number.
    */
    auto nodesBottomFrontLeftCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    \return Node number.
    */
    auto nodesBottomLeftFrontCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    \return Node number.
    */
    auto nodesLeftFrontBottomCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    \return Node number.
    */
    auto nodesLeftBottomFrontCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    \return Node number.
    */
    auto nodesFrontRightBottomCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    \return Node number.
    */
    auto nodesBottomFrontRightCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    \return Node number.
    */
    auto nodesBottomRightFrontCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    \return Node number.
    */
    auto nodesRightFrontBottomCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    \return Node number.
    */
    auto nodesRightBottomFrontCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    \return Node number.
    */
    auto nodesFrontLeftTopCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    \return Node number.
    */
    auto nodesTopFrontLeftCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    \return Node number.
    */
    auto nodesTopLeftFrontCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    \return Node number.
    */
    auto nodesLeftFrontTopCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    \return Node number.
    */
    auto nodesLeftTopFrontCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    \return Node number.
    */
    auto nodesFrontRightTopCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    \return Node number.
    */
    auto nodesTopFrontRightCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    \return Node number.
    */
    auto nodesTopRightFrontCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    \return Node number.
    */
    auto nodesRightFrontTopCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    \return Node number.
    */
    auto nodesRightTopFrontCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    \return Node number.
    */
    auto nodesBackLeftBottomCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    \return Node number.
    */
    auto nodesBottomBackLeftCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    \return Node number.
    */
    auto nodesBottomLeftBackCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    \return Node number.
    */
    auto nodesLeftBackBottomCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    \return Node number.
    */
    auto nodesLeftBottomBackCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    \return Node number.
    */
    auto nodesBackRightBottomCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    \return Node number.
    */
    auto nodesBottomBackRightCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    \return Node number.
    */
    auto nodesBottomRightBackCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    \return Node number.
    */
    auto nodesRightBackBottomCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    \return Node number.
    */
    auto nodesRightBottomBackCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    \return Node number.
    */
    auto nodesBackLeftTopCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    \return Node number.
    */
    auto nodesTopBackLeftCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    \return Node number.
    */
    auto nodesTopLeftBackCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    \return Node number.
    */
    auto nodesLeftBackTopCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    \return Node number.
    */
    auto nodesLeftTopBackCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    \return Node number.
    */
    auto nodesBackRightTopCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    \return Node number.
    */
    auto nodesTopBackRightCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    \return Node number.
    */
    auto nodesTopRightBackCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    \return Node number.
    */
    auto nodesRightBackTopCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    \return Node number.
    */
    auto nodesRightTopBackCorner() const;

private:
    auto derived_cast() -> derived_type&;
    auto derived_cast() const -> const derived_type&;

    friend class RegularBase<D>;

    xt::xtensor<size_t, 2> nodesPeriodic_impl() const;
    auto nodesOrigin_impl() const;
};

/**
Find overlapping nodes. The output has the following structure:

    [[nodes_from_mesh_a],
     [nodes_from_mesh_b]]

\param coor_a Nodal coordinates of mesh "a" [nnode, ndim].
\param coor_b Nodal coordinates of mesh "b" [nnode, ndim].
\param rtol Relative tolerance for position match.
\param atol Absolute tolerance for position match.
\return Overlapping nodes.
*/
template <class S, class T>
inline xt::xtensor<size_t, 2>
overlapping(const S& coor_a, const T& coor_b, double rtol = 1e-5, double atol = 1e-8);

/**
Stitch two mesh objects, specifying overlapping nodes by hand.
*/
class ManualStitch {
public:
    ManualStitch() = default;

    /**
    \param coor_a Nodal coordinates of mesh "a"  [nnode, ndim].
    \param conn_a Connectivity of mesh "a" [nelem, nne].
    \param overlapping_nodes_a Node-numbers of mesh "a" that overlap with mesh "b" [n].
    \param coor_b Nodal coordinates of mesh "b"  [nnode, ndim].
    \param conn_b Connectivity of mesh "b" [nelem, nne].
    \param overlapping_nodes_b Node-numbers of mesh "b" that overlap with mesh "a" [n].
    \param check_position If ``true`` the nodes are checked for position overlap.
    \param rtol Relative tolerance for check on position overlap.
    \param atol Absolute tolerance for check on position overlap.
    */
    template <class CA, class EA, class NA, class CB, class EB, class NB>
    ManualStitch(
        const CA& coor_a,
        const EA& conn_a,
        const NA& overlapping_nodes_a,
        const CB& coor_b,
        const EB& conn_b,
        const NB& overlapping_nodes_b,
        bool check_position = true,
        double rtol = 1e-5,
        double atol = 1e-8);

    /**
    Number of sub meshes == 2.
    \return unsigned int
    */
    size_t nmesh() const;

    /**
    Number of elements.
    \return unsigned int
    */
    size_t nelem() const;

    /**
    Number of nodes.
    \return unsigned int
    */
    size_t nnode() const;

    /**
    Number of nodes-per-element.
    \return unsigned int
    */
    size_t nne() const;

    /**
    Number of dimensions.
    \return unsigned int
    */
    size_t ndim() const;

    /**
    Nodal coordinates [#nnode, #ndim].
    \return coordinates per node
    */
    xt::xtensor<double, 2> coor() const;

    /**
    Connectivity [#nelem, #nne].
    \return nodes per element
    */
    xt::xtensor<size_t, 2> conn() const;

    /**
    DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
    \return DOFs per node
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    Node-map per sub-mesh.
    \return nodes per mesh
    */
    std::vector<xt::xtensor<size_t, 1>> nodemap() const;

    /**
    Element-map per sub-mesh.
    \return elements per mesh
    */
    std::vector<xt::xtensor<size_t, 1>> elemmap() const;

    /**
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return Node-map for a given mesh.
    */
    xt::xtensor<size_t, 1> nodemap(size_t mesh_index) const;

    /**
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return Element-map for a given mesh.
    */
    xt::xtensor<size_t, 1> elemmap(size_t mesh_index) const;

    /**
    Convert set of node numbers for an original mesh to the stitched mesh.

    \param set List of node numbers.
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return List of node numbers for the stitched mesh.
    */
    template <class T>
    T nodeset(const T& set, size_t mesh_index) const;

    /**
    Convert set of element numbers for an original mesh to the stitched mesh.

    \param set List of element numbers.
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return List of element numbers for the stitched mesh.
    */
    template <class T>
    T elemset(const T& set, size_t mesh_index) const;

private:
    xt::xtensor<double, 2> m_coor;
    xt::xtensor<size_t, 2> m_conn;
    xt::xtensor<size_t, 1> m_map_b;
    size_t m_nnd_a;
    size_t m_nel_a;
    size_t m_nel_b;
};

/**
Stitch mesh objects, automatically searching for overlapping nodes.
*/
class Stitch {
public:
    /**
    \param rtol Relative tolerance for position match.
    \param atol Absolute tolerance for position match.
    */
    Stitch(double rtol = 1e-5, double atol = 1e-8);

    /**
    Add mesh to be stitched.

    \param coor Nodal coordinates [nnode, ndim].
    \param conn Connectivity [nelem, nne].
    */
    template <class C, class E>
    void push_back(const C& coor, const E& conn);

    /**
    Number of sub meshes.
    \return unsigned int
    */
    size_t nmesh() const;

    /**
    Number of elements.
    \return unsigned int
    */
    size_t nelem() const;

    /**
    Number of nodes.
    \return unsigned int
    */
    size_t nnode() const;

    /**
    Number of nodes-per-element.
    \return unsigned int
    */
    size_t nne() const;

    /**
    Number of dimensions.
    \return unsigned int
    */
    size_t ndim() const;

    /**
    Nodal coordinates [#nnode, #ndim].
    \return coordinates per node
    */
    xt::xtensor<double, 2> coor() const;

    /**
    Connectivity [#nelem, #nne].
    \return nodes per element
    */
    xt::xtensor<size_t, 2> conn() const;

    /**
    DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
    \return DOFs per node
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    Node-map per sub-mesh.
    \return nodes per mesh
    */
    std::vector<xt::xtensor<size_t, 1>> nodemap() const;

    /**
    Element-map per sub-mesh.
    \return elements per mesh
    */
    std::vector<xt::xtensor<size_t, 1>> elemmap() const;

    /**
    The node numbers in the stitched mesh that are coming from a specific sub-mesh.

    \param mesh_index Index of the sub-mesh.
    \return List of node numbers.
    */
    xt::xtensor<size_t, 1> nodemap(size_t mesh_index) const;

    /**
    The element numbers in the stitched mesh that are coming from a specific sub-mesh.

    \param mesh_index Index of the sub-mesh.
    \return List of element numbers.
    */
    xt::xtensor<size_t, 1> elemmap(size_t mesh_index) const;

    /**
    Convert set of node-numbers for a sub-mesh to the stitched mesh.

    \param set List of node numbers.
    \param mesh_index Index of the sub-mesh.
    \return List of node numbers for the stitched mesh.
    */
    template <class T>
    T nodeset(const T& set, size_t mesh_index) const;

    /**
    Convert set of element-numbers for a sub-mesh to the stitched mesh.

    \param set List of element numbers.
    \param mesh_index Index of the sub-mesh.
    \return List of element numbers for the stitched mesh.
    */
    template <class T>
    T elemset(const T& set, size_t mesh_index) const;

    /**
    Combine set of node numbers for an original to the final mesh (removes duplicates).

    \param set List of node numbers per mesh.
    \return List of node numbers for the stitched mesh.
    */
    template <class T>
    T nodeset(const std::vector<T>& set) const;

    /** \copydoc nodeset(const std::vector<T>&) const */
    template <class T>
    T nodeset(std::initializer_list<T> set) const;

    /**
    Combine set of element numbers for an original to the final mesh.

    \param set List of element numbers per mesh.
    \return List of element numbers for the stitched mesh.
    */
    template <class T>
    T elemset(const std::vector<T>& set) const;

    /** \copydoc elemset(const std::vector<T>&) const */
    template <class T>
    T elemset(std::initializer_list<T> set) const;

protected:
    xt::xtensor<double, 2> m_coor; ///< Nodal coordinates [#nnode, #ndim]
    xt::xtensor<size_t, 2> m_conn; ///< Connectivity [#nelem, #nne]
    std::vector<xt::xtensor<size_t, 1>> m_map; ///< See nodemap(size_t)
    std::vector<size_t> m_nel; ///< Number of elements per sub-mesh.
    std::vector<size_t> m_el_offset; ///< First element of every sub-mesh.
    double m_rtol; ///< Relative tolerance to find overlapping nodes.
    double m_atol; ///< Absolute tolerance to find overlapping nodes.
};

/**
Vertically stack meshes.
*/
class Vstack : public Stitch {
public:
    /**
    \param check_overlap Check if nodes are overlapping when adding a mesh.
    \param rtol Relative tolerance for position match.
    \param atol Absolute tolerance for position match.
    */
    Vstack(bool check_overlap = true, double rtol = 1e-5, double atol = 1e-8);

    /**
    Add a mesh to the top of the current stack.
    Each time the current `nodes_bottom` are stitched with the then highest `nodes_top`.

    \param coor Nodal coordinates [nnode, ndim].
    \param conn Connectivity [nelem, nne].
    \param nodes_bottom Nodes along the bottom edge [n].
    \param nodes_top Nodes along the top edge [n].
    */
    template <class C, class E, class N>
    void push_back(const C& coor, const E& conn, const N& nodes_bottom, const N& nodes_top);

private:
    std::vector<xt::xtensor<size_t, 1>> m_nodes_bot; ///< Bottom nodes of each mesh (renumbered).
    std::vector<xt::xtensor<size_t, 1>> m_nodes_top; ///< Top nodes of each mesh (renumbered).
    bool m_check_overlap; ///< Check if nodes are overlapping when adding a mesh.
};

/**
Renumber indices to lowest possible index. For example:

\f$ \begin{bmatrix} 0 & 1 \\ 5 & 4 \end{bmatrix} \f$

is renumbered to

\f$ \begin{bmatrix} 0 & 1 \\ 3 & 2 \end{bmatrix} \f$

Or, in pseudo-code, the result of this function is that:

    dofs = renumber(dofs)
    sort(unique(dofs[:])) == range(max(dofs+1))

\note One can use the wrapper function renumber(). This class gives more advanced features.
*/
class Renumber {
public:
    Renumber() = default;

    /**
    \param dofs DOF-numbers.
    */
    template <class T>
    Renumber(const T& dofs);

    /**
    Apply renumbering to other set.

    \param list List of (DOF-)numbers.
    \return Renumbered list of (DOF-)numbers.
    */
    template <class T>
    T apply(const T& list) const;

    /**
    Get the list needed to renumber, e.g.:

        dofs_renumbered(i, j) = index(dofs(i, j))

    \return Renumber-index.
    */
    xt::xtensor<size_t, 1> index() const;

private:
    xt::xtensor<size_t, 1> m_renum;
};

/**
Renumber to lowest possible index (see GooseFEM::Mesh::Renumber).

\param dofs DOF-numbers [nnode, ndim].
\return Renumbered DOF-numbers.
*/
template <class T>
inline T renumber(const T& dofs);

/**
Reorder to lowest possible index, in specific order.

For example for ``Reorder({iiu, iip})`` after reordering:

    iiu = xt::range<size_t>(nnu);
    iip = xt::range<size_t>(nnp) + nnu;
*/
class Reorder {
public:
    Reorder() = default;

    /**
    \param args List of (DOF-)numbers.
    */
    template <class T>
    Reorder(const std::initializer_list<T> args);

    /**
    Apply reordering to other set.

    \param list List of (DOF-)numbers.
    \return Reordered list of (DOF-)numbers.
    */
    template <class T>
    T apply(const T& list) const;

    /**
    Get the list needed to reorder, e.g.:

        dofs_reordered(i, j) = index(dofs(i, j))

    \return Reorder-index.
    */
    xt::xtensor<size_t, 1> index() const;

private:
    xt::xtensor<size_t, 1> m_renum;
};

/**
List with DOF-numbers in sequential order.
The output is a sequential list of DOF-numbers for each vector-component of each node.
For example for 3 nodes in 2 dimensions the output is

\f$ \begin{bmatrix} 0 & 1 \\ 2 & 3 \\ 4 & 5 \end{bmatrix} \f$

\param nnode Number of nodes.
\param ndim Number of dimensions.
\return DOF-numbers.
*/
inline xt::xtensor<size_t, 2> dofs(size_t nnode, size_t ndim);

/**
Number of elements connected to each node.

\param conn Connectivity [nelem, nne].
\return Coordination per node.
*/
template <class E>
inline xt::xtensor<size_t, 1> coordination(const E& conn);

/**
Elements connected to each node.

\param conn Connectivity.
\param sorted If ``true`` the output is sorted.
\return Elements per node.
*/
template <class E>
inline std::vector<std::vector<size_t>> elem2node(const E& conn, bool sorted = true);

/**
Return size of each element edge.

\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Edge-sizes per element.
*/
template <class C, class E>
inline xt::xtensor<double, 2> edgesize(const C& coor, const E& conn, ElementType type);

/**
Return size of each element edge.
The element-type is automatically determined, see defaultElementType().

\param coor Nodal coordinates.
\param conn Connectivity.
\return Edge-sizes per element.
*/
template <class C, class E>
inline xt::xtensor<double, 2> edgesize(const C& coor, const E& conn);

/**
Coordinates of the center of each element.

\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Center of each element.
*/
template <class C, class E>
inline xt::xtensor<double, 2> centers(const C& coor, const E& conn, ElementType type);

/**
Coordinates of the center of each element.
The element-type is automatically determined, see defaultElementType().

\param coor Nodal coordinates.
\param conn Connectivity.
\return Center of each element.
*/
template <class C, class E>
inline xt::xtensor<double, 2> centers(const C& coor, const E& conn);

/**
Convert an element-map to a node-map.

\param elem_map Element-map such that ``new_elvar = elvar[elem_map]``.
\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Node-map such that ``new_nodevar = nodevar[node_map]``
*/
template <class T, class C, class E>
inline xt::xtensor<size_t, 1>
elemmap2nodemap(const T& elem_map, const C& coor, const E& conn, ElementType type);

/**
Convert an element-map to a node-map.
The element-type is automatically determined, see defaultElementType().

\param elem_map Element-map such that ``new_elvar = elvar[elem_map]``.
\param coor Nodal coordinates.
\param conn Connectivity.
\return Node-map such that ``new_nodevar = nodevar[node_map]``
*/
template <class T, class C, class E>
inline xt::xtensor<size_t, 1> elemmap2nodemap(const T& elem_map, const C& coor, const E& conn);

/**
Compute the center of gravity of a mesh.

\tparam C e.g. `xt::xtensor<double, 2>`
\tparam E e.g. `xt::xtensor<size_t, 2>`
\param coor Nodal coordinates `[nnode, ndim]`.
\param conn Connectivity `[nelem, nne]`.
\param type ElementType.
\return Center of gravity `[ndim]`.
*/
template <class C, class E>
inline xt::xtensor<double, 1> center_of_gravity(const C& coor, const E& conn, ElementType type);

/**
Compute the center of gravity of a mesh.

\tparam C e.g. `xt::xtensor<double, 2>`
\tparam E e.g. `xt::xtensor<size_t, 2>`
\param coor Nodal coordinates `[nnode, ndim]`.
\param conn Connectivity `[nelem, nne]`.
\return Center of gravity `[ndim]`.
*/
template <class C, class E>
inline xt::xtensor<double, 1> center_of_gravity(const C& coor, const E& conn);

} // namespace Mesh
} // namespace GooseFEM

#include "Mesh.hpp"

#endif

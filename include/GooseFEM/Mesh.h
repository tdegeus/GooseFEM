/**
Generic mesh operations.

\file Mesh.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESH_H
#define GOOSEFEM_MESH_H

#include "config.h"

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

\param coor Nodal coordinates.
\param conn Connectivity.
\return ElementType().
*/
inline ElementType defaultElementType(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn);

/**
Abstract class for regular meshes in 2-d.
This class does not have a specific element-type in mind, it is used mostly internally
to derive from such that common methods do not have to be reimplementation.
*/
class RegularBase2d {
public:

    RegularBase2d() = default;

    virtual ~RegularBase2d() = default;

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
    Number of nodes-per-element == 4.
    \return unsigned int
    */
    size_t nne() const;

    /**
    Number of dimensions == 2.
    \return unsigned int
    */
    size_t ndim() const;

    /**
    Number of elements in x-direction == width of the mesh in units of #h.
    \return unsigned int
    */
    virtual size_t nelx() const;

    /**
    Number of elements in y-direction == height of the mesh, in units of #h,
    \return unsigned int
    */
    virtual size_t nely() const;

    /**
    Linear edge size of one 'block'.
    \return double
    */
    double h() const;

    /**
    The ElementType().
    \return element type
    */
    virtual ElementType getElementType() const = 0;

    /**
    Nodal coordinates [#nnode, #ndim].
    \return coordinates per node
    */
    virtual xt::xtensor<double, 2> coor() const = 0;

    /**
    Connectivity [#nelem, #nne].
    \return nodes per element
    */
    virtual xt::xtensor<size_t, 2> conn() const = 0;

    /**
    DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
    \return DOFs per node
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    Nodes along the bottom edge (y = 0), in order of increasing x.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomEdge() const = 0;

    /**
    Nodes along the top edge (y = #nely * #h), in order of increasing x.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopEdge() const = 0;

    /**
    Nodes along the left edge (x = 0), in order of increasing y.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesLeftEdge() const = 0;

    /**
    Nodes along the right edge (x = #nelx * #h), in order of increasing y.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesRightEdge() const = 0;

    /**
    Nodes along the bottom edge (y = 0), without the corners (at x = 0 and x = #nelx * #h).
    Same as: nodesBottomEdge()[1: -1].

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomOpenEdge() const = 0;

    /**
    Nodes along the top edge (y = #nely * #h), without the corners (at x = 0 and x = #nelx * #h).
    Same as: nodesTopEdge()[1: -1].

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopOpenEdge() const = 0;

    /**
    Nodes along the left edge (x = 0), without the corners (at y = 0 and y = #nely * #h).
    Same as: nodesLeftEdge()[1: -1].

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesLeftOpenEdge() const = 0;

    /**
    Nodes along the right edge (x = #nelx * #h), without the corners (at y = 0 and y = #nely * #h).
    Same as: nodesRightEdge()[1: -1].

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesRightOpenEdge() const = 0;

    /**
    The bottom-left corner node (at x = 0, y = 0).
    Same as nodesBottomEdge()[0] and nodesLeftEdge()[0].

    \return Node number.
    */
    virtual size_t nodesBottomLeftCorner() const = 0;

    /**
    The bottom-right corner node (at x = #nelx * #h, y = 0).
    Same as nodesBottomEdge()[-1] and nodesRightEdge()[0].

    \return Node number.
    */
    virtual size_t nodesBottomRightCorner() const = 0;

    /**
    The top-left corner node (at x = 0, y = #nely * #h).
    Same as nodesTopEdge()[0] and nodesRightEdge()[-1].

    \return Node number.
    */
    virtual size_t nodesTopLeftCorner() const = 0;

    /**
    The top-right corner node (at x = #nelx * #h, y = #nely * #h).
    Same as nodesTopEdge()[-1] and nodesRightEdge()[-1].

    \return Node number.
    */
    virtual size_t nodesTopRightCorner() const = 0;

    /**
    Alias of nodesBottomLeftCorner().
    */
    size_t nodesLeftBottomCorner() const;

    /**
    Alias of nodesTopLeftCorner().
    */
    size_t nodesLeftTopCorner() const;

    /**
    Alias of nodesBottomRightCorner().
    */
    size_t nodesRightBottomCorner() const;

    /**
    Alias of nodesTopRightCorner().
    */
    size_t nodesRightTopCorner() const;

    /**
    DOF-numbers for the case that the periodicity if fully eliminated. Such that:

        dofs[nodesRightOpenEdge(), :] = dofs[nodesLeftOpenEdge(), :]
        dofs[nodesTopOpenEdge(), :] = dofs[nodesBottomOpenEdge(), :]
        dofs[nodesBottomRightCorner(), :] = dofs[nodesBottomLeftCorner(), :]
        dofs[nodesTopRightCorner(), :] = dofs[nodesBottomLeftCorner(), :]
        dofs[nodesTopLeftCorner(), :] = dofs[nodesBottomLeftCorner(), :]

    \return DOF numbers for each node [#nnode, #ndim].
    */
    xt::xtensor<size_t, 2> dofsPeriodic() const;

    /**
    Periodic node pairs, in two columns: (independent, dependent).

    -   nodesRightOpenEdge() are tied to nodesLeftOpenEdge().
    -   nodesTopOpenEdge() are tied to nodesBottomOpenEdge().
    -   nodesBottomRightCorner() are tied to nodesBottomLeftCorner().
    -   nodesTopRightCorner() are tied to nodesBottomLeftCorner().
    -   nodesTopLeftCorner() are tied to nodesBottomLeftCorner().

    \return [ntyings, #ndim].
    */
    xt::xtensor<size_t, 2> nodesPeriodic() const;

    /**
    Reference node to use for periodicity, because all corners are tied to it.

    \return Alias of nodesBottomLeftCorner().
    */
    size_t nodesOrigin() const;

protected:
    double m_h;     ///< See h()
    size_t m_nelx;  ///< See nelx()
    size_t m_nely;  ///< See nely()
    size_t m_nelem; ///< See nelem()
    size_t m_nnode; ///< See nnode()
    size_t m_nne;   ///< See nne()
    size_t m_ndim;  ///< See ndim()
};

/**
Abstract class for regular meshes in 3-d.
This class does not have a specific element-type in mind, it is used mostly internally
to derive from such that common methods do not have to be reimplementation.
*/
class RegularBase3d {
public:

    RegularBase3d() = default;

    virtual ~RegularBase3d() = default;

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
    Number of nodes-per-element == 4.
    \return unsigned int
    */
    size_t nne() const;

    /**
    Number of dimensions == 2.
    \return unsigned int
    */
    size_t ndim() const;

    /**
    Number of elements in x-direction == width of the mesh in units of #h.
    \return unsigned int
    */
    virtual size_t nelx() const;

    /**
    Number of elements in y-direction == height of the mesh, in units of #h,
    \return unsigned int
    */
    virtual size_t nely() const;

    /**
    Number of elements in y-direction == height of the mesh, in units of #h,
    \return unsigned int
    */
    virtual size_t nelz() const;

    /**
    Linear edge size of one 'block'.
    \return double
    */
    double h() const;

    /**
    The ElementType().
    \return element type
    */
    virtual ElementType getElementType() const = 0;

    /**
    Nodal coordinates [#nnode, #ndim].
    \return coordinates per node
    */
    virtual xt::xtensor<double, 2> coor() const = 0;

    /**
    Connectivity [#nelem, #nne].
    \return nodes per element
    */
    virtual xt::xtensor<size_t, 2> conn() const = 0;

    /**
    DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
    \return DOFs per node
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    Nodes along the bottom face (y = 0).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottom() const = 0;

    /**
    Nodes along the top face (y = #nely * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTop() const = 0;

    /**
    Nodes along the left face (x = 0).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesLeft() const = 0;

    /**
    Nodes along the right face (x = #nelx * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesRight() const = 0;

    /**
    Nodes along the front face (z = 0).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFront() const = 0;

    /**
    Nodes along the back face (z = #nelz * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBack() const = 0;

    /**
    Nodes along the edge at the intersection of the front and bottom faces
    (z = 0 and y = 0).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontBottomEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the front and top faces
    (z = 0 and y = #nely * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontTopEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the front and left faces
    (z = 0 and x = 0).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontLeftEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the front and right faces
    (z = 0 and x = #nelx * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontRightEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the back and bottom faces
    (z = #nelz * #h and y = #nely * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackBottomEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the back and top faces
    (z = #nelz * #h and x = 0).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackTopEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the back and left faces
    (z = #nelz * #h and x = #nelx * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackLeftEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the back and right faces
    (? = #nelz * #h and ? = ?).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackRightEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the bottom and left faces
    (y = 0 and x = 0).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomLeftEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the bottom and right faces
    (y = 0 and x = #nelx * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomRightEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the top and left faces
    (y = 0 and x = #nelx * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopLeftEdge() const = 0;

    /**
    Nodes along the edge at the intersection of the top and right faces
    (y = #nely * #h and x = #nelx * #h).

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopRightEdge() const = 0;

    /**
    Alias of nodesFrontBottomEdge()
    */
    xt::xtensor<size_t, 1> nodesBottomFrontEdge() const;

    /**
    Alias of nodesBackBottomEdge()
    */
    xt::xtensor<size_t, 1> nodesBottomBackEdge() const;

    /**
    Alias of nodesFrontTopEdge()
    */
    xt::xtensor<size_t, 1> nodesTopFrontEdge() const;

    /**
    Alias of nodesBackTopEdge()
    */
    xt::xtensor<size_t, 1> nodesTopBackEdge() const;

    /**
    Alias of nodesBottomLeftEdge()
    */
    xt::xtensor<size_t, 1> nodesLeftBottomEdge() const;

    /**
    Alias of nodesFrontLeftEdge()
    */
    xt::xtensor<size_t, 1> nodesLeftFrontEdge() const;

    /**
    Alias of nodesBackLeftEdge()
    */
    xt::xtensor<size_t, 1> nodesLeftBackEdge() const;

    /**
    Alias of nodesTopLeftEdge()
    */
    xt::xtensor<size_t, 1> nodesLeftTopEdge() const;

    /**
    Alias of nodesBottomRightEdge()
    */
    xt::xtensor<size_t, 1> nodesRightBottomEdge() const;

    /**
    Alias of nodesTopRightEdge()
    */
    xt::xtensor<size_t, 1> nodesRightTopEdge() const;

    /**
    Alias of nodesFrontRightEdge()
    */
    xt::xtensor<size_t, 1> nodesRightFrontEdge() const;

    /**
    Alias of nodesBackRightEdge()
    */
    xt::xtensor<size_t, 1> nodesRightBackEdge() const;


    /**
    Nodes along the front face excluding edges.
    Same as different between nodesFront() and
    [nodesFrontBottomEdge(), nodesFrontTopEdge(), nodesFrontLeftEdge(), nodesFrontRightEdge()]

    \return list of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontFace() const = 0;

    /**
    Nodes along the back face excluding edges.
    Same as different between nodesBack() and
    [nodesBackBottomEdge(), nodesBackTopEdge(), nodesBackLeftEdge(), nodesBackRightEdge()]

    \return list of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackFace() const = 0;

    /**
    Nodes along the left face excluding edges.
    Same as different between nodesLeft() and
    [nodesFrontLeftEdge(), nodesBackLeftEdge(), nodesBottomLeftEdge(), nodesTopLeftEdge()]

    \return list of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesLeftFace() const = 0;

    /**
    Nodes along the right face excluding edges.
    Same as different between nodesRight() and
    [nodesFrontRightEdge(), nodesBackRightEdge(), nodesBottomRightEdge(), nodesTopRightEdge()]

    \return list of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesRightFace() const = 0;

    /**
    Nodes along the bottom face excluding edges.
    Same as different between nodesBottom() and
    [nodesBackBottomEdge(), nodesBackTopEdge(), nodesBackLeftEdge(), nodesBackRightEdge()]

    \return list of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomFace() const = 0;

    /**
    Nodes along the top face excluding edges.
    Same as different between nodesTop() and
    [nodesFrontBottomEdge(), nodesFrontTopEdge(), nodesFrontLeftEdge(), nodesFrontRightEdge()]

    \return list of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopFace() const = 0;

    /**
    Same as nodesFrontBottomEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontBottomOpenEdge() const = 0;

    /**
    Same as nodesFrontTopEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontTopOpenEdge() const = 0;

    /**
    Same as nodesFrontLeftEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontLeftOpenEdge() const = 0;

    /**
    Same as nodesFrontRightEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesFrontRightOpenEdge() const = 0;

    /**
    Same as nodesBackBottomEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackBottomOpenEdge() const = 0;

    /**
    Same as nodesBackTopEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackTopOpenEdge() const = 0;

    /**
    Same as nodesBackLeftEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackLeftOpenEdge() const = 0;

    /**
    Same as nodesBackRightEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBackRightOpenEdge() const = 0;

    /**
    Same as nodesBottomLeftEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomLeftOpenEdge() const = 0;

    /**
    Same as nodesBottomRightEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomRightOpenEdge() const = 0;

    /**
    Same as nodesTopLeftEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopLeftOpenEdge() const = 0;

    /**
    Same as nodesTopRightEdge() but without corners.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopRightOpenEdge() const = 0;

    /**
    Alias of nodesFrontBottomOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesBottomFrontOpenEdge() const;

    /**
    Alias of nodesBackBottomOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesBottomBackOpenEdge() const;

    /**
    Alias of nodesFrontTopOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesTopFrontOpenEdge() const;

    /**
    Alias of nodesBackTopOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesTopBackOpenEdge() const;

    /**
    Alias of nodesBottomLeftOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesLeftBottomOpenEdge() const;

    /**
    Alias of nodesFrontLeftOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesLeftFrontOpenEdge() const;

    /**
    Alias of nodesBackLeftOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesLeftBackOpenEdge() const;

    /**
    Alias of nodesTopLeftOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesLeftTopOpenEdge() const;

    /**
    Alias of nodesBottomRightOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesRightBottomOpenEdge() const;

    /**
    Alias of nodesTopRightOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesRightTopOpenEdge() const;

    /**
    Alias of nodesFrontRightOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesRightFrontOpenEdge() const;

    /**
    Alias of nodesBackRightOpenEdge().
    */
    xt::xtensor<size_t, 1> nodesRightBackOpenEdge() const;

    /**
    Front-Bottom-Left corner node.

    \return Node number.
    */
    virtual size_t nodesFrontBottomLeftCorner() const = 0;

    /**
    Front-Bottom-Right corner node.

    \return Node number.
    */
    virtual size_t nodesFrontBottomRightCorner() const = 0;

    /**
    Front-Top-Left corner node.

    \return Node number.
    */
    virtual size_t nodesFrontTopLeftCorner() const = 0;

    /**
    Front-Top-Right corner node.

    \return Node number.
    */
    virtual size_t nodesFrontTopRightCorner() const = 0;

    /**
    Back-Bottom-Left corner node.

    \return Node number.
    */
    virtual size_t nodesBackBottomLeftCorner() const = 0;

    /**
    Back-Bottom-Right corner node.

    \return Node number.
    */
    virtual size_t nodesBackBottomRightCorner() const = 0;

    /**
    Back-Top-Left corner node.

    \return Node number.
    */
    virtual size_t nodesBackTopLeftCorner() const = 0;

    /**
    Back-Top-Right corner node.

    \return Node number.
    */
    virtual size_t nodesBackTopRightCorner() const = 0;

    /**
    Alias of nodesFrontBottomLeftCorner().
    */
    size_t nodesFrontLeftBottomCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    */
    size_t nodesBottomFrontLeftCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    */
    size_t nodesBottomLeftFrontCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    */
    size_t nodesLeftFrontBottomCorner() const;

    /**
    Alias of nodesFrontBottomLeftCorner().
    */
    size_t nodesLeftBottomFrontCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    */
    size_t nodesFrontRightBottomCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    */
    size_t nodesBottomFrontRightCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    */
    size_t nodesBottomRightFrontCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    */
    size_t nodesRightFrontBottomCorner() const;

    /**
    Alias of nodesFrontBottomRightCorner().
    */
    size_t nodesRightBottomFrontCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    */
    size_t nodesFrontLeftTopCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    */
    size_t nodesTopFrontLeftCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    */
    size_t nodesTopLeftFrontCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    */
    size_t nodesLeftFrontTopCorner() const;

    /**
    Alias of nodesFrontTopLeftCorner().
    */
    size_t nodesLeftTopFrontCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    */
    size_t nodesFrontRightTopCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    */
    size_t nodesTopFrontRightCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    */
    size_t nodesTopRightFrontCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    */
    size_t nodesRightFrontTopCorner() const;

    /**
    Alias of nodesFrontTopRightCorner().
    */
    size_t nodesRightTopFrontCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    */
    size_t nodesBackLeftBottomCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    */
    size_t nodesBottomBackLeftCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    */
    size_t nodesBottomLeftBackCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    */
    size_t nodesLeftBackBottomCorner() const;

    /**
    Alias of nodesBackBottomLeftCorner().
    */
    size_t nodesLeftBottomBackCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    */
    size_t nodesBackRightBottomCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    */
    size_t nodesBottomBackRightCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    */
    size_t nodesBottomRightBackCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    */
    size_t nodesRightBackBottomCorner() const;

    /**
    Alias of nodesBackBottomRightCorner().
    */
    size_t nodesRightBottomBackCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    */
    size_t nodesBackLeftTopCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    */
    size_t nodesTopBackLeftCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    */
    size_t nodesTopLeftBackCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    */
    size_t nodesLeftBackTopCorner() const;

    /**
    Alias of nodesBackTopLeftCorner().
    */
    size_t nodesLeftTopBackCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    */
    size_t nodesBackRightTopCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    */
    size_t nodesTopBackRightCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    */
    size_t nodesTopRightBackCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    */
    size_t nodesRightBackTopCorner() const;

    /**
    Alias of nodesBackTopRightCorner().
    */
    size_t nodesRightTopBackCorner() const;

    /**
    DOF-numbers for the case that the periodicity if fully eliminated. Such that:

        dofs[nodesBackFace(), :] = dofs[nodesFrontFace(), :]
        dofs[nodesRightFace(), :] = dofs[nodesLeftFace(), :]
        dofs[nodesTopFace(), :] = dofs[nodesBottomFace(), :]
        dofs[nodesBackBottomOpenEdge(), :] = dofs[nodesFrontBottomOpenEdge(), :]
        dofs[nodesBackTopOpenEdge(), :] = dofs[nodesFrontBottomOpenEdge(), :]
        dofs[nodesFrontTopOpenEdge(), :] = dofs[nodesFrontBottomOpenEdge(), :]
        dofs[nodesBottomRightOpenEdge(), :] = dofs[nodesBottomLeftOpenEdge(), :]
        dofs[nodesTopRightOpenEdge(), :] = dofs[nodesBottomLeftOpenEdge(), :]
        dofs[nodesTopLeftOpenEdge(), :] = dofs[nodesBottomLeftOpenEdge(), :]
        dofs[nodesFrontRightOpenEdge(), :] = dofs[nodesFrontLeftOpenEdge(), :]
        dofs[nodesBackRightOpenEdge(), :] = dofs[nodesFrontLeftOpenEdge(), :]
        dofs[nodesBackLeftOpenEdge(), :] = dofs[nodesFrontLeftOpenEdge(), :]
        dofs[nodesFrontBottomRightCorner(), :] = dofs[nodesFrontBottomLeftCorner(), :]
        dofs[nodesBackBottomRightCorner(), :] = dofs[nodesFrontBottomLeftCorner(), :]
        dofs[nodesBackBottomLeftCorner(), :] = dofs[nodesFrontBottomLeftCorner(), :]
        dofs[nodesFrontTopLeftCorner(), :] = dofs[nodesFrontBottomLeftCorner(), :]
        dofs[nodesFrontTopRightCorner(), :] = dofs[nodesFrontBottomLeftCorner(), :]
        dofs[nodesBackTopRightCorner(), :] = dofs[nodesFrontBottomLeftCorner(), :]
        dofs[nodesBackTopLeftCorner(), :] = dofs[nodesFrontBottomLeftCorner(), :]

    \return DOF numbers for each node [#nnode, #ndim].
    */
    xt::xtensor<size_t, 2> dofsPeriodic() const;

    /**
    Periodic node pairs, in two columns: (independent, dependent).

    -   nodesBackFace() are tied to nodesFrontFace().
    -   nodesRightFace() are tied to nodesLeftFace().
    -   nodesTopFace() are tied to nodesBottomFace().
    -   nodesBackBottomOpenEdge() are tied to nodesFrontBottomOpenEdge().
    -   nodesBackTopOpenEdge() are tied to nodesFrontBottomOpenEdge().
    -   nodesFrontTopOpenEdge() are tied to nodesFrontBottomOpenEdge().
    -   nodesBottomRightOpenEdge() are tied to nodesBottomLeftOpenEdge().
    -   nodesTopRightOpenEdge() are tied to nodesBottomLeftOpenEdge().
    -   nodesTopLeftOpenEdge() are tied to nodesBottomLeftOpenEdge().
    -   nodesFrontRightOpenEdge() are tied to nodesFrontLeftOpenEdge().
    -   nodesBackRightOpenEdge() are tied to nodesFrontLeftOpenEdge().
    -   nodesBackLeftOpenEdge() are tied to nodesFrontLeftOpenEdge().
    -   nodesFrontBottomRightCorner() are tied to nodesFrontBottomLeftCorner().
    -   nodesBackBottomRightCorner() are tied to nodesFrontBottomLeftCorner().
    -   nodesBackBottomLeftCorner() are tied to nodesFrontBottomLeftCorner().
    -   nodesFrontTopLeftCorner() are tied to nodesFrontBottomLeftCorner().
    -   nodesFrontTopRightCorner() are tied to nodesFrontBottomLeftCorner().
    -   nodesBackTopRightCorner() are tied to nodesFrontBottomLeftCorner().
    -   nodesBackTopLeftCorner() are tied to nodesFrontBottomLeftCorner().

    \return [ntyings, #ndim].
    */
    xt::xtensor<size_t, 2> nodesPeriodic() const;

    /**
    Reference node to use for periodicity, because all corners are tied to it.

    \return Alias of nodesFrontBottomLeftCorner().
    */
    size_t nodesOrigin() const;

protected:
    double m_h;     ///< See h()
    size_t m_nelx;  ///< See nelx()
    size_t m_nely;  ///< See nely()
    size_t m_nelz;  ///< See nely()
    size_t m_nelem; ///< See nelem()
    size_t m_nnode; ///< See nnode()
    size_t m_nne;   ///< See nne()
    size_t m_ndim;  ///< See ndim()
};

/**
Find overlapping nodes. The output has the following structure:

    [[nodes_from_mesh_a],
     [nodes_from_mesh_b]]

\param coor_a Nodal coordinates of mesh "a".
\param coor_b Nodal coordinates of mesh "b".
\param rtol Relative tolerance for position match.
\param atol Absolute tolerance for position match.
\return Overlapping nodes.
*/
inline xt::xtensor<size_t, 2> overlapping(
    const xt::xtensor<double, 2>& coor_a,
    const xt::xtensor<double, 2>& coor_b,
    double rtol = 1e-5,
    double atol = 1e-8);

/**
Stitch two mesh objects, specifying overlapping nodes by hand.
*/
class ManualStitch {
public:
    ManualStitch() = default;

    /**
    \param coor_a Nodal coordinates of mesh "a".
    \param conn_a Connectivity of mesh "a".
    \param overlapping_nodes_a Node-numbers of mesh "a" that overlap with mesh "b".
    \param coor_b Nodal coordinates of mesh "b".
    \param conn_b Connectivity of mesh "b".
    \param overlapping_nodes_b Node-numbers of mesh "b" that overlap with mesh "a".
    \param check_position If ``true`` the nodes are checked for position overlap.
    \param rtol Relative tolerance for check on position overlap.
    \param atol Absolute tolerance for check on position overlap.
    */
    ManualStitch(
        const xt::xtensor<double, 2>& coor_a,
        const xt::xtensor<size_t, 2>& conn_a,
        const xt::xtensor<size_t, 1>& overlapping_nodes_a,
        const xt::xtensor<double, 2>& coor_b,
        const xt::xtensor<size_t, 2>& conn_b,
        const xt::xtensor<size_t, 1>& overlapping_nodes_b,
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
    xt::xtensor<size_t, 1> nodeset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const;

    /**
    Convert set of element numbers for an original mesh to the stitched mesh.

    \param set List of element numbers.
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return List of element numbers for the stitched mesh.
    */
    xt::xtensor<size_t, 1> elemset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const;

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

    \param coor Nodal coordinates.
    \param conn Connectivity.
    */
    void push_back(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn);

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
    xt::xtensor<size_t, 1> nodeset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const;

    /**
    Convert set of element-numbers for a sub-mesh to the stitched mesh.

    \param set List of element numbers.
    \param mesh_index Index of the sub-mesh.
    \return List of element numbers for the stitched mesh.
    */
    xt::xtensor<size_t, 1> elemset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const;

    /**
    Combine set of node numbers for an original to the final mesh (removes duplicates).

    \param set List of node numbers per mesh.
    \return List of node numbers for the stitched mesh.
    */
    xt::xtensor<size_t, 1> nodeset(const std::vector<xt::xtensor<size_t, 1>>& set) const;

    /**
    Combine set of element numbers for an original to the final mesh.

    \param set List of element numbers per mesh.
    \return List of element numbers for the stitched mesh.
    */
    xt::xtensor<size_t, 1> elemset(const std::vector<xt::xtensor<size_t, 1>>& set) const;

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

    \param coor Nodal coordinates.
    \param conn Connectivity.
    \param nodes_bottom Nodes along the bottom edge.
    \param nodes_top Nodes along the top edge.
    */
    void push_back(
        const xt::xtensor<double, 2>& coor,
        const xt::xtensor<size_t, 2>& conn,
        const xt::xtensor<size_t, 1>& nodes_bottom,
        const xt::xtensor<size_t, 1>& nodes_top);

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
    Get renumbered DOFs (same as ``Renumber::apply(dofs)``).

    \param dofs List of (DOF-)numbers.
    \return Renumbered list of (DOF-)numbers.
    */
    [[deprecated]]
    xt::xtensor<size_t, 2> get(const xt::xtensor<size_t, 2>& dofs) const;

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

\param dofs DOF-numbers.
\return Renumbered DOF-numbers.
*/
inline xt::xtensor<size_t, 2> renumber(const xt::xtensor<size_t, 2>& dofs);

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
    Reorder(const std::initializer_list<xt::xtensor<size_t, 1>> args);

    /**
    Get reordered DOFs (same as ``Reorder::apply(dofs)``).

    \param dofs List of (DOF-)numbers.
    \return Reordered list of (DOF-)numbers.
    */
    [[deprecated]]
    xt::xtensor<size_t, 2> get(const xt::xtensor<size_t, 2>& dofs) const;

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

\param conn Connectivity.
\return Coordination per node.
*/
inline xt::xtensor<size_t, 1> coordination(const xt::xtensor<size_t, 2>& conn);

/**
Elements connected to each node.

\param conn Connectivity.
\param sorted If ``true`` the output is sorted.
\return Elements per node.
*/
inline std::vector<std::vector<size_t>> elem2node(
    const xt::xtensor<size_t, 2>& conn,
    bool sorted = true);

/**
Return size of each element edge.

\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Edge-sizes per element.
*/
inline xt::xtensor<double, 2> edgesize(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn,
    ElementType type);

/**
Return size of each element edge.
The element-type is automatically determined, see defaultElementType().

\param coor Nodal coordinates.
\param conn Connectivity.
\return Edge-sizes per element.
*/
inline xt::xtensor<double, 2> edgesize(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn);

/**
Coordinates of the center of each element.

\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Center of each element.
*/
inline xt::xtensor<double, 2> centers(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn,
    ElementType type);

/**
Coordinates of the center of each element.
The element-type is automatically determined, see defaultElementType().

\param coor Nodal coordinates.
\param conn Connectivity.
\return Center of each element.
*/
inline xt::xtensor<double, 2> centers(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn);

/**
Convert an element-map to a node-map.

\param elem_map Element-map such that ``new_elvar = elvar[elem_map]``.
\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Node-map such that ``new_nodevar = nodevar[node_map]``
*/
inline xt::xtensor<size_t, 1> elemmap2nodemap(
    const xt::xtensor<size_t, 1>& elem_map,
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn,
    ElementType type);

/**
Convert an element-map to a node-map.
The element-type is automatically determined, see defaultElementType().

\param elem_map Element-map such that ``new_elvar = elvar[elem_map]``.
\param coor Nodal coordinates.
\param conn Connectivity.
\return Node-map such that ``new_nodevar = nodevar[node_map]``
*/
inline xt::xtensor<size_t, 1> elemmap2nodemap(
    const xt::xtensor<size_t, 1>& elem_map,
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn);

} // namespace Mesh
} // namespace GooseFEM

#include "Mesh.hpp"

#endif

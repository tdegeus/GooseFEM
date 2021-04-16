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
Base-class for regular meshes in 2-d.
This class does not have a specific element-type in mind, it is used mostly internally
to derive from such that common methods do not have to be reimplementation.
*/
class RegularBase2d {
public:

    RegularBase2d() = default;

    /**
    \return Number of elements.
    */
    size_t nelem() const;

    /**
    \return Number of nodes.
    */
    size_t nnode() const;

    /**
    \return Number of nodes-per-element == 4.
    */
    size_t nne() const;

    /**
    \return Number of dimensions == 2.
    */
    size_t ndim() const;

    /**
    \return Number of elements in x-direction == width of the mesh in units of #h.
    */
    virtual size_t nelx() const;

    /**
    \return Number of elements in y-direction == height of the mesh, in units of #h,
    */
    virtual size_t nely() const;

    /**
    \return Linear edge size of one 'block'.
    */
    double h() const;

    /**
    \return The ElementType().
    */
    virtual ElementType getElementType() const;

    /**
    \return Nodal coordinates [#nnode, #ndim].
    */
    virtual xt::xtensor<double, 2> coor() const;

    /**
    \return Connectivity [#nelem, #nne].
    */
    virtual xt::xtensor<size_t, 2> conn() const;

    /**
    \return DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
    */
    virtual xt::xtensor<size_t, 2> dofs() const;

    /**
    Nodes along the bottom edge (y = 0), in order of increasing x.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomEdge() const;

    /**
    Nodes along the top edge (y = #nely * #h), in order of increasing x.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopEdge() const;

    /**
    Nodes along the left edge (x = 0), in order of increasing y.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesLeftEdge() const;

    /**
    Nodes along the right edge (x = #nelx * #h), in order of increasing y.

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesRightEdge() const;

    /**
    Nodes along the bottom edge (y = 0), without the corners (at x = 0 and x = #nelx * #h).
    Same as: nodesBottomEdge()[1: -1].

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesBottomOpenEdge() const;

    /**
    Nodes along the top edge (y = #nely * #h), without the corners (at x = 0 and x = #nelx * #h).
    Same as: nodesTopEdge()[1: -1].

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesTopOpenEdge() const;

    /**
    Nodes along the left edge (x = 0), without the corners (at y = 0 and y = #nely * #h).
    Same as: nodesLeftEdge()[1: -1].

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesLeftOpenEdge() const;

    /**
    Nodes along the right edge (x = #nelx * #h), without the corners (at y = 0 and y = #nely * #h).
    Same as: nodesRightEdge()[1: -1].

    \return List of node numbers.
    */
    virtual xt::xtensor<size_t, 1> nodesRightOpenEdge() const;

    /**
    The bottom-left corner node (at x = 0, y = 0).
    Same as nodesBottomEdge()[0] and nodesLeftEdge()[0].

    \return Node number.
    */
    virtual size_t nodesBottomLeftCorner() const;

    /**
    The bottom-right corner node (at x = #nelx * #h, y = 0).
    Same as nodesBottomEdge()[-1] and nodesRightEdge()[0].

    \return Node number.
    */
    virtual size_t nodesBottomRightCorner() const;

    /**
    The top-left corner node (at x = 0, y = #nely * #h).
    Same as nodesTopEdge()[0] and nodesRightEdge()[-1].

    \return Node number.
    */
    virtual size_t nodesTopLeftCorner() const;

    /**
    The top-right corner node (at x = #nelx * #h, y = #nely * #h).
    Same as nodesTopEdge()[-1] and nodesRightEdge()[-1].

    \return Node number.
    */
    virtual size_t nodesTopRightCorner() const;

    /**
    \return Alias of nodesBottomLeftCorner().
    */
    size_t nodesLeftBottomCorner() const;

    /**
    \return Alias of nodesTopLeftCorner().
    */
    size_t nodesLeftTopCorner() const;

    /**
    \return Alias of nodesBottomRightCorner().
    */
    size_t nodesRightBottomCorner() const;

    /**
    \return Alias of nodesTopRightCorner().
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
    Periodic node pairs, in two columns: (independent, dependent)

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
    \return Number of sub meshes == 2.
    */
    size_t nmesh() const;

    /**
    \return Number of elements.
    */
    size_t nelem() const;

    /**
    \return Number of nodes.
    */
    size_t nnode() const;

    /**
    \return Number of nodes-per-element.
    */
    size_t nne() const;

    /**
    \return Number of dimensions.
    */
    size_t ndim() const;

    /**
    \return Nodal coordinates [#nnode, #ndim].
    */
    xt::xtensor<double, 2> coor() const;

    /**
    \return Connectivity [#nelem, #nne].
    */
    xt::xtensor<size_t, 2> conn() const;

    /**
    \return DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    \return Node-map per sub-mesh.
    */
    std::vector<xt::xtensor<size_t, 1>> nodemap() const;

    /**
    \return Element-map per sub-mesh.
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
    \return Number of sub meshes.
    */
    size_t nmesh() const;

    /**
    \return Number of elements.
    */
    size_t nelem() const;

    /**
    \return Number of nodes.
    */
    size_t nnode() const;

    /**
    \return Number of nodes-per-element.
    */
    size_t nne() const;

    /**
    \return Number of dimensions.
    */
    size_t ndim() const;

    /**
    \return Nodal coordinates [#nnode, #ndim].
    */
    xt::xtensor<double, 2> coor() const;

    /**
    \return Connectivity [#nelem, #nne].
    */
    xt::xtensor<size_t, 2> conn() const;

    /**
    \return DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    \return Node-map per sub-mesh.
    */
    std::vector<xt::xtensor<size_t, 1>> nodemap() const;

    /**
    \return Element-map per sub-mesh.
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

private:
    xt::xtensor<double, 2> m_coor;
    xt::xtensor<size_t, 2> m_conn;
    std::vector<xt::xtensor<size_t, 1>> m_map;
    std::vector<size_t> m_nel; ///< Number of elements per sub-mesh.
    std::vector<size_t> m_el_offset;
    double m_rtol;
    double m_atol;
};

/**
\rst
Renumber indices to lowest possible index. For example:

.. math::

    \begin{bmatrix}
        0 & 1 \\
        5 & 4
    \end{bmatrix}

is renumbered to

.. math::

    \begin{bmatrix}
        0 & 1 \\
        3 & 2
    \end{bmatrix}

Or, in pseudo-code, the result of this function is that:

.. code-block:: python

    dofs = renumber(dofs)

    sort(unique(dofs[:])) == range(max(dofs+1))

.. tip::

    One can use the wrapper function :cpp:func:`GooseFEM::Mesh::renumber`.
    This class gives more advanced features.
\endrst
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

\rst
.. math::

    \begin{bmatrix}
        0 & 1 \\
        2 & 3 \\
        4 & 5
    \end{bmatrix}
\endrst

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

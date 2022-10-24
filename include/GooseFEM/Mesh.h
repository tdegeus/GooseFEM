/**
 * Generic mesh operations.
 *
 * @file Mesh.h
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#ifndef GOOSEFEM_MESH_H
#define GOOSEFEM_MESH_H

#include "ElementQuad4.h"
#include "MatrixDiagonal.h"
#include "Vector.h"
#include "assertions.h"
#include "config.h"

namespace GooseFEM {

/**
 * Generic mesh operations, and simple mesh definitions.
 */
namespace Mesh {

template <class D>
inline std::vector<std::vector<size_t>> nodaltyings(const D& dofs);

/**
 * Enumerator for element-types
 */
enum class ElementType {
    Unknown, ///< Unknown element-type
    Quad4, ///< Quadrilateral: 4-noded element in 2-d
    Hex8, ///< Hexahedron: 8-noded element in 3-d
    Tri3 ///< Triangle: 3-noded element in 2-d
};

/**
 * Extract the element type based on the connectivity.
 *
 * @param coor Nodal coordinates [nnode, ndim].
 * @param conn Connectivity [nelem, nne].
 * @return ElementType().
 */
template <class S, class T>
inline ElementType defaultElementType(const S& coor, const T& conn)
{
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);

    if (coor.shape(1) == 2ul && conn.shape(1) == 3ul) {
        return ElementType::Tri3;
    }
    if (coor.shape(1) == 2ul && conn.shape(1) == 4ul) {
        return ElementType::Quad4;
    }
    if (coor.shape(1) == 3ul && conn.shape(1) == 8ul) {
        return ElementType::Hex8;
    }

    throw std::runtime_error("Element-type not implemented");
}

namespace detail {

template <class T, class R>
inline T renum(const T& arg, const R& mapping)
{
    T ret = T::from_shape(arg.shape());

    auto jt = ret.begin();

    for (auto it = arg.begin(); it != arg.end(); ++it, ++jt) {
        *jt = mapping(*it);
    }

    return ret;
}

} // namespace detail

/**
 * List with DOF-numbers in sequential order.
 * The output is a sequential list of DOF-numbers for each vector-component of each node.
 * For example for 3 nodes in 2 dimensions the output is
 *
 * \f$ \begin{bmatrix} 0 & 1 \\ 2 & 3 \\ 4 & 5 \end{bmatrix} \f$
 *
 * @param nnode Number of nodes.
 * @param ndim Number of dimensions.
 * @return DOF-numbers.
 */
inline array_type::tensor<size_t, 2> dofs(size_t nnode, size_t ndim)
{
    return xt::reshape_view(xt::arange<size_t>(nnode * ndim), {nnode, ndim});
}

/**
 * Renumber indices to lowest possible index. For example:
 *
 * \f$ \begin{bmatrix} 0 & 1 \\ 5 & 4 \end{bmatrix} \f$
 *
 * is renumbered to
 *
 * \f$ \begin{bmatrix} 0 & 1 \\ 3 & 2 \end{bmatrix} \f$
 *
 * Or, in pseudo-code, the result of this function is that:
 *
 *      dofs = renumber(dofs)
 *      sort(unique(dofs[:])) == range(max(dofs+1))
 *
 * \note One can use the wrapper function renumber(). This class gives more advanced features.
 */
class Renumber {
public:
    Renumber() = default;

    /**
     * @param dofs DOF-numbers.
     */
    template <class T>
    Renumber(const T& dofs)
    {
        size_t n = xt::amax(dofs)() + 1;
        size_t i = 0;

        array_type::tensor<size_t, 1> unique = xt::unique(dofs);

        m_renum = xt::empty<size_t>({n});

        for (auto& j : unique) {
            m_renum(j) = i;
            ++i;
        }
    }

    /**
     * Apply renumbering to other set.
     *
     * @param list List of (DOF-)numbers.
     * @return Renumbered list of (DOF-)numbers.
     */
    template <class T>
    T apply(const T& list) const
    {
        return detail::renum(list, m_renum);
    }

    /**
     * Get the list needed to renumber, e.g.:
     *
     *      dofs_renumbered(i, j) = index(dofs(i, j))
     *
     * @return Renumber-index.
     */
    const array_type::tensor<size_t, 1>& index() const
    {
        return m_renum;
    }

private:
    array_type::tensor<size_t, 1> m_renum;
};

/**
 * Renumber to lowest possible index (see GooseFEM::Mesh::Renumber).
 *
 * @param dofs DOF-numbers [nnode, ndim].
 * @return Renumbered DOF-numbers.
 */
template <class T>
inline T renumber(const T& dofs)
{
    return Renumber(dofs).apply(dofs);
}

/**
 * CRTP base class for regular meshes.
 */
template <class D>
class RegularBase {
public:
    /**
     * Underlying type.
     */
    using derived_type = D;

    /**
     * Number of elements.
     * @return unsigned int
     */
    auto nelem() const
    {
        return derived_cast().m_nelem;
    }

    /**
     * Number of nodes.
     * @return unsigned int
     */
    auto nnode() const
    {
        return derived_cast().m_nnode;
    }

    /**
     * Number of nodes-per-element == 4.
     * @return unsigned int
     */
    auto nne() const
    {
        return derived_cast().m_nne;
    }

    /**
     * Number of dimensions == 2.
     * @return unsigned int
     */
    auto ndim() const
    {
        return derived_cast().m_ndim;
    }

    /**
     * Number of elements in x-direction == width of the mesh in units of #h.
     * @return unsigned int
     */
    auto nelx() const
    {
        return derived_cast().nelx_impl();
    }

    /**
     * Number of elements in y-direction == height of the mesh, in units of #h,
     * @return unsigned int
     */
    auto nely() const
    {
        return derived_cast().nely_impl();
    }

    /**
     * Linear edge size of one 'block'.
     * @return double
     */
    auto h() const
    {
        return derived_cast().m_h;
    }

    /**
     * The ElementType().
     * @return element type
     */
    auto getElementType() const
    {
        return derived_cast().getElementType_impl();
    }

    /**
     * Nodal coordinates [#nnode, #ndim].
     * @return coordinates per node
     */
    auto coor() const
    {
        return derived_cast().coor_impl();
    }

    /**
     * Connectivity [#nelem, #nne].
     * @return nodes per element
     */
    auto conn() const
    {
        return derived_cast().conn_impl();
    }

    /**
     * DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
     * @return DOFs per node
     */
    auto dofs() const
    {
        return GooseFEM::Mesh::dofs(this->nnode(), this->ndim());
    }

    /**
     * DOF-numbers for the case that the periodicity if fully eliminated.
     * @return DOF numbers for each node [#nnode, #ndim].
     */
    auto dofsPeriodic() const
    {
        array_type::tensor<size_t, 2> ret = this->dofs();
        array_type::tensor<size_t, 2> nodePer = this->nodesPeriodic();
        array_type::tensor<size_t, 1> independent = xt::view(nodePer, xt::all(), 0);
        array_type::tensor<size_t, 1> dependent = xt::view(nodePer, xt::all(), 1);

        for (size_t j = 0; j < this->ndim(); ++j) {
            xt::view(ret, xt::keep(dependent), j) = xt::view(ret, xt::keep(independent), j);
        }

        return GooseFEM::Mesh::renumber(ret);
    }

    /**
     * Periodic node pairs, in two columns: (independent, dependent).
     * @return [ntyings, #ndim].
     */
    auto nodesPeriodic() const
    {
        return derived_cast().nodesPeriodic_impl();
    }

    /**
     * Reference node to use for periodicity, because all corners are tied to it.
     * @return Node number.
     */
    auto nodesOrigin() const
    {
        return derived_cast().nodesOrigin_impl();
    }

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }
};

/**
 * CRTP base class for regular meshes in 2d.
 */
template <class D>
class RegularBase2d : public RegularBase<D> {
public:
    /**
     * Underlying type.
     */
    using derived_type = D;

    /**
     * Nodes along the bottom edge (y = 0), in order of increasing x.
     * @return List of node numbers.
     */
    auto nodesBottomEdge() const
    {
        return derived_cast().nodesBottomEdge_impl();
    }

    /**
     * Nodes along the top edge (y = #nely * #h), in order of increasing x.
     * @return List of node numbers.
     */
    auto nodesTopEdge() const
    {
        return derived_cast().nodesTopEdge_impl();
    }

    /**
     * Nodes along the left edge (x = 0), in order of increasing y.
     * @return List of node numbers.
     */
    auto nodesLeftEdge() const
    {
        return derived_cast().nodesLeftEdge_impl();
    }

    /**
     * Nodes along the right edge (x = #nelx * #h), in order of increasing y.
     * @return List of node numbers.
     */
    auto nodesRightEdge() const
    {
        return derived_cast().nodesRightEdge_impl();
    }

    /**
     * Nodes along the bottom edge (y = 0), without the corners (at x = 0 and x = #nelx * #h).
     * Same as: nodesBottomEdge()[1: -1].
     * @return List of node numbers.
     */
    auto nodesBottomOpenEdge() const
    {
        return derived_cast().nodesBottomOpenEdge_impl();
    }

    /**
     * Nodes along the top edge (y = #nely * #h), without the corners (at x = 0 and x = #nelx * #h).
     * Same as: nodesTopEdge()[1: -1].
     * @return List of node numbers.
     */
    auto nodesTopOpenEdge() const
    {
        return derived_cast().nodesTopOpenEdge_impl();
    }

    /**
     * Nodes along the left edge (x = 0), without the corners (at y = 0 and y = #nely * #h).
     * Same as: nodesLeftEdge()[1: -1].
     * @return List of node numbers.
     */
    auto nodesLeftOpenEdge() const
    {
        return derived_cast().nodesLeftOpenEdge_impl();
    }

    /**
     * Nodes along the right edge (x = #nelx * #h), without the corners (at y = 0 and y = #nely *
     * #h). Same as: nodesRightEdge()[1: -1].
     * @return List of node numbers.
     */
    auto nodesRightOpenEdge() const
    {
        return derived_cast().nodesRightOpenEdge_impl();
    }

    /**
     * The bottom-left corner node (at x = 0, y = 0).
     * Same as nodesBottomEdge()[0] and nodesLeftEdge()[0].
     * @return Node number.
     */
    auto nodesBottomLeftCorner() const
    {
        return derived_cast().nodesBottomLeftCorner_impl();
    }

    /**
     * The bottom-right corner node (at x = #nelx * #h, y = 0).
     * Same as nodesBottomEdge()[-1] and nodesRightEdge()[0].
     * @return Node number.
     */
    auto nodesBottomRightCorner() const
    {
        return derived_cast().nodesBottomRightCorner_impl();
    }

    /**
     * The top-left corner node (at x = 0, y = #nely * #h).
     * Same as nodesTopEdge()[0] and nodesRightEdge()[-1].
     * @return Node number.
     */
    auto nodesTopLeftCorner() const
    {
        return derived_cast().nodesTopLeftCorner_impl();
    }

    /**
     * The top-right corner node (at x = #nelx * #h, y = #nely * #h).
     * Same as nodesTopEdge()[-1] and nodesRightEdge()[-1].
     * @return Node number.
     */
    auto nodesTopRightCorner() const
    {
        return derived_cast().nodesTopRightCorner_impl();
    }

    /**
     * Alias of nodesBottomLeftCorner().
     * @return Node number.
     */
    auto nodesLeftBottomCorner() const
    {
        return derived_cast().nodesBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesTopLeftCorner().
     * @return Node number.
     */
    auto nodesLeftTopCorner() const
    {
        return derived_cast().nodesTopLeftCorner_impl();
    }

    /**
     * Alias of nodesBottomRightCorner().
     * @return Node number.
     */
    auto nodesRightBottomCorner() const
    {
        return derived_cast().nodesBottomRightCorner_impl();
    }

    /**
     * Alias of nodesTopRightCorner().
     * @return Node number.
     */
    auto nodesRightTopCorner() const
    {
        return derived_cast().nodesTopRightCorner_impl();
    }

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

    friend class RegularBase<D>;

    array_type::tensor<size_t, 2> nodesPeriodic_impl() const
    {
        array_type::tensor<size_t, 1> bot = derived_cast().nodesBottomOpenEdge_impl();
        array_type::tensor<size_t, 1> top = derived_cast().nodesTopOpenEdge_impl();
        array_type::tensor<size_t, 1> lft = derived_cast().nodesLeftOpenEdge_impl();
        array_type::tensor<size_t, 1> rgt = derived_cast().nodesRightOpenEdge_impl();
        std::array<size_t, 2> shape = {bot.size() + lft.size() + size_t(3), size_t(2)};
        auto ret = array_type::tensor<size_t, 2>::from_shape(shape);

        ret(0, 0) = derived_cast().nodesBottomLeftCorner_impl();
        ret(0, 1) = derived_cast().nodesBottomRightCorner_impl();

        ret(1, 0) = derived_cast().nodesBottomLeftCorner_impl();
        ret(1, 1) = derived_cast().nodesTopRightCorner_impl();

        ret(2, 0) = derived_cast().nodesBottomLeftCorner_impl();
        ret(2, 1) = derived_cast().nodesTopLeftCorner_impl();

        size_t i = 3;

        xt::view(ret, xt::range(i, i + bot.size()), 0) = bot;
        xt::view(ret, xt::range(i, i + bot.size()), 1) = top;

        i += bot.size();

        xt::view(ret, xt::range(i, i + lft.size()), 0) = lft;
        xt::view(ret, xt::range(i, i + lft.size()), 1) = rgt;

        return ret;
    }

    auto nodesOrigin_impl() const
    {
        return derived_cast().nodesBottomLeftCorner_impl();
    }
};

/**
 * CRTP base class for regular meshes in 3d.
 */
template <class D>
class RegularBase3d : public RegularBase<D> {
public:
    /**
     * Underlying type.
     */
    using derived_type = D;

    /**
     * Number of elements in y-direction == height of the mesh, in units of #h,
     * @return unsigned int
     */
    auto nelz() const
    {
        return derived_cast().nelz_impl();
    }

    /**
     * Nodes along the bottom face (y = 0).
     * @return List of node numbers.
     */
    auto nodesBottom() const
    {
        return derived_cast().nodesBottom_impl();
    }

    /**
     * Nodes along the top face (y = #nely * #h).
     * @return List of node numbers.
     */
    auto nodesTop() const
    {
        return derived_cast().nodesTop_impl();
    }

    /**
     * Nodes along the left face (x = 0).
     * @return List of node numbers.
     */
    auto nodesLeft() const
    {
        return derived_cast().nodesLeft_impl();
    }

    /**
     * Nodes along the right face (x = #nelx * #h).
     * @return List of node numbers.
     */
    auto nodesRight() const
    {
        return derived_cast().nodesRight_impl();
    }

    /**
     * Nodes along the front face (z = 0).
     * @return List of node numbers.
     */
    auto nodesFront() const
    {
        return derived_cast().nodesFront_impl();
    }

    /**
     * Nodes along the back face (z = #nelz * #h).
     * @return List of node numbers.
     */
    auto nodesBack() const
    {
        return derived_cast().nodesBack_impl();
    }

    /**
     * Nodes along the edge at the intersection of the front and bottom faces
     * (z = 0 and y = 0).
     * @return List of node numbers.
     */
    auto nodesFrontBottomEdge() const
    {
        return derived_cast().nodesFrontBottomEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the front and top faces
     * (z = 0 and y = #nely * #h).
     * @return List of node numbers.
     */
    auto nodesFrontTopEdge() const
    {
        return derived_cast().nodesFrontTopEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the front and left faces
     * (z = 0 and x = 0).
     * @return List of node numbers.
     */
    auto nodesFrontLeftEdge() const
    {
        return derived_cast().nodesFrontLeftEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the front and right faces
     * (z = 0 and x = #nelx * #h).
     * @return List of node numbers.
     */
    auto nodesFrontRightEdge() const
    {
        return derived_cast().nodesFrontRightEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the back and bottom faces
     * (z = #nelz * #h and y = #nely * #h).
     * @return List of node numbers.
     */
    auto nodesBackBottomEdge() const
    {
        return derived_cast().nodesBackBottomEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the back and top faces
     * (z = #nelz * #h and x = 0).
     * @return List of node numbers.
     */
    auto nodesBackTopEdge() const
    {
        return derived_cast().nodesBackTopEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the back and left faces
     * (z = #nelz * #h and x = #nelx * #h).
     * @return List of node numbers.
     */
    auto nodesBackLeftEdge() const
    {
        return derived_cast().nodesBackLeftEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the back and right faces
     * (? = #nelz * #h and ? = ?).
     * @return List of node numbers.
     */
    auto nodesBackRightEdge() const
    {
        return derived_cast().nodesBackRightEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the bottom and left faces
     * (y = 0 and x = 0).
     * @return List of node numbers.
     */
    auto nodesBottomLeftEdge() const
    {
        return derived_cast().nodesBottomLeftEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the bottom and right faces
     * (y = 0 and x = #nelx * #h).
     * @return List of node numbers.
     */
    auto nodesBottomRightEdge() const
    {
        return derived_cast().nodesBottomRightEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the top and left faces
     * (y = 0 and x = #nelx * #h).
     * @return List of node numbers.
     */
    auto nodesTopLeftEdge() const
    {
        return derived_cast().nodesTopLeftEdge_impl();
    }

    /**
     * Nodes along the edge at the intersection of the top and right faces
     * (y = #nely * #h and x = #nelx * #h).
     * @return List of node numbers.
     */
    auto nodesTopRightEdge() const
    {
        return derived_cast().nodesTopRightEdge_impl();
    }

    /**
     * Alias of nodesFrontBottomEdge()
     * @return List of node numbers.
     */
    auto nodesBottomFrontEdge() const
    {
        return derived_cast().nodesFrontBottomEdge_impl();
    }

    /**
     * Alias of nodesBackBottomEdge()
     * @return List of node numbers.
     */
    auto nodesBottomBackEdge() const
    {
        return derived_cast().nodesBackBottomEdge_impl();
    }

    /**
     * Alias of nodesFrontTopEdge()
     * @return List of node numbers.
     */
    auto nodesTopFrontEdge() const
    {
        return derived_cast().nodesFrontTopEdge_impl();
    }

    /**
     * Alias of nodesBackTopEdge()
     * @return List of node numbers.
     */
    auto nodesTopBackEdge() const
    {
        return derived_cast().nodesBackTopEdge_impl();
    }

    /**
     * Alias of nodesBottomLeftEdge()
     * @return List of node numbers.
     */
    auto nodesLeftBottomEdge() const
    {
        return derived_cast().nodesBottomLeftEdge_impl();
    }

    /**
     * Alias of nodesFrontLeftEdge()
     * @return List of node numbers.
     */
    auto nodesLeftFrontEdge() const
    {
        return derived_cast().nodesFrontLeftEdge_impl();
    }

    /**
     * Alias of nodesBackLeftEdge()
     * @return List of node numbers.
     */
    auto nodesLeftBackEdge() const
    {
        return derived_cast().nodesBackLeftEdge_impl();
    }

    /**
     * Alias of nodesTopLeftEdge()
     * @return List of node numbers.
     */
    auto nodesLeftTopEdge() const
    {
        return derived_cast().nodesTopLeftEdge_impl();
    }

    /**
     * Alias of nodesBottomRightEdge()
     * @return List of node numbers.
     */
    auto nodesRightBottomEdge() const
    {
        return derived_cast().nodesBottomRightEdge_impl();
    }

    /**
     * Alias of nodesTopRightEdge()
     * @return List of node numbers.
     */
    auto nodesRightTopEdge() const
    {
        return derived_cast().nodesTopRightEdge_impl();
    }

    /**
     * Alias of nodesFrontRightEdge()
     * @return List of node numbers.
     */
    auto nodesRightFrontEdge() const
    {
        return derived_cast().nodesFrontRightEdge_impl();
    }

    /**
     * Alias of nodesBackRightEdge()
     * @return List of node numbers.
     */
    auto nodesRightBackEdge() const
    {
        return derived_cast().nodesBackRightEdge_impl();
    }

    /**
     * Nodes along the front face excluding edges.
     * Same as different between nodesFront() and
     * [nodesFrontBottomEdge(), nodesFrontTopEdge(), nodesFrontLeftEdge(), nodesFrontRightEdge()]
     * @return list of node numbers.
     */
    auto nodesFrontFace() const
    {
        return derived_cast().nodesFrontFace_impl();
    }

    /**
     * Nodes along the back face excluding edges.
     * Same as different between nodesBack() and
     * [nodesBackBottomEdge(), nodesBackTopEdge(), nodesBackLeftEdge(), nodesBackRightEdge()]
     * @return list of node numbers.
     */
    auto nodesBackFace() const
    {
        return derived_cast().nodesBackFace_impl();
    }

    /**
     * Nodes along the left face excluding edges.
     * Same as different between nodesLeft() and
     * [nodesFrontLeftEdge(), nodesBackLeftEdge(), nodesBottomLeftEdge(), nodesTopLeftEdge()]
     * @return list of node numbers.
     */
    auto nodesLeftFace() const
    {
        return derived_cast().nodesLeftFace_impl();
    }

    /**
     * Nodes along the right face excluding edges.
     * Same as different between nodesRight() and
     * [nodesFrontRightEdge(), nodesBackRightEdge(), nodesBottomRightEdge(), nodesTopRightEdge()]
     * @return list of node numbers.
     */
    auto nodesRightFace() const
    {
        return derived_cast().nodesRightFace_impl();
    }

    /**
     * Nodes along the bottom face excluding edges.
     * Same as different between nodesBottom() and
     * [nodesBackBottomEdge(), nodesBackTopEdge(), nodesBackLeftEdge(), nodesBackRightEdge()]
     * @return list of node numbers.
     */
    auto nodesBottomFace() const
    {
        return derived_cast().nodesBottomFace_impl();
    }

    /**
     * Nodes along the top face excluding edges.
     * Same as different between nodesTop() and
     * [nodesFrontBottomEdge(), nodesFrontTopEdge(), nodesFrontLeftEdge(), nodesFrontRightEdge()]
     * @return list of node numbers.
     */
    auto nodesTopFace() const
    {
        return derived_cast().nodesTopFace_impl();
    }

    /**
     * Same as nodesFrontBottomEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesFrontBottomOpenEdge() const
    {
        return derived_cast().nodesFrontBottomOpenEdge_impl();
    }

    /**
     * Same as nodesFrontTopEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesFrontTopOpenEdge() const
    {
        return derived_cast().nodesFrontTopOpenEdge_impl();
    }

    /**
     * Same as nodesFrontLeftEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesFrontLeftOpenEdge() const
    {
        return derived_cast().nodesFrontLeftOpenEdge_impl();
    }

    /**
     * Same as nodesFrontRightEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesFrontRightOpenEdge() const
    {
        return derived_cast().nodesFrontRightOpenEdge_impl();
    }

    /**
     * Same as nodesBackBottomEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesBackBottomOpenEdge() const
    {
        return derived_cast().nodesBackBottomOpenEdge_impl();
    }

    /**
     * Same as nodesBackTopEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesBackTopOpenEdge() const
    {
        return derived_cast().nodesBackTopOpenEdge_impl();
    }

    /**
     * Same as nodesBackLeftEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesBackLeftOpenEdge() const
    {
        return derived_cast().nodesBackLeftOpenEdge_impl();
    }

    /**
     * Same as nodesBackRightEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesBackRightOpenEdge() const
    {
        return derived_cast().nodesBackRightOpenEdge_impl();
    }

    /**
     * Same as nodesBottomLeftEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesBottomLeftOpenEdge() const
    {
        return derived_cast().nodesBottomLeftOpenEdge_impl();
    }

    /**
     * Same as nodesBottomRightEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesBottomRightOpenEdge() const
    {
        return derived_cast().nodesBottomRightOpenEdge_impl();
    }

    /**
     * Same as nodesTopLeftEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesTopLeftOpenEdge() const
    {
        return derived_cast().nodesTopLeftOpenEdge_impl();
    }

    /**
     * Same as nodesTopRightEdge() but without corners.
     * @return List of node numbers.
     */
    auto nodesTopRightOpenEdge() const
    {
        return derived_cast().nodesTopRightOpenEdge_impl();
    }

    /**
     * Alias of nodesFrontBottomOpenEdge().
     * @return List of node numbers.
     */
    auto nodesBottomFrontOpenEdge() const
    {
        return derived_cast().nodesFrontBottomOpenEdge_impl();
    }

    /**
     * Alias of nodesBackBottomOpenEdge().
     * @return List of node numbers.
     */
    auto nodesBottomBackOpenEdge() const
    {
        return derived_cast().nodesBackBottomOpenEdge_impl();
    }

    /**
     * Alias of nodesFrontTopOpenEdge().
     * @return List of node numbers.
     */
    auto nodesTopFrontOpenEdge() const
    {
        return derived_cast().nodesFrontTopOpenEdge_impl();
    }

    /**
     * Alias of nodesBackTopOpenEdge().
     * @return List of node numbers.
     */
    auto nodesTopBackOpenEdge() const
    {
        return derived_cast().nodesBackTopOpenEdge_impl();
    }

    /**
     * Alias of nodesBottomLeftOpenEdge().
     * @return List of node numbers.
     */
    auto nodesLeftBottomOpenEdge() const
    {
        return derived_cast().nodesBottomLeftOpenEdge_impl();
    }

    /**
     * Alias of nodesFrontLeftOpenEdge().
     * @return List of node numbers.
     */
    auto nodesLeftFrontOpenEdge() const
    {
        return derived_cast().nodesFrontLeftOpenEdge_impl();
    }

    /**
     * Alias of nodesBackLeftOpenEdge().
     * @return List of node numbers.
     */
    auto nodesLeftBackOpenEdge() const
    {
        return derived_cast().nodesBackLeftOpenEdge_impl();
    }

    /**
     * Alias of nodesTopLeftOpenEdge().
     * @return List of node numbers.
     */
    auto nodesLeftTopOpenEdge() const
    {
        return derived_cast().nodesTopLeftOpenEdge_impl();
    }

    /**
     * Alias of nodesBottomRightOpenEdge().
     * @return List of node numbers.
     */
    auto nodesRightBottomOpenEdge() const
    {
        return derived_cast().nodesBottomRightOpenEdge_impl();
    }

    /**
     * Alias of nodesTopRightOpenEdge().
     * @return List of node numbers.
     */
    auto nodesRightTopOpenEdge() const
    {
        return derived_cast().nodesTopRightOpenEdge_impl();
    }

    /**
     * Alias of nodesFrontRightOpenEdge().
     * @return List of node numbers.
     */
    auto nodesRightFrontOpenEdge() const
    {
        return derived_cast().nodesFrontRightOpenEdge_impl();
    }

    /**
     * Alias of nodesBackRightOpenEdge().
     * @return List of node numbers.
     */
    auto nodesRightBackOpenEdge() const
    {
        return derived_cast().nodesBackRightOpenEdge_impl();
    }

    /**
     * Front-Bottom-Left corner node.
     * @return Node number.
     */
    auto nodesFrontBottomLeftCorner() const
    {
        return derived_cast().nodesFrontBottomLeftCorner_impl();
    }

    /**
     * Front-Bottom-Right corner node.
     * @return Node number.
     */
    auto nodesFrontBottomRightCorner() const
    {
        return derived_cast().nodesFrontBottomRightCorner_impl();
    }

    /**
     * Front-Top-Left corner node.
     * @return Node number.
     */
    auto nodesFrontTopLeftCorner() const
    {
        return derived_cast().nodesFrontTopLeftCorner_impl();
    }

    /**
     * Front-Top-Right corner node.
     * @return Node number.
     */
    auto nodesFrontTopRightCorner() const
    {
        return derived_cast().nodesFrontTopRightCorner_impl();
    }

    /**
     * Back-Bottom-Left corner node.
     * @return Node number.
     */
    auto nodesBackBottomLeftCorner() const
    {
        return derived_cast().nodesBackBottomLeftCorner_impl();
    }

    /**
     * Back-Bottom-Right corner node.
     * @return Node number.
     */
    auto nodesBackBottomRightCorner() const
    {
        return derived_cast().nodesBackBottomRightCorner_impl();
    }

    /**
     * Back-Top-Left corner node.
     * @return Node number.
     */
    auto nodesBackTopLeftCorner() const
    {
        return derived_cast().nodesBackTopLeftCorner_impl();
    }

    /**
     * Back-Top-Right corner node.
     * @return Node number.
     */
    auto nodesBackTopRightCorner() const
    {
        return derived_cast().nodesBackTopRightCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomLeftCorner().
     * @return Node number.
     */
    auto nodesFrontLeftBottomCorner() const
    {
        return derived_cast().nodesFrontBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomLeftCorner().
     * @return Node number.
     */
    auto nodesBottomFrontLeftCorner() const
    {
        return derived_cast().nodesFrontBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomLeftCorner().
     * @return Node number.
     */
    auto nodesBottomLeftFrontCorner() const
    {
        return derived_cast().nodesFrontBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomLeftCorner().
     * @return Node number.
     */
    auto nodesLeftFrontBottomCorner() const
    {
        return derived_cast().nodesFrontBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomLeftCorner().
     * @return Node number.
     */
    auto nodesLeftBottomFrontCorner() const
    {
        return derived_cast().nodesFrontBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomRightCorner().
     * @return Node number.
     */
    auto nodesFrontRightBottomCorner() const
    {
        return derived_cast().nodesFrontBottomRightCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomRightCorner().
     * @return Node number.
     */
    auto nodesBottomFrontRightCorner() const
    {
        return derived_cast().nodesFrontBottomRightCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomRightCorner().
     * @return Node number.
     */
    auto nodesBottomRightFrontCorner() const
    {
        return derived_cast().nodesFrontBottomRightCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomRightCorner().
     * @return Node number.
     */
    auto nodesRightFrontBottomCorner() const
    {
        return derived_cast().nodesFrontBottomRightCorner_impl();
    }

    /**
     * Alias of nodesFrontBottomRightCorner().
     * @return Node number.
     */
    auto nodesRightBottomFrontCorner() const
    {
        return derived_cast().nodesFrontBottomRightCorner_impl();
    }

    /**
     * Alias of nodesFrontTopLeftCorner().
     * @return Node number.
     */
    auto nodesFrontLeftTopCorner() const
    {
        return derived_cast().nodesFrontTopLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontTopLeftCorner().
     * @return Node number.
     */
    auto nodesTopFrontLeftCorner() const
    {
        return derived_cast().nodesFrontTopLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontTopLeftCorner().
     * @return Node number.
     */
    auto nodesTopLeftFrontCorner() const
    {
        return derived_cast().nodesFrontTopLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontTopLeftCorner().
     * @return Node number.
     */
    auto nodesLeftFrontTopCorner() const
    {
        return derived_cast().nodesFrontTopLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontTopLeftCorner().
     * @return Node number.
     */
    auto nodesLeftTopFrontCorner() const
    {
        return derived_cast().nodesFrontTopLeftCorner_impl();
    }

    /**
     * Alias of nodesFrontTopRightCorner().
     * @return Node number.
     */
    auto nodesFrontRightTopCorner() const
    {
        return derived_cast().nodesFrontTopRightCorner_impl();
    }

    /**
     * Alias of nodesFrontTopRightCorner().
     * @return Node number.
     */
    auto nodesTopFrontRightCorner() const
    {
        return derived_cast().nodesFrontTopRightCorner_impl();
    }

    /**
     * Alias of nodesFrontTopRightCorner().
     * @return Node number.
     */
    auto nodesTopRightFrontCorner() const
    {
        return derived_cast().nodesFrontTopRightCorner_impl();
    }

    /**
     * Alias of nodesFrontTopRightCorner().
     * @return Node number.
     */
    auto nodesRightFrontTopCorner() const
    {
        return derived_cast().nodesFrontTopRightCorner_impl();
    }

    /**
     * Alias of nodesFrontTopRightCorner().
     * @return Node number.
     */
    auto nodesRightTopFrontCorner() const
    {
        return derived_cast().nodesFrontTopRightCorner_impl();
    }

    /**
     * Alias of nodesBackBottomLeftCorner().
     * @return Node number.
     */
    auto nodesBackLeftBottomCorner() const
    {
        return derived_cast().nodesBackBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesBackBottomLeftCorner().
     * @return Node number.
     */
    auto nodesBottomBackLeftCorner() const
    {
        return derived_cast().nodesBackBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesBackBottomLeftCorner().
     * @return Node number.
     */
    auto nodesBottomLeftBackCorner() const
    {
        return derived_cast().nodesBackBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesBackBottomLeftCorner().
     * @return Node number.
     */
    auto nodesLeftBackBottomCorner() const
    {
        return derived_cast().nodesBackBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesBackBottomLeftCorner().
     * @return Node number.
     */
    auto nodesLeftBottomBackCorner() const
    {
        return derived_cast().nodesBackBottomLeftCorner_impl();
    }

    /**
     * Alias of nodesBackBottomRightCorner().
     * @return Node number.
     */
    auto nodesBackRightBottomCorner() const
    {
        return derived_cast().nodesBackBottomRightCorner_impl();
    }

    /**
     * Alias of nodesBackBottomRightCorner().
     * @return Node number.
     */
    auto nodesBottomBackRightCorner() const
    {
        return derived_cast().nodesBackBottomRightCorner_impl();
    }

    /**
     * Alias of nodesBackBottomRightCorner().
     * @return Node number.
     */
    auto nodesBottomRightBackCorner() const
    {
        return derived_cast().nodesBackBottomRightCorner_impl();
    }

    /**
     * Alias of nodesBackBottomRightCorner().
     * @return Node number.
     */
    auto nodesRightBackBottomCorner() const
    {
        return derived_cast().nodesBackBottomRightCorner_impl();
    }

    /**
     * Alias of nodesBackBottomRightCorner().
     * @return Node number.
     */
    auto nodesRightBottomBackCorner() const
    {
        return derived_cast().nodesBackBottomRightCorner_impl();
    }

    /**
     * Alias of nodesBackTopLeftCorner().
     * @return Node number.
     */
    auto nodesBackLeftTopCorner() const
    {
        return derived_cast().nodesBackTopLeftCorner_impl();
    }

    /**
     * Alias of nodesBackTopLeftCorner().
     * @return Node number.
     */
    auto nodesTopBackLeftCorner() const
    {
        return derived_cast().nodesBackTopLeftCorner_impl();
    }

    /**
     * Alias of nodesBackTopLeftCorner().
     * @return Node number.
     */
    auto nodesTopLeftBackCorner() const
    {
        return derived_cast().nodesBackTopLeftCorner_impl();
    }

    /**
     * Alias of nodesBackTopLeftCorner().
     * @return Node number.
     */
    auto nodesLeftBackTopCorner() const
    {
        return derived_cast().nodesBackTopLeftCorner_impl();
    }

    /**
     * Alias of nodesBackTopLeftCorner().
     * @return Node number.
     */
    auto nodesLeftTopBackCorner() const
    {
        return derived_cast().nodesBackTopLeftCorner_impl();
    }

    /**
     * Alias of nodesBackTopRightCorner().
     * @return Node number.
     */
    auto nodesBackRightTopCorner() const
    {
        return derived_cast().nodesBackTopRightCorner_impl();
    }

    /**
     * Alias of nodesBackTopRightCorner().
     * @return Node number.
     */
    auto nodesTopBackRightCorner() const
    {
        return derived_cast().nodesBackTopRightCorner_impl();
    }

    /**
     * Alias of nodesBackTopRightCorner().
     * @return Node number.
     */
    auto nodesTopRightBackCorner() const
    {
        return derived_cast().nodesBackTopRightCorner_impl();
    }

    /**
     * Alias of nodesBackTopRightCorner().
     * @return Node number.
     */
    auto nodesRightBackTopCorner() const
    {
        return derived_cast().nodesBackTopRightCorner_impl();
    }

    /**
     * Alias of nodesBackTopRightCorner().
     * @return Node number.
     */
    auto nodesRightTopBackCorner() const
    {
        return derived_cast().nodesBackTopRightCorner_impl();
    }

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

    friend class RegularBase<D>;

    array_type::tensor<size_t, 2> nodesPeriodic_impl() const
    {
        array_type::tensor<size_t, 1> fro = derived_cast().nodesFrontFace_impl();
        array_type::tensor<size_t, 1> bck = derived_cast().nodesBackFace_impl();
        array_type::tensor<size_t, 1> lft = derived_cast().nodesLeftFace_impl();
        array_type::tensor<size_t, 1> rgt = derived_cast().nodesRightFace_impl();
        array_type::tensor<size_t, 1> bot = derived_cast().nodesBottomFace_impl();
        array_type::tensor<size_t, 1> top = derived_cast().nodesTopFace_impl();

        array_type::tensor<size_t, 1> froBot = derived_cast().nodesFrontBottomOpenEdge_impl();
        array_type::tensor<size_t, 1> froTop = derived_cast().nodesFrontTopOpenEdge_impl();
        array_type::tensor<size_t, 1> froLft = derived_cast().nodesFrontLeftOpenEdge_impl();
        array_type::tensor<size_t, 1> froRgt = derived_cast().nodesFrontRightOpenEdge_impl();
        array_type::tensor<size_t, 1> bckBot = derived_cast().nodesBackBottomOpenEdge_impl();
        array_type::tensor<size_t, 1> bckTop = derived_cast().nodesBackTopOpenEdge_impl();
        array_type::tensor<size_t, 1> bckLft = derived_cast().nodesBackLeftOpenEdge_impl();
        array_type::tensor<size_t, 1> bckRgt = derived_cast().nodesBackRightOpenEdge_impl();
        array_type::tensor<size_t, 1> botLft = derived_cast().nodesBottomLeftOpenEdge_impl();
        array_type::tensor<size_t, 1> botRgt = derived_cast().nodesBottomRightOpenEdge_impl();
        array_type::tensor<size_t, 1> topLft = derived_cast().nodesTopLeftOpenEdge_impl();
        array_type::tensor<size_t, 1> topRgt = derived_cast().nodesTopRightOpenEdge_impl();

        size_t tface = fro.size() + lft.size() + bot.size();
        size_t tedge = 3 * froBot.size() + 3 * froLft.size() + 3 * botLft.size();
        size_t tnode = 7;
        array_type::tensor<size_t, 2> ret =
            xt::empty<size_t>({tface + tedge + tnode, std::size_t(2)});

        size_t i = 0;

        ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
        ret(i, 1) = derived_cast().nodesFrontBottomRightCorner_impl();
        ++i;
        ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
        ret(i, 1) = derived_cast().nodesBackBottomRightCorner_impl();
        ++i;
        ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
        ret(i, 1) = derived_cast().nodesBackBottomLeftCorner_impl();
        ++i;
        ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
        ret(i, 1) = derived_cast().nodesFrontTopLeftCorner_impl();
        ++i;
        ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
        ret(i, 1) = derived_cast().nodesFrontTopRightCorner_impl();
        ++i;
        ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
        ret(i, 1) = derived_cast().nodesBackTopRightCorner_impl();
        ++i;
        ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
        ret(i, 1) = derived_cast().nodesBackTopLeftCorner_impl();
        ++i;

        for (size_t j = 0; j < froBot.size(); ++j) {
            ret(i, 0) = froBot(j);
            ret(i, 1) = bckBot(j);
            ++i;
        }
        for (size_t j = 0; j < froBot.size(); ++j) {
            ret(i, 0) = froBot(j);
            ret(i, 1) = bckTop(j);
            ++i;
        }
        for (size_t j = 0; j < froBot.size(); ++j) {
            ret(i, 0) = froBot(j);
            ret(i, 1) = froTop(j);
            ++i;
        }
        for (size_t j = 0; j < botLft.size(); ++j) {
            ret(i, 0) = botLft(j);
            ret(i, 1) = botRgt(j);
            ++i;
        }
        for (size_t j = 0; j < botLft.size(); ++j) {
            ret(i, 0) = botLft(j);
            ret(i, 1) = topRgt(j);
            ++i;
        }
        for (size_t j = 0; j < botLft.size(); ++j) {
            ret(i, 0) = botLft(j);
            ret(i, 1) = topLft(j);
            ++i;
        }
        for (size_t j = 0; j < froLft.size(); ++j) {
            ret(i, 0) = froLft(j);
            ret(i, 1) = froRgt(j);
            ++i;
        }
        for (size_t j = 0; j < froLft.size(); ++j) {
            ret(i, 0) = froLft(j);
            ret(i, 1) = bckRgt(j);
            ++i;
        }
        for (size_t j = 0; j < froLft.size(); ++j) {
            ret(i, 0) = froLft(j);
            ret(i, 1) = bckLft(j);
            ++i;
        }

        for (size_t j = 0; j < fro.size(); ++j) {
            ret(i, 0) = fro(j);
            ret(i, 1) = bck(j);
            ++i;
        }
        for (size_t j = 0; j < lft.size(); ++j) {
            ret(i, 0) = lft(j);
            ret(i, 1) = rgt(j);
            ++i;
        }
        for (size_t j = 0; j < bot.size(); ++j) {
            ret(i, 0) = bot(j);
            ret(i, 1) = top(j);
            ++i;
        }

        return ret;
    }

    auto nodesOrigin_impl() const
    {
        return derived_cast().nodesFrontBottomLeftCorner_impl();
    }
};

/**
 * Find overlapping nodes. The output has the following structure:
 *
 *      [[nodes_from_mesh_a],
 *      [nodes_from_mesh_b]]
 *
 * @param coor_a Nodal coordinates of mesh "a" [nnode, ndim].
 * @param coor_b Nodal coordinates of mesh "b" [nnode, ndim].
 * @param rtol Relative tolerance for position match.
 * @param atol Absolute tolerance for position match.
 * @return Overlapping nodes.
 */
template <class S, class T>
inline array_type::tensor<size_t, 2>
overlapping(const S& coor_a, const T& coor_b, double rtol = 1e-5, double atol = 1e-8)
{
    GOOSEFEM_ASSERT(coor_a.dimension() == 2);
    GOOSEFEM_ASSERT(coor_b.dimension() == 2);
    GOOSEFEM_ASSERT(coor_a.shape(1) == coor_b.shape(1));

    std::vector<size_t> ret_a;
    std::vector<size_t> ret_b;

    for (size_t i = 0; i < coor_a.shape(0); ++i) {

        auto idx = xt::flatten_indices(xt::argwhere(
            xt::prod(xt::isclose(coor_b, xt::view(coor_a, i, xt::all()), rtol, atol), 1)));

        for (auto& j : idx) {
            ret_a.push_back(i);
            ret_b.push_back(j);
        }
    }

    array_type::tensor<size_t, 2> ret = xt::empty<size_t>({size_t(2), ret_a.size()});
    for (size_t i = 0; i < ret_a.size(); ++i) {
        ret(0, i) = ret_a[i];
        ret(1, i) = ret_b[i];
    }

    return ret;
}

/**
 * Stitch two mesh objects, specifying overlapping nodes by hand.
 */
class ManualStitch {
public:
    ManualStitch() = default;

    /**
     * @param coor_a Nodal coordinates of mesh "a"  [nnode, ndim].
     * @param conn_a Connectivity of mesh "a" [nelem, nne].
     * @param overlapping_nodes_a Node-numbers of mesh "a" that overlap with mesh "b" [n].
     * @param coor_b Nodal coordinates of mesh "b"  [nnode, ndim].
     * @param conn_b Connectivity of mesh "b" [nelem, nne].
     * @param overlapping_nodes_b Node-numbers of mesh "b" that overlap with mesh "a" [n].
     * @param check_position If `true` the nodes are checked for position overlap.
     * @param rtol Relative tolerance for check on position overlap.
     * @param atol Absolute tolerance for check on position overlap.
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
        double atol = 1e-8)
    {
        UNUSED(rtol);
        UNUSED(atol);

        GOOSEFEM_ASSERT(coor_a.dimension() == 2);
        GOOSEFEM_ASSERT(conn_a.dimension() == 2);
        GOOSEFEM_ASSERT(overlapping_nodes_a.dimension() == 1);
        GOOSEFEM_ASSERT(coor_b.dimension() == 2);
        GOOSEFEM_ASSERT(conn_b.dimension() == 2);
        GOOSEFEM_ASSERT(overlapping_nodes_b.dimension() == 1);
        GOOSEFEM_ASSERT(xt::has_shape(overlapping_nodes_a, overlapping_nodes_b.shape()));
        GOOSEFEM_ASSERT(coor_a.shape(1) == coor_b.shape(1));
        GOOSEFEM_ASSERT(conn_a.shape(1) == conn_b.shape(1));

        if (check_position) {
            GOOSEFEM_CHECK(xt::allclose(
                xt::view(coor_a, xt::keep(overlapping_nodes_a), xt::all()),
                xt::view(coor_b, xt::keep(overlapping_nodes_b), xt::all()),
                rtol,
                atol));
        }

        size_t nnda = coor_a.shape(0);
        size_t nndb = coor_b.shape(0);
        size_t ndim = coor_a.shape(1);
        size_t nelim = overlapping_nodes_a.size();

        size_t nela = conn_a.shape(0);
        size_t nelb = conn_b.shape(0);
        size_t nne = conn_a.shape(1);

        m_nel_a = nela;
        m_nel_b = nelb;
        m_nnd_a = nnda;

        array_type::tensor<size_t, 1> keep_b =
            xt::setdiff1d(xt::arange<size_t>(nndb), overlapping_nodes_b);

        m_map_b = xt::empty<size_t>({nndb});
        xt::view(m_map_b, xt::keep(overlapping_nodes_b)) = overlapping_nodes_a;
        xt::view(m_map_b, xt::keep(keep_b)) = xt::arange<size_t>(keep_b.size()) + nnda;

        m_conn = xt::empty<size_t>({nela + nelb, nne});
        xt::view(m_conn, xt::range(0, nela), xt::all()) = conn_a;
        xt::view(m_conn, xt::range(nela, nela + nelb), xt::all()) = detail::renum(conn_b, m_map_b);

        m_coor = xt::empty<size_t>({nnda + nndb - nelim, ndim});
        xt::view(m_coor, xt::range(0, nnda), xt::all()) = coor_a;
        xt::view(m_coor, xt::range(nnda, nnda + nndb - nelim), xt::all()) =
            xt::view(coor_b, xt::keep(keep_b), xt::all());
    }

    /**
     * Number of sub meshes == 2.
     * @return unsigned int
     */
    size_t nmesh() const
    {
        return 2;
    }

    /**
     * Number of elements.
     * @return unsigned int
     */
    size_t nelem() const
    {
        return m_conn.shape(0);
    }

    /**
     * Number of nodes.
     * @return unsigned int
     */
    size_t nnode() const
    {
        return m_coor.shape(0);
    }

    /**
     * Number of nodes-per-element.
     * @return unsigned int
     */
    size_t nne() const
    {
        return m_conn.shape(1);
    }

    /**
     * Number of dimensions.
     * @return unsigned int
     */
    size_t ndim() const
    {
        return m_coor.shape(1);
    }

    /**
     * Nodal coordinates [#nnode, #ndim].
     * @return coordinates per node
     */
    const array_type::tensor<double, 2>& coor() const
    {
        return m_coor;
    }

    /**
     * Connectivity [#nelem, #nne].
     * @return nodes per element
     */
    const array_type::tensor<size_t, 2>& conn() const
    {
        return m_conn;
    }

    /**
     * DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
     * @return DOFs per node
     */
    array_type::tensor<size_t, 2> dofs() const
    {
        size_t nnode = this->nnode();
        size_t ndim = this->ndim();
        return xt::reshape_view(xt::arange<size_t>(nnode * ndim), {nnode, ndim});
    }

    /**
     * Node-map per sub-mesh.
     * @return nodes per mesh
     */
    std::vector<array_type::tensor<size_t, 1>> nodemap() const
    {
        std::vector<array_type::tensor<size_t, 1>> ret(this->nmesh());
        for (size_t i = 0; i < this->nmesh(); ++i) {
            ret[i] = this->nodemap(i);
        }
        return ret;
    }

    /**
     * Element-map per sub-mesh.
     * @return elements per mesh
     */
    std::vector<array_type::tensor<size_t, 1>> elemmap() const
    {
        std::vector<array_type::tensor<size_t, 1>> ret(this->nmesh());
        for (size_t i = 0; i < this->nmesh(); ++i) {
            ret[i] = this->elemmap(i);
        }
        return ret;
    }

    /**
     * @param mesh_index Index of the mesh ("a" = 1, "b" = 1).
     * @return Node-map for a given mesh.
     */
    array_type::tensor<size_t, 1> nodemap(size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index <= 1);

        if (mesh_index == 0) {
            return xt::arange<size_t>(m_nnd_a);
        }

        return m_map_b;
    }

    /**
     * @param mesh_index Index of the mesh ("a" = 1, "b" = 1).
     * @return Element-map for a given mesh.
     */
    array_type::tensor<size_t, 1> elemmap(size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index <= 1);

        if (mesh_index == 0) {
            return xt::arange<size_t>(m_nel_a);
        }

        return xt::arange<size_t>(m_nel_b) + m_nel_a;
    }

    /**
     * Convert set of node numbers for an original mesh to the stitched mesh.
     *
     * @param set List of node numbers.
     * @param mesh_index Index of the mesh ("a" = 1, "b" = 1).
     * @return List of node numbers for the stitched mesh.
     */
    template <class T>
    T nodeset(const T& set, size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index <= 1);

        if (mesh_index == 0) {
            GOOSEFEM_ASSERT(xt::amax(set)() < m_nnd_a);
            return set;
        }

        GOOSEFEM_ASSERT(xt::amax(set)() < m_map_b.size());
        return detail::renum(set, m_map_b);
    }

    /**
     * Convert set of element numbers for an original mesh to the stitched mesh.
     *
     * @param set List of element numbers.
     * @param mesh_index Index of the mesh ("a" = 1, "b" = 1).
     * @return List of element numbers for the stitched mesh.
     */
    template <class T>
    T elemset(const T& set, size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index <= 1);

        if (mesh_index == 0) {
            GOOSEFEM_ASSERT(xt::amax(set)() < m_nel_a);
            return set;
        }

        GOOSEFEM_ASSERT(xt::amax(set)() < m_nel_b);
        return set + m_nel_a;
    }

private:
    array_type::tensor<double, 2> m_coor;
    array_type::tensor<size_t, 2> m_conn;
    array_type::tensor<size_t, 1> m_map_b;
    size_t m_nnd_a;
    size_t m_nel_a;
    size_t m_nel_b;
};

/**
 * Stitch mesh objects, automatically searching for overlapping nodes.
 */
class Stitch {
public:
    /**
     * @param rtol Relative tolerance for position match.
     * @param atol Absolute tolerance for position match.
     */
    Stitch(double rtol = 1e-5, double atol = 1e-8)
    {
        m_rtol = rtol;
        m_atol = atol;
    }

    /**
     * Add mesh to be stitched.
     *
     * @param coor Nodal coordinates [nnode, ndim].
     * @param conn Connectivity [nelem, nne].
     */
    template <class C, class E>
    void push_back(const C& coor, const E& conn)
    {
        GOOSEFEM_ASSERT(coor.dimension() == 2);
        GOOSEFEM_ASSERT(conn.dimension() == 2);

        if (m_map.size() == 0) {
            m_coor = coor;
            m_conn = conn;
            m_map.push_back(xt::eval(xt::arange<size_t>(coor.shape(0))));
            m_nel.push_back(conn.shape(0));
            m_el_offset.push_back(0);
            return;
        }

        auto overlap = overlapping(m_coor, coor, m_rtol, m_atol);
        size_t index = m_map.size();

        ManualStitch stitch(
            m_coor,
            m_conn,
            xt::eval(xt::view(overlap, 0, xt::all())),
            coor,
            conn,
            xt::eval(xt::view(overlap, 1, xt::all())),
            false);

        m_coor = stitch.coor();
        m_conn = stitch.conn();
        m_map.push_back(stitch.nodemap(1));
        m_nel.push_back(conn.shape(0));
        m_el_offset.push_back(m_el_offset[index - 1] + m_nel[index - 1]);
    }

    /**
     * Number of sub meshes.
     * @return unsigned int
     */
    size_t nmesh() const
    {
        return m_map.size();
    }

    /**
     * Number of elements.
     * @return unsigned int
     */
    size_t nelem() const
    {
        return m_conn.shape(0);
    }

    /**
     * Number of nodes.
     * @return unsigned int
     */
    size_t nnode() const
    {
        return m_coor.shape(0);
    }

    /**
     * Number of nodes-per-element.
     * @return unsigned int
     */
    size_t nne() const
    {
        return m_conn.shape(1);
    }

    /**
     * Number of dimensions.
     * @return unsigned int
     */
    size_t ndim() const
    {
        return m_coor.shape(1);
    }

    /**
     * Nodal coordinates [#nnode, #ndim].
     * @return coordinates per node
     */
    const array_type::tensor<double, 2>& coor() const
    {
        return m_coor;
    }

    /**
     * Connectivity [#nelem, #nne].
     * @return nodes per element
     */
    const array_type::tensor<size_t, 2>& conn() const
    {
        return m_conn;
    }

    /**
     * DOF numbers for each node (numbered sequentially) [#nnode, #ndim].
     * @return DOFs per node
     */
    array_type::tensor<size_t, 2> dofs() const
    {
        size_t nnode = this->nnode();
        size_t ndim = this->ndim();
        return xt::reshape_view(xt::arange<size_t>(nnode * ndim), {nnode, ndim});
    }

    /**
     * Node-map per sub-mesh.
     * @return nodes per mesh
     */
    std::vector<array_type::tensor<size_t, 1>> nodemap() const
    {
        std::vector<array_type::tensor<size_t, 1>> ret(this->nmesh());
        for (size_t i = 0; i < this->nmesh(); ++i) {
            ret[i] = this->nodemap(i);
        }
        return ret;
    }

    /**
     * Element-map per sub-mesh.
     * @return elements per mesh
     */
    std::vector<array_type::tensor<size_t, 1>> elemmap() const
    {
        std::vector<array_type::tensor<size_t, 1>> ret(this->nmesh());
        for (size_t i = 0; i < this->nmesh(); ++i) {
            ret[i] = this->elemmap(i);
        }
        return ret;
    }

    /**
     * The node numbers in the stitched mesh that are coming from a specific sub-mesh.
     *
     * @param mesh_index Index of the sub-mesh.
     * @return List of node numbers.
     */
    array_type::tensor<size_t, 1> nodemap(size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index < m_map.size());
        return m_map[mesh_index];
    }

    /**
     * The element numbers in the stitched mesh that are coming from a specific sub-mesh.
     *
     * @param mesh_index Index of the sub-mesh.
     * @return List of element numbers.
     */
    array_type::tensor<size_t, 1> elemmap(size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index < m_map.size());
        return xt::arange<size_t>(m_nel[mesh_index]) + m_el_offset[mesh_index];
    }

    /**
     * Convert set of node-numbers for a sub-mesh to the stitched mesh.
     *
     * @param set List of node numbers.
     * @param mesh_index Index of the sub-mesh.
     * @return List of node numbers for the stitched mesh.
     */
    template <class T>
    T nodeset(const T& set, size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index < m_map.size());
        GOOSEFEM_ASSERT(xt::amax(set)() < m_map[mesh_index].size());
        return detail::renum(set, m_map[mesh_index]);
    }

    /**
     * Convert set of element-numbers for a sub-mesh to the stitched mesh.
     *
     * @param set List of element numbers.
     * @param mesh_index Index of the sub-mesh.
     * @return List of element numbers for the stitched mesh.
     */
    template <class T>
    T elemset(const T& set, size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index < m_map.size());
        GOOSEFEM_ASSERT(xt::amax(set)() < m_nel[mesh_index]);
        return set + m_el_offset[mesh_index];
    }

    /**
     * Combine set of node numbers for an original to the final mesh (removes duplicates).
     *
     * @param set List of node numbers per mesh.
     * @return List of node numbers for the stitched mesh.
     */
    template <class T>
    T nodeset(const std::vector<T>& set) const
    {
        GOOSEFEM_ASSERT(set.size() == m_map.size());

        size_t n = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            n += set[i].size();
        }

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        n = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            xt::view(ret, xt::range(n, n + set[i].size())) = this->nodeset(set[i], i);
            n += set[i].size();
        }

        return xt::unique(ret);
    }

    /**
     * @copydoc nodeset(const std::vector<T>&) const
     */
    template <class T>
    T nodeset(std::initializer_list<T> set) const
    {
        return this->nodeset(std::vector<T>(set));
    }

    /**
     * Combine set of element numbers for an original to the final mesh.
     *
     * @param set List of element numbers per mesh.
     * @return List of element numbers for the stitched mesh.
     */
    template <class T>
    T elemset(const std::vector<T>& set) const
    {
        GOOSEFEM_ASSERT(set.size() == m_map.size());

        size_t n = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            n += set[i].size();
        }

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        n = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            xt::view(ret, xt::range(n, n + set[i].size())) = this->elemset(set[i], i);
            n += set[i].size();
        }

        return ret;
    }

    /**
     * @copydoc elemset(const std::vector<T>&) const
     */
    template <class T>
    T elemset(std::initializer_list<T> set) const
    {
        return this->elemset(std::vector<T>(set));
    }

protected:
    array_type::tensor<double, 2> m_coor; ///< Nodal coordinates [#nnode, #ndim]
    array_type::tensor<size_t, 2> m_conn; ///< Connectivity [#nelem, #nne]
    std::vector<array_type::tensor<size_t, 1>> m_map; ///< See nodemap(size_t)
    std::vector<size_t> m_nel; ///< Number of elements per sub-mesh.
    std::vector<size_t> m_el_offset; ///< First element of every sub-mesh.
    double m_rtol; ///< Relative tolerance to find overlapping nodes.
    double m_atol; ///< Absolute tolerance to find overlapping nodes.
};

/**
 * Vertically stack meshes.
 */
class Vstack : public Stitch {
public:
    /**
     * @param check_overlap Check if nodes are overlapping when adding a mesh.
     * @param rtol Relative tolerance for position match.
     * @param atol Absolute tolerance for position match.
     */
    Vstack(bool check_overlap = true, double rtol = 1e-5, double atol = 1e-8)
    {
        m_check_overlap = check_overlap;
        m_rtol = rtol;
        m_atol = atol;
    }

    /**
     * Add a mesh to the top of the current stack.
     * Each time the current `nodes_bot` are stitched with the then highest `nodes_top`.
     *
     * @param coor Nodal coordinates [nnode, ndim].
     * @param conn Connectivity [nelem, nne].
     * @param nodes_bot Nodes along the bottom edge [n].
     * @param nodes_top Nodes along the top edge [n].
     */
    template <class C, class E, class N>
    void push_back(const C& coor, const E& conn, const N& nodes_bot, const N& nodes_top)
    {
        if (m_map.size() == 0) {
            m_coor = coor;
            m_conn = conn;
            m_map.push_back(xt::eval(xt::arange<size_t>(coor.shape(0))));
            m_nel.push_back(conn.shape(0));
            m_el_offset.push_back(0);
            m_nodes_bot.push_back(nodes_bot);
            m_nodes_top.push_back(nodes_top);
            return;
        }

        GOOSEFEM_ASSERT(nodes_bot.size() == m_nodes_top.back().size());

        size_t index = m_map.size();

        double shift = xt::amax(xt::view(m_coor, xt::all(), 1))();
        auto x = coor;
        xt::view(x, xt::all(), 1) += shift;

        ManualStitch stitch(
            m_coor,
            m_conn,
            m_nodes_top.back(),
            x,
            conn,
            nodes_bot,
            m_check_overlap,
            m_rtol,
            m_atol);

        m_nodes_bot.push_back(stitch.nodeset(nodes_bot, 1));
        m_nodes_top.push_back(stitch.nodeset(nodes_top, 1));

        m_coor = stitch.coor();
        m_conn = stitch.conn();
        m_map.push_back(stitch.nodemap(1));
        m_nel.push_back(conn.shape(0));
        m_el_offset.push_back(m_el_offset[index - 1] + m_nel[index - 1]);
    }

private:
    std::vector<array_type::tensor<size_t, 1>>
        m_nodes_bot; ///< Bottom nodes of each mesh (renumbered).
    std::vector<array_type::tensor<size_t, 1>>
        m_nodes_top; ///< Top nodes of each mesh (renumbered).
    bool m_check_overlap; ///< Check if nodes are overlapping when adding a mesh.
};

/**
 * Reorder to lowest possible index, in specific order.
 *
 * For example for ``Reorder({iiu, iip})`` after reordering:
 *
 *      iiu = xt::range<size_t>(nnu);
 *      iip = xt::range<size_t>(nnp) + nnu;
 */
class Reorder {
public:
    Reorder() = default;

    /**
     * @param args List of (DOF-)numbers.
     */
    template <class T>
    Reorder(const std::initializer_list<T> args)
    {
        size_t n = 0;
        size_t i = 0;

        for (auto& arg : args) {
            if (arg.size() == 0) {
                continue;
            }
            n = std::max(n, xt::amax(arg)() + 1);
        }

#ifdef GOOSEFEM_ENABLE_ASSERT
        for (auto& arg : args) {
            GOOSEFEM_ASSERT(is_unique(arg));
        }
#endif

        m_renum = xt::empty<size_t>({n});

        for (auto& arg : args) {
            for (auto& j : arg) {
                m_renum(j) = i;
                ++i;
            }
        }
    }

    /**
     * Apply reordering to other set.
     *
     * @param list List of (DOF-)numbers.
     * @return Reordered list of (DOF-)numbers.
     */
    template <class T>
    T apply(const T& list) const
    {
        T ret = T::from_shape(list.shape());

        auto jt = ret.begin();

        for (auto it = list.begin(); it != list.end(); ++it, ++jt) {
            *jt = m_renum(*it);
        }

        return ret;
    }

    /**
     * Get the list needed to reorder, e.g.:
     *
     *      dofs_reordered(i, j) = index(dofs(i, j))
     *
     * @return Reorder-index.
     */
    const array_type::tensor<size_t, 1>& index() const
    {
        return m_renum;
    }

private:
    array_type::tensor<size_t, 1> m_renum;
};

/**
 * Number of elements connected to each node.
 *
 * @param conn Connectivity [nelem, nne].
 * @return Coordination per node.
 */
template <class E>
inline array_type::tensor<size_t, 1> coordination(const E& conn)
{
    GOOSEFEM_ASSERT(conn.dimension() == 2);

    size_t nnode = xt::amax(conn)() + 1;

    array_type::tensor<size_t, 1> N = xt::zeros<size_t>({nnode});

    for (auto it = conn.begin(); it != conn.end(); ++it) {
        N(*it) += 1;
    }

    return N;
}

/**
 * Elements connected to each node.
 *
 * @param conn Connectivity [nelem, nne].
 * @param sorted If `true` the list of elements for each node is sorted.
 * @return Elements per node [nnode, ...].
 */
template <class E>
inline std::vector<std::vector<size_t>> elem2node(const E& conn, bool sorted = true)
{
    auto N = coordination(conn);
    auto nnode = N.size();

    std::vector<std::vector<size_t>> ret(nnode);
    for (size_t i = 0; i < nnode; ++i) {
        ret[i].reserve(N(i));
    }

    for (size_t e = 0; e < conn.shape(0); ++e) {
        for (size_t m = 0; m < conn.shape(1); ++m) {
            ret[conn(e, m)].push_back(e);
        }
    }

    if (sorted) {
        for (auto& row : ret) {
            std::sort(row.begin(), row.end());
        }
    }

    return ret;
}

/**
 * @copydoc elem2node(const E&, bool)
 *
 * @param dofs DOFs per node, allowing accounting for periodicity [nnode, ndim].
 */
template <class E, class D>
inline std::vector<std::vector<size_t>> elem2node(const E& conn, const D& dofs, bool sorted = true)
{
    size_t nnode = dofs.shape(0);
    auto ret = elem2node(conn, sorted);
    auto nties = nodaltyings(dofs);

    for (size_t m = 0; m < nnode; ++m) {
        if (nties[m].size() <= 1) {
            continue;
        }
        if (m > nties[m][0]) {
            continue;
        }
        size_t n = ret[m].size();
        for (size_t j = 1; j < nties[m].size(); ++j) {
            n += ret[nties[m][j]].size();
        }
        ret[m].reserve(n);
        for (size_t j = 1; j < nties[m].size(); ++j) {
            ret[m].insert(ret[m].end(), ret[nties[m][j]].begin(), ret[nties[m][j]].end());
        }
        if (sorted) {
            std::sort(ret[m].begin(), ret[m].end());
        }
        for (size_t j = 1; j < nties[m].size(); ++j) {
            ret[nties[m][j]] = ret[m];
        }
    }

    return ret;
}

/**
 * Nodes connected to each DOF.
 *
 * @param dofs DOFs per node [nnode, ndim].
 * @param sorted If `true` the list of nodes for each DOF is sorted.
 * @return Nodes per DOF [ndof, ...].
 */
template <class D>
inline std::vector<std::vector<size_t>> node2dof(const D& dofs, bool sorted = true)
{
    return elem2node(dofs, sorted);
}

/**
 * Return size of each element edge.
 *
 * @param coor Nodal coordinates.
 * @param conn Connectivity.
 * @param type ElementType.
 * @return Edge-sizes per element.
 */
template <class C, class E>
inline array_type::tensor<double, 2> edgesize(const C& coor, const E& conn, ElementType type)
{
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);
    GOOSEFEM_ASSERT(xt::amax(conn)() < coor.shape(0));

    if (type == ElementType::Quad4) {
        GOOSEFEM_ASSERT(coor.shape(1) == 2ul);
        GOOSEFEM_ASSERT(conn.shape(1) == 4ul);
        array_type::tensor<size_t, 1> n0 = xt::view(conn, xt::all(), 0);
        array_type::tensor<size_t, 1> n1 = xt::view(conn, xt::all(), 1);
        array_type::tensor<size_t, 1> n2 = xt::view(conn, xt::all(), 2);
        array_type::tensor<size_t, 1> n3 = xt::view(conn, xt::all(), 3);
        array_type::tensor<double, 1> x0 = xt::view(coor, xt::keep(n0), 0);
        array_type::tensor<double, 1> x1 = xt::view(coor, xt::keep(n1), 0);
        array_type::tensor<double, 1> x2 = xt::view(coor, xt::keep(n2), 0);
        array_type::tensor<double, 1> x3 = xt::view(coor, xt::keep(n3), 0);
        array_type::tensor<double, 1> y0 = xt::view(coor, xt::keep(n0), 1);
        array_type::tensor<double, 1> y1 = xt::view(coor, xt::keep(n1), 1);
        array_type::tensor<double, 1> y2 = xt::view(coor, xt::keep(n2), 1);
        array_type::tensor<double, 1> y3 = xt::view(coor, xt::keep(n3), 1);
        array_type::tensor<double, 2> ret = xt::empty<double>(conn.shape());
        xt::view(ret, xt::all(), 0) = xt::sqrt(xt::pow(x1 - x0, 2.0) + xt::pow(y1 - y0, 2.0));
        xt::view(ret, xt::all(), 1) = xt::sqrt(xt::pow(x2 - x1, 2.0) + xt::pow(y2 - y1, 2.0));
        xt::view(ret, xt::all(), 2) = xt::sqrt(xt::pow(x3 - x2, 2.0) + xt::pow(y3 - y2, 2.0));
        xt::view(ret, xt::all(), 3) = xt::sqrt(xt::pow(x0 - x3, 2.0) + xt::pow(y0 - y3, 2.0));
        return ret;
    }

    throw std::runtime_error("Element-type not implemented");
}

/**
 * Return size of each element edge.
 * The element-type is automatically determined, see defaultElementType().
 *
 * @param coor Nodal coordinates.
 * @param conn Connectivity.
 * @return Edge-sizes per element.
 */
template <class C, class E>
inline array_type::tensor<double, 2> edgesize(const C& coor, const E& conn)
{
    return edgesize(coor, conn, defaultElementType(coor, conn));
}

/**
 * Coordinates of the center of each element.
 *
 * @param coor Nodal coordinates.
 * @param conn Connectivity.
 * @param type ElementType.
 * @return Center of each element.
 */
template <class C, class E>
inline array_type::tensor<double, 2> centers(const C& coor, const E& conn, ElementType type)
{
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);
    GOOSEFEM_ASSERT(xt::amax(conn)() < coor.shape(0));
    array_type::tensor<double, 2> ret = xt::zeros<double>({conn.shape(0), coor.shape(1)});

    if (type == ElementType::Quad4) {
        GOOSEFEM_ASSERT(coor.shape(1) == 2);
        GOOSEFEM_ASSERT(conn.shape(1) == 4);
        for (size_t i = 0; i < 4; ++i) {
            auto n = xt::view(conn, xt::all(), i);
            ret += xt::view(coor, xt::keep(n), xt::all());
        }
        ret /= 4.0;
        return ret;
    }

    throw std::runtime_error("Element-type not implemented");
}

/**
 * Coordinates of the center of each element.
 * The element-type is automatically determined, see defaultElementType().
 *
 * @param coor Nodal coordinates.
 * @param conn Connectivity.
 * @return Center of each element.
 */
template <class C, class E>
inline array_type::tensor<double, 2> centers(const C& coor, const E& conn)
{
    return centers(coor, conn, defaultElementType(coor, conn));
}

/**
 * Convert an element-map to a node-map.
 *
 * @param elem_map Element-map such that `new_elvar = elvar[elem_map]`.
 * @param coor Nodal coordinates.
 * @param conn Connectivity.
 * @param type ElementType.
 * @return Node-map such that `new_nodevar = nodevar[node_map]`
 */
template <class T, class C, class E>
inline array_type::tensor<size_t, 1>
elemmap2nodemap(const T& elem_map, const C& coor, const E& conn, ElementType type)
{
    GOOSEFEM_ASSERT(elem_map.dimension() == 1);
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);
    GOOSEFEM_ASSERT(xt::amax(conn)() < coor.shape(0));
    GOOSEFEM_ASSERT(elem_map.size() == conn.shape(0));
    size_t N = coor.shape(0);

    array_type::tensor<size_t, 1> ret = N * xt::ones<size_t>({N});

    if (type == ElementType::Quad4) {
        GOOSEFEM_ASSERT(coor.shape(1) == 2);
        GOOSEFEM_ASSERT(conn.shape(1) == 4);

        for (size_t i = 0; i < 4; ++i) {
            array_type::tensor<size_t, 1> t = N * xt::ones<size_t>({N});
            auto old_nd = xt::view(conn, xt::all(), i);
            auto new_nd = xt::view(conn, xt::keep(elem_map), i);
            xt::view(t, xt::keep(old_nd)) = new_nd;
            ret = xt::where(xt::equal(ret, N), t, ret);
        }

        return ret;
    }

    throw std::runtime_error("Element-type not implemented");
}

/**
 * Convert an element-map to a node-map.
 * The element-type is automatically determined, see defaultElementType().
 *
 * @param elem_map Element-map such that `new_elvar = elvar[elem_map]`.
 * @param coor Nodal coordinates.
 * @param conn Connectivity.
 * @return Node-map such that `new_nodevar = nodevar[node_map]`
 */
template <class T, class C, class E>
inline array_type::tensor<size_t, 1>
elemmap2nodemap(const T& elem_map, const C& coor, const E& conn)
{
    return elemmap2nodemap(elem_map, coor, conn, defaultElementType(coor, conn));
}

/**
 * Compute the mass of each node in the mesh.
 * If nodes are not part of the connectivity the mass is set to zero,
 * such that the center of gravity is simply::
 *
 *      average(coor, GooseFEM.Mesh.nodal_mass(coor, conn, type), axis=0);
 *
 * @tparam C e.g. `array_type::tensor<double, 2>`
 * @tparam E e.g. `array_type::tensor<size_t, 2>`
 * @param coor Nodal coordinates `[nnode, ndim]`.
 * @param conn Connectivity `[nelem, nne]`.
 * @param type ElementType.
 * @return Center of gravity `[ndim]`.
 */
template <class C, class E>
inline array_type::tensor<double, 2> nodal_mass(const C& coor, const E& conn, ElementType type)
{
    auto dof = dofs(coor.shape(0), coor.shape(1));
    GooseFEM::MatrixDiagonal M(conn, dof);
    GooseFEM::Vector vector(conn, dof);
    array_type::tensor<double, 2> rho = xt::ones<double>(conn.shape());

    if (type == ElementType::Quad4) {
        GooseFEM::Element::Quad4::Quadrature quad(
            vector.AsElement(coor),
            GooseFEM::Element::Quad4::Nodal::xi(),
            GooseFEM::Element::Quad4::Nodal::w());
        M.assemble(quad.Int_N_scalar_NT_dV(rho));
    }
    else {
        throw std::runtime_error("Element-type not implemented");
    }

    return vector.AsNode(M.data());
}

/**
 * Compute the mass of each node in the mesh.
 * If nodes are not part of the connectivity the mass is set to zero,
 * such that the center of gravity is simply::
 *
 *      average(coor, GooseFEM.Mesh.nodal_mass(coor, conn), axis=0);
 *
 * @tparam C e.g. `array_type::tensor<double, 2>`
 * @tparam E e.g. `array_type::tensor<size_t, 2>`
 * @param coor Nodal coordinates `[nnode, ndim]`.
 * @param conn Connectivity `[nelem, nne]`.
 * @return Center of gravity `[ndim]`.
 */
template <class C, class E>
inline array_type::tensor<double, 2> nodal_mass(const C& coor, const E& conn)
{
    return nodal_mass(coor, conn, defaultElementType(coor, conn));
}

namespace detail {

// todo: remove after upstream fix
template <class T>
array_type::tensor<double, 1> average_axis_0(const T& data, const T& weights)
{
    GOOSEFEM_ASSERT(data.dimension() == 2);
    GOOSEFEM_ASSERT(xt::has_shape(data, weights.shape()));

    array_type::tensor<double, 1> ret = xt::zeros<double>({weights.shape(1)});

    for (size_t j = 0; j < weights.shape(1); ++j) {
        double norm = 0.0;
        for (size_t i = 0; i < weights.shape(0); ++i) {
            ret(j) += data(i, j) * weights(i, j);
            norm += weights(i, j);
        }
        ret(j) /= norm;
    }
    return ret;
}

} // namespace detail

/**
 * Compute the center of gravity of a mesh.
 *
 * @tparam C e.g. `array_type::tensor<double, 2>`
 * @tparam E e.g. `array_type::tensor<size_t, 2>`
 * @param coor Nodal coordinates `[nnode, ndim]`.
 * @param conn Connectivity `[nelem, nne]`.
 * @param type ElementType.
 * @return Center of gravity `[ndim]`.
 */
template <class C, class E>
inline array_type::tensor<double, 1>
center_of_gravity(const C& coor, const E& conn, ElementType type)
{
    // todo: remove after upstream fix
    return detail::average_axis_0(coor, nodal_mass(coor, conn, type));
    // return xt::average(coor, nodal_mass(coor, conn, type), 0);
}

/**
 * Compute the center of gravity of a mesh.
 *
 * @tparam C e.g. `array_type::tensor<double, 2>`
 * @tparam E e.g. `array_type::tensor<size_t, 2>`
 * @param coor Nodal coordinates `[nnode, ndim]`.
 * @param conn Connectivity `[nelem, nne]`.
 * @return Center of gravity `[ndim]`.
 */
template <class C, class E>
inline array_type::tensor<double, 1> center_of_gravity(const C& coor, const E& conn)
{
    // todo: remove after upstream fix
    return detail::average_axis_0(coor, nodal_mass(coor, conn, defaultElementType(coor, conn)));
    // return xt::average(coor, nodal_mass(coor, conn, defaultElementType(coor, conn)), 0);
}

/**
 * List nodal tyings based on DOF-numbers per node.
 *
 * @param dofs DOFs per node [nnode, ndim].
 * @return Nodes to which the nodes is connected (sorted) [nnode, ...]
 */
template <class D>
inline std::vector<std::vector<size_t>> nodaltyings(const D& dofs)
{
    size_t nnode = dofs.shape(0);
    size_t ndim = dofs.shape(1);
    auto nodemap = node2dof(dofs);
    std::vector<std::vector<size_t>> ret(nnode);

    for (size_t m = 0; m < nnode; ++m) {
        auto r = nodemap[dofs(m, 0)];
        std::sort(r.begin(), r.end());
        ret[m] = r;
#ifdef GOOSEFEM_ENABLE_ASSERT
        for (size_t i = 1; i < ndim; ++i) {
            auto t = nodemap[dofs(m, i)];
            std::sort(t.begin(), t.end());
            GOOSEFEM_ASSERT(std::equal(r.begin(), r.end(), t.begin()));
        }
#endif
    }

    return ret;
}

} // namespace Mesh
} // namespace GooseFEM

#endif

/**
Generate simple meshes of 3-noded triangular elements in 2d (GooseFEM::Mesh::ElementType::Tri3).

\file MeshTri3.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESHTRI3_H
#define GOOSEFEM_MESHTRI3_H

#include "Mesh.h"
#include "config.h"

namespace GooseFEM {
namespace Mesh {

/**
Simple meshes of and mesh operations for triangular elements of type ElementType::Tri3.
*/
namespace Tri3 {

/**
Regular grid of squares, with each square cut into two triangular elements.
*/
class Regular : public RegularBase2d<Regular> {
public:
    /**
    Constructor.

    \param nelx Number of elements in x-direction.
    \param nely Number of elements in y-direction.
    \param h Edge-size (of the squares, and thus of two of three element-edges).
    */
    Regular(size_t nelx, size_t nely, double h = 1.0);

private:
    friend class RegularBase<Regular>;
    friend class RegularBase2d<Regular>;

    size_t nelx_impl() const;
    size_t nely_impl() const;
    ElementType getElementType_impl() const;
    array_type::tensor<double, 2> coor_impl() const;
    array_type::tensor<size_t, 2> conn_impl() const;
    array_type::tensor<size_t, 1> nodesBottomEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopEdge_impl() const;
    array_type::tensor<size_t, 1> nodesLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesRightOpenEdge_impl() const;
    size_t nodesBottomLeftCorner_impl() const;
    size_t nodesBottomRightCorner_impl() const;
    size_t nodesTopLeftCorner_impl() const;
    size_t nodesTopRightCorner_impl() const;

    double m_h; ///< See h()
    size_t m_nelx; ///< See nelx()
    size_t m_nely; ///< See nely()
    size_t m_nelem; ///< See nelem()
    size_t m_nnode; ///< See nnode()
    size_t m_nne; ///< See nne()
    size_t m_ndim; ///< See ndim()
};

// read / set the orientation (-1/+1) of all triangles

/**
Read the orientation of a mesh of triangular elements of type ElementType::Tri3.

\param coor Nodal coordinates [nnode, ndim].
\param conn Connectivity [nelem, nne].
\return Orientation (-1 or +1) per element [nelem].
*/
inline array_type::tensor<int, 1> getOrientation(
    const array_type::tensor<double, 2>& coor,
    const array_type::tensor<size_t, 2>& conn);

/**
Set the orientation of a mesh of triangular elements of type ElementType::Tri3.

\param coor Nodal coordinates [nnode, ndim].
\param conn Connectivity [nelem, nne].
\param orientation Target orientation (applied to all elements).
\return Connectivity (order of nodes-per-element may have changed) [nelem, nne].
*/
inline array_type::tensor<size_t, 2> setOrientation(
    const array_type::tensor<double, 2>& coor,
    const array_type::tensor<size_t, 2>& conn,
    int orientation = -1);

/**
Set the orientation of a mesh of triangular elements of type ElementType::Tri3.
For efficiency this function reuses the output of getOrientation().

\param coor Nodal coordinates [nnode, ndim].
\param conn Connectivity [nelem, nne].
\param current Current orientation per element (output of getOrientation()) [nelem].
\param orientation Target orientation (applied to all elements).
\return Connectivity (order of nodes-per-element may have changed) [nelem, nne].
*/
inline array_type::tensor<size_t, 2> setOrientation(
    const array_type::tensor<double, 2>& coor,
    const array_type::tensor<size_t, 2>& conn,
    const array_type::tensor<int, 1>& current,
    int orientation = -1);

} // namespace Tri3
} // namespace Mesh
} // namespace GooseFEM

#include "MeshTri3.hpp"

#endif

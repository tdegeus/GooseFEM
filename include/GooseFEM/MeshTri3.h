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
    Regular(size_t nelx, size_t nely, double h = 1.0)
    {
        m_h = h;
        m_nelx = nelx;
        m_nely = nely;
        m_ndim = 2;
        m_nne = 3;

        GOOSEFEM_ASSERT(m_nelx >= 1);
        GOOSEFEM_ASSERT(m_nely >= 1);

        m_nnode = (m_nelx + 1) * (m_nely + 1);
        m_nelem = m_nelx * m_nely * 2;
    }

private:
    friend class RegularBase<Regular>;
    friend class RegularBase2d<Regular>;

    size_t nelx_impl() const
    {
        return m_nelx;
    }

    size_t nely_impl() const
    {
        return m_nely;
    }

    ElementType getElementType_impl() const
    {
        return ElementType::Tri3;
    }

    array_type::tensor<double, 2> coor_impl() const
    {
        array_type::tensor<double, 2> ret = xt::empty<double>({m_nnode, m_ndim});

        array_type::tensor<double, 1> x =
            xt::linspace<double>(0.0, m_h * static_cast<double>(m_nelx), m_nelx + 1);
        array_type::tensor<double, 1> y =
            xt::linspace<double>(0.0, m_h * static_cast<double>(m_nely), m_nely + 1);

        size_t inode = 0;

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
                ret(inode, 0) = x(ix);
                ret(inode, 1) = y(iy);
                ++inode;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 2> conn_impl() const
    {
        array_type::tensor<size_t, 2> ret = xt::empty<size_t>({m_nelem, m_nne});

        size_t ielem = 0;

        for (size_t iy = 0; iy < m_nely; ++iy) {
            for (size_t ix = 0; ix < m_nelx; ++ix) {
                ret(ielem, 0) = (iy) * (m_nelx + 1) + (ix);
                ret(ielem, 1) = (iy) * (m_nelx + 1) + (ix + 1);
                ret(ielem, 2) = (iy + 1) * (m_nelx + 1) + (ix);
                ++ielem;
                ret(ielem, 0) = (iy) * (m_nelx + 1) + (ix + 1);
                ret(ielem, 1) = (iy + 1) * (m_nelx + 1) + (ix + 1);
                ret(ielem, 2) = (iy + 1) * (m_nelx + 1) + (ix);
                ++ielem;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx + 1});

        for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
            ret(ix) = ix;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx + 1});

        for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
            ret(ix) = ix + m_nely * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesLeftEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely + 1});

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            ret(iy) = iy * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesRightEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely + 1});

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            ret(iy) = iy * (m_nelx + 1) + m_nelx;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx - 1});

        for (size_t ix = 1; ix < m_nelx; ++ix) {
            ret(ix - 1) = ix;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx - 1});

        for (size_t ix = 1; ix < m_nelx; ++ix) {
            ret(ix - 1) = ix + m_nely * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesLeftOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely - 1});

        for (size_t iy = 1; iy < m_nely; ++iy) {
            ret(iy - 1) = iy * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesRightOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely - 1});

        for (size_t iy = 1; iy < m_nely; ++iy) {
            ret(iy - 1) = iy * (m_nelx + 1) + m_nelx;
        }

        return ret;
    }

    size_t nodesBottomLeftCorner_impl() const
    {
        return 0;
    }

    size_t nodesBottomRightCorner_impl() const
    {
        return m_nelx;
    }

    size_t nodesTopLeftCorner_impl() const
    {
        return m_nely * (m_nelx + 1);
    }

    size_t nodesTopRightCorner_impl() const
    {
        return m_nely * (m_nelx + 1) + m_nelx;
    }

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
inline array_type::tensor<int, 1>
getOrientation(const array_type::tensor<double, 2>& coor, const array_type::tensor<size_t, 2>& conn)
{
    GOOSEFEM_ASSERT(conn.shape(1) == 3ul);
    GOOSEFEM_ASSERT(coor.shape(1) == 2ul);

    double k;
    size_t nelem = conn.shape(0);

    array_type::tensor<int, 1> ret = xt::empty<int>({nelem});

    for (size_t ielem = 0; ielem < nelem; ++ielem) {
        auto v1 =
            xt::view(coor, conn(ielem, 0), xt::all()) - xt::view(coor, conn(ielem, 1), xt::all());
        auto v2 =
            xt::view(coor, conn(ielem, 2), xt::all()) - xt::view(coor, conn(ielem, 1), xt::all());

        k = v1(0) * v2(1) - v2(0) * v1(1);

        if (k < 0) {
            ret(ielem) = -1;
        }
        else {
            ret(ielem) = +1;
        }
    }

    return ret;
}

/**
Set the orientation of a mesh of triangular elements of type ElementType::Tri3.
For efficiency this function reuses the output of getOrientation().

\param coor Nodal coordinates [nnode, ndim].
\param conn Connectivity [nelem, nne].
\param val Current orientation per element (output of getOrientation()) [nelem].
\param orientation Target orientation (applied to all elements).
\return Connectivity (order of nodes-per-element may have changed) [nelem, nne].
*/
inline array_type::tensor<size_t, 2> setOrientation(
    const array_type::tensor<double, 2>& coor,
    const array_type::tensor<size_t, 2>& conn,
    const array_type::tensor<int, 1>& val,
    int orientation = -1)
{
    GOOSEFEM_ASSERT(conn.shape(1) == 3ul);
    GOOSEFEM_ASSERT(coor.shape(1) == 2ul);
    GOOSEFEM_ASSERT(conn.shape(0) == val.size());
    GOOSEFEM_ASSERT(orientation == -1 || orientation == +1);

    UNUSED(coor);

    size_t nelem = conn.shape(0);
    array_type::tensor<size_t, 2> ret = conn;

    for (size_t ielem = 0; ielem < nelem; ++ielem) {
        if ((orientation == -1 && val(ielem) > 0) || (orientation == +1 && val(ielem) < 0)) {
            std::swap(ret(ielem, 2), ret(ielem, 1));
        }
    }

    return ret;
}

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
    int orientation = -1)
{
    GOOSEFEM_ASSERT(conn.shape(1) == 3ul);
    GOOSEFEM_ASSERT(coor.shape(1) == 2ul);
    GOOSEFEM_ASSERT(orientation == -1 || orientation == +1);

    array_type::tensor<int, 1> val = getOrientation(coor, conn);

    return setOrientation(coor, conn, val, orientation);
}

} // namespace Tri3
} // namespace Mesh
} // namespace GooseFEM

#endif

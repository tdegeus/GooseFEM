/**
Implementation of MeshTri3.h

\file MeshTri3.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESHTRI3_HPP
#define GOOSEFEM_MESHTRI3_HPP

#include "MeshTri3.h"

namespace GooseFEM {
namespace Mesh {
namespace Tri3 {

inline Regular::Regular(size_t nelx, size_t nely, double h)
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

inline ElementType Regular::getElementType() const
{
    return ElementType::Tri3;
}

inline xt::xtensor<double, 2> Regular::coor() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({m_nnode, m_ndim});

    xt::xtensor<double, 1> x =
        xt::linspace<double>(0.0, m_h * static_cast<double>(m_nelx), m_nelx + 1);
    xt::xtensor<double, 1> y =
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

inline xt::xtensor<size_t, 2> Regular::conn() const
{
    xt::xtensor<size_t, 2> ret = xt::empty<size_t>({m_nelem, m_nne});

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

inline xt::xtensor<size_t, 1> Regular::nodesBottomEdge() const
{
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nelx + 1});

    for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
        ret(ix) = ix;
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Regular::nodesTopEdge() const
{
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nelx + 1});

    for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
        ret(ix) = ix + m_nely * (m_nelx + 1);
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Regular::nodesLeftEdge() const
{
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nely + 1});

    for (size_t iy = 0; iy < m_nely + 1; ++iy) {
        ret(iy) = iy * (m_nelx + 1);
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Regular::nodesRightEdge() const
{
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nely + 1});

    for (size_t iy = 0; iy < m_nely + 1; ++iy) {
        ret(iy) = iy * (m_nelx + 1) + m_nelx;
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Regular::nodesBottomOpenEdge() const
{
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nelx - 1});

    for (size_t ix = 1; ix < m_nelx; ++ix) {
        ret(ix - 1) = ix;
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Regular::nodesTopOpenEdge() const
{
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nelx - 1});

    for (size_t ix = 1; ix < m_nelx; ++ix) {
        ret(ix - 1) = ix + m_nely * (m_nelx + 1);
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Regular::nodesLeftOpenEdge() const
{
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nely - 1});

    for (size_t iy = 1; iy < m_nely; ++iy) {
        ret(iy - 1) = iy * (m_nelx + 1);
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Regular::nodesRightOpenEdge() const
{
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nely - 1});

    for (size_t iy = 1; iy < m_nely; ++iy) {
        ret(iy - 1) = iy * (m_nelx + 1) + m_nelx;
    }

    return ret;
}

inline size_t Regular::nodesBottomLeftCorner() const
{
    return 0;
}

inline size_t Regular::nodesBottomRightCorner() const
{
    return m_nelx;
}

inline size_t Regular::nodesTopLeftCorner() const
{
    return m_nely * (m_nelx + 1);
}

inline size_t Regular::nodesTopRightCorner() const
{
    return m_nely * (m_nelx + 1) + m_nelx;
}

inline xt::xtensor<int, 1>
getOrientation(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn)
{
    GOOSEFEM_ASSERT(conn.shape(1) == 3ul);
    GOOSEFEM_ASSERT(coor.shape(1) == 2ul);

    double k;
    size_t nelem = conn.shape(0);

    xt::xtensor<int, 1> ret = xt::empty<int>({nelem});

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

inline xt::xtensor<size_t, 2> setOrientation(
    const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn, int orientation)
{
    GOOSEFEM_ASSERT(conn.shape(1) == 3ul);
    GOOSEFEM_ASSERT(coor.shape(1) == 2ul);
    GOOSEFEM_ASSERT(orientation == -1 || orientation == +1);

    xt::xtensor<int, 1> val = getOrientation(coor, conn);

    return setOrientation(coor, conn, val, orientation);
}

inline xt::xtensor<size_t, 2> setOrientation(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn,
    const xt::xtensor<int, 1>& val,
    int orientation)
{
    GOOSEFEM_ASSERT(conn.shape(1) == 3ul);
    GOOSEFEM_ASSERT(coor.shape(1) == 2ul);
    GOOSEFEM_ASSERT(conn.shape(0) == val.size());
    GOOSEFEM_ASSERT(orientation == -1 || orientation == +1);

    UNUSED(coor);

    size_t nelem = conn.shape(0);
    xt::xtensor<size_t, 2> ret = conn;

    for (size_t ielem = 0; ielem < nelem; ++ielem) {
        if ((orientation == -1 && val(ielem) > 0) || (orientation == +1 && val(ielem) < 0)) {
            std::swap(ret(ielem, 2), ret(ielem, 1));
        }
    }

    return ret;
}

} // namespace Tri3
} // namespace Mesh
} // namespace GooseFEM

#endif

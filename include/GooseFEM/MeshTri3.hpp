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

inline Regular::Regular(size_t nelx, size_t nely, double h) : m_h(h), m_nelx(nelx), m_nely(nely)
{
    GOOSEFEM_ASSERT(m_nelx >= 1ul);
    GOOSEFEM_ASSERT(m_nely >= 1ul);

    m_nnode = (m_nelx + 1) * (m_nely + 1);
    m_nelem = m_nelx * m_nely * 2;
}

inline size_t Regular::nelem() const
{
    return m_nelem;
}

inline size_t Regular::nnode() const
{
    return m_nnode;
}

inline size_t Regular::nne() const
{
    return m_nne;
}

inline size_t Regular::ndim() const
{
    return m_ndim;
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

inline size_t Regular::nodesLeftBottomCorner() const
{
    return nodesBottomLeftCorner();
}

inline size_t Regular::nodesLeftTopCorner() const
{
    return nodesTopLeftCorner();
}

inline size_t Regular::nodesRightBottomCorner() const
{
    return nodesBottomRightCorner();
}

inline size_t Regular::nodesRightTopCorner() const
{
    return nodesTopRightCorner();
}

inline xt::xtensor<size_t, 2> Regular::nodesPeriodic() const
{
    xt::xtensor<size_t, 1> bot = nodesBottomOpenEdge();
    xt::xtensor<size_t, 1> top = nodesTopOpenEdge();
    xt::xtensor<size_t, 1> lft = nodesLeftOpenEdge();
    xt::xtensor<size_t, 1> rgt = nodesRightOpenEdge();

    size_t tedge = bot.size() + lft.size();
    size_t tnode = 3;
    xt::xtensor<size_t, 2> ret = xt::empty<size_t>({tedge + tnode, std::size_t(2)});

    size_t i = 0;

    ret(i, 0) = nodesBottomLeftCorner();
    ret(i, 1) = nodesBottomRightCorner();
    ++i;
    ret(i, 0) = nodesBottomLeftCorner();
    ret(i, 1) = nodesTopRightCorner();
    ++i;
    ret(i, 0) = nodesBottomLeftCorner();
    ret(i, 1) = nodesTopLeftCorner();
    ++i;

    for (size_t j = 0; j < bot.size(); ++j) {
        ret(i, 0) = bot(j);
        ret(i, 1) = top(j);
        ++i;
    }
    for (size_t j = 0; j < lft.size(); ++j) {
        ret(i, 0) = lft(j);
        ret(i, 1) = rgt(j);
        ++i;
    }

    return ret;
}

inline size_t Regular::nodesOrigin() const
{
    return nodesBottomLeftCorner();
}

inline xt::xtensor<size_t, 2> Regular::dofs() const
{
    return GooseFEM::Mesh::dofs(m_nnode, m_ndim);
}

inline xt::xtensor<size_t, 2> Regular::dofsPeriodic() const
{
    xt::xtensor<size_t, 2> ret = GooseFEM::Mesh::dofs(m_nnode, m_ndim);
    xt::xtensor<size_t, 2> nodePer = nodesPeriodic();

    for (size_t i = 0; i < nodePer.shape(0); ++i) {
        for (size_t j = 0; j < m_ndim; ++j) {
            ret(nodePer(i, 1), j) = ret(nodePer(i, 0), j);
        }
    }

    return GooseFEM::Mesh::renumber(ret);
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

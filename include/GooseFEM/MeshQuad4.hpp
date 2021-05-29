/**
Implementation of MeshQuad4.h

\file MeshQuad4.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESHQUAD4_HPP
#define GOOSEFEM_MESHQUAD4_HPP

#include "MeshQuad4.h"

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {

inline Regular::Regular(size_t nelx, size_t nely, double h)
{
    m_h = h;
    m_nelx = nelx;
    m_nely = nely;
    m_ndim = 2;
    m_nne = 4;

    GOOSEFEM_ASSERT(m_nelx >= 1);
    GOOSEFEM_ASSERT(m_nely >= 1);

    m_nnode = (m_nelx + 1) * (m_nely + 1);
    m_nelem = m_nelx * m_nely;
}

inline ElementType Regular::getElementType_impl() const
{
    return ElementType::Quad4;
}

inline size_t Regular::nelx_impl() const
{
    return m_nelx;
}

inline size_t Regular::nely_impl() const
{
    return m_nely;
}

inline xt::xtensor<double, 2> Regular::coor_impl() const
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

inline xt::xtensor<size_t, 2> Regular::conn_impl() const
{
    xt::xtensor<size_t, 2> ret = xt::empty<size_t>({m_nelem, m_nne});

    size_t ielem = 0;

    for (size_t iy = 0; iy < m_nely; ++iy) {
        for (size_t ix = 0; ix < m_nelx; ++ix) {
            ret(ielem, 0) = (iy) * (m_nelx + 1) + (ix);
            ret(ielem, 1) = (iy) * (m_nelx + 1) + (ix + 1);
            ret(ielem, 3) = (iy + 1) * (m_nelx + 1) + (ix);
            ret(ielem, 2) = (iy + 1) * (m_nelx + 1) + (ix + 1);
            ++ielem;
        }
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Regular::nodesBottomEdge_impl() const
{
    return xt::arange<size_t>(m_nelx + 1);
}

inline xt::xtensor<size_t, 1> Regular::nodesTopEdge_impl() const
{
    return xt::arange<size_t>(m_nelx + 1) + m_nely * (m_nelx + 1);
}

inline xt::xtensor<size_t, 1> Regular::nodesLeftEdge_impl() const
{
    return xt::arange<size_t>(m_nely + 1) * (m_nelx + 1);
}

inline xt::xtensor<size_t, 1> Regular::nodesRightEdge_impl() const
{
    return xt::arange<size_t>(m_nely + 1) * (m_nelx + 1) + m_nelx;
}

inline xt::xtensor<size_t, 1> Regular::nodesBottomOpenEdge_impl() const
{
    return xt::arange<size_t>(1, m_nelx);
}

inline xt::xtensor<size_t, 1> Regular::nodesTopOpenEdge_impl() const
{
    return xt::arange<size_t>(1, m_nelx) + m_nely * (m_nelx + 1);
}

inline xt::xtensor<size_t, 1> Regular::nodesLeftOpenEdge_impl() const
{
    return xt::arange<size_t>(1, m_nely) * (m_nelx + 1);
}

inline xt::xtensor<size_t, 1> Regular::nodesRightOpenEdge_impl() const
{
    return xt::arange<size_t>(1, m_nely) * (m_nelx + 1) + m_nelx;
}

inline size_t Regular::nodesBottomLeftCorner_impl() const
{
    return 0;
}

inline size_t Regular::nodesBottomRightCorner_impl() const
{
    return m_nelx;
}

inline size_t Regular::nodesTopLeftCorner_impl() const
{
    return m_nely * (m_nelx + 1);
}

inline size_t Regular::nodesTopRightCorner_impl() const
{
    return m_nely * (m_nelx + 1) + m_nelx;
}

inline xt::xtensor<size_t, 2> Regular::elementgrid() const
{
    return xt::arange<size_t>(m_nelem).reshape({m_nely, m_nelx});
}

inline FineLayer::FineLayer(size_t nelx, size_t nely, double h, size_t nfine)
{
    this->init(nelx, nely, h, nfine);
}

inline FineLayer::FineLayer(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn)
{
    this->map(coor, conn);
}

inline void FineLayer::init(size_t nelx, size_t nely, double h, size_t nfine)
{
    GOOSEFEM_ASSERT(nelx >= 1ul);
    GOOSEFEM_ASSERT(nely >= 1ul);

    m_h = h;
    m_ndim = 2;
    m_nne = 4;
    m_Lx = m_h * static_cast<double>(nelx);

    // compute element size in y-direction (use symmetry, compute upper half)

    // temporary variables
    size_t nmin, ntot;
    xt::xtensor<size_t, 1> nhx = xt::ones<size_t>({nely});
    xt::xtensor<size_t, 1> nhy = xt::ones<size_t>({nely});
    xt::xtensor<int, 1> refine = -1 * xt::ones<int>({nely});

    // minimum height in y-direction (half of the height because of symmetry)
    if (nely % 2 == 0) {
        nmin = nely / 2;
    }
    else {
        nmin = (nely + 1) / 2;
    }

    // minimum number of fine layers in y-direction (minimum 1, middle layer part of this half)
    if (nfine % 2 == 0) {
        nfine = nfine / 2 + 1;
    }
    else {
        nfine = (nfine + 1) / 2;
    }

    if (nfine < 1) {
        nfine = 1;
    }

    if (nfine > nmin) {
        nfine = nmin;
    }

    // loop over element layers in y-direction, try to coarsen using these rules:
    // (1) element size in y-direction <= distance to origin in y-direction
    // (2) element size in x-direction should fit the total number of elements in x-direction
    // (3) a certain number of layers have the minimum size "1" (are fine)
    for (size_t iy = nfine;;) {
        // initialize current size in y-direction
        if (iy == nfine) {
            ntot = nfine;
        }
        // check to stop
        if (iy >= nely || ntot >= nmin) {
            nely = iy;
            break;
        }

        // rules (1,2) satisfied: coarsen in x-direction
        if (3 * nhy(iy) <= ntot && nelx % (3 * nhx(iy)) == 0 && ntot + nhy(iy) < nmin) {
            refine(iy) = 0;
            nhy(iy) *= 2;
            auto vnhy = xt::view(nhy, xt::range(iy + 1, _));
            auto vnhx = xt::view(nhx, xt::range(iy, _));
            vnhy *= 3;
            vnhx *= 3;
        }

        // update the number of elements in y-direction
        ntot += nhy(iy);
        // proceed to next element layer in y-direction
        ++iy;
        // check to stop
        if (iy >= nely || ntot >= nmin) {
            nely = iy;
            break;
        }
    }

    // symmetrize, compute full information

    // allocate mesh constructor parameters
    m_nhx = xt::empty<size_t>({nely * 2 - 1});
    m_nhy = xt::empty<size_t>({nely * 2 - 1});
    m_refine = xt::empty<int>({nely * 2 - 1});
    m_layer_nelx = xt::empty<size_t>({nely * 2 - 1});
    m_nnd = xt::empty<size_t>({nely * 2});
    m_startElem = xt::empty<size_t>({nely * 2 - 1});
    m_startNode = xt::empty<size_t>({nely * 2});

    // fill
    // - lower half
    for (size_t iy = 0; iy < nely; ++iy) {
        m_nhx(iy) = nhx(nely - iy - 1);
        m_nhy(iy) = nhy(nely - iy - 1);
        m_refine(iy) = refine(nely - iy - 1);
    }
    // - upper half
    for (size_t iy = 0; iy < nely - 1; ++iy) {
        m_nhx(iy + nely) = nhx(iy + 1);
        m_nhy(iy + nely) = nhy(iy + 1);
        m_refine(iy + nely) = refine(iy + 1);
    }

    // update size
    nely = m_nhx.size();

    // compute the number of elements per element layer in y-direction
    for (size_t iy = 0; iy < nely; ++iy) {
        m_layer_nelx(iy) = nelx / m_nhx(iy);
    }

    // compute the number of nodes per node layer in y-direction
    for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
        m_nnd(iy) = m_layer_nelx(iy) + 1;
    }
    for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
        m_nnd(iy + 1) = m_layer_nelx(iy) + 1;
    }

    // compute mesh dimensions

    // initialize
    m_nnode = 0;
    m_nelem = 0;
    m_startNode(0) = 0;

    // loop over element layers (bottom -> middle, elements become finer)
    for (size_t i = 0; i < (nely - 1) / 2; ++i) {
        // - store the first element of the layer
        m_startElem(i) = m_nelem;
        // - add the nodes of this layer
        if (m_refine(i) == 0) {
            m_nnode += (3 * m_layer_nelx(i) + 1);
        }
        else {
            m_nnode += (m_layer_nelx(i) + 1);
        }
        // - add the elements of this layer
        if (m_refine(i) == 0) {
            m_nelem += (4 * m_layer_nelx(i));
        }
        else {
            m_nelem += (m_layer_nelx(i));
        }
        // - store the starting node of the next layer
        m_startNode(i + 1) = m_nnode;
    }

    // loop over element layers (middle -> top, elements become coarser)
    for (size_t i = (nely - 1) / 2; i < nely; ++i) {
        // - store the first element of the layer
        m_startElem(i) = m_nelem;
        // - add the nodes of this layer
        if (m_refine(i) == 0) {
            m_nnode += (5 * m_layer_nelx(i) + 1);
        }
        else {
            m_nnode += (m_layer_nelx(i) + 1);
        }
        // - add the elements of this layer
        if (m_refine(i) == 0) {
            m_nelem += (4 * m_layer_nelx(i));
        }
        else {
            m_nelem += (m_layer_nelx(i));
        }
        // - store the starting node of the next layer
        m_startNode(i + 1) = m_nnode;
    }
    // - add the top row of nodes
    m_nnode += m_layer_nelx(nely - 1) + 1;
}

inline size_t FineLayer::nelx_impl() const
{
    return xt::amax(m_layer_nelx)();
}

inline size_t FineLayer::nely_impl() const
{
    return xt::sum(m_nhy)();
}

inline xt::xtensor<size_t, 1> FineLayer::elemrow_nhx() const
{
    return m_nhx;
}

inline xt::xtensor<size_t, 1> FineLayer::elemrow_nhy() const
{
    return m_nhy;
}

inline xt::xtensor<size_t, 1> FineLayer::elemrow_nelem() const
{
    return m_layer_nelx;
}

inline ElementType FineLayer::getElementType_impl() const
{
    return ElementType::Quad4;
}

inline xt::xtensor<double, 2> FineLayer::coor_impl() const
{
    // allocate output
    xt::xtensor<double, 2> ret = xt::empty<double>({m_nnode, m_ndim});

    // current node, number of element layers
    size_t inode = 0;
    size_t nely = static_cast<size_t>(m_nhy.size());

    // y-position of each main node layer (i.e. excluding node layers for refinement/coarsening)
    // - allocate
    xt::xtensor<double, 1> y = xt::empty<double>({nely + 1});
    // - initialize
    y(0) = 0.0;
    // - compute
    for (size_t iy = 1; iy < nely + 1; ++iy) {
        y(iy) = y(iy - 1) + m_nhy(iy - 1) * m_h;
    }

    // loop over element layers (bottom -> middle) : add bottom layer (+ refinement layer) of nodes

    for (size_t iy = 0;; ++iy) {
        // get positions along the x- and z-axis
        xt::xtensor<double, 1> x = xt::linspace<double>(0.0, m_Lx, m_layer_nelx(iy) + 1);

        // add nodes of the bottom layer of this element
        for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
            ret(inode, 0) = x(ix);
            ret(inode, 1) = y(iy);
            ++inode;
        }

        // stop at middle layer
        if (iy == (nely - 1) / 2) {
            break;
        }

        // add extra nodes of the intermediate layer, for refinement in x-direction
        if (m_refine(iy) == 0) {
            // - get position offset in x- and y-direction
            double dx = m_h * static_cast<double>(m_nhx(iy) / 3);
            double dy = m_h * static_cast<double>(m_nhy(iy) / 2);
            // - add nodes of the intermediate layer
            for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                for (size_t j = 0; j < 2; ++j) {
                    ret(inode, 0) = x(ix) + dx * static_cast<double>(j + 1);
                    ret(inode, 1) = y(iy) + dy;
                    ++inode;
                }
            }
        }
    }

    // loop over element layers (middle -> top) : add (refinement layer +) top layer of nodes

    for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
        // get positions along the x- and z-axis
        xt::xtensor<double, 1> x = xt::linspace<double>(0.0, m_Lx, m_layer_nelx(iy) + 1);

        // add extra nodes of the intermediate layer, for refinement in x-direction
        if (m_refine(iy) == 0) {
            // - get position offset in x- and y-direction
            double dx = m_h * static_cast<double>(m_nhx(iy) / 3);
            double dy = m_h * static_cast<double>(m_nhy(iy) / 2);
            // - add nodes of the intermediate layer
            for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                for (size_t j = 0; j < 2; ++j) {
                    ret(inode, 0) = x(ix) + dx * static_cast<double>(j + 1);
                    ret(inode, 1) = y(iy) + dy;
                    ++inode;
                }
            }
        }

        // add nodes of the top layer of this element
        for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
            ret(inode, 0) = x(ix);
            ret(inode, 1) = y(iy + 1);
            ++inode;
        }
    }

    return ret;
}

inline xt::xtensor<size_t, 2> FineLayer::conn_impl() const
{
    // allocate output
    xt::xtensor<size_t, 2> ret = xt::empty<size_t>({m_nelem, m_nne});

    // current element, number of element layers, starting nodes of each node layer
    size_t ielem = 0;
    size_t nely = static_cast<size_t>(m_nhy.size());
    size_t bot, mid, top;

    // loop over all element layers
    for (size_t iy = 0; iy < nely; ++iy) {
        // - get: starting nodes of bottom(, middle) and top layer
        bot = m_startNode(iy);
        mid = m_startNode(iy) + m_nnd(iy);
        top = m_startNode(iy + 1);

        // - define connectivity: no coarsening/refinement
        if (m_refine(iy) == -1) {
            for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                ret(ielem, 0) = bot + (ix);
                ret(ielem, 1) = bot + (ix + 1);
                ret(ielem, 2) = top + (ix + 1);
                ret(ielem, 3) = top + (ix);
                ielem++;
            }
        }

        // - define connectivity: refinement along the x-direction (below the middle layer)
        else if (m_refine(iy) == 0 && iy <= (nely - 1) / 2) {
            for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                // -- bottom element
                ret(ielem, 0) = bot + (ix);
                ret(ielem, 1) = bot + (ix + 1);
                ret(ielem, 2) = mid + (2 * ix + 1);
                ret(ielem, 3) = mid + (2 * ix);
                ielem++;
                // -- top-right element
                ret(ielem, 0) = bot + (ix + 1);
                ret(ielem, 1) = top + (3 * ix + 3);
                ret(ielem, 2) = top + (3 * ix + 2);
                ret(ielem, 3) = mid + (2 * ix + 1);
                ielem++;
                // -- top-center element
                ret(ielem, 0) = mid + (2 * ix);
                ret(ielem, 1) = mid + (2 * ix + 1);
                ret(ielem, 2) = top + (3 * ix + 2);
                ret(ielem, 3) = top + (3 * ix + 1);
                ielem++;
                // -- top-left element
                ret(ielem, 0) = bot + (ix);
                ret(ielem, 1) = mid + (2 * ix);
                ret(ielem, 2) = top + (3 * ix + 1);
                ret(ielem, 3) = top + (3 * ix);
                ielem++;
            }
        }

        // - define connectivity: coarsening along the x-direction (above the middle layer)
        else if (m_refine(iy) == 0 && iy > (nely - 1) / 2) {
            for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                // -- lower-left element
                ret(ielem, 0) = bot + (3 * ix);
                ret(ielem, 1) = bot + (3 * ix + 1);
                ret(ielem, 2) = mid + (2 * ix);
                ret(ielem, 3) = top + (ix);
                ielem++;
                // -- lower-center element
                ret(ielem, 0) = bot + (3 * ix + 1);
                ret(ielem, 1) = bot + (3 * ix + 2);
                ret(ielem, 2) = mid + (2 * ix + 1);
                ret(ielem, 3) = mid + (2 * ix);
                ielem++;
                // -- lower-right element
                ret(ielem, 0) = bot + (3 * ix + 2);
                ret(ielem, 1) = bot + (3 * ix + 3);
                ret(ielem, 2) = top + (ix + 1);
                ret(ielem, 3) = mid + (2 * ix + 1);
                ielem++;
                // -- upper element
                ret(ielem, 0) = mid + (2 * ix);
                ret(ielem, 1) = mid + (2 * ix + 1);
                ret(ielem, 2) = top + (ix + 1);
                ret(ielem, 3) = top + (ix);
                ielem++;
            }
        }
    }

    return ret;
}

inline xt::xtensor<size_t, 1> FineLayer::elementsMiddleLayer() const
{
    size_t nely = m_nhy.size();
    size_t iy = (nely - 1) / 2;
    return m_startElem(iy) + xt::arange<size_t>(m_layer_nelx(iy));
}

inline xt::xtensor<size_t, 1> FineLayer::elementsLayer(size_t iy) const
{
    GOOSEFEM_ASSERT(iy < m_layer_nelx.size());
    size_t n = m_layer_nelx(iy);
    if (m_refine(iy) != -1) {
        n *= 4;
    }
    return m_startElem(iy) + xt::arange<size_t>(n);
}

inline xt::xtensor<size_t, 1> FineLayer::elementgrid_ravel(
    std::vector<size_t> start_stop_rows,
    std::vector<size_t> start_stop_cols) const
{
    GOOSEFEM_ASSERT(start_stop_rows.size() == 0 || start_stop_rows.size() == 2);
    GOOSEFEM_ASSERT(start_stop_cols.size() == 0 || start_stop_cols.size() == 2);

    std::array<size_t, 2> rows;
    std::array<size_t, 2> cols;

    if (start_stop_rows.size() == 2) {
        std::copy(start_stop_rows.begin(), start_stop_rows.end(), rows.begin());
        GOOSEFEM_ASSERT(rows[0] <= this->nely());
        GOOSEFEM_ASSERT(rows[1] <= this->nely());
    }
    else {
        rows[0] = 0;
        rows[1] = this->nely();
    }

    if (start_stop_cols.size() == 2) {
        std::copy(start_stop_cols.begin(), start_stop_cols.end(), cols.begin());
        GOOSEFEM_ASSERT(cols[0] <= this->nelx());
        GOOSEFEM_ASSERT(cols[1] <= this->nelx());
    }
    else {
        cols[0] = 0;
        cols[1] = this->nelx();
    }

    if (rows[0] == rows[1] || cols[0] == cols[1]) {
        xt::xtensor<size_t, 1> ret = xt::empty<size_t>({0});
        return ret;
    }

    // Compute dimensions

    auto H = xt::cumsum(m_nhy);
    size_t yl = 0;
    if (rows[0] > 0) {
        yl = xt::argmax(H > rows[0])();
    }
    size_t yu = xt::argmax(H >= rows[1])();
    size_t hx = std::max(m_nhx(yl), m_nhx(yu));
    size_t xl = (cols[0] - cols[0] % hx) / hx;
    size_t xu = (cols[1] - cols[1] % hx) / hx;

    // Allocate output

    size_t N = 0;

    for (size_t i = yl; i <= yu; ++i) {
        // no refinement
        size_t n = (xu - xl) * hx / m_nhx(i);
        // refinement
        if (m_refine(i) != -1) {
            n *= 4;
        }
        N += n;
    }

    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({N});

    // Write output

    N = 0;

    for (size_t i = yl; i <= yu; ++i) {
        // no refinement
        size_t n = (xu - xl) * hx / m_nhx(i);
        size_t h = hx;
        // refinement
        if (m_refine(i) != -1) {
            n *= 4;
            h *= 4;
        }
        xt::xtensor<size_t, 1> e = m_startElem(i) + xl * h / m_nhx(i) + xt::arange<size_t>(n);
        xt::view(ret, xt::range(N, N + n)) = e;
        N += n;
    }

    return ret;
}

inline xt::xtensor<size_t, 1> FineLayer::elementgrid_around_ravel(
    size_t e,
    size_t size,
    bool periodic)
{
    GOOSEFEM_WIP_ASSERT(periodic == true);

    size_t iy = xt::argmin(m_startElem <= e)() - 1;
    size_t nel = m_layer_nelx(iy);

    GOOSEFEM_WIP_ASSERT(iy == (m_nhy.size() - 1) / 2);

    auto H = xt::cumsum(m_nhy);

    if (2 * size >= H(H.size() - 1)) {
        return xt::arange<size_t>(this->nelem());
    }

    size_t hy = H(iy);
    size_t l = xt::argmax(H > (hy - size - 1))();
    size_t u = xt::argmax(H >= (hy + size))();
    size_t lh = 0;
    if (l > 0) {
        lh = H(l - 1);
    }
    size_t uh = H(u);

    size_t step = xt::amax(m_nhx)();
    size_t relx = (e - m_startElem(iy)) % step;
    size_t mid = (nel / step - (nel / step) % 2) / 2 * step + relx;
    size_t nroll = (nel - (nel - mid + e - m_startElem(iy)) % nel) / step;
    size_t dx = m_nhx(u);
    size_t xl = 0;
    size_t xu = nel;
    if (mid >= size) {
        xl = mid - size;
    }
    if (mid + size < nel) {
        xu = mid + size + 1;
    }
    xl = xl - xl % dx;
    xu = xu - xu % dx;
    if (mid - xl < size) {
        if (xl < dx) {
            xl = 0;
        }
        else {
            xl -= dx;
        }
    }
    if (xu - mid < size) {
        if (xu > nel - dx) {
            xu = nel;
        }
        else {
            xu += dx;
        }
    }

    auto ret = this->elementgrid_ravel({lh, uh}, {xl, xu});
    auto map = this->roll(nroll);
    return xt::view(map, xt::keep(ret));
}

inline xt::xtensor<size_t, 1> FineLayer::elementgrid_leftright(
    size_t e,
    size_t left,
    size_t right,
    bool periodic)
{
    GOOSEFEM_WIP_ASSERT(periodic == true);

    size_t iy = xt::argmin(m_startElem <= e)() - 1;
    size_t nel = m_layer_nelx(iy);

    GOOSEFEM_WIP_ASSERT(iy == (m_nhy.size() - 1) / 2);

    size_t step = xt::amax(m_nhx)();
    size_t relx = (e - m_startElem(iy)) % step;
    size_t mid = (nel / step - (nel / step) % 2) / 2 * step + relx;
    size_t nroll = (nel - (nel - mid + e - m_startElem(iy)) % nel) / step;
    size_t dx = m_nhx(iy);
    size_t xl = 0;
    size_t xu = nel;
    if (mid >= left) {
        xl = mid - left;
    }
    if (mid + right < nel) {
        xu = mid + right + 1;
    }
    xl = xl - xl % dx;
    xu = xu - xu % dx;
    if (mid - xl < left) {
        if (xl < dx) {
            xl = 0;
        }
        else {
            xl -= dx;
        }
    }
    if (xu - mid < right) {
        if (xu > nel - dx) {
            xu = nel;
        }
        else {
            xu += dx;
        }
    }

    auto H = xt::cumsum(m_nhy);
    size_t yl = 0;
    if (iy > 0) {
        yl = H(iy - 1);
    }
    auto ret = this->elementgrid_ravel({yl, H(iy)}, {xl, xu});
    auto map = this->roll(nroll);
    return xt::view(map, xt::keep(ret));
}

inline xt::xtensor<size_t, 1> FineLayer::nodesBottomEdge_impl() const
{
    return m_startNode(0) + xt::arange<size_t>(m_layer_nelx(0) + 1);
}

inline xt::xtensor<size_t, 1> FineLayer::nodesTopEdge_impl() const
{
    size_t nely = m_nhy.size();
    return m_startNode(nely) + xt::arange<size_t>(m_layer_nelx(nely - 1) + 1);
}

inline xt::xtensor<size_t, 1> FineLayer::nodesLeftEdge_impl() const
{
    size_t nely = m_nhy.size();

    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({nely + 1});

    size_t i = 0;
    size_t j = (nely + 1) / 2;
    size_t k = (nely - 1) / 2;
    size_t l = nely;

    xt::view(ret, xt::range(i, j)) = xt::view(m_startNode, xt::range(i, j));
    xt::view(ret, xt::range(k + 1, l + 1)) = xt::view(m_startNode, xt::range(k + 1, l + 1));

    return ret;
}

inline xt::xtensor<size_t, 1> FineLayer::nodesRightEdge_impl() const
{
    size_t nely = m_nhy.size();

    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({nely + 1});

    size_t i = 0;
    size_t j = (nely + 1) / 2;
    size_t k = (nely - 1) / 2;
    size_t l = nely;

    xt::view(ret, xt::range(i, j)) =
        xt::view(m_startNode, xt::range(i, j)) + xt::view(m_layer_nelx, xt::range(i, j));

    xt::view(ret, xt::range(k + 1, l + 1)) =
        xt::view(m_startNode, xt::range(k + 1, l + 1)) + xt::view(m_layer_nelx, xt::range(k, l));

    return ret;
}

inline xt::xtensor<size_t, 1> FineLayer::nodesBottomOpenEdge_impl() const
{
    return m_startNode(0) + xt::arange<size_t>(1, m_layer_nelx(0));
}

inline xt::xtensor<size_t, 1> FineLayer::nodesTopOpenEdge_impl() const
{
    size_t nely = m_nhy.size();

    return m_startNode(nely) + xt::arange<size_t>(1, m_layer_nelx(nely - 1));
}

inline xt::xtensor<size_t, 1> FineLayer::nodesLeftOpenEdge_impl() const
{
    size_t nely = m_nhy.size();

    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({nely - 1});

    size_t i = 0;
    size_t j = (nely + 1) / 2;
    size_t k = (nely - 1) / 2;
    size_t l = nely;

    xt::view(ret, xt::range(i, j - 1)) = xt::view(m_startNode, xt::range(i + 1, j));
    xt::view(ret, xt::range(k, l - 1)) = xt::view(m_startNode, xt::range(k + 1, l));

    return ret;
}

inline xt::xtensor<size_t, 1> FineLayer::nodesRightOpenEdge_impl() const
{
    size_t nely = m_nhy.size();

    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({nely - 1});

    size_t i = 0;
    size_t j = (nely + 1) / 2;
    size_t k = (nely - 1) / 2;
    size_t l = nely;

    xt::view(ret, xt::range(i, j - 1)) =
        xt::view(m_startNode, xt::range(i + 1, j)) + xt::view(m_layer_nelx, xt::range(i + 1, j));

    xt::view(ret, xt::range(k, l - 1)) =
        xt::view(m_startNode, xt::range(k + 1, l)) + xt::view(m_layer_nelx, xt::range(k, l - 1));

    return ret;
}

inline size_t FineLayer::nodesBottomLeftCorner_impl() const
{
    return m_startNode(0);
}

inline size_t FineLayer::nodesBottomRightCorner_impl() const
{
    return m_startNode(0) + m_layer_nelx(0);
}

inline size_t FineLayer::nodesTopLeftCorner_impl() const
{
    size_t nely = m_nhy.size();

    return m_startNode(nely);
}

inline size_t FineLayer::nodesTopRightCorner_impl() const
{
    size_t nely = m_nhy.size();

    return m_startNode(nely) + m_layer_nelx(nely - 1);
}

inline xt::xtensor<size_t, 1> FineLayer::roll(size_t n)
{
    auto conn = this->conn();
    size_t nely = static_cast<size_t>(m_nhy.size());
    xt::xtensor<size_t, 1> ret = xt::empty<size_t>({m_nelem});

    // loop over all element layers
    for (size_t iy = 0; iy < nely; ++iy) {

        // no refinement
        size_t shift = n * (m_layer_nelx(iy) / m_layer_nelx(0));
        size_t nel = m_layer_nelx(iy);

        // refinement
        if (m_refine(iy) != -1) {
            shift = n * (m_layer_nelx(iy) / m_layer_nelx(0)) * 4;
            nel = m_layer_nelx(iy) * 4;
        }

        // element numbers of the layer, and roll them
        auto e = m_startElem(iy) + xt::arange<size_t>(nel);
        xt::view(ret, xt::range(m_startElem(iy), m_startElem(iy) + nel)) = xt::roll(e, shift);
    }

    return ret;
}

inline void FineLayer::map(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn)
{
    GOOSEFEM_ASSERT(coor.shape(1) == 2);
    GOOSEFEM_ASSERT(conn.shape(1) == 4);
    GOOSEFEM_ASSERT(conn.shape(0) > 0);
    GOOSEFEM_ASSERT(coor.shape(0) >= 4);
    GOOSEFEM_ASSERT(xt::amax(conn)() < coor.shape(0));

    if (conn.shape(0) == 1) {
        this->init(1, 1, coor(conn(0, 1), 0) - coor(conn(0, 0), 0));
        GOOSEFEM_CHECK(xt::all(xt::equal(this->conn(), conn)));
        GOOSEFEM_CHECK(xt::allclose(this->coor(), coor));
        return;
    }

    // Identify the middle layer

    size_t emid = (conn.shape(0) - conn.shape(0) % 2) / 2;
    size_t eleft = emid;
    size_t eright = emid;

    while (conn(eleft, 0) == conn(eleft - 1, 1) && eleft > 0) {
        eleft--;
    }

    while (conn(eright, 1) == conn(eright + 1, 0) && eright < conn.shape(0) - 1) {
        eright++;
    }

    GOOSEFEM_CHECK(xt::allclose(coor(conn(eleft, 0), 0), 0.0));

    // Get element sizes along the middle layer

    auto n0 = xt::view(conn, xt::range(eleft, eright + 1), 0);
    auto n1 = xt::view(conn, xt::range(eleft, eright + 1), 1);
    auto n2 = xt::view(conn, xt::range(eleft, eright + 1), 2);
    auto dx = xt::view(coor, xt::keep(n1), 0) - xt::view(coor, xt::keep(n0), 0);
    auto dy = xt::view(coor, xt::keep(n2), 1) - xt::view(coor, xt::keep(n1), 1);
    auto hx = xt::amin(dx)();
    auto hy = xt::amin(dy)();

    GOOSEFEM_CHECK(xt::allclose(hx, hy));
    GOOSEFEM_CHECK(xt::allclose(dx, hx));
    GOOSEFEM_CHECK(xt::allclose(dy, hy));

    // Extract shape and initialise

    size_t nelx = eright - eleft + 1;
    size_t nely = static_cast<size_t>((coor(coor.shape(0) - 1, 1) - coor(0, 1)) / hx);
    this->init(nelx, nely, hx);
    GOOSEFEM_CHECK(xt::all(xt::equal(this->conn(), conn)));
    GOOSEFEM_CHECK(xt::allclose(this->coor(), coor));
    GOOSEFEM_CHECK(xt::all(xt::equal(this->elementsMiddleLayer(), eleft + xt::arange<size_t>(nelx))));
}

namespace Map {

inline RefineRegular::RefineRegular(
    const GooseFEM::Mesh::Quad4::Regular& mesh, size_t nx, size_t ny)
    : m_coarse(mesh), m_nx(nx), m_ny(ny)
{
    m_fine = Regular(nx * m_coarse.nelx(), ny * m_coarse.nely(), m_coarse.h());

    xt::xtensor<size_t, 2> elmat_coarse = m_coarse.elementgrid();
    xt::xtensor<size_t, 2> elmat_fine = m_fine.elementgrid();

    m_coarse2fine = xt::empty<size_t>({m_coarse.nelem(), nx * ny});

    for (size_t i = 0; i < elmat_coarse.shape(0); ++i) {
        for (size_t j = 0; j < elmat_coarse.shape(1); ++j) {
            xt::view(m_coarse2fine, elmat_coarse(i, j), xt::all()) = xt::flatten(xt::view(
                elmat_fine, xt::range(i * ny, (i + 1) * ny), xt::range(j * nx, (j + 1) * nx)));
        }
    }
}

inline size_t RefineRegular::nx() const
{
    return m_nx;
}

inline size_t RefineRegular::ny() const
{
    return m_ny;
}

inline GooseFEM::Mesh::Quad4::Regular RefineRegular::getCoarseMesh() const
{
    return m_coarse;
}

inline GooseFEM::Mesh::Quad4::Regular RefineRegular::getFineMesh() const
{
    return m_fine;
}

inline xt::xtensor<size_t, 2> RefineRegular::getMap() const
{
    return m_coarse2fine;
}

template <class T, size_t rank>
inline xt::xtensor<T, rank> RefineRegular::meanToCoarse(const xt::xtensor<T, rank>& data) const
{
    GOOSEFEM_ASSERT(data.shape(0) == m_coarse2fine.size());

    std::array<size_t, rank> shape;
    std::copy(data.shape().cbegin(), data.shape().cend(), &shape[0]);
    shape[0] = m_coarse2fine.shape(0);

    xt::xtensor<T, rank> ret = xt::empty<T>(shape);

    for (size_t i = 0; i < m_coarse2fine.shape(0); ++i) {
        auto e = xt::view(m_coarse2fine, i, xt::all());
        auto d = xt::view(data, xt::keep(e));
        xt::view(ret, i) = xt::mean(d, 0);
    }

    return ret;
}

template <class T, size_t rank, class S>
inline xt::xtensor<T, rank> RefineRegular::averageToCoarse(
    const xt::xtensor<T, rank>& data,
    const xt::xtensor<S, rank>& weights) const
{
    GOOSEFEM_ASSERT(data.shape(0) == m_coarse2fine.size());

    std::array<size_t, rank> shape;
    std::copy(data.shape().cbegin(), data.shape().cend(), &shape[0]);
    shape[0] = m_coarse2fine.shape(0);

    xt::xtensor<T, rank> ret = xt::empty<T>(shape);

    for (size_t i = 0; i < m_coarse2fine.shape(0); ++i) {
        auto e = xt::view(m_coarse2fine, i, xt::all());
        xt::xtensor<T, rank> d = xt::view(data, xt::keep(e));
        xt::xtensor<T, rank> w = xt::view(weights, xt::keep(e));
        xt::view(ret, i) = xt::average(d, w, {0});
    }

    return ret;
}

template <class T, size_t rank>
inline xt::xtensor<T, rank> RefineRegular::mapToFine(const xt::xtensor<T, rank>& data) const
{
    GOOSEFEM_ASSERT(data.shape(0) == m_coarse2fine.shape(0));

    std::array<size_t, rank> shape;
    std::copy(data.shape().cbegin(), data.shape().cend(), &shape[0]);
    shape[0] = m_coarse2fine.size();

    xt::xtensor<T, rank> ret = xt::empty<T>(shape);

    for (size_t e = 0; e < m_coarse2fine.shape(0); ++e) {
        for (size_t i = 0; i < m_coarse2fine.shape(1); ++i) {
            xt::view(ret, m_coarse2fine(e, i)) = xt::view(data, e);
        }
    }

    return ret;
}

inline FineLayer2Regular::FineLayer2Regular(const GooseFEM::Mesh::Quad4::FineLayer& mesh)
    : m_finelayer(mesh)
{
    // ------------
    // Regular-mesh
    // ------------

    m_regular = GooseFEM::Mesh::Quad4::Regular(
        xt::amax(m_finelayer.m_layer_nelx)(), xt::sum(m_finelayer.m_nhy)(), m_finelayer.m_h);

    // -------
    // mapping
    // -------

    // allocate mapping
    m_elem_regular.resize(m_finelayer.m_nelem);
    m_frac_regular.resize(m_finelayer.m_nelem);

    // alias
    xt::xtensor<size_t, 1> nhx = m_finelayer.m_nhx;
    xt::xtensor<size_t, 1> nhy = m_finelayer.m_nhy;
    xt::xtensor<size_t, 1> nelx = m_finelayer.m_layer_nelx;
    xt::xtensor<size_t, 1> start = m_finelayer.m_startElem;

    // 'matrix' of element numbers of the Regular-mesh
    xt::xtensor<size_t, 2> elementgrid = m_regular.elementgrid();

    // cumulative number of element-rows of the Regular-mesh per layer of the FineLayer-mesh
    xt::xtensor<size_t, 1> cum_nhy =
        xt::concatenate(xt::xtuple(xt::zeros<size_t>({1}), xt::cumsum(nhy)));

    // number of element layers in y-direction of the FineLayer-mesh
    size_t nely = nhy.size();

    // loop over layers of the FineLayer-mesh
    for (size_t iy = 0; iy < nely; ++iy) {
        // element numbers of the Regular-mesh along this layer of the FineLayer-mesh
        auto el_new = xt::view(elementgrid, xt::range(cum_nhy(iy), cum_nhy(iy + 1)), xt::all());

        // no coarsening/refinement
        // ------------------------

        if (m_finelayer.m_refine(iy) == -1) {
            // element numbers of the FineLayer-mesh along this layer
            xt::xtensor<size_t, 1> el_old = start(iy) + xt::arange<size_t>(nelx(iy));

            // loop along this layer of the FineLayer-mesh
            for (size_t ix = 0; ix < nelx(iy); ++ix) {
                // get the element numbers of the Regular-mesh for this element of the
                // FineLayer-mesh
                auto block =
                    xt::view(el_new, xt::all(), xt::range(ix * nhx(iy), (ix + 1) * nhx(iy)));

                // write to mapping
                for (auto& i : block) {
                    m_elem_regular[el_old(ix)].push_back(i);
                    m_frac_regular[el_old(ix)].push_back(1.0);
                }
            }
        }

        // refinement along the x-direction (below the middle layer)
        // ---------------------------------------------------------

        else if (m_finelayer.m_refine(iy) == 0 && iy <= (nely - 1) / 2) {
            // element numbers of the FineLayer-mesh along this layer
            // rows: coarse block, columns element numbers per block
            xt::xtensor<size_t, 2> el_old =
                start(iy) + xt::arange<size_t>(nelx(iy) * 4ul).reshape({-1, 4});

            // loop along this layer of the FineLayer-mesh
            for (size_t ix = 0; ix < nelx(iy); ++ix) {
                // get the element numbers of the Regular-mesh for this block of the FineLayer-mesh
                auto block =
                    xt::view(el_new, xt::all(), xt::range(ix * nhx(iy), (ix + 1) * nhx(iy)));

                // bottom: wide-to-narrow
                {
                    for (size_t j = 0; j < nhy(iy) / 2; ++j) {
                        auto e = xt::view(block, j, xt::range(j, nhx(iy) - j));

                        m_elem_regular[el_old(ix, 0)].push_back(e(0));
                        m_frac_regular[el_old(ix, 0)].push_back(0.5);

                        for (size_t k = 1; k < e.size() - 1; ++k) {
                            m_elem_regular[el_old(ix, 0)].push_back(e(k));
                            m_frac_regular[el_old(ix, 0)].push_back(1.0);
                        }

                        m_elem_regular[el_old(ix, 0)].push_back(e(e.size() - 1));
                        m_frac_regular[el_old(ix, 0)].push_back(0.5);
                    }
                }

                // top: regular small element
                {
                    auto e = xt::view(
                        block,
                        xt::range(nhy(iy) / 2, nhy(iy)),
                        xt::range(1 * nhx(iy) / 3, 2 * nhx(iy) / 3));

                    for (auto& i : e) {
                        m_elem_regular[el_old(ix, 2)].push_back(i);
                        m_frac_regular[el_old(ix, 2)].push_back(1.0);
                    }
                }

                // left
                {
                    // left-bottom: narrow-to-wide
                    for (size_t j = 0; j < nhy(iy) / 2; ++j) {
                        auto e = xt::view(block, j, xt::range(0, j + 1));

                        for (size_t k = 0; k < e.size() - 1; ++k) {
                            m_elem_regular[el_old(ix, 3)].push_back(e(k));
                            m_frac_regular[el_old(ix, 3)].push_back(1.0);
                        }

                        m_elem_regular[el_old(ix, 3)].push_back(e(e.size() - 1));
                        m_frac_regular[el_old(ix, 3)].push_back(0.5);
                    }

                    // left-top: regular
                    {
                        auto e = xt::view(
                            block,
                            xt::range(nhy(iy) / 2, nhy(iy)),
                            xt::range(0 * nhx(iy) / 3, 1 * nhx(iy) / 3));

                        for (auto& i : e) {
                            m_elem_regular[el_old(ix, 3)].push_back(i);
                            m_frac_regular[el_old(ix, 3)].push_back(1.0);
                        }
                    }
                }

                // right
                {
                    // right-bottom: narrow-to-wide
                    for (size_t j = 0; j < nhy(iy) / 2; ++j) {
                        auto e = xt::view(block, j, xt::range(nhx(iy) - j - 1, nhx(iy)));

                        m_elem_regular[el_old(ix, 1)].push_back(e(0));
                        m_frac_regular[el_old(ix, 1)].push_back(0.5);

                        for (size_t k = 1; k < e.size(); ++k) {
                            m_elem_regular[el_old(ix, 1)].push_back(e(k));
                            m_frac_regular[el_old(ix, 1)].push_back(1.0);
                        }
                    }

                    // right-top: regular
                    {
                        auto e = xt::view(
                            block,
                            xt::range(nhy(iy) / 2, nhy(iy)),
                            xt::range(2 * nhx(iy) / 3, 3 * nhx(iy) / 3));

                        for (auto& i : e) {
                            m_elem_regular[el_old(ix, 1)].push_back(i);
                            m_frac_regular[el_old(ix, 1)].push_back(1.0);
                        }
                    }
                }
            }
        }

        // coarsening along the x-direction (above the middle layer)
        else if (m_finelayer.m_refine(iy) == 0 && iy > (nely - 1) / 2) {
            // element numbers of the FineLayer-mesh along this layer
            // rows: coarse block, columns element numbers per block
            xt::xtensor<size_t, 2> el_old =
                start(iy) + xt::arange<size_t>(nelx(iy) * 4ul).reshape({-1, 4});

            // loop along this layer of the FineLayer-mesh
            for (size_t ix = 0; ix < nelx(iy); ++ix) {
                // get the element numbers of the Regular-mesh for this block of the FineLayer-mesh
                auto block =
                    xt::view(el_new, xt::all(), xt::range(ix * nhx(iy), (ix + 1) * nhx(iy)));

                // top: narrow-to-wide
                {
                    for (size_t j = 0; j < nhy(iy) / 2; ++j) {
                        auto e = xt::view(
                            block,
                            nhy(iy) / 2 + j,
                            xt::range(1 * nhx(iy) / 3 - j - 1, 2 * nhx(iy) / 3 + j + 1));

                        m_elem_regular[el_old(ix, 3)].push_back(e(0));
                        m_frac_regular[el_old(ix, 3)].push_back(0.5);

                        for (size_t k = 1; k < e.size() - 1; ++k) {
                            m_elem_regular[el_old(ix, 3)].push_back(e(k));
                            m_frac_regular[el_old(ix, 3)].push_back(1.0);
                        }

                        m_elem_regular[el_old(ix, 3)].push_back(e(e.size() - 1));
                        m_frac_regular[el_old(ix, 3)].push_back(0.5);
                    }
                }

                // bottom: regular small element
                {
                    auto e = xt::view(
                        block,
                        xt::range(0, nhy(iy) / 2),
                        xt::range(1 * nhx(iy) / 3, 2 * nhx(iy) / 3));

                    for (auto& i : e) {
                        m_elem_regular[el_old(ix, 1)].push_back(i);
                        m_frac_regular[el_old(ix, 1)].push_back(1.0);
                    }
                }

                // left
                {
                    // left-bottom: regular
                    {
                        auto e = xt::view(
                            block,
                            xt::range(0, nhy(iy) / 2),
                            xt::range(0 * nhx(iy) / 3, 1 * nhx(iy) / 3));

                        for (auto& i : e) {
                            m_elem_regular[el_old(ix, 0)].push_back(i);
                            m_frac_regular[el_old(ix, 0)].push_back(1.0);
                        }
                    }

                    // left-top: narrow-to-wide
                    for (size_t j = 0; j < nhy(iy) / 2; ++j) {
                        auto e =
                            xt::view(block, nhy(iy) / 2 + j, xt::range(0, 1 * nhx(iy) / 3 - j));

                        for (size_t k = 0; k < e.size() - 1; ++k) {
                            m_elem_regular[el_old(ix, 0)].push_back(e(k));
                            m_frac_regular[el_old(ix, 0)].push_back(1.0);
                        }

                        m_elem_regular[el_old(ix, 0)].push_back(e(e.size() - 1));
                        m_frac_regular[el_old(ix, 0)].push_back(0.5);
                    }
                }

                // right
                {
                    // right-bottom: regular
                    {
                        auto e = xt::view(
                            block,
                            xt::range(0, nhy(iy) / 2),
                            xt::range(2 * nhx(iy) / 3, 3 * nhx(iy) / 3));

                        for (auto& i : e) {
                            m_elem_regular[el_old(ix, 2)].push_back(i);
                            m_frac_regular[el_old(ix, 2)].push_back(1.0);
                        }
                    }

                    // right-top: narrow-to-wide
                    for (size_t j = 0; j < nhy(iy) / 2; ++j) {
                        auto e = xt::view(
                            block, nhy(iy) / 2 + j, xt::range(2 * nhx(iy) / 3 + j, nhx(iy)));

                        m_elem_regular[el_old(ix, 2)].push_back(e(0));
                        m_frac_regular[el_old(ix, 2)].push_back(0.5);

                        for (size_t k = 1; k < e.size(); ++k) {
                            m_elem_regular[el_old(ix, 2)].push_back(e(k));
                            m_frac_regular[el_old(ix, 2)].push_back(1.0);
                        }
                    }
                }
            }
        }
    }
}

inline GooseFEM::Mesh::Quad4::Regular FineLayer2Regular::getRegularMesh() const
{
    return m_regular;
}

inline GooseFEM::Mesh::Quad4::FineLayer FineLayer2Regular::getFineLayerMesh() const
{
    return m_finelayer;
}

inline std::vector<std::vector<size_t>> FineLayer2Regular::getMap() const
{
    return m_elem_regular;
}

inline std::vector<std::vector<double>> FineLayer2Regular::getMapFraction() const
{
    return m_frac_regular;
}

template <class T, size_t rank>
inline xt::xtensor<T, rank> FineLayer2Regular::mapToRegular(const xt::xtensor<T, rank>& data) const
{
    GOOSEFEM_ASSERT(data.shape(0) == m_finelayer.nelem());

    std::array<size_t, rank> shape;
    std::copy(data.shape().cbegin(), data.shape().cend(), &shape[0]);
    shape[0] = m_regular.nelem();

    xt::xtensor<T, rank> ret = xt::zeros<T>(shape);

    for (size_t e = 0; e < m_finelayer.nelem(); ++e) {
        for (size_t i = 0; i < m_elem_regular[e].size(); ++i) {
            xt::view(ret, m_elem_regular[e][i]) += m_frac_regular[e][i] * xt::view(data, e);
        }
    }

    return ret;
}

} // namespace Map

} // namespace Quad4
} // namespace Mesh
} // namespace GooseFEM

#endif

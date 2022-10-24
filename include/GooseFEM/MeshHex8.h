/**
 * Generate simple meshes of 8-noded hexahedral elements in 3d (GooseFEM::Mesh::ElementType::Hex8).
 *
 * @file MeshHex8.h
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#ifndef GOOSEFEM_MESHHEX8_H
#define GOOSEFEM_MESHHEX8_H

#include "Mesh.h"
#include "config.h"

namespace GooseFEM {
namespace Mesh {

/**
 * Simple meshes of 8-noded hexahedral elements in 3d (ElementType::Hex8).
 */
namespace Hex8 {

/**
 * Regular mesh: equi-sized elements.
 */
class Regular : public RegularBase3d<Regular> {
public:
    /**
     * Constructor.
     *
     * @param nelx Number of elements in horizontal (x) direction.
     * @param nely Number of elements in vertical (y) direction.
     * @param nelz Number of elements in vertical (z) direction.
     * @param h Edge size (width == height == depth).
     */
    Regular(size_t nelx, size_t nely, size_t nelz, double h = 1.0)
    {
        m_h = h;
        m_nelx = nelx;
        m_nely = nely;
        m_nelz = nelz;
        m_nne = 8;
        m_ndim = 3;

        GOOSEFEM_ASSERT(m_nelx >= 1ul);
        GOOSEFEM_ASSERT(m_nely >= 1ul);
        GOOSEFEM_ASSERT(m_nelz >= 1ul);

        m_nnode = (m_nelx + 1) * (m_nely + 1) * (m_nelz + 1);
        m_nelem = m_nelx * m_nely * m_nelz;
    }

private:
    friend class RegularBase<Regular>;
    friend class RegularBase3d<Regular>;

    size_t nelx_impl() const
    {
        return m_nelx;
    }

    size_t nely_impl() const
    {
        return m_nely;
    }

    size_t nelz_impl() const
    {
        return m_nelz;
    }

    ElementType getElementType_impl() const
    {
        return ElementType::Hex8;
    }

    array_type::tensor<double, 2> coor_impl() const
    {
        array_type::tensor<double, 2> ret = xt::empty<double>({m_nnode, m_ndim});

        array_type::tensor<double, 1> x =
            xt::linspace<double>(0.0, m_h * static_cast<double>(m_nelx), m_nelx + 1);
        array_type::tensor<double, 1> y =
            xt::linspace<double>(0.0, m_h * static_cast<double>(m_nely), m_nely + 1);
        array_type::tensor<double, 1> z =
            xt::linspace<double>(0.0, m_h * static_cast<double>(m_nelz), m_nelz + 1);

        size_t inode = 0;

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            for (size_t iy = 0; iy < m_nely + 1; ++iy) {
                for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
                    ret(inode, 0) = x(ix);
                    ret(inode, 1) = y(iy);
                    ret(inode, 2) = z(iz);
                    ++inode;
                }
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 2> conn_impl() const
    {
        array_type::tensor<size_t, 2> ret = xt::empty<size_t>({m_nelem, m_nne});

        size_t ielem = 0;

        for (size_t iz = 0; iz < m_nelz; ++iz) {
            for (size_t iy = 0; iy < m_nely; ++iy) {
                for (size_t ix = 0; ix < m_nelx; ++ix) {
                    ret(ielem, 0) = iy * (m_nelx + 1) + ix + iz * (m_nely + 1) * (m_nelx + 1);
                    ret(ielem, 1) = iy * (m_nelx + 1) + (ix + 1) + iz * (m_nely + 1) * (m_nelx + 1);
                    ret(ielem, 3) = (iy + 1) * (m_nelx + 1) + ix + iz * (m_nely + 1) * (m_nelx + 1);
                    ret(ielem, 2) =
                        (iy + 1) * (m_nelx + 1) + (ix + 1) + iz * (m_nely + 1) * (m_nelx + 1);
                    ret(ielem, 4) = iy * (m_nelx + 1) + ix + (iz + 1) * (m_nely + 1) * (m_nelx + 1);
                    ret(ielem, 5) =
                        (iy) * (m_nelx + 1) + (ix + 1) + (iz + 1) * (m_nely + 1) * (m_nelx + 1);
                    ret(ielem, 7) =
                        (iy + 1) * (m_nelx + 1) + ix + (iz + 1) * (m_nely + 1) * (m_nelx + 1);
                    ret(ielem, 6) =
                        (iy + 1) * (m_nelx + 1) + (ix + 1) + (iz + 1) * (m_nely + 1) * (m_nelx + 1);
                    ++ielem;
                }
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFront_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nelx + 1) * (m_nely + 1)});

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
                ret(iy * (m_nelx + 1) + ix) = iy * (m_nelx + 1) + ix;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBack_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nelx + 1) * (m_nely + 1)});

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
                ret(iy * (m_nelx + 1) + ix) =
                    iy * (m_nelx + 1) + ix + m_nelz * (m_nely + 1) * (m_nelx + 1);
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesLeft_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nely + 1) * (m_nelz + 1)});

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            for (size_t iy = 0; iy < m_nely + 1; ++iy) {
                ret(iz * (m_nely + 1) + iy) = iy * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1);
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesRight_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nely + 1) * (m_nelz + 1)});

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            for (size_t iy = 0; iy < m_nely + 1; ++iy) {
                ret(iz * (m_nely + 1) + iy) =
                    iy * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1) + m_nelx;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottom_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nelx + 1) * (m_nelz + 1)});

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
                ret(iz * (m_nelx + 1) + ix) = ix + iz * (m_nelx + 1) * (m_nely + 1);
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTop_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nelx + 1) * (m_nelz + 1)});

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
                ret(iz * (m_nelx + 1) + ix) =
                    ix + m_nely * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1);
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontFace_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nelx - 1) * (m_nely - 1)});

        for (size_t iy = 1; iy < m_nely; ++iy) {
            for (size_t ix = 1; ix < m_nelx; ++ix) {
                ret((iy - 1) * (m_nelx - 1) + (ix - 1)) = iy * (m_nelx + 1) + ix;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackFace_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nelx - 1) * (m_nely - 1)});

        for (size_t iy = 1; iy < m_nely; ++iy) {
            for (size_t ix = 1; ix < m_nelx; ++ix) {
                ret((iy - 1) * (m_nelx - 1) + (ix - 1)) =
                    iy * (m_nelx + 1) + ix + m_nelz * (m_nely + 1) * (m_nelx + 1);
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesLeftFace_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nely - 1) * (m_nelz - 1)});

        for (size_t iz = 1; iz < m_nelz; ++iz) {
            for (size_t iy = 1; iy < m_nely; ++iy) {
                ret((iz - 1) * (m_nely - 1) + (iy - 1)) =
                    iy * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1);
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesRightFace_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nely - 1) * (m_nelz - 1)});

        for (size_t iz = 1; iz < m_nelz; ++iz) {
            for (size_t iy = 1; iy < m_nely; ++iy) {
                ret((iz - 1) * (m_nely - 1) + (iy - 1)) =
                    iy * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1) + m_nelx;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomFace_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nelx - 1) * (m_nelz - 1)});

        for (size_t iz = 1; iz < m_nelz; ++iz) {
            for (size_t ix = 1; ix < m_nelx; ++ix) {
                ret((iz - 1) * (m_nelx - 1) + (ix - 1)) = ix + iz * (m_nelx + 1) * (m_nely + 1);
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopFace_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({(m_nelx - 1) * (m_nelz - 1)});

        for (size_t iz = 1; iz < m_nelz; ++iz) {
            for (size_t ix = 1; ix < m_nelx; ++ix) {
                ret((iz - 1) * (m_nelx - 1) + (ix - 1)) =
                    ix + m_nely * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1);
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontBottomEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx + 1});

        for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
            ret(ix) = ix;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontTopEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx + 1});

        for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
            ret(ix) = ix + m_nely * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontLeftEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely + 1});

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            ret(iy) = iy * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontRightEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely + 1});

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            ret(iy) = iy * (m_nelx + 1) + m_nelx;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackBottomEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx + 1});

        for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
            ret(ix) = ix + m_nelz * (m_nely + 1) * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackTopEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx + 1});

        for (size_t ix = 0; ix < m_nelx + 1; ++ix) {
            ret(ix) = m_nely * (m_nelx + 1) + ix + m_nelz * (m_nely + 1) * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackLeftEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely + 1});

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            ret(iy) = iy * (m_nelx + 1) + m_nelz * (m_nelx + 1) * (m_nely + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackRightEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely + 1});

        for (size_t iy = 0; iy < m_nely + 1; ++iy) {
            ret(iy) = iy * (m_nelx + 1) + m_nelz * (m_nelx + 1) * (m_nely + 1) + m_nelx;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomLeftEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelz + 1});

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            ret(iz) = iz * (m_nelx + 1) * (m_nely + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomRightEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelz + 1});

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            ret(iz) = iz * (m_nelx + 1) * (m_nely + 1) + m_nelx;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopLeftEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelz + 1});

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            ret(iz) = m_nely * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopRightEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelz + 1});

        for (size_t iz = 0; iz < m_nelz + 1; ++iz) {
            ret(iz) = m_nely * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1) + m_nelx;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontBottomOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx - 1});

        for (size_t ix = 1; ix < m_nelx; ++ix) {
            ret(ix - 1) = ix;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontTopOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx - 1});

        for (size_t ix = 1; ix < m_nelx; ++ix) {
            ret(ix - 1) = ix + m_nely * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontLeftOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely - 1});

        for (size_t iy = 1; iy < m_nely; ++iy) {
            ret(iy - 1) = iy * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontRightOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely - 1});

        for (size_t iy = 1; iy < m_nely; ++iy) {
            ret(iy - 1) = iy * (m_nelx + 1) + m_nelx;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackBottomOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx - 1});

        for (size_t ix = 1; ix < m_nelx; ++ix) {
            ret(ix - 1) = ix + m_nelz * (m_nely + 1) * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackTopOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelx - 1});

        for (size_t ix = 1; ix < m_nelx; ++ix) {
            ret(ix - 1) = m_nely * (m_nelx + 1) + ix + m_nelz * (m_nely + 1) * (m_nelx + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackLeftOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely - 1});

        for (size_t iy = 1; iy < m_nely; ++iy) {
            ret(iy - 1) = iy * (m_nelx + 1) + m_nelz * (m_nelx + 1) * (m_nely + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackRightOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nely - 1});

        for (size_t iy = 1; iy < m_nely; ++iy) {
            ret(iy - 1) = iy * (m_nelx + 1) + m_nelz * (m_nelx + 1) * (m_nely + 1) + m_nelx;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomLeftOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelz - 1});

        for (size_t iz = 1; iz < m_nelz; ++iz) {
            ret(iz - 1) = iz * (m_nelx + 1) * (m_nely + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomRightOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelz - 1});

        for (size_t iz = 1; iz < m_nelz; ++iz) {
            ret(iz - 1) = iz * (m_nelx + 1) * (m_nely + 1) + m_nelx;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopLeftOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelz - 1});

        for (size_t iz = 1; iz < m_nelz; ++iz) {
            ret(iz - 1) = m_nely * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopRightOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nelz - 1});

        for (size_t iz = 1; iz < m_nelz; ++iz) {
            ret(iz - 1) = m_nely * (m_nelx + 1) + iz * (m_nelx + 1) * (m_nely + 1) + m_nelx;
        }

        return ret;
    }

    size_t nodesFrontBottomLeftCorner_impl() const
    {
        return 0;
    }

    size_t nodesFrontBottomRightCorner_impl() const
    {
        return m_nelx;
    }

    size_t nodesFrontTopLeftCorner_impl() const
    {
        return m_nely * (m_nelx + 1);
    }

    size_t nodesFrontTopRightCorner_impl() const
    {
        return m_nely * (m_nelx + 1) + m_nelx;
    }

    size_t nodesBackBottomLeftCorner_impl() const
    {
        return m_nelz * (m_nely + 1) * (m_nelx + 1);
    }

    size_t nodesBackBottomRightCorner_impl() const
    {
        return m_nelx + m_nelz * (m_nely + 1) * (m_nelx + 1);
    }

    size_t nodesBackTopLeftCorner_impl() const
    {
        return m_nely * (m_nelx + 1) + m_nelz * (m_nely + 1) * (m_nelx + 1);
    }

    size_t nodesBackTopRightCorner_impl() const
    {
        return m_nely * (m_nelx + 1) + m_nelx + m_nelz * (m_nely + 1) * (m_nelx + 1);
    }

    double m_h; ///< See h()
    size_t m_nelx; ///< See nelx()
    size_t m_nely; ///< See nely()
    size_t m_nelz; ///< See nely()
    size_t m_nelem; ///< See nelem()
    size_t m_nnode; ///< See nnode()
    size_t m_nne; ///< See nne()
    size_t m_ndim; ///< See ndim()
};

/**
 * Mesh with fine middle layer, and coarser elements towards the top and bottom.
 */
class FineLayer : public RegularBase3d<FineLayer> {
public:
    /**
     * Constructor.
     *
     * @param nelx Number of elements (along the middle layer) in horizontal (x) direction.
     * @param nely Approximate equivalent number of elements in vertical (y) direction.
     * @param nelz Number of elements (along the middle layer) in depth (z) direction.
     * @param h Edge size (width == height == depth) of elements along the weak layer.
     *
     * @param nfine
     *      Extra number of fine layers around the middle layer.
     *      By default the element size is kept smaller than the distance to the middle layer.
     */
    FineLayer(size_t nelx, size_t nely, size_t nelz, double h = 1.0, size_t nfine = 1)
    {
        m_h = h;
        m_nne = 8;
        m_ndim = 3;

        // basic assumptions
        GOOSEFEM_ASSERT(nelx >= 1ul);
        GOOSEFEM_ASSERT(nely >= 1ul);
        GOOSEFEM_ASSERT(nelz >= 1ul);

        // store basic info
        m_Lx = m_h * static_cast<double>(nelx);
        m_Lz = m_h * static_cast<double>(nelz);

        // compute element size in y-direction (use symmetry, compute upper half)

        // temporary variables
        size_t nmin, ntot;
        array_type::tensor<size_t, 1> nhx = xt::ones<size_t>({nely});
        array_type::tensor<size_t, 1> nhy = xt::ones<size_t>({nely});
        array_type::tensor<size_t, 1> nhz = xt::ones<size_t>({nely});
        array_type::tensor<int, 1> refine = -1 * xt::ones<int>({nely});

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
        // (2) element size in x-(z-)direction should fit the total number of elements in
        // x-(z-)direction (3) a certain number of layers have the minimum size "1" (are fine)
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

            // rules (1,2) satisfied: coarsen in x-direction (and z-direction)
            if (3 * nhy(iy) <= ntot && nelx % (3 * nhx(iy)) == 0 && ntot + nhy(iy) < nmin) {
                // - process refinement in x-direction
                refine(iy) = 0;
                nhy(iy) *= 2;
                auto vnhy = xt::view(nhy, xt::range(iy + 1, _));
                auto vnhx = xt::view(nhx, xt::range(iy, _));
                vnhy *= 3;
                vnhx *= 3;

                // - rule (2) satisfied: coarsen next element layer in z-direction
                if (iy + 1 < nely && ntot + 2 * nhy(iy) < nmin) {
                    if (nelz % (3 * nhz(iy + 1)) == 0) {
                        // - update the number of elements in y-direction
                        ntot += nhy(iy);
                        // - proceed to next element layer in y-direction
                        ++iy;
                        // - process refinement in z-direction
                        refine(iy) = 2;
                        nhy(iy) = nhy(iy - 1);
                        auto vnhz = xt::view(nhz, xt::range(iy, _));
                        vnhz *= 3;
                    }
                }
            }

            // rules (1,2) satisfied: coarse in z-direction
            else if (3 * nhy(iy) <= ntot && nelz % (3 * nhz(iy)) == 0 && ntot + nhy(iy) < nmin) {
                // - process refinement in z-direction
                refine(iy) = 2;
                nhy(iy) *= 2;
                auto vnhy = xt::view(nhy, xt::range(iy + 1, _));
                auto vnhz = xt::view(nhz, xt::range(iy, _));
                vnhy *= 3;
                vnhz *= 3;
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
        m_nhz = xt::empty<size_t>({nely * 2 - 1});
        m_refine = xt::empty<int>({nely * 2 - 1});
        m_layer_nelx = xt::empty<size_t>({nely * 2 - 1});
        m_layer_nelz = xt::empty<size_t>({nely * 2 - 1});
        m_nnd = xt::empty<size_t>({nely * 2});
        m_startElem = xt::empty<size_t>({nely * 2 - 1});
        m_startNode = xt::empty<size_t>({nely * 2});

        // fill
        // - lower half
        for (size_t iy = 0; iy < nely; ++iy) {
            m_nhx(iy) = nhx(nely - iy - 1);
            m_nhy(iy) = nhy(nely - iy - 1);
            m_nhz(iy) = nhz(nely - iy - 1);
            m_refine(iy) = refine(nely - iy - 1);
        }
        // - upper half
        for (size_t iy = 0; iy < nely - 1; ++iy) {
            m_nhx(iy + nely) = nhx(iy + 1);
            m_nhy(iy + nely) = nhy(iy + 1);
            m_nhz(iy + nely) = nhz(iy + 1);
            m_refine(iy + nely) = refine(iy + 1);
        }

        // update size
        nely = m_nhx.size();

        // compute the number of elements per element layer in y-direction
        for (size_t iy = 0; iy < nely; ++iy) {
            m_layer_nelx(iy) = nelx / m_nhx(iy);
            m_layer_nelz(iy) = nelz / m_nhz(iy);
        }

        // compute the number of nodes per node layer in y-direction
        // - bottom half
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy)
            m_nnd(iy) = (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy) + 1);
        // - top half
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy)
            m_nnd(iy + 1) = (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy) + 1);

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
                m_nnode += (3 * m_layer_nelx(i) + 1) * (m_layer_nelz(i) + 1);
            }
            else if (m_refine(i) == 2) {
                m_nnode += (m_layer_nelx(i) + 1) * (3 * m_layer_nelz(i) + 1);
            }
            else {
                m_nnode += (m_layer_nelx(i) + 1) * (m_layer_nelz(i) + 1);
            }
            // - add the elements of this layer
            if (m_refine(i) == 0) {
                m_nelem += (4 * m_layer_nelx(i)) * (m_layer_nelz(i));
            }
            else if (m_refine(i) == 2) {
                m_nelem += (m_layer_nelx(i)) * (4 * m_layer_nelz(i));
            }
            else {
                m_nelem += (m_layer_nelx(i)) * (m_layer_nelz(i));
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
                m_nnode += (5 * m_layer_nelx(i) + 1) * (m_layer_nelz(i) + 1);
            }
            else if (m_refine(i) == 2) {
                m_nnode += (m_layer_nelx(i) + 1) * (5 * m_layer_nelz(i) + 1);
            }
            else {
                m_nnode += (m_layer_nelx(i) + 1) * (m_layer_nelz(i) + 1);
            }
            // - add the elements of this layer
            if (m_refine(i) == 0) {
                m_nelem += (4 * m_layer_nelx(i)) * (m_layer_nelz(i));
            }
            else if (m_refine(i) == 2) {
                m_nelem += (m_layer_nelx(i)) * (4 * m_layer_nelz(i));
            }
            else {
                m_nelem += (m_layer_nelx(i)) * (m_layer_nelz(i));
            }
            // - store the starting node of the next layer
            m_startNode(i + 1) = m_nnode;
        }
        // - add the top row of nodes
        m_nnode += (m_layer_nelx(nely - 1) + 1) * (m_layer_nelz(nely - 1) + 1);
    }

    /**
     * Elements in the middle (fine) layer.
     * @return List of element numbers (copy, involves computation).
     */
    array_type::tensor<size_t, 1> elementsMiddleLayer() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());
        size_t y = (nely - 1) / 2;

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(y) * m_layer_nelz(y)});

        for (size_t x = 0; x < m_layer_nelx(y); ++x) {
            for (size_t z = 0; z < m_layer_nelz(y); ++z) {
                ret(x + z * m_layer_nelx(y)) = m_startElem(y) + x + z * m_layer_nelx(y);
            }
        }

        return ret;
    }

private:
    friend class RegularBase<FineLayer>;
    friend class RegularBase3d<FineLayer>;

    size_t nelx_impl() const
    {
        return xt::amax(m_layer_nelx)();
    }

    size_t nely_impl() const
    {
        return xt::sum(m_nhy)();
    }

    size_t nelz_impl() const
    {
        return xt::amax(m_layer_nelz)();
    }
    ElementType getElementType_impl() const
    {
        return ElementType::Hex8;
    }

    array_type::tensor<double, 2> coor_impl() const
    {
        // allocate output
        array_type::tensor<double, 2> ret = xt::empty<double>({m_nnode, m_ndim});

        // current node, number of element layers
        size_t inode = 0;
        size_t nely = static_cast<size_t>(m_nhy.size());

        // y-position of each main node layer (i.e. excluding node layers for refinement/coarsening)
        // - allocate
        array_type::tensor<double, 1> y = xt::empty<double>({nely + 1});
        // - initialize
        y(0) = 0.0;
        // - compute
        for (size_t iy = 1; iy < nely + 1; ++iy) {
            y(iy) = y(iy - 1) + m_nhy(iy - 1) * m_h;
        }

        // loop over element layers (bottom -> middle) : add bottom layer (+ refinement layer) of
        // nodes

        for (size_t iy = 0;; ++iy) {
            // get positions along the x- and z-axis
            array_type::tensor<double, 1> x = xt::linspace<double>(0.0, m_Lx, m_layer_nelx(iy) + 1);
            array_type::tensor<double, 1> z = xt::linspace<double>(0.0, m_Lz, m_layer_nelz(iy) + 1);

            // add nodes of the bottom layer of this element
            for (size_t iz = 0; iz < m_layer_nelz(iy) + 1; ++iz) {
                for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
                    ret(inode, 0) = x(ix);
                    ret(inode, 1) = y(iy);
                    ret(inode, 2) = z(iz);
                    ++inode;
                }
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
                for (size_t iz = 0; iz < m_layer_nelz(iy) + 1; ++iz) {
                    for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                        for (size_t j = 0; j < 2; ++j) {
                            ret(inode, 0) = x(ix) + dx * static_cast<double>(j + 1);
                            ret(inode, 1) = y(iy) + dy;
                            ret(inode, 2) = z(iz);
                            ++inode;
                        }
                    }
                }
            }

            // add extra nodes of the intermediate layer, for refinement in z-direction
            else if (m_refine(iy) == 2) {
                // - get position offset in y- and z-direction
                double dz = m_h * static_cast<double>(m_nhz(iy) / 3);
                double dy = m_h * static_cast<double>(m_nhy(iy) / 2);
                // - add nodes of the intermediate layer
                for (size_t iz = 0; iz < m_layer_nelz(iy); ++iz) {
                    for (size_t j = 0; j < 2; ++j) {
                        for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
                            ret(inode, 0) = x(ix);
                            ret(inode, 1) = y(iy) + dy;
                            ret(inode, 2) = z(iz) + dz * static_cast<double>(j + 1);
                            ++inode;
                        }
                    }
                }
            }
        }

        // loop over element layers (middle -> top) : add (refinement layer +) top layer of nodes

        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            // get positions along the x- and z-axis
            array_type::tensor<double, 1> x = xt::linspace<double>(0.0, m_Lx, m_layer_nelx(iy) + 1);
            array_type::tensor<double, 1> z = xt::linspace<double>(0.0, m_Lz, m_layer_nelz(iy) + 1);

            // add extra nodes of the intermediate layer, for refinement in x-direction
            if (m_refine(iy) == 0) {
                // - get position offset in x- and y-direction
                double dx = m_h * static_cast<double>(m_nhx(iy) / 3);
                double dy = m_h * static_cast<double>(m_nhy(iy) / 2);
                // - add nodes of the intermediate layer
                for (size_t iz = 0; iz < m_layer_nelz(iy) + 1; ++iz) {
                    for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                        for (size_t j = 0; j < 2; ++j) {
                            ret(inode, 0) = x(ix) + dx * static_cast<double>(j + 1);
                            ret(inode, 1) = y(iy) + dy;
                            ret(inode, 2) = z(iz);
                            ++inode;
                        }
                    }
                }
            }

            // add extra nodes of the intermediate layer, for refinement in z-direction
            else if (m_refine(iy) == 2) {
                // - get position offset in y- and z-direction
                double dz = m_h * static_cast<double>(m_nhz(iy) / 3);
                double dy = m_h * static_cast<double>(m_nhy(iy) / 2);
                // - add nodes of the intermediate layer
                for (size_t iz = 0; iz < m_layer_nelz(iy); ++iz) {
                    for (size_t j = 0; j < 2; ++j) {
                        for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
                            ret(inode, 0) = x(ix);
                            ret(inode, 1) = y(iy) + dy;
                            ret(inode, 2) = z(iz) + dz * static_cast<double>(j + 1);
                            ++inode;
                        }
                    }
                }
            }

            // add nodes of the top layer of this element
            for (size_t iz = 0; iz < m_layer_nelz(iy) + 1; ++iz) {
                for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
                    ret(inode, 0) = x(ix);
                    ret(inode, 1) = y(iy + 1);
                    ret(inode, 2) = z(iz);
                    ++inode;
                }
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 2> conn_impl() const
    {
        // allocate output
        array_type::tensor<size_t, 2> ret = xt::empty<size_t>({m_nelem, m_nne});

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
                for (size_t iz = 0; iz < m_layer_nelz(iy); ++iz) {
                    for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                        ret(ielem, 0) = bot + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = top + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = bot + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = top + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ielem++;
                    }
                }
            }

            // - define connectivity: refinement along the x-direction (below the middle layer)
            else if (m_refine(iy) == 0 && iy <= (nely - 1) / 2) {
                for (size_t iz = 0; iz < m_layer_nelz(iy); ++iz) {
                    for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                        // -- bottom element
                        ret(ielem, 0) = bot + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = mid + (2 * ix + 1) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 3) = mid + (2 * ix) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 4) = bot + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = mid + (2 * ix + 1) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ret(ielem, 7) = mid + (2 * ix) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ielem++;
                        // -- top-right element
                        ret(ielem, 0) = bot + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = top + (3 * ix + 3) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = top + (3 * ix + 2) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = mid + (2 * ix + 1) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 4) = bot + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = top + (3 * ix + 3) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = top + (3 * ix + 2) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = mid + (2 * ix + 1) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ielem++;
                        // -- top-center element
                        ret(ielem, 0) = mid + (2 * ix) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 1) = mid + (2 * ix + 1) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 2) = top + (3 * ix + 2) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + (3 * ix + 1) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = mid + (2 * ix) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ret(ielem, 5) = mid + (2 * ix + 1) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ret(ielem, 6) = top + (3 * ix + 2) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + (3 * ix + 1) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ielem++;
                        // -- top-left element
                        ret(ielem, 0) = bot + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = mid + (2 * ix) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 2) = top + (3 * ix + 1) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + (3 * ix) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = bot + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = mid + (2 * ix) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ret(ielem, 6) = top + (3 * ix + 1) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + (3 * ix) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ielem++;
                    }
                }
            }

            // - define connectivity: coarsening along the x-direction (above the middle layer)
            else if (m_refine(iy) == 0 && iy > (nely - 1) / 2) {
                for (size_t iz = 0; iz < m_layer_nelz(iy); ++iz) {
                    for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                        // -- lower-left element
                        ret(ielem, 0) = bot + (3 * ix) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (3 * ix + 1) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = mid + (2 * ix) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 3) = top + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = bot + (3 * ix) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (3 * ix + 1) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = mid + (2 * ix) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ret(ielem, 7) = top + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ielem++;
                        // -- lower-center element
                        ret(ielem, 0) = bot + (3 * ix + 1) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (3 * ix + 2) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = mid + (2 * ix + 1) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 3) = mid + (2 * ix) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 4) = bot + (3 * ix + 1) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (3 * ix + 2) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = mid + (2 * ix + 1) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ret(ielem, 7) = mid + (2 * ix) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ielem++;
                        // -- lower-right element
                        ret(ielem, 0) = bot + (3 * ix + 2) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (3 * ix + 3) + iz * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = top + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = mid + (2 * ix + 1) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 4) = bot + (3 * ix + 2) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (3 * ix + 3) + (iz + 1) * (3 * m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = top + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = mid + (2 * ix + 1) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ielem++;
                        // -- upper element
                        ret(ielem, 0) = mid + (2 * ix) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 1) = mid + (2 * ix + 1) + iz * (2 * m_layer_nelx(iy));
                        ret(ielem, 2) = top + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = mid + (2 * ix) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ret(ielem, 5) = mid + (2 * ix + 1) + (iz + 1) * (2 * m_layer_nelx(iy));
                        ret(ielem, 6) = top + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ielem++;
                    }
                }
            }

            // - define connectivity: refinement along the z-direction (below the middle layer)
            else if (m_refine(iy) == 2 && iy <= (nely - 1) / 2) {
                for (size_t iz = 0; iz < m_layer_nelz(iy); ++iz) {
                    for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                        // -- bottom element
                        ret(ielem, 0) = bot + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = bot + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = bot + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = mid + ix + 2 * iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = mid + ix + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = mid + (ix + 1) + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = mid + (ix + 1) + 2 * iz * (m_layer_nelx(iy) + 1);
                        ielem++;
                        // -- top-back element
                        ret(ielem, 0) = mid + ix + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = mid + (ix + 1) + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = top + (ix + 1) + (3 * iz + 2) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + ix + (3 * iz + 2) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = bot + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = top + (ix + 1) + (3 * iz + 3) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + ix + (3 * iz + 3) * (m_layer_nelx(iy) + 1);
                        ielem++;
                        // -- top-center element
                        ret(ielem, 0) = mid + ix + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = mid + (ix + 1) + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = top + (ix + 1) + (3 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + ix + (3 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = mid + ix + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = mid + (ix + 1) + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = top + (ix + 1) + (3 * iz + 2) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + ix + (3 * iz + 2) * (m_layer_nelx(iy) + 1);
                        ielem++;
                        // -- top-front element
                        ret(ielem, 0) = bot + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = top + (ix + 1) + (3 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + ix + (3 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = mid + ix + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = mid + (ix + 1) + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = top + (ix + 1) + (3 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + ix + (3 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ielem++;
                    }
                }
            }

            // - define connectivity: coarsening along the z-direction (above the middle layer)
            else if (m_refine(iy) == 2 && iy > (nely - 1) / 2) {
                for (size_t iz = 0; iz < m_layer_nelz(iy); ++iz) {
                    for (size_t ix = 0; ix < m_layer_nelx(iy); ++ix) {
                        // -- bottom-front element
                        ret(ielem, 0) = bot + ix + (3 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (ix + 1) + (3 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = top + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = bot + ix + (3 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (ix + 1) + (3 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = mid + (ix + 1) + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = mid + ix + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ielem++;
                        // -- bottom-center element
                        ret(ielem, 0) = bot + ix + (3 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (ix + 1) + (3 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = mid + (ix + 1) + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = mid + ix + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = bot + ix + (3 * iz + 2) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (ix + 1) + (3 * iz + 2) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = mid + (ix + 1) + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = mid + ix + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ielem++;
                        // -- bottom-back element
                        ret(ielem, 0) = bot + ix + (3 * iz + 2) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = bot + (ix + 1) + (3 * iz + 2) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = mid + (ix + 1) + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = mid + ix + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = bot + ix + (3 * iz + 3) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = bot + (ix + 1) + (3 * iz + 3) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = top + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ielem++;
                        // -- top element
                        ret(ielem, 0) = mid + ix + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 1) = mid + (ix + 1) + (2 * iz) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 2) = top + (ix + 1) + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 3) = top + ix + iz * (m_layer_nelx(iy) + 1);
                        ret(ielem, 4) = mid + ix + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 5) = mid + (ix + 1) + (2 * iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 6) = top + (ix + 1) + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ret(ielem, 7) = top + ix + (iz + 1) * (m_layer_nelx(iy) + 1);
                        ielem++;
                    }
                }
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFront_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // number of boundary nodes
        // - initialize
        size_t n = 0;
        // - bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            if (m_refine(iy) == 0) {
                n += m_layer_nelx(iy) * 3 + 1;
            }
            else {
                n += m_layer_nelx(iy) + 1;
            }
        }
        // - top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            if (m_refine(iy) == 0) {
                n += m_layer_nelx(iy) * 3 + 1;
            }
            else {
                n += m_layer_nelx(iy) + 1;
            }
        }

        // allocate node-list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        // initialize counter: current index in the node-list "ret"
        size_t j = 0;

        // bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            // -- bottom node layer
            for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
                ret(j) = m_startNode(iy) + ix;
                ++j;
            }
            // -- refinement layer
            if (m_refine(iy) == 0) {
                for (size_t ix = 0; ix < 2 * m_layer_nelx(iy); ++ix) {
                    ret(j) = m_startNode(iy) + ix + m_nnd(iy);
                    ++j;
                }
            }
        }

        // top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            // -- refinement layer
            if (m_refine(iy) == 0) {
                for (size_t ix = 0; ix < 2 * m_layer_nelx(iy); ++ix) {
                    ret(j) = m_startNode(iy) + ix + m_nnd(iy);
                    ++j;
                }
            }
            // -- top node layer
            for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
                ret(j) = m_startNode(iy + 1) + ix;
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBack_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // number of boundary nodes
        // - initialize
        size_t n = 0;
        // - bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            if (m_refine(iy) == 0) {
                n += m_layer_nelx(iy) * 3 + 1;
            }
            else {
                n += m_layer_nelx(iy) + 1;
            }
        }
        // - top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            if (m_refine(iy) == 0) {
                n += m_layer_nelx(iy) * 3 + 1;
            }
            else {
                n += m_layer_nelx(iy) + 1;
            }
        }

        // allocate node-list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        // initialize counter: current index in the node-list "ret"
        size_t j = 0;

        // bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            // -- bottom node layer
            for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
                ret(j) = m_startNode(iy) + ix + (m_layer_nelx(iy) + 1) * m_layer_nelz(iy);
                ++j;
            }
            // -- refinement layer
            if (m_refine(iy) == 0) {
                for (size_t ix = 0; ix < 2 * m_layer_nelx(iy); ++ix) {
                    ret(j) =
                        m_startNode(iy) + ix + m_nnd(iy) + 2 * m_layer_nelx(iy) * m_layer_nelz(iy);
                    ++j;
                }
            }
        }

        // top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            // -- refinement layer
            if (m_refine(iy) == 0) {
                for (size_t ix = 0; ix < 2 * m_layer_nelx(iy); ++ix) {
                    ret(j) =
                        m_startNode(iy) + ix + m_nnd(iy) + 2 * m_layer_nelx(iy) * m_layer_nelz(iy);
                    ++j;
                }
            }
            // -- top node layer
            for (size_t ix = 0; ix < m_layer_nelx(iy) + 1; ++ix) {
                ret(j) = m_startNode(iy + 1) + ix + (m_layer_nelx(iy) + 1) * m_layer_nelz(iy);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesLeft_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // number of boundary nodes
        // - initialize
        size_t n = 0;
        // - bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            if (m_refine(iy) == 2) {
                n += m_layer_nelz(iy) * 3 + 1;
            }
            else {
                n += m_layer_nelz(iy) + 1;
            }
        }
        // - top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            if (m_refine(iy) == 2) {
                n += m_layer_nelz(iy) * 3 + 1;
            }
            else {
                n += m_layer_nelz(iy) + 1;
            }
        }

        // allocate node-list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        // initialize counter: current index in the node-list "ret"
        size_t j = 0;

        // bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            // -- bottom node layer
            for (size_t iz = 0; iz < m_layer_nelz(iy) + 1; ++iz) {
                ret(j) = m_startNode(iy) + iz * (m_layer_nelx(iy) + 1);
                ++j;
            }
            // -- refinement layer
            if (m_refine(iy) == 2) {
                for (size_t iz = 0; iz < 2 * m_layer_nelz(iy); ++iz) {
                    ret(j) = m_startNode(iy) + iz * (m_layer_nelx(iy) + 1) + m_nnd(iy);
                    ++j;
                }
            }
        }

        // top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            // -- refinement layer
            if (m_refine(iy) == 2) {
                for (size_t iz = 0; iz < 2 * m_layer_nelz(iy); ++iz) {
                    ret(j) = m_startNode(iy) + iz * (m_layer_nelx(iy) + 1) + m_nnd(iy);
                    ++j;
                }
            }
            // -- top node layer
            for (size_t iz = 0; iz < m_layer_nelz(iy) + 1; ++iz) {
                ret(j) = m_startNode(iy + 1) + iz * (m_layer_nelx(iy) + 1);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesRight_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // number of boundary nodes
        // - initialize
        size_t n = 0;
        // - bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            if (m_refine(iy) == 2)
                n += m_layer_nelz(iy) * 3 + 1;
            else
                n += m_layer_nelz(iy) + 1;
        }
        // - top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            if (m_refine(iy) == 2)
                n += m_layer_nelz(iy) * 3 + 1;
            else
                n += m_layer_nelz(iy) + 1;
        }

        // allocate node-list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        // initialize counter: current index in the node-list "ret"
        size_t j = 0;

        // bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            // -- bottom node layer
            for (size_t iz = 0; iz < m_layer_nelz(iy) + 1; ++iz) {
                ret(j) = m_startNode(iy) + iz * (m_layer_nelx(iy) + 1) + m_layer_nelx(iy);
                ++j;
            }
            // -- refinement layer
            if (m_refine(iy) == 2) {
                for (size_t iz = 0; iz < 2 * m_layer_nelz(iy); ++iz) {
                    ret(j) = m_startNode(iy) + m_nnd(iy) + iz * (m_layer_nelx(iy) + 1) +
                             m_layer_nelx(iy);
                    ++j;
                }
            }
        }

        // top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            // -- refinement layer
            if (m_refine(iy) == 2) {
                for (size_t iz = 0; iz < 2 * m_layer_nelz(iy); ++iz) {
                    ret(j) = m_startNode(iy) + m_nnd(iy) + iz * (m_layer_nelx(iy) + 1) +
                             m_layer_nelx(iy);
                    ++j;
                }
            }
            // -- top node layer
            for (size_t iz = 0; iz < m_layer_nelz(iy) + 1; ++iz) {
                ret(j) = m_startNode(iy + 1) + iz * (m_layer_nelx(iy) + 1) + m_layer_nelx(iy);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottom_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // allocate node list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nnd(nely)});

        // counter
        size_t j = 0;

        // fill node list
        for (size_t ix = 0; ix < m_layer_nelx(0) + 1; ++ix) {
            for (size_t iz = 0; iz < m_layer_nelz(0) + 1; ++iz) {
                ret(j) = m_startNode(0) + ix + iz * (m_layer_nelx(0) + 1);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTop_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // allocate node list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_nnd(nely)});

        // counter
        size_t j = 0;

        // fill node list
        for (size_t ix = 0; ix < m_layer_nelx(nely - 1) + 1; ++ix) {
            for (size_t iz = 0; iz < m_layer_nelz(nely - 1) + 1; ++iz) {
                ret(j) = m_startNode(nely) + ix + iz * (m_layer_nelx(nely - 1) + 1);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontFace_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // number of boundary nodes
        // - initialize
        size_t n = 0;
        // - bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            if (m_refine(iy) == 0) {
                n += m_layer_nelx(iy) * 3 - 1;
            }
            else {
                n += m_layer_nelx(iy) - 1;
            }
        }
        // - top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            if (m_refine(iy) == 0) {
                n += m_layer_nelx(iy) * 3 - 1;
            }
            else {
                n += m_layer_nelx(iy) - 1;
            }
        }

        // allocate node-list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        // initialize counter: current index in the node-list "ret"
        size_t j = 0;

        // bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            // -- bottom node layer
            for (size_t ix = 1; ix < m_layer_nelx(iy); ++ix) {
                ret(j) = m_startNode(iy) + ix;
                ++j;
            }
            // -- refinement layer
            if (m_refine(iy) == 0) {
                for (size_t ix = 0; ix < 2 * m_layer_nelx(iy); ++ix) {
                    ret(j) = m_startNode(iy) + ix + m_nnd(iy);
                    ++j;
                }
            }
        }

        // top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            // -- refinement layer
            if (m_refine(iy) == 0) {
                for (size_t ix = 0; ix < 2 * m_layer_nelx(iy); ++ix) {
                    ret(j) = m_startNode(iy) + ix + m_nnd(iy);
                    ++j;
                }
            }
            // -- top node layer
            for (size_t ix = 1; ix < m_layer_nelx(iy); ++ix) {
                ret(j) = m_startNode(iy + 1) + ix;
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackFace_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // number of boundary nodes
        // - initialize
        size_t n = 0;
        // - bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            if (m_refine(iy) == 0) {
                n += m_layer_nelx(iy) * 3 - 1;
            }
            else {
                n += m_layer_nelx(iy) - 1;
            }
        }
        // - top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            if (m_refine(iy) == 0) {
                n += m_layer_nelx(iy) * 3 - 1;
            }
            else {
                n += m_layer_nelx(iy) - 1;
            }
        }

        // allocate node-list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        // initialize counter: current index in the node-list "ret"
        size_t j = 0;

        // bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            // -- bottom node layer
            for (size_t ix = 1; ix < m_layer_nelx(iy); ++ix) {
                ret(j) = m_startNode(iy) + ix + (m_layer_nelx(iy) + 1) * m_layer_nelz(iy);
                ++j;
            }
            // -- refinement layer
            if (m_refine(iy) == 0) {
                for (size_t ix = 0; ix < 2 * m_layer_nelx(iy); ++ix) {
                    ret(j) =
                        m_startNode(iy) + ix + m_nnd(iy) + 2 * m_layer_nelx(iy) * m_layer_nelz(iy);
                    ++j;
                }
            }
        }

        // top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            // -- refinement layer
            if (m_refine(iy) == 0) {
                for (size_t ix = 0; ix < 2 * m_layer_nelx(iy); ++ix) {
                    ret(j) =
                        m_startNode(iy) + ix + m_nnd(iy) + 2 * m_layer_nelx(iy) * m_layer_nelz(iy);
                    ++j;
                }
            }
            // -- top node layer
            for (size_t ix = 1; ix < m_layer_nelx(iy); ++ix) {
                ret(j) = m_startNode(iy + 1) + ix + (m_layer_nelx(iy) + 1) * m_layer_nelz(iy);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesLeftFace_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // number of boundary nodes
        // - initialize
        size_t n = 0;
        // - bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            if (m_refine(iy) == 2) {
                n += m_layer_nelz(iy) * 3 - 1;
            }
            else {
                n += m_layer_nelz(iy) - 1;
            }
        }
        // - top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            if (m_refine(iy) == 2) {
                n += m_layer_nelz(iy) * 3 - 1;
            }
            else {
                n += m_layer_nelz(iy) - 1;
            }
        }

        // allocate node-list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        // initialize counter: current index in the node-list "ret"
        size_t j = 0;

        // bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            // -- bottom node layer
            for (size_t iz = 1; iz < m_layer_nelz(iy); ++iz) {
                ret(j) = m_startNode(iy) + iz * (m_layer_nelx(iy) + 1);
                ++j;
            }
            // -- refinement layer
            if (m_refine(iy) == 2) {
                for (size_t iz = 0; iz < 2 * m_layer_nelz(iy); ++iz) {
                    ret(j) = m_startNode(iy) + iz * (m_layer_nelx(iy) + 1) + m_nnd(iy);
                    ++j;
                }
            }
        }

        // top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            // -- refinement layer
            if (m_refine(iy) == 2) {
                for (size_t iz = 0; iz < 2 * m_layer_nelz(iy); ++iz) {
                    ret(j) = m_startNode(iy) + iz * (m_layer_nelx(iy) + 1) + m_nnd(iy);
                    ++j;
                }
            }
            // -- top node layer
            for (size_t iz = 1; iz < m_layer_nelz(iy); ++iz) {
                ret(j) = m_startNode(iy + 1) + iz * (m_layer_nelx(iy) + 1);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesRightFace_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // number of boundary nodes
        // - initialize
        size_t n = 0;
        // - bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            if (m_refine(iy) == 2) {
                n += m_layer_nelz(iy) * 3 - 1;
            }
            else {
                n += m_layer_nelz(iy) - 1;
            }
        }
        // - top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            if (m_refine(iy) == 2) {
                n += m_layer_nelz(iy) * 3 - 1;
            }
            else {
                n += m_layer_nelz(iy) - 1;
            }
        }

        // allocate node-list
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({n});

        // initialize counter: current index in the node-list "ret"
        size_t j = 0;

        // bottom half: bottom node layer (+ middle node layer)
        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            // -- bottom node layer
            for (size_t iz = 1; iz < m_layer_nelz(iy); ++iz) {
                ret(j) = m_startNode(iy) + iz * (m_layer_nelx(iy) + 1) + m_layer_nelx(iy);
                ++j;
            }
            // -- refinement layer
            if (m_refine(iy) == 2) {
                for (size_t iz = 0; iz < 2 * m_layer_nelz(iy); ++iz) {
                    ret(j) = m_startNode(iy) + m_nnd(iy) + iz * (m_layer_nelx(iy) + 1) +
                             m_layer_nelx(iy);
                    ++j;
                }
            }
        }

        // top half: (middle node layer +) top node layer
        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            // -- refinement layer
            if (m_refine(iy) == 2) {
                for (size_t iz = 0; iz < 2 * m_layer_nelz(iy); ++iz) {
                    ret(j) = m_startNode(iy) + m_nnd(iy) + iz * (m_layer_nelx(iy) + 1) +
                             m_layer_nelx(iy);
                    ++j;
                }
            }
            // -- top node layer
            for (size_t iz = 1; iz < m_layer_nelz(iy); ++iz) {
                ret(j) = m_startNode(iy + 1) + iz * (m_layer_nelx(iy) + 1) + m_layer_nelx(iy);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomFace_impl() const
    {
        // allocate node list
        array_type::tensor<size_t, 1> ret =
            xt::empty<size_t>({(m_layer_nelx(0) - 1) * (m_layer_nelz(0) - 1)});

        // counter
        size_t j = 0;

        // fill node list
        for (size_t ix = 1; ix < m_layer_nelx(0); ++ix) {
            for (size_t iz = 1; iz < m_layer_nelz(0); ++iz) {
                ret(j) = m_startNode(0) + ix + iz * (m_layer_nelx(0) + 1);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopFace_impl() const
    {
        // number of element layers in y-direction
        size_t nely = static_cast<size_t>(m_nhy.size());

        // allocate node list
        array_type::tensor<size_t, 1> ret =
            xt::empty<size_t>({(m_layer_nelx(nely - 1) - 1) * (m_layer_nelz(nely - 1) - 1)});

        // counter
        size_t j = 0;

        // fill node list
        for (size_t ix = 1; ix < m_layer_nelx(nely - 1); ++ix) {
            for (size_t iz = 1; iz < m_layer_nelz(nely - 1); ++iz) {
                ret(j) = m_startNode(nely) + ix + iz * (m_layer_nelx(nely - 1) + 1);
                ++j;
            }
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontBottomEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(0) + 1});

        for (size_t ix = 0; ix < m_layer_nelx(0) + 1; ++ix) {
            ret(ix) = m_startNode(0) + ix;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontTopEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(nely - 1) + 1});

        for (size_t ix = 0; ix < m_layer_nelx(nely - 1) + 1; ++ix) {
            ret(ix) = m_startNode(nely) + ix;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontLeftEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({nely + 1});

        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            ret(iy) = m_startNode(iy);
        }

        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            ret(iy + 1) = m_startNode(iy + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontRightEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({nely + 1});

        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            ret(iy) = m_startNode(iy) + m_layer_nelx(iy);
        }

        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            ret(iy + 1) = m_startNode(iy + 1) + m_layer_nelx(iy);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackBottomEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(0) + 1});

        for (size_t ix = 0; ix < m_layer_nelx(0) + 1; ++ix) {
            ret(ix) = m_startNode(0) + ix + (m_layer_nelx(0) + 1) * (m_layer_nelz(0));
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackTopEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(nely - 1) + 1});

        for (size_t ix = 0; ix < m_layer_nelx(nely - 1) + 1; ++ix) {
            ret(ix) =
                m_startNode(nely) + ix + (m_layer_nelx(nely - 1) + 1) * (m_layer_nelz(nely - 1));
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackLeftEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({nely + 1});

        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            ret(iy) = m_startNode(iy) + (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy));
        }

        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            ret(iy + 1) = m_startNode(iy + 1) + (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy));
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackRightEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({nely + 1});

        for (size_t iy = 0; iy < (nely + 1) / 2; ++iy) {
            ret(iy) =
                m_startNode(iy) + m_layer_nelx(iy) + (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy));
        }

        for (size_t iy = (nely - 1) / 2; iy < nely; ++iy) {
            ret(iy + 1) = m_startNode(iy + 1) + m_layer_nelx(iy) +
                          (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy));
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomLeftEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelz(0) + 1});

        for (size_t iz = 0; iz < m_layer_nelz(0) + 1; ++iz) {
            ret(iz) = m_startNode(0) + iz * (m_layer_nelx(0) + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomRightEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelz(0) + 1});

        for (size_t iz = 0; iz < m_layer_nelz(0) + 1; ++iz) {
            ret(iz) = m_startNode(0) + m_layer_nelx(0) + iz * (m_layer_nelx(0) + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopLeftEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelz(nely - 1) + 1});

        for (size_t iz = 0; iz < m_layer_nelz(nely - 1) + 1; ++iz) {
            ret(iz) = m_startNode(nely) + iz * (m_layer_nelx(nely - 1) + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopRightEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelz(nely - 1) + 1});

        for (size_t iz = 0; iz < m_layer_nelz(nely - 1) + 1; ++iz) {
            ret(iz) =
                m_startNode(nely) + m_layer_nelx(nely - 1) + iz * (m_layer_nelx(nely - 1) + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontBottomOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(0) - 1});

        for (size_t ix = 1; ix < m_layer_nelx(0); ++ix) {
            ret(ix - 1) = m_startNode(0) + ix;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontTopOpenEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(nely - 1) - 1});

        for (size_t ix = 1; ix < m_layer_nelx(nely - 1); ++ix) {
            ret(ix - 1) = m_startNode(nely) + ix;
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontLeftOpenEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({nely - 1});

        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            ret(iy - 1) = m_startNode(iy);
        }

        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            ret(iy) = m_startNode(iy + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesFrontRightOpenEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({nely - 1});

        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            ret(iy - 1) = m_startNode(iy) + m_layer_nelx(iy);
        }

        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            ret(iy) = m_startNode(iy + 1) + m_layer_nelx(iy);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackBottomOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(0) - 1});

        for (size_t ix = 1; ix < m_layer_nelx(0); ++ix) {
            ret(ix - 1) = m_startNode(0) + ix + (m_layer_nelx(0) + 1) * (m_layer_nelz(0));
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackTopOpenEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelx(nely - 1) - 1});

        for (size_t ix = 1; ix < m_layer_nelx(nely - 1); ++ix) {
            ret(ix - 1) =
                m_startNode(nely) + ix + (m_layer_nelx(nely - 1) + 1) * (m_layer_nelz(nely - 1));
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackLeftOpenEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({nely - 1});

        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            ret(iy - 1) = m_startNode(iy) + (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy));
        }

        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            ret(iy) = m_startNode(iy + 1) + (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy));
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBackRightOpenEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({nely - 1});

        for (size_t iy = 1; iy < (nely + 1) / 2; ++iy) {
            ret(iy - 1) =
                m_startNode(iy) + m_layer_nelx(iy) + (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy));
        }

        for (size_t iy = (nely - 1) / 2; iy < nely - 1; ++iy) {
            ret(iy) = m_startNode(iy + 1) + m_layer_nelx(iy) +
                      (m_layer_nelx(iy) + 1) * (m_layer_nelz(iy));
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomLeftOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelz(0) - 1});

        for (size_t iz = 1; iz < m_layer_nelz(0); ++iz) {
            ret(iz - 1) = m_startNode(0) + iz * (m_layer_nelx(0) + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesBottomRightOpenEdge_impl() const
    {
        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelz(0) - 1});

        for (size_t iz = 1; iz < m_layer_nelz(0); ++iz) {
            ret(iz - 1) = m_startNode(0) + m_layer_nelx(0) + iz * (m_layer_nelx(0) + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopLeftOpenEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelz(nely - 1) - 1});

        for (size_t iz = 1; iz < m_layer_nelz(nely - 1); ++iz) {
            ret(iz - 1) = m_startNode(nely) + iz * (m_layer_nelx(nely - 1) + 1);
        }

        return ret;
    }

    array_type::tensor<size_t, 1> nodesTopRightOpenEdge_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());

        array_type::tensor<size_t, 1> ret = xt::empty<size_t>({m_layer_nelz(nely - 1) - 1});

        for (size_t iz = 1; iz < m_layer_nelz(nely - 1); ++iz) {
            ret(iz - 1) =
                m_startNode(nely) + m_layer_nelx(nely - 1) + iz * (m_layer_nelx(nely - 1) + 1);
        }

        return ret;
    }

    size_t nodesFrontBottomLeftCorner_impl() const
    {
        return m_startNode(0);
    }

    size_t nodesFrontBottomRightCorner_impl() const
    {
        return m_startNode(0) + m_layer_nelx(0);
    }

    size_t nodesFrontTopLeftCorner_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());
        return m_startNode(nely);
    }

    size_t nodesFrontTopRightCorner_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());
        return m_startNode(nely) + m_layer_nelx(nely - 1);
    }

    size_t nodesBackBottomLeftCorner_impl() const
    {
        return m_startNode(0) + (m_layer_nelx(0) + 1) * (m_layer_nelz(0));
    }

    size_t nodesBackBottomRightCorner_impl() const
    {
        return m_startNode(0) + m_layer_nelx(0) + (m_layer_nelx(0) + 1) * (m_layer_nelz(0));
    }

    size_t nodesBackTopLeftCorner_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());
        return m_startNode(nely) + (m_layer_nelx(nely - 1) + 1) * (m_layer_nelz(nely - 1));
    }

    size_t nodesBackTopRightCorner_impl() const
    {
        size_t nely = static_cast<size_t>(m_nhy.size());
        return m_startNode(nely) + m_layer_nelx(nely - 1) +
               (m_layer_nelx(nely - 1) + 1) * (m_layer_nelz(nely - 1));
    }

    double m_h; ///< See h()
    size_t m_nelx; ///< See nelx()
    size_t m_nely; ///< See nely()
    size_t m_nelz; ///< See nely()
    size_t m_nelem; ///< See nelem()
    size_t m_nnode; ///< See nnode()
    size_t m_nne; ///< See nne()
    size_t m_ndim; ///< See ndim()
    double m_Lx; ///< mesh size in "x"
    double m_Lz; ///< mesh size in "z"
    array_type::tensor<size_t, 1> m_layer_nelx; ///< number of elem in "x" per element layer in "y"
    array_type::tensor<size_t, 1> m_layer_nelz; ///< number of elem in "z" per element layer in "y"
    array_type::tensor<size_t, 1> m_nnd; ///< num nodes in the main node layer per node layer in "y"
    array_type::tensor<size_t, 1> m_nhx; ///< element size in x-direction per element layer in "y"
    array_type::tensor<size_t, 1> m_nhy; ///< element size in y-direction per element layer in "y"
    array_type::tensor<size_t, 1> m_nhz; ///< element size in z-direction per element layer in "y"
    array_type::tensor<int, 1>
        m_refine; ///< refine direction (-1:no refine, 0:"x", 2:"z") per element layer in "y"
    array_type::tensor<size_t, 1> m_startElem; ///< start element per element layer in "y"
    array_type::tensor<size_t, 1> m_startNode; ///< start node per node layer in "y"
};

} // namespace Hex8
} // namespace Mesh
} // namespace GooseFEM

#endif

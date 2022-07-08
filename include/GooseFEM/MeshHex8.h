/**
Generate simple meshes of 8-noded hexahedral elements in 3d (GooseFEM::Mesh::ElementType::Hex8).

\file MeshHex8.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESHHEX8_H
#define GOOSEFEM_MESHHEX8_H

#include "Mesh.h"
#include "config.h"

namespace GooseFEM {
namespace Mesh {

/**
Simple meshes of 8-noded hexahedral elements in 3d (ElementType::Hex8).
*/
namespace Hex8 {

/**
Regular mesh: equi-sized elements.
*/
class Regular : public RegularBase3d<Regular> {
public:
    /**
    Constructor.

    \param nelx Number of elements in horizontal (x) direction.
    \param nely Number of elements in vertical (y) direction.
    \param nelz Number of elements in vertical (z) direction.
    \param h Edge size (width == height == depth).
    */
    Regular(size_t nelx, size_t nely, size_t nelz, double h = 1.0);

private:
    friend class RegularBase<Regular>;
    friend class RegularBase3d<Regular>;

    size_t nelx_impl() const;
    size_t nely_impl() const;
    size_t nelz_impl() const;
    ElementType getElementType_impl() const;
    array_type::tensor<double, 2> coor_impl() const;
    array_type::tensor<size_t, 2> conn_impl() const;
    array_type::tensor<size_t, 1> nodesFront_impl() const;
    array_type::tensor<size_t, 1> nodesBack_impl() const;
    array_type::tensor<size_t, 1> nodesLeft_impl() const;
    array_type::tensor<size_t, 1> nodesRight_impl() const;
    array_type::tensor<size_t, 1> nodesBottom_impl() const;
    array_type::tensor<size_t, 1> nodesTop_impl() const;
    array_type::tensor<size_t, 1> nodesFrontFace_impl() const;
    array_type::tensor<size_t, 1> nodesBackFace_impl() const;
    array_type::tensor<size_t, 1> nodesLeftFace_impl() const;
    array_type::tensor<size_t, 1> nodesRightFace_impl() const;
    array_type::tensor<size_t, 1> nodesBottomFace_impl() const;
    array_type::tensor<size_t, 1> nodesTopFace_impl() const;
    array_type::tensor<size_t, 1> nodesFrontBottomEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontTopEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackBottomEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackTopEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontBottomOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontTopOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontRightOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackBottomOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackTopOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackRightOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomRightOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopRightOpenEdge_impl() const;
    size_t nodesFrontBottomLeftCorner_impl() const;
    size_t nodesFrontBottomRightCorner_impl() const;
    size_t nodesFrontTopLeftCorner_impl() const;
    size_t nodesFrontTopRightCorner_impl() const;
    size_t nodesBackBottomLeftCorner_impl() const;
    size_t nodesBackBottomRightCorner_impl() const;
    size_t nodesBackTopLeftCorner_impl() const;
    size_t nodesBackTopRightCorner_impl() const;

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
Mesh with fine middle layer, and coarser elements towards the top and bottom.
*/
class FineLayer : public RegularBase3d<FineLayer> {
public:
    /**
    Constructor.

    \param nelx Number of elements (along the middle layer) in horizontal (x) direction.
    \param nely Approximate equivalent number of elements in vertical (y) direction.
    \param nelz Number of elements (along the middle layer) in depth (z) direction.
    \param h Edge size (width == height == depth) of elements along the weak layer.

    \param nfine
        Extra number of fine layers around the middle layer.
        By default the element size is kept smaller than the distance to the middle layer.
    */
    FineLayer(size_t nelx, size_t nely, size_t nelz, double h = 1.0, size_t nfine = 1);

    /**
    Elements in the middle (fine) layer.

    \return List of element numbers.
    */
    array_type::tensor<size_t, 1> elementsMiddleLayer() const;

private:
    friend class RegularBase<FineLayer>;
    friend class RegularBase3d<FineLayer>;

    size_t nelx_impl() const;
    size_t nely_impl() const;
    size_t nelz_impl() const;
    ElementType getElementType_impl() const;
    array_type::tensor<double, 2> coor_impl() const;
    array_type::tensor<size_t, 2> conn_impl() const;
    array_type::tensor<size_t, 1> nodesFront_impl() const;
    array_type::tensor<size_t, 1> nodesBack_impl() const;
    array_type::tensor<size_t, 1> nodesLeft_impl() const;
    array_type::tensor<size_t, 1> nodesRight_impl() const;
    array_type::tensor<size_t, 1> nodesBottom_impl() const;
    array_type::tensor<size_t, 1> nodesTop_impl() const;
    array_type::tensor<size_t, 1> nodesFrontFace_impl() const;
    array_type::tensor<size_t, 1> nodesBackFace_impl() const;
    array_type::tensor<size_t, 1> nodesLeftFace_impl() const;
    array_type::tensor<size_t, 1> nodesRightFace_impl() const;
    array_type::tensor<size_t, 1> nodesBottomFace_impl() const;
    array_type::tensor<size_t, 1> nodesTopFace_impl() const;
    array_type::tensor<size_t, 1> nodesFrontBottomEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontTopEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackBottomEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackTopEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopLeftEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopRightEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontBottomOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontTopOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesFrontRightOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackBottomOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackTopOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBackRightOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesBottomRightOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopLeftOpenEdge_impl() const;
    array_type::tensor<size_t, 1> nodesTopRightOpenEdge_impl() const;
    size_t nodesFrontBottomLeftCorner_impl() const;
    size_t nodesFrontBottomRightCorner_impl() const;
    size_t nodesFrontTopLeftCorner_impl() const;
    size_t nodesFrontTopRightCorner_impl() const;
    size_t nodesBackBottomLeftCorner_impl() const;
    size_t nodesBackBottomRightCorner_impl() const;
    size_t nodesBackTopLeftCorner_impl() const;
    size_t nodesBackTopRightCorner_impl() const;

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

#include "MeshHex8.hpp"

#endif

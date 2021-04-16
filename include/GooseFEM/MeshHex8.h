/**
Generate simple meshes of 8-noded hexahedral elements in 3d (GooseFEM::Mesh::ElementType::Hex8).

\file MeshHex8.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESHHEX8_H
#define GOOSEFEM_MESHHEX8_H

#include "config.h"
#include "Mesh.h"

namespace GooseFEM {
namespace Mesh {
namespace Hex8 {

/**
Regular mesh: equi-sized elements.
*/
class Regular : public RegularBase3d {
public:

    /**
    Constructor.

    \param nelx Number of elements in horizontal (x) direction.
    \param nely Number of elements in vertical (y) direction.
    \param nelz Number of elements in vertical (z) direction.
    \param h Edge size (width == height == depth).
    */
    Regular(size_t nelx, size_t nely, size_t nelz, double h = 1.0);

    ElementType getElementType() const override;
    xt::xtensor<double, 2> coor() const override;
    xt::xtensor<size_t, 2> conn() const override;
    xt::xtensor<size_t, 1> nodesFront() const override;
    xt::xtensor<size_t, 1> nodesBack() const override;
    xt::xtensor<size_t, 1> nodesLeft() const override;
    xt::xtensor<size_t, 1> nodesRight() const override;
    xt::xtensor<size_t, 1> nodesBottom() const override;
    xt::xtensor<size_t, 1> nodesTop() const override;
    xt::xtensor<size_t, 1> nodesFrontFace() const override;
    xt::xtensor<size_t, 1> nodesBackFace() const override;
    xt::xtensor<size_t, 1> nodesLeftFace() const override;
    xt::xtensor<size_t, 1> nodesRightFace() const override;
    xt::xtensor<size_t, 1> nodesBottomFace() const override;
    xt::xtensor<size_t, 1> nodesTopFace() const override;
    xt::xtensor<size_t, 1> nodesFrontBottomEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontTopEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontLeftEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontRightEdge() const override;
    xt::xtensor<size_t, 1> nodesBackBottomEdge() const override;
    xt::xtensor<size_t, 1> nodesBackTopEdge() const override;
    xt::xtensor<size_t, 1> nodesBackLeftEdge() const override;
    xt::xtensor<size_t, 1> nodesBackRightEdge() const override;
    xt::xtensor<size_t, 1> nodesBottomLeftEdge() const override;
    xt::xtensor<size_t, 1> nodesBottomRightEdge() const override;
    xt::xtensor<size_t, 1> nodesTopLeftEdge() const override;
    xt::xtensor<size_t, 1> nodesTopRightEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontBottomOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontTopOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontLeftOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontRightOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBackBottomOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBackTopOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBackLeftOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBackRightOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBottomLeftOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBottomRightOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesTopLeftOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesTopRightOpenEdge() const override;
    size_t nodesFrontBottomLeftCorner() const override;
    size_t nodesFrontBottomRightCorner() const override;
    size_t nodesFrontTopLeftCorner() const override;
    size_t nodesFrontTopRightCorner() const override;
    size_t nodesBackBottomLeftCorner() const override;
    size_t nodesBackBottomRightCorner() const override;
    size_t nodesBackTopLeftCorner() const override;
    size_t nodesBackTopRightCorner() const override;
};

/**
Mesh with fine middle layer, and coarser elements towards the top and bottom.
*/
class FineLayer : public RegularBase3d {
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

    ElementType getElementType() const override;
    xt::xtensor<double, 2> coor() const override;
    xt::xtensor<size_t, 2> conn() const override;
    xt::xtensor<size_t, 1> nodesFront() const override;
    xt::xtensor<size_t, 1> nodesBack() const override;
    xt::xtensor<size_t, 1> nodesLeft() const override;
    xt::xtensor<size_t, 1> nodesRight() const override;
    xt::xtensor<size_t, 1> nodesBottom() const override;
    xt::xtensor<size_t, 1> nodesTop() const override;
    xt::xtensor<size_t, 1> nodesFrontFace() const override;
    xt::xtensor<size_t, 1> nodesBackFace() const override;
    xt::xtensor<size_t, 1> nodesLeftFace() const override;
    xt::xtensor<size_t, 1> nodesRightFace() const override;
    xt::xtensor<size_t, 1> nodesBottomFace() const override;
    xt::xtensor<size_t, 1> nodesTopFace() const override;
    xt::xtensor<size_t, 1> nodesFrontBottomEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontTopEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontLeftEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontRightEdge() const override;
    xt::xtensor<size_t, 1> nodesBackBottomEdge() const override;
    xt::xtensor<size_t, 1> nodesBackTopEdge() const override;
    xt::xtensor<size_t, 1> nodesBackLeftEdge() const override;
    xt::xtensor<size_t, 1> nodesBackRightEdge() const override;
    xt::xtensor<size_t, 1> nodesBottomLeftEdge() const override;
    xt::xtensor<size_t, 1> nodesBottomRightEdge() const override;
    xt::xtensor<size_t, 1> nodesTopLeftEdge() const override;
    xt::xtensor<size_t, 1> nodesTopRightEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontBottomOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontTopOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontLeftOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesFrontRightOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBackBottomOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBackTopOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBackLeftOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBackRightOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBottomLeftOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesBottomRightOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesTopLeftOpenEdge() const override;
    xt::xtensor<size_t, 1> nodesTopRightOpenEdge() const override;
    size_t nodesFrontBottomLeftCorner() const override;
    size_t nodesFrontBottomRightCorner() const override;
    size_t nodesFrontTopLeftCorner() const override;
    size_t nodesFrontTopRightCorner() const override;
    size_t nodesBackBottomLeftCorner() const override;
    size_t nodesBackBottomRightCorner() const override;
    size_t nodesBackTopLeftCorner() const override;
    size_t nodesBackTopRightCorner() const override;

    size_t nelx() const override;
    size_t nely() const override;
    size_t nelz() const override;

    /**
    Elements in the middle (fine) layer.

    \return List of element numbers.
    */
    xt::xtensor<size_t, 1> elementsMiddleLayer() const;

private:
    double m_Lx;                         ///< mesh size in "x"
    double m_Lz;                         ///< mesh size in "z"
    xt::xtensor<size_t, 1> m_layer_nelx; ///< number of elements in "x" per element layer in "y"
    xt::xtensor<size_t, 1> m_layer_nelz; ///< number of elements in "z" per element layer in "y"
    xt::xtensor<size_t, 1> m_nnd;        ///< number of nodes in the main node layer per node layer in "y"
    xt::xtensor<size_t, 1> m_nhx;        ///< element size in x-direction per element layer in "y"
    xt::xtensor<size_t, 1> m_nhy;        ///< element size in y-direction per element layer in "y"
    xt::xtensor<size_t, 1> m_nhz;        ///< element size in z-direction per element layer in "y"
    xt::xtensor<int, 1> m_refine;        ///< refine direction (-1:no refine, 0:"x", 2:"z") per element layer in "y"
    xt::xtensor<size_t, 1> m_startElem;  ///< start element per element layer in "y"
    xt::xtensor<size_t, 1> m_startNode;  ///< start node per node layer in "y"
};

} // namespace Hex8
} // namespace Mesh
} // namespace GooseFEM

#include "MeshHex8.hpp"

#endif

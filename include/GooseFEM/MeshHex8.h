/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_MESHHEX8_H
#define GOOSEFEM_MESHHEX8_H

#include "config.h"

namespace GooseFEM {
namespace Mesh {
namespace Hex8 {

class Regular {
public:
    Regular(size_t nelx, size_t nely, size_t nelz, double h = 1.);

    // size
    size_t nelem() const; // number of elements
    size_t nnode() const; // number of nodes
    size_t nne() const;   // number of nodes-per-element
    size_t ndim() const;  // number of dimensions

    // type
    ElementType getElementType() const;

    // mesh
    xt::xtensor<double, 2> coor() const; // nodal positions [nnode, ndim]
    xt::xtensor<size_t, 2> conn() const; // connectivity [nelem, nne]

    // boundary nodes: planes
    xt::xtensor<size_t, 1> nodesFront() const;
    xt::xtensor<size_t, 1> nodesBack() const;
    xt::xtensor<size_t, 1> nodesLeft() const;
    xt::xtensor<size_t, 1> nodesRight() const;
    xt::xtensor<size_t, 1> nodesBottom() const;
    xt::xtensor<size_t, 1> nodesTop() const;

    // boundary nodes: faces
    xt::xtensor<size_t, 1> nodesFrontFace() const;
    xt::xtensor<size_t, 1> nodesBackFace() const;
    xt::xtensor<size_t, 1> nodesLeftFace() const;
    xt::xtensor<size_t, 1> nodesRightFace() const;
    xt::xtensor<size_t, 1> nodesBottomFace() const;
    xt::xtensor<size_t, 1> nodesTopFace() const;

    // boundary nodes: edges
    xt::xtensor<size_t, 1> nodesFrontBottomEdge() const;
    xt::xtensor<size_t, 1> nodesFrontTopEdge() const;
    xt::xtensor<size_t, 1> nodesFrontLeftEdge() const;
    xt::xtensor<size_t, 1> nodesFrontRightEdge() const;
    xt::xtensor<size_t, 1> nodesBackBottomEdge() const;
    xt::xtensor<size_t, 1> nodesBackTopEdge() const;
    xt::xtensor<size_t, 1> nodesBackLeftEdge() const;
    xt::xtensor<size_t, 1> nodesBackRightEdge() const;
    xt::xtensor<size_t, 1> nodesBottomLeftEdge() const;
    xt::xtensor<size_t, 1> nodesBottomRightEdge() const;
    xt::xtensor<size_t, 1> nodesTopLeftEdge() const;
    xt::xtensor<size_t, 1> nodesTopRightEdge() const;

    // boundary nodes: edges (aliases)
    xt::xtensor<size_t, 1> nodesBottomFrontEdge() const;
    xt::xtensor<size_t, 1> nodesBottomBackEdge() const;
    xt::xtensor<size_t, 1> nodesTopFrontEdge() const;
    xt::xtensor<size_t, 1> nodesTopBackEdge() const;
    xt::xtensor<size_t, 1> nodesLeftBottomEdge() const;
    xt::xtensor<size_t, 1> nodesLeftFrontEdge() const;
    xt::xtensor<size_t, 1> nodesLeftBackEdge() const;
    xt::xtensor<size_t, 1> nodesLeftTopEdge() const;
    xt::xtensor<size_t, 1> nodesRightBottomEdge() const;
    xt::xtensor<size_t, 1> nodesRightTopEdge() const;
    xt::xtensor<size_t, 1> nodesRightFrontEdge() const;
    xt::xtensor<size_t, 1> nodesRightBackEdge() const;

    // boundary nodes: edges, without corners
    xt::xtensor<size_t, 1> nodesFrontBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesFrontTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesFrontLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesFrontRightOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBackBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBackTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBackLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBackRightOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBottomLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBottomRightOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopRightOpenEdge() const;

    // boundary nodes: edges, without corners (aliases)
    xt::xtensor<size_t, 1> nodesBottomFrontOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBottomBackOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopFrontOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopBackOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftFrontOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftBackOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightFrontOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightBackOpenEdge() const;

    // boundary nodes: corners
    size_t nodesFrontBottomLeftCorner() const;
    size_t nodesFrontBottomRightCorner() const;
    size_t nodesFrontTopLeftCorner() const;
    size_t nodesFrontTopRightCorner() const;
    size_t nodesBackBottomLeftCorner() const;
    size_t nodesBackBottomRightCorner() const;
    size_t nodesBackTopLeftCorner() const;
    size_t nodesBackTopRightCorner() const;

    // boundary nodes: corners (aliases)
    size_t nodesFrontLeftBottomCorner() const;
    size_t nodesBottomFrontLeftCorner() const;
    size_t nodesBottomLeftFrontCorner() const;
    size_t nodesLeftFrontBottomCorner() const;
    size_t nodesLeftBottomFrontCorner() const;
    size_t nodesFrontRightBottomCorner() const;
    size_t nodesBottomFrontRightCorner() const;
    size_t nodesBottomRightFrontCorner() const;
    size_t nodesRightFrontBottomCorner() const;
    size_t nodesRightBottomFrontCorner() const;
    size_t nodesFrontLeftTopCorner() const;
    size_t nodesTopFrontLeftCorner() const;
    size_t nodesTopLeftFrontCorner() const;
    size_t nodesLeftFrontTopCorner() const;
    size_t nodesLeftTopFrontCorner() const;
    size_t nodesFrontRightTopCorner() const;
    size_t nodesTopFrontRightCorner() const;
    size_t nodesTopRightFrontCorner() const;
    size_t nodesRightFrontTopCorner() const;
    size_t nodesRightTopFrontCorner() const;
    size_t nodesBackLeftBottomCorner() const;
    size_t nodesBottomBackLeftCorner() const;
    size_t nodesBottomLeftBackCorner() const;
    size_t nodesLeftBackBottomCorner() const;
    size_t nodesLeftBottomBackCorner() const;
    size_t nodesBackRightBottomCorner() const;
    size_t nodesBottomBackRightCorner() const;
    size_t nodesBottomRightBackCorner() const;
    size_t nodesRightBackBottomCorner() const;
    size_t nodesRightBottomBackCorner() const;
    size_t nodesBackLeftTopCorner() const;
    size_t nodesTopBackLeftCorner() const;
    size_t nodesTopLeftBackCorner() const;
    size_t nodesLeftBackTopCorner() const;
    size_t nodesLeftTopBackCorner() const;
    size_t nodesBackRightTopCorner() const;
    size_t nodesTopBackRightCorner() const;
    size_t nodesTopRightBackCorner() const;
    size_t nodesRightBackTopCorner() const;
    size_t nodesRightTopBackCorner() const;

    // DOF-numbers for each component of each node (sequential)
    xt::xtensor<size_t, 2> dofs() const;

    // DOF-numbers for the case that the periodicity if fully eliminated
    xt::xtensor<size_t, 2> dofsPeriodic() const;

    // periodic node pairs [:,2]: (independent, dependent)
    xt::xtensor<size_t, 2> nodesPeriodic() const;

    // front-bottom-left node, used as reference for periodicity
    size_t nodesOrigin() const;

private:
    double m_h;                     // elementary element edge-size (in all directions)
    size_t m_nelx;                  // number of elements in x-direction (length == "m_nelx * m_h")
    size_t m_nely;                  // number of elements in y-direction (length == "m_nely * m_h")
    size_t m_nelz;                  // number of elements in z-direction (length == "m_nely * m_h")
    size_t m_nelem;                 // number of elements
    size_t m_nnode;                 // number of nodes
    static const size_t m_nne = 8;  // number of nodes-per-element
    static const size_t m_ndim = 3; // number of dimensions
};

class FineLayer {
public:
    FineLayer(size_t nelx, size_t nely, size_t nelz, double h = 1.0, size_t nfine = 1);

    // size
    size_t nelem() const; // number of elements
    size_t nnode() const; // number of nodes
    size_t nne() const;   // number of nodes-per-element
    size_t ndim() const;  // number of dimensions
    size_t nelx() const;  // number of elements in x-direction
    size_t nely() const;  // number of elements in y-direction
    size_t nelz() const;  // number of elements in y-direction

    // type
    ElementType getElementType() const;

    // mesh
    xt::xtensor<double, 2> coor() const; // nodal positions [nnode, ndim]
    xt::xtensor<size_t, 2> conn() const; // connectivity [nelem, nne]

    // element sets
    xt::xtensor<size_t, 1> elementsMiddleLayer() const; // elements in the middle (fine) layer

    // boundary nodes: planes
    xt::xtensor<size_t, 1> nodesFront() const;
    xt::xtensor<size_t, 1> nodesBack() const;
    xt::xtensor<size_t, 1> nodesLeft() const;
    xt::xtensor<size_t, 1> nodesRight() const;
    xt::xtensor<size_t, 1> nodesBottom() const;
    xt::xtensor<size_t, 1> nodesTop() const;

    // boundary nodes: faces
    xt::xtensor<size_t, 1> nodesFrontFace() const;
    xt::xtensor<size_t, 1> nodesBackFace() const;
    xt::xtensor<size_t, 1> nodesLeftFace() const;
    xt::xtensor<size_t, 1> nodesRightFace() const;
    xt::xtensor<size_t, 1> nodesBottomFace() const;
    xt::xtensor<size_t, 1> nodesTopFace() const;

    // boundary nodes: edges
    xt::xtensor<size_t, 1> nodesFrontBottomEdge() const;
    xt::xtensor<size_t, 1> nodesFrontTopEdge() const;
    xt::xtensor<size_t, 1> nodesFrontLeftEdge() const;
    xt::xtensor<size_t, 1> nodesFrontRightEdge() const;
    xt::xtensor<size_t, 1> nodesBackBottomEdge() const;
    xt::xtensor<size_t, 1> nodesBackTopEdge() const;
    xt::xtensor<size_t, 1> nodesBackLeftEdge() const;
    xt::xtensor<size_t, 1> nodesBackRightEdge() const;
    xt::xtensor<size_t, 1> nodesBottomLeftEdge() const;
    xt::xtensor<size_t, 1> nodesBottomRightEdge() const;
    xt::xtensor<size_t, 1> nodesTopLeftEdge() const;
    xt::xtensor<size_t, 1> nodesTopRightEdge() const;

    // boundary nodes: edges (aliases)
    xt::xtensor<size_t, 1> nodesBottomFrontEdge() const;
    xt::xtensor<size_t, 1> nodesBottomBackEdge() const;
    xt::xtensor<size_t, 1> nodesTopFrontEdge() const;
    xt::xtensor<size_t, 1> nodesTopBackEdge() const;
    xt::xtensor<size_t, 1> nodesLeftBottomEdge() const;
    xt::xtensor<size_t, 1> nodesLeftFrontEdge() const;
    xt::xtensor<size_t, 1> nodesLeftBackEdge() const;
    xt::xtensor<size_t, 1> nodesLeftTopEdge() const;
    xt::xtensor<size_t, 1> nodesRightBottomEdge() const;
    xt::xtensor<size_t, 1> nodesRightTopEdge() const;
    xt::xtensor<size_t, 1> nodesRightFrontEdge() const;
    xt::xtensor<size_t, 1> nodesRightBackEdge() const;

    // boundary nodes: edges, without corners
    xt::xtensor<size_t, 1> nodesFrontBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesFrontTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesFrontLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesFrontRightOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBackBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBackTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBackLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBackRightOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBottomLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBottomRightOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopRightOpenEdge() const;

    // boundary nodes: edges, without corners (aliases)
    xt::xtensor<size_t, 1> nodesBottomFrontOpenEdge() const;
    xt::xtensor<size_t, 1> nodesBottomBackOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopFrontOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopBackOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftFrontOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftBackOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightFrontOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightBackOpenEdge() const;

    // boundary nodes: corners
    size_t nodesFrontBottomLeftCorner() const;
    size_t nodesFrontBottomRightCorner() const;
    size_t nodesFrontTopLeftCorner() const;
    size_t nodesFrontTopRightCorner() const;
    size_t nodesBackBottomLeftCorner() const;
    size_t nodesBackBottomRightCorner() const;
    size_t nodesBackTopLeftCorner() const;
    size_t nodesBackTopRightCorner() const;

    // boundary nodes: corners (aliases)
    size_t nodesFrontLeftBottomCorner() const;
    size_t nodesBottomFrontLeftCorner() const;
    size_t nodesBottomLeftFrontCorner() const;
    size_t nodesLeftFrontBottomCorner() const;
    size_t nodesLeftBottomFrontCorner() const;
    size_t nodesFrontRightBottomCorner() const;
    size_t nodesBottomFrontRightCorner() const;
    size_t nodesBottomRightFrontCorner() const;
    size_t nodesRightFrontBottomCorner() const;
    size_t nodesRightBottomFrontCorner() const;
    size_t nodesFrontLeftTopCorner() const;
    size_t nodesTopFrontLeftCorner() const;
    size_t nodesTopLeftFrontCorner() const;
    size_t nodesLeftFrontTopCorner() const;
    size_t nodesLeftTopFrontCorner() const;
    size_t nodesFrontRightTopCorner() const;
    size_t nodesTopFrontRightCorner() const;
    size_t nodesTopRightFrontCorner() const;
    size_t nodesRightFrontTopCorner() const;
    size_t nodesRightTopFrontCorner() const;
    size_t nodesBackLeftBottomCorner() const;
    size_t nodesBottomBackLeftCorner() const;
    size_t nodesBottomLeftBackCorner() const;
    size_t nodesLeftBackBottomCorner() const;
    size_t nodesLeftBottomBackCorner() const;
    size_t nodesBackRightBottomCorner() const;
    size_t nodesBottomBackRightCorner() const;
    size_t nodesBottomRightBackCorner() const;
    size_t nodesRightBackBottomCorner() const;
    size_t nodesRightBottomBackCorner() const;
    size_t nodesBackLeftTopCorner() const;
    size_t nodesTopBackLeftCorner() const;
    size_t nodesTopLeftBackCorner() const;
    size_t nodesLeftBackTopCorner() const;
    size_t nodesLeftTopBackCorner() const;
    size_t nodesBackRightTopCorner() const;
    size_t nodesTopBackRightCorner() const;
    size_t nodesTopRightBackCorner() const;
    size_t nodesRightBackTopCorner() const;
    size_t nodesRightTopBackCorner() const;

    // DOF-numbers for each component of each node (sequential)
    xt::xtensor<size_t, 2> dofs() const;

    // DOF-numbers for the case that the periodicity if fully eliminated
    xt::xtensor<size_t, 2> dofsPeriodic() const;

    // periodic node pairs [:,2]: (independent, dependent)
    xt::xtensor<size_t, 2> nodesPeriodic() const;

    // front-bottom-left node, used as reference for periodicity
    size_t nodesOrigin() const;

private:
    double m_h;                         // elementary element edge-size (in all directions)
    double m_Lx;                        // mesh size in "x"
    double m_Lz;                        // mesh size in "z"
    size_t m_nelem;                     // number of elements
    size_t m_nnode;                     // number of nodes
    static const size_t m_nne = 8;      // number of nodes-per-element
    static const size_t m_ndim = 3;     // number of dimensions
    xt::xtensor<size_t, 1> m_nelx;      // number of elements in "x" (*)
    xt::xtensor<size_t, 1> m_nelz;      // number of elements in "z" (*)
    xt::xtensor<size_t, 1> m_nnd;       // number of nodes in the main node layer (**)
    xt::xtensor<size_t, 1> m_nhx;       // element size in x-direction (*)
    xt::xtensor<size_t, 1> m_nhy;       // element size in y-direction (*)
    xt::xtensor<size_t, 1> m_nhz;       // element size in z-direction (*)
    xt::xtensor<int, 1> m_refine;       // refine direction (-1:no refine, 0:"x", 2:"z") (*)
    xt::xtensor<size_t, 1> m_startElem; // start element (*)
    xt::xtensor<size_t, 1> m_startNode; // start node (**)
    // (*) per element layer in "y"
    // (**) per node layer in "y"
};

} // namespace Hex8
} // namespace Mesh
} // namespace GooseFEM

#include "MeshHex8.hpp"

#endif

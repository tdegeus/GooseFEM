/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHHEX8_H
#define GOOSEFEM_MESHHEX8_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ===================================== GooseFEM::Mesh::Hex8 ======================================

namespace GooseFEM {
namespace Mesh {
namespace Hex8 {

// ----------------------------------------- regular mesh ------------------------------------------

class Regular
{
private:
  double m_h;      // elementary element edge-size (in all directions)
  size_t m_nelx;   // number of elements in x-direction (length == "m_nelx * m_h")
  size_t m_nely;   // number of elements in y-direction (length == "m_nely * m_h")
  size_t m_nelz;   // number of elements in z-direction (length == "m_nely * m_h")
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=8;  // number of nodes-per-element
  size_t m_ndim=3; // number of dimensions

public:
  // mesh with "nelx*nely*nelz" 'elements' of edge size "h"
  Regular(size_t nelx, size_t nely, size_t nelz, double h=1.);
  // sizes
  size_t nelem() const;                       // number of elements
  size_t nnode() const;                       // number of nodes
  size_t nne() const;                         // number of nodes-per-element
  size_t ndim() const;                        // number of dimensions
  // mesh
  MatD   coor() const;                        // nodal positions [nnode ,ndim]
  MatS   conn() const;                        // connectivity    [nelem ,nne ]
  // boundary nodes: planes
  ColS   nodesFront() const;                  // node-numbers along the front  plane
  ColS   nodesBack() const;                   // node-numbers along the back   plane
  ColS   nodesLeft() const;                   // node-numbers along the left   plane
  ColS   nodesRight() const;                  // node-numbers along the right  plane
  ColS   nodesBottom() const;                 // node-numbers along the bottom plane
  ColS   nodesTop() const;                    // node-numbers along the top    plane
  // boundary nodes: faces
  ColS   nodesFrontFace() const;              // node-numbers along the front  face
  ColS   nodesBackFace() const;               // node-numbers along the back   face
  ColS   nodesLeftFace() const;               // node-numbers along the left   face
  ColS   nodesRightFace() const;              // node-numbers along the right  face
  ColS   nodesBottomFace() const;             // node-numbers along the bottom face
  ColS   nodesTopFace() const;                // node-numbers along the top    face
  // boundary nodes: edges
  ColS   nodesFrontBottomEdge() const;        // node-numbers along the front  - bottom edge
  ColS   nodesFrontTopEdge() const;           // node-numbers along the front  - top    edge
  ColS   nodesFrontLeftEdge() const;          // node-numbers along the front  - left   edge
  ColS   nodesFrontRightEdge() const;         // node-numbers along the front  - right  edge
  ColS   nodesBackBottomEdge() const;         // node-numbers along the back   - bottom edge
  ColS   nodesBackTopEdge() const;            // node-numbers along the back   - top    edge
  ColS   nodesBackLeftEdge() const;           // node-numbers along the back   - left   edge
  ColS   nodesBackRightEdge() const;          // node-numbers along the back   - right  edge
  ColS   nodesBottomLeftEdge() const;         // node-numbers along the bottom - left   edge
  ColS   nodesBottomRightEdge() const;        // node-numbers along the bottom - right  edge
  ColS   nodesTopLeftEdge() const;            // node-numbers along the top    - left   edge
  ColS   nodesTopRightEdge() const;           // node-numbers along the top    - right  edge
  // boundary nodes: edges (aliases)
  ColS   nodesBottomFrontEdge() const;        // alias, see above: nodesFrontBottomEdge
  ColS   nodesBottomBackEdge() const;         // alias, see above: nodesBackBottomEdge
  ColS   nodesTopFrontEdge() const;           // alias, see above: nodesFrontTopEdge
  ColS   nodesTopBackEdge() const;            // alias, see above: nodesBackTopEdge
  ColS   nodesLeftBottomEdge() const;         // alias, see above: nodesBottomLeftEdge
  ColS   nodesLeftFrontEdge() const;          // alias, see above: nodesFrontLeftEdge
  ColS   nodesLeftBackEdge() const;           // alias, see above: nodesBackLeftEdge
  ColS   nodesLeftTopEdge() const;            // alias, see above: nodesTopLeftEdge
  ColS   nodesRightBottomEdge() const;        // alias, see above: nodesBottomRightEdge
  ColS   nodesRightTopEdge() const;           // alias, see above: nodesTopRightEdge
  ColS   nodesRightFrontEdge() const;         // alias, see above: nodesFrontRightEdge
  ColS   nodesRightBackEdge() const;          // alias, see above: nodesBackRightEdge
  // boundary nodes: edges, without corners
  ColS   nodesFrontBottomOpenEdge() const;    // node-numbers along the front  - bottom edge
  ColS   nodesFrontTopOpenEdge() const;       // node-numbers along the front  - top    edge
  ColS   nodesFrontLeftOpenEdge() const;      // node-numbers along the front  - left   edge
  ColS   nodesFrontRightOpenEdge() const;     // node-numbers along the front  - right  edge
  ColS   nodesBackBottomOpenEdge() const;     // node-numbers along the back   - bottom edge
  ColS   nodesBackTopOpenEdge() const;        // node-numbers along the back   - top    edge
  ColS   nodesBackLeftOpenEdge() const;       // node-numbers along the back   - left   edge
  ColS   nodesBackRightOpenEdge() const;      // node-numbers along the back   - right  edge
  ColS   nodesBottomLeftOpenEdge() const;     // node-numbers along the bottom - left   edge
  ColS   nodesBottomRightOpenEdge() const;    // node-numbers along the bottom - right  edge
  ColS   nodesTopLeftOpenEdge() const;        // node-numbers along the top    - left   edge
  ColS   nodesTopRightOpenEdge() const;       // node-numbers along the top    - right  edge
  // boundary nodes: edges, without corners (aliases)
  ColS   nodesBottomFrontOpenEdge() const;    // alias, see above: nodesFrontBottomOpenEdge
  ColS   nodesBottomBackOpenEdge() const;     // alias, see above: nodesBackBottomOpenEdge
  ColS   nodesTopFrontOpenEdge() const;       // alias, see above: nodesFrontTopOpenEdge
  ColS   nodesTopBackOpenEdge() const;        // alias, see above: nodesBackTopOpenEdge
  ColS   nodesLeftBottomOpenEdge() const;     // alias, see above: nodesBottomLeftOpenEdge
  ColS   nodesLeftFrontOpenEdge() const;      // alias, see above: nodesFrontLeftOpenEdge
  ColS   nodesLeftBackOpenEdge() const;       // alias, see above: nodesBackLeftOpenEdge
  ColS   nodesLeftTopOpenEdge() const;        // alias, see above: nodesTopLeftOpenEdge
  ColS   nodesRightBottomOpenEdge() const;    // alias, see above: nodesBottomRightOpenEdge
  ColS   nodesRightTopOpenEdge() const;       // alias, see above: nodesTopRightOpenEdge
  ColS   nodesRightFrontOpenEdge() const;     // alias, see above: nodesFrontRightOpenEdge
  ColS   nodesRightBackOpenEdge() const;      // alias, see above: nodesBackRightOpenEdge
  // boundary nodes: corners
  size_t nodesFrontBottomLeftCorner() const;  // node-number of the front - bottom - left  corner
  size_t nodesFrontBottomRightCorner() const; // node-number of the front - bottom - right corner
  size_t nodesFrontTopLeftCorner() const;     // node-number of the front - top    - left  corner
  size_t nodesFrontTopRightCorner() const;    // node-number of the front - top    - right corner
  size_t nodesBackBottomLeftCorner() const;   // node-number of the back  - bottom - left  corner
  size_t nodesBackBottomRightCorner() const;  // node-number of the back  - bottom - right corner
  size_t nodesBackTopLeftCorner() const;      // node-number of the back  - top    - left  corner
  size_t nodesBackTopRightCorner() const;     // node-number of the back  - top    - right corner
  // boundary nodes: corners (aliases)
  size_t nodesFrontLeftBottomCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesBottomFrontLeftCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesBottomLeftFrontCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesLeftFrontBottomCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesLeftBottomFrontCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesFrontRightBottomCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesBottomFrontRightCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesBottomRightFrontCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesRightFrontBottomCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesRightBottomFrontCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesFrontLeftTopCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesTopFrontLeftCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesTopLeftFrontCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesLeftFrontTopCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesLeftTopFrontCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesFrontRightTopCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesTopFrontRightCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesTopRightFrontCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesRightFrontTopCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesRightTopFrontCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesBackLeftBottomCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBottomBackLeftCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBottomLeftBackCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesLeftBackBottomCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesLeftBottomBackCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBackRightBottomCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBottomBackRightCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBottomRightBackCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesRightBackBottomCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesRightBottomBackCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBackLeftTopCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesTopBackLeftCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesTopLeftBackCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesLeftBackTopCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesLeftTopBackCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesBackRightTopCorner() const;     // alias, see above: nodesBackTopRightCorner
  size_t nodesTopBackRightCorner() const;     // alias, see above: nodesBackTopRightCorner
  size_t nodesTopRightBackCorner() const;     // alias, see above: nodesBackTopRightCorner
  size_t nodesRightBackTopCorner() const;     // alias, see above: nodesBackTopRightCorner
  size_t nodesRightTopBackCorner() const;     // alias, see above: nodesBackTopRightCorner
  // periodicity
  MatS   nodesPeriodic() const;               // periodic node pairs [:,2]: (independent, dependent)
  size_t nodesOrigin() const;                 // front-bottom-left node, used as reference for periodicity
  MatS   dofs() const;                        // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic() const;                // ,, for the case that the periodicity if fully eliminated
};

// ---------------------- mesh with a fine layer that exponentially coarsens -----------------------

class FineLayer
{
private:
  double m_h;                 // elementary element edge-size (in all directions)
  double m_Lx, m_Lz;          // mesh size in "x" and "z"
  ColS   m_nelx, m_nelz;      // number of elements in "x" and "z"             (per el.layer in "y")
  ColS   m_nnd;               // total number of nodes in the main node layer  (per nd.layer in "y")
  ColS   m_nhx, m_nhy, m_nhz; // element size in each direction                (per el.layer in "y")
  ColI   m_refine;            // refine direction (-1:no refine, 0:"x", 2:"z") (per el.layer in "y")
  ColS   m_startElem;         // start element                                 (per el.layer in "y")
  ColS   m_startNode;         // start node                                    (per nd.layer in "y")
  size_t m_nelem;             // number of elements
  size_t m_nnode;             // number of nodes
  size_t m_nne=8;             // number of nodes-per-element
  size_t m_ndim=3;            // number of dimensions

public:
  // mesh with "nelx*nely*nelz" elements of edge size "h"; elements are coarsened in "y"-direction
  FineLayer(size_t nelx, size_t nely, size_t nelz, double h=1., size_t nfine=1);
  // sizes
  size_t nelem() const;                       // number of elements
  size_t nnode() const;                       // number of nodes
  size_t nne() const;                         // number of nodes-per-element
  size_t ndim() const;                        // number of dimensions
  size_t shape(size_t i) const;               // actual shape in a certain direction
  // mesh
  MatD   coor() const;                        // nodal positions [nnode ,ndim]
  MatS   conn() const;                        // connectivity    [nelem ,nne ]
  // element sets
  ColS   elementsMiddleLayer() const;         // elements in the middle, fine, layer
  // boundary nodes: planes
  ColS   nodesFront() const;                  // node-numbers along the front  plane
  ColS   nodesBack() const;                   // node-numbers along the back   plane
  ColS   nodesLeft() const;                   // node-numbers along the left   plane
  ColS   nodesRight() const;                  // node-numbers along the right  plane
  ColS   nodesBottom() const;                 // node-numbers along the bottom plane
  ColS   nodesTop() const;                    // node-numbers along the top    plane
  // boundary nodes: faces
  ColS   nodesFrontFace() const;              // node-numbers along the front  face
  ColS   nodesBackFace() const;               // node-numbers along the back   face
  ColS   nodesLeftFace() const;               // node-numbers along the left   face
  ColS   nodesRightFace() const;              // node-numbers along the right  face
  ColS   nodesBottomFace() const;             // node-numbers along the bottom face
  ColS   nodesTopFace() const;                // node-numbers along the top    face
  // boundary nodes: edges
  ColS   nodesFrontBottomEdge() const;        // node-numbers along the front  - bottom edge
  ColS   nodesFrontTopEdge() const;           // node-numbers along the front  - top    edge
  ColS   nodesFrontLeftEdge() const;          // node-numbers along the front  - left   edge
  ColS   nodesFrontRightEdge() const;         // node-numbers along the front  - right  edge
  ColS   nodesBackBottomEdge() const;         // node-numbers along the back   - bottom edge
  ColS   nodesBackTopEdge() const;            // node-numbers along the back   - top    edge
  ColS   nodesBackLeftEdge() const;           // node-numbers along the back   - left   edge
  ColS   nodesBackRightEdge() const;          // node-numbers along the back   - right  edge
  ColS   nodesBottomLeftEdge() const;         // node-numbers along the bottom - left   edge
  ColS   nodesBottomRightEdge() const;        // node-numbers along the bottom - right  edge
  ColS   nodesTopLeftEdge() const;            // node-numbers along the top    - left   edge
  ColS   nodesTopRightEdge() const;           // node-numbers along the top    - right  edge
  // boundary nodes: edges (aliases)
  ColS   nodesBottomFrontEdge() const;        // alias, see above: nodesFrontBottomEdge
  ColS   nodesBottomBackEdge() const;         // alias, see above: nodesBackBottomEdge
  ColS   nodesTopFrontEdge() const;           // alias, see above: nodesFrontTopEdge
  ColS   nodesTopBackEdge() const;            // alias, see above: nodesBackTopEdge
  ColS   nodesLeftBottomEdge() const;         // alias, see above: nodesBottomLeftEdge
  ColS   nodesLeftFrontEdge() const;          // alias, see above: nodesFrontLeftEdge
  ColS   nodesLeftBackEdge() const;           // alias, see above: nodesBackLeftEdge
  ColS   nodesLeftTopEdge() const;            // alias, see above: nodesTopLeftEdge
  ColS   nodesRightBottomEdge() const;        // alias, see above: nodesBottomRightEdge
  ColS   nodesRightTopEdge() const;           // alias, see above: nodesTopRightEdge
  ColS   nodesRightFrontEdge() const;         // alias, see above: nodesFrontRightEdge
  ColS   nodesRightBackEdge() const;          // alias, see above: nodesBackRightEdge
  // boundary nodes: edges, without corners
  ColS   nodesFrontBottomOpenEdge() const;    // node-numbers along the front  - bottom edge
  ColS   nodesFrontTopOpenEdge() const;       // node-numbers along the front  - top    edge
  ColS   nodesFrontLeftOpenEdge() const;      // node-numbers along the front  - left   edge
  ColS   nodesFrontRightOpenEdge() const;     // node-numbers along the front  - right  edge
  ColS   nodesBackBottomOpenEdge() const;     // node-numbers along the back   - bottom edge
  ColS   nodesBackTopOpenEdge() const;        // node-numbers along the back   - top    edge
  ColS   nodesBackLeftOpenEdge() const;       // node-numbers along the back   - left   edge
  ColS   nodesBackRightOpenEdge() const;      // node-numbers along the back   - right  edge
  ColS   nodesBottomLeftOpenEdge() const;     // node-numbers along the bottom - left   edge
  ColS   nodesBottomRightOpenEdge() const;    // node-numbers along the bottom - right  edge
  ColS   nodesTopLeftOpenEdge() const;        // node-numbers along the top    - left   edge
  ColS   nodesTopRightOpenEdge() const;       // node-numbers along the top    - right  edge
  // boundary nodes: edges, without corners (aliases)
  ColS   nodesBottomFrontOpenEdge() const;    // alias, see above: nodesFrontBottomOpenEdge
  ColS   nodesBottomBackOpenEdge() const;     // alias, see above: nodesBackBottomOpenEdge
  ColS   nodesTopFrontOpenEdge() const;       // alias, see above: nodesFrontTopOpenEdge
  ColS   nodesTopBackOpenEdge() const;        // alias, see above: nodesBackTopOpenEdge
  ColS   nodesLeftBottomOpenEdge() const;     // alias, see above: nodesBottomLeftOpenEdge
  ColS   nodesLeftFrontOpenEdge() const;      // alias, see above: nodesFrontLeftOpenEdge
  ColS   nodesLeftBackOpenEdge() const;       // alias, see above: nodesBackLeftOpenEdge
  ColS   nodesLeftTopOpenEdge() const;        // alias, see above: nodesTopLeftOpenEdge
  ColS   nodesRightBottomOpenEdge() const;    // alias, see above: nodesBottomRightOpenEdge
  ColS   nodesRightTopOpenEdge() const;       // alias, see above: nodesTopRightOpenEdge
  ColS   nodesRightFrontOpenEdge() const;     // alias, see above: nodesFrontRightOpenEdge
  ColS   nodesRightBackOpenEdge() const;      // alias, see above: nodesBackRightOpenEdge
  // boundary nodes: corners
  size_t nodesFrontBottomLeftCorner() const;  // node-number of the front - bottom - left  corner
  size_t nodesFrontBottomRightCorner() const; // node-number of the front - bottom - right corner
  size_t nodesFrontTopLeftCorner() const;     // node-number of the front - top    - left  corner
  size_t nodesFrontTopRightCorner() const;    // node-number of the front - top    - right corner
  size_t nodesBackBottomLeftCorner() const;   // node-number of the back  - bottom - left  corner
  size_t nodesBackBottomRightCorner() const;  // node-number of the back  - bottom - right corner
  size_t nodesBackTopLeftCorner() const;      // node-number of the back  - top    - left  corner
  size_t nodesBackTopRightCorner() const;     // node-number of the back  - top    - right corner
  // boundary nodes: corners (aliases)
  size_t nodesFrontLeftBottomCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesBottomFrontLeftCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesBottomLeftFrontCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesLeftFrontBottomCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesLeftBottomFrontCorner() const;  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesFrontRightBottomCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesBottomFrontRightCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesBottomRightFrontCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesRightFrontBottomCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesRightBottomFrontCorner() const; // alias, see above: nodesFrontBottomRightCorner
  size_t nodesFrontLeftTopCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesTopFrontLeftCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesTopLeftFrontCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesLeftFrontTopCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesLeftTopFrontCorner() const;     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesFrontRightTopCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesTopFrontRightCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesTopRightFrontCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesRightFrontTopCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesRightTopFrontCorner() const;    // alias, see above: nodesFrontTopRightCorner
  size_t nodesBackLeftBottomCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBottomBackLeftCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBottomLeftBackCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesLeftBackBottomCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesLeftBottomBackCorner() const;   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBackRightBottomCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBottomBackRightCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBottomRightBackCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesRightBackBottomCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesRightBottomBackCorner() const;  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBackLeftTopCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesTopBackLeftCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesTopLeftBackCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesLeftBackTopCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesLeftTopBackCorner() const;      // alias, see above: nodesBackTopLeftCorner
  size_t nodesBackRightTopCorner() const;     // alias, see above: nodesBackTopRightCorner
  size_t nodesTopBackRightCorner() const;     // alias, see above: nodesBackTopRightCorner
  size_t nodesTopRightBackCorner() const;     // alias, see above: nodesBackTopRightCorner
  size_t nodesRightBackTopCorner() const;     // alias, see above: nodesBackTopRightCorner
  size_t nodesRightTopBackCorner() const;     // alias, see above: nodesBackTopRightCorner
  // periodicity
  MatS   nodesPeriodic() const; // periodic node pairs [:,2]: (independent, dependent)
  size_t nodesOrigin() const;   // front-bottom-left node, used as reference for periodicity
  MatS   dofs() const;          // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic() const;  // ,, for the case that the periodicity if fully eliminated
};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

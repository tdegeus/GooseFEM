/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_MESHHEX8_H
#define XGOOSEFEM_MESHHEX8_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ===================================== xGooseFEM::Mesh::Hex8 ======================================

namespace xGooseFEM {
namespace Mesh {
namespace Hex8 {

// ----------------------------------------- regular mesh ------------------------------------------

class Regular
{
private:
  double m_h;                   // elementary element edge-size (in all directions)
  size_t m_nelx;                // number of elements in x-direction (length == "m_nelx * m_h")
  size_t m_nely;                // number of elements in y-direction (length == "m_nely * m_h")
  size_t m_nelz;                // number of elements in z-direction (length == "m_nely * m_h")
  size_t m_nelem;               // number of elements
  size_t m_nnode;               // number of nodes
  static const size_t m_nne=8;  // number of nodes-per-element
  static const size_t m_ndim=3; // number of dimensions

public:
  // mesh with "nelx*nely*nelz" 'elements' of edge size "h"
  Regular(size_t nelx, size_t nely, size_t nelz, double h=1.);
  // sizes
  size_t nelem() const;                       // number of elements
  size_t nnode() const;                       // number of nodes
  size_t nne() const;                         // number of nodes-per-element
  size_t ndim() const;                        // number of dimensions
  // mesh
  xt::xtensor<double,2> coor() const;                        // nodal positions [nnode ,ndim]
  xt::xtensor<size_t,2> conn() const;                        // connectivity    [nelem ,nne ]
  // boundary nodes: planes
  xt::xtensor<size_t,1> nodesFront() const;                  // node-numbers along the front  plane
  xt::xtensor<size_t,1> nodesBack() const;                   // node-numbers along the back   plane
  xt::xtensor<size_t,1> nodesLeft() const;                   // node-numbers along the left   plane
  xt::xtensor<size_t,1> nodesRight() const;                  // node-numbers along the right  plane
  xt::xtensor<size_t,1> nodesBottom() const;                 // node-numbers along the bottom plane
  xt::xtensor<size_t,1> nodesTop() const;                    // node-numbers along the top    plane
  // boundary nodes: faces
  xt::xtensor<size_t,1> nodesFrontFace() const;              // node-numbers along the front  face
  xt::xtensor<size_t,1> nodesBackFace() const;               // node-numbers along the back   face
  xt::xtensor<size_t,1> nodesLeftFace() const;               // node-numbers along the left   face
  xt::xtensor<size_t,1> nodesRightFace() const;              // node-numbers along the right  face
  xt::xtensor<size_t,1> nodesBottomFace() const;             // node-numbers along the bottom face
  xt::xtensor<size_t,1> nodesTopFace() const;                // node-numbers along the top    face
  // boundary nodes: edges
  xt::xtensor<size_t,1> nodesFrontBottomEdge() const;        // node-numbers along the front  - bottom edge
  xt::xtensor<size_t,1> nodesFrontTopEdge() const;           // node-numbers along the front  - top    edge
  xt::xtensor<size_t,1> nodesFrontLeftEdge() const;          // node-numbers along the front  - left   edge
  xt::xtensor<size_t,1> nodesFrontRightEdge() const;         // node-numbers along the front  - right  edge
  xt::xtensor<size_t,1> nodesBackBottomEdge() const;         // node-numbers along the back   - bottom edge
  xt::xtensor<size_t,1> nodesBackTopEdge() const;            // node-numbers along the back   - top    edge
  xt::xtensor<size_t,1> nodesBackLeftEdge() const;           // node-numbers along the back   - left   edge
  xt::xtensor<size_t,1> nodesBackRightEdge() const;          // node-numbers along the back   - right  edge
  xt::xtensor<size_t,1> nodesBottomLeftEdge() const;         // node-numbers along the bottom - left   edge
  xt::xtensor<size_t,1> nodesBottomRightEdge() const;        // node-numbers along the bottom - right  edge
  xt::xtensor<size_t,1> nodesTopLeftEdge() const;            // node-numbers along the top    - left   edge
  xt::xtensor<size_t,1> nodesTopRightEdge() const;           // node-numbers along the top    - right  edge
  // boundary nodes: edges (aliases)
  xt::xtensor<size_t,1> nodesBottomFrontEdge() const;        // alias, see above: nodesFrontBottomEdge
  xt::xtensor<size_t,1> nodesBottomBackEdge() const;         // alias, see above: nodesBackBottomEdge
  xt::xtensor<size_t,1> nodesTopFrontEdge() const;           // alias, see above: nodesFrontTopEdge
  xt::xtensor<size_t,1> nodesTopBackEdge() const;            // alias, see above: nodesBackTopEdge
  xt::xtensor<size_t,1> nodesLeftBottomEdge() const;         // alias, see above: nodesBottomLeftEdge
  xt::xtensor<size_t,1> nodesLeftFrontEdge() const;          // alias, see above: nodesFrontLeftEdge
  xt::xtensor<size_t,1> nodesLeftBackEdge() const;           // alias, see above: nodesBackLeftEdge
  xt::xtensor<size_t,1> nodesLeftTopEdge() const;            // alias, see above: nodesTopLeftEdge
  xt::xtensor<size_t,1> nodesRightBottomEdge() const;        // alias, see above: nodesBottomRightEdge
  xt::xtensor<size_t,1> nodesRightTopEdge() const;           // alias, see above: nodesTopRightEdge
  xt::xtensor<size_t,1> nodesRightFrontEdge() const;         // alias, see above: nodesFrontRightEdge
  xt::xtensor<size_t,1> nodesRightBackEdge() const;          // alias, see above: nodesBackRightEdge
  // boundary nodes: edges, without corners
  xt::xtensor<size_t,1> nodesFrontBottomOpenEdge() const;    // node-numbers along the front  - bottom edge
  xt::xtensor<size_t,1> nodesFrontTopOpenEdge() const;       // node-numbers along the front  - top    edge
  xt::xtensor<size_t,1> nodesFrontLeftOpenEdge() const;      // node-numbers along the front  - left   edge
  xt::xtensor<size_t,1> nodesFrontRightOpenEdge() const;     // node-numbers along the front  - right  edge
  xt::xtensor<size_t,1> nodesBackBottomOpenEdge() const;     // node-numbers along the back   - bottom edge
  xt::xtensor<size_t,1> nodesBackTopOpenEdge() const;        // node-numbers along the back   - top    edge
  xt::xtensor<size_t,1> nodesBackLeftOpenEdge() const;       // node-numbers along the back   - left   edge
  xt::xtensor<size_t,1> nodesBackRightOpenEdge() const;      // node-numbers along the back   - right  edge
  xt::xtensor<size_t,1> nodesBottomLeftOpenEdge() const;     // node-numbers along the bottom - left   edge
  xt::xtensor<size_t,1> nodesBottomRightOpenEdge() const;    // node-numbers along the bottom - right  edge
  xt::xtensor<size_t,1> nodesTopLeftOpenEdge() const;        // node-numbers along the top    - left   edge
  xt::xtensor<size_t,1> nodesTopRightOpenEdge() const;       // node-numbers along the top    - right  edge
  // boundary nodes: edges, without corners (aliases)
  xt::xtensor<size_t,1> nodesBottomFrontOpenEdge() const;    // alias, see above: nodesFrontBottomOpenEdge
  xt::xtensor<size_t,1> nodesBottomBackOpenEdge() const;     // alias, see above: nodesBackBottomOpenEdge
  xt::xtensor<size_t,1> nodesTopFrontOpenEdge() const;       // alias, see above: nodesFrontTopOpenEdge
  xt::xtensor<size_t,1> nodesTopBackOpenEdge() const;        // alias, see above: nodesBackTopOpenEdge
  xt::xtensor<size_t,1> nodesLeftBottomOpenEdge() const;     // alias, see above: nodesBottomLeftOpenEdge
  xt::xtensor<size_t,1> nodesLeftFrontOpenEdge() const;      // alias, see above: nodesFrontLeftOpenEdge
  xt::xtensor<size_t,1> nodesLeftBackOpenEdge() const;       // alias, see above: nodesBackLeftOpenEdge
  xt::xtensor<size_t,1> nodesLeftTopOpenEdge() const;        // alias, see above: nodesTopLeftOpenEdge
  xt::xtensor<size_t,1> nodesRightBottomOpenEdge() const;    // alias, see above: nodesBottomRightOpenEdge
  xt::xtensor<size_t,1> nodesRightTopOpenEdge() const;       // alias, see above: nodesTopRightOpenEdge
  xt::xtensor<size_t,1> nodesRightFrontOpenEdge() const;     // alias, see above: nodesFrontRightOpenEdge
  xt::xtensor<size_t,1> nodesRightBackOpenEdge() const;      // alias, see above: nodesBackRightOpenEdge
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
  xt::xtensor<size_t,2> nodesPeriodic() const; // periodic node pairs [:,2]: (independent, dependent)
  size_t                nodesOrigin() const;   // front-bottom-left node, used as reference for periodicity
  xt::xtensor<size_t,2> dofs() const;          // DOF-numbers for each component of each node (sequential)
  xt::xtensor<size_t,2> dofsPeriodic() const;  // ,, for the case that the periodicity if fully eliminated
};

// ---------------------- mesh with a fine layer that exponentially coarsens -----------------------

class FineLayer
{
private:
  double m_h;                                // elementary element edge-size (in all directions)
  double m_Lx, m_Lz;                         // mesh size in "x" and "z"
  size_t m_nelem;                            // number of elements
  size_t m_nnode;                            // number of nodes
  static const size_t m_nne=8;               // number of nodes-per-element
  static const size_t m_ndim=3;              // number of dimensions
  xt::xtensor<size_t,1> m_nelx, m_nelz;      // number of elements in "x" and "z"             (per el.layer in "y")
  xt::xtensor<size_t,1> m_nnd;               // total number of nodes in the main node layer  (per nd.layer in "y")
  xt::xtensor<size_t,1> m_nhx, m_nhy, m_nhz; // element size in each direction                (per el.layer in "y")
  xt::xtensor<int   ,1> m_refine;            // refine direction (-1:no refine, 0:"x", 2:"z") (per el.layer in "y")
  xt::xtensor<size_t,1> m_startElem;         // start element                                 (per el.layer in "y")
  xt::xtensor<size_t,1> m_startNode;         // start node                                    (per nd.layer in "y")

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
  xt::xtensor<double,2> coor() const;                        // nodal positions [nnode ,ndim]
  xt::xtensor<size_t,2> conn() const;                        // connectivity    [nelem ,nne ]
  // element sets
  xt::xtensor<size_t,1> elementsMiddleLayer() const;         // elements in the middle, fine, layer
  // boundary nodes: planes
  xt::xtensor<size_t,1> nodesFront() const;                  // node-numbers along the front  plane
  xt::xtensor<size_t,1> nodesBack() const;                   // node-numbers along the back   plane
  xt::xtensor<size_t,1> nodesLeft() const;                   // node-numbers along the left   plane
  xt::xtensor<size_t,1> nodesRight() const;                  // node-numbers along the right  plane
  xt::xtensor<size_t,1> nodesBottom() const;                 // node-numbers along the bottom plane
  xt::xtensor<size_t,1> nodesTop() const;                    // node-numbers along the top    plane
  // boundary nodes: faces
  xt::xtensor<size_t,1> nodesFrontFace() const;              // node-numbers along the front  face
  xt::xtensor<size_t,1> nodesBackFace() const;               // node-numbers along the back   face
  xt::xtensor<size_t,1> nodesLeftFace() const;               // node-numbers along the left   face
  xt::xtensor<size_t,1> nodesRightFace() const;              // node-numbers along the right  face
  xt::xtensor<size_t,1> nodesBottomFace() const;             // node-numbers along the bottom face
  xt::xtensor<size_t,1> nodesTopFace() const;                // node-numbers along the top    face
  // boundary nodes: edges
  xt::xtensor<size_t,1> nodesFrontBottomEdge() const;        // node-numbers along the front  - bottom edge
  xt::xtensor<size_t,1> nodesFrontTopEdge() const;           // node-numbers along the front  - top    edge
  xt::xtensor<size_t,1> nodesFrontLeftEdge() const;          // node-numbers along the front  - left   edge
  xt::xtensor<size_t,1> nodesFrontRightEdge() const;         // node-numbers along the front  - right  edge
  xt::xtensor<size_t,1> nodesBackBottomEdge() const;         // node-numbers along the back   - bottom edge
  xt::xtensor<size_t,1> nodesBackTopEdge() const;            // node-numbers along the back   - top    edge
  xt::xtensor<size_t,1> nodesBackLeftEdge() const;           // node-numbers along the back   - left   edge
  xt::xtensor<size_t,1> nodesBackRightEdge() const;          // node-numbers along the back   - right  edge
  xt::xtensor<size_t,1> nodesBottomLeftEdge() const;         // node-numbers along the bottom - left   edge
  xt::xtensor<size_t,1> nodesBottomRightEdge() const;        // node-numbers along the bottom - right  edge
  xt::xtensor<size_t,1> nodesTopLeftEdge() const;            // node-numbers along the top    - left   edge
  xt::xtensor<size_t,1> nodesTopRightEdge() const;           // node-numbers along the top    - right  edge
  // boundary nodes: edges (aliases)
  xt::xtensor<size_t,1> nodesBottomFrontEdge() const;        // alias, see above: nodesFrontBottomEdge
  xt::xtensor<size_t,1> nodesBottomBackEdge() const;         // alias, see above: nodesBackBottomEdge
  xt::xtensor<size_t,1> nodesTopFrontEdge() const;           // alias, see above: nodesFrontTopEdge
  xt::xtensor<size_t,1> nodesTopBackEdge() const;            // alias, see above: nodesBackTopEdge
  xt::xtensor<size_t,1> nodesLeftBottomEdge() const;         // alias, see above: nodesBottomLeftEdge
  xt::xtensor<size_t,1> nodesLeftFrontEdge() const;          // alias, see above: nodesFrontLeftEdge
  xt::xtensor<size_t,1> nodesLeftBackEdge() const;           // alias, see above: nodesBackLeftEdge
  xt::xtensor<size_t,1> nodesLeftTopEdge() const;            // alias, see above: nodesTopLeftEdge
  xt::xtensor<size_t,1> nodesRightBottomEdge() const;        // alias, see above: nodesBottomRightEdge
  xt::xtensor<size_t,1> nodesRightTopEdge() const;           // alias, see above: nodesTopRightEdge
  xt::xtensor<size_t,1> nodesRightFrontEdge() const;         // alias, see above: nodesFrontRightEdge
  xt::xtensor<size_t,1> nodesRightBackEdge() const;          // alias, see above: nodesBackRightEdge
  // boundary nodes: edges, without corners
  xt::xtensor<size_t,1> nodesFrontBottomOpenEdge() const;    // node-numbers along the front  - bottom edge
  xt::xtensor<size_t,1> nodesFrontTopOpenEdge() const;       // node-numbers along the front  - top    edge
  xt::xtensor<size_t,1> nodesFrontLeftOpenEdge() const;      // node-numbers along the front  - left   edge
  xt::xtensor<size_t,1> nodesFrontRightOpenEdge() const;     // node-numbers along the front  - right  edge
  xt::xtensor<size_t,1> nodesBackBottomOpenEdge() const;     // node-numbers along the back   - bottom edge
  xt::xtensor<size_t,1> nodesBackTopOpenEdge() const;        // node-numbers along the back   - top    edge
  xt::xtensor<size_t,1> nodesBackLeftOpenEdge() const;       // node-numbers along the back   - left   edge
  xt::xtensor<size_t,1> nodesBackRightOpenEdge() const;      // node-numbers along the back   - right  edge
  xt::xtensor<size_t,1> nodesBottomLeftOpenEdge() const;     // node-numbers along the bottom - left   edge
  xt::xtensor<size_t,1> nodesBottomRightOpenEdge() const;    // node-numbers along the bottom - right  edge
  xt::xtensor<size_t,1> nodesTopLeftOpenEdge() const;        // node-numbers along the top    - left   edge
  xt::xtensor<size_t,1> nodesTopRightOpenEdge() const;       // node-numbers along the top    - right  edge
  // boundary nodes: edges, without corners (aliases)
  xt::xtensor<size_t,1> nodesBottomFrontOpenEdge() const;    // alias, see above: nodesFrontBottomOpenEdge
  xt::xtensor<size_t,1> nodesBottomBackOpenEdge() const;     // alias, see above: nodesBackBottomOpenEdge
  xt::xtensor<size_t,1> nodesTopFrontOpenEdge() const;       // alias, see above: nodesFrontTopOpenEdge
  xt::xtensor<size_t,1> nodesTopBackOpenEdge() const;        // alias, see above: nodesBackTopOpenEdge
  xt::xtensor<size_t,1> nodesLeftBottomOpenEdge() const;     // alias, see above: nodesBottomLeftOpenEdge
  xt::xtensor<size_t,1> nodesLeftFrontOpenEdge() const;      // alias, see above: nodesFrontLeftOpenEdge
  xt::xtensor<size_t,1> nodesLeftBackOpenEdge() const;       // alias, see above: nodesBackLeftOpenEdge
  xt::xtensor<size_t,1> nodesLeftTopOpenEdge() const;        // alias, see above: nodesTopLeftOpenEdge
  xt::xtensor<size_t,1> nodesRightBottomOpenEdge() const;    // alias, see above: nodesBottomRightOpenEdge
  xt::xtensor<size_t,1> nodesRightTopOpenEdge() const;       // alias, see above: nodesTopRightOpenEdge
  xt::xtensor<size_t,1> nodesRightFrontOpenEdge() const;     // alias, see above: nodesFrontRightOpenEdge
  xt::xtensor<size_t,1> nodesRightBackOpenEdge() const;      // alias, see above: nodesBackRightOpenEdge
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
  xt::xtensor<size_t,2> nodesPeriodic() const; // periodic node pairs [:,2]: (independent, dependent)
  size_t                nodesOrigin() const;   // front-bottom-left node, used as reference for periodicity
  xt::xtensor<size_t,2> dofs() const;          // DOF-numbers for each component of each node (sequential)
  xt::xtensor<size_t,2> dofsPeriodic() const;  // ,, for the case that the periodicity if fully eliminated
};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

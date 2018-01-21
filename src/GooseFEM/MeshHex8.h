/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHHEX8_H
#define GOOSEFEM_MESHHEX8_H

// -------------------------------------------------------------------------------------------------

#include "Mesh.h"

// ===================================== GooseFEM::Mesh::Hex8 ======================================

namespace GooseFEM {
namespace Mesh {
namespace Hex8 {

// ========================================== REGULAR MESH =========================================

class Regular
{
private:
  size_t m_nelx;     // number of 'pixels' x-direction (length == "m_nelx * m_h")
  size_t m_nely;     // number of 'pixels' y-direction (length == "m_nely * m_h")
  size_t m_nelz;     // number of 'pixels' z-direction (length == "m_nely * m_h")
  double m_h;      // size of the element edge (equal in both directions)
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=8;  // number of nodes-per-element
  size_t m_ndim=3; // number of dimensions

public:
  // mesh with "nelx*nely*nelz" 'pixels' and edge size "h"
  Regular(size_t nelx, size_t nely, size_t nelz, double h=1.);
  // sizes
  size_t nelem();                       // number of elements
  size_t nnode();                       // number of nodes
  size_t nne();                         // number of nodes-per-element
  size_t ndim();                        // number of dimensions
  // mesh
  MatD   coor();                        // nodal positions [ nnode , ndim ]
  MatS   conn();                        // connectivity    [ nelem , nne  ]
  // boundary nodes: planes
  ColS   nodesFront();                  // node-numbers along the front  plane
  ColS   nodesBack();                   // node-numbers along the back   plane
  ColS   nodesLeft();                   // node-numbers along the left   plane
  ColS   nodesRight();                  // node-numbers along the right  plane
  ColS   nodesBottom();                 // node-numbers along the bottom plane
  ColS   nodesTop();                    // node-numbers along the top    plane
  // boundary nodes: faces
  ColS   nodesFrontFace();              // node-numbers along the front  face
  ColS   nodesBackFace();               // node-numbers along the back   face
  ColS   nodesLeftFace();               // node-numbers along the left   face
  ColS   nodesRightFace();              // node-numbers along the right  face
  ColS   nodesBottomFace();             // node-numbers along the bottom face
  ColS   nodesTopFace();                // node-numbers along the top    face
  // boundary nodes: edges
  ColS   nodesFrontBottomEdge();        // node-numbers along the front  - bottom edge
  ColS   nodesFrontTopEdge();           // node-numbers along the front  - top    edge
  ColS   nodesFrontLeftEdge();          // node-numbers along the front  - left   edge
  ColS   nodesFrontRightEdge();         // node-numbers along the front  - right  edge
  ColS   nodesBackBottomEdge();         // node-numbers along the back   - bottom edge
  ColS   nodesBackTopEdge();            // node-numbers along the back   - top    edge
  ColS   nodesBackLeftEdge();           // node-numbers along the back   - left   edge
  ColS   nodesBackRightEdge();          // node-numbers along the back   - right  edge
  ColS   nodesBottomLeftEdge();         // node-numbers along the bottom - left   edge
  ColS   nodesBottomRightEdge();        // node-numbers along the bottom - right  edge
  ColS   nodesTopLeftEdge();            // node-numbers along the top    - left   edge
  ColS   nodesTopRightEdge();           // node-numbers along the top    - right  edge
  // boundary nodes: edges (aliases)
  ColS   nodesBottomFrontEdge();        // alias, see above: nodesFrontBottomEdge
  ColS   nodesBottomBackEdge();         // alias, see above: nodesBackBottomEdge
  ColS   nodesTopFrontEdge();           // alias, see above: nodesFrontTopEdge
  ColS   nodesTopBackEdge();            // alias, see above: nodesBackTopEdge
  ColS   nodesLeftBottomEdge();         // alias, see above: nodesBottomLeftEdge
  ColS   nodesLeftFrontEdge();          // alias, see above: nodesFrontLeftEdge
  ColS   nodesLeftBackEdge();           // alias, see above: nodesBackLeftEdge
  ColS   nodesLeftTopEdge();            // alias, see above: nodesTopLeftEdge
  ColS   nodesRightBottomEdge();        // alias, see above: nodesBottomRightEdge
  ColS   nodesRightTopEdge();           // alias, see above: nodesTopRightEdge
  ColS   nodesRightFrontEdge();         // alias, see above: nodesFrontRightEdge
  ColS   nodesRightBackEdge();          // alias, see above: nodesBackRightEdge
  // boundary nodes: corners
  size_t nodesFrontBottomLeftCorner();  // node-number of the front - bottom - left  corner
  size_t nodesFrontBottomRightCorner(); // node-number of the front - bottom - right corner
  size_t nodesFrontTopLeftCorner();     // node-number of the front - top    - left  corner
  size_t nodesFrontTopRightCorner();    // node-number of the front - top    - right corner
  size_t nodesBackBottomLeftCorner();   // node-number of the back  - bottom - left  corner
  size_t nodesBackBottomRightCorner();  // node-number of the back  - bottom - right corner
  size_t nodesBackTopLeftCorner();      // node-number of the back  - top    - left  corner
  size_t nodesBackTopRightCorner();     // node-number of the back  - top    - right corner
  // boundary nodes: corners (aliases)
  size_t nodesFrontLeftBottomCorner();  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesBottomFrontLeftCorner();  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesBottomLeftFrontCorner();  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesLeftFrontBottomCorner();  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesLeftBottomFrontCorner();  // alias, see above: nodesFrontBottomLeftCorner
  size_t nodesFrontRightBottomCorner(); // alias, see above: nodesFrontBottomRightCorner
  size_t nodesBottomFrontRightCorner(); // alias, see above: nodesFrontBottomRightCorner
  size_t nodesBottomRightFrontCorner(); // alias, see above: nodesFrontBottomRightCorner
  size_t nodesRightFrontBottomCorner(); // alias, see above: nodesFrontBottomRightCorner
  size_t nodesRightBottomFrontCorner(); // alias, see above: nodesFrontBottomRightCorner
  size_t nodesFrontLeftTopCorner();     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesTopFrontLeftCorner();     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesTopLeftFrontCorner();     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesLeftFrontTopCorner();     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesLeftTopFrontCorner();     // alias, see above: nodesFrontTopLeftCorner
  size_t nodesFrontRightTopCorner();    // alias, see above: nodesFrontTopRightCorner
  size_t nodesTopFrontRightCorner();    // alias, see above: nodesFrontTopRightCorner
  size_t nodesTopRightFrontCorner();    // alias, see above: nodesFrontTopRightCorner
  size_t nodesRightFrontTopCorner();    // alias, see above: nodesFrontTopRightCorner
  size_t nodesRightTopFrontCorner();    // alias, see above: nodesFrontTopRightCorner
  size_t nodesBackLeftBottomCorner();   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBottomBackLeftCorner();   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBottomLeftBackCorner();   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesLeftBackBottomCorner();   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesLeftBottomBackCorner();   // alias, see above: nodesBackBottomLeftCorner
  size_t nodesBackRightBottomCorner();  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBottomBackRightCorner();  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBottomRightBackCorner();  // alias, see above: nodesBackBottomRightCorner
  size_t nodesRightBackBottomCorner();  // alias, see above: nodesBackBottomRightCorner
  size_t nodesRightBottomBackCorner();  // alias, see above: nodesBackBottomRightCorner
  size_t nodesBackLeftTopCorner();      // alias, see above: nodesBackTopLeftCorner
  size_t nodesTopBackLeftCorner();      // alias, see above: nodesBackTopLeftCorner
  size_t nodesTopLeftBackCorner();      // alias, see above: nodesBackTopLeftCorner
  size_t nodesLeftBackTopCorner();      // alias, see above: nodesBackTopLeftCorner
  size_t nodesLeftTopBackCorner();      // alias, see above: nodesBackTopLeftCorner
  size_t nodesBackRightTopCorner();     // alias, see above: nodesBackTopRightCorner
  size_t nodesTopBackRightCorner();     // alias, see above: nodesBackTopRightCorner
  size_t nodesTopRightBackCorner();     // alias, see above: nodesBackTopRightCorner
  size_t nodesRightBackTopCorner();     // alias, see above: nodesBackTopRightCorner
  size_t nodesRightTopBackCorner();     // alias, see above: nodesBackTopRightCorner
  // periodicity
  MatS   nodesPeriodic();               // periodic node pairs [:,2]: (independent, dependent)
  size_t nodesOrigin();                 // front-bottom-left node, used as reference for periodicity
  MatS   dofs();                        // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic();                // ,, for the case that the periodicity if fully eliminated
};

// ====================== MESH WITH A FINE LAYER THAT EXPONENTIALLY COARSENS =======================

class FineLayer
{
private:
  double m_h;                 // elementary element edge-size (in all directions)
  double m_Lx, m_Lz;          // mesh size in "x" and "z"
  ColS   m_nelx, m_nelz;      // number of elements in "x" and "z"             (per el.layer in "y")
  ColS   m_nnd;               // total number of nodes in the main node layer  (per nd.layer in "y")
  ColS   m_nhx, m_nhy, m_nhz; // element size in each direction                (per el.layer in "y")
  ColI   m_refine;            // refine direction (-1: no refine, 0: "x", ...) (per el.layer in "y")
  ColS   m_startElem;         // start element                                 (per el.layer in "y")
  ColS   m_startNode;         // start node                                    (per nd.layer in "y")
  size_t m_nelem;             // number of elements
  size_t m_nnode;             // number of nodes
  size_t m_nne=8;             // number of nodes-per-element
  size_t m_ndim=3;            // number of dimensions

public:
  // mesh with "nelx*nely*nelz" 'pixels' of edge size "h"; elements are coarsened in "y"-direction
  FineLayer(size_t nelx, size_t nely, size_t nelz, double h=1., size_t nfine=1);
  // sizes
  size_t nelem();                   // number of elements
  size_t nnode();                   // number of nodes
  size_t nne();                     // number of nodes-per-element
  size_t ndim();                    // number of dimensions
  size_t shape(size_t i);           // actual shape in a certain direction
  // mesh
  MatD   coor();                    // nodal positions [nnode ,ndim]
  MatS   conn();                    // connectivity    [nelem ,nne ]
  // boundary nodes: planes
  ColS   nodesFront();                  // node-numbers along the front  plane
  ColS   nodesBack();                   // node-numbers along the back   plane
  ColS   nodesLeft();                   // node-numbers along the left   plane
  ColS   nodesRight();                  // node-numbers along the right  plane
  ColS   nodesBottom();                 // node-numbers along the bottom plane
  ColS   nodesTop();                    // node-numbers along the top    plane
  // boundary nodes: faces
  ColS   nodesFrontFace();              // node-numbers along the front  face
  ColS   nodesBackFace();               // node-numbers along the back   face
  ColS   nodesLeftFace();               // node-numbers along the left   face
  ColS   nodesRightFace();              // node-numbers along the right  face
  ColS   nodesBottomFace();             // node-numbers along the bottom face
  ColS   nodesTopFace();                // node-numbers along the top    face
  // boundary nodes: edges
  ColS   nodesFrontBottomEdge();        // node-numbers along the front  - bottom edge
  ColS   nodesFrontTopEdge();           // node-numbers along the front  - top    edge
  ColS   nodesFrontLeftEdge();          // node-numbers along the front  - left   edge
  ColS   nodesFrontRightEdge();         // node-numbers along the front  - right  edge
  ColS   nodesBackBottomEdge();         // node-numbers along the back   - bottom edge
  ColS   nodesBackTopEdge();            // node-numbers along the back   - top    edge
  ColS   nodesBackLeftEdge();           // node-numbers along the back   - left   edge
  ColS   nodesBackRightEdge();          // node-numbers along the back   - right  edge
  ColS   nodesBottomLeftEdge();         // node-numbers along the bottom - left   edge
  ColS   nodesBottomRightEdge();        // node-numbers along the bottom - right  edge
  ColS   nodesTopLeftEdge();            // node-numbers along the top    - left   edge
  ColS   nodesTopRightEdge();           // node-numbers along the top    - right  edge

  // boundary nodes: edges
  ColS   nodesFrontBottomOpenEdge();        // node-numbers along the front  - bottom edge
  ColS   nodesFrontTopOpenEdge();           // node-numbers along the front  - top    edge
  ColS   nodesFrontLeftOpenEdge();          // node-numbers along the front  - left   edge
  ColS   nodesFrontRightOpenEdge();         // node-numbers along the front  - right  edge
  ColS   nodesBackBottomOpenEdge();         // node-numbers along the back   - bottom edge
  ColS   nodesBackTopOpenEdge();            // node-numbers along the back   - top    edge
  ColS   nodesBackLeftOpenEdge();           // node-numbers along the back   - left   edge
  ColS   nodesBackRightOpenEdge();          // node-numbers along the back   - right  edge
  ColS   nodesBottomLeftOpenEdge();         // node-numbers along the bottom - left   edge
  ColS   nodesBottomRightOpenEdge();        // node-numbers along the bottom - right  edge
  ColS   nodesTopLeftOpenEdge();            // node-numbers along the top    - left   edge
  ColS   nodesTopRightOpenEdge();           // node-numbers along the top    - right  edge
  // // boundary nodes: edges (aliases)
  // ColS   nodesBottomFrontEdge();        // alias, see above: nodesFrontBottomEdge
  // ColS   nodesBottomBackEdge();         // alias, see above: nodesBackBottomEdge
  // ColS   nodesTopFrontEdge();           // alias, see above: nodesFrontTopEdge
  // ColS   nodesTopBackEdge();            // alias, see above: nodesBackTopEdge
  // ColS   nodesLeftBottomEdge();         // alias, see above: nodesBottomLeftEdge
  // ColS   nodesLeftFrontEdge();          // alias, see above: nodesFrontLeftEdge
  // ColS   nodesLeftBackEdge();           // alias, see above: nodesBackLeftEdge
  // ColS   nodesLeftTopEdge();            // alias, see above: nodesTopLeftEdge
  // ColS   nodesRightBottomEdge();        // alias, see above: nodesBottomRightEdge
  // ColS   nodesRightTopEdge();           // alias, see above: nodesTopRightEdge
  // ColS   nodesRightFrontEdge();         // alias, see above: nodesFrontRightEdge
  // ColS   nodesRightBackEdge();          // alias, see above: nodesBackRightEdge
  // boundary nodes: corners
  size_t nodesFrontBottomLeftCorner();  // node-number of the front - bottom - left  corner
  size_t nodesFrontBottomRightCorner(); // node-number of the front - bottom - right corner
  size_t nodesFrontTopLeftCorner();     // node-number of the front - top    - left  corner
  size_t nodesFrontTopRightCorner();    // node-number of the front - top    - right corner
  size_t nodesBackBottomLeftCorner();   // node-number of the back  - bottom - left  corner
  size_t nodesBackBottomRightCorner();  // node-number of the back  - bottom - right corner
  size_t nodesBackTopLeftCorner();      // node-number of the back  - top    - left  corner
  size_t nodesBackTopRightCorner();     // node-number of the back  - top    - right corner
  // // boundary nodes: corners (aliases)
  // size_t nodesFrontLeftBottomCorner();  // alias, see above: nodesFrontBottomLeftCorner
  // size_t nodesBottomFrontLeftCorner();  // alias, see above: nodesFrontBottomLeftCorner
  // size_t nodesBottomLeftFrontCorner();  // alias, see above: nodesFrontBottomLeftCorner
  // size_t nodesLeftFrontBottomCorner();  // alias, see above: nodesFrontBottomLeftCorner
  // size_t nodesLeftBottomFrontCorner();  // alias, see above: nodesFrontBottomLeftCorner
  // size_t nodesFrontRightBottomCorner(); // alias, see above: nodesFrontBottomRightCorner
  // size_t nodesBottomFrontRightCorner(); // alias, see above: nodesFrontBottomRightCorner
  // size_t nodesBottomRightFrontCorner(); // alias, see above: nodesFrontBottomRightCorner
  // size_t nodesRightFrontBottomCorner(); // alias, see above: nodesFrontBottomRightCorner
  // size_t nodesRightBottomFrontCorner(); // alias, see above: nodesFrontBottomRightCorner
  // size_t nodesFrontLeftTopCorner();     // alias, see above: nodesFrontTopLeftCorner
  // size_t nodesTopFrontLeftCorner();     // alias, see above: nodesFrontTopLeftCorner
  // size_t nodesTopLeftFrontCorner();     // alias, see above: nodesFrontTopLeftCorner
  // size_t nodesLeftFrontTopCorner();     // alias, see above: nodesFrontTopLeftCorner
  // size_t nodesLeftTopFrontCorner();     // alias, see above: nodesFrontTopLeftCorner
  // size_t nodesFrontRightTopCorner();    // alias, see above: nodesFrontTopRightCorner
  // size_t nodesTopFrontRightCorner();    // alias, see above: nodesFrontTopRightCorner
  // size_t nodesTopRightFrontCorner();    // alias, see above: nodesFrontTopRightCorner
  // size_t nodesRightFrontTopCorner();    // alias, see above: nodesFrontTopRightCorner
  // size_t nodesRightTopFrontCorner();    // alias, see above: nodesFrontTopRightCorner
  // size_t nodesBackLeftBottomCorner();   // alias, see above: nodesBackBottomLeftCorner
  // size_t nodesBottomBackLeftCorner();   // alias, see above: nodesBackBottomLeftCorner
  // size_t nodesBottomLeftBackCorner();   // alias, see above: nodesBackBottomLeftCorner
  // size_t nodesLeftBackBottomCorner();   // alias, see above: nodesBackBottomLeftCorner
  // size_t nodesLeftBottomBackCorner();   // alias, see above: nodesBackBottomLeftCorner
  // size_t nodesBackRightBottomCorner();  // alias, see above: nodesBackBottomRightCorner
  // size_t nodesBottomBackRightCorner();  // alias, see above: nodesBackBottomRightCorner
  // size_t nodesBottomRightBackCorner();  // alias, see above: nodesBackBottomRightCorner
  // size_t nodesRightBackBottomCorner();  // alias, see above: nodesBackBottomRightCorner
  // size_t nodesRightBottomBackCorner();  // alias, see above: nodesBackBottomRightCorner
  // size_t nodesBackLeftTopCorner();      // alias, see above: nodesBackTopLeftCorner
  // size_t nodesTopBackLeftCorner();      // alias, see above: nodesBackTopLeftCorner
  // size_t nodesTopLeftBackCorner();      // alias, see above: nodesBackTopLeftCorner
  // size_t nodesLeftBackTopCorner();      // alias, see above: nodesBackTopLeftCorner
  // size_t nodesLeftTopBackCorner();      // alias, see above: nodesBackTopLeftCorner
  // size_t nodesBackRightTopCorner();     // alias, see above: nodesBackTopRightCorner
  // size_t nodesTopBackRightCorner();     // alias, see above: nodesBackTopRightCorner
  // size_t nodesTopRightBackCorner();     // alias, see above: nodesBackTopRightCorner
  // size_t nodesRightBackTopCorner();     // alias, see above: nodesBackTopRightCorner
  // size_t nodesRightTopBackCorner();     // alias, see above: nodesBackTopRightCorner
  // // periodicity
  // MatS   nodesPeriodic();               // periodic node pairs [:,2]: (independent, dependent)
  // size_t nodesOrigin();                 // front-bottom-left node, used as reference for periodicity
  // MatS   dofs();                        // DOF-numbers for each component of each node (sequential)
  // MatS   dofsPeriodic();                // ,, for the case that the periodicity if fully eliminated
};

// =================================================================================================

}}} // namespace ...

// =================================================================================================

#endif

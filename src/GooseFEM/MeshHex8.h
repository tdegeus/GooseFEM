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
  size_t m_nx;     // number of 'pixels' x-direction (length == "m_nx * m_h")
  size_t m_ny;     // number of 'pixels' y-direction (length == "m_ny * m_h")
  size_t m_nz;     // number of 'pixels' z-direction (length == "m_ny * m_h")
  double m_h;      // size of the element edge (equal in both directions)
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=8;  // number of nodes-per-element
  size_t m_ndim=3; // number of dimensions

public:
  // mesh with "nx*ny*nz" 'pixels' and edge size "h"
  Regular(size_t nx, size_t ny, size_t nz, double h=1.);
  // sizes
  size_t nelem();                       // number of elements
  size_t nnode();                       // number of nodes
  size_t nne();                         // number of nodes-per-element
  size_t ndim();                        // number of dimensions
  // mesh
  MatD   coor();                        // nodal positions [ nnode , ndim ]
  MatS   conn();                        // connectivity    [ nelem , nne  ]
  // boundary nodes: planes
  ColS   nodesBottom();                 // node-numbers along the bottom plane
  ColS   nodesTop();                    // node-numbers along the top    plane
  ColS   nodesLeft();                   // node-numbers along the left   plane
  ColS   nodesRight();                  // node-numbers along the right  plane
  ColS   nodesFront();                  // node-numbers along the front  plane
  ColS   nodesBack();                   // node-numbers along the back   plane
  // boundary nodes: faces
  ColS   nodesBottomFace();             // node-numbers along the bottom face
  ColS   nodesTopFace();                // node-numbers along the top    face
  ColS   nodesLeftFace();               // node-numbers along the left   face
  ColS   nodesRightFace();              // node-numbers along the right  face
  ColS   nodesFrontFace();              // node-numbers along the front  face
  ColS   nodesBackFace();               // node-numbers along the back   face
  // boundary nodes: edges
  ColS   nodesBottomFrontEdge();        // node-numbers along the bottom - front edge
  ColS   nodesBottomBackEdge();         // node-numbers along the bottom - back  edge
  ColS   nodesBottomLeftEdge();         // node-numbers along the bottom - left  edge
  ColS   nodesBottomRightEdge();        // node-numbers along the bottom - right edge
  ColS   nodesTopFrontEdge();           // node-numbers along the top    - front edge
  ColS   nodesTopBackEdge();            // node-numbers along the top    - back  edge
  ColS   nodesTopLeftEdge();            // node-numbers along the top    - left  edge
  ColS   nodesTopRightEdge();           // node-numbers along the top    - right edge
  ColS   nodesFrontLeftEdge();          // node-numbers along the front  - left  edge
  ColS   nodesFrontRightEdge();         // node-numbers along the front  - right edge
  ColS   nodesBackLeftEdge();           // node-numbers along the back   - left  edge
  ColS   nodesBackRightEdge();          // node-numbers along the back   - right edge
  // boundary nodes: edges (aliases)
  ColS   nodesFrontBottomEdge();        // alias, see above: nodesBottomFrontEdge
  ColS   nodesFrontTopEdge();           // alias, see above: nodesTopFrontEdge
  ColS   nodesBackBottomEdge();         // alias, see above: nodesBottomBackEdge
  ColS   nodesBackTopEdge();            // alias, see above: nodesTopBackEdge
  ColS   nodesLeftFrontEdge();          // alias, see above: nodesFrontLeftEdge
  ColS   nodesLeftBottomEdge();         // alias, see above: nodesBottomLeftEdge
  ColS   nodesLeftTopEdge();            // alias, see above: nodesTopLeftEdge
  ColS   nodesLeftBackEdge();           // alias, see above: nodesBackLeftEdge
  ColS   nodesRightFrontEdge();         // alias, see above: nodesFrontRightEdge
  ColS   nodesRightBackEdge();          // alias, see above: nodesBackRightEdge
  ColS   nodesRightBottomEdge();        // alias, see above: nodesBottomRightEdge
  ColS   nodesRightTopEdge();           // alias, see above: nodesTopRightEdge
  // boundary nodes: corners
  size_t nodesBottomFrontLeftCorner();  // node-number of the bottom - front - left  corner
  size_t nodesBottomFrontRightCorner(); // node-number of the bottom - front - right corner
  size_t nodesBottomBackLeftCorner();   // node-number of the bottom - back  - left  corner
  size_t nodesBottomBackRightCorner();  // node-number of the bottom - back  - right corner
  size_t nodesTopFrontLeftCorner();     // node-number of the top    - front - left  corner
  size_t nodesTopFrontRightCorner();    // node-number of the top    - front - right corner
  size_t nodesTopBackLeftCorner();      // node-number of the top    - back  - left  corner
  size_t nodesTopBackRightCorner();     // node-number of the top    - back  - right corner
  // boundary nodes: corners (aliases)
  size_t nodesBottomLeftFrontCorner();  // alias, see above: nodesBottomFrontLeftCorner
  size_t nodesFrontBottomLeftCorner();  // alias, see above: nodesBottomFrontLeftCorner
  size_t nodesFrontLeftBottomCorner();  // alias, see above: nodesBottomFrontLeftCorner
  size_t nodesLeftBottomFrontCorner();  // alias, see above: nodesBottomFrontLeftCorner
  size_t nodesLeftFrontBottomCorner();  // alias, see above: nodesBottomFrontLeftCorner
  size_t nodesBottomRightFrontCorner(); // alias, see above: nodesBottomFrontRightCorner
  size_t nodesFrontBottomRightCorner(); // alias, see above: nodesBottomFrontRightCorner
  size_t nodesFrontRightBottomCorner(); // alias, see above: nodesBottomFrontRightCorner
  size_t nodesRightBottomFrontCorner(); // alias, see above: nodesBottomFrontRightCorner
  size_t nodesRightFrontBottomCorner(); // alias, see above: nodesBottomFrontRightCorner
  size_t nodesBottomLeftBackCorner();   // alias, see above: nodesBottomBackLeftCorner
  size_t nodesBackBottomLeftCorner();   // alias, see above: nodesBottomBackLeftCorner
  size_t nodesBackLeftBottomCorner();   // alias, see above: nodesBottomBackLeftCorner
  size_t nodesLeftBottomBackCorner();   // alias, see above: nodesBottomBackLeftCorner
  size_t nodesLeftBackBottomCorner();   // alias, see above: nodesBottomBackLeftCorner
  size_t nodesBottomRightBackCorner();  // alias, see above: nodesBottomBackRightCorner
  size_t nodesBackBottomRightCorner();  // alias, see above: nodesBottomBackRightCorner
  size_t nodesBackRightBottomCorner();  // alias, see above: nodesBottomBackRightCorner
  size_t nodesRightBottomBackCorner();  // alias, see above: nodesBottomBackRightCorner
  size_t nodesRightBackBottomCorner();  // alias, see above: nodesBottomBackRightCorner
  size_t nodesTopLeftFrontCorner();     // alias, see above: nodesTopFrontLeftCorner
  size_t nodesFrontTopLeftCorner();     // alias, see above: nodesTopFrontLeftCorner
  size_t nodesFrontLeftTopCorner();     // alias, see above: nodesTopFrontLeftCorner
  size_t nodesLeftTopFrontCorner();     // alias, see above: nodesTopFrontLeftCorner
  size_t nodesLeftFrontTopCorner();     // alias, see above: nodesTopFrontLeftCorner
  size_t nodesTopRightFrontCorner();    // alias, see above: nodesTopFrontRightCorner
  size_t nodesFrontTopRightCorner();    // alias, see above: nodesTopFrontRightCorner
  size_t nodesFrontRightTopCorner();    // alias, see above: nodesTopFrontRightCorner
  size_t nodesRightTopFrontCorner();    // alias, see above: nodesTopFrontRightCorner
  size_t nodesRightFrontTopCorner();    // alias, see above: nodesTopFrontRightCorner
  size_t nodesTopLeftBackCorner();      // alias, see above: nodesTopBackLeftCorner
  size_t nodesBackTopLeftCorner();      // alias, see above: nodesTopBackLeftCorner
  size_t nodesBackLeftTopCorner();      // alias, see above: nodesTopBackLeftCorner
  size_t nodesLeftTopBackCorner();      // alias, see above: nodesTopBackLeftCorner
  size_t nodesLeftBackTopCorner();      // alias, see above: nodesTopBackLeftCorner
  size_t nodesTopRightBackCorner();     // alias, see above: nodesTopBackRightCorner
  size_t nodesBackTopRightCorner();     // alias, see above: nodesTopBackRightCorner
  size_t nodesBackRightTopCorner();     // alias, see above: nodesTopBackRightCorner
  size_t nodesRightTopBackCorner();     // alias, see above: nodesTopBackRightCorner
  size_t nodesRightBackTopCorner();     // alias, see above: nodesTopBackRightCorner
  // periodicity
  MatS   nodesPeriodic();   // periodic node pairs [ : , 2 ]: ( independent , dependent )
  size_t nodesOrigin();     // bottom-left node, to be used as reference for periodicity
  MatS   dofs();            // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic();    // DOF-numbers for each component of each node (sequential)
};


}}} // namespace ...

#endif

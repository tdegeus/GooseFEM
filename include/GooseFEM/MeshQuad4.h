/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHQUAD4_H
#define GOOSEFEM_MESHQUAD4_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// ===================================== GooseFEM::Mesh::Quad4 =====================================

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {

// ----------------------------------------- regular mesh ------------------------------------------

class Regular
{
private:
  double m_h;      // elementary element edge-size (in both directions)
  size_t m_nelx;   // number of elements in x-direction (length == "m_nelx * m_h")
  size_t m_nely;   // number of elements in y-direction (length == "m_nely * m_h")
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=4;  // number of nodes-per-element
  size_t m_ndim=2; // number of dimensions

public:
  // mesh with "nelx*nely" 'elements' of edge size "h"
  Regular(size_t nelx, size_t nely, double h=1.);
  // sizes
  size_t nelem() const;                   // number of elements
  size_t nne() const;                     // number of nodes-per-element
  size_t ndim() const;                    // number of dimensions
  size_t nnode() const;                   // number of nodes
  size_t nnodePeriodic() const;           // number of nodes, after elimination of periodicity
  // mesh
  MatD   coor() const;                    // nodal positions [nnode ,ndim]
  MatS   conn() const;                    // connectivity    [nelem ,nne ]
  // boundary nodes: edges
  ColS   nodesBottomEdge() const;         // node-numbers along the bottom edge
  ColS   nodesTopEdge() const;            // node-numbers along the top    edge
  ColS   nodesLeftEdge() const;           // node-numbers along the left   edge
  ColS   nodesRightEdge() const;          // node-numbers along the right  edge
  // boundary nodes: edges, without corners
  ColS   nodesBottomOpenEdge() const;     // node-numbers along the bottom edge
  ColS   nodesTopOpenEdge() const;        // node-numbers along the top    edge
  ColS   nodesLeftOpenEdge() const;       // node-numbers along the left   edge
  ColS   nodesRightOpenEdge() const;      // node-numbers along the right  edge
  // boundary nodes: corners
  size_t nodesBottomLeftCorner() const;   // node-number of the bottom - left  corner
  size_t nodesBottomRightCorner() const;  // node-number of the bottom - right corner
  size_t nodesTopLeftCorner() const;      // node-number of the top    - left  corner
  size_t nodesTopRightCorner() const;     // node-number of the top    - right corner
  // boundary nodes: corners (aliases)
  size_t nodesLeftBottomCorner() const;   // alias, see above: nodesBottomLeftCorner
  size_t nodesLeftTopCorner() const;      // alias, see above: nodesBottomRightCorner
  size_t nodesRightBottomCorner() const;  // alias, see above: nodesTopLeftCorner
  size_t nodesRightTopCorner() const;     // alias, see above: nodesTopRightCorner
  // periodicity
  MatS   nodesPeriodic() const; // periodic node pairs [:,2]: (independent, dependent)
  size_t nodesOrigin() const;   // bottom-left node, used as reference for periodicity
  MatS   dofs() const;          // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic() const;  // ,, for the case that the periodicity if fully eliminated
};

// ---------------------- mesh with a fine layer that exponentially coarsens -----------------------

class FineLayer
{
private:
  double m_h;          // elementary element edge-size (in all directions)
  double m_Lx;         // mesh size in "x"
  ColS   m_nelx;       // number of elements in "x"                     (per el.layer in "y")
  ColS   m_nnd;        // total number of nodes in the main node layer  (per nd.layer in "y")
  ColS   m_nhx, m_nhy; // element size in each direction                (per el.layer in "y")
  ColI   m_refine;     // refine direction (-1:no refine, 0:"x")        (per el.layer in "y")
  ColS   m_startElem;  // start element                                 (per el.layer in "y")
  ColS   m_startNode;  // start node                                    (per nd.layer in "y")
  size_t m_nelem;      // number of elements
  size_t m_nnode;      // number of nodes
  size_t m_nne=4;      // number of nodes-per-element
  size_t m_ndim=2;     // number of dimensions

public:
  // mesh with "nelx*nely" elements of edge size "h"; elements are coarsened in "y"-direction
  FineLayer(size_t nelx, size_t nely, double h=1., size_t nfine=1);
  // sizes
  size_t nelem() const;                   // number of elements
  size_t nnode() const;                   // number of nodes
  size_t nne() const;                     // number of nodes-per-element
  size_t ndim() const;                    // number of dimensions
  size_t shape(size_t i) const;           // actual shape in a certain direction
  // mesh
  MatD   coor() const;                    // nodal positions [nnode ,ndim]
  MatS   conn() const;                    // connectivity    [nelem ,nne ]
  // element sets
  ColS   elementsMiddleLayer() const;     // elements in the middle, fine, layer
  // boundary nodes: edges
  ColS   nodesBottomEdge() const;         // node-numbers along the bottom edge
  ColS   nodesTopEdge() const;            // node-numbers along the top    edge
  ColS   nodesLeftEdge() const;           // node-numbers along the left   edge
  ColS   nodesRightEdge() const;          // node-numbers along the right  edge
  // boundary nodes: edges, without corners
  ColS   nodesBottomOpenEdge() const;     // node-numbers along the bottom edge
  ColS   nodesTopOpenEdge() const;        // node-numbers along the top    edge
  ColS   nodesLeftOpenEdge() const;       // node-numbers along the left   edge
  ColS   nodesRightOpenEdge() const;      // node-numbers along the right  edge
  // boundary nodes: corners
  size_t nodesBottomLeftCorner() const;   // node-number of the bottom - left  corner
  size_t nodesBottomRightCorner() const;  // node-number of the bottom - right corner
  size_t nodesTopLeftCorner() const;      // node-number of the top    - left  corner
  size_t nodesTopRightCorner() const;     // node-number of the top    - right corner
  // boundary nodes: corners (aliases)
  size_t nodesLeftBottomCorner() const;   // alias, see above: nodesBottomLeftCorner
  size_t nodesLeftTopCorner() const;      // alias, see above: nodesBottomRightCorner
  size_t nodesRightBottomCorner() const;  // alias, see above: nodesTopLeftCorner
  size_t nodesRightTopCorner() const;     // alias, see above: nodesTopRightCorner
  // periodicity
  MatS   nodesPeriodic() const; // periodic node pairs [:,2]: (independent, dependent)
  size_t nodesOrigin() const;   // bottom-left node, used as reference for periodicity
  MatS   dofs() const;          // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic() const;  // ,, for the case that the periodicity if fully eliminated
};

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

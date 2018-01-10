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

// ========================================== REGULAR MESH =========================================

class Regular
{
private:
  size_t m_nx;              // number of 'pixels' horizontal direction (length == "m_nx * m_h")
  size_t m_ny;              // number of 'pixels' vertical direction   (length == "m_ny * m_h")
  double m_h;               // size of the element edge (equal in both directions)
  size_t m_nelem;           // number of elements
  size_t m_nnode;           // number of nodes
  size_t m_nne=4;           // number of nodes-per-element
  size_t m_ndim=2;          // number of dimensions

public:
  // mesh with "nx" pixels in horizontal direction, "ny" in vertical direction and "h" the edge size
  Regular(size_t nx, size_t ny, double h=1.);

  size_t nelem();           // number of elements
  size_t nnode();           // number of nodes
  size_t nne();             // number of nodes-per-element
  size_t ndim();            // number of dimensions
  MatD   coor();            // nodal positions [ nnode , ndim ]
  MatS   conn();            // connectivity    [ nelem , nne  ]
  ColS   nodesBottom();     // nodes along the bottom edge
  ColS   nodesTop();        // nodes along the top    edge
  ColS   nodesLeft();       // nodes along the left   edge
  ColS   nodesRight();      // nodes along the right  edge
  MatS   nodesPeriodic();   // periodic node pairs [ : , 2 ]: ( independent , dependent )
  size_t nodeOrigin();      // bottom-left node, to be used as reference for periodicity
  MatS   dofs();            // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic();    // DOF-numbers for each component of each node (sequential)
};

// ====================== MESH WITH A FINE LAYER THAT EXPONENTIALLY COARSENS =======================

class FineLayer
{
private:
  double m_h;               // base size of the element edge (equal in both directions)
  size_t m_nx;              // number of elements in vertical direction
  ColS   m_nh;              // element size in vertical direction (number of time "h")
  ColS   m_startNode;       // start node of each row
  ColS   m_startElem;       // start element of each row
  size_t m_nelem;           // number of elements
  size_t m_nnode;           // number of nodes
  size_t m_nne=4;           // number of nodes-per-element
  size_t m_ndim=2;          // number of dimensions

public:
  // mesh with "nx" pixels in horizontal direction, "ny" in vertical direction and "h" the edge size
  FineLayer(size_t nx, size_t ny, double h=1., size_t nfine=0, size_t nskip=0);

  size_t nelem();           // number of elements
  size_t nnode();           // number of nodes
  size_t nne();             // number of nodes-per-element
  size_t ndim();            // number of dimensions
  size_t shape(size_t i);   // actual shape in horizontal and vertical direction
  MatD   coor();            // nodal positions [ nnode , ndim ]
  MatS   conn();            // connectivity    [ nelem , nne  ]
  ColS   elementsFine();    // elements in the middle, fine, layer
  ColS   nodesBottom();     // nodes along the bottom edge
  ColS   nodesTop();        // nodes along the top    edge
  ColS   nodesLeft();       // nodes along the left   edge
  ColS   nodesRight();      // nodes along the right  edge
  MatS   nodesPeriodic();   // periodic node pairs [ : , 2 ]: ( independent , dependent )
  size_t nodeOrigin();      // bottom-left node, to be used as reference for periodicity
  MatS   dofs();            // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic();    // DOF-numbers for each component of each node (sequential)
};

// =================================================================================================

}}} // namespace ...

// =================================================================================================

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHTRI3_H
#define GOOSEFEM_MESHTRI3_H

// -------------------------------------------------------------------------------------------------

#include "Mesh.h"

// ===================================== GooseFEM::Mesh::Tri3 =====================================

namespace GooseFEM {
namespace Mesh {
namespace Tri3 {

// ========================================== REGULAR MESH =========================================

class Regular
{
private:
  size_t m_nx;              // number of 'pixels' horizontal direction (length == "m_nx * m_h")
  size_t m_ny;              // number of 'pixels' vertical direction   (length == "m_ny * m_h")
  double m_h;               // size of the element edge (equal in both directions)
  size_t m_nelem;           // number of elements
  size_t m_nnode;           // number of nodes
  size_t m_nne=3;           // number of nodes-per-element
  size_t m_ndim=2;          // number of dimensions

public:
  // mesh with "nx" pixels in horizontal direction, "ny" in vertical direction and "h" the edge size
  Regular(size_t nx, size_t ny, double h=1.);
  // sizes
  size_t nelem();           // number of elements
  size_t nnode();           // number of nodes
  size_t nne();             // number of nodes-per-element
  size_t ndim();            // number of dimensions
  // mesh
  MatD   coor();            // nodal positions [ nnode , ndim ]
  MatS   conn();            // connectivity    [ nelem , nne  ]
  // boundary nodes: edges
  ColS   nodesBottomEdge(); // nodes along the bottom edge
  ColS   nodesTopEdge();    // nodes along the top    edge
  ColS   nodesLeftEdge();   // nodes along the left   edge
  ColS   nodesRightEdge();  // nodes along the right  edge
  // boundary nodes: corners
  size_t nodesBottomLeftCorner();   // bottom - left  corner node
  size_t nodesBottomRightCorner();  // bottom - right corner node
  size_t nodesTopLeftCorner();      // top    - left  corner node
  size_t nodesTopRightCorner();     // top    - right corner node
  // boundary nodes: corners (aliases)
  size_t nodesLeftBottomCorner();   // bottom - left  corner node
  size_t nodesLeftTopCorner();      // top    - left  corner node
  size_t nodesRightBottomCorner();  // bottom - right corner node
  size_t nodesRightTopCorner();     // top    - right corner node
  // periodicity
  MatS   nodesPeriodic();   // periodic node pairs [ : , 2 ]: (independent, dependent)
  size_t nodesOrigin();     // bottom-left node, to be used as reference for periodicity
  MatS   dofs();            // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic();    // DOF-numbers for each component of each node (sequential)
};

// ========================================= MESH ANALYSIS =========================================

// read / set the orientation (-1 / +1) of all triangles
inline ColI getOrientation(const MatD &coor, const MatS &conn                                     );
inline MatS setOrientation(const MatD &coor, const MatS &conn,                  int orientation=-1);
inline MatS setOrientation(const MatD &coor, const MatS &conn, const ColI &val, int orientation=-1);

// ======================================= RE-TRIANGULATION ========================================

// simple interface to compute the full re-triangulation; it uses, depending on the input mesh:
// (1) the minimal evasive "TriUpdate"
// (2) the more rigorous "TriCompute"
inline MatS retriangulate(const MatD &coor, const MatS &conn, int orientation=-1);

// ------------------------- support class - update existing triangulation -------------------------

// minimal evasive re-triangulation which only flips edges of the existing connectivity
class TriUpdate
{
private:
  MatS   m_edge;    // the element that neighbors along each edge (m_nelem: no neighbor)
  MatS   m_conn;    // connectivity (updated)
  MatD   m_coor;    // nodal positions (does not change)
  size_t m_nelem;   // #elements
  size_t m_nnode;   // #nodes
  size_t m_nne;     // #nodes-per-element
  size_t m_ndim;    // #dimensions
  ColS   m_elem;    // the two elements involved in the last element change (see below)
  ColS   m_node;    // the four nodes   involved in the last element change (see below)
  // old: m_elem(0) = [ m_node(0) , m_node(1) , m_node(2) ]
  //      m_elem(1) = [ m_node(1) , m_node(3) , m_node(2) ]
  // new: m_elem(0) = [ m_node(0) , m_node(3) , m_node(2) ]
  //      m_elem(1) = [ m_node(0) , m_node(1) , m_node(3) ]

  // compute neighbors per edge of all elements
  void edge();

  // update edges around renumbered elements
  void chedge(size_t edge, size_t old_elem, size_t new_elem);

public:
  TriUpdate(){};
  TriUpdate(const MatD &coor, const MatS &conn);

  bool eval     ();                     // re-triangulate the full mesh (returns "true" if changed)
  bool increment();                     // one re-triangulation step    (returns "true" if changed)
  MatS conn     () { return m_conn; };  // return (new) connectivity
  MatS ch_elem  () { return m_elem; };  // return element involved in last element change
  MatS ch_node  () { return m_node; };  // return nodes   involved in last element change
};

// ================================ SUPPORT CLASS - EDGE DEFINITION ================================

// support class to allow the storage of a list of edges
class Edge {
public:
  size_t n1  ; // node 1 (edge from node 1-2)
  size_t n2  ; // node 2 (edge from node 1-2)
  size_t elem; // element to which the edge belong
  size_t edge; // edge index within the element (e.g. edge==1 -> n1=conn(0,elem), n2=conn(1,elem))

  Edge(){};
  Edge(size_t i, size_t j, size_t el, size_t ed, bool sort=false);
};

// -------------------------------------------------------------------------------------------------

// compare edges
inline bool Edge_cmp (Edge a, Edge b); // check equality
inline bool Edge_sort(Edge a, Edge b); // check if "a" is smaller than "b" in terms of node-numbers

// =================================================================================================

}}} // namespace ...

#endif

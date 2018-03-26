/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHTRI3_H
#define GOOSEFEM_MESHTRI3_H

// -------------------------------------------------------------------------------------------------

#include "Mesh.h"

// ===================================== GooseFEM::Mesh::Tri3 ======================================

namespace GooseFEM {
namespace Mesh {
namespace Tri3 {

// ------------------------------------------ regular mesh -----------------------------------------

class Regular
{
private:
  double m_h;      // elementary element edge-size (in both directions)
  size_t m_nelx;   // number of elements in x-direction (length == "m_nelx * m_h")
  size_t m_nely;   // number of elements in y-direction (length == "m_nely * m_h")
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=3;  // number of nodes-per-element
  size_t m_ndim=2; // number of dimensions

public:
  // mesh with "2*nelx*nely" 'elements' of edge size "h"
  Regular(size_t nelx, size_t nely, double h=1.);
  // sizes
  size_t nelem();                   // number of elements
  size_t nnode();                   // number of nodes
  size_t nne();                     // number of nodes-per-element
  size_t ndim();                    // number of dimensions
  // mesh
  MatD   coor();                    // nodal positions [nnode ,ndim]
  MatS   conn();                    // connectivity    [nelem ,nne ]
  // boundary nodes: edges
  ColS   nodesBottomEdge();         // node-numbers along the bottom edge
  ColS   nodesTopEdge();            // node-numbers along the top    edge
  ColS   nodesLeftEdge();           // node-numbers along the left   edge
  ColS   nodesRightEdge();          // node-numbers along the right  edge
  // boundary nodes: edges, without corners
  ColS   nodesBottomOpenEdge();     // node-numbers along the bottom edge
  ColS   nodesTopOpenEdge();        // node-numbers along the top    edge
  ColS   nodesLeftOpenEdge();       // node-numbers along the left   edge
  ColS   nodesRightOpenEdge();      // node-numbers along the right  edge
  // boundary nodes: corners
  size_t nodesBottomLeftCorner();   // node-number of the bottom - left  corner
  size_t nodesBottomRightCorner();  // node-number of the bottom - right corner
  size_t nodesTopLeftCorner();      // node-number of the top    - left  corner
  size_t nodesTopRightCorner();     // node-number of the top    - right corner
  // boundary nodes: corners (aliases)
  size_t nodesLeftBottomCorner();   // alias, see above: nodesBottomLeftCorner
  size_t nodesLeftTopCorner();      // alias, see above: nodesBottomRightCorner
  size_t nodesRightBottomCorner();  // alias, see above: nodesTopLeftCorner
  size_t nodesRightTopCorner();     // alias, see above: nodesTopRightCorner
  // periodicity
  MatS   nodesPeriodic();           // periodic node pairs [:,2]: (independent, dependent)
  size_t nodesOrigin();             // bottom-left node, used as reference for periodicity
  MatS   dofs();                    // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic();            // ,, for the case that the periodicity if fully eliminated
};

// ----------------------------------------- mesh analysis -----------------------------------------

// read / set the orientation (-1 / +1) of all triangles
inline ColI getOrientation(const MatD &coor, const MatS &conn                                     );
inline MatS setOrientation(const MatD &coor, const MatS &conn,                  int orientation=-1);
inline MatS setOrientation(const MatD &coor, const MatS &conn, const ColI &val, int orientation=-1);

// --------------------------------------- re-triangulation ----------------------------------------

// simple interface to compute the full re-triangulation; it uses, depending on the input mesh:
// (1) the minimal evasive "TriUpdate"
// (2) the more rigorous "TriCompute"
inline MatS retriangulate(const MatD &coor, const MatS &conn, int orientation=-1);



// ================================= GooseFEM::Mesh::Tri3::Private =================================

namespace Private {

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

// -------------------------------- support class - edge definition --------------------------------

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

// -------------------------------------------------------------------------------------------------

} // namespace Private

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHTRI3_H
#define GOOSEFEM_MESHTRI3_H

// -------------------------------------------------------------------------------------------------

#include "config.h"

// =================================================================================================

namespace GooseFEM {
namespace Mesh {
namespace Tri3 {

// -------------------------------------------------------------------------------------------------

class Regular
{
private:
  double m_h;                   // elementary element edge-size (in both directions)
  size_t m_nelx;                // number of elements in x-direction (length == "m_nelx * m_h")
  size_t m_nely;                // number of elements in y-direction (length == "m_nely * m_h")
  size_t m_nelem;               // number of elements
  size_t m_nnode;               // number of nodes
  static const size_t m_nne=3;  // number of nodes-per-element
  static const size_t m_ndim=2; // number of dimensions

public:
  // mesh with "2*nelx*nely" 'elements' of edge size "h"
  Regular(size_t nelx, size_t nely, double h=1.);
  // sizes
  size_t nelem() const; // number of elements
  size_t nnode() const; // number of nodes
  size_t nne() const;   // number of nodes-per-element
  size_t ndim() const;  // number of dimensions
  // mesh
  xt::xtensor<double,2> coor() const;                    // nodal positions [nnode ,ndim]
  xt::xtensor<size_t,2> conn() const;                    // connectivity    [nelem ,nne ]
  // boundary nodes: edges
  xt::xtensor<size_t,1> nodesBottomEdge() const;         // node-numbers along the bottom edge
  xt::xtensor<size_t,1> nodesTopEdge() const;            // node-numbers along the top    edge
  xt::xtensor<size_t,1> nodesLeftEdge() const;           // node-numbers along the left   edge
  xt::xtensor<size_t,1> nodesRightEdge() const;          // node-numbers along the right  edge
  // boundary nodes: edges, without corners
  xt::xtensor<size_t,1> nodesBottomOpenEdge() const;     // node-numbers along the bottom edge
  xt::xtensor<size_t,1> nodesTopOpenEdge() const;        // node-numbers along the top    edge
  xt::xtensor<size_t,1> nodesLeftOpenEdge() const;       // node-numbers along the left   edge
  xt::xtensor<size_t,1> nodesRightOpenEdge() const;      // node-numbers along the right  edge
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
  xt::xtensor<size_t,2> nodesPeriodic() const; // periodic node pairs [:,2]: (independent, dependent)
  size_t                nodesOrigin() const;   // bottom-left node, used as reference for periodicity
  xt::xtensor<size_t,2> dofs() const;          // DOF-numbers for each component of each node (sequential)
  xt::xtensor<size_t,2> dofsPeriodic() const;  // ,, for the case that the periodicity if fully eliminated
};

// -------------------------------------------------------------------------------------------------

// read / set the orientation (-1 / +1) of all triangles
inline xt::xtensor<int   ,1> getOrientation(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn                                                   );
inline xt::xtensor<size_t,2> setOrientation(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn,                                int orientation=-1);
inline xt::xtensor<size_t,2> setOrientation(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn, const xt::xtensor<int,1> &val, int orientation=-1);

// -------------------------------------------------------------------------------------------------

// simple interface to compute the full re-triangulation; it uses, depending on the input mesh:
// (1) the minimal evasive "TriUpdate"
// (2) the more rigorous "TriCompute"
inline xt::xtensor<size_t,2> retriangulate(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn, int orientation=-1);

// =================================================================================================

namespace Private {

// -------------------------------------------------------------------------------------------------

// minimal evasive re-triangulation which only flips edges of the existing connectivity
class TriUpdate
{
private:
  size_t                m_nelem; // #elements
  size_t                m_nnode; // #nodes
  size_t                m_nne;   // #nodes-per-element
  size_t                m_ndim;  // #dimensions
  xt::xtensor<size_t,2> m_edge;  // the element that neighbors along each edge (m_nelem: no neighbor)
  xt::xtensor<size_t,2> m_conn;  // connectivity (updated)
  xt::xtensor<double,2> m_coor;  // nodal positions (does not change)
  xt::xtensor_fixed<size_t,xt::xshape<2>> m_elem; // the two elements in the last element change
  xt::xtensor_fixed<size_t,xt::xshape<4>> m_node; // the four nodes   in the last element change
  // old: m_elem(0) = [ m_node(0) , m_node(1) , m_node(2) ]
  //      m_elem(1) = [ m_node(1) , m_node(3) , m_node(2) ]
  // new: m_elem(0) = [ m_node(0) , m_node(3) , m_node(2) ]
  //      m_elem(1) = [ m_node(0) , m_node(1) , m_node(3) ]

  // compute neighbors per edge of all elements
  void edge();

  // update edges around renumbered elements
  void chedge(size_t edge, size_t old_elem, size_t new_elem);

public:
  TriUpdate() = default;
  TriUpdate(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn);

  bool eval     (); // re-triangulate the full mesh (returns "true" if changed)
  bool increment(); // one re-triangulation step    (returns "true" if changed)
  xt::xtensor<size_t,2> conn()    { return m_conn; };  // (new) connectivity
  xt::xtensor<size_t,2> ch_elem() { return m_elem; };  // element involved in last element change
  xt::xtensor<size_t,2> ch_node() { return m_node; };  // nodes   involved in last element change
};

// -------------------------------------------------------------------------------------------------

// support class to allow the storage of a list of edges
class Edge {
public:
  size_t n1  ; // node 1 (edge from node 1-2)
  size_t n2  ; // node 2 (edge from node 1-2)
  size_t elem; // element to which the edge belong
  size_t edge; // edge index within the element (e.g. edge==1 -> n1=conn(0,elem), n2=conn(1,elem))

  Edge() = default;
  Edge(size_t i, size_t j, size_t el, size_t ed, bool sort=false);
};

// -------------------------------------------------------------------------------------------------

// compare edges
inline bool Edge_cmp (Edge a, Edge b); // check equality
inline bool Edge_sort(Edge a, Edge b); // check if "a" is smaller than "b" in terms of node-numbers

// -------------------------------------------------------------------------------------------------

} // namespace Private

// =================================================================================================

}}} // namespace ...

// =================================================================================================

#include "MeshTri3.hpp"

// =================================================================================================

#endif

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
  // type
  ElementType getElementType() const;
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

// read / set the orientation (-1/+1) of all triangles
inline xt::xtensor<int,1> getOrientation(
  const xt::xtensor<double,2> &coor,  // nodal coordinates
  const xt::xtensor<size_t,2> &conn); // connectivity

inline xt::xtensor<size_t,2> setOrientation(
  const xt::xtensor<double,2> &coor,  // nodal coordinates
  const xt::xtensor<size_t,2> &conn,  // connectivity
  int orientation=-1);                // target orientation (-1/+1)

inline xt::xtensor<size_t,2> setOrientation(
  const xt::xtensor<double,2> &coor,  // nodal coordinates
  const xt::xtensor<size_t,2> &conn,  // connectivity
  const xt::xtensor<int,1> &val,      // current orientations (== output of "getOrientation")
  int orientation=-1);                // target orientation (-1/+1)

// =================================================================================================

}}} // namespace ...

// =================================================================================================

#include "MeshTri3.hpp"

// =================================================================================================

#endif

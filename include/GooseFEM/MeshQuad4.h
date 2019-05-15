/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHQUAD4_H
#define GOOSEFEM_MESHQUAD4_H

// -------------------------------------------------------------------------------------------------

#include "config.h"

// =================================================================================================

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {
namespace Map {

class FineLayer2Regular;

}}}}

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {

class Regular
{
public:

  Regular() = default;

  // mesh with "nelx*nely" 'elements' of edge size "h"
  Regular(size_t nelx, size_t nely, double h=1.);

  // size
  size_t nelem() const; // number of elements
  size_t nnode() const; // number of nodes
  size_t nne() const;   // number of nodes-per-element
  size_t ndim() const;  // number of dimensions
  size_t nelx() const;  // number of elements in x-direction
  size_t nely() const;  // number of elements in y-direction
  double h() const;     // edge size

  // type
  ElementType getElementType() const;

  // mesh
  xt::xtensor<double,2> coor() const; // nodal positions [nnode ,ndim]
  xt::xtensor<size_t,2> conn() const; // connectivity    [nelem ,nne ]

  // boundary nodes: edges
  xt::xtensor<size_t,1> nodesBottomEdge() const; // node-numbers along the bottom edge
  xt::xtensor<size_t,1> nodesTopEdge() const;    // node-numbers along the top    edge
  xt::xtensor<size_t,1> nodesLeftEdge() const;   // node-numbers along the left   edge
  xt::xtensor<size_t,1> nodesRightEdge() const;  // node-numbers along the right  edge

  // boundary nodes: edges, without corners
  xt::xtensor<size_t,1> nodesBottomOpenEdge() const; // node-numbers along the bottom edge
  xt::xtensor<size_t,1> nodesTopOpenEdge() const;    // node-numbers along the top    edge
  xt::xtensor<size_t,1> nodesLeftOpenEdge() const;   // node-numbers along the left   edge
  xt::xtensor<size_t,1> nodesRightOpenEdge() const;  // node-numbers along the right  edge

  // boundary nodes: corners (including aliases)
  size_t nodesBottomLeftCorner() const;  // node-number of the bottom - left  corner
  size_t nodesBottomRightCorner() const; // node-number of the bottom - right corner
  size_t nodesTopLeftCorner() const;     // node-number of the top    - left  corner
  size_t nodesTopRightCorner() const;    // node-number of the top    - right corner
  size_t nodesLeftBottomCorner() const;  // alias, see above: nodesBottomLeftCorner
  size_t nodesLeftTopCorner() const;     // alias, see above: nodesBottomRightCorner
  size_t nodesRightBottomCorner() const; // alias, see above: nodesTopLeftCorner
  size_t nodesRightTopCorner() const;    // alias, see above: nodesTopRightCorner

  // DOFs
  xt::xtensor<size_t,2> dofs() const; // DOF-numbers for each component of each node (sequential)

  // periodicity
  xt::xtensor<size_t,2> nodesPeriodic() const; // periodic node pairs [:,2]: (independent, dependent)
  size_t                nodesOrigin() const;   // bottom-left node, used as reference for periodicity
  xt::xtensor<size_t,2> dofsPeriodic() const;  // DOF-numbers when the periodicity if fully eliminated

  // element numbers as matrix
  xt::xtensor<size_t,2> elementMatrix() const;


private:

  double m_h;                   // elementary element edge-size (in all directions)
  size_t m_nelx;                // number of elements in x-direction (length == "m_nelx * m_h")
  size_t m_nely;                // number of elements in y-direction (length == "m_nely * m_h")
  size_t m_nelem;               // number of elements
  size_t m_nnode;               // number of nodes
  static const size_t m_nne=4;  // number of nodes-per-element
  static const size_t m_ndim=2; // number of dimensions
};

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {

class FineLayer
{
public:

  FineLayer() = default;

  // mesh with "nelx*nely" elements of edge size "h"; elements are coarsened in "y"-direction
  FineLayer(size_t nelx, size_t nely, double h=1., size_t nfine=1);

  // sizes
  size_t nelem() const; // number of elements
  size_t nnode() const; // number of nodes
  size_t nne() const;   // number of nodes-per-element
  size_t ndim() const;  // number of dimensions
  size_t nelx() const;  // number of elements in x-direction
  size_t nely() const;  // number of elements in y-direction
  double h() const;     // edge size

  // type
  ElementType getElementType() const;

  // mesh
  xt::xtensor<double,2> coor() const; // nodal positions [nnode ,ndim]
  xt::xtensor<size_t,2> conn() const; // connectivity    [nelem ,nne ]

  // element sets
  xt::xtensor<size_t,1> elementsMiddleLayer() const; // elements in the middle, fine, layer

  // boundary nodes: edges
  xt::xtensor<size_t,1> nodesBottomEdge() const; // node-numbers along the bottom edge
  xt::xtensor<size_t,1> nodesTopEdge() const;    // node-numbers along the top    edge
  xt::xtensor<size_t,1> nodesLeftEdge() const;   // node-numbers along the left   edge
  xt::xtensor<size_t,1> nodesRightEdge() const;  // node-numbers along the right  edge

  // boundary nodes: edges, without corners
  xt::xtensor<size_t,1> nodesBottomOpenEdge() const; // node-numbers along the bottom edge
  xt::xtensor<size_t,1> nodesTopOpenEdge() const;    // node-numbers along the top    edge
  xt::xtensor<size_t,1> nodesLeftOpenEdge() const;   // node-numbers along the left   edge
  xt::xtensor<size_t,1> nodesRightOpenEdge() const;  // node-numbers along the right  edge

  // boundary nodes: corners (including aliases)
  size_t nodesBottomLeftCorner() const;  // node-number of the bottom - left  corner
  size_t nodesBottomRightCorner() const; // node-number of the bottom - right corner
  size_t nodesTopLeftCorner() const;     // node-number of the top    - left  corner
  size_t nodesTopRightCorner() const;    // node-number of the top    - right corner
  size_t nodesLeftBottomCorner() const;  // alias, see above: nodesBottomLeftCorner
  size_t nodesLeftTopCorner() const;     // alias, see above: nodesBottomRightCorner
  size_t nodesRightBottomCorner() const; // alias, see above: nodesTopLeftCorner
  size_t nodesRightTopCorner() const;    // alias, see above: nodesTopRightCorner

  // DOFs
  xt::xtensor<size_t,2> dofs() const; // DOF-numbers for each component of each node (sequential)

  // periodicity
  xt::xtensor<size_t,2> nodesPeriodic() const; // periodic node pairs [:,2]: (independent, dependent)
  size_t                nodesOrigin() const;   // bottom-left node, used as reference for periodicity
  xt::xtensor<size_t,2> dofsPeriodic() const;  // DOF-numbers when the periodicity if fully eliminated

private:

  double m_h;                         // elementary element edge-size (in all directions)
  double m_Lx;                        // mesh size in "x"
  size_t m_nelem;                     // number of elements
  size_t m_nnode;                     // number of nodes
  static const size_t m_nne=4;        // number of nodes-per-element
  static const size_t m_ndim=2;       // number of dimensions
  xt::xtensor<size_t,1> m_nelx;       // number of elements in "x"                    (per el.layer in "y")
  xt::xtensor<size_t,1> m_nnd;        // total number of nodes in the main node layer (per nd.layer in "y")
  xt::xtensor<size_t,1> m_nhx, m_nhy; // element size in each direction               (per el.layer in "y")
  xt::xtensor<int   ,1> m_refine;     // refine direction (-1:no refine, 0:"x")       (per el.layer in "y")
  xt::xtensor<size_t,1> m_startElem;  // start element                                (per el.layer in "y")
  xt::xtensor<size_t,1> m_startNode;  // start node                                   (per nd.layer in "y")

  friend class GooseFEM::Mesh::Quad4::Map::FineLayer2Regular;
};

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {
namespace Map {

class RefineRegular
{
public:

  // constructor
  RefineRegular() = default;
  RefineRegular(const GooseFEM::Mesh::Quad4::Regular &mesh, size_t nx, size_t ny);

  // return the one of the two meshes
  GooseFEM::Mesh::Quad4::Regular getCoarseMesh() const;
  GooseFEM::Mesh::Quad4::Regular getFineMesh() const;

  // elements of the Fine mesh per element of the Coarse mesh
  xt::xtensor<size_t,2> getMap() const;

  // map field
  xt::xtensor<double,2> mapToCoarse(const xt::xtensor<double,1> &data) const; // scalar per element
  xt::xtensor<double,2> mapToCoarse(const xt::xtensor<double,2> &data) const; // scalar per int.pnt.
  xt::xtensor<double,4> mapToCoarse(const xt::xtensor<double,4> &data) const; // tensor per int.pnt.

  // map field
  xt::xtensor<double,1> mapToFine(const xt::xtensor<double,1> &data) const; // scalar per element
  xt::xtensor<double,2> mapToFine(const xt::xtensor<double,2> &data) const; // scalar per int.pnt.
  xt::xtensor<double,4> mapToFine(const xt::xtensor<double,4> &data) const; // tensor per int.pnt.

private:

  // the meshes
  GooseFEM::Mesh::Quad4::Regular m_coarse;
  GooseFEM::Mesh::Quad4::Regular m_fine;

  // mapping
  xt::xtensor<size_t,1> m_fine2coarse;
  xt::xtensor<size_t,1> m_fine2coarse_index;
  xt::xtensor<size_t,2> m_coarse2fine;

};

}}}} // namespace ...

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {
namespace Map {

class FineLayer2Regular
{
public:

  // constructor
  FineLayer2Regular() = default;
  FineLayer2Regular(const GooseFEM::Mesh::Quad4::FineLayer &mesh);

  // return either of the meshes
  GooseFEM::Mesh::Quad4::Regular   getRegularMesh()   const;
  GooseFEM::Mesh::Quad4::FineLayer getFineLayerMesh() const;

  // elements of the Regular mesh per element of the FineLayer mesh
  // and the fraction by which the overlap is
  std::vector<std::vector<size_t>> getMap() const;
  std::vector<std::vector<double>> getMapFraction() const;

  // map field
  xt::xtensor<double,1> mapToRegular(const xt::xtensor<double,1> &data) const; // scalar per element
  xt::xtensor<double,2> mapToRegular(const xt::xtensor<double,2> &data) const; // scalar per int.pnt.
  xt::xtensor<double,4> mapToRegular(const xt::xtensor<double,4> &data) const; // tensor per int.pnt.

private:

  // the "FineLayer" mesh to map
  GooseFEM::Mesh::Quad4::FineLayer m_finelayer;

  // the new "Regular" mesh to which to map
  GooseFEM::Mesh::Quad4::Regular m_regular;

  // mapping
  std::vector<std::vector<size_t>> m_elem_regular;
  std::vector<std::vector<double>> m_frac_regular;

};

}}}} // namespace ...

// -------------------------------------------------------------------------------------------------

// =================================================================================================

#include "MeshQuad4.hpp"

// =================================================================================================

#endif

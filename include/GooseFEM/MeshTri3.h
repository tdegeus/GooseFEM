/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_MESHTRI3_H
#define GOOSEFEM_MESHTRI3_H

#include "config.h"

namespace GooseFEM {
namespace Mesh {
namespace Tri3 {

class Regular {
public:
    Regular(size_t nelx, size_t nely, double h = 1.);

    // size
    size_t nelem() const; // number of elements
    size_t nnode() const; // number of nodes
    size_t nne() const;   // number of nodes-per-element
    size_t ndim() const;  // number of dimensions

    // type
    ElementType getElementType() const;

    // mesh
    xt::xtensor<double,2> coor() const; // nodal positions [nnode, ndim]
    xt::xtensor<size_t,2> conn() const; // connectivity [nelem, nne]

    // boundary nodes: edges
    xt::xtensor<size_t,1> nodesBottomEdge() const;
    xt::xtensor<size_t,1> nodesTopEdge() const;
    xt::xtensor<size_t,1> nodesLeftEdge() const;
    xt::xtensor<size_t,1> nodesRightEdge() const;

    // boundary nodes: edges, without corners
    xt::xtensor<size_t,1> nodesBottomOpenEdge() const;
    xt::xtensor<size_t,1> nodesTopOpenEdge() const;
    xt::xtensor<size_t,1> nodesLeftOpenEdge() const;
    xt::xtensor<size_t,1> nodesRightOpenEdge() const;

    // boundary nodes: corners (including aliases)
    size_t nodesBottomLeftCorner() const;
    size_t nodesBottomRightCorner() const;
    size_t nodesTopLeftCorner() const;
    size_t nodesTopRightCorner() const;
    size_t nodesLeftBottomCorner() const;
    size_t nodesLeftTopCorner() const;
    size_t nodesRightBottomCorner() const;
    size_t nodesRightTopCorner() const;

    // DOF-numbers for each component of each node (sequential)
    xt::xtensor<size_t,2> dofs() const;

    // DOF-numbers for the case that the periodicity if fully eliminated
    xt::xtensor<size_t,2> dofsPeriodic() const;

    // periodic node pairs [:,2]: (independent, dependent)
    xt::xtensor<size_t,2> nodesPeriodic() const;

    // front-bottom-left node, used as reference for periodicity
    size_t nodesOrigin() const;

private:
    double m_h;                     // elementary element edge-size (in all directions)
    size_t m_nelx;                  // number of elements in x-direction (length == "m_nelx * m_h")
    size_t m_nely;                  // number of elements in y-direction (length == "m_nely * m_h")
    size_t m_nelem;                 // number of elements
    size_t m_nnode;                 // number of nodes
    static const size_t m_nne = 3;  // number of nodes-per-element
    static const size_t m_ndim = 2; // number of dimensions
};

// read / set the orientation (-1/+1) of all triangles
inline xt::xtensor<int,1> getOrientation(
    const xt::xtensor<double,2>& coor,
    const xt::xtensor<size_t,2>& conn);

inline xt::xtensor<size_t,2> setOrientation(
    const xt::xtensor<double,2>& coor,
    const xt::xtensor<size_t,2>& conn,
    int orientation = -1);

inline xt::xtensor<size_t,2> setOrientation(
    const xt::xtensor<double,2>& coor,
    const xt::xtensor<size_t,2>& conn,
    const xt::xtensor<int,1>& current, // (output of "getOrientation")
    int orientation = -1);

} // namespace Tri3
} // namespace Mesh
} // namespace GooseFEM

#include "MeshTri3.hpp"

#endif

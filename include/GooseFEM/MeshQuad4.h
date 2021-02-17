/**
Generate mesh with 4-noded quadrilateral elements.

\file MeshQuad4.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESHQUAD4_H
#define GOOSEFEM_MESHQUAD4_H

#include "config.h"

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {

// pre-allocation

namespace Map {
    class FineLayer2Regular;
}

// Regular mesh: equi-sized elements

class Regular {
public:
    Regular() = default;
    Regular(size_t nelx, size_t nely, double h = 1.0);

    // size
    size_t nelem() const; // number of elements
    size_t nnode() const; // number of nodes
    size_t nne() const;   // number of nodes-per-element
    size_t ndim() const;  // number of dimensions
    size_t nelx() const;  // number of elements in x-direction
    size_t nely() const;  // number of elements in y-direction
    double h() const;     // edge size

    /**
    Get the element type.

    \return GooseFEM::Mesh::ElementType().
    */
    ElementType getElementType() const;

    // mesh
    xt::xtensor<double, 2> coor() const; // nodal positions [nnode, ndim]
    xt::xtensor<size_t, 2> conn() const; // connectivity [nelem, nne]

    // boundary nodes: edges
    xt::xtensor<size_t, 1> nodesBottomEdge() const;
    xt::xtensor<size_t, 1> nodesTopEdge() const;
    xt::xtensor<size_t, 1> nodesLeftEdge() const;
    xt::xtensor<size_t, 1> nodesRightEdge() const;

    // boundary nodes: edges, without corners
    xt::xtensor<size_t, 1> nodesBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightOpenEdge() const;

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
    xt::xtensor<size_t, 2> dofs() const;

    // DOF-numbers for the case that the periodicity if fully eliminated
    xt::xtensor<size_t, 2> dofsPeriodic() const;

    // periodic node pairs [:,2]: (independent, dependent)
    xt::xtensor<size_t, 2> nodesPeriodic() const;

    // front-bottom-left node, used as reference for periodicity
    size_t nodesOrigin() const;

    // element numbers as matrix
    xt::xtensor<size_t, 2> elementgrid() const;

private:
    double m_h;                     // elementary element edge-size (in all directions)
    size_t m_nelx;                  // number of elements in x-direction (length == "m_nelx * m_h")
    size_t m_nely;                  // number of elements in y-direction (length == "m_nely * m_h")
    size_t m_nelem;                 // number of elements
    size_t m_nnode;                 // number of nodes
    static const size_t m_nne = 4;  // number of nodes-per-element
    static const size_t m_ndim = 2; // number of dimensions
};

// Mesh with fine middle layer, and coarser elements towards the top and bottom

class FineLayer {
public:
    FineLayer() = default;
    FineLayer(size_t nelx, size_t nely, double h = 1.0, size_t nfine = 1);

    // Reconstruct class for given coordinates / connectivity
    FineLayer(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn);

    // size
    size_t nelem() const; // number of elements
    size_t nnode() const; // number of nodes
    size_t nne() const;   // number of nodes-per-element
    size_t ndim() const;  // number of dimensions
    size_t nelx() const;  // number of elements in x-direction
    size_t nely() const;  // number of elements in y-direction
    double h() const;     // edge size

    // edge size, per row of elements (in units of "h")
    xt::xtensor<size_t, 1> elemrow_nhx() const;
    xt::xtensor<size_t, 1> elemrow_nhy() const;
    xt::xtensor<size_t, 1> elemrow_nelem() const;

    // type
    ElementType getElementType() const;

    // mesh
    xt::xtensor<double, 2> coor() const; // nodal positions [nnode, ndim]
    xt::xtensor<size_t, 2> conn() const; // connectivity [nelem, nne]

    // elements in the middle (fine) layer
    xt::xtensor<size_t, 1> elementsMiddleLayer() const;

    // extract elements along a layer
    xt::xtensor<size_t, 1> elementsLayer(size_t layer) const;

    // select region of elements from 'matrix' of element numbers
    xt::xtensor<size_t, 1> elementgrid_ravel(
        std::vector<size_t> rows_start_stop,
        std::vector<size_t> cols_start_stop) const;

    /**
    Select region of elements from 'matrix' of element numbers around an element:
    square box with edge-size ``(2 * size + 1) * h``, around ``element``.

    \param element The element around which to select elements.
    \param size Edge size of the square box encapsulating the selected element.
    \param periodic Assume the mesh periodic.
    \returns List of elements.
    */
    xt::xtensor<size_t, 1> elementgrid_around_ravel(
        size_t element,
        size_t size,
        bool periodic = true);

    /**
    Select region of elements from 'matrix' of element numbers around an element:
    left/right from ``element`` (on the same layer).

    \param element The element around which to select elements.
    \param left Number of elements to select to the left.
    \param right Number of elements to select to the right.
    \param periodic Assume the mesh periodic.
    \returns List of elements.
    */
    // -
    xt::xtensor<size_t, 1> elementgrid_leftright(
        size_t element,
        size_t left,
        size_t right,
        bool periodic = true);

    // boundary nodes: edges
    xt::xtensor<size_t, 1> nodesBottomEdge() const;
    xt::xtensor<size_t, 1> nodesTopEdge() const;
    xt::xtensor<size_t, 1> nodesLeftEdge() const;
    xt::xtensor<size_t, 1> nodesRightEdge() const;

    // boundary nodes: edges, without corners
    xt::xtensor<size_t, 1> nodesBottomOpenEdge() const;
    xt::xtensor<size_t, 1> nodesTopOpenEdge() const;
    xt::xtensor<size_t, 1> nodesLeftOpenEdge() const;
    xt::xtensor<size_t, 1> nodesRightOpenEdge() const;

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
    xt::xtensor<size_t, 2> dofs() const;

    // DOF-numbers for the case that the periodicity if fully eliminated
    xt::xtensor<size_t, 2> dofsPeriodic() const;

    // periodic node pairs [:,2]: (independent, dependent)
    xt::xtensor<size_t, 2> nodesPeriodic() const;

    // front-bottom-left node, used as reference for periodicity
    size_t nodesOrigin() const;

    // mapping to 'roll' periodically in the x-direction,
    // returns element mapping, such that: new_elemvar = elemvar[elem_map]
    xt::xtensor<size_t, 1> roll(size_t n);

private:
    double m_h;                         // elementary element edge-size (in all directions)
    double m_Lx;                        // mesh size in "x"
    size_t m_nelem;                     // number of elements
    size_t m_nnode;                     // number of nodes
    static const size_t m_nne = 4;      // number of nodes-per-element
    static const size_t m_ndim = 2;     // number of dimensions
    xt::xtensor<size_t, 1> m_nelx;      // number of elements in "x" (*)
    xt::xtensor<size_t, 1> m_nnd;       // total number of nodes in the main node layer (**)
    xt::xtensor<size_t, 1> m_nhx;       // element size in x-direction (*)
    xt::xtensor<size_t, 1> m_nhy;       // element size in y-direction (*)
    xt::xtensor<int, 1> m_refine;       // refine direction (-1:no refine, 0:"x" (*)
    xt::xtensor<size_t, 1> m_startElem; // start element (*)
    xt::xtensor<size_t, 1> m_startNode; // start node (**)
    // (*) per element layer in "y"
    // (**) per node layer in "y"

    void init(size_t nelx, size_t nely, double h, size_t nfine = 1);
    void map(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn);

    friend class GooseFEM::Mesh::Quad4::Map::FineLayer2Regular;
};

// Mesh mappings

namespace Map {

    // Return "FineLayer"-class responsible for generating a connectivity
    // Throws if conversion is not possible

    GooseFEM::Mesh::Quad4::FineLayer FineLayer(
        const xt::xtensor<double, 2>& coor,
        const xt::xtensor<size_t, 2>& conn);

    // Refine a regular mesh: sub-divide elements in several smaller elements

    class RefineRegular {
    public:
        // Constructors
        RefineRegular() = default;
        RefineRegular(const GooseFEM::Mesh::Quad4::Regular& mesh, size_t nx, size_t ny);

        // return the coarse or the fine mesh objects
        GooseFEM::Mesh::Quad4::Regular getCoarseMesh() const;
        GooseFEM::Mesh::Quad4::Regular getFineMesh() const;

        // elements of the fine mesh per element of the coarse mesh
        xt::xtensor<size_t, 2> getMap() const;

        // map field
        xt::xtensor<double, 2> mapToCoarse(const xt::xtensor<double, 1>& data) const; // scalar per el
        xt::xtensor<double, 2> mapToCoarse(const xt::xtensor<double, 2>& data) const; // scalar per intpnt
        xt::xtensor<double, 4> mapToCoarse(const xt::xtensor<double, 4>& data) const; // tensor per intpnt

        // map field
        xt::xtensor<double, 1> mapToFine(const xt::xtensor<double, 1>& data) const; // scalar per el
        xt::xtensor<double, 2> mapToFine(const xt::xtensor<double, 2>& data) const; // scalar per intpnt
        xt::xtensor<double, 4> mapToFine(const xt::xtensor<double, 4>& data) const; // tensor per intpnt

    private:
        // the meshes
        GooseFEM::Mesh::Quad4::Regular m_coarse;
        GooseFEM::Mesh::Quad4::Regular m_fine;

        // mapping
        xt::xtensor<size_t, 1> m_fine2coarse;
        xt::xtensor<size_t, 1> m_fine2coarse_index;
        xt::xtensor<size_t, 2> m_coarse2fine;
    };

    class FineLayer2Regular {
    public:
        // constructor
        FineLayer2Regular() = default;
        FineLayer2Regular(const GooseFEM::Mesh::Quad4::FineLayer& mesh);

        // return either of the meshes
        GooseFEM::Mesh::Quad4::Regular getRegularMesh() const;
        GooseFEM::Mesh::Quad4::FineLayer getFineLayerMesh() const;

        // elements of the Regular mesh per element of the FineLayer mesh
        // and the fraction by which the overlap is
        std::vector<std::vector<size_t>> getMap() const;
        std::vector<std::vector<double>> getMapFraction() const;

        // map field
        xt::xtensor<double, 1> mapToRegular(const xt::xtensor<double, 1>& data) const; // scalar per el
        xt::xtensor<double, 2> mapToRegular(const xt::xtensor<double, 2>& data) const; // scalar per intpnt
        xt::xtensor<double, 4> mapToRegular(const xt::xtensor<double, 4>& data) const; // tensor per intpnt

    private:
        // the "FineLayer" mesh to map
        GooseFEM::Mesh::Quad4::FineLayer m_finelayer;

        // the new "Regular" mesh to which to map
        GooseFEM::Mesh::Quad4::Regular m_regular;

        // mapping
        std::vector<std::vector<size_t>> m_elem_regular;
        std::vector<std::vector<double>> m_frac_regular;
    };

} // namespace Map

} // namespace Quad4
} // namespace Mesh
} // namespace GooseFEM

#include "MeshQuad4.hpp"

#endif

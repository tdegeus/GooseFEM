/**
Generate simple meshes of 4-noded quadrilateral elements in 2d (GooseFEM::Mesh::ElementType::Quad4).

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

/**
Regular mesh: equi-sized elements.
*/
class Regular {
public:

    Regular() = default;

    /**
    Constructor.

    \param nelx Number of elements in horizontal (x) direction.
    \param nely Number of elements in vertical (y) direction.
    \param h Edge size (width == height).
    */
    Regular(size_t nelx, size_t nely, double h = 1.0);

    /**
    Number of elements.

    \return unsigned int.
    */
    size_t nelem() const;

    /**
    Number of nodes.

    \return unsigned int.
    */
    size_t nnode() const;

    /**
    Number of nodes-per-element.

    \return unsigned int.
    */
    size_t nne() const;

    /**
    Number of dimensions.

    \return unsigned int.
    */
    size_t ndim() const;

    /**
    Number of elements in x-direction == width of the mesh in units of h().

    \return unsigned int.
    */
    size_t nelx() const;

    /**
    Number of elements in y-direction == height of the mesh, in units of h(),

    \return unsigned int.
    */
    size_t nely() const;

    /**
    Edge size of one element.

    \return double.
    */
    double h() const;

    /**
    Element type.

    \return GooseFEM::Mesh::ElementType().
    */
    ElementType getElementType() const;

    /**
    Nodal coordinates.

    \return ``[nnode, ndim]``.
    */
    xt::xtensor<double, 2> coor() const;

    /**
    Connectivity.

    \return ``[nelem, nne]``.
    */
    xt::xtensor<size_t, 2> conn() const;

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

    /**
    Reconstruct class for given coordinates / connectivity.

    \param coor Nodal coordinates ``[nnode, ndim]`` with ``ndim == 2``.
    \param conn Connectivity ``[nne, nne]`` with ``nne == 4``.
    \throw GOOSEFEM_CHECK()
    */
    FineLayer(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn);

    /**
    Number of elements.

    \return unsigned int.
    */
    size_t nelem() const;

    /**
    Number of nodes.

    \return unsigned int.
    */
    size_t nnode() const;

    /**
    Number of nodes-per-element.

    \return unsigned int.
    */
    size_t nne() const;

    /**
    Number of dimensions.

    \return unsigned int.
    */
    size_t ndim() const;

    /**
    Number of elements in x-direction along the middle layer == width of the mesh in units of h().

    \return unsigned int.
    */
    size_t nelx() const;

    /**
    Height of the mesh, in units of h()

    \return unsigned int.
    */
    size_t nely() const;

    /**
    Edge size of the smallest elements (along the middle layer).

    \return double.
    */
    double h() const;

    // edge size, per row of elements (in units of "h")
    xt::xtensor<size_t, 1> elemrow_nhx() const;
    xt::xtensor<size_t, 1> elemrow_nhy() const;
    xt::xtensor<size_t, 1> elemrow_nelem() const;

    /**
    Element type.

    \return GooseFEM::Mesh::ElementType().
    */
    ElementType getElementType() const;

    /**
    Nodal coordinates.

    \return ``[nnode, ndim]``.
    */
    xt::xtensor<double, 2> coor() const;

    /**
    Connectivity.

    \return ``[nelem, nne]``.
    */
    xt::xtensor<size_t, 2> conn() const;

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

    [[ deprecated ]]
    GooseFEM::Mesh::Quad4::FineLayer FineLayer(
        const xt::xtensor<double, 2>& coor,
        const xt::xtensor<size_t, 2>& conn);

    /**
    Refine a Regular mesh: subdivide elements in several smaller elements.
    */
    class RefineRegular {
    public:

        RefineRegular() = default;

        /**
        Constructor.

        \param mesh the coarse mesh.
        \param nx for each coarse element: number of fine elements in x-direction.
        \param ny for each coarse element: number of fine elements in y-direction.
        */
        RefineRegular(const GooseFEM::Mesh::Quad4::Regular& mesh, size_t nx, size_t ny);

        /**
        For each coarse element: number of fine elements in x-direction.

        \return unsigned int (same as used in constructor)
        */
        size_t nx() const;

        /**
        For each coarse element: number of fine elements in y-direction.

        \return unsigned int (same as used in constructor)
        */
        size_t ny() const;

        /**
        Obtain the coarse mesh (copy of the mesh passed to the constructor).

        \return mesh
        */
        GooseFEM::Mesh::Quad4::Regular getCoarseMesh() const;

        /**
        Obtain the fine mesh.

        \return mesh
        */
        GooseFEM::Mesh::Quad4::Regular getFineMesh() const;

        /**
        Get element-mapping: elements of the fine mesh per element of the coarse mesh.

        \return [nelem_coarse, nx() * ny()]
        */
        xt::xtensor<size_t, 2> getMap() const;

        // map field
        [[ deprecated ]]
        xt::xtensor<double, 2> mapToCoarse(const xt::xtensor<double, 1>& data) const; // scalar per el

        [[ deprecated ]]
        xt::xtensor<double, 2> mapToCoarse(const xt::xtensor<double, 2>& data) const; // scalar per intpnt

        [[ deprecated ]]
        xt::xtensor<double, 4> mapToCoarse(const xt::xtensor<double, 4>& data) const; // tensor per intpnt

        /**
        Compute the mean of the quantity define on the fine mesh when mapped on the coarse mesh.

        \tparam T type of the data (e.g. ``double``).
        \tparam rank rank of the data.
        \param data the data [nelem_fine, ...]
        \return the average data of the coarse mesh [nelem_coarse, ...]
        */
        template <class T, size_t rank>
        xt::xtensor<T, rank> meanToCoarse(const xt::xtensor<T, rank>& data) const;

        /**
        Compute the average of the quantity define on the fine mesh when mapped on the coarse mesh.

        \tparam T type of the data (e.g. ``double``).
        \tparam rank rank of the data.
        \tparam S type of the weights (e.g. ``double``).
        \param data the data [nelem_fine, ...]
        \param weights the weights [nelem_fine, ...]
        \return the average data of the coarse mesh [nelem_coarse, ...]
        */
        template <class T, size_t rank, class S>
        xt::xtensor<T, rank> averageToCoarse(
            const xt::xtensor<T, rank>& data,
            const xt::xtensor<S, rank>& weights) const;

        /**
        Map element quantities to the fine mesh.
        The mapping is a bit simplistic: no interpolation is involved.
        The mapping is such that::

            ret[e_fine, ...] <- data[e_coarse, ...]

        \tparam T type of the data (e.g. ``double``).
        \tparam rank rank of the data.
        \param data the data.
        \return mapped data.
        */
        template <class T, size_t rank>
        xt::xtensor<T, rank> mapToFine(const xt::xtensor<T, rank>& data) const;

    private:

        GooseFEM::Mesh::Quad4::Regular m_coarse; ///< the coarse mesh
        GooseFEM::Mesh::Quad4::Regular m_fine; ///< the fine mesh
        size_t m_nx; ///< see nx()
        size_t m_ny; ///< see ny()
        xt::xtensor<size_t, 2> m_coarse2fine; ///< see getMap()
    };

    /**
    Map a FineLayer mesh to a Regular mesh.
    The element size of the Regular corresponds to the smallest elements of the FineLayer mesh
    (along the middle layer).
    */
    class FineLayer2Regular {
    public:

        FineLayer2Regular() = default;

        /**
        Constructors.

        \param mesh The FineLayer mesh.
        */
        FineLayer2Regular(const GooseFEM::Mesh::Quad4::FineLayer& mesh);

        /**
        Obtain the Regular mesh.

        \return mesh.
        */
        GooseFEM::Mesh::Quad4::Regular getRegularMesh() const;

        /**
        Obtain the FineLayer mesh (copy of the mesh passed to the constructor).

        \return mesh.
        */
        GooseFEM::Mesh::Quad4::FineLayer getFineLayerMesh() const;

        // elements of the Regular mesh per element of the FineLayer mesh
        // and the fraction by which the overlap is

        /**
        Get element-mapping: elements of the Regular mesh per element of the FineLayer mesh.
        The number of Regular elements varies between elements of the FineLayer mesh.

        \return [nelem_finelayer, ?]
        */
        std::vector<std::vector<size_t>> getMap() const;

        /**
        To overlap fraction for each item in the mapping in getMap().

        \return [nelem_finelayer, ?]
        */
        std::vector<std::vector<double>> getMapFraction() const;

        /**
        Map element quantities to Regular.
        The mapping is a bit simplistic: no interpolation is involved, the function just
        accounts the fraction of overlap between the FineLayer element and the Regular element.
        The mapping is such that::

            ret[e_regular, ...] <- arg[e_finelayer, ...]

        \tparam T type of the data (e.g. ``double``).
        \tparam rank rank of the data.
        \param arg data.
        \return mapped data.
        */
        template <class T, size_t rank>
        xt::xtensor<T, rank> mapToRegular(const xt::xtensor<T, rank>& arg) const;

    private:

        GooseFEM::Mesh::Quad4::FineLayer m_finelayer; ///< the FineLayer mesh to map
        GooseFEM::Mesh::Quad4::Regular m_regular; ///< the new Regular mesh to which to map
        std::vector<std::vector<size_t>> m_elem_regular; ///< see getMap()
        std::vector<std::vector<double>> m_frac_regular; ///< see getMapFraction()
    };

} // namespace Map

} // namespace Quad4
} // namespace Mesh
} // namespace GooseFEM

#include "MeshQuad4.hpp"

#endif

/**
Generate simple meshes of 4-noded quadrilateral elements in 2d (GooseFEM::Mesh::ElementType::Quad4).

\file MeshQuad4.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESHQUAD4_H
#define GOOSEFEM_MESHQUAD4_H

#include "Mesh.h"
#include "config.h"

namespace GooseFEM {
namespace Mesh {

/**
Simple meshes of 3-noded quadrilateral elements in 2d (ElementType::Quad4).
*/
namespace Quad4 {

// pre-allocation
namespace Map {
class FineLayer2Regular;
}

/**
Regular mesh: equi-sized elements.
*/
class Regular : public RegularBase2d<Regular> {
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
    Element numbers as 'matrix'.

    \return [#nely, #nelx].
    */
    xt::xtensor<size_t, 2> elementgrid() const;

private:
    friend class RegularBase<Regular>;
    friend class RegularBase2d<Regular>;

    size_t nelx_impl() const;
    size_t nely_impl() const;
    ElementType getElementType_impl() const;
    xt::xtensor<double, 2> coor_impl() const;
    xt::xtensor<size_t, 2> conn_impl() const;
    xt::xtensor<size_t, 1> nodesBottomEdge_impl() const;
    xt::xtensor<size_t, 1> nodesTopEdge_impl() const;
    xt::xtensor<size_t, 1> nodesLeftEdge_impl() const;
    xt::xtensor<size_t, 1> nodesRightEdge_impl() const;
    xt::xtensor<size_t, 1> nodesBottomOpenEdge_impl() const;
    xt::xtensor<size_t, 1> nodesTopOpenEdge_impl() const;
    xt::xtensor<size_t, 1> nodesLeftOpenEdge_impl() const;
    xt::xtensor<size_t, 1> nodesRightOpenEdge_impl() const;
    size_t nodesBottomLeftCorner_impl() const;
    size_t nodesBottomRightCorner_impl() const;
    size_t nodesTopLeftCorner_impl() const;
    size_t nodesTopRightCorner_impl() const;

    double m_h; ///< See h()
    size_t m_nelx; ///< See nelx()
    size_t m_nely; ///< See nely()
    size_t m_nelem; ///< See nelem()
    size_t m_nnode; ///< See nnode()
    size_t m_nne; ///< See nne()
    size_t m_ndim; ///< See ndim()
};

/**
Mesh with fine middle layer, and coarser elements towards the top and bottom.
*/
class FineLayer : public RegularBase2d<FineLayer> {
public:
    FineLayer() = default;

    /**
    Constructor.

    \param nelx Number of elements (along the middle layer) in horizontal (x) direction.
    \param nely Approximate equivalent number of elements in vertical (y) direction.
    \param h Edge size (width == height) of elements along the weak layer.

    \param nfine
        Extra number of fine layers around the middle layer.
        By default the element size is kept smaller than the distance to the middle layer.
    */
    FineLayer(size_t nelx, size_t nely, double h = 1.0, size_t nfine = 1);

    /**
    Reconstruct class for given coordinates / connectivity.

    \tparam C e.g. `xt::xtensor<double, 2>`
    \tparam E e.g. `xt::xtensor<size_t, 2>`
    \param coor Nodal coordinates ``[nnode, ndim]`` with ``ndim == 2``.
    \param conn Connectivity ``[nne, nne]`` with ``nne == 4``.
    \throw GOOSEFEM_CHECK()
    */
    template <class C, class E, std::enable_if_t<xt::is_xexpression<C>::value, bool> = true>
    FineLayer(const C& coor, const E& conn);

    /**
    Edge size in x-direction of a block, in units of #h, per row of blocks.
    Note that a block is equal to an element except in refinement layers
    where it contains three elements.

    \return List of size equal to the number of rows of blocks.
    */
    xt::xtensor<size_t, 1> elemrow_nhx() const;

    /**
    Edge size in y-direction of a block, in units of #h, per row of blocks.
    Note that a block is equal to an element except in refinement layers
    where it contains three elements.

    \return List of size equal to the number of rows of blocks.
    */
    xt::xtensor<size_t, 1> elemrow_nhy() const;

    /**
    Per row of blocks:
    *   `-1`: normal layer
    *   `0`: transition layer to match coarse and finer element on the previous/next row.

    \return List of size equal to the number of rows of blocks.
    */
    xt::xtensor<int, 1> elemrow_type() const;

    /**
    Number of elements per row of blocks.
    Note that a block is equal to an element except in refinement layers
    where it contains three elements.

    \return List of size equal to the number of rows of blocks.
    */
    xt::xtensor<size_t, 1> elemrow_nelem() const;

    /**
    Elements in the middle (fine) layer.

    \return List of element numbers.
    */
    xt::xtensor<size_t, 1> elementsMiddleLayer() const;

    /**
    Elements along a layer.

    \return List of element numbers.
    */
    xt::xtensor<size_t, 1> elementsLayer(size_t layer) const;

    /**
    Select region of elements from 'matrix' of element numbers.

    \return List of element numbers.
    */
    xt::xtensor<size_t, 1> elementgrid_ravel(
        std::vector<size_t> rows_start_stop,
        std::vector<size_t> cols_start_stop) const;

    /**
    Select region of elements from 'matrix' of element numbers around an element:
    square box with edge-size ``(2 * size + 1) * h``, around ``element``.

    \param element The element around which to select elements.
    \param size Edge size of the square box encapsulating the selected element.
    \param periodic Assume the mesh periodic.
    \return List of elements.
    */
    xt::xtensor<size_t, 1>
    elementgrid_around_ravel(size_t element, size_t size, bool periodic = true);

    /**
    Select region of elements from 'matrix' of element numbers around an element:
    left/right from ``element`` (on the same layer).

    \param element The element around which to select elements.
    \param left Number of elements to select to the left.
    \param right Number of elements to select to the right.
    \param periodic Assume the mesh periodic.
    \return List of elements.
    */
    // -
    xt::xtensor<size_t, 1>
    elementgrid_leftright(size_t element, size_t left, size_t right, bool periodic = true);

    /**
    Mapping to 'roll' periodically in the x-direction,

    \return element mapping, such that: new_elemvar = elemvar[elem_map]
    */
    xt::xtensor<size_t, 1> roll(size_t n);

private:
    friend class RegularBase<FineLayer>;
    friend class RegularBase2d<FineLayer>;
    friend class GooseFEM::Mesh::Quad4::Map::FineLayer2Regular;

    size_t nelx_impl() const;
    size_t nely_impl() const;
    ElementType getElementType_impl() const;
    xt::xtensor<double, 2> coor_impl() const;
    xt::xtensor<size_t, 2> conn_impl() const;
    xt::xtensor<size_t, 1> nodesBottomEdge_impl() const;
    xt::xtensor<size_t, 1> nodesTopEdge_impl() const;
    xt::xtensor<size_t, 1> nodesLeftEdge_impl() const;
    xt::xtensor<size_t, 1> nodesRightEdge_impl() const;
    xt::xtensor<size_t, 1> nodesBottomOpenEdge_impl() const;
    xt::xtensor<size_t, 1> nodesTopOpenEdge_impl() const;
    xt::xtensor<size_t, 1> nodesLeftOpenEdge_impl() const;
    xt::xtensor<size_t, 1> nodesRightOpenEdge_impl() const;
    size_t nodesBottomLeftCorner_impl() const;
    size_t nodesBottomRightCorner_impl() const;
    size_t nodesTopLeftCorner_impl() const;
    size_t nodesTopRightCorner_impl() const;

    double m_h; ///< See h()
    size_t m_nelem; ///< See nelem()
    size_t m_nnode; ///< See nnode()
    size_t m_nne; ///< See nne()
    size_t m_ndim; ///< See ndim()
    double m_Lx; ///< Mesh size in x-direction.
    xt::xtensor<size_t, 1> m_layer_nelx; ///< See elemrow_nelem().
    xt::xtensor<size_t, 1> m_nhx; ///< See elemrow_nhx().
    xt::xtensor<size_t, 1> m_nhy; ///< See elemrow_nhy().
    xt::xtensor<size_t, 1> m_nnd; ///< total #nodes in the main node layer per node layer in "y"
    xt::xtensor<int, 1> m_refine; ///< See elemrow_type().
    xt::xtensor<size_t, 1> m_startElem; ///< start element per element layer in "y"
    xt::xtensor<size_t, 1> m_startNode; ///< start node per node layer in "y"

    /**
    \copydoc FineLayer::FineLayer(size_t, size_t, double, size_t)
    */
    void init(size_t nelx, size_t nely, double h, size_t nfine = 1);

    /**
    \copydoc FineLayer::FineLayer(const C&, const E&)
    */
    template <class C, class E>
    void map(const C& coor, const E& conn);
};

/**
Mesh mappings.
*/
namespace Map {

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
    xt::xtensor<T, rank>
    averageToCoarse(const xt::xtensor<T, rank>& data, const xt::xtensor<S, rank>& weights) const;

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

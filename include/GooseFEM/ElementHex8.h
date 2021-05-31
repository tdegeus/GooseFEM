/**
Quadrature for 8-noded hexahedral element in 3d (GooseFEM::Mesh::ElementType::Hex8),
in a Cartesian coordinate system.

\file ElementHex8.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENTHEX8_H
#define GOOSEFEM_ELEMENTHEX8_H

#include "config.h"

namespace GooseFEM {
namespace Element {

/**
8-noded hexahedral element in 3d (GooseFEM::Mesh::ElementType::Hex8).
*/
namespace Hex8 {

/**
gauss quadrature: quadrature points such that integration is exact for these bi-linear elements::
*/
namespace Gauss {

    /**
    Number of integration points::

        nip = nne = 8

    \return unsigned int
    */
    inline size_t nip();

    /**
    Integration point coordinates (local coordinates).

    \return Coordinates [#nip, ndim], with `ndim = 3`.
    */
    inline xt::xtensor<double, 2> xi();

    /**
    Integration point weights.

    \return Coordinates [#nip].
    */
    inline xt::xtensor<double, 1> w();

} // namespace Gauss

/**
nodal quadrature: quadrature points coincide with the nodes.
The order is the same as in the connectivity.
*/
namespace Nodal {

    /**
    Number of integration points::

        nip = nne = 8

    \return unsigned int
    */
    inline size_t nip();

    /**
    Integration point coordinates (local coordinates).

    \return Coordinates [#nip, `ndim`], with ``ndim = 3``.
    */
    inline xt::xtensor<double, 2> xi();

    /**
    Integration point weights.

    \return Coordinates [#nip].
    */
    inline xt::xtensor<double, 1> w();

} // namespace Nodal

/**
Interpolation and quadrature.

Fixed dimensions:
-   ``ndim = 3``: number of dimensions.
-   ``nne = 8``: number of nodes per element.

Naming convention:
-    ``elemmat``:  matrices stored per element, [#nelem, #nne * #ndim, #nne * #ndim]
-    ``elemvec``:  nodal vectors stored per element, [#nelem, #nne, #ndim]
-    ``qtensor``:  integration point tensor, [#nelem, #nip, #ndim, #ndim]
-    ``qscalar``:  integration point scalar, [#nelem, #nip]
*/
class Quadrature : public QuadratureBaseCartesian<Quadrature> {
public:

    Quadrature() = default;

    /**
    Constructor: use default Gauss integration.
    The following is pre-computed during construction:
    -   the shape functions,
    -   the shape function gradients (in local and global) coordinates,
    -   the integration points volumes.
    They can be reused without any cost.
    They only have to be recomputed when the nodal position changes
    (note that they are assumed to be constant under a small-strain assumption).
    In that case use update_x() to update the nodal positions and
    to recompute the above listed quantities.

    \param x nodal coordinates (``elemvec``).
    */
    template <class T>
    Quadrature(const T& x);

    /**
    Constructor with custom integration.
    The following is pre-computed during construction:
    -   the shape functions,
    -   the shape function gradients (in local and global) coordinates,
    -   the integration points volumes.
    They can be reused without any cost.
    They only have to be recomputed when the nodal position changes
    (note that they are assumed to be constant under a small-strain assumption).
    In that case use update_x() to update the nodal positions and
    to recompute the above listed quantities.

    \param x nodal coordinates (``elemvec``).
    \param xi Integration point coordinates (local coordinates) [#nip].
    \param w Integration point weights [#nip].
    */
    template <class T, class X, class W>
    Quadrature(const T& x, const X& xi, const W& w);

private:

    friend QuadratureBase<Quadrature>;
    friend QuadratureBaseCartesian<Quadrature>;

    template <class T, class R>
    void int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const;

    template <class T, class R>
    void int_gradN_dot_tensor2_dV_impl(const T& qtensor, R& elemvec) const;

    constexpr static size_t s_nne = 8; ///< Number of nodes per element.
    constexpr static size_t s_ndim = 3; ///< Number of dimensions for nodal vectors.
    constexpr static size_t s_tdim = 3; ///< Number of dimensions for tensors.
    size_t m_tdim = 3; ///< Dynamic alias of s_tdim (remove in C++17)
    size_t m_nelem; ///< Number of elements.
    size_t m_nip;   ///< Number of integration points per element.
    xt::xtensor<double, 3> m_x;    ///< nodal positions stored per element [#nelem, #nne, #ndim]
    xt::xtensor<double, 1> m_w;    ///< weight of each integration point [nip]
    xt::xtensor<double, 2> m_xi;   ///< local coordinate of each integration point [#nip, #ndim]
    xt::xtensor<double, 2> m_N;    ///< shape functions [#nip, #nne]
    xt::xtensor<double, 3> m_dNxi; ///< shape function grad. wrt local  coor. [#nip, #nne, #ndim]
    xt::xtensor<double, 4> m_dNx;  ///< shape function grad. wrt global coor. [#nelem, #nip, #nne, #ndim]
    xt::xtensor<double, 2> m_vol;  ///< integration point volume [#nelem, #nip]
};

} // namespace Hex8
} // namespace Element
} // namespace GooseFEM

#include "ElementHex8.hpp"

#endif

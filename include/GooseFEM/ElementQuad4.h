/**
Quadrature for 4-noded quadrilateral element in 2d (GooseFEM::Mesh::ElementType::Quad4),
in a Cartesian coordinate system.

\file ElementQuad4.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENTQUAD4_H
#define GOOSEFEM_ELEMENTQUAD4_H

#include "config.h"
#include "Element.h"

namespace GooseFEM {
namespace Element {

/**
4-noded quadrilateral element in 2d (GooseFEM::Mesh::ElementType::Quad4).
*/
namespace Quad4 {

/**
Gauss quadrature: quadrature points such that integration is exact for this bi-linear element::

    + ----------- +
    |             |
    |   3     2   |
    |             |
    |   0     1   |
    |             |
    + ----------- +
*/
namespace Gauss {

    /**
    Number of integration points:

        nip = nne = 4

    \return unsigned int
    */
    inline size_t nip();

    /**
    Integration point coordinates (local coordinates).

    \return Coordinates [#nip, `ndim`], with `ndim = 2`.
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
The order is the same as in the connectivity:

    3 -- 2
    |    |
    0 -- 1
*/
namespace Nodal {

    /**
    Number of integration points::

        nip = nne = 4

    \return unsigned int
    */
    inline size_t nip();

    /**
    Integration point coordinates (local coordinates).

    \return Coordinates [#nip, `ndim`], with ``ndim = 2``.
    */
    inline xt::xtensor<double, 2> xi();

    /**
    Integration point weights.

    \return Coordinates [#nip].
    */
    inline xt::xtensor<double, 1> w();

} // namespace Nodal

/**
midpoint quadrature: quadrature points in the middle of the element::

    + ------- +
    |         |
    |    0    |
    |         |
    + ------- +
*/
namespace MidPoint {

    /**
    Number of integration points::

        nip = 1

    \return unsigned int
    */
    inline size_t nip();

    /**
    Integration point coordinates (local coordinates).

    \return Coordinates [#nip, ``ndim``], with ``ndim = 2``.
    */
    inline xt::xtensor<double, 2> xi();

    /**
    Integration point weights.

    \return Coordinates [#nip].
    */
    inline xt::xtensor<double, 1> w();

} // namespace MidPoint

/**
Interpolation and quadrature.

Fixed dimensions:
-   ``ndim = 2``: number of dimensions.
-   ``nne = 4``: number of nodes per element.

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
    void interpQuad_vector_impl(const T&, R& qvector) const;

    template <class T, class R>
    void gradN_vector_impl(const T&, R& qtensor) const;

    template <class T, class R>
    void gradN_vector_T_impl(const T&, R& qtensor) const;

    template <class T, class R>
    void symGradN_vector_impl(const T&, R& qtensor) const;

    template <class T, class R>
    void int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const;

    template <class T, class R>
    void int_gradN_dot_tensor2_dV_impl(const T& qtensor, R&) const;

    void compute_dN_impl();

    constexpr static size_t s_nne = 4;  ///< Number of nodes per element.
    constexpr static size_t s_ndim = 2; ///< Number of dimensions for nodal vectors.
    constexpr static size_t s_tdim = 2; ///< Number of dimensions for tensors.
    size_t m_tdim = 2; ///< Dynamic alias of s_tdim (remove in C++17)
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

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#include "ElementQuad4.hpp"

#endif

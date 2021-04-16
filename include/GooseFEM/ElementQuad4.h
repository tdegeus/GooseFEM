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
    Number of integration points::

        nip = nne = 4

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

} // namespace Gauss

/**
nodal quadrature: quadrature points coincide with the nodes.
The order is the same as in the connectivity::

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

    \return Coordinates [#nip, ``ndim``], with ``ndim = 2``.
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
class Quadrature : public GooseFEM::Element::QuadratureBaseCartesian<4, 2, 2> {
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
    Quadrature(const xt::xtensor<double, 3>& x);

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
    Quadrature(
        const xt::xtensor<double, 3>& x,
        const xt::xtensor<double, 2>& xi,
        const xt::xtensor<double, 1>& w);

    void interpQuad_vector(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 3>& qvector) const override;

    void gradN_vector(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const override;

    void gradN_vector_T(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const override;

    void symGradN_vector(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const override;

    void int_N_scalar_NT_dV(
        const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 3>& elemmat) const override;

    void int_gradN_dot_tensor2_dV(
        const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const override;

protected:
    void compute_dN() override;
};

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#include "ElementQuad4.hpp"

#endif

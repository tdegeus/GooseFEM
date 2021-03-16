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

    \return Coordinates [nip(), ``ndim``], with ``ndim = 3``.
    */
    inline xt::xtensor<double, 2> xi();

    /**
    Integration point weights.

    \return Coordinates [nip()].
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

    \return Coordinates [nip(), ``ndim``], with ``ndim = 3``.
    */
    inline xt::xtensor<double, 2> xi();

    /**
    Integration point weights.

    \return Coordinates [nip()].
    */
    inline xt::xtensor<double, 1> w();

} // namespace Nodal


/**
Interpolation and quadrature.

Fixed dimensions:
-   ``ndim = 3``: number of dimensions.
-   ``nne = 8``: number of nodes per element.

Naming convention:
-    ``elemmat``:  matrices stored per element, [nelem(), nne() * ndim(), nne() * ndim()]
-    ``elemvec``:  nodal vectors stored per element, [nelem(), nne(), ndim()]
-    ``qtensor``:  integration point tensor, [nelem(), nip(), ndim(), ndim()]
-    ``qscalar``:  integration point scalar, [nelem(), nip()]
*/
class Quadrature : public GooseFEM::Element::QuadratureBaseCartesian<8, 3, 3> {
public:

    Quadrature() = default;

    /**
    Constructor: use default Gauss integration.
    During construction the values of the shape functions and the shape function gradients
    (in local and global) coordinates are computed. They can be reused without any cost.
    They only have to be recomputed when the nodal position changes, taken care of in update_x().

    \param x nodal coordinates (``elemvec``).
    */
    Quadrature(const xt::xtensor<double, 3>& x);

    /**
    Constructor with custom integration.
    During construction the values of the shape functions and the shape function gradients
    (in local and global) coordinates are computed. They can be reused without any cost.
    They only have to be recomputed when the nodal position changes, taken care of in update_x().

    \param x nodal coordinates (``elemvec``).
    \param xi Integration point coordinates (local coordinates) [nip()].
    \param w Integration point weights [nip()].
    */
    Quadrature(
        const xt::xtensor<double, 3>& x,
        const xt::xtensor<double, 2>& xi,
        const xt::xtensor<double, 1>& w);

    void int_N_scalar_NT_dV(
        const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 3>& elemmat) const override;

    void int_gradN_dot_tensor2_dV(
        const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const override;
};

} // namespace Hex8
} // namespace Element
} // namespace GooseFEM

#include "ElementHex8.hpp"

#endif

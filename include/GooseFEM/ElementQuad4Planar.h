/**
Quadrature for 4-noded quadrilateral element in 2d (GooseFEM::Mesh::ElementType::Quad4),
in a Cartesian coordinate system.
The different with ElementQuad4.h is that here the tensors live in 3d and are assumed plane strain.

\file ElementQuad4Planar.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENTQUAD4PLANAR_H
#define GOOSEFEM_ELEMENTQUAD4PLANAR_H

#include "config.h"

namespace GooseFEM {
namespace Element {
namespace Quad4 {

/**
Interpolation and quadrature.
Similar to Element::Quad4::Quadrature with the only different that quadrature point tensors
are 3d ("plane strain") while the mesh is 2d.

Fixed dimensions:
-   ``ndim = 2``: number of dimensions.
-   ``tdim = 3``: number of dimensions or tensor.
-   ``nne = 4``: number of nodes per element.

Naming convention:
-    ``elemmat``:  matrices stored per element, [#nelem, #nne * #ndim, #nne * #ndim]
-    ``elemvec``:  nodal vectors stored per element, [#nelem, #nne, #ndim]
-    ``qtensor``:  integration point tensor, [#nelem, #nip, #tdim, #tdim]
-    ``qscalar``:  integration point scalar, [#nelem, #nip]
*/
class QuadraturePlanar : public GooseFEM::Element::QuadratureBaseCartesian<4, 2, 3> {
public:

    QuadraturePlanar() = default;

    /**
    Constructor: use default Gauss integration.
    During construction the values of the shape functions and the shape function gradients
    (in local and global) coordinates are computed. They can be reused without any cost.
    They only have to be recomputed when the nodal position changes.
    In that case use update_x() to update the nodal positions and to recompute the
    shape functions and their gradients.

    \param x nodal coordinates (``elemvec``).
    \param thick out-of-plane thickness (incorporated in the element volume).
    */
    QuadraturePlanar(const xt::xtensor<double, 3>& x, double thick = 1.0);

    /**
    Constructor with custom integration.
    During construction the values of the shape functions and the shape function gradients
    (in local and global) coordinates are computed. They can be reused without any cost.
    They only have to be recomputed when the nodal position changes.
    In that case use update_x() to update the nodal positions and to recompute the
    shape functions and their gradients.

    \param x nodal coordinates (``elemvec``).
    \param xi Integration point coordinates (local coordinates) [#nip].
    \param w Integration point weights [#nip].
    \param thick out-of-plane thickness (incorporated in the element volume).
    */
    QuadraturePlanar(
        const xt::xtensor<double, 3>& x,
        const xt::xtensor<double, 2>& xi,
        const xt::xtensor<double, 1>& w,
        double thick = 1.0);

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

private:
    double m_thick; ///< out-of-plane thickness
};

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#include "ElementQuad4Planar.hpp"

#endif

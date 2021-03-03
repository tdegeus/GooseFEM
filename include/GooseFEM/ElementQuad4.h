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
namespace Quad4 {

/**
gauss quadrature: quadrature points such that integration is exact for these bi-linear elements::

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

    \return Coordinates [nip(), ``ndim``], with ``ndim = 2``.
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

    \return Coordinates [nip(), ``ndim``], with ``ndim = 2``.
    */
    inline xt::xtensor<double, 2> xi();

    /**
    Integration point weights.

    \return Coordinates [nip()].
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

    \return Coordinates [nip(), ``ndim``], with ``ndim = 2``.
    */
    inline xt::xtensor<double, 2> xi();

    /**
    Integration point weights.

    \return Coordinates [nip()].
    */
    inline xt::xtensor<double, 1> w();

} // namespace MidPoint

/**
Interpolation and quadrature.

Fixed dimensions:
-   ``ndim = 2``: number of dimensions.
-   ``nne  = 4``: number of nodes per element.

Naming convention:
-    ``elemmat``:  matrices stored per element, [nelem(), nne() * ndim(), nne() * ndim()]
-    ``elemvec``:  nodal vectors stored per element, [nelem(), nne(), ndim()]
-    ``qtensor``:  integration point tensor, [nelem(), nip(), ndim(), ndim()]
-    ``qscalar``:  integration point scalar, [nelem(), nip()]
*/
class Quadrature : public GooseFEM::Element::QuadratureBase<4, 2, 2> {
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

    /**
    Update the nodal positions.
    This recomputes the values of the shape functions and the shape function gradients
    (in local and global) coordinates.
    Under the small deformations assumption this should not be called.

    \param x nodal coordinates (``elemvec``). Shape should match the earlier definition.
    */
    void update_x(const xt::xtensor<double, 3>& x);

    /**
    Shape function gradients (in global coordinates).

    \return ``dN / dx`` stored per element, per integration point [nelem(), nip(), nne(), ndim()].
    */
    xt::xtensor<double, 4> GradN() const;

    /**
    Integration volume.

    \return volume stored per element, per integration point [nelem(), nip()].
    */
    xt::xtensor<double, 2> dV() const;

    /**
    Interpolate element vector.

    \param elemvec nodal vector stored per element (shape: [nelem(), nne(), ndim()]).
    \param qvector integration point vector, overwritten (shape: [nelem(), nip(), ndim()]).
    */
    template <class T>
    void interp_N_vector(const xt::xtensor<T, 3>& elemvec, xt::xtensor<T, 3>& qvector) const;

    /**
    Same as interp_N_vector(), but returns auto-allocated data.

    \param elemvec nodal vector stored per element (shape: [nelem(), nne(), ndim()]).
    \return integration point vector (shape: [nelem(), nip(), ndim()]).
    */
    template <class T>
    xt::xtensor<T, 3> Interp_N_vector(const xt::xtensor<T, 3>& elemvec) const;

    /**
    Element-by-element: dyadic product of the shape function gradients (GradN()) and a nodal vector::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    qtensor(e, q, i, j) += dNdx(e, q, m, i) * elemvec(e, m, j)

    \param elemvec [nelem(), nne(), ndim()]
    \param qtensor overwritten [nelem(), nip(), ndim(), ndim()]
    */
    void gradN_vector(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;

    /**
    Element-by-element: dyadic product of the shape function gradients (GradN()) and a nodal vector.
    The transpose of the resulting tensor is stored::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    qtensor(e, q, j, i) += dNdx(e, q, m, i) * elemvec(e, m, j)

    \param elemvec [nelem(), nne(), ndim()]
    \param qtensor overwritten [nelem(), nip(), ndim(), ndim()]
    */
    void gradN_vector_T(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;

    /**
    Element-by-element: dyadic product of the shape function gradients (GradN()) and a nodal vector.
    The symmetric part of the resulting tensor is stored::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    qtensor(e, q, i, j) += 0.5 * dNdx(e, q, m, i) * elemvec(e, m, j)
                    qtensor(e, q, j, i) += 0.5 * dNdx(e, q, m, i) * elemvec(e, m, j)

    \param elemvec [nelem(), nne(), ndim()]
    \param qtensor overwritten [nelem(), nip(), ndim(), ndim()]
    */
    void symGradN_vector(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;

    /**
    Element-by-element: integral of the scalar product.
    Within one one element::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    for n in range(nne):
                        elemmat(e, m * ndim + i, n * ndim + i) += N(e, q, m) * qscalar(e, q) * N(e, q, n) * dV(e, q)

    \param qscalar [nelem(), nip()]
    \param elemmat overwritten [nelem(), nne() * ndim(), nne() * ndim()]
    */
    void int_N_scalar_NT_dV(
        const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 3>& elemmat) const;

    // Integral of the dot product
    // elemvec(m,j) += dNdx(m,i) * qtensor(i,j) * dV
    void int_gradN_dot_tensor2_dV(
        const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const;

    // Integral of the dot product
    // elemmat(m*2+j, n*2+k) += dNdx(m,i) * qtensor(i,j,k,l) * dNdx(n,l) * dV
    void int_gradN_dot_tensor4_dot_gradNT_dV(
        const xt::xtensor<double, 6>& qtensor, xt::xtensor<double, 3>& elemmat) const;

    // Auto-allocation of the functions above
    xt::xtensor<double, 4> GradN_vector(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 4> GradN_vector_T(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 4> SymGradN_vector(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 3> Int_N_scalar_NT_dV(const xt::xtensor<double, 2>& qscalar) const;
    xt::xtensor<double, 3> Int_gradN_dot_tensor2_dV(const xt::xtensor<double, 4>& qtensor) const;
    xt::xtensor<double, 3> Int_gradN_dot_tensor4_dot_gradNT_dV(const xt::xtensor<double, 6>& qtensor) const;


private:
    /**
    Compute m_vol and m_dNdx based on current m_x.
    */
    void compute_dN();

private:
    xt::xtensor<double, 3> m_x;    ///< nodal positions stored per element [nelem(), nne(), ndim()]
    xt::xtensor<double, 1> m_w;    ///< weight of each integration point [nip]
    xt::xtensor<double, 2> m_xi;   ///< local coordinate of each integration point [nip(), ndim()]
    xt::xtensor<double, 2> m_N;    ///< shape functions [nip(), nne()]
    xt::xtensor<double, 3> m_dNxi; ///< shape function grad. wrt local  coor. [nip(), nne(), ndim()]
    xt::xtensor<double, 4> m_dNx;  ///< shape function grad. wrt global coor. [nelem(), nip(), nne(), ndim()]
    xt::xtensor<double, 2> m_vol;  ///< integration point volume [nelem(), nip()]
};

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#include "ElementQuad4.hpp"

#endif

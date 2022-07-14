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
#include "detail.h"

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
class QuadraturePlanar : public QuadratureBaseCartesian<QuadraturePlanar> {
public:
    QuadraturePlanar() = default;

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
    \param thick out-of-plane thickness (incorporated in the element volume).
    */
    template <class T>
    QuadraturePlanar(const T& x, double thick = 1.0)
        : QuadraturePlanar(x, Gauss::xi(), Gauss::w(), thick)
    {
    }

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
    \param thick out-of-plane thickness (incorporated in the element volume).
    */
    template <class T, class X, class W>
    QuadraturePlanar(const T& x, const X& xi, const W& w, double thick = 1.0)
    {
        m_x = x;
        m_w = w;
        m_xi = xi;
        m_nip = w.size();
        m_nelem = m_x.shape(0);
        m_N = xt::empty<double>({m_nip, s_nne});
        m_dNxi = xt::empty<double>({m_nip, s_nne, s_ndim});
        m_thick = thick;

        for (size_t q = 0; q < m_nip; ++q) {
            m_N(q, 0) = 0.25 * (1.0 - xi(q, 0)) * (1.0 - xi(q, 1));
            m_N(q, 1) = 0.25 * (1.0 + xi(q, 0)) * (1.0 - xi(q, 1));
            m_N(q, 2) = 0.25 * (1.0 + xi(q, 0)) * (1.0 + xi(q, 1));
            m_N(q, 3) = 0.25 * (1.0 - xi(q, 0)) * (1.0 + xi(q, 1));
        }

        for (size_t q = 0; q < m_nip; ++q) {
            // - dN / dxi_0
            m_dNxi(q, 0, 0) = -0.25 * (1.0 - xi(q, 1));
            m_dNxi(q, 1, 0) = +0.25 * (1.0 - xi(q, 1));
            m_dNxi(q, 2, 0) = +0.25 * (1.0 + xi(q, 1));
            m_dNxi(q, 3, 0) = -0.25 * (1.0 + xi(q, 1));
            // - dN / dxi_1
            m_dNxi(q, 0, 1) = -0.25 * (1.0 - xi(q, 0));
            m_dNxi(q, 1, 1) = -0.25 * (1.0 + xi(q, 0));
            m_dNxi(q, 2, 1) = +0.25 * (1.0 + xi(q, 0));
            m_dNxi(q, 3, 1) = +0.25 * (1.0 - xi(q, 0));
        }

        GOOSEFEM_ASSERT(m_x.shape(1) == s_nne);
        GOOSEFEM_ASSERT(m_x.shape(2) == s_ndim);
        GOOSEFEM_ASSERT(xt::has_shape(m_xi, {m_nip, s_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(m_w, {m_nip}));
        GOOSEFEM_ASSERT(xt::has_shape(m_N, {m_nip, s_nne}));
        GOOSEFEM_ASSERT(xt::has_shape(m_dNxi, {m_nip, s_nne, s_ndim}));

        m_dNx = xt::empty<double>({m_nelem, m_nip, s_nne, s_ndim});
        m_vol = xt::empty<double>(this->shape_qscalar());

        this->compute_dN_impl();
    }

private:
    friend QuadratureBase<QuadraturePlanar>;
    friend QuadratureBaseCartesian<QuadraturePlanar>;

    template <class T, class R>
    void gradN_vector_impl(const T& elemvec, R& qtensor) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));

        qtensor.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<s_nne, s_ndim>());
                auto gradu = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_tdim, s_tdim>());

                // gradu(i,j) += dNx(m,i) * u(m,j)
                gradu(0, 0) = dNx(0, 0) * u(0, 0) + dNx(1, 0) * u(1, 0) + dNx(2, 0) * u(2, 0) +
                              dNx(3, 0) * u(3, 0);
                gradu(0, 1) = dNx(0, 0) * u(0, 1) + dNx(1, 0) * u(1, 1) + dNx(2, 0) * u(2, 1) +
                              dNx(3, 0) * u(3, 1);
                gradu(1, 0) = dNx(0, 1) * u(0, 0) + dNx(1, 1) * u(1, 0) + dNx(2, 1) * u(2, 0) +
                              dNx(3, 1) * u(3, 0);
                gradu(1, 1) = dNx(0, 1) * u(0, 1) + dNx(1, 1) * u(1, 1) + dNx(2, 1) * u(2, 1) +
                              dNx(3, 1) * u(3, 1);
            }
        }
    }

    template <class T, class R>
    void gradN_vector_T_impl(const T& elemvec, R& qtensor) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));

        qtensor.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<s_nne, s_ndim>());
                auto gradu = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_tdim, s_tdim>());

                // gradu(j,i) += dNx(m,i) * u(m,j)
                gradu(0, 0) = dNx(0, 0) * u(0, 0) + dNx(1, 0) * u(1, 0) + dNx(2, 0) * u(2, 0) +
                              dNx(3, 0) * u(3, 0);
                gradu(1, 0) = dNx(0, 0) * u(0, 1) + dNx(1, 0) * u(1, 1) + dNx(2, 0) * u(2, 1) +
                              dNx(3, 0) * u(3, 1);
                gradu(0, 1) = dNx(0, 1) * u(0, 0) + dNx(1, 1) * u(1, 0) + dNx(2, 1) * u(2, 0) +
                              dNx(3, 1) * u(3, 0);
                gradu(1, 1) = dNx(0, 1) * u(0, 1) + dNx(1, 1) * u(1, 1) + dNx(2, 1) * u(2, 1) +
                              dNx(3, 1) * u(3, 1);
            }
        }
    }

    template <class T, class R>
    void symGradN_vector_impl(const T& elemvec, R& qtensor) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));

        qtensor.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto u = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<s_nne, s_ndim>());
                auto eps = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_tdim, s_tdim>());

                // gradu(i,j) += dNx(m,i) * u(m,j)
                // eps(j,i) = 0.5 * (gradu(i,j) + gradu(j,i))
                eps(0, 0) = dNx(0, 0) * u(0, 0) + dNx(1, 0) * u(1, 0) + dNx(2, 0) * u(2, 0) +
                            dNx(3, 0) * u(3, 0);
                eps(1, 1) = dNx(0, 1) * u(0, 1) + dNx(1, 1) * u(1, 1) + dNx(2, 1) * u(2, 1) +
                            dNx(3, 1) * u(3, 1);
                eps(0, 1) = 0.5 * (dNx(0, 0) * u(0, 1) + dNx(1, 0) * u(1, 1) + dNx(2, 0) * u(2, 1) +
                                   dNx(3, 0) * u(3, 1) + dNx(0, 1) * u(0, 0) + dNx(1, 1) * u(1, 0) +
                                   dNx(2, 1) * u(2, 0) + dNx(3, 1) * u(3, 0));
                eps(1, 0) = eps(0, 1);
            }
        }
    }

    template <class T, class R>
    void int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(qscalar, this->shape_qscalar()));
        GOOSEFEM_ASSERT(xt::has_shape(elemmat, this->shape_elemmat()));

        elemmat.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto M = xt::adapt(&elemmat(e, 0, 0), xt::xshape<s_nne * s_ndim, s_nne * s_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto N = xt::adapt(&m_N(q, 0), xt::xshape<s_nne>());
                auto& vol = m_vol(e, q);
                auto& rho = qscalar(e, q);

                // M(m*ndim+i,n*ndim+i) += N(m) * scalar * N(n) * dV
                for (size_t m = 0; m < s_nne; ++m) {
                    for (size_t n = 0; n < s_nne; ++n) {
                        M(m * s_ndim + 0, n * s_ndim + 0) += N(m) * rho * N(n) * vol;
                        M(m * s_ndim + 1, n * s_ndim + 1) += N(m) * rho * N(n) * vol;
                    }
                }
            }
        }
    }

    template <class T, class R>
    void int_gradN_dot_tensor2_dV_impl(const T& qtensor, R& elemvec) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->shape_qtensor<2>()));
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));

        elemvec.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {

            auto f = xt::adapt(&elemvec(e, 0, 0), xt::xshape<s_nne, s_ndim>());

            for (size_t q = 0; q < m_nip; ++q) {

                auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<s_nne, s_ndim>());
                auto sig = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<s_tdim, s_tdim>());
                auto& vol = m_vol(e, q);

                for (size_t m = 0; m < s_nne; ++m) {
                    f(m, 0) += (dNx(m, 0) * sig(0, 0) + dNx(m, 1) * sig(1, 0)) * vol;
                    f(m, 1) += (dNx(m, 0) * sig(0, 1) + dNx(m, 1) * sig(1, 1)) * vol;
                }
            }
        }
    }

    void compute_dN_impl()
    {
#pragma omp parallel
        {
            array_type::tensor<double, 2> J = xt::empty<double>({2, 2});
            array_type::tensor<double, 2> Jinv = xt::empty<double>({2, 2});

#pragma omp for
            for (size_t e = 0; e < m_nelem; ++e) {

                auto x = xt::adapt(&m_x(e, 0, 0), xt::xshape<s_nne, s_ndim>());

                for (size_t q = 0; q < m_nip; ++q) {

                    auto dNxi = xt::adapt(&m_dNxi(q, 0, 0), xt::xshape<s_nne, s_ndim>());
                    auto dNx = xt::adapt(&m_dNx(e, q, 0, 0), xt::xshape<s_nne, s_ndim>());

                    // J(i,j) += dNxi(m,i) * x(m,j);
                    J(0, 0) = dNxi(0, 0) * x(0, 0) + dNxi(1, 0) * x(1, 0) + dNxi(2, 0) * x(2, 0) +
                              dNxi(3, 0) * x(3, 0);
                    J(0, 1) = dNxi(0, 0) * x(0, 1) + dNxi(1, 0) * x(1, 1) + dNxi(2, 0) * x(2, 1) +
                              dNxi(3, 0) * x(3, 1);
                    J(1, 0) = dNxi(0, 1) * x(0, 0) + dNxi(1, 1) * x(1, 0) + dNxi(2, 1) * x(2, 0) +
                              dNxi(3, 1) * x(3, 0);
                    J(1, 1) = dNxi(0, 1) * x(0, 1) + dNxi(1, 1) * x(1, 1) + dNxi(2, 1) * x(2, 1) +
                              dNxi(3, 1) * x(3, 1);

                    double Jdet = detail::tensor<2>::inv(J, Jinv);

                    // dNx(m,i) += Jinv(i,j) * dNxi(m,j);
                    for (size_t m = 0; m < s_nne; ++m) {
                        dNx(m, 0) = Jinv(0, 0) * dNxi(m, 0) + Jinv(0, 1) * dNxi(m, 1);
                        dNx(m, 1) = Jinv(1, 0) * dNxi(m, 0) + Jinv(1, 1) * dNxi(m, 1);
                    }

                    m_vol(e, q) = m_w(q) * Jdet * m_thick;
                }
            }
        }
    }

    constexpr static size_t s_nne = 4; ///< Number of nodes per element.
    constexpr static size_t s_ndim = 2; ///< Number of dimensions for nodal vectors.
    constexpr static size_t s_tdim = 3; ///< Dynamic alias of s_tdim (remove in C++17)
    size_t m_tdim = 3; ///< Number of dimensions for tensors.
    size_t m_nelem; ///< Number of elements.
    size_t m_nip; ///< Number of integration points per element.
    array_type::tensor<double, 3> m_x; ///< nodal positions stored per element [#nelem, #nne, #ndim]
    array_type::tensor<double, 1> m_w; ///< weight of each integration point [nip]
    array_type::tensor<double, 2> m_xi; ///< local coordinate per integration point [#nip, #ndim]
    array_type::tensor<double, 2> m_N; ///< shape functions [#nip, #nne]
    array_type::tensor<double, 3> m_dNxi; ///< local shape func grad [#nip, #nne, #ndim]
    array_type::tensor<double, 4> m_dNx; ///< global shape func grad [#nelem, #nip, #nne, #ndim]
    array_type::tensor<double, 2> m_vol; ///< integration point volume [#nelem, #nip]
    double m_thick; ///< out-of-plane thickness
};

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#endif

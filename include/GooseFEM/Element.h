/**
Convenience methods for integration point data.

\file Element.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENT_H
#define GOOSEFEM_ELEMENT_H

#include "Allocate.h"
#include "config.h"
#include "detail.hpp"

namespace GooseFEM {

/**
Element quadrature and interpolation.
*/
namespace Element {

/**
Convert nodal vector with ("nodevec", shape:``[nnode, ndim]``) to nodal vector stored per element
("elemvec", shape: ``[nelem, nne, ndim]``).

\param conn Connectivity.
\param nodevec "nodevec".
\return "elemvec".
*/
inline array_type::tensor<double, 3> asElementVector(
    const array_type::tensor<size_t, 2>& conn,
    const array_type::tensor<double, 2>& nodevec)
{
    size_t nelem = conn.shape(0);
    size_t nne = conn.shape(1);
    size_t ndim = nodevec.shape(1);

    array_type::tensor<double, 3> elemvec = xt::empty<double>({nelem, nne, ndim});

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {
        for (size_t m = 0; m < nne; ++m) {
            for (size_t i = 0; i < ndim; ++i) {
                elemvec(e, m, i) = nodevec(conn(e, m), i);
            }
        }
    }

    return elemvec;
}

/**
Assemble nodal vector stored per element ("elemvec", shape ``[nelem, nne, ndim]``) to nodal vector
("nodevec", shape ``[nnode, ndim]``).

\param conn Connectivity.
\param elemvec "elemvec".
\return "nodevec".
*/
inline array_type::tensor<double, 2> assembleNodeVector(
    const array_type::tensor<size_t, 2>& conn,
    const array_type::tensor<double, 3>& elemvec)
{
    size_t nelem = conn.shape(0);
    size_t nne = conn.shape(1);
    size_t ndim = elemvec.shape(2);
    size_t nnode = xt::amax(conn)() + 1;

    GOOSEFEM_ASSERT(elemvec.shape(0) == nelem);
    GOOSEFEM_ASSERT(elemvec.shape(1) == nne);

    array_type::tensor<double, 2> nodevec = xt::zeros<double>({nnode, ndim});

    for (size_t e = 0; e < nelem; ++e) {
        for (size_t m = 0; m < nne; ++m) {
            for (size_t i = 0; i < ndim; ++i) {
                nodevec(conn(e, m), i) += elemvec(e, m, i);
            }
        }
    }

    return nodevec;
}

/**
Check that DOFs leave no holes.

\param dofs DOFs ("nodevec")
\return ``true`` if there are no holds.
*/
template <class E>
inline bool isSequential(const E& dofs)
{
    size_t ndof = xt::amax(dofs)() + 1;

    array_type::tensor<int, 1> exists = xt::zeros<int>({ndof});

    for (auto& i : dofs) {
        exists[i]++;
    }

    for (auto& i : dofs) {
        if (exists[i] == 0) {
            return false;
        }
    }

    return true;
}

/**
Check that all of the matrices stored per elemmat (shape: ``[nelem, nne * ndim, nne * ndim]``)
are diagonal.

\param elemmat Element-vectors ("elemmat")
\return ``true`` if all element matrices are diagonal.
*/
bool isDiagonal(const array_type::tensor<double, 3>& elemmat)
{
    GOOSEFEM_ASSERT(elemmat.shape(1) == elemmat.shape(2));

    size_t nelem = elemmat.shape(0);
    size_t N = elemmat.shape(1);

    double eps = std::numeric_limits<double>::epsilon();

#pragma omp parallel for
    for (size_t e = 0; e < nelem; ++e) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                if (i != j) {
                    if (std::abs(elemmat(e, i, j)) > eps) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

/**
CRTP base class for quadrature.
*/
template <class D>
class QuadratureBase {
public:
    /**
    Underlying type.
    */
    using derived_type = D;

    /**
    Number of elements.

    \return Scalar.
    */
    auto nelem() const
    {
        return derived_cast().m_nelem;
    }

    /**
    Number of nodes per element.

    \return Scalar.
    */
    auto nne() const
    {
        return D::s_nne;
    }

    /**
    Number of dimensions for node vectors.

    \return Scalar.
    */
    auto ndim() const
    {
        return D::s_ndim;
    }

    /**
    Number of dimensions for integration point tensors.

    \return Scalar.
    */
    auto tdim() const
    {
        return D::s_tdim;
    }

    /**
    Number of integration points.

    \return Scalar.
    */
    auto nip() const
    {
        return derived_cast().m_nip;
    }

    /**
    Convert "qscalar" to "qtensor" of certain rank.
    Fully allocated output passed as reference, use AsTensor to allocate and return data.

    \param arg A "qscalar".
    \param ret A "qtensor".
    */
    template <class T, class R>
    void asTensor(const T& arg, R& ret) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_qscalar()));
        GooseFEM::asTensor(arg, ret);
    }

    /**
    Convert "qscalar" to "qtensor" of certain rank.

    \param arg A "qscalar".
    \return "qtensor".
    */
    template <size_t rank, class T>
    auto AsTensor(const T& arg) const
    {
        return GooseFEM::AsTensor<rank>(arg, derived_cast().m_tdim);
    }

    /**
    Convert "qscalar" to "qtensor" of certain rank.

    \param rank Tensor rank.
    \param arg A "qscalar".
    \return "qtensor".
    */
    template <class T>
    auto AsTensor(size_t rank, const T& arg) const
    {
        return GooseFEM::AsTensor(rank, arg, derived_cast().m_tdim);
    }

    /**
    Get the shape of an "elemvec".

    \returns [#nelem, #nne, #ndim].
    */
    auto shape_elemvec() const -> std::array<size_t, 3>
    {
        return std::array<size_t, 3>{derived_cast().m_nelem, D::s_nne, D::s_ndim};
    }

    /**
    Get the shape of an "elemvec".

    \param arg The vector dimension.
    \returns [#nelem, #nne, tdim].
    */
    auto shape_elemvec(size_t arg) const -> std::array<size_t, 3>
    {
        return std::array<size_t, 3>{derived_cast().m_nelem, D::s_nne, arg};
    }

    /**
    Get the shape of an "elemmat".

    \returns [#nelem, #nne * #ndim, #nne * #ndim].
    */
    auto shape_elemmat() const -> std::array<size_t, 3>
    {
        return std::array<size_t, 3>{
            derived_cast().m_nelem, D::s_nne * D::s_ndim, D::s_nne * D::s_ndim};
    }

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \tparam rank The rank of the tensor.
        Since this function is templated, the output is fixed-size of type `std::array<size_t, n>`.

    \returns [#nelem, #nip, #tdim, ...].
    */
    template <size_t rank = 0>
    auto shape_qtensor() const -> std::array<size_t, 2 + rank>
    {
        std::array<size_t, 2 + rank> shape;
        shape[0] = derived_cast().m_nelem;
        shape[1] = derived_cast().m_nip;
        std::fill(shape.begin() + 2, shape.end(), derived_cast().m_tdim);
        return shape;
    }

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \param rank The tensor rank.
    \returns [#nelem, #nip, #tdim, ...].
    */
    auto shape_qtensor(size_t rank) const -> std::vector<size_t>
    {
        std::vector<size_t> shape(2 + rank);
        shape[0] = derived_cast().m_nelem;
        shape[1] = derived_cast().m_nip;
        std::fill(shape.begin() + 2, shape.end(), derived_cast().m_tdim);
        return shape;
    }

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \tparam rank The rank of the tensor.
        Since this function is templated, the output is fixed-size of type `std::array<size_t, n>`.

    \param rank The tensor rank.
        Effectively useless, but is there to distinguish from the dynamic-sized overloads.
    \param arg The tensor dimension.
    \returns [#nelem, #nip, tdim, ...].
    */
    template <size_t trank>
    auto shape_qtensor(size_t rank, size_t arg) const -> std::array<size_t, 2 + trank>
    {
        GOOSEFEM_ASSERT(trank == rank);
        std::array<size_t, 2 + trank> shape;
        shape[0] = derived_cast().m_nelem;
        shape[1] = derived_cast().m_nip;
        std::fill(shape.begin() + 2, shape.end(), arg);
        return shape;
    }

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \param rank The tensor rank.
    \param arg The tensor dimension.
    \returns [#nelem, #nip, tdim, ...].
    */
    auto shape_qtensor(size_t rank, size_t arg) const -> std::vector<size_t>
    {
        std::vector<size_t> shape(2 + rank);
        shape[0] = derived_cast().m_nelem;
        shape[1] = derived_cast().m_nip;
        std::fill(shape.begin() + 2, shape.end(), arg);
        return shape;
    }

    /**
    Get the shape of a "qscalar" (a "qtensor" of rank 0)
    \returns [#nelem, #nip].
    */
    auto shape_qscalar() const -> std::array<size_t, 2>
    {
        return std::array<size_t, 2>{derived_cast().m_nelem, derived_cast().m_nip};
    }

    /**
    Get the shape of a "qvector" (a "qtensor" of rank 1)
    \returns [#nelem, #nip, #tdim].
    */
    auto shape_qvector() const -> std::array<size_t, 3>
    {
        return std::array<size_t, 3>{derived_cast().m_nelem, derived_cast().m_nip, D::s_tdim};
    }

    /**
    Get the shape of a "qvector" (a "qtensor" of rank 1)
    \param arg Tensor dimension.
    \returns [#nelem, #nip, #tdim].
    */
    auto shape_qvector(size_t arg) const -> std::array<size_t, 3>
    {
        return std::array<size_t, 3>{derived_cast().m_nelem, derived_cast().m_nip, arg};
    }

    /**
    Get an allocated `array_type::tensor` to store a "elemvec".
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class R>
    auto allocate_elemvec() const
    {
        return array_type::tensor<R, 3>::from_shape(this->shape_elemvec());
    }

    /**
    Get an allocated and initialised `xt::xarray` to store a "elemvec".

    \tparam R value-type of the array, e.g. `double`.
    \param val The value to which to initialise all items.
    \returns `array_type::tensor` container of the correct shape.
    */
    template <class R>
    auto allocate_elemvec(R val) const
    {
        auto ret = array_type::tensor<R, 3>::from_shape(this->shape_elemvec());
        ret.fill(val);
        return ret;
    }

    /**
    Get an allocated `array_type::tensor` to store a "elemmat".
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class R>
    auto allocate_elemmat() const
    {
        return array_type::tensor<R, 3>::from_shape(this->shape_elemmat());
    }

    /**
    Get an allocated and initialised `xt::xarray` to store a "elemmat".

    \tparam R value-type of the array, e.g. `double`.
    \param val The value to which to initialise all items.
    \returns `array_type::tensor` container of the correct shape.
    */
    template <class R>
    auto allocate_elemmat(R val) const
    {
        auto ret = array_type::tensor<R, 3>::from_shape(this->shape_elemmat());
        ret.fill(val);
        return ret;
    }

    /**
    Get an allocated `array_type::tensor` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \returns [#nelem, #nip].
    */
    template <size_t rank = 0, class R>
    auto allocate_qtensor() const
    {
        return array_type::tensor<R, 2 + rank>::from_shape(this->shape_qtensor<rank>());
    }

    /**
    Get an allocated and initialised `array_type::tensor` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.

    \tparam R value-type of the array, e.g. `double`.
    \param val The value to which to initialise all items.
    \returns `array_type::tensor` container of the correct shape (and rank).
    */
    template <size_t rank = 0, class R>
    auto allocate_qtensor(R val) const
    {
        auto ret = array_type::tensor<R, 2 + rank>::from_shape(this->shape_qtensor<rank>());
        ret.fill(val);
        return ret;
    }

    /**
    Get an allocated `xt::xarray` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \param rank The tensor rank.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class R>
    auto allocate_qtensor(size_t rank) const
    {
        return xt::xarray<R>::from_shape(this->shape_qtensor(rank));
    }

    /**
    Get an allocated and initialised `xt::xarray` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \tparam R value-type of the array, e.g. `double`.
    \param rank The tensor rank.
    \param val The value to which to initialise all items.
    \returns `array_type::tensor` container of the correct shape (and rank).
    */
    template <class R>
    auto allocate_qtensor(size_t rank, R val) const
    {
        auto ret = xt::xarray<R>::from_shape(this->shape_qtensor(rank));
        ret.fill(val);
        return ret;
    }

    /**
    Get an allocated `array_type::tensor` to store a "qscalar" (a "qtensor" of rank 0).
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class R>
    auto allocate_qscalar() const
    {
        return this->allocate_qtensor<0, R>();
    }

    /**
    Get an allocated and initialised `xt::xarray` to store a "qscalar" (a "qtensor" of rank 0).

    \tparam R value-type of the array, e.g. `double`.
    \param val The value to which to initialise all items.
    \returns `array_type::tensor` container of the correct shape (and rank).
    */
    template <class R>
    auto allocate_qscalar(R val) const
    {
        return this->allocate_qtensor<0, R>(val);
    }

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }
};

/**
CRTP base class for interpolation and quadrature for a generic element in Cartesian coordinates.

Naming convention:
-    ``elemmat``:  matrices stored per element, [#nelem, #nne * #ndim, #nne * #ndim]
-    ``elemvec``:  nodal vectors stored per element, [#nelem, #nne, #ndim]
-    ``qtensor``:  integration point tensor, [#nelem, #nip, #tdim, #tdim]
-    ``qscalar``:  integration point scalar, [#nelem, #nip]
*/
template <class D>
class QuadratureBaseCartesian : public QuadratureBase<D> {
public:
    /**
    Underlying type.
    */
    using derived_type = D;

    /**
    Update the nodal positions.
    This recomputes:
    -   the shape functions,
    -   the shape function gradients (in local and global) coordinates,
    -   the integration points volumes.
    Under the small deformations assumption this function should not be called.

    \param x nodal coordinates (``elemvec``). Shape should match the earlier definition.
    */
    template <class T>
    void update_x(const T& x)
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, derived_cast().m_x.shape()));
        xt::noalias(derived_cast().m_x) = x;
        derived_cast().compute_dN_impl();
    }

    /**
    Shape function gradients (in global coordinates).
    \return ``gradN`` stored per element, per integration point [#nelem, #nip, #nne, #ndim].
    */
    auto GradN() const -> const array_type::tensor<double, 4>&
    {
        return derived_cast().m_dNx;
    }

    /**
    Integration volume.
    \return volume stored per element, per integration point [#nelem, #nip].
    */
    auto dV() const -> const array_type::tensor<double, 2>&
    {
        return derived_cast().m_vol;
    }

    /**
    Interpolate element vector and evaluate at each quadrature point.

    \f$ \vec{u}(\vec{x}_q) = N_i^e(\vec{x}) \vec{u}_i^e \f$

    \param elemvec nodal vector stored per element [#nelem, #nne, #ndim].
    \return qvector [#nelem, #nip, #ndim].
    */
    template <class T>
    auto InterpQuad_vector(const T& elemvec) const -> array_type::tensor<double, 3>
    {
        size_t n = elemvec.shape(2);
        auto qvector = array_type::tensor<double, 3>::from_shape(this->shape_qvector(n));
        derived_cast().interpQuad_vector_impl(elemvec, qvector);
        return qvector;
    }

    /**
    Same as InterpQuad_vector(), but writing to preallocated return.

    \param elemvec nodal vector stored per element [#nelem, #nne, #ndim].
    \param qvector [#nelem, #nip, #ndim].
    */
    template <class T, class R>
    void interpQuad_vector(const T& elemvec, R& qvector) const
    {
        derived_cast().interpQuad_vector_impl(elemvec, qvector);
    }

    /**
    Element-by-element: dyadic product of the shape function gradients and a nodal vector.
    Typical input: nodal displacements. Typical output: quadrature point strains.
    Within one element::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    qtensor(e, q, i, j) += dNdx(e, q, m, i) * elemvec(e, m, j)

    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \param elemvec [#nelem, #nne, #ndim]
    \return qtensor [#nelem, #nip, #tdim, #tdim]
    */
    template <class T>
    auto GradN_vector(const T& elemvec) const -> array_type::tensor<double, 4>
    {
        auto qtensor = array_type::tensor<double, 4>::from_shape(this->template shape_qtensor<2>());
        derived_cast().gradN_vector_impl(elemvec, qtensor);
        return qtensor;
    }

    /**
    Same as GradN_vector(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    template <class T, class R>
    void gradN_vector(const T& elemvec, R& qtensor) const
    {
        derived_cast().gradN_vector_impl(elemvec, qtensor);
    }

    /**
    The transposed output of GradN_vector().
    Within one element::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    qtensor(e, q, j, i) += dNdx(e, q, m, i) * elemvec(e, m, j)

    \param elemvec [#nelem, #nne, #ndim]
    \return qtensor [#nelem, #nip, #tdim, #tdim]
    */
    template <class T>
    auto GradN_vector_T(const T& elemvec) const -> array_type::tensor<double, 4>
    {
        auto qtensor = array_type::tensor<double, 4>::from_shape(this->template shape_qtensor<2>());
        derived_cast().gradN_vector_T_impl(elemvec, qtensor);
        return qtensor;
    }

    /**
    Same as GradN_vector_T(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    template <class T, class R>
    void gradN_vector_T(const T& elemvec, R& qtensor) const
    {
        derived_cast().gradN_vector_T_impl(elemvec, qtensor);
    }

    /**
    The symmetric output of GradN_vector().
    Without one element::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    qtensor(e, q, i, j) += 0.5 * dNdx(e, q, m, i) * elemvec(e, m, j)
                    qtensor(e, q, j, i) += 0.5 * dNdx(e, q, m, i) * elemvec(e, m, j)

    \param elemvec [#nelem, #nne, #ndim]
    \return qtensor [#nelem, #nip, #tdim, #tdim]
    */
    template <class T>
    auto SymGradN_vector(const T& elemvec) const -> array_type::tensor<double, 4>
    {
        auto qtensor = array_type::tensor<double, 4>::from_shape(this->template shape_qtensor<2>());
        derived_cast().symGradN_vector_impl(elemvec, qtensor);
        return qtensor;
    }

    /**
    Same as SymGradN_vector(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    template <class T, class R>
    void symGradN_vector(const T& elemvec, R& qtensor) const
    {
        derived_cast().symGradN_vector_impl(elemvec, qtensor);
    }

    /**
    Element-by-element: integral of a continuous vector-field.

    \f$ \vec{f}_i^e = \int N_i^e(\vec{x}) \vec{f}(\vec{x}) d\Omega_e \f$

    which is integration numerically as follows

    \f$ \vec{f}_i^e = \sum\limits_q N_i^e(\vec{x}_q) \vec{f}(\vec{x}_q) \f$

    \param qvector [#nelem, #nip. #ndim]
    \return elemvec [#nelem, #nne. #ndim]
    */
    template <class T>
    auto Int_N_vector_dV(const T& qvector) const -> array_type::tensor<double, 3>
    {
        size_t n = qvector.shape(2);
        auto elemvec = array_type::tensor<double, 3>::from_shape(this->shape_elemvec(n));
        derived_cast().int_N_vector_dV_impl(qvector, elemvec);
        return elemvec;
    }

    /**
    Same as Int_N_vector_dV(), but writing to preallocated return.

    \param qvector [#nelem, #nip. #ndim]
    \param elemvec overwritten [#nelem, #nne. #ndim]
    */
    template <class T, class R>
    void int_N_vector_dV(const T& qvector, R& elemvec) const
    {
        derived_cast().int_N_vector_dV_impl(qvector, elemvec);
    }

    /**
    Element-by-element: integral of the scalar product of the shape function with a scalar.
    Within one one element::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    for n in range(nne):
                        elemmat(e, m * ndim + i, n * ndim + i) +=
                            N(e, q, m) * qscalar(e, q) * N(e, q, n) * dV(e, q)

    with ``i`` a tensor dimension.
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \param qscalar [#nelem, #nip]
    \return elemmat [#nelem, #nne * #ndim, #nne * #ndim]
    */
    template <class T>
    auto Int_N_scalar_NT_dV(const T& qscalar) const -> array_type::tensor<double, 3>
    {
        auto elemmat = array_type::tensor<double, 3>::from_shape(this->shape_elemmat());
        derived_cast().int_N_scalar_NT_dV_impl(qscalar, elemmat);
        return elemmat;
    }

    /**
    Same as Int_N_scalar_NT_dV(), but writing to preallocated return.

    \param qscalar [#nelem, #nip]
    \param elemmat overwritten [#nelem, #nne * #ndim, #nne * #ndim]
    */
    template <class T, class R>
    void int_N_scalar_NT_dV(const T& qscalar, R& elemmat) const
    {
        derived_cast().int_N_scalar_NT_dV_impl(qscalar, elemmat);
    }

    /**
    Element-by-element: integral of the dot product of the shape function gradients with
    a second order tensor. Typical input: stress. Typical output: nodal force.
    Within one one element::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    elemvec(e, m, j) += dNdx(e, q, m, i) * qtensor(e, q, i, j) * dV(e, q)

    with ``i`` and ``j`` tensor dimensions.
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \param qtensor [#nelem, #nip, #ndim, #ndim]
    \return elemvec [#nelem, #nne. #ndim]
    */
    template <class T>
    auto Int_gradN_dot_tensor2_dV(const T& qtensor) const -> array_type::tensor<double, 3>
    {
        auto elemvec = array_type::tensor<double, 3>::from_shape(this->shape_elemvec());
        derived_cast().int_gradN_dot_tensor2_dV_impl(qtensor, elemvec);
        return elemvec;
    }

    /**
    Same as Int_gradN_dot_tensor2_dV(), but writing to preallocated return.

    \param qtensor [#nelem, #nip, #ndim, #ndim]
    \param elemvec overwritten [#nelem, #nne. #ndim]
    */
    template <class T, class R>
    void int_gradN_dot_tensor2_dV(const T& qtensor, R& elemvec) const
    {
        derived_cast().int_gradN_dot_tensor2_dV_impl(qtensor, elemvec);
    }

    // Integral of the dot product
    // elemmat(m*2+j, n*2+k) += dNdx(m,i) * qtensor(i,j,k,l) * dNdx(n,l) * dV

    /**
    Element-by-element: integral of the dot products of the shape function gradients with
    a fourth order tensor. Typical input: stiffness tensor. Typical output: stiffness matrix.
    Within one one element::

        for e in range(nelem):
            for q in range(nip):
                for m in range(nne):
                    for n in range(nne):
                        elemmat(e, m * ndim + j, n * ndim + k) +=
                            dNdx(e,q,m,i) * qtensor(e,q,i,j,k,l) * dNdx(e,q,n,l) * dV(e,q)

    with ``i``, ``j``, ``k``, and ``l`` tensor dimensions.
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \param qtensor [#nelem, #nip, #ndim, #ndim, #ndim, #ndim]
    \return elemmat [#nelem, #nne * #ndim, #nne * #ndim]
    */
    template <class T>
    auto Int_gradN_dot_tensor4_dot_gradNT_dV(const T& qtensor) const
        -> array_type::tensor<double, 3>
    {
        auto elemmat = array_type::tensor<double, 3>::from_shape(this->shape_elemmat());
        derived_cast().int_gradN_dot_tensor4_dot_gradNT_dV_impl(qtensor, elemmat);
        return elemmat;
    }

    /**
    Same as Int_gradN_dot_tensor4_dot_gradNT_dV(), but writing to preallocated return.

    \param qtensor [#nelem, #nip, #ndim, #ndim, #ndim, #ndim]
    \param elemmat overwritten [#nelem, #nne * #ndim, #nne * #ndim]
    */
    template <class T, class R>
    void int_gradN_dot_tensor4_dot_gradNT_dV(const T& qtensor, R& elemmat) const
    {
        derived_cast().int_gradN_dot_tensor4_dot_gradNT_dV_impl(qtensor, elemmat);
    }

protected:
    /**
    Update the shape function gradients (called when the nodal positions are updated).
    */
    void compute_dN()
    {
        derived_cast().compute_dN_impl();
    }

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

    friend class QuadratureBase<D>;

    template <class T, class R>
    void interpQuad_vector_impl(const T& elemvec, R& qvector) const
    {
        size_t n = elemvec.shape(2);
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec(n)));
        GOOSEFEM_ASSERT(xt::has_shape(qvector, this->shape_qvector(n)));

        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& N = derived_cast().m_N;

        qvector.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < nelem; ++e) {

            auto fq = &elemvec(e, 0, 0);

            for (size_t q = 0; q < nip; ++q) {

                auto Nq = &N(q, 0);
                auto tq = &qvector(e, q, 0);

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t i = 0; i < n; ++i) {
                        tq[i] += Nq[m] * fq[m * n + i];
                    }
                }
            }
        }
    }

    template <class T, class R>
    void gradN_vector_impl(const T& elemvec, R& qtensor) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<2>()));

        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& dNx = derived_cast().m_dNx;

        qtensor.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < nelem; ++e) {

            auto ue = xt::adapt(&elemvec(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

            for (size_t q = 0; q < nip; ++q) {

                auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
                auto graduq = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<D::s_tdim, D::s_tdim>());

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t i = 0; i < D::s_ndim; ++i) {
                        for (size_t j = 0; j < D::s_ndim; ++j) {
                            graduq(i, j) += dNxq(m, i) * ue(m, j);
                        }
                    }
                }
            }
        }
    }

    template <class T, class R>
    void gradN_vector_T_impl(const T& elemvec, R& qtensor) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<2>()));

        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& dNx = derived_cast().m_dNx;

        qtensor.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < nelem; ++e) {

            auto ue = xt::adapt(&elemvec(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

            for (size_t q = 0; q < nip; ++q) {

                auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
                auto graduq = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<D::s_tdim, D::s_tdim>());

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t i = 0; i < D::s_ndim; ++i) {
                        for (size_t j = 0; j < D::s_ndim; ++j) {
                            graduq(j, i) += dNxq(m, i) * ue(m, j);
                        }
                    }
                }
            }
        }
    }

    template <class T, class R>
    void symGradN_vector_impl(const T& elemvec, R& qtensor) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<2>()));

        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& dNx = derived_cast().m_dNx;

        qtensor.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < nelem; ++e) {

            auto ue = xt::adapt(&elemvec(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

            for (size_t q = 0; q < nip; ++q) {

                auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
                auto epsq = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<D::s_tdim, D::s_tdim>());

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t i = 0; i < D::s_ndim; ++i) {
                        for (size_t j = 0; j < D::s_ndim; ++j) {
                            epsq(i, j) += 0.5 * dNxq(m, i) * ue(m, j);
                            epsq(j, i) += 0.5 * dNxq(m, i) * ue(m, j);
                        }
                    }
                }
            }
        }
    }

    template <class T, class R>
    void int_N_vector_dV_impl(const T& qvector, R& elemvec) const
    {
        size_t n = qvector.shape(2);
        GOOSEFEM_ASSERT(xt::has_shape(qvector, this->shape_qvector(n)));
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec(n)));

        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& N = derived_cast().m_N;
        auto& vol = derived_cast().m_vol;

        elemvec.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < nelem; ++e) {

            auto f = &elemvec(e, 0, 0);

            for (size_t q = 0; q < nip; ++q) {

                auto Ne = &N(q, 0);
                auto tq = &qvector(e, q, 0);
                auto& volq = vol(e, q);

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t i = 0; i < n; ++i) {
                        f[m * n + i] += Ne[m] * tq[i] * volq;
                    }
                }
            }
        }
    }

    template <class T, class R>
    void int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(qscalar, this->shape_qscalar()));
        GOOSEFEM_ASSERT(xt::has_shape(elemmat, this->shape_elemmat()));

        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& N = derived_cast().m_N;
        auto& vol = derived_cast().m_vol;

        elemmat.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < nelem; ++e) {

            auto Me = xt::adapt(
                &elemmat(e, 0, 0), xt::xshape<D::s_nne * D::s_ndim, D::s_nne * D::s_ndim>());

            for (size_t q = 0; q < nip; ++q) {

                auto Ne = xt::adapt(&N(q, 0), xt::xshape<D::s_nne>());
                auto& volq = vol(e, q);
                auto& rho = qscalar(e, q);

                // M(m * D::s_ndim + i, n * D::s_ndim + i) += N(m) * scalar * N(n) * dV
                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t n = 0; n < D::s_nne; ++n) {
                        for (size_t i = 0; i < D::s_ndim; ++i) {
                            Me(m * D::s_ndim + i, n * D::s_ndim + i) += Ne(m) * rho * Ne(n) * volq;
                        }
                    }
                }
            }
        }
    }

    template <class T, class R>
    void int_gradN_dot_tensor2_dV_impl(const T& qtensor, R& elemvec) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<2>()));
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, this->shape_elemvec()));

        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& dNx = derived_cast().m_m_dNx;
        auto& vol = derived_cast().m_vol;

        elemvec.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < nelem; ++e) {

            auto fe = xt::adapt(&elemvec(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

            for (size_t q = 0; q < nip; ++q) {

                auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
                auto sigq = xt::adapt(&qtensor(e, q, 0, 0), xt::xshape<D::s_tdim, D::s_tdim>());
                auto& volq = vol(e, q);

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t i = 0; i < D::s_ndim; ++i) {
                        for (size_t j = 0; j < D::s_ndim; ++j) {
                            fe(m, j) += dNxq(m, i) * sigq(i, j) * volq;
                        }
                    }
                }
            }
        }
    }

    template <class T, class R>
    void int_gradN_dot_tensor4_dot_gradNT_dV_impl(const T& qtensor, R& elemmat) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(qtensor, this->template shape_qtensor<4>()));
        GOOSEFEM_ASSERT(xt::has_shape(elemmat, this->shape_elemmat()));

        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& dNx = derived_cast().m_dNx;
        auto& vol = derived_cast().m_vol;

        elemmat.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < nelem; ++e) {

            auto K = xt::adapt(
                &elemmat(e, 0, 0), xt::xshape<D::s_nne * D::s_ndim, D::s_nne * D::s_ndim>());

            for (size_t q = 0; q < nip; ++q) {

                auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
                auto Cq = xt::adapt(
                    &qtensor(e, q, 0, 0, 0, 0),
                    xt::xshape<D::s_tdim, D::s_tdim, D::s_tdim, D::s_tdim>());
                auto& volq = vol(e, q);

                for (size_t m = 0; m < D::s_nne; ++m) {
                    for (size_t n = 0; n < D::s_nne; ++n) {
                        for (size_t i = 0; i < D::s_ndim; ++i) {
                            for (size_t j = 0; j < D::s_ndim; ++j) {
                                for (size_t k = 0; k < D::s_ndim; ++k) {
                                    for (size_t l = 0; l < D::s_ndim; ++l) {
                                        K(m * D::s_ndim + j, n * D::s_ndim + k) +=
                                            dNxq(m, i) * Cq(i, j, k, l) * dNxq(n, l) * volq;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void compute_dN_impl()
    {
        auto nelem = derived_cast().m_nelem;
        auto nip = derived_cast().m_nip;
        auto& vol = derived_cast().m_vol;
        auto& w = derived_cast().m_w;
        auto& dNxi = derived_cast().m_dNxi;
        auto& dNx = derived_cast().m_dNx;
        auto& x = derived_cast().m_x;

        dNx.fill(0.0);

#pragma omp parallel
        {
            auto J = array_type::tensor<double, 2>::from_shape({D::s_ndim, D::s_ndim});
            auto Jinv = array_type::tensor<double, 2>::from_shape({D::s_ndim, D::s_ndim});

#pragma omp for
            for (size_t e = 0; e < nelem; ++e) {

                auto xe = xt::adapt(&x(e, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

                for (size_t q = 0; q < nip; ++q) {

                    auto dNxiq = xt::adapt(&dNxi(q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());
                    auto dNxq = xt::adapt(&dNx(e, q, 0, 0), xt::xshape<D::s_nne, D::s_ndim>());

                    J.fill(0.0);

                    for (size_t m = 0; m < D::s_nne; ++m) {
                        for (size_t i = 0; i < D::s_ndim; ++i) {
                            for (size_t j = 0; j < D::s_ndim; ++j) {
                                J(i, j) += dNxiq(m, i) * xe(m, j);
                            }
                        }
                    }

                    double Jdet = detail::tensor<D::s_ndim>::inv(J, Jinv);

                    for (size_t m = 0; m < D::s_nne; ++m) {
                        for (size_t i = 0; i < D::s_ndim; ++i) {
                            for (size_t j = 0; j < D::s_ndim; ++j) {
                                dNxq(m, i) += Jinv(i, j) * dNxiq(m, i);
                            }
                        }
                    }

                    vol(e, q) = w(q) * Jdet;
                }
            }
        }
    }
};

} // namespace Element
} // namespace GooseFEM

#endif

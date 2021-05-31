/**
Convenience methods for integration point data.

\file Element.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ELEMENT_H
#define GOOSEFEM_ELEMENT_H

#include "config.h"
#include "Allocate.h"

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
inline xt::xtensor<double, 3> asElementVector(
    const xt::xtensor<size_t, 2>& conn, const xt::xtensor<double, 2>& nodevec);

/**
Assemble nodal vector stored per element ("elemvec", shape ``[nelem, nne, ndim]``) to nodal vector
("nodevec", shape ``[nnode, ndim]``).

\param conn Connectivity.
\param elemvec "elemvec".
\return "nodevec".
*/
inline xt::xtensor<double, 2> assembleNodeVector(
    const xt::xtensor<size_t, 2>& conn, const xt::xtensor<double, 3>& elemvec);

/**
Check that DOFs leave no holes.

\param dofs DOFs ("nodevec")
\return ``true`` if there are no holds.
*/
template <class E>
inline bool isSequential(const E& dofs);


/**
Check that all of the matrices stored per elemmat (shape: ``[nelem, nne * ndim, nne * ndim]``)
are diagonal.

\param elemmat Element-vectors ("elemmat")
\return ``true`` if all element matrices are diagonal.
*/
bool isDiagonal(const xt::xtensor<double, 3>& elemmat);

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
    auto nelem() const;

    /**
    Number of nodes per element.

    \return Scalar.
    */
    auto nne() const;

    /**
    Number of dimensions for node vectors.

    \return Scalar.
    */
    auto ndim() const;

    /**
    Number of dimensions for integration point tensors.

    \return Scalar.
    */
    auto tdim() const;

    /**
    Number of integration points.

    \return Scalar.
    */
    auto nip() const;

    /**
    Convert "qscalar" to "qtensor" of certain rank.
    Fully allocated output passed as reference, use AsTensor to allocate and return data.

    \param qscalar A "qscalar".
    \param qtensor A "qtensor".
    */
    template <class T, class R>
    void asTensor(const T& qscalar, R& qtensor) const;

    /**
    Convert "qscalar" to "qtensor" of certain rank.

    \param "qscalar".
    \return "qtensor".
    */
    template <size_t rank, class T>
    auto AsTensor(const T& qscalar) const;

    /**
    Convert "qscalar" to "qtensor" of certain rank.

    \param rank Tensor rank.
    \param qscalar A "qscalar".
    \return "qtensor".
    */
    template <class T>
    auto AsTensor(size_t rank, const T& qscalar) const;

    /**
    Get the shape of an "elemvec".

    \returns [#nelem, #nne, #ndim].
    */
    auto shape_elemvec() const -> std::array<size_t, 3>;

    /**
    Get the shape of an "elemvec".

    \param tdim The vector dimension.
    \returns [#nelem, #nne, tdim].
    */
    auto shape_elemvec(size_t tdim) const -> std::array<size_t, 3>;

    /**
    Get the shape of an "elemmat".

    \returns [#nelem, #nne * #ndim, #nne * #ndim].
    */
    auto shape_elemmat() const -> std::array<size_t, 3>;

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \tparam rank The rank of the tensor.
        Since this function is templated, the output is fixed-size of type `std::array<size_t, n>`.

    \returns [#nelem, #nip, #tdim, ...].
    */
    template <size_t rank = 0>
    auto shape_qtensor() const -> std::array<size_t, 2 + rank>;

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \param rank The tensor rank.
    \returns [#nelem, #nip, #tdim, ...].
    */
    auto shape_qtensor(size_t rank) const -> std::vector<size_t>;

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \tparam rank The rank of the tensor.
        Since this function is templated, the output is fixed-size of type `std::array<size_t, n>`.

    \param rank The tensor rank.
        Effectively useless, but is there to distinguish from the dynamic-sized overloads.
    \param tdim The tensor dimension.
    \returns [#nelem, #nip, tdim, ...].
    */
    template <size_t trank>
    auto shape_qtensor(size_t rank, size_t tdim) const -> std::array<size_t, 2 + trank>;

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \param rank The tensor rank.
    \param tdim The tensor dimension.
    \returns [#nelem, #nip, tdim, ...].
    */
    auto shape_qtensor(size_t rank, size_t tdim) const -> std::vector<size_t>;

    /**
    Get the shape of a "qscalar" (a "qtensor" of rank 0)
    \returns [#nelem, #nip].
    */
    auto shape_qscalar() const -> std::array<size_t, 2>;

    /**
    Get the shape of a "qvector" (a "qtensor" of rank 1)
    \returns [#nelem, #nip, #tdim].
    */
    auto shape_qvector() const -> std::array<size_t, 3>;

    /**
    Get the shape of a "qvector" (a "qtensor" of rank 1)
    \param tdim Tensor dimension.
    \returns [#nelem, #nip, #tdim].
    */
    auto shape_qvector(size_t tdim) const -> std::array<size_t, 3>;

    /**
    Get an allocated `xt::xtensor` to store a "elemvec".
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class R>
    auto allocate_elemvec() const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "elemvec".

    \tparam R value-type of the array, e.g. `double`.
    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape.
    */
    template <class R>
    auto allocate_elemvec(R val) const;

    /**
    Get an allocated `xt::xtensor` to store a "elemmat".
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class R>
    auto allocate_elemmat() const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "elemmat".

    \tparam R value-type of the array, e.g. `double`.
    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape.
    */
    template <class R>
    auto allocate_elemmat(R val) const;

    /**
    Get an allocated `xt::xtensor` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \returns [#nelem, #nip].
    */
    template <size_t rank = 0, class R>
    auto allocate_qtensor() const;

    /**
    Get an allocated and initialised `xt::xtensor` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.

    \tparam R value-type of the array, e.g. `double`.
    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    template <size_t rank = 0, class R>
    auto allocate_qtensor(R val) const;

    /**
    Get an allocated `xt::xarray` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \param rank The tensor rank.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class R>
    auto allocate_qtensor(size_t rank) const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \tparam R value-type of the array, e.g. `double`.
    \param rank The tensor rank.
    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    template <class R>
    auto allocate_qtensor(size_t rank, R val) const;

    /**
    Get an allocated `xt::xtensor` to store a "qscalar" (a "qtensor" of rank 0).
    Note: the container is not (zero-)initialised.

    \tparam R value-type of the array, e.g. `double`.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class R>
    auto allocate_qscalar() const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "qscalar" (a "qtensor" of rank 0).

    \tparam R value-type of the array, e.g. `double`.
    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    template <class R>
    auto allocate_qscalar(R val) const;

private:

    auto derived_cast() -> derived_type&;
    auto derived_cast() const -> const derived_type&;
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
    void update_x(const T& x);

    /**
    Get the shape function gradients (in global coordinates).

    \return ``gradN`` stored per element, per integration point [#nelem, #nip, #nne, #ndim].
    */
    auto GradN() const -> xt::xtensor<double, 4>;

    /**
    Get the integration volume.

    \return volume stored per element, per integration point [#nelem, #nip].
    */
    auto dV() const -> xt::xtensor<double, 2>;

    /**
    Interpolate element vector and evaluate at each quadrature point.

    \f$ \vec{u}(\vec{x}_q) = N_i^e(\vec{x}) \vec{u}_i^e \f$

    \param elemvec nodal vector stored per element [#nelem, #nne, #ndim].
    \return qvector [#nelem, #nip, #ndim].
    */
    template <class T>
    auto InterpQuad_vector(const T& elemvec) const -> xt::xtensor<double, 3>;

    /**
    Same as InterpQuad_vector(), but writing to preallocated return.

    \param elemvec nodal vector stored per element [#nelem, #nne, #ndim].
    \param qvector [#nelem, #nip, #ndim].
    */
    template <class T, class R>
    void interpQuad_vector(const T& elemvec, R& qvector) const;

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
    auto GradN_vector(const T& elemvec) const -> xt::xtensor<double, 4>;

    /**
    Same as GradN_vector(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    template <class T, class R>
    void gradN_vector(const T& elemvec, R& qtensor) const;

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
    auto GradN_vector_T(const T& elemvec) const -> xt::xtensor<double, 4>;

    /**
    Same as GradN_vector_T(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    template <class T, class R>
    void gradN_vector_T(const T& elemvec, R& qtensor) const;

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
    auto SymGradN_vector(const T& elemvec) const -> xt::xtensor<double, 4>;

    /**
    Same as SymGradN_vector(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    template <class T, class R>
    void symGradN_vector(const T& elemvec, R& qtensor) const;

    /**
    Element-by-element: integral of a continuous vector-field.

    \f$ \vec{f}_i^e = \int N_i^e(\vec{x}) \vec{f}(\vec{x}) d\Omega_e \f$

    which is integration numerically as follows

    \f$ \vec{f}_i^e = \sum\limits_q N_i^e(\vec{x}_q) \vec{f}(\vec{x}_q) \f$

    \param qvector [#nelem, #nip. #ndim]
    \return elemvec [#nelem, #nne. #ndim]
    */
    template <class T>
    auto Int_N_vector_dV(const T& qvector) const -> xt::xtensor<double, 3>;

    /**
    Same as Int_N_vector_dV(), but writing to preallocated return.

    \param qvector [#nelem, #nip. #ndim]
    \param elemvec overwritten [#nelem, #nne. #ndim]
    */
    template <class T, class R>
    void int_N_vector_dV(const T& qvector, R& elemvec) const;

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
    auto Int_N_scalar_NT_dV(const T& qscalar) const -> xt::xtensor<double, 3>;

    /**
    Same as Int_N_scalar_NT_dV(), but writing to preallocated return.

    \param qscalar [#nelem, #nip]
    \param elemmat overwritten [#nelem, #nne * #ndim, #nne * #ndim]
    */
    template <class T, class R>
    void int_N_scalar_NT_dV(const T& qscalar, R& elemmat) const;

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
    auto Int_gradN_dot_tensor2_dV(const T& qtensor) const -> xt::xtensor<double, 3>;

    /**
    Same as Int_gradN_dot_tensor2_dV(), but writing to preallocated return.

    \param qtensor [#nelem, #nip, #ndim, #ndim]
    \param elemvec overwritten [#nelem, #nne. #ndim]
    */
    template <class T, class R>
    void int_gradN_dot_tensor2_dV(const T& qtensor, R& elemvec) const;

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
                            dNdx(e, q, m, i) * qtensor(e, q, i, j, k, l) * dNdx(e, q, n, l) * dV(e, q)

    with ``i``, ``j``, ``k``, and ``l`` tensor dimensions.
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \param qtensor [#nelem, #nip, #ndim, #ndim, #ndim, #ndim]
    \return elemmat [#nelem, #nne * #ndim, #nne * #ndim]
    */
    template <class T>
    auto Int_gradN_dot_tensor4_dot_gradNT_dV(const T& qtensor) const -> xt::xtensor<double, 3>;

    /**
    Same as Int_gradN_dot_tensor4_dot_gradNT_dV(), but writing to preallocated return.

    \param qtensor [#nelem, #nip, #ndim, #ndim, #ndim, #ndim]
    \param elemmat overwritten [#nelem, #nne * #ndim, #nne * #ndim]
    */
    template <class T, class R>
    void int_gradN_dot_tensor4_dot_gradNT_dV(const T& qtensor, R& elemmat) const;

protected:

    /**
    Update the shape function gradients (called when the nodal positions are updated).
    */
    void compute_dN();

private:

    auto derived_cast() -> derived_type&;
    auto derived_cast() const -> const derived_type&;

    friend class QuadratureBase<D>;

    template <class T, class R>
    void interpQuad_vector_impl(const T& elemvec, R& qvector) const;

    template <class T, class R>
    void gradN_vector_impl(const T& elemvec, R& qtensor) const;

    template <class T, class R>
    void gradN_vector_T_impl(const T& elemvec, R& qtensor) const;

    template <class T, class R>
    void symGradN_vector_impl(const T& elemvec, R& qtensor) const;

    template <class T, class R>
    void int_N_vector_dV_impl(const T& qvector, R& elemvec) const;

    template <class T, class R>
    void int_N_scalar_NT_dV_impl(const T& qscalar, R& elemmat) const;

    template <class T, class R>
    void int_gradN_dot_tensor2_dV_impl(const T& qtensor, R& elemvec) const;

    template <class T, class R>
    void int_gradN_dot_tensor4_dot_gradNT_dV_impl(const T& qtensor, R& elemmat) const;

    void compute_dN_impl();
};

} // namespace Element
} // namespace GooseFEM

#include "Element.hpp"

#endif

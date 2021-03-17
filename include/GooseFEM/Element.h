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
Base quadrature-class.
This class does not have a specific element-type in mind, it is used mostly internally
to derive from such that common methods do not have to be reimplementation.

\tparam ne Number of nodes per element.
\tparam nd Number of dimensions for node vectors.
\tparam td Number of dimensions for integration point tensors.
*/
template <size_t ne, size_t nd, size_t td>
class QuadratureBase {
public:
    QuadratureBase() = default;

    /**
    Constructor
    */
    QuadratureBase(size_t nelem, size_t nip);

    /**
    Number of elements.

    \return Scalar.
    */
    size_t nelem() const;

    /**
    Number of nodes per element.

    \return Scalar.
    */
    size_t nne() const;

    /**
    Number of dimensions for node vectors.

    \return Scalar.
    */
    size_t ndim() const;

    /**
    Number of dimensions for integration point tensors.

    \return Scalar.
    */
    size_t tdim() const;

    /**
    Number of integration points.

    \return Scalar.
    */
    size_t nip() const;

    /**
    Convert "qscalar" to "qtensor" of certain rank.
    Fully allocated output passed as reference, use AsTensor to allocate and return data.

    \param qscalar A "qscalar".
    \param qtensor A "qtensor".
    */
    template <size_t rank = 0, class T>
    void asTensor(const xt::xtensor<T, 2>& qscalar, xt::xtensor<T, 2 + rank>& qtensor) const;

    /**
    Convert "qscalar" to "qtensor" of certain rank.

    \param "qscalar".
    \return "qtensor".
    */
    template <size_t rank = 0, class T>
    xt::xtensor<T, 2 + rank> AsTensor(const xt::xtensor<T, 2>& qscalar) const;

    /**
    Convert "qscalar" to "qtensor" of certain rank.

    \param rank Tensor rank.
    \param qscalar A "qscalar".
    \return "qtensor".
    */
    template <class T>
    xt::xarray<T> AsTensor(size_t rank, const xt::xtensor<T, 2>& qscalar) const;

    /**
    Get the shape of an "elemvec".

    \returns [#nelem, #nne, #ndim].
    */
    std::array<size_t, 3> shape_elemvec() const;

    /**
    Get the shape of an "elemmat".

    \returns [#nelem, #nne * #ndim, #nne * #ndim].
    */
    std::array<size_t, 3> shape_elemmat() const;

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.

    \returns [#nelem, #nip, #tdim, ...].
    */
    template <size_t rank = 0>
    std::array<size_t, 2 + rank> shape_qtensor() const;

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \param rank The tensor rank.
    \returns [#nelem, #nip, #tdim, ...].
    */
    std::vector<size_t> shape_qtensor(size_t rank) const;

    /**
    Get the shape of a "qscalar" (a "qtensor" of rank 0)

    \returns [#nelem, #nip].
    */
    std::vector<size_t> shape_qscalar() const;

    /**
    Get an allocated `xt::xtensor` to store a "elemvec".
    Note: the container is not (zero-)initialised.

    \returns `xt::xarray` container of the correct shape.
    */
    template <class T>
    xt::xtensor<T, 3> allocate_elemvec() const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "elemvec".

    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape.
    */
    template <class T>
    xt::xtensor<T, 3> allocate_elemvec(T val) const;

    /**
    Get an allocated `xt::xtensor` to store a "elemmat".
    Note: the container is not (zero-)initialised.

    \returns `xt::xarray` container of the correct shape.
    */
    template <class T>
    xt::xtensor<T, 3> allocate_elemmat() const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "elemmat".

    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape.
    */
    template <class T>
    xt::xtensor<T, 3> allocate_elemmat(T val) const;

    /**
    Get an allocated `xt::xtensor` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.
    Note: the container is not (zero-)initialised.

    \returns [#nelem, #nip].
    */
    template <size_t rank = 0, class T>
    xt::xtensor<T, 2 + rank> allocate_qtensor() const;

    /**
    Get an allocated and initialised `xt::xtensor` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.

    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    template <size_t rank = 0, class T>
    xt::xtensor<T, 2 + rank> allocate_qtensor(T val) const;

    /**
    Get an allocated `xt::xarray` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Note: the container is not (zero-)initialised.

    \param rank The tensor rank.
    \returns `xt::xarray` container of the correct shape.
    */
    template <class T>
    xt::xarray<T> allocate_qtensor(size_t rank) const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \param rank The tensor rank.
    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    template <class T>
    xt::xarray<T> allocate_qtensor(size_t rank, T val) const;

    /**
    Get an allocated `xt::xtensor` to store a "qscalar" (a "qtensor" of rank 0).
    Note: the container is not (zero-)initialised.

    \returns `xt::xarray` container of the correct shape.
    */
    template <class T>
    xt::xtensor<T, 2> allocate_qscalar() const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "qscalar" (a "qtensor" of rank 0).

    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    template <class T>
    xt::xtensor<T, 2> allocate_qscalar(T val) const;

    /**
    \cond
    */
    template <size_t rank = 0>
    [[ deprecated ]]
    std::array<size_t, 2 + rank> ShapeQtensor() const;

    [[ deprecated ]]
    std::vector<size_t> ShapeQtensor(size_t rank) const;

    template <size_t rank = 0, class T>
    [[ deprecated ]]
    xt::xtensor<T, 2 + rank> AllocateQtensor() const;

    [[ deprecated ]]
    std::vector<size_t> ShapeQscalar() const;

    template <size_t rank = 0, class T>
    [[ deprecated ]]
    xt::xtensor<T, 2 + rank> AllocateQtensor(T val) const;

    template <class T>
    [[ deprecated ]]
    xt::xarray<T> AllocateQtensor(size_t rank) const;

    template <class T>
    [[ deprecated ]]
    xt::xarray<T> AllocateQtensor(size_t rank, T val) const;

    template <class T>
    [[ deprecated ]]
    xt::xtensor<T, 2> AllocateQscalar() const;

    template <class T>
    [[ deprecated ]]
    xt::xtensor<T, 2> AllocateQscalar(T val) const;
    /**
    \endcond
    */

protected:
    /**
    Wrapper of constructor, for derived classes.
    */
    void initQuadratureBase(size_t nelem, size_t nip);

protected:
    size_t m_nelem; ///< Number of elements.
    size_t m_nip; ///< Number of integration points per element.
    constexpr static size_t m_nne = ne; ///< Number of nodes per element.
    constexpr static size_t m_ndim = nd; ///< Number of dimensions for nodal vectors.
    constexpr static size_t m_tdim = td; ///< Number of dimensions for integration point tensors.
};

/**
Interpolation and quadrature for a generic element in Cartesian coordinates.

Naming convention:
-    ``elemmat``:  matrices stored per element, [#nelem, #nne * #ndim, #nne * #ndim]
-    ``elemvec``:  nodal vectors stored per element, [#nelem, #nne, #ndim]
-    ``qtensor``:  integration point tensor, [#nelem, #nip, #tdim, #tdim]
-    ``qscalar``:  integration point scalar, [#nelem, #nip]
*/
template <size_t ne, size_t nd, size_t td>
class QuadratureBaseCartesian : public QuadratureBase<ne, nd, td> {
public:

    QuadratureBaseCartesian() = default;

    virtual ~QuadratureBaseCartesian(){};

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
    \param N Shape functions in local coordinates [#nip, #nne].
    \param dNdxi Shape function gradient w.r.t. local coordinates [#nip, #nne, #ndim].
    */
    QuadratureBaseCartesian(
        const xt::xtensor<double, 3>& x,
        const xt::xtensor<double, 2>& xi,
        const xt::xtensor<double, 1>& w,
        const xt::xtensor<double, 2>& N,
        const xt::xtensor<double, 3>& dNdxi);

    /**
    Update the nodal positions.
    This recomputes the values of the shape functions and the shape function gradients
    (in local and global) coordinates.
    Under the small deformations assumption this should not be called.

    \param x nodal coordinates (``elemvec``). Shape should match the earlier definition.
    */
    void update_x(const xt::xtensor<double, 3>& x);

    /**
    Get the shape function gradients (in global coordinates).
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \return ``gradN`` stored per element, per integration point [#nelem, #nip, #nne, #ndim].
    */
    virtual xt::xtensor<double, 4> GradN() const;

    /**
    Get the integration volume.
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \return volume stored per element, per integration point [#nelem, #nip].
    */
    xt::xtensor<double, 2> dV() const;

    /**
    Interpolate element vector.
    Note that the functions and their gradients are precomputed upon construction,
    or updated when calling update_x().

    \param elemvec nodal vector stored per element (shape: [#nelem, #nne, #ndim]).
    \return integration point vector (shape: [#nelem, #nip, #ndim]).
    */
    xt::xtensor<double, 3> Interp_N_vector(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Same as Interp_N_vector(), but writing to preallocated return.

    \param elemvec nodal vector stored per element (shape: [#nelem, #nne, #ndim]).
    \param qvector integration point vector, overwritten (shape: [#nelem, #nip, #ndim]).
    */
    virtual void interp_N_vector(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 3>& qvector) const;

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
    xt::xtensor<double, 4> GradN_vector(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Same as GradN_vector(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    virtual void gradN_vector(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;

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
    xt::xtensor<double, 4> GradN_vector_T(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Same as GradN_vector_T(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    virtual void gradN_vector_T(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;

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
    xt::xtensor<double, 4> SymGradN_vector(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Same as SymGradN_vector(), but writing to preallocated return.

    \param elemvec [#nelem, #nne, #ndim]
    \param qtensor overwritten [#nelem, #nip, #tdim, #tdim]
    */
    virtual void symGradN_vector(
        const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;

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
    xt::xtensor<double, 3> Int_N_scalar_NT_dV(const xt::xtensor<double, 2>& qscalar) const;

    /**
    Same as Int_N_scalar_NT_dV(), but writing to preallocated return.

    \param qscalar [#nelem, #nip]
    \param elemmat overwritten [#nelem, #nne * #ndim, #nne * #ndim]
    */
    virtual void int_N_scalar_NT_dV(
        const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 3>& elemmat) const;

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
    xt::xtensor<double, 3> Int_gradN_dot_tensor2_dV(const xt::xtensor<double, 4>& qtensor) const;

    /**
    Same as Int_gradN_dot_tensor2_dV(), but writing to preallocated return.

    \param qtensor [#nelem, #nip, #ndim, #ndim]
    \param elemvec overwritten [#nelem, #nne. #ndim]
    */
    virtual void int_gradN_dot_tensor2_dV(
        const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const;

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
    xt::xtensor<double, 3> Int_gradN_dot_tensor4_dot_gradNT_dV(
        const xt::xtensor<double, 6>& qtensor) const;

    /**
    Same as Int_gradN_dot_tensor4_dot_gradNT_dV(), but writing to preallocated return.

    \param qtensor [#nelem, #nip, #ndim, #ndim, #ndim, #ndim]
    \param elemmat overwritten [#nelem, #nne * #ndim, #nne * #ndim]
    */
    virtual void int_gradN_dot_tensor4_dot_gradNT_dV(
        const xt::xtensor<double, 6>& qtensor, xt::xtensor<double, 3>& elemmat) const;

protected:
    /**
    Constructor alias.
    \param x nodal coordinates (``elemvec``).
    \param xi Integration point coordinates (local coordinates) [#nip].
    \param w Integration point weights [#nip].
    \param N Shape functions in the integration points [#nip, #nne].
    \param dNdxi Shape function gradient w.r.t. local coordinates [#nip, #nne, #ndim].
    */
    void initQuadratureBaseCartesian(
        const xt::xtensor<double, 3>& x,
        const xt::xtensor<double, 2>& xi,
        const xt::xtensor<double, 1>& w,
        const xt::xtensor<double, 2>& N,
        const xt::xtensor<double, 3>& dNdxi);

    /**
    Compute m_vol and m_dNdx based on current m_x.
    */
    virtual void compute_dN();

protected:
    using QuadratureBase<ne, nd, td>::m_nelem;
    using QuadratureBase<ne, nd, td>::m_nip;
    using QuadratureBase<ne, nd, td>::m_nne;
    using QuadratureBase<ne, nd, td>::m_ndim;
    using QuadratureBase<ne, nd, td>::m_tdim;

    xt::xtensor<double, 3> m_x;    ///< nodal positions stored per element [#nelem, #nne, #ndim]
    xt::xtensor<double, 1> m_w;    ///< weight of each integration point [nip]
    xt::xtensor<double, 2> m_xi;   ///< local coordinate of each integration point [#nip, #ndim]
    xt::xtensor<double, 2> m_N;    ///< shape functions [#nip, #nne]
    xt::xtensor<double, 3> m_dNxi; ///< shape function grad. wrt local  coor. [#nip, #nne, #ndim]
    xt::xtensor<double, 4> m_dNx;  ///< shape function grad. wrt global coor. [#nelem, #nip, #nne, #ndim]
    xt::xtensor<double, 2> m_vol;  ///< integration point volume [#nelem, #nip]
};

} // namespace Element
} // namespace GooseFEM

#include "Element.hpp"

#endif

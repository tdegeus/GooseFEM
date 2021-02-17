/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_ELEMENT_H
#define GOOSEFEM_ELEMENT_H

#include "config.h"
#include "Allocate.h"

namespace GooseFEM {
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
Check that all of the matrices stored per element (shape: ``[nelem, nne * ndim, nne * ndim]``)
are diagonal.

\param element Element-vectors ("element")
\return ``true`` if all element matrices are diagonal.
*/
bool isDiagonal(const xt::xtensor<double, 3>& elemmat);

/**
Base quadrature-class.
This class does not have a specific element-type in mind, it is used mostly internally
to derive from such that common methods do not have to be reimplementation.
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
    Number of dimensions for node-vectors.

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
    template <size_t rank = 0>
    void asTensor(const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 2 + rank>& qtensor) const;

    /**
    Convert "qscalar" to "qtensor" of certain rank.

    \param "qscalar".
    \return "qtensor".
    */
    template <size_t rank = 0>
    xt::xtensor<double, 2 + rank> AsTensor(const xt::xtensor<double, 2>& qscalar) const;

    /**
    Convert "qscalar" to "qtensor" of certain rank.

    \param rank Tensor rank.
    \param qscalar A "qscalar".
    \return "qtensor".
    */
    xt::xarray<double> AsTensor(size_t rank, const xt::xtensor<double, 2>& qscalar) const;

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.

    \returns Shape as `std::array`.
    */
    template <size_t rank = 0>
    std::array<size_t, 2 + rank> ShapeQtensor() const;

    /**
    Get the shape of a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \param rank The tensor rank.
    \returns Shape as `std::vector`.
    */
    std::vector<size_t> ShapeQtensor(size_t rank) const;

    /**
    Get the shape of a "qscalar" (a "qtensor" of rank 0)

    \returns Shape as `std::vector`.
    */
    std::vector<size_t> ShapeQscalar() const;

    /**
    Get an allocated `xt::xtensor` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.
    Note: the container is not (zero-)initialised.

    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    template <size_t rank = 0>
    xt::xtensor<double, 2 + rank> AllocateQtensor() const;

    /**
    Get an allocated and initialised `xt::xtensor` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Default: rank = 0, a.k.a. scalar.

    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    template <size_t rank = 0>
    xt::xtensor<double, 2 + rank> AllocateQtensor(double val) const;

    /**
    Get an allocated `xt::xarray` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).
    Note: the container is not (zero-)initialised.

    \param rank The tensor rank.
    \returns `xt::xarray` container of the correct shape.
    */
    xt::xarray<double> AllocateQtensor(size_t rank) const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "qtensor" of a certain rank
    (0 = scalar, 1, vector, 2 = 2nd-order tensor, etc.).

    \param rank The tensor rank.
    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    xt::xarray<double> AllocateQtensor(size_t rank, double val) const;

    /**
    Get an allocated `xt::xtensor` to store a "qscalar" (a "qtensor" of rank 0).
    Note: the container is not (zero-)initialised.

    \returns `xt::xarray` container of the correct shape.
    */
    xt::xtensor<double, 2> AllocateQscalar() const;

    /**
    Get an allocated and initialised `xt::xarray` to store a "qscalar" (a "qtensor" of rank 0).

    \param val The value to which to initialise all items.
    \returns `xt::xtensor` container of the correct shape (and rank).
    */
    xt::xtensor<double, 2> AllocateQscalar(double val) const;

protected:
    void initQuadratureBase(size_t nelem, size_t nip);

protected:
    size_t m_nelem; ///< Number of elements.
    size_t m_nip; ///< Number of integration points per element.
    constexpr static size_t m_nne = ne; ///< Number of nodes per element.
    constexpr static size_t m_ndim = nd; ///< Number of dimensions for nodal vectors.
    constexpr static size_t m_tdim = td; ///< Number of dimensions for integration point tensors.
};

} // namespace Element
} // namespace GooseFEM

#include "Element.hpp"

#endif

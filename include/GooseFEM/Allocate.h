/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_ALLOCATE_H
#define GOOSEFEM_ALLOCATE_H

#include "config.h"

namespace GooseFEM {

/**
"Broadcast" a scalar stored in an array (e.g.: ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).

\param arg A "qscalar".
\param ret Corresponding "qtensor", fully allocate to correct shape.
*/
template <size_t dim, size_t rank>
inline void asTensor(const xt::xtensor<double, dim>& arg, xt::xtensor<double, dim + rank>& ret);

/**
"Broadcast" a scalar stored in an array (e.g.: ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).

\param arg A "qscalar".
\param shape The shape of the added dimensions (e.g.: ``[i, j]``).
\return Corresponding "qtensor".
*/
template <size_t dim, size_t rank>
inline xt::xtensor<double, dim + rank> AsTensor(
    const xt::xtensor<double, dim>& arg,
    const std::array<size_t, rank>& shape);

/**
"Broadcast" a scalar stored in an array (e.g.: ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, n, n]``).

\param arg A "qscalar".
\param n The shape along each of the added dimensions.
\return Corresponding "qtensor".
*/
template <size_t dim, size_t rank>
inline xt::xtensor<double, dim + rank> AsTensor(const xt::xtensor<double, dim>& arg, size_t n);

/**
"Broadcast" a scalar stored in an array (e.g.: ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).

\param rank Rank of the tensor (== number of dimensions to add).
\param arg A "qscalar".
\param shape The shape of the added dimensions (e.g.: ``[i, j]``).
\return Corresponding "qtensor".
*/
template <class T>
inline xt::xarray<double> AsTensor(size_t rank, const T& arg, const std::vector<size_t>& shape);

/**
"Broadcast" a scalar stored in an array (e.g.: ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, n, n]``).

\param rank Rank of the tensor (== number of dimensions to add).
\param arg A "qscalar".
\param n The shape along each of the added dimensions.
\return Corresponding "qtensor".
*/
template <class T>
inline xt::xarray<double> AsTensor(size_t rank, const T& arg, size_t n);

/**
Zero-pad columns to a matrix until is that shape ``[m, 3]``.

\param arg A "nodevec" (``arg.shape(1) <= 3``).
\return Corresponding "nodevec" in 3-d (``ret.shape(1) == 3``)
*/
inline xt::xtensor<double, 2> as3d(const xt::xtensor<double, 2>& arg);

} // namespace GooseFEM

#include "Allocate.hpp"

#endif

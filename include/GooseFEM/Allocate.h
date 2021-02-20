/**
Common allocation methods.

\file Allocate.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ALLOCATE_H
#define GOOSEFEM_ALLOCATE_H

#include "config.h"

namespace GooseFEM {

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).

\tparam dim Number of dimensions of scalar array (rank of the input).
\tparam rank Number of tensor dimensions (number of dimensions to add to the input).
\tparam T Type of the data.
\param arg An array with scalars.
\param ret Corresponding array with tensors.
*/
template <size_t dim, size_t rank, class T>
inline void asTensor(const xt::xtensor<T, dim>& arg, xt::xtensor<T, dim + rank>& ret);

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).

\tparam dim Number of dimensions of scalar array (rank of the input).
\tparam rank Number of tensor dimensions (number of dimensions to add to the input).
\tparam T Type of the data.
\param arg An array with scalars.
\param shape The shape of the added tensor dimensions (e.g.: ``[i, j]``).
\return Corresponding array with tensors.
*/
template <size_t dim, size_t rank, class T>
inline xt::xtensor<T, dim + rank> AsTensor(
    const xt::xtensor<T, dim>& arg,
    const std::array<size_t, rank>& shape);

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, n, n]``).

\tparam dim Number of dimensions of scalar array (rank of the input).
\tparam rank Number of tensor dimensions (number of dimensions to add to the input).
\param arg An array with scalars.
\param n The shape along each of the added dimensions.
\return Corresponding array with tensors.
*/
template <size_t dim, size_t rank, class T>
inline xt::xtensor<T, dim + rank> AsTensor(const xt::xtensor<T, dim>& arg, size_t n);

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).

\tparam T Type of the data.
\param rank Number of tensor dimensions (number of dimensions to add to the input).
\param arg An array with scalars.
\param shape The shape of the added dimensions (e.g.: ``[i, j]``).
\return Corresponding array with tensors.
*/
template <class T>
inline xt::xarray<typename T::value_type> AsTensor(
    size_t rank,
    const T& arg,
    const std::vector<size_t>& shape);

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, n, n]``).

\param rank Number of tensor dimensions (number of dimensions to add to the input).
\param arg An array with scalars.
\param n The shape along each of the added dimensions.
\return Corresponding array with tensors.
*/
template <class T>
inline xt::xarray<typename T::value_type> AsTensor(size_t rank, const T& arg, size_t n);

/**
Zero-pad columns to a matrix until is that shape ``[m, 3]``.

\tparam T Type of the data.
\param arg A "nodevec" (``arg.shape(1) <= 3``).
\return Corresponding "nodevec" in 3-d (``ret.shape(1) == 3``)
*/
template <class T>
inline xt::xtensor<T, 2> as3d(const xt::xtensor<T, 2>& arg);

} // namespace GooseFEM

#include "Allocate.hpp"

#endif

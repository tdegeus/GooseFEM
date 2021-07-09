/**
 *  Common allocation methods.
 *
 *  \file Allocate.h
 *  \copyright Copyright 2017. Tom de Geus. All rights reserved.
 *  \license This project is released under the GNU Public License (GPLv3).
 */

#ifndef GOOSEFEM_ALLOCATE_H
#define GOOSEFEM_ALLOCATE_H

#include "config.h"

namespace GooseFEM {

/**
 *  "Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
 *  tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).
 *
 *  \param arg An array with scalars.
 *  \param ret Corresponding array with tensors.
 */
template <class T, class R>
inline void asTensor(const T& arg, R& ret);

/**
 *  "Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
 *  tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).
 *
 *  \param arg An array with scalars.
 *  \param shape The shape of the added tensor dimensions (e.g.: ``[i, j]``).
 *  \return Corresponding array with tensors.
 */
template <class T, class S>
inline auto AsTensor(const T& arg, const S& shape);

/**
 *  \copydoc AsTensor(const T& arg, const S& shape)
 */
template <class T, class I, size_t L>
inline auto AsTensor(const T& arg, const I (&shape)[L]);

/**
 *  "Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
 *  tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, n, n]``).
 *
 *  \tparam rank Number of tensor dimensions (number of dimensions to add to the input).
 *  \param arg An array with scalars.
 *  \param n The shape along each of the added dimensions.
 *  \return Corresponding array with tensors.
 */
template <size_t rank, class T>
inline auto AsTensor(const T& arg, size_t n);

/**
 *  "Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
 *  tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, n, n]``).
 *
 *  \param rank Number of tensor dimensions (number of dimensions to add to the input).
 *  \param arg An array with scalars.
 *  \param n The shape along each of the added dimensions.
 *  \return Corresponding array with tensors.
 */
template <class T>
inline auto AsTensor(size_t rank, const T& arg, size_t n);

/**
 *  Zero-pad columns to a matrix until is that shape ``[m, 3]``.
 *
 *  \param arg A "nodevec" (``arg.shape(1) <= 3``).
 *  \return Corresponding "nodevec" in 3-d (``ret.shape(1) == 3``)
 */
template <class T>
inline T as3d(const T& arg);

} // namespace GooseFEM

#include "Allocate.hpp"

#endif

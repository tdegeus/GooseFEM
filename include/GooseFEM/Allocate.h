/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_ALLOCATE_H
#define GOOSEFEM_ALLOCATE_H

#include "config.h"

namespace GooseFEM {

// "Broadcast"

template <size_t dim, size_t rank>
inline void asTensor(const xt::xtensor<double, dim>& arg, xt::xtensor<double, dim + rank>& ret);

template <size_t dim, size_t rank>
inline xt::xtensor<double, dim + rank>
AsTensor(const xt::xtensor<double, dim>& arg, const std::array<size_t, rank>& shape);

template <size_t dim, size_t rank>
inline xt::xtensor<double, dim + rank> AsTensor(const xt::xtensor<double, dim>& arg, size_t n);

template <class T>
inline xt::xarray<double> AsTensor(size_t rank, const T& arg, const std::vector<size_t>& shape);

template <class T>
inline xt::xarray<double> AsTensor(size_t rank, const T& arg, size_t n);

} // namespace GooseFEM

#include "Allocate.hpp"

#endif

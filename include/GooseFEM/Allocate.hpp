/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_ALLOCATE_HPP
#define GOOSEFEM_ALLOCATE_HPP

#include "Allocate.h"

namespace GooseFEM {

namespace detail {

template <class T, class S>
inline bool has_shape_begin(const T& t, const S& s)
{
    return s.dimension() >= t.dimension() &&
           std::equal(t.shape().cbegin(), t.shape().cend(), s.shape().begin());
}

} // namespace detail

template <size_t dim, size_t rank>
inline void asTensor(const xt::xtensor<double, dim>& arg, xt::xtensor<double, dim + rank>& ret)
{
    using strides_type = typename xt::xtensor<double, dim>::strides_type::value_type;
    GOOSEFEM_ASSERT(detail::has_shape_begin(arg, ret));
    std::array<strides_type, dim + rank> ret_strides;
    std::copy(arg.strides().begin(), arg.strides().end(), ret_strides.begin());
    std::fill(ret_strides.begin() + dim, ret_strides.end(), 0);
    ret = xt::strided_view(arg, ret.shape(), std::move(ret_strides), 0ul, xt::layout_type::dynamic);
}

template <size_t dim, size_t rank>
inline xt::xtensor<double, dim + rank> AsTensor(
    const xt::xtensor<double, dim>& arg,
    const std::array<size_t, rank>& shape)
{
    std::array<size_t, dim + rank> ret_shape;
    std::copy(arg.shape().begin(), arg.shape().end(), ret_shape.begin());
    std::copy(shape.begin(), shape.end(), ret_shape.begin() + dim);
    xt::xtensor<double, dim + rank> ret = xt::empty<double>(ret_shape);
    GooseFEM::asTensor<dim, rank>(arg, ret);
    return ret;
}

template <size_t dim, size_t rank>
inline xt::xtensor<double, dim + rank> AsTensor(const xt::xtensor<double, dim>& arg, size_t n)
{
    std::array<size_t, dim + rank> ret_shape;
    std::copy(arg.shape().begin(), arg.shape().end(), ret_shape.begin());
    std::fill(ret_shape.begin() + dim, ret_shape.end(), n);
    xt::xtensor<double, dim + rank> ret = xt::empty<double>(ret_shape);
    GooseFEM::asTensor<dim, rank>(arg, ret);
    return ret;
}

template <class T>
inline xt::xarray<double> AsTensor(size_t rank, const T& arg, const std::vector<size_t>& shape)
{
    GOOSEFEM_ASSERT(rank == shape.size());
    size_t dim = arg.dimension();
    std::vector<size_t> ret_shape(dim + rank);
    xt::dynamic_shape<ptrdiff_t> ret_strides(dim + rank);
    std::copy(arg.shape().begin(), arg.shape().end(), ret_shape.begin());
    std::copy(arg.strides().begin(), arg.strides().end(), ret_strides.begin());
    std::copy(shape.begin(), shape.end(), ret_shape.begin() + dim);
    std::fill(ret_strides.begin() + dim, ret_strides.end(), 0);
    xt::xarray<double> ret = xt::empty<double>(ret_shape);
    ret = xt::strided_view(arg, ret.shape(), std::move(ret_strides), 0ul, xt::layout_type::dynamic);
    return ret;
}

template <class T>
inline xt::xarray<double> AsTensor(size_t rank, const T& arg, size_t n)
{
    size_t dim = arg.dimension();
    using strides_type = typename T::strides_type::value_type;
    std::vector<size_t> ret_shape(dim + rank);
    xt::svector<strides_type> ret_strides(dim + rank);
    std::copy(arg.shape().begin(), arg.shape().end(), ret_shape.begin());
    std::copy(arg.strides().begin(), arg.strides().end(), ret_strides.begin());
    std::fill(ret_shape.begin() + dim, ret_shape.end(), n);
    std::fill(ret_strides.begin() + dim, ret_strides.end(), 0);
    xt::xarray<double> ret = xt::empty<double>(ret_shape);
    ret = xt::strided_view(arg, ret.shape(), std::move(ret_strides), 0ul, xt::layout_type::dynamic);
    return ret;
}

inline xt::xtensor<double, 2> as3d(const xt::xtensor<double, 2>& data)
{
    GOOSEFEM_ASSERT(data.shape(1) > 0 && data.shape(1) < 4)

    if (data.shape(1) == 3ul) {
        return data;
    }

    xt::xtensor<double, 2> ret = xt::zeros<double>(std::array<size_t, 2>{data.shape(0), 3ul});

    if (data.shape(1) == 2ul) {
        xt::view(ret, xt::all(), xt::keep(0, 1)) = data;
    }

    if (data.shape(1) == 1ul) {
        xt::view(ret, xt::all(), xt::keep(0)) = data;
    }

    return ret;
}

} // namespace GooseFEM

#endif

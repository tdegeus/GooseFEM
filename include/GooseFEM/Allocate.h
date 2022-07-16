/**
Common allocation methods.

\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ALLOCATE_H
#define GOOSEFEM_ALLOCATE_H

#include "config.h"

namespace GooseFEM {

template <class T, class R>
inline void asTensor(const T& arg, R& ret);

namespace detail {

/**
Check that two shapes partly overlap. If `s` has more dimensions that `t` the excess dimensions
of `s` are ignored and the first `t.size()` dimensions are checked for equality.
*/
template <class T, class S>
inline bool has_shape_begin(const T& t, const S& s)
{
    return s.dimension() >= t.dimension() &&
           std::equal(t.shape().cbegin(), t.shape().cend(), s.shape().begin());
}

/**
Static identification of an std::array
*/
template <class T>
struct is_std_array : std::false_type {
};

template <class T, size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {
};

/**
Helper for std_array_size
*/
template <class T, std::size_t N>
auto std_array_size_impl(const std::array<T, N>&) -> std::integral_constant<std::size_t, N>;

/**
Get the size of an std:array (T::size is not static)
*/
template <class T>
using std_array_size = decltype(std_array_size_impl(std::declval<const T&>()));

/**
Return as std::array.
*/
template <class I, std::size_t L>
std::array<I, L> to_std_array(const I (&shape)[L])
{
    std::array<I, L> r;
    std::copy(&shape[0], &shape[0] + L, r.begin());
    return r;
}

/**
asTensor for array_type::array.
*/
template <class T, class R, typename = void>
struct asTensor_write {
    static void impl(const T& arg, R& ret)
    {
        GOOSEFEM_ASSERT(arg.dimension() <= ret.dimension());
        GOOSEFEM_ASSERT(detail::has_shape_begin(arg, ret));
        using strides_type = typename T::strides_type::value_type;
        std::vector<strides_type> ret_strides(ret.dimension());
        std::copy(arg.strides().begin(), arg.strides().end(), ret_strides.begin());
        std::fill(ret_strides.begin() + arg.dimension(), ret_strides.end(), 0);
        ret = xt::strided_view(
            arg, ret.shape(), std::move(ret_strides), 0ul, xt::layout_type::dynamic);
    }
};

/**
asTensor for xt::tensor.
*/
template <class T, class R>
struct asTensor_write<
    T,
    R,
    typename std::enable_if_t<xt::has_fixed_rank_t<T>::value && xt::has_fixed_rank_t<R>::value>> {
    static void impl(const T& arg, R& ret)
    {
        static_assert(T::rank <= R::rank, "Return must be fixed rank too");
        GOOSEFEM_ASSERT(detail::has_shape_begin(arg, ret));
        using strides_type = typename T::strides_type::value_type;
        std::array<strides_type, R::rank> ret_strides;
        std::copy(arg.strides().begin(), arg.strides().end(), ret_strides.begin());
        std::fill(ret_strides.begin() + T::rank, ret_strides.end(), 0);
        ret = xt::strided_view(
            arg, ret.shape(), std::move(ret_strides), 0ul, xt::layout_type::dynamic);
    }
};

/**
AsTensor for array_type::array.
*/
template <class T, class S, typename = void>
struct asTensor_allocate {
    static auto impl(const T& arg, const S& shape)
    {
        using value_type = typename T::value_type;
        size_t dim = arg.dimension();
        size_t rank = shape.size();
        std::vector<size_t> ret_shape(dim + rank);
        std::copy(arg.shape().begin(), arg.shape().end(), ret_shape.begin());
        std::copy(shape.begin(), shape.end(), ret_shape.begin() + dim);
        xt::xarray<value_type> ret(ret_shape);
        GooseFEM::asTensor(arg, ret);
        return ret;
    }
};

/**
AsTensor for xt::tensor.
*/
template <class T, class S>
struct asTensor_allocate<T, S, typename std::enable_if_t<detail::is_std_array<S>::value>> {
    static auto impl(const T& arg, const S& shape)
    {
        using value_type = typename T::value_type;
        static constexpr size_t dim = T::rank;
        static constexpr size_t rank = std_array_size<S>::value;
        std::array<size_t, dim + rank> ret_shape;
        std::copy(arg.shape().begin(), arg.shape().end(), ret_shape.begin());
        std::copy(shape.begin(), shape.end(), ret_shape.begin() + dim);
        array_type::tensor<value_type, dim + rank> ret = xt::empty<value_type>(ret_shape);
        detail::asTensor_write<std::decay_t<T>, decltype(ret)>::impl(arg, ret);
        return ret;
    }
};

} // namespace detail

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).

\param arg An array with scalars.
\param ret Corresponding array with tensors.
*/
template <class T, class R>
inline void asTensor(const T& arg, R& ret)
{
    detail::asTensor_write<std::decay_t<T>, std::decay_t<R>>::impl(arg, ret);
}

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, i, j]``).

\param arg An array with scalars.
\param shape The shape of the added tensor dimensions (e.g.: ``[i, j]``).
\return Corresponding array with tensors.
*/
template <class T, class S>
inline auto AsTensor(const T& arg, const S& shape)
{
    return detail::asTensor_allocate<std::decay_t<T>, std::decay_t<S>>::impl(arg, shape);
}

/**
\copydoc AsTensor(const T& arg, const S& shape)
*/
template <class T, class I, size_t L>
inline auto AsTensor(const T& arg, const I (&shape)[L])
{
    auto s = detail::to_std_array(shape);
    return detail::asTensor_allocate<std::decay_t<T>, decltype(s)>::impl(arg, s);
}

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, n, n]``).

\tparam rank Number of tensor dimensions (number of dimensions to add to the input).
\param arg An array with scalars.
\param n The shape along each of the added dimensions.
\return Corresponding array with tensors.
*/
template <size_t rank, class T>
inline auto AsTensor(const T& arg, size_t n)
{
    std::array<size_t, rank> shape;
    std::fill(shape.begin(), shape.end(), n);
    return detail::asTensor_allocate<std::decay_t<T>, decltype(shape)>::impl(arg, shape);
}

/**
"Broadcast" a scalar stored in an array (e.g. ``[r, s]``) to the same scalar of all
tensor components of a tensor of certain rank (e.g. for rank 2: ``[r, s, n, n]``).

\param rank Number of tensor dimensions (number of dimensions to add to the input).
\param arg An array with scalars.
\param n The shape along each of the added dimensions.
\return Corresponding array with tensors.
*/
template <class T>
inline auto AsTensor(size_t rank, const T& arg, size_t n)
{
    std::vector<size_t> shape(rank);
    std::fill(shape.begin(), shape.end(), n);
    return detail::asTensor_allocate<std::decay_t<T>, decltype(shape)>::impl(arg, shape);
}

/**
Zero-pad columns to a matrix until is that shape ``[m, 3]``.

\param arg A "nodevec" (``arg.shape(1) <= 3``).
\return Corresponding "nodevec" in 3-d (``ret.shape(1) == 3``)
*/
template <class T>
inline T as3d(const T& arg)
{
    GOOSEFEM_ASSERT(arg.dimension() == 2);
    GOOSEFEM_ASSERT(arg.shape(1) > 0 && arg.shape(1) < 4);

    if (arg.shape(1) == 3ul) {
        return arg;
    }

    T ret = xt::zeros<typename T::value_type>(std::array<size_t, 2>{arg.shape(0), 3ul});

    if (arg.shape(1) == 2ul) {
        xt::view(ret, xt::all(), xt::keep(0, 1)) = arg;
    }

    if (arg.shape(1) == 1ul) {
        xt::view(ret, xt::all(), xt::keep(0)) = arg;
    }

    return ret;
}

} // namespace GooseFEM

#endif

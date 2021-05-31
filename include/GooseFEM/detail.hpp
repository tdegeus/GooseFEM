/**
Support functions.
Not part of public API.

\file detail.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_DETAIL_HPP
#define GOOSEFEM_DETAIL_HPP

namespace GooseFEM {
namespace detail {

    template <size_t nd, typename = void>
    struct tensor
    {
    };

    template <size_t nd>
    struct tensor<nd, typename std::enable_if_t<nd == 2>>
    {
        /**
        Inverse of a 2nd order tensor.

        \param The tensor.
        \param The inverse (overwritten).
        \return Determinant.
        */
        template <class T>
        static double inv(const T& A, T& Ainv)
        {
            double det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

            Ainv(0, 0) = A(1, 1) / det;
            Ainv(0, 1) = -1.0 * A(0, 1) / det;
            Ainv(1, 0) = -1.0 * A(1, 0) / det;
            Ainv(1, 1) = A(0, 0) / det;

            return det;
        }
    };

    template <size_t nd>
    struct tensor<nd, typename std::enable_if_t<nd == 3>>
    {
        /**
        Inverse of a 2nd order tensor (shape: [3, 3]).

        \param The tensor.
        \param The inverse (overwritten).
        \return Determinant.
        */
        template <class T>
        static double inv(const T& A, T& Ainv)
        {
            double det =
                (A(0, 0) * A(1, 1) * A(2, 2) + A(0, 1) * A(1, 2) * A(2, 0) + A(0, 2) * A(1, 0) * A(2, 1)) -
                (A(0, 2) * A(1, 1) * A(2, 0) + A(0, 1) * A(1, 0) * A(2, 2) + A(0, 0) * A(1, 2) * A(2, 1));

            Ainv(0, 0) = (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) / det;
            Ainv(0, 1) = (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2)) / det;
            Ainv(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) / det;

            Ainv(1, 0) = (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)) / det;
            Ainv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) / det;
            Ainv(1, 2) = (A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2)) / det;

            Ainv(2, 0) = (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) / det;
            Ainv(2, 1) = (A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1)) / det;
            Ainv(2, 2) = (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0)) / det;

            return det;
        }
    };

} // namespace detail
} // namespace GooseFEM

#endif

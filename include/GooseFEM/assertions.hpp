/**
Implementation of assertions.h

\file assertions.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ASSERTIONS_HPP
#define GOOSEFEM_ASSERTIONS_HPP

#include "assertions.h"

namespace GooseFEM {

template <class T>
inline bool is_unique(const T& arg)
{
    xt::xtensor<typename T::value_type, 1> tmp = xt::ravel(arg);
    return xt::unique(tmp) == xt::sort(tmp);
}

} // namespace GooseFEM

#endif

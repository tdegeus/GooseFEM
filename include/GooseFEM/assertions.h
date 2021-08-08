/**
 *  Common assertions with verbosity.
 *
 *  \file assertions.h
 *  \copyright Copyright 2017. Tom de Geus. All rights reserved.
 *  \license This project is released under the GNU Public License (GPLv3).
 */

#ifndef GOOSEFEM_ASSERTIONS_H
#define GOOSEFEM_ASSERTIONS_H

#include "config.h"

namespace GooseFEM {

/**
 *  Returns true is a list is unique (has not duplicate items).
 *  \param a Array-like.
 *  \return `true` if `a` is unique.
 */
template <class T>
inline bool is_unique(const T& a);

} // namespace GooseFEM

#include "assertions.hpp"

#endif

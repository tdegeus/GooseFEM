/**
Implementation of Iterate.h

\file Iterate.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ITERATE_HPP
#define GOOSEFEM_ITERATE_HPP

#include "Iterate.h"

namespace GooseFEM {
namespace Iterate {

inline StopList::StopList(size_t n)
{
    m_res.resize(n);
    reset();
}

inline void StopList::reset()
{
    std::fill(m_res.begin(), m_res.end(), std::numeric_limits<double>::infinity());
}

inline void StopList::reset(size_t n)
{
    m_res.resize(n);
    reset();
}

inline void StopList::roll_insert(double res)
{
    std::rotate(m_res.begin(), m_res.begin() + 1, m_res.end());
    m_res.back() = res;
}

inline bool StopList::descending() const
{
    return std::is_sorted(m_res.cbegin(), m_res.cend(), std::greater<double>());
}

inline bool StopList::all_less(double tol) const
{
    return !std::any_of(m_res.cbegin(), m_res.cend(), [=](const auto& i) { return i >= tol; });
}

inline auto StopList::get() const
{
    return m_res;
}

} // namespace Iterate
} // namespace GooseFEM

#endif

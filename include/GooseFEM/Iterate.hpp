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

inline bool StopList::stop_simple(double res, double tol)
{
    std::rotate(m_res.begin(), m_res.begin() + 1, m_res.end());
    m_res.back() = res;
    return !std::any_of(m_res.cbegin(), m_res.cend(), [=](const auto& i) { return i >= tol; });
}

inline bool StopList::stop(double res, double tol)
{
    // move residual one place back and add new residual to the end
    std::rotate(m_res.begin(), m_res.begin() + 1, m_res.end());
    m_res.back() = res;

    // check for convergence: all residuals should be below the tolerance
    for (size_t i = 0; i < m_res.size(); ++i) {
        if (m_res[i] > tol) {
            return false;
        }
    }

    // check for convergence: all residuals should be decreasing
    for (size_t i = 1; i < m_res.size(); ++i) {
        if (m_res[i] > m_res[i - 1]) {
            return false;
        }
    }

    // all checks passed: signal convergence
    return true;
}

inline auto StopList::get() const
{
    return m_res;
}

} // namespace Iterate
} // namespace GooseFEM

#endif

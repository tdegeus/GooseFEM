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
    for (auto& i : m_res) {
        i = std::numeric_limits<double>::infinity();
    }
}

inline void StopList::reset(size_t n)
{
    m_res.resize(n);
    reset();
}

inline bool StopList::stop(double res, double tol)
{
    // move residual one place back (forgetting the first)
    for (size_t i = 1; i < m_res.size(); ++i) {
        m_res[i - 1] = m_res[i];
    }

    // add new residual to the end
    m_res[m_res.size() - 1] = res;

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

} // namespace Iterate
} // namespace GooseFEM

#endif

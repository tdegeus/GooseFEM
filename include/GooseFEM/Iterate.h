/**
Support function for iterations.

\file Iterate.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_ITERATE_H
#define GOOSEFEM_ITERATE_H

#include "config.h"

namespace GooseFEM {

/**
Support function for iterations in end-user programs.
*/
namespace Iterate {

/**
Class to perform a residual check based on the last "n" iterations.
A typical usage is in dynamic simulations where equilibrium is checked based on a force residual.
Fluctuations could however be responsible for this criterion to be triggered too early.
By checking several time-steps such case can be avoided.
*/
class StopList {
public:

    /**
    Constructor.

    \param n Number of consecutive iterations to consider.
    */
    StopList(size_t n = 1);

    /**
    Reset all residuals to infinity.
    */
    void reset();

    /**
    Reset all residuals to infinity, and change the number of residuals to check.

    \param n Number of consecutive iterations to consider.
    */
    void reset(size_t n);

    /**
    Update list of residuals, return `true` if all residuals are below the tolerance.

    \param res Current residual.
    \param tol Tolerance below which all last "n" iterations must lie.
    */
    bool stop(double res, double tol);

private:
    std::vector<double> m_res; ///< List with residuals.
};

} // namespace Iterate
} // namespace GooseFEM

#include "Iterate.hpp"

#endif

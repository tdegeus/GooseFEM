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
     *  Roll the list with the residuals, and add a new residual to the end.
     *  In Python code this function corresponds to::
     *
     *      residuals = residuals[1:] + [new_residual]
     *
     *  I.e. the residual of `n` iterations ago will be forgotten.
     *
     *  \param res New residual to add to the list of residuals.
     */
    void roll_insert(double res);

    /**
     *  Check of the sequence of `n` residuals is in descending order.
     *
     *  \return `true` if the `n` residuals are in descending order.
     */
    bool descending() const;

    /**
     *  Check of the sequence of `n` residuals are all below a tolerance.
     *
     *  \param tol Tolerance.
     *  \return `true` if all `n` residuals are less than the tolerance.
     */
    bool all_less(double tol) const;

    /**
    Update list of residuals, return `true` if all residuals are below the tolerance.

    \param res Current residual.
    \param tol Tolerance below which all last "n" iterations must lie.
    */
    [[deprecated]] bool stop_simple(double res, double tol);

    /**
    Update list of residuals, return `true` if all residuals are below the tolerance,
    and all residuals are in a descending order.

    \param res Current residual.
    \param tol Tolerance below which all last "n" iterations must lie.
    */
    [[deprecated]] bool stop(double res, double tol);

    /**
    Get the historic residuals.
    */
    auto get() const;

private:
    std::vector<double> m_res; ///< List with residuals.
};

} // namespace Iterate
} // namespace GooseFEM

#include "Iterate.hpp"

#endif

/**
Methods to store and apply DOF tyings.

\file TyingsPeriodic.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_TYINGSPERIODIC_H
#define GOOSEFEM_TYINGSPERIODIC_H

#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace GooseFEM {
namespace Tyings {

class Periodic {
public:
    // Constructors
    Periodic() = default;

    Periodic(
        const xt::xtensor<double, 2>& coor,
        const xt::xtensor<size_t, 2>& dofs,
        const xt::xtensor<size_t, 2>& control_dofs,
        const xt::xtensor<size_t, 2>& nodal_tyings); // (independent, dependent)

    Periodic(
        const xt::xtensor<double, 2>& coor,
        const xt::xtensor<size_t, 2>& dofs,
        const xt::xtensor<size_t, 2>& control_dofs,
        const xt::xtensor<size_t, 2>& nodal_tyings, // (independent, dependent)
        const xt::xtensor<size_t, 1>& iip);

    // Dimensions
    size_t nnd() const; // dependent DOFs
    size_t nni() const; // independent DOFs
    size_t nnu() const; // independent, unknown DOFs
    size_t nnp() const; // independent, prescribed DOFs

    // DOF lists
    xt::xtensor<size_t, 2> dofs() const;    // DOFs
    xt::xtensor<size_t, 2> control() const; // control DOFs
    xt::xtensor<size_t, 1> iid() const;     // dependent DOFs
    xt::xtensor<size_t, 1> iii() const;     // independent DOFs
    xt::xtensor<size_t, 1> iiu() const;     // independent, unknown DOFs
    xt::xtensor<size_t, 1> iip() const;     // independent, prescribed DOFs

    // Return the tying matrix
    // u_d = C_di * u_i
    // u_d = [C_du, C_dp]^T * [u_u, u_p] = C_du * u_u + C_dp * u_p
    Eigen::SparseMatrix<double> Cdi() const;
    Eigen::SparseMatrix<double> Cdu() const;
    Eigen::SparseMatrix<double> Cdp() const;

private:
    size_t m_nnu;
    size_t m_nnp;
    size_t m_nni;
    size_t m_nnd;
    size_t m_ndim;
    size_t m_nties; // number of nodal ties
    xt::xtensor<size_t, 2> m_dofs;
    xt::xtensor<size_t, 2> m_control;
    xt::xtensor<size_t, 2> m_tyings; // nodal ties: (independent, dependent)
    xt::xtensor<double, 2> m_coor;
};

class Control {
public:
    // Constructors
    Control() = default;
    Control(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& dofs);

    // Extract new lists
    xt::xtensor<double, 2> coor() const;
    xt::xtensor<size_t, 2> dofs() const;
    xt::xtensor<size_t, 2> controlDofs() const;
    xt::xtensor<size_t, 1> controlNodes() const;

private:
    xt::xtensor<double, 2> m_coor;
    xt::xtensor<size_t, 2> m_dofs;
    xt::xtensor<size_t, 2> m_control_dofs;
    xt::xtensor<size_t, 1> m_control_nodes;
};

} // namespace Tyings
} // namespace GooseFEM

#include "TyingsPeriodic.hpp"

#endif

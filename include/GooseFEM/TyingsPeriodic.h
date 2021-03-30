/**
Tools to store and apply nodal/DOF tyings.

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

/**
Tools to store and apply nodal/DOF tyings.
*/
namespace Tyings {

/**
Nodal tyings per periodic boundary conditions.
The idea is that the displacement of all DOFs of a node are tied to another node
and to the average displacement gradient.
The latter is applied/measured using 'virtual' control nodes.

Consider the DOF list \f$ u \f$ renumbered such that it is split up in
independent and dependent DOFs as follows

\f$ u = \begin{bmatrix} u_i \\ u_d \end{bmatrix}\f$

whereby the independent DOFs are furthermore split up in unknown and prescribed nodes as follows

\f$ u_i = \begin{bmatrix} u_u \\ u_p \end{bmatrix}\f$

such that

\f$ u = \begin{bmatrix} u_u \\ u_p \\ u_d \end{bmatrix}\f$

\todo Document how the DOFs are tied to the control nodes, and what the has to do with the mean.
*/
class Periodic {
public:

    Periodic() = default;

    /**
    Constructor.

    \param coor Nodal coordinates [nnode, ndim].
    \param dofs DOF-numbers per node [nnode, ndim].
    \param control_dofs DOF-numbers per control node [ndim, ndim].
    \param nodal_tyings List of nodal tyings, see nodal_tyings(). [ntyings, 2].
    */
    Periodic(
        const xt::xtensor<double, 2>& coor,
        const xt::xtensor<size_t, 2>& dofs,
        const xt::xtensor<size_t, 2>& control_dofs,
        const xt::xtensor<size_t, 2>& nodal_tyings);

    /**
    Constructor.

    \param coor Nodal coordinates [nnode, ndim].
    \param dofs DOF-numbers per node [nnode, ndim].
    \param control_dofs DOF-numbers per control node [ndim, ndim].
    \param nodal_tyings List of nodal tyings, see nodal_tyings(). [ntyings, 2].
    \param iip List of prescribed DOF-numbers.
    */
    Periodic(
        const xt::xtensor<double, 2>& coor,
        const xt::xtensor<size_t, 2>& dofs,
        const xt::xtensor<size_t, 2>& control_dofs,
        const xt::xtensor<size_t, 2>& nodal_tyings,
        const xt::xtensor<size_t, 1>& iip);

    /**
    Number of dependent DOFs.

    \return unsigned int
    */
    size_t nnd() const;

    /**
    Number of independent DOFs.

    \return unsigned int
    */
    size_t nni() const;

    /**
    Number of independent unknown DOFs.

    \return unsigned int
    */
    size_t nnu() const;

    /**
    Number of independent prescribed DOFs.

    \return unsigned int
    */
    size_t nnp() const;

    /**
    Return the DOF-numbers per node, as used internally (after renumbering).

    \return [nnode, ndim].
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    Return the DOF-numbers for each control node, as used internally (after renumbering).

    \return [ndim, ndim].
    */
    xt::xtensor<size_t, 2> control() const;

    /**
    Return the applied nodal tyings.
    Per tying (row) two node numbers are specified,
    according to the convention (independent, dependent).

    \return [ntyings, 2].
    */
    xt::xtensor<size_t, 2> nodal_tyings() const;

    /**
    Dependent DOFs.

    \return List of DOF numbers.
    */
    xt::xtensor<size_t, 1> iid() const;

    /**
    Independent DOFs.

    \return List of DOF numbers.
    */
    xt::xtensor<size_t, 1> iii() const;

    /**
    Independent unknown DOFs.

    \return List of DOF numbers.
    */
    xt::xtensor<size_t, 1> iiu() const;

    /**
    Independent prescribed DOFs.

    \return List of DOF numbers.
    */
    xt::xtensor<size_t, 1> iip() const;

    /**
    Return tying matrix such as to get the dependent DOFs \f$ u_d \f$ from
    the independent DOFs \f$ u_i \f$ as follows

    \f$ u_d = C_{di} u_i \f$

    Note that this can be further partitioned in

    \f$ u_d = C_{du} u_u + C_{dp} u_p \f$

    See Cdu() and Cdp().

    \return Sparse matrix.
    */
    Eigen::SparseMatrix<double> Cdi() const;

    /**
    Unknown part of the partitioned tying matrix, see Cdi().

    \return Sparse matrix.
    */
    Eigen::SparseMatrix<double> Cdu() const;

    /**
    Prescribed part of the partitioned tying matrix, see Cdi().

    \return Sparse matrix.
    */
    Eigen::SparseMatrix<double> Cdp() const;

private:
    size_t m_nnu; ///< See nnu().
    size_t m_nnp; ///< See nnp().
    size_t m_nni; ///< See nni().
    size_t m_nnd; ///< See nnd().
    size_t m_ndim; ///< Number of dimensions.
    size_t m_nties; ///< Number of nodal ties.
    xt::xtensor<size_t, 2> m_dofs; ///< See dofs().
    xt::xtensor<size_t, 2> m_control; ///< See control().
    xt::xtensor<size_t, 2> m_tyings; ///< See nodal_tyings().
    xt::xtensor<double, 2> m_coor; ///< Nodal coordinates [nnode, ndim].
};

/**
Add control nodes to an existing system.
*/
class Control {
public:

    Control() = default;

    /**
    Constructor.

    \param coor Nodal coordinates [nnode, ndim].
    \param dofs DOF-numbers per node [nnode, ndim].
    */
    Control(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& dofs);

    /**
    Nodal coordinates, for the system with control nodes added to it.

    \param [nnode + ndim, ndim], with nnode the number of nodes of the original system.
    */
    xt::xtensor<double, 2> coor() const;

    /**
    DOF-numbers per node, for the system with control nodes added to it.

    \param [nnode + ndim, ndim], with nnode the number of nodes of the original system.
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    DOF-numbers of each control node.

    \param [ndim, ndim].
    */
    xt::xtensor<size_t, 2> controlDofs() const;

    /**
    Node-numbers of the control nodes.

    \param [ndim].
    */
    xt::xtensor<size_t, 1> controlNodes() const;

private:
    xt::xtensor<double, 2> m_coor; ///< See coor().
    xt::xtensor<size_t, 2> m_dofs; ///< See dofs().
    xt::xtensor<size_t, 2> m_control_dofs;  ///< See controlDofs().
    xt::xtensor<size_t, 1> m_control_nodes; ///< See controlNodes().
};

} // namespace Tyings
} // namespace GooseFEM

#include "TyingsPeriodic.hpp"

#endif

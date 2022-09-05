/**
Tools to store and apply nodal/DOF tyings.

\file TyingsPeriodic.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_TYINGSPERIODIC_H
#define GOOSEFEM_TYINGSPERIODIC_H

#include "Mesh.h"
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

    \tparam C array_type::tensor<double, 2>
    \tparam D array_type::tensor<size_t, 2>
    \tparam S array_type::tensor<size_t, 2>
    \tparam T array_type::tensor<size_t, 2>
    \param coor Nodal coordinates [nnode, ndim].
    \param dofs DOF-numbers per node [nnode, ndim].
    \param control_dofs DOF-numbers per control node [ndim, ndim].
    \param nodal_tyings List of nodal tyings, see nodal_tyings(). [ntyings, 2].
    */
    template <class C, class D, class S, class T>
    Periodic(const C& coor, const D& dofs, const S& control_dofs, const T& nodal_tyings)
        : Periodic(coor, dofs, control_dofs, nodal_tyings, xt::eval(xt::empty<size_t>({0})))
    {
    }

    /**
    Constructor.

    \tparam C array_type::tensor<double, 2>
    \tparam D array_type::tensor<size_t, 2>
    \tparam S array_type::tensor<size_t, 2>
    \tparam T array_type::tensor<size_t, 2>
    \tparam U array_type::tensor<size_t, 1>
    \param coor Nodal coordinates [nnode, ndim].
    \param dofs DOF-numbers per node [nnode, ndim].
    \param control_dofs DOF-numbers per control node [ndim, ndim].
    \param nodal_tyings List of nodal tyings, see nodal_tyings(). [ntyings, 2].
    \param iip List of prescribed DOF-numbers.
    */
    template <class C, class D, class S, class T, class U>
    Periodic(
        const C& coor,
        const D& dofs,
        const S& control_dofs,
        const T& nodal_tyings,
        const U& iip)
    {
        m_tyings = nodal_tyings;
        m_coor = coor;
        m_ndim = m_coor.shape(1);
        m_nties = m_tyings.shape(0);

        GOOSEFEM_ASSERT(xt::has_shape(m_tyings, {m_nties, size_t(2)}));
        GOOSEFEM_ASSERT(xt::has_shape(control_dofs, {m_ndim, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(dofs, m_coor.shape()));
        GOOSEFEM_ASSERT(xt::amax(control_dofs)() <= xt::amax(dofs)());
        GOOSEFEM_ASSERT(xt::amin(control_dofs)() >= xt::amin(dofs)());
        GOOSEFEM_ASSERT(xt::amax(iip)() <= xt::amax(dofs)());
        GOOSEFEM_ASSERT(xt::amin(iip)() >= xt::amin(dofs)());
        GOOSEFEM_ASSERT(xt::amax(nodal_tyings)() < m_coor.shape(0));

        array_type::tensor<size_t, 1> dependent = xt::view(m_tyings, xt::all(), 1);
        array_type::tensor<size_t, 2> dependent_dofs =
            xt::view(dofs, xt::keep(dependent), xt::all());
        U iid = xt::flatten(dependent_dofs);
        U iii = xt::setdiff1d(dofs, iid);
        U iiu = xt::setdiff1d(iii, iip);

        m_nnu = iiu.size();
        m_nnp = iip.size();
        m_nni = iii.size();
        m_nnd = iid.size();

        GooseFEM::Mesh::Reorder reorder({iiu, iip, iid});

        m_dofs = reorder.apply(dofs);
        m_control = reorder.apply(control_dofs);
    }

    /**
    \return Number of dependent DOFs.
    */
    size_t nnd() const
    {
        return m_nnd;
    }

    /**
    \return Number of independent DOFs.
    */
    size_t nni() const
    {
        return m_nni;
    }

    /**
    \return Number of independent unknown DOFs.
    */
    size_t nnu() const
    {
        return m_nnu;
    }

    /**
    \return Number of independent prescribed DOFs.
    */
    size_t nnp() const
    {
        return m_nnp;
    }

    /**
    \return DOF-numbers per node, as used internally (after renumbering), [nnode, ndim].
    */
    const array_type::tensor<size_t, 2>& dofs() const
    {
        return m_dofs;
    }

    /**
    \return DOF-numbers for each control node, as used internally (after renumbering), [ndim, ndim].
    */
    const array_type::tensor<size_t, 2>& control() const
    {
        return m_control;
    }

    /**
    Return the applied nodal tyings.
    Per tying (row) two node numbers are specified,
    according to the convention (independent, dependent).

    \return [ntyings, 2].
    */
    const array_type::tensor<size_t, 2>& nodal_tyings() const
    {
        return m_tyings;
    }

    /**
    Dependent DOFs.

    \return List of DOF numbers.
    */
    array_type::tensor<size_t, 1> iid() const
    {
        return xt::arange<size_t>(m_nni, m_nni + m_nnd);
    }

    /**
    Independent DOFs.

    \return List of DOF numbers.
    */
    array_type::tensor<size_t, 1> iii() const
    {
        return xt::arange<size_t>(m_nni);
    }

    /**
    Independent unknown DOFs.

    \return List of DOF numbers.
    */
    array_type::tensor<size_t, 1> iiu() const
    {
        return xt::arange<size_t>(m_nnu);
    }

    /**
    Independent prescribed DOFs.

    \return List of DOF numbers.
    */
    array_type::tensor<size_t, 1> iip() const
    {
        return xt::arange<size_t>(m_nnp) + m_nnu;
    }

    /**
    Return tying matrix such as to get the dependent DOFs \f$ u_d \f$ from
    the independent DOFs \f$ u_i \f$ as follows

    \f$ u_d = C_{di} u_i \f$

    Note that this can be further partitioned in

    \f$ u_d = C_{du} u_u + C_{dp} u_p \f$

    See Cdu() and Cdp().

    \return Sparse matrix.
    */
    Eigen::SparseMatrix<double> Cdi() const
    {
        std::vector<Eigen::Triplet<double>> data;

        data.reserve(m_nties * m_ndim * (m_ndim + 1));

        for (size_t i = 0; i < m_nties; ++i) {
            for (size_t j = 0; j < m_ndim; ++j) {

                size_t ni = m_tyings(i, 0);
                size_t nd = m_tyings(i, 1);

                data.push_back(Eigen::Triplet<double>(i * m_ndim + j, m_dofs(ni, j), +1.0));

                for (size_t k = 0; k < m_ndim; ++k) {
                    data.push_back(Eigen::Triplet<double>(
                        i * m_ndim + j, m_control(j, k), m_coor(nd, k) - m_coor(ni, k)));
                }
            }
        }

        Eigen::SparseMatrix<double> Cdi;
        Cdi.resize(m_nnd, m_nni);
        Cdi.setFromTriplets(data.begin(), data.end());

        return Cdi;
    }

    /**
    Unknown part of the partitioned tying matrix, see Cdi().

    \return Sparse matrix.
    */
    Eigen::SparseMatrix<double> Cdu() const
    {
        std::vector<Eigen::Triplet<double>> data;

        data.reserve(m_nties * m_ndim * (m_ndim + 1));

        for (size_t i = 0; i < m_nties; ++i) {
            for (size_t j = 0; j < m_ndim; ++j) {

                size_t ni = m_tyings(i, 0);
                size_t nd = m_tyings(i, 1);

                if (m_dofs(ni, j) < m_nnu) {
                    data.push_back(Eigen::Triplet<double>(i * m_ndim + j, m_dofs(ni, j), +1.0));
                }

                for (size_t k = 0; k < m_ndim; ++k) {
                    if (m_control(j, k) < m_nnu) {
                        data.push_back(Eigen::Triplet<double>(
                            i * m_ndim + j, m_control(j, k), m_coor(nd, k) - m_coor(ni, k)));
                    }
                }
            }
        }

        Eigen::SparseMatrix<double> Cdu;
        Cdu.resize(m_nnd, m_nnu);
        Cdu.setFromTriplets(data.begin(), data.end());

        return Cdu;
    }

    /**
    Prescribed part of the partitioned tying matrix, see Cdi().

    \return Sparse matrix.
    */
    Eigen::SparseMatrix<double> Cdp() const
    {
        std::vector<Eigen::Triplet<double>> data;

        data.reserve(m_nties * m_ndim * (m_ndim + 1));

        for (size_t i = 0; i < m_nties; ++i) {
            for (size_t j = 0; j < m_ndim; ++j) {

                size_t ni = m_tyings(i, 0);
                size_t nd = m_tyings(i, 1);

                if (m_dofs(ni, j) >= m_nnu) {
                    data.push_back(
                        Eigen::Triplet<double>(i * m_ndim + j, m_dofs(ni, j) - m_nnu, +1.0));
                }

                for (size_t k = 0; k < m_ndim; ++k) {
                    if (m_control(j, k) >= m_nnu) {
                        data.push_back(Eigen::Triplet<double>(
                            i * m_ndim + j,
                            m_control(j, k) - m_nnu,
                            m_coor(nd, k) - m_coor(ni, k)));
                    }
                }
            }
        }

        Eigen::SparseMatrix<double> Cdp;
        Cdp.resize(m_nnd, m_nnp);
        Cdp.setFromTriplets(data.begin(), data.end());

        return Cdp;
    }

private:
    size_t m_nnu; ///< See nnu().
    size_t m_nnp; ///< See nnp().
    size_t m_nni; ///< See nni().
    size_t m_nnd; ///< See nnd().
    size_t m_ndim; ///< Number of dimensions.
    size_t m_nties; ///< Number of nodal ties.
    array_type::tensor<size_t, 2> m_dofs; ///< See dofs().
    array_type::tensor<size_t, 2> m_control; ///< See control().
    array_type::tensor<size_t, 2> m_tyings; ///< See nodal_tyings().
    array_type::tensor<double, 2> m_coor; ///< Nodal coordinates [nnode, ndim].
};

/**
Add control nodes to an existing system.
*/
class Control {
public:
    Control() = default;

    /**
    Constructor.

    \tparam C array_type::tensor<double, 2>
    \tparam N array_type::tensor<size_t, 2>
    \param coor Nodal coordinates [nnode, ndim].
    \param dofs DOF-numbers per node [nnode, ndim].
    */
    template <class C, class N>
    Control(const C& coor, const N& dofs)
    {
        GOOSEFEM_ASSERT(coor.shape().size() == 2);
        GOOSEFEM_ASSERT(coor.shape() == dofs.shape());

        m_coor = coor;
        m_dofs = dofs;

        size_t nnode = coor.shape(0);
        size_t ndim = coor.shape(1);

        m_control_dofs = xt::arange<size_t>(ndim * ndim).reshape({ndim, ndim});
        m_control_dofs += xt::amax(dofs)() + 1;

        m_control_nodes = nnode + xt::arange<size_t>(ndim);

        m_coor = xt::concatenate(xt::xtuple(coor, xt::zeros<double>({ndim, ndim})));
        m_dofs = xt::concatenate(xt::xtuple(dofs, m_control_dofs));
    }

    /**
    Nodal coordinates, for the system with control nodes added to it.

    \return [nnode + ndim, ndim], with nnode the number of nodes of the original system.
    */
    const array_type::tensor<double, 2>& coor() const
    {
        return m_coor;
    }

    /**
    DOF-numbers per node, for the system with control nodes added to it.

    \return [nnode + ndim, ndim], with nnode the number of nodes of the original system.
    */
    const array_type::tensor<size_t, 2>& dofs() const
    {
        return m_dofs;
    }

    /**
    DOF-numbers of each control node.

    \return [ndim, ndim].
    */
    const array_type::tensor<size_t, 2>& controlDofs() const
    {
        return m_control_dofs;
    }

    /**
    Node-numbers of the control nodes.

    \return [ndim].
    */
    const array_type::tensor<size_t, 1>& controlNodes() const
    {
        return m_control_nodes;
    }

private:
    array_type::tensor<double, 2> m_coor; ///< See coor().
    array_type::tensor<size_t, 2> m_dofs; ///< See dofs().
    array_type::tensor<size_t, 2> m_control_dofs; ///< See controlDofs().
    array_type::tensor<size_t, 1> m_control_nodes; ///< See controlNodes().
};

} // namespace Tyings
} // namespace GooseFEM

#endif

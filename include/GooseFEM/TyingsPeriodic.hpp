/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_TYINGSPERIODIC_HPP
#define GOOSEFEM_TYINGSPERIODIC_HPP

#include "Mesh.h"
#include "TyingsPeriodic.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace GooseFEM {
namespace Tyings {

inline Periodic::Periodic(
    const xt::xtensor<double,2>& coor,
    const xt::xtensor<size_t,2>& dofs,
    const xt::xtensor<size_t,2>& control,
    const xt::xtensor<size_t,2>& nodal_tyings)
    : Periodic(coor, dofs, control, nodal_tyings, xt::empty<size_t>({0}))
{
}

inline Periodic::Periodic(
    const xt::xtensor<double,2>& coor,
    const xt::xtensor<size_t,2>& dofs,
    const xt::xtensor<size_t,2>& control,
    const xt::xtensor<size_t,2>& nodal_tyings,
    const xt::xtensor<size_t,1>& iip)
    : m_tyings(nodal_tyings), m_coor(coor)
{
    m_ndim = m_coor.shape(1);
    m_nties = m_tyings.shape(0);

    xt::xtensor<size_t,1> dependent = xt::view(m_tyings, xt::all(), 1);
    xt::xtensor<size_t,2> dependent_dofs = xt::view(dofs, xt::keep(dependent), xt::all());
    xt::xtensor<size_t,1> iid = xt::flatten(dependent_dofs);
    xt::xtensor<size_t,1> iii = xt::setdiff1d(dofs, iid);
    xt::xtensor<size_t,1> iiu = xt::setdiff1d(iii, iip);

    m_nnu = iiu.size();
    m_nnp = iip.size();
    m_nni = iii.size();
    m_nnd = iid.size();

    GooseFEM::Mesh::Reorder reorder({iiu, iip, iid});

    m_dofs = reorder.apply(dofs);
    m_control = reorder.apply(control);
}

inline xt::xtensor<size_t,2> Periodic::dofs() const
{
    return m_dofs;
}

inline xt::xtensor<size_t,2> Periodic::control() const
{
    return m_control;
}

inline size_t Periodic::nnu() const
{
    return m_nnu;
}

inline size_t Periodic::nnp() const
{
    return m_nnp;
}

inline size_t Periodic::nni() const
{
    return m_nni;
}

inline size_t Periodic::nnd() const
{
    return m_nnd;
}

inline xt::xtensor<size_t,1> Periodic::iiu() const
{
    return xt::arange<size_t>(m_nnu);
}

inline xt::xtensor<size_t,1> Periodic::iip() const
{
    return xt::arange<size_t>(m_nnp) + m_nnu;
}

inline xt::xtensor<size_t,1> Periodic::iii() const
{
    return xt::arange<size_t>(m_nni);
}

inline xt::xtensor<size_t,1> Periodic::iid() const
{
    return xt::arange<size_t>(m_nni, m_nni + m_nnd);
}

inline Eigen::SparseMatrix<double> Periodic::Cdi() const
{
    std::vector<Eigen::Triplet<double>> data;

    data.reserve(m_nties * m_ndim * (m_ndim + 1));

    for (size_t i = 0; i < m_nties; ++i) {
        for (size_t j = 0; j < m_ndim; ++j) {

            size_t ni = m_tyings(i, 0);
            size_t nd = m_tyings(i, 1);

            data.push_back(Eigen::Triplet<double>(i * m_ndim + j, m_dofs(ni, j), +1.));

            for (size_t k = 0; k < m_ndim; ++k) {
                data.push_back(Eigen::Triplet<double>(
                    i * m_ndim + j,
                    m_control(j, k),
                    m_coor(nd, k) - m_coor(ni, k)));
            }
        }
    }

    Eigen::SparseMatrix<double> Cdi;
    Cdi.resize(m_nnd, m_nni);
    Cdi.setFromTriplets(data.begin(), data.end());

    return Cdi;
}

inline Eigen::SparseMatrix<double> Periodic::Cdu() const
{
    std::vector<Eigen::Triplet<double>> data;

    data.reserve(m_nties * m_ndim * (m_ndim + 1));

    for (size_t i = 0; i < m_nties; ++i) {
        for (size_t j = 0; j < m_ndim; ++j) {

            size_t ni = m_tyings(i, 0);
            size_t nd = m_tyings(i, 1);

            if (m_dofs(ni, j) < m_nnu) {
                data.push_back(Eigen::Triplet<double>(i * m_ndim + j, m_dofs(ni, j), +1.));
            }

            for (size_t k = 0; k < m_ndim; ++k) {
                if (m_control(j, k) < m_nnu) {
                    data.push_back(Eigen::Triplet<double>(
                        i * m_ndim + j,
                        m_control(j, k),
                        m_coor(nd, k) - m_coor(ni, k)));
                }
            }
        }
    }

    Eigen::SparseMatrix<double> Cdu;
    Cdu.resize(m_nnd, m_nnu);
    Cdu.setFromTriplets(data.begin(), data.end());

    return Cdu;
}

inline Eigen::SparseMatrix<double> Periodic::Cdp() const
{
    std::vector<Eigen::Triplet<double>> data;

    data.reserve(m_nties * m_ndim * (m_ndim + 1));

    for (size_t i = 0; i < m_nties; ++i) {
        for (size_t j = 0; j < m_ndim; ++j) {

            size_t ni = m_tyings(i, 0);
            size_t nd = m_tyings(i, 1);

            if (m_dofs(ni, j) >= m_nnu) {
                data.push_back(Eigen::Triplet<double>(i * m_ndim + j, m_dofs(ni, j) - m_nnu, +1.));
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

inline Control::Control(const xt::xtensor<double,2>& coor, const xt::xtensor<size_t,2>& dofs)
    : m_coor(coor), m_dofs(dofs)
{
    GOOSEFEM_ASSERT(coor.shape().size() == 2);
    GOOSEFEM_ASSERT(coor.shape() == dofs.shape());

    size_t nnode = coor.shape(0);
    size_t ndim = coor.shape(1);

    m_control_dofs = xt::arange<size_t>(ndim * ndim).reshape({ndim, ndim});
    m_control_dofs += xt::amax(dofs)[0] + 1;

    m_control_nodes = nnode + xt::arange<size_t>(ndim);

    m_coor = xt::concatenate(xt::xtuple(coor, xt::zeros<double>({ndim, ndim})));
    m_dofs = xt::concatenate(xt::xtuple(dofs, m_control_dofs));
}

inline xt::xtensor<double,2> Control::coor() const
{
    return m_coor;
}

inline xt::xtensor<size_t,2> Control::dofs() const
{
    return m_dofs;
}

inline xt::xtensor<size_t,2> Control::controlDofs() const
{
    return m_control_dofs;
}

inline xt::xtensor<size_t,1> Control::controlNodes() const
{
    return m_control_nodes;
}

} // namespace Tyings
} // namespace GooseFEM

#endif

/**
Methods to switch between storage types based on a mesh and DOFs that are partitioned in:
-   unknown DOFs
-   prescribed DOFs
-   dependent DOFs

\file VectorPartitionedTyings.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VECTORPARTITIONEDTYINGS_H
#define GOOSEFEM_VECTORPARTITIONEDTYINGS_H

#include "Vector.h"
#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace GooseFEM {

/**
Class to switch between storage types. In particular:

-   "nodevec": nodal vectors [#nnode, #ndim].
-   "elemvec": nodal vectors stored per element [nelem, #nne, #ndim].
-   "dofval": DOF values [#ndof].
-   "dofval_u": DOF values (Unknown), `== dofval[iiu]`, [#nnu].
-   "dofval_p": DOF values (Prescribed), `== dofval[iiu]`,  [#nnp].
*/
class VectorPartitionedTyings : public Vector {
private:
    array_type::tensor<size_t, 1> m_iiu; ///< See iiu().
    array_type::tensor<size_t, 1> m_iip; ///< See iip().
    array_type::tensor<size_t, 1> m_iid; ///< See iid().
    size_t m_nnu; ///< See nnu().
    size_t m_nnp; ///< See nnp().
    size_t m_nni; ///< See nni().
    size_t m_nnd; ///< See nnd().
    Eigen::SparseMatrix<double> m_Cdu; ///< Tying matrix, see Tyings::Periodic::Cdu().
    Eigen::SparseMatrix<double> m_Cdp; ///< Tying matrix, see Tyings::Periodic::Cdp().
    Eigen::SparseMatrix<double> m_Cdi; ///< Tying matrix, see Tyings::Periodic::Cdi().
    Eigen::SparseMatrix<double> m_Cud; ///< Transpose of "m_Cdu".
    Eigen::SparseMatrix<double> m_Cpd; ///< Transpose of "m_Cdp".
    Eigen::SparseMatrix<double> m_Cid; ///< Transpose of "m_Cdi".

private:
    /**
    Convert to "dofval" (overwrite entries that occur more than once).
    Only the dependent DOFs are retained.

    \param nodevec nodal vectors [#nnode, #ndim].
    \return dofval[iid()] [#nnd].
    */
    template <class T>
    Eigen::VectorXd Eigen_asDofs_d(const T& nodevec) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

        Eigen::VectorXd dofval_d(m_nnd, 1);

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_dofs(m, i) >= m_nni) {
                    dofval_d(m_dofs(m, i) - m_nni) = nodevec(m, i);
                }
            }
        }

        return dofval_d;
    }

public:
    VectorPartitionedTyings() = default;

    /**
    Constructor.

    \tparam E e.g. `array_type::tensor<size_t, 2>`
    \tparam M e.g. `Eigen::SparseMatrix<double>`
    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    \param Cdu See Tyings::Periodic::Cdu().
    \param Cdp See Tyings::Periodic::Cdp().
    \param Cdi See Tyings::Periodic::Cdi().
    */
    template <class E, class M>
    VectorPartitionedTyings(const E& conn, const E& dofs, const M& Cdu, const M& Cdp, const M& Cdi)
        : Vector(conn, dofs), m_Cdu(Cdu), m_Cdp(Cdp), m_Cdi(Cdi)
    {
        GOOSEFEM_ASSERT(Cdu.rows() == Cdp.rows());
        GOOSEFEM_ASSERT(Cdi.rows() == Cdp.rows());

        m_nnu = static_cast<size_t>(m_Cdu.cols());
        m_nnp = static_cast<size_t>(m_Cdp.cols());
        m_nnd = static_cast<size_t>(m_Cdp.rows());
        m_nni = m_nnu + m_nnp;
        m_iiu = xt::arange<size_t>(m_nnu);
        m_iip = xt::arange<size_t>(m_nnu, m_nnu + m_nnp);
        m_iid = xt::arange<size_t>(m_nni, m_nni + m_nnd);
        m_Cud = m_Cdu.transpose();
        m_Cpd = m_Cdp.transpose();
        m_Cid = m_Cdi.transpose();

        GOOSEFEM_ASSERT(static_cast<size_t>(m_Cdi.cols()) == m_nni);
        GOOSEFEM_ASSERT(m_ndof == xt::amax(m_dofs)() + 1);
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
    Dependent DOFs.

    \return List of DOF numbers.
    */
    array_type::tensor<size_t, 1> iid() const
    {
        return m_iid;
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
        return m_iiu;
    }

    /**
    Independent prescribed DOFs.

    \return List of DOF numbers.
    */
    array_type::tensor<size_t, 1> iip() const
    {
        return m_iip;
    }

    /**
    Copy (part of) "dofval" to another "dofval": dofval_dest[iip()] = dofval_src[iip()].

    \param dofval_src DOF values, iip() updated, [#ndof].
    \param dofval_dest DOF values, iip() updated, [#ndof].
    */
    template <class T>
    void copy_p(const T& dofval_src, T& dofval_dest) const
    {
        GOOSEFEM_ASSERT(dofval_src.dimension() == 1);
        GOOSEFEM_ASSERT(dofval_dest.dimension() == 1);
        GOOSEFEM_ASSERT(dofval_src.size() == m_ndof || dofval_src.size() == m_nni);
        GOOSEFEM_ASSERT(dofval_dest.size() == m_ndof || dofval_dest.size() == m_nni);

#pragma omp parallel for
        for (size_t i = m_nnu; i < m_nni; ++i) {
            dofval_dest(i) = dofval_src(i);
        }
    }

    /**
    Convert to "dofval" (overwrite entries that occur more than once).
    Only the independent DOFs are retained.

    \param nodevec nodal vectors [#nnode, #ndim].
    \return dofval[iii()] [#nni].
    */
    template <class T>
    array_type::tensor<double, 1> AsDofs_i(const T& nodevec) const
    {
        array_type::tensor<double, 1> dofval = xt::empty<double>({m_nni});
        this->asDofs_i(nodevec, dofval);
        return dofval;
    }

    /**
    Same as InterpQuad_vector(), but writing to preallocated return.

    \param nodevec [#nnode, #ndim].
    \param dofval_i [#nni].
    \param apply_tyings If `true` the dependent DOFs are eliminated.
    */
    template <class T, class R>
    void asDofs_i(const T& nodevec, R& dofval_i, bool apply_tyings = true) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(dofval_i.size() == m_nni);

        dofval_i.fill(0.0);

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_dofs(m, i) < m_nni) {
                    dofval_i(m_dofs(m, i)) = nodevec(m, i);
                }
            }
        }

        if (!apply_tyings) {
            return;
        }

        Eigen::VectorXd Dofval_d = this->Eigen_asDofs_d(nodevec);
        Eigen::VectorXd Dofval_i = m_Cid * Dofval_d;

#pragma omp parallel for
        for (size_t i = 0; i < m_nni; ++i) {
            dofval_i(i) += Dofval_i(i);
        }
    }
};

} // namespace GooseFEM

#endif

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

#include "config.h"
#include "Vector.h"

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
public:

    VectorPartitionedTyings() = default;

    /**
    Constructor.

    \tparam E e.g. `xt::xtensor<size_t, 2>`
    \tparam M e.g. `Eigen::SparseMatrix<double>`
    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    \param Cdu See Tyings::Periodic::Cdu().
    \param Cdp See Tyings::Periodic::Cdp().
    \param Cdi See Tyings::Periodic::Cdi().
    */
    template <class E, class M>
    VectorPartitionedTyings(const E& conn, const E& dofs, const M& Cdu, const M& Cdp, const M& Cdi);

    /**
    \return Number of dependent DOFs.
    */
    size_t nnd() const;

    /**
    \return Number of independent DOFs.
    */
    size_t nni() const;

    /**
    \return Number of independent unknown DOFs.
    */
    size_t nnu() const;

    /**
    \return Number of independent prescribed DOFs.
    */
    size_t nnp() const;

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
    Copy (part of) "dofval" to another "dofval": dofval_dest[iip()] = dofval_src[iip()].

    \param dofval_src DOF values, iip() updated, [#ndof].
    \param dofval_dest DOF values, iip() updated, [#ndof].
    */
    template <class T>
    void copy_p(const T& dofval_src, T& dofval_dest) const;

    /**
    Convert to "dofval" (overwrite entries that occur more than once).
    Only the independent DOFs are retained.

    \param nodevec nodal vectors [#nnode, #ndim].
    \return dofval[iii()] [#nni].
    */
    template <class T>
    xt::xtensor<double, 1> AsDofs_i(const T& nodevec) const;

    /**
    Same as InterpQuad_vector(), but writing to preallocated return.

    \param nodevec [#nnode, #ndim].
    \param dofval_i [#nni].
    \param apply_tyings If `true` the dependent DOFs are eliminated.
    */
    template <class T, class R>
    void asDofs_i(const T& nodevec, R& dofval_i, bool apply_tyings = true) const;

private:
    xt::xtensor<size_t, 1> m_iiu; ///< See iiu().
    xt::xtensor<size_t, 1> m_iip; ///< See iip().
    xt::xtensor<size_t, 1> m_iid; ///< See iid().
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

    /**
    Convert to "dofval" (overwrite entries that occur more than once).
    Only the dependent DOFs are retained.

    \param nodevec nodal vectors [#nnode, #ndim].
    \return dofval[iid()] [#nnd].
    */
    template <class T>
    Eigen::VectorXd Eigen_asDofs_d(const T& nodevec) const;
};

} // namespace GooseFEM

#include "VectorPartitionedTyings.hpp"

#endif

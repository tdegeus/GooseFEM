/**
\file
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIX_H
#define GOOSEFEM_MATRIX_H

#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace GooseFEM {

// forward declaration
template <class>
class MatrixSolver;

/**
CRTP base class for a solver class.
*/
template <class D>
class MatrixSolverBase {
public:
    /**
    Underlying type.
    */
    using derived_type = D;

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

public:
    /**
    Solve \f$ x = A^{-1} b \f$.

    \param A GooseFEM (sparse) matrix, see e.g. GooseFEM::Matrix().
    \param b nodevec [nelem, ndim].
    \return x nodevec [nelem, ndim].
    */
    template <class M>
    array_type::tensor<double, 2> Solve(M& A, const array_type::tensor<double, 2>& b)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b, A.shape_nodevec()));
        array_type::tensor<double, 2> x = xt::empty_like(b);
        derived_cast().solve_nodevec_impl(A, b, x);
        return x;
    }

    /**
    Solve \f$ x = A^{-1} b \f$.

    \param A GooseFEM (sparse) matrix, see e.g. GooseFEM::Matrix().
    \param b dofval [ndof].
    \return x dofval [ndof].
    */
    template <class M>
    array_type::tensor<double, 1> Solve(M& A, const array_type::tensor<double, 1>& b)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b, A.shape_dofval()));
        array_type::tensor<double, 1> x = xt::empty_like(b);
        derived_cast().solve_dofval_impl(A, b, x);
        return x;
    }

    /**
    Solve \f$ x = A^{-1} b \f$.

    \param A GooseFEM (sparse) matrix, see e.g. GooseFEM::Matrix().
    \param b nodevec [nelem, ndim].
    \param x (overwritten) nodevec [nelem, ndim].
    */
    template <class M>
    void solve(M& A, const array_type::tensor<double, 2>& b, array_type::tensor<double, 2>& x)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b, A.shape_nodevec()));
        GOOSEFEM_ASSERT(xt::has_shape(x, A.shape_nodevec()));
        derived_cast().solve_nodevec_impl(A, b, x);
    }

    /**
    Solve \f$ x = A^{-1} b \f$.

    \param A GooseFEM (sparse) matrix, see e.g. GooseFEM::Matrix().
    \param b dofval [ndof].
    \param x (overwritten) dofval [ndof].
    */
    template <class M>
    void solve(M& A, const array_type::tensor<double, 1>& b, array_type::tensor<double, 1>& x)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b, A.shape_dofval()));
        GOOSEFEM_ASSERT(xt::has_shape(x, A.shape_dofval()));
        derived_cast().solve_dofval_impl(A, b, x);
    }
};

/**
CRTP base class for a extra functions for a partitioned solver class.
*/
template <class D>
class MatrixSolverPartitionedBase {
public:
    /**
    Underlying type.
    */
    using derived_type = D;

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

public:
    /**
    Solve \f$ x = A^{-1} b \f$.

    \param A GooseFEM (sparse) matrix, see e.g. GooseFEM::MatrixPartitioned().
    \param b_u unknown dofval [nnu].
    \param x_p prescribed dofval [nnp]
    \return x_u unknown dofval [nnu].
    */
    template <class M>
    array_type::tensor<double, 1> Solve_u(
        M& A,
        const array_type::tensor<double, 1>& b_u,
        const array_type::tensor<double, 1>& x_p)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b_u, {A.nnu()}));
        GOOSEFEM_ASSERT(xt::has_shape(x_p, {A.nnp()}));
        array_type::tensor<double, 1> x_u = xt::empty_like(b_u);
        derived_cast().solve_u_impl(A, b_u, x_p, x_u);
        return x_u;
    }

    /**
    Same as
    Solve \f$ x = A^{-1} b \f$.

    \param A GooseFEM (sparse) matrix, see e.g. GooseFEM::MatrixPartitioned().
    \param b_u unknown dofval [nnu].
    \param x_p prescribed dofval [nnp]
    \param x_u (overwritten) unknown dofval [nnu].
    */
    template <class M>
    void solve_u(
        M& A,
        const array_type::tensor<double, 1>& b_u,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& x_u)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b_u, {A.nnu()}));
        GOOSEFEM_ASSERT(xt::has_shape(x_p, {A.nnp()}));
        GOOSEFEM_ASSERT(xt::has_shape(x_u, {A.nnu()}));
        derived_cast().solve_u_impl(A, b_u, x_p, x_u);
    }
};

/**
CRTP base class for a matrix.
*/
template <class D>
class MatrixBase {
protected:
    array_type::tensor<size_t, 2> m_conn; ///< Connectivity [#nelem, #nne].
    array_type::tensor<size_t, 2> m_dofs; ///< DOF-numbers per node [#nnode, #ndim].

    size_t m_nelem; ///< See nelem().
    size_t m_nne; ///< See nne().
    size_t m_nnode; ///< See nnode().
    size_t m_ndim; ///< See ndim().
    size_t m_ndof; ///< See ndof().

    bool m_changed = true; ///< Signal changes to data.

public:
    /**
    Underlying type.
    */
    using derived_type = D;

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

public:
    /**
    Number of elements.
    \return Unsigned integer.
    */
    size_t nelem() const
    {
        return derived_cast().m_nelem;
    }

    /**
    Number of nodes per element.
    \return Unsigned integer.
    */
    size_t nne() const
    {
        return derived_cast().m_nne;
    }

    /**
    Number of nodes.
    \return Unsigned integer.
    */
    size_t nnode() const
    {
        return derived_cast().m_nnode;
    }

    /**
    Number of dimensions.
    \return Unsigned integer.
    */
    size_t ndim() const
    {
        return derived_cast().m_ndim;
    }

    /**
    Number of DOFs.
    \return Unsigned integer.
    */
    size_t ndof() const
    {
        return derived_cast().m_ndof;
    }

    /**
    DOFs per node.
    \return [#nnode, #ndim].
    */
    const array_type::tensor<size_t, 2>& dofs() const
    {
        return derived_cast().m_dofs;
    }

    /**
    Connectivity.
    \return [#nelem, #nne].
    */
    const array_type::tensor<size_t, 2>& conn() const
    {
        return derived_cast().m_conn;
    }

    /**
    Shape of "dofval".
    \return [#ndof].
    */
    std::array<size_t, 1> shape_dofval() const
    {
        return std::array<size_t, 1>{derived_cast().m_ndof};
    }

    /**
    Shape of "nodevec".
    \return [#nnode, #ndim].
    */
    std::array<size_t, 2> shape_nodevec() const
    {
        return std::array<size_t, 2>{derived_cast().m_nnode, derived_cast().m_ndim};
    }

    /**
    Shape of "elemmat".
    \return [#nelem, #nne * #ndim, #nne * #ndim].
    */
    std::array<size_t, 3> shape_elemmat() const
    {
        return std::array<size_t, 3>{
            derived_cast().m_nelem,
            derived_cast().m_nne * derived_cast().m_ndim,
            derived_cast().m_nne * derived_cast().m_ndim};
    }

    /**
    Assemble from "elemmat".
    \param elemmat [#nelem, #nne * #ndim, #nne * #ndim].
    */
    template <class T>
    void assemble(const T& elemmat)
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemmat, this->shape_elemmat()));
        derived_cast().assemble_impl(elemmat);
    }

    /**
    Copy as dense matrix.
    \return [#ndof, #ndof].
    */
    array_type::tensor<double, 2> Todense() const
    {
        size_t ndof = derived_cast().m_ndof;
        array_type::tensor<double, 2> ret = xt::empty<double>({ndof, ndof});
        derived_cast().todense_impl(ret);
        return ret;
    }

    /**
    Copy to dense matrix.
    \param ret overwritten [#ndof, #ndof].
    */
    template <class T>
    void todense(T& ret) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(ret, {derived_cast().m_ndof, derived_cast().m_ndof}));
        derived_cast().todense_impl(ret);
    }

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x nodevec [#nelem, #ndim].
    \return b nodevec [#nelem, #ndim].
    */
    array_type::tensor<double, 2> Dot(const array_type::tensor<double, 2>& x) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, this->shape_nodevec()));
        array_type::tensor<double, 2> b = xt::empty_like(x);
        derived_cast().dot_nodevec_impl(x, b);
        return b;
    }

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x dofval [#ndof].
    \return b dofval [#ndof].
    */
    array_type::tensor<double, 1> Dot(const array_type::tensor<double, 1>& x) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, this->shape_dofval()));
        array_type::tensor<double, 1> b = xt::empty_like(x);
        derived_cast().dot_dofval_impl(x, b);
        return b;
    }

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x nodevec [#nelem, #ndim].
    \param b (overwritten) nodevec [#nelem, #ndim].
    */
    void dot(const array_type::tensor<double, 2>& x, array_type::tensor<double, 2>& b) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, this->shape_nodevec()));
        GOOSEFEM_ASSERT(xt::has_shape(b, this->shape_nodevec()));
        derived_cast().dot_nodevec_impl(x, b);
    }

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x dofval [#ndof].
    \param b (overwritten) dofval [#ndof].
    */
    void dot(const array_type::tensor<double, 1>& x, array_type::tensor<double, 1>& b) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, this->shape_dofval()));
        GOOSEFEM_ASSERT(xt::has_shape(b, this->shape_dofval()));
        derived_cast().dot_dofval_impl(x, b);
    }
};

/**
CRTP base class for a partitioned matrix.
*/
template <class D>
class MatrixPartitionedBase : public MatrixBase<D> {
protected:
    array_type::tensor<size_t, 1> m_iiu; ///< See iiu()
    array_type::tensor<size_t, 1> m_iip; ///< See iip()

    size_t m_nnu; ///< See #nnu
    size_t m_nnp; ///< See #nnp

public:
    /**
    Underlying type.
    */
    using derived_type = D;

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

public:
    /**
    Number of unknown DOFs.
    \return Unsigned integer.
    */
    size_t nnu() const
    {
        return derived_cast().m_nnu;
    }

    /**
    Number of prescribed DOFs.
    \return Unsigned integer.
    */
    size_t nnp() const
    {
        return derived_cast().m_nnp;
    }

    /**
    Unknown DOFs.
    \return [#nnu].
    */
    const array_type::tensor<size_t, 1>& iiu() const
    {
        return derived_cast().m_iiu;
    }

    /**
    Prescribed DOFs.
    \return [#nnp].
    */
    const array_type::tensor<size_t, 1>& iip() const
    {
        return derived_cast().m_iip;
    }

    /**
    Right-hand-size for corresponding to the prescribed DOFs:

    \f$ b_p = A_{pu} * x_u + A_{pp} * x_p \f$

    and assemble them to the appropriate places in "nodevec".

    \param x "nodevec" [#nnode, #ndim].
    \param b "nodevec" [#nnode, #ndim].
    \return Copy of `b` with \f$ b_p \f$ overwritten.
    */
    array_type::tensor<double, 2>
    Reaction(const array_type::tensor<double, 2>& x, const array_type::tensor<double, 2>& b) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, this->shape_nodevec()));
        GOOSEFEM_ASSERT(xt::has_shape(b, this->shape_nodevec()));
        array_type::tensor<double, 2> ret = b;
        derived_cast().reaction_nodevec_impl(x, ret);
        return ret;
    }

    /**
    Same as Reaction(const array_type::tensor<double, 2>&, const array_type::tensor<double, 2>&),
    but of "dofval" input and output.

    \param x "dofval" [#ndof].
    \param b "dofval" [#ndof].
    \return Copy of `b` with \f$ b_p \f$ overwritten.
    */
    array_type::tensor<double, 1>
    Reaction(const array_type::tensor<double, 1>& x, const array_type::tensor<double, 1>& b) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, this->shape_dofval()));
        GOOSEFEM_ASSERT(xt::has_shape(b, this->shape_dofval()));
        array_type::tensor<double, 1> ret = b;
        derived_cast().reaction_dofval_impl(x, ret);
        return ret;
    }

    /**
    Same as Reaction(const array_type::tensor<double, 2>&, const array_type::tensor<double, 2>&),
    but inserting in-place.

    \param x "nodevec" [#nnode, #ndim].
    \param b "nodevec" [#nnode, #ndim], \f$ b_p \f$ overwritten.
    */
    void reaction(const array_type::tensor<double, 2>& x, array_type::tensor<double, 2>& b) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, this->shape_nodevec()));
        GOOSEFEM_ASSERT(xt::has_shape(b, this->shape_nodevec()));
        derived_cast().reaction_nodevec_impl(x, b);
    }

    /**
    Same as Reaction(const array_type::tensor<double, 1>&, const array_type::tensor<double, 1>&),
    but inserting in-place.

    \param x "dofval" [#ndof].
    \param b "dofval" [#ndof], \f$ b_p \f$ overwritten.
    */
    void reaction(const array_type::tensor<double, 1>& x, array_type::tensor<double, 1>& b) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, this->shape_dofval()));
        GOOSEFEM_ASSERT(xt::has_shape(b, this->shape_dofval()));
        derived_cast().reaction_dofval_impl(x, b);
    }

    /**
    Same as Reaction(const array_type::tensor<double, 1>&, const array_type::tensor<double, 1>&),
    but with partitioned input and output.

    \param x_u unknown "dofval" [#nnu].
    \param x_p prescribed "dofval" [#nnp].
    \return b_p prescribed "dofval" [#nnp].
    */
    array_type::tensor<double, 1> Reaction_p(
        const array_type::tensor<double, 1>& x_u,
        const array_type::tensor<double, 1>& x_p) const
    {
        array_type::tensor<double, 1> b_p = xt::empty<double>({m_nnp});
        derived_cast().reaction_p_impl(x_u, x_p, b_p);
        return b_p;
    }

    /**
    Same as Reaction_p(const array_type::tensor<double, 1>&, const array_type::tensor<double, 1>&),
    but writing to preallocated output.

    \param x_u unknown "dofval" [#nnu].
    \param x_p prescribed "dofval" [#nnp].
    \param b_p (overwritten) prescribed "dofval" [#nnp].
    */
    void reaction_p(
        const array_type::tensor<double, 1>& x_u,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& b_p) const
    {
        derived_cast().reaction_p_impl(x_u, x_p, b_p);
    }
};

/**
CRTP base class for a partitioned matrix with tying.
*/
template <class D>
class MatrixPartitionedTyingsBase : public MatrixPartitionedBase<D> {
protected:
    array_type::tensor<size_t, 1> m_iii; ///< See iii()
    array_type::tensor<size_t, 1> m_iid; ///< See iid()

    size_t m_nni; ///< See #nni
    size_t m_nnd; ///< See #nnd

public:
    /**
    Underlying type.
    */
    using derived_type = D;

private:
    auto derived_cast() -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }

    auto derived_cast() const -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

public:
    /**
    Number of independent DOFs.
    \return Unsigned integer.
    */
    size_t nni() const
    {
        return derived_cast().m_nni;
    }

    /**
    Number of dependent DOFs.
    \return Unsigned integer.
    */
    size_t nnd() const
    {
        return derived_cast().m_nnd;
    }

    /**
    Independent DOFs.
    \return [#nnu].
    */
    const array_type::tensor<size_t, 1>& iii() const
    {
        return derived_cast().m_iii;
    }

    /**
    Dependent DOFs.
    \return [#nnp].
    */
    const array_type::tensor<size_t, 1>& iid() const
    {
        return derived_cast().m_iid;
    }
};

/**
Sparse matrix.

See GooseFEM::Vector() for bookkeeping definitions.
*/
class Matrix : public MatrixBase<Matrix> {
private:
    friend MatrixBase<Matrix>;

private:
    bool m_changed = true; ///< Signal changes to data.

    Eigen::SparseMatrix<double> m_A; ///< The matrix.

    std::vector<Eigen::Triplet<double>> m_T; ///< Matrix entries.

    /**
    Class to solve the system (allowing single factorisation for multiple right-hand-sides).
    */
    template <class>
    friend class MatrixSolver;

public:
    Matrix() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    */
    Matrix(const array_type::tensor<size_t, 2>& conn, const array_type::tensor<size_t, 2>& dofs)
    {
        m_conn = conn;
        m_dofs = dofs;
        m_nelem = m_conn.shape(0);
        m_nne = m_conn.shape(1);
        m_nnode = m_dofs.shape(0);
        m_ndim = m_dofs.shape(1);
        m_ndof = xt::amax(m_dofs)() + 1;
        m_T.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
        m_A.resize(m_ndof, m_ndof);

        GOOSEFEM_ASSERT(xt::amax(m_conn)() + 1 <= m_nnode);
        GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
    }

private:
    template <class T>
    void assemble_impl(const T& elemmat)
    {
        m_T.clear();

        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    for (size_t n = 0; n < m_nne; ++n) {
                        for (size_t j = 0; j < m_ndim; ++j) {
                            m_T.push_back(Eigen::Triplet<double>(
                                m_dofs(m_conn(e, m), i),
                                m_dofs(m_conn(e, n), j),
                                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
                        }
                    }
                }
            }
        }

        m_A.setFromTriplets(m_T.begin(), m_T.end());
        m_changed = true;
    }

public:
    /**
    Overwrite matrix.

    \param rows Row numbers [m].
    \param cols Column numbers [n].
    \param matrix Data entries `matrix(i, j)` for `rows(i), cols(j)` [m, n].
    */
    void
    set(const array_type::tensor<size_t, 1>& rows,
        const array_type::tensor<size_t, 1>& cols,
        const array_type::tensor<double, 2>& matrix)
    {
        GOOSEFEM_ASSERT(rows.size() == matrix.shape(0));
        GOOSEFEM_ASSERT(cols.size() == matrix.shape(1));
        GOOSEFEM_ASSERT(xt::amax(cols)() < m_ndof);
        GOOSEFEM_ASSERT(xt::amax(rows)() < m_ndof);

        std::vector<Eigen::Triplet<double>> T;

        for (size_t i = 0; i < rows.size(); ++i) {
            for (size_t j = 0; j < cols.size(); ++j) {
                T.push_back(Eigen::Triplet<double>(rows(i), cols(j), matrix(i, j)));
            }
        }

        m_A.setFromTriplets(T.begin(), T.end());
        m_changed = true;
    }

    /**
    Add matrix.

    \param rows Row numbers [m].
    \param cols Column numbers [n].
    \param matrix Data entries `matrix(i, j)` for `rows(i), cols(j)` [m, n].
    */
    void
    add(const array_type::tensor<size_t, 1>& rows,
        const array_type::tensor<size_t, 1>& cols,
        const array_type::tensor<double, 2>& matrix)
    {
        GOOSEFEM_ASSERT(rows.size() == matrix.shape(0));
        GOOSEFEM_ASSERT(cols.size() == matrix.shape(1));
        GOOSEFEM_ASSERT(xt::amax(cols)() < m_ndof);
        GOOSEFEM_ASSERT(xt::amax(rows)() < m_ndof);

        std::vector<Eigen::Triplet<double>> T;

        Eigen::SparseMatrix<double> A(m_ndof, m_ndof);

        for (size_t i = 0; i < rows.size(); ++i) {
            for (size_t j = 0; j < cols.size(); ++j) {
                T.push_back(Eigen::Triplet<double>(rows(i), cols(j), matrix(i, j)));
            }
        }

        A.setFromTriplets(T.begin(), T.end());
        m_A += A;
        m_changed = true;
    }

private:
    template <class T>
    void todense_impl(T& ret) const
    {
        ret.fill(0.0);

        for (int k = 0; k < m_A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_A, k); it; ++it) {
                ret(it.row(), it.col()) = it.value();
            }
        }
    }

    template <class T>
    void dot_nodevec_impl(const T& x, T& b) const
    {
        this->Eigen_asNode_dofval_nodevec(m_A * this->Eigen_AsDofs_nodevec(x), b);
    }

    template <class T>
    void dot_dofval_impl(const T& x, T& b) const
    {
        Eigen::Map<Eigen::VectorXd>(b.data(), b.size()).noalias() =
            m_A * Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
    }

private:
    /**
    Convert "nodevec" to "dofval" (overwrite entries that occur more than once).

    \param nodevec input [#nnode, #ndim]
    \return dofval output [#ndof]
    */
    template <class T>
    Eigen::VectorXd Eigen_AsDofs_nodevec(const T& nodevec) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

        Eigen::VectorXd dofval = Eigen::VectorXd::Zero(m_ndof, 1);

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                dofval(m_dofs(m, i)) = nodevec(m, i);
            }
        }

        return dofval;
    }

    /**
    Convert "dofval" to "nodevec" (overwrite entries that occur more than once).

    \param dofval input [#ndof]
    \param nodevec output [#nnode, #ndim]
    */
    template <class T>
    void Eigen_asNode_dofval_nodevec(const Eigen::VectorXd& dofval, T& nodevec) const
    {
        GOOSEFEM_ASSERT(static_cast<size_t>(dofval.size()) == m_ndof);
        GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                nodevec(m, i) = dofval(m_dofs(m, i));
            }
        }
    }
};

/**
Solve \f$ x = A^{-1} b \f$, for `A` of the GooseFEM::Matrix() class.
You can solve for multiple right-hand-sides using one factorisation.
*/
template <class Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
class MatrixSolver : public MatrixSolverBase<MatrixSolver<Solver>> {
private:
    friend MatrixSolverBase<MatrixSolver<Solver>>;

public:
    MatrixSolver() = default;

private:
    template <class T>
    void solve_nodevec_impl(Matrix& A, const T& b, T& x)
    {
        this->factorize(A);
        Eigen::VectorXd X = m_solver.solve(A.Eigen_AsDofs_nodevec(b));
        A.Eigen_asNode_dofval_nodevec(X, x);
    }

    template <class T>
    void solve_dofval_impl(Matrix& A, const T& b, T& x)
    {
        this->factorize(A);
        Eigen::Map<Eigen::VectorXd>(x.data(), x.size()).noalias() =
            m_solver.solve(Eigen::Map<const Eigen::VectorXd>(b.data(), A.m_ndof));
    }

private:
    Solver m_solver; ///< Solver.
    bool m_factor = true; ///< Signal to force factorization.

    /**
    Compute inverse (evaluated by "solve").
    */
    void factorize(Matrix& A)
    {
        if (!A.m_changed && !m_factor) {
            return;
        }
        m_solver.compute(A.m_A);
        m_factor = false;
        A.m_changed = false;
    }
};

} // namespace GooseFEM

#endif

/**
Diagonal matrix that is partitioned in:
-   unknown DOFs
-   prescribed DOFs

\file MatrixDiagonalPartitioned.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXDIAGONALPARTITIONED_H
#define GOOSEFEM_MATRIXDIAGONALPARTITIONED_H

#include "MatrixDiagonal.h"
#include "Mesh.h"
#include "config.h"

namespace GooseFEM {

/**
Diagonal and partitioned matrix.

See Vector() for bookkeeping definitions.
*/
class MatrixDiagonalPartitioned : public MatrixDiagonal {
public:
    MatrixDiagonalPartitioned() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    \param iip prescribed DOFs [#nnp].
    */
    MatrixDiagonalPartitioned(
        const array_type::tensor<size_t, 2>& conn,
        const array_type::tensor<size_t, 2>& dofs,
        const array_type::tensor<size_t, 1>& iip)
    {
        m_conn = conn;
        m_dofs = dofs;
        m_nelem = m_conn.shape(0);
        m_nne = m_conn.shape(1);
        m_nnode = m_dofs.shape(0);
        m_ndim = m_dofs.shape(1);
        m_ndof = xt::amax(m_dofs)() + 1;

        GOOSEFEM_ASSERT(is_unique(iip));
        GOOSEFEM_ASSERT(xt::amax(conn)() + 1 <= m_nnode);
        GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
        GOOSEFEM_ASSERT(xt::amax(iip)() <= xt::amax(dofs)());

        m_iip = iip;
        m_iiu = xt::setdiff1d(dofs, iip);
        m_nnp = m_iip.size();
        m_nnu = m_iiu.size();
        m_part = Mesh::Reorder({m_iiu, m_iip}).apply(m_dofs);
        m_Auu = xt::empty<double>({m_nnu});
        m_App = xt::empty<double>({m_nnp});
        m_inv_uu = xt::empty<double>({m_nnu});
    }

    /**
    Number of unknown DOFs.
    */
    size_t nnu() const
    {
        return m_nnu;
    }

    /**
    Number of prescribed DOFs.
    */
    size_t nnp() const
    {
        return m_nnp;
    }

    /**
    Unknown DOFs [#nnu].
    */
    array_type::tensor<size_t, 1> iiu() const
    {
        return m_iiu;
    }

    /**
    Prescribed DOFs [#nnp].
    */
    array_type::tensor<size_t, 1> iip() const
    {
        return m_iip;
    }

    void assemble(const array_type::tensor<double, 3>& elemmat) override
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));
        GOOSEFEM_ASSERT(Element::isDiagonal(elemmat));

        m_Auu.fill(0.0);
        m_App.fill(0.0);

        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {

                    size_t d = m_part(m_conn(e, m), i);

                    if (d < m_nnu) {
                        m_Auu(d) += elemmat(e, m * m_ndim + i, m * m_ndim + i);
                    }
                    else {
                        m_App(d - m_nnu) += elemmat(e, m * m_ndim + i, m * m_ndim + i);
                    }
                }
            }
        }

        m_changed = true;
    }

    void set(const array_type::tensor<double, 1>& A) override
    {
        GOOSEFEM_ASSERT(A.size() == m_ndof);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            m_Auu(d) = A(m_iiu(d));
        }

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            m_App(d) = A(m_iip(d));
        }

        m_changed = true;
    }

    array_type::tensor<double, 1> Todiagonal() const override
    {
        array_type::tensor<double, 1> ret = xt::zeros<double>({m_ndof});

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            ret(m_iiu(d)) = m_Auu(d);
        }

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            ret(m_iip(d)) = m_App(d);
        }

        return ret;
    }

    void
    dot(const array_type::tensor<double, 2>& x, array_type::tensor<double, 2>& b) const override
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {

                size_t d = m_part(m, i);

                if (d < m_nnu) {
                    b(m, i) = m_Auu(d) * x(m, i);
                }
                else {
                    b(m, i) = m_App(d - m_nnu) * x(m, i);
                }
            }
        }
    }

    void
    dot(const array_type::tensor<double, 1>& x, array_type::tensor<double, 1>& b) const override
    {
        GOOSEFEM_ASSERT(x.size() == m_ndof);
        GOOSEFEM_ASSERT(b.size() == m_ndof);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            b(m_iiu(d)) = m_Auu(d) * x(m_iiu(d));
        }

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            b(m_iip(d)) = m_App(d) * x(m_iip(d));
        }
    }

    /**
    \param x_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \return b_u dofval [#nnu].
    */
    array_type::tensor<double, 1>
    Dot_u(const array_type::tensor<double, 1>& x_u, const array_type::tensor<double, 1>& x_p) const
    {
        array_type::tensor<double, 1> b_u = xt::empty<double>({m_nnu});
        this->dot_u(x_u, x_p, b_u);
        return b_u;
    }

    /**
    \param x_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \param b_u (overwritten) dofval [#nnu].
    */
    void dot_u(
        const array_type::tensor<double, 1>& x_u,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& b_u) const
    {
        UNUSED(x_p);

        GOOSEFEM_ASSERT(x_u.size() == m_nnu);
        GOOSEFEM_ASSERT(x_p.size() == m_nnp);
        GOOSEFEM_ASSERT(b_u.size() == m_nnu);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            b_u(d) = m_Auu(d) * x_u(d);
        }
    }

    /**
    \param x_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \return b_p dofval [#nnp].
    */
    array_type::tensor<double, 1>
    Dot_p(const array_type::tensor<double, 1>& x_u, const array_type::tensor<double, 1>& x_p) const
    {
        array_type::tensor<double, 1> b_p = xt::empty<double>({m_nnp});
        this->dot_p(x_u, x_p, b_p);
        return b_p;
    }

    /**
    \param x_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \param b_p (overwritten) dofval [#nnp].
    */
    void dot_p(
        const array_type::tensor<double, 1>& x_u,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& b_p) const
    {
        UNUSED(x_u);

        GOOSEFEM_ASSERT(x_u.size() == m_nnu);
        GOOSEFEM_ASSERT(x_p.size() == m_nnp);
        GOOSEFEM_ASSERT(b_p.size() == m_nnp);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            b_p(d) = m_App(d) * x_p(d);
        }
    }

    /**
    Solve \f$ x_u = A_{uu}^{-1} (b_u - A_{up} * x_p) \equiv A_{uu}^{-1} b_u \f$.

    \param b nodevec [#nelem, #ndim].
    \param x nodevec, modified with `x_u` [#nelem, #ndim].
    */
    void solve(const array_type::tensor<double, 2>& b, array_type::tensor<double, 2>& x) override
    {
        GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));

        this->factorize();

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) < m_nnu) {
                    x(m, i) = m_inv_uu(m_part(m, i)) * b(m, i);
                }
            }
        }
    }

    /**
    Solve \f$ x_u = A_{uu}^{-1} (b_u - A_{up} * x_p) \equiv A_{uu}^{-1} b_u \f$.

    \param b dofval [#ndof].
    \param x dofval, modified with `x_u` [#ndof].
    */
    void solve(const array_type::tensor<double, 1>& b, array_type::tensor<double, 1>& x) override
    {
        GOOSEFEM_ASSERT(b.size() == m_ndof);
        GOOSEFEM_ASSERT(x.size() == m_ndof);

        this->factorize();

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            x(m_iiu(d)) = m_inv_uu(d) * b(m_iiu(d));
        }
    }

    /**
    \param b_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \return x_u dofval [#nnu].
    */
    array_type::tensor<double, 1>
    Solve_u(const array_type::tensor<double, 1>& b_u, const array_type::tensor<double, 1>& x_p)
    {
        array_type::tensor<double, 1> x_u = xt::empty<double>({m_nnu});
        this->solve_u(b_u, x_p, x_u);
        return x_u;
    }

    /**
    \param b_u dofval [#nnu].
    \param x_p dofval [#nnp].
    \param x_u (overwritten) dofval [#nnu].
    */
    void solve_u(
        const array_type::tensor<double, 1>& b_u,
        const array_type::tensor<double, 1>& x_p,
        array_type::tensor<double, 1>& x_u)
    {
        UNUSED(x_p);

        GOOSEFEM_ASSERT(b_u.size() == m_nnu);
        GOOSEFEM_ASSERT(x_p.size() == m_nnp);
        GOOSEFEM_ASSERT(x_u.size() == m_nnu);

        this->factorize();

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            x_u(d) = m_inv_uu(d) * b_u(d);
        }
    }

    /**
    Get right-hand-size for corresponding to the prescribed DOFs.

    \f$ b_p = A_{pu} * x_u + A_{pp} * x_p = A_{pp} * x_p \equiv A_{pp} * x_p \f$

    and assemble them to the appropriate places in "nodevec".

    \param x "nodevec" [#nnode, #ndim].
    \param b "nodevec" [#nnode, #ndim].
    \return Copy of `b` with \f$ b_p \f$ overwritten.
    */
    array_type::tensor<double, 2>
    Reaction(const array_type::tensor<double, 2>& x, const array_type::tensor<double, 2>& b) const
    {
        array_type::tensor<double, 2> ret = b;
        this->reaction(x, ret);
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
        GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) >= m_nnu) {
                    b(m, i) = m_App(m_part(m, i) - m_nnu) * x(m, i);
                }
            }
        }
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
        array_type::tensor<double, 1> ret = b;
        this->reaction(x, ret);
        return ret;
    }

    /**
    Same as Reaction(const array_type::tensor<double, 1>&, const array_type::tensor<double, 1>&),
    but inserting in-place.

    \param x "dofval" [#ndof].
    \param b "dofval" [#ndof], \f$ b_p \f$ overwritten.
    */
    void reaction(const array_type::tensor<double, 1>& x, array_type::tensor<double, 1>& b) const
    {
        GOOSEFEM_ASSERT(x.size() == m_ndof);
        GOOSEFEM_ASSERT(b.size() == m_ndof);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            b(m_iip(d)) = m_App(d) * x(m_iip(d));
        }
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
        this->reaction_p(x_u, x_p, b_p);
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
        UNUSED(x_u);

        GOOSEFEM_ASSERT(x_u.size() == m_nnu);
        GOOSEFEM_ASSERT(x_p.size() == m_nnp);
        GOOSEFEM_ASSERT(b_p.size() == m_nnp);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            b_p(d) = m_App(d) * x_p(d);
        }
    }

private:
    // The diagonal matrix, and its inverse (re-used to solve different RHS)
    array_type::tensor<double, 1> m_Auu;
    array_type::tensor<double, 1> m_App;
    array_type::tensor<double, 1> m_inv_uu;

    // Bookkeeping
    array_type::tensor<size_t, 2> m_part; // DOF-numbers per node, renumbered  [nnode, ndim]
    array_type::tensor<size_t, 1> m_iiu; // DOF-numbers that are unknown      [nnu]
    array_type::tensor<size_t, 1> m_iip; // DOF-numbers that are prescribed   [nnp]

    // Dimensions
    size_t m_nnu; // number of unknown DOFs
    size_t m_nnp; // number of prescribed DOFs

    // Compute inverse (automatically evaluated by "solve")
    void factorize()
    {
        if (!m_changed) {
            return;
        }

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            m_inv_uu(d) = 1.0 / m_Auu(d);
        }

        m_changed = false;
    }
};

} // namespace GooseFEM

#endif

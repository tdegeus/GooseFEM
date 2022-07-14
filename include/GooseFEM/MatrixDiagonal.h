/**
Diagonal matrix.

\file MatrixDiagonal.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MATRIXDIAGONAL_H
#define GOOSEFEM_MATRIXDIAGONAL_H

#include "Element.h"
#include "config.h"

namespace GooseFEM {

/**
Diagonal matrix.

See Vector() for bookkeeping definitions.
*/
class MatrixDiagonal {
public:
    MatrixDiagonal() = default;

    virtual ~MatrixDiagonal() = default;

    /**
    Constructor.

    \tparam C e.g. `array_type::tensor<size_t, 2>`
    \tparam D e.g. `array_type::tensor<size_t, 2>`
    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    */
    template <class C, class D>
    MatrixDiagonal(const C& conn, const D& dofs) : m_conn(conn), m_dofs(dofs)
    {
        m_nelem = m_conn.shape(0);
        m_nne = m_conn.shape(1);
        m_nnode = m_dofs.shape(0);
        m_ndim = m_dofs.shape(1);
        m_ndof = xt::amax(m_dofs)() + 1;
        m_A = xt::empty<double>({m_ndof});
        m_inv = xt::empty<double>({m_ndof});

        GOOSEFEM_ASSERT(xt::amax(m_conn)() + 1 <= m_nnode);
        GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
    }

    /**
    \return Number of elements.
    */
    size_t nelem() const
    {
        return m_nelem;
    }

    /**
    \return Number of nodes per element.
    */
    size_t nne() const
    {
        return m_nne;
    }

    /**
    \return Number of nodes.
    */
    size_t nnode() const
    {
        return m_nnode;
    }

    /**
    \return Number of dimensions.
    */
    size_t ndim() const
    {
        return m_ndim;
    }

    /**
    \return Number of DOFs.
    */
    size_t ndof() const
    {
        return m_ndof;
    }

    /**
    \return DOFs per node [#nnode, #ndim]
    */
    array_type::tensor<size_t, 2> dofs() const
    {
        return m_dofs;
    }

    /**
    Assemble from matrices stored per element.
    \warning Ignores any off-diagonal terms.

    \param elemmat [#nelem, #nne * #ndim, #nne * #ndim].
    */
    virtual void assemble(const array_type::tensor<double, 3>& elemmat)
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemmat, {m_nelem, m_nne * m_ndim, m_nne * m_ndim}));
        GOOSEFEM_ASSERT(Element::isDiagonal(elemmat));

        m_A.fill(0.0);

        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    m_A(m_dofs(m_conn(e, m), i)) += elemmat(e, m * m_ndim + i, m * m_ndim + i);
                }
            }
        }

        m_changed = true;
    }

    /**
    Set all (diagonal) matrix components.

    \param A The matrix [#ndof].
    */
    virtual void set(const array_type::tensor<double, 1>& A)
    {
        GOOSEFEM_ASSERT(A.size() == m_ndof);
        std::copy(A.begin(), A.end(), m_A.begin());
        m_changed = true;
    }

    /**
    Return matrix as diagonal matrix.

    \param [#ndof].
    */
    virtual array_type::tensor<double, 1> Todiagonal() const
    {
        return m_A;
    }

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x nodevec [#nelem, #ndim].
    \return b nodevec overwritten [#nelem, #ndim].
    */
    array_type::tensor<double, 2> Dot(const array_type::tensor<double, 2>& x) const
    {
        array_type::tensor<double, 2> b = xt::empty<double>({m_nnode, m_ndim});
        this->dot(x, b);
        return b;
    }

    /**
    Same as Dot(const array_type::tensor<double, 2>&, array_type::tensor<double, 2>& b)
    but writing to preallocated data.

    \param x nodevec [#nelem, #ndim].
    \param b nodevec overwritten [#nelem, #ndim].
    */
    virtual void dot(const array_type::tensor<double, 2>& x, array_type::tensor<double, 2>& b) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                b(m, i) = m_A(m_dofs(m, i)) * x(m, i);
            }
        }
    }

    /**
    Dot-product \f$ b_i = A_{ij} x_j \f$.

    \param x dofval [#ndof].
    \return b dofval overwritten [#ndof].
    */
    array_type::tensor<double, 1> Dot(const array_type::tensor<double, 1>& x) const
    {
        array_type::tensor<double, 1> b = xt::empty<double>({m_ndof});
        this->dot(x, b);
        return b;
    }

    /**
    Same as Dot(const array_type::tensor<double, 1>&, array_type::tensor<double, 1>& b)
    but writing to preallocated data.

    \param x dofval [#ndof].
    \param b dofval overwritten [#ndof].
    */
    virtual void dot(const array_type::tensor<double, 1>& x, array_type::tensor<double, 1>& b) const
    {
        GOOSEFEM_ASSERT(x.size() == m_ndof);
        GOOSEFEM_ASSERT(b.size() == m_ndof);

        xt::noalias(b) = m_A * x;
    }

    /**
    Solve \f$ x = A^{-1} b \f$.

    \param b nodevec [nelem, ndim].
    \return x nodevec [nelem, ndim].
    */
    array_type::tensor<double, 2> Solve(const array_type::tensor<double, 2>& b)
    {
        array_type::tensor<double, 2> x = xt::empty<double>({m_nnode, m_ndim});
        this->solve(b, x);
        return x;
    }

    /**
    Same as Solve(const array_type::tensor<double, 2>&)
    but writing to preallocated data.

    \param b nodevec [nelem, ndim].
    \param x nodevec overwritten [nelem, ndim].
    */
    virtual void solve(const array_type::tensor<double, 2>& b, array_type::tensor<double, 2>& x)
    {
        GOOSEFEM_ASSERT(xt::has_shape(b, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(x, {m_nnode, m_ndim}));

        this->factorize();

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                x(m, i) = m_inv(m_dofs(m, i)) * b(m, i);
            }
        }
    }

    /**
    Same as Solve(const array_type::tensor<double, 2>&)
    but for "dofval" input and output.

    \param b dofval [ndof].
    \return x dofval [ndof].
    */
    array_type::tensor<double, 1> Solve(const array_type::tensor<double, 1>& b)
    {
        array_type::tensor<double, 1> x = xt::empty<double>({m_ndof});
        this->solve(b, x);
        return x;
    }

    /**
    Same as Solve(const array_type::tensor<double, 1>&)
    but writing to preallocated data.

    \param b dofval [ndof].
    \param x dofval overwritten [ndof].
    */
    virtual void solve(const array_type::tensor<double, 1>& b, array_type::tensor<double, 1>& x)
    {
        GOOSEFEM_ASSERT(b.size() == m_ndof);
        GOOSEFEM_ASSERT(x.size() == m_ndof);
        this->factorize();
        xt::noalias(x) = m_inv * b;
    }

protected:
    array_type::tensor<size_t, 2> m_conn; ///< Connectivity [#nelem, #nne].
    array_type::tensor<size_t, 2> m_dofs; ///< DOF-numbers per node [#nnode, #ndim].
    size_t m_nelem; ///< See nelem().
    size_t m_nne; ///< See nne().
    size_t m_nnode; ///< See nnode().
    size_t m_ndim; ///< See ndim().
    size_t m_ndof; ///< See ndof().
    bool m_changed = true; ///< Signal changes to data.

private:
    array_type::tensor<double, 1> m_A; ///< The matrix.
    array_type::tensor<double, 1> m_inv; /// Inverse of the matrix.

    /**
    Compute inverse (automatically evaluated by "solve").
    */
    void factorize()
    {
        if (!m_changed) {
            return;
        }

#pragma omp parallel for
        for (size_t d = 0; d < m_ndof; ++d) {
            m_inv(d) = 1.0 / m_A(d);
        }

        m_changed = false;
    }
};

} // namespace GooseFEM

#endif

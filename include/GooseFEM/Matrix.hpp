/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_MATRIX_HPP
#define GOOSEFEM_MATRIX_HPP

#include "Matrix.h"

namespace GooseFEM {

template <class Solver>
inline Matrix<Solver>::Matrix(
    const xt::xtensor<size_t,2>& conn, const xt::xtensor<size_t,2>& dofs)
    : m_conn(conn), m_dofs(dofs)
{
    m_nelem = m_conn.shape(0);
    m_nne = m_conn.shape(1);
    m_nnode = m_dofs.shape(0);
    m_ndim = m_dofs.shape(1);
    m_ndof = xt::amax(m_dofs)[0] + 1;
    m_T.reserve(m_nelem * m_nne * m_ndim * m_nne * m_ndim);
    m_A.resize(m_ndof, m_ndof);

    GOOSEFEM_ASSERT(xt::amax(m_conn)[0] + 1 == m_nnode);
    GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
}

template <class Solver>
inline size_t Matrix<Solver>::nelem() const
{
    return m_nelem;
}

template <class Solver>
inline size_t Matrix<Solver>::nne() const
{
    return m_nne;
}

template <class Solver>
inline size_t Matrix<Solver>::nnode() const
{
    return m_nnode;
}

template <class Solver>
inline size_t Matrix<Solver>::ndim() const
{
    return m_ndim;
}

template <class Solver>
inline size_t Matrix<Solver>::ndof() const
{
    return m_ndof;
}

template <class Solver>
inline xt::xtensor<size_t,2> Matrix<Solver>::dofs() const
{
    return m_dofs;
}

template <class Solver>
inline void Matrix<Solver>::factorize()
{
    if (!m_factor) {
        return;
    }

    m_solver.compute(m_A);
    m_factor = false;
}

template <class Solver>
inline void Matrix<Solver>::assemble(const xt::xtensor<double,3>& elemmat)
{
    GOOSEFEM_ASSERT(
        elemmat.shape() ==
        std::decay_t<decltype(elemmat)>::shape_type({m_nelem, m_nne * m_ndim, m_nne * m_ndim}));

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

    m_factor = true;
}

template <class Solver>
inline void Matrix<Solver>::solve(const xt::xtensor<double,2>& b, xt::xtensor<double,2>& x)
{
    GOOSEFEM_ASSERT(b.shape() == std::decay_t<decltype(b)>::shape_type({m_nnode, m_ndim}));
    GOOSEFEM_ASSERT(x.shape() == std::decay_t<decltype(x)>::shape_type({m_nnode, m_ndim}));

    this->factorize();
    Eigen::VectorXd B = this->asDofs(b);
    Eigen::VectorXd X = m_solver.solve(B);
    this->asNode(X, x);
}

template <class Solver>
inline void Matrix<Solver>::solve(const xt::xtensor<double,1>& b, xt::xtensor<double,1>& x)
{
    GOOSEFEM_ASSERT(b.size() == m_ndof);
    GOOSEFEM_ASSERT(x.size() == m_ndof);

    this->factorize();
    Eigen::VectorXd B(m_ndof, 1);
    std::copy(b.begin(), b.end(), B.data());
    Eigen::VectorXd X = m_solver.solve(B);
    std::copy(X.data(), X.data() + m_ndof, x.begin());
}

template <class Solver>
inline xt::xtensor<double,2> Matrix<Solver>::Solve(const xt::xtensor<double,2>& b)
{
    xt::xtensor<double,2> x = xt::empty<double>({m_nnode, m_ndim});
    this->solve(b, x);
    return x;
}

template <class Solver>
inline xt::xtensor<double,1> Matrix<Solver>::Solve(const xt::xtensor<double,1>& b)
{
    xt::xtensor<double,1> x = xt::empty<double>({m_ndof});
    this->solve(b, x);
    return x;
}

template <class Solver>
inline Eigen::VectorXd Matrix<Solver>::asDofs(const xt::xtensor<double,2>& nodevec) const
{
    assert(nodevec.shape() == std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

    Eigen::VectorXd dofval(m_ndof, 1);

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            dofval(m_dofs(m, i)) = nodevec(m, i);
        }
    }

    return dofval;
}

template <class Solver>
inline void
Matrix<Solver>::asNode(const Eigen::VectorXd& dofval, xt::xtensor<double,2>& nodevec) const
{
    assert(static_cast<size_t>(dofval.size()) == m_ndof);
    assert(nodevec.shape() == std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

    #pragma omp parallel for
    for (size_t m = 0; m < m_nnode; ++m) {
        for (size_t i = 0; i < m_ndim; ++i) {
            nodevec(m, i) = dofval(m_dofs(m, i));
        }
    }
}

} // namespace GooseFEM

#endif

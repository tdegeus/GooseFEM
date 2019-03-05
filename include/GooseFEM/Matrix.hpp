/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIX_HPP
#define GOOSEFEM_MATRIX_HPP

// -------------------------------------------------------------------------------------------------

#include "Matrix.h"

// =================================================================================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

inline Matrix::Matrix(
  const xt::xtensor<size_t,2> &conn,
  const xt::xtensor<size_t,2> &dofs) :
  m_conn(conn), m_dofs(dofs)
{
  m_nelem = m_conn.shape()[0];
  m_nne   = m_conn.shape()[1];
  m_nnode = m_dofs.shape()[0];
  m_ndim  = m_dofs.shape()[1];

  m_ndof  = xt::amax(m_dofs)[0] + 1;

  m_T.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);

  m_A.resize(m_ndof,m_ndof);

  assert(xt::amax(m_conn)[0] + 1 == m_nnode);
  assert(m_ndof <= m_nnode * m_ndim);
}

// -------------------------------------------------------------------------------------------------

inline size_t Matrix::nelem() const
{ return m_nelem; }

inline size_t Matrix::nne() const
{ return m_nne; }

inline size_t Matrix::nnode() const
{ return m_nnode; }

inline size_t Matrix::ndim() const
{ return m_ndim; }

inline size_t Matrix::ndof() const
{ return m_ndof; }

inline xt::xtensor<size_t,2> Matrix::dofs() const
{ return m_dofs; }

// -------------------------------------------------------------------------------------------------

inline void Matrix::factorize()
{
  if ( ! m_factor ) return;

  m_solver.compute(m_A);

  m_factor = false;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::assemble(const xt::xtensor<double,3> &elemmat)
{
  assert(elemmat.shape()[0] == m_nelem     );
  assert(elemmat.shape()[1] == m_nne*m_ndim);
  assert(elemmat.shape()[2] == m_nne*m_ndim);

  m_T.clear();

  for (size_t e = 0 ; e < m_nelem ; ++e)
    for (size_t m = 0 ; m < m_nne ; ++m)
      for (size_t i = 0 ; i < m_ndim ; ++i)
        for (size_t n = 0 ; n < m_nne ; ++n)
          for (size_t j = 0 ; j < m_ndim ; ++j)
            m_T.push_back(Eigen::Triplet<double>(
              m_dofs(m_conn(e,m),i), m_dofs(m_conn(e,n),j), elemmat(e,m*m_ndim+i,n*m_ndim+j)
            ));

  m_A.setFromTriplets(m_T.begin(), m_T.end());

  m_factor = true;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::solve(
  const xt::xtensor<double,2> &b,
        xt::xtensor<double,2> &x)
{
  assert(b.shape()[0] == m_nnode);
  assert(b.shape()[1] == m_ndim );
  assert(x.shape()[0] == m_nnode);
  assert(x.shape()[1] == m_ndim );

  this->factorize();

  Eigen::VectorXd B = this->asDofs(b);

  Eigen::VectorXd X = m_solver.solve(B);

  this->asNode(X, x);
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::solve(
  const xt::xtensor<double,1> &b,
        xt::xtensor<double,1> &x)
{
  assert(b.size() == m_ndof);
  assert(x.size() == m_ndof);

  this->factorize();

  Eigen::VectorXd B(m_ndof,1);

  std::copy(b.begin(), b.end(), B.data());

  Eigen::VectorXd X = m_solver.solve(B);

  std::copy(X.data(), X.data()+m_ndof, x.begin());
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Matrix::Solve(
  const xt::xtensor<double,2> &b)
{
  xt::xtensor<double,2> x = xt::empty<double>({m_nnode, m_ndim});

  this->solve(b, x);

  return x;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Matrix::Solve(
  const xt::xtensor<double,1> &b)
{
  xt::xtensor<double,1> x = xt::empty<double>({m_ndof});

  this->solve(b, x);

  return x;
}

// -------------------------------------------------------------------------------------------------

inline Eigen::VectorXd Matrix::asDofs(const xt::xtensor<double,2> &nodevec) const
{
  assert(nodevec.shape()[0] == m_nnode);
  assert(nodevec.shape()[1] == m_ndim );

  Eigen::VectorXd dofval(m_ndof,1);

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      dofval(m_dofs(m,i)) = nodevec(m,i);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::asNode(const Eigen::VectorXd &dofval, xt::xtensor<double,2> &nodevec) const
{
  assert(static_cast<size_t>(dofval.size()) == m_ndof);
  assert(nodevec.shape()[0] == m_nnode);
  assert(nodevec.shape()[1] == m_ndim );

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      nodevec(m,i) = dofval(m_dofs(m,i));
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

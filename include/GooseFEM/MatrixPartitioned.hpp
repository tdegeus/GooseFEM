/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIXPARTITIONED_HPP
#define GOOSEFEM_MATRIXPARTITIONED_HPP

// -------------------------------------------------------------------------------------------------

#include "MatrixPartitioned.h"
#include "Mesh.h"

// =================================================================================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline MatrixPartitioned<Solver>::MatrixPartitioned(
  const xt::xtensor<size_t,2> &conn,
  const xt::xtensor<size_t,2> &dofs,
  const xt::xtensor<size_t,1> &iip) :
  m_conn(conn), m_dofs(dofs), m_iip(iip)
{
  m_nelem = m_conn.shape()[0];
  m_nne   = m_conn.shape()[1];
  m_nnode = m_dofs.shape()[0];
  m_ndim  = m_dofs.shape()[1];

  m_iiu   = xt::setdiff1d(dofs, iip);

  m_ndof  = xt::amax(m_dofs)[0] + 1;
  m_nnp   = m_iip.size();
  m_nnu   = m_iiu.size();

  m_part  = Mesh::Reorder({m_iiu, m_iip}).get(m_dofs);

  m_Tuu.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tup.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tpu.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tpp.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);

  m_Auu.resize(m_nnu,m_nnu);
  m_Aup.resize(m_nnu,m_nnp);
  m_Apu.resize(m_nnp,m_nnu);
  m_App.resize(m_nnp,m_nnp);

  GOOSEFEM_ASSERT(xt::amax(m_conn)[0] + 1 == m_nnode);
  GOOSEFEM_ASSERT(xt::amax(m_iip)[0] <= xt::amax(m_dofs)[0]);
  GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitioned<Solver>::nelem() const
{ return m_nelem; }

template <class Solver>
inline size_t MatrixPartitioned<Solver>::nne() const
{ return m_nne; }

template <class Solver>
inline size_t MatrixPartitioned<Solver>::nnode() const
{ return m_nnode; }

template <class Solver>
inline size_t MatrixPartitioned<Solver>::ndim() const
{ return m_ndim; }

template <class Solver>
inline size_t MatrixPartitioned<Solver>::ndof() const
{ return m_ndof; }

template <class Solver>
inline size_t MatrixPartitioned<Solver>::nnu() const
{ return m_nnu; }

template <class Solver>
inline size_t MatrixPartitioned<Solver>::nnp() const
{ return m_nnp; }

template <class Solver>
inline xt::xtensor<size_t,2> MatrixPartitioned<Solver>::dofs() const
{ return m_dofs; }

template <class Solver>
inline xt::xtensor<size_t,1> MatrixPartitioned<Solver>::iiu() const
{ return m_iiu; }

template <class Solver>
inline xt::xtensor<size_t,1> MatrixPartitioned<Solver>::iip() const
{ return m_iip; }

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitioned<Solver>::factorize()
{
  if (!m_factor) return;


  m_solver.compute(m_Auu);

  m_factor = false;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitioned<Solver>::assemble(const xt::xtensor<double,3> &elemmat)
{
  GOOSEFEM_ASSERT(elemmat.shape() ==\
    std::decay_t<decltype(elemmat)>::shape_type({m_nelem, m_nne*m_ndim, m_nne*m_ndim}));

  m_Tuu.clear();
  m_Tup.clear();
  m_Tpu.clear();
  m_Tpp.clear();

  for (size_t e = 0 ; e < m_nelem ; ++e) {
    for (size_t m = 0 ; m < m_nne ; ++m) {
      for (size_t i = 0 ; i < m_ndim ; ++i) {

        size_t di = m_part(m_conn(e,m),i);

        for (size_t n = 0 ; n < m_nne ; ++n) {
          for (size_t j = 0 ; j < m_ndim ; ++j) {

            size_t dj = m_part(m_conn(e,n),j);

            if      (di < m_nnu and dj < m_nnu)
              m_Tuu.push_back(Eigen::Triplet<double>(di      ,dj      ,elemmat(e,m*m_ndim+i,n*m_ndim+j)));
            else if (di < m_nnu)
              m_Tup.push_back(Eigen::Triplet<double>(di      ,dj-m_nnu,elemmat(e,m*m_ndim+i,n*m_ndim+j)));
            else if (dj < m_nnu)
              m_Tpu.push_back(Eigen::Triplet<double>(di-m_nnu,dj      ,elemmat(e,m*m_ndim+i,n*m_ndim+j)));
            else
              m_Tpp.push_back(Eigen::Triplet<double>(di-m_nnu,dj-m_nnu,elemmat(e,m*m_ndim+i,n*m_ndim+j)));
          }
        }
      }
    }
  }

  m_Auu.setFromTriplets(m_Tuu.begin(), m_Tuu.end());
  m_Aup.setFromTriplets(m_Tup.begin(), m_Tup.end());
  m_Apu.setFromTriplets(m_Tpu.begin(), m_Tpu.end());
  m_App.setFromTriplets(m_Tpp.begin(), m_Tpp.end());

  m_factor = true;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitioned<Solver>::solve(
  const xt::xtensor<double,2> &b,
        xt::xtensor<double,2> &x)
{
  GOOSEFEM_ASSERT(b.shape() ==\
    std::decay_t<decltype(b)>::shape_type({m_nnode, m_ndim}));
  GOOSEFEM_ASSERT(x.shape() ==\
    std::decay_t<decltype(x)>::shape_type({m_nnode, m_ndim}));

  this->factorize();

  Eigen::VectorXd B_u = this->asDofs_u(b);
  Eigen::VectorXd X_p = this->asDofs_p(x);

  Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - m_Aup * X_p));

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      if ( m_part(m,i) < m_nnu )
        x(m,i) = X_u(m_part(m,i));
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitioned<Solver>::solve(
  const xt::xtensor<double,1> &b,
        xt::xtensor<double,1> &x)
{
  GOOSEFEM_ASSERT(b.size() == m_ndof);
  GOOSEFEM_ASSERT(x.size() == m_ndof);

  this->factorize();

  Eigen::VectorXd B_u = this->asDofs_u(b);
  Eigen::VectorXd X_p = this->asDofs_p(x);

  Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - m_Aup * X_p));

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnu ; ++d)
    x(m_iiu(d)) = X_u(d);
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitioned<Solver>::solve_u(
  const xt::xtensor<double,1> &b_u,
  const xt::xtensor<double,1> &x_p,
        xt::xtensor<double,1> &x_u)
{
  GOOSEFEM_ASSERT(b_u.size() == m_nnu);
  GOOSEFEM_ASSERT(x_p.size() == m_nnp);
  GOOSEFEM_ASSERT(x_u.size() == m_nnu);

  this->factorize();

  Eigen::VectorXd B_u(m_nnu,1);
  Eigen::VectorXd X_p(m_nnp,1);

  std::copy(b_u.begin(), b_u.end(), B_u.data());
  std::copy(x_p.begin(), x_p.end(), X_p.data());

  Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - m_Aup * X_p));

  std::copy(X_u.data(), X_u.data()+m_nnu, x_u.begin());
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitioned<Solver>::reaction(
  const xt::xtensor<double,2> &x,
        xt::xtensor<double,2> &b) const
{
  GOOSEFEM_ASSERT(x.shape() ==\
    std::decay_t<decltype(x)>::shape_type({m_nnode, m_ndim}));
  GOOSEFEM_ASSERT(b.shape() ==\
    std::decay_t<decltype(b)>::shape_type({m_nnode, m_ndim}));

  Eigen::VectorXd X_u = this->asDofs_u(x);
  Eigen::VectorXd X_p = this->asDofs_p(x);

  Eigen::VectorXd B_p = m_Apu * X_u + m_App * X_p;

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      if (m_part(m,i) >= m_nnu)
        b(m,i) = B_p(m_part(m,i)-m_nnu);
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitioned<Solver>::reaction(
  const xt::xtensor<double,1> &x,
        xt::xtensor<double,1> &b) const
{
  GOOSEFEM_ASSERT(x.size() == m_ndof);
  GOOSEFEM_ASSERT(b.size() == m_ndof);

  Eigen::VectorXd X_u = this->asDofs_u(x);
  Eigen::VectorXd X_p = this->asDofs_p(x);

  Eigen::VectorXd B_p = m_Apu * X_u + m_App * X_p;

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnp ; ++d)
    b(m_iip(d)) = B_p(d);
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitioned<Solver>::reaction_p(
  const xt::xtensor<double,1> &x_u,
  const xt::xtensor<double,1> &x_p,
        xt::xtensor<double,1> &b_p) const
{
  GOOSEFEM_ASSERT(x_u.size() == m_nnu);
  GOOSEFEM_ASSERT(x_p.size() == m_nnp);
  GOOSEFEM_ASSERT(b_p.size() == m_nnp);

  Eigen::VectorXd X_u(m_nnu,1);
  Eigen::VectorXd X_p(m_nnp,1);

  std::copy(x_u.begin(), x_u.end(), X_u.data());
  std::copy(x_p.begin(), x_p.end(), X_p.data());

  Eigen::VectorXd B_p = m_Apu * X_u + m_App * X_p;

  std::copy(B_p.data(), B_p.data()+m_nnp, b_p.begin());
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,2> MatrixPartitioned<Solver>::Solve(
  const xt::xtensor<double,2> &b,
  const xt::xtensor<double,2> &x)
{
  xt::xtensor<double,2> out = x;
  this->solve(b, out);
  return out;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,1> MatrixPartitioned<Solver>::Solve(
  const xt::xtensor<double,1> &b,
  const xt::xtensor<double,1> &x)
{
  xt::xtensor<double,1> out = x;
  this->solve(b, out);
  return out;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,1> MatrixPartitioned<Solver>::Solve_u(
  const xt::xtensor<double,1> &b_u,
  const xt::xtensor<double,1> &x_p)
{
  xt::xtensor<double,1> x_u = xt::empty<double>({m_nnu});
  this->solve_u(b_u, x_p, x_u);
  return x_u;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,2> MatrixPartitioned<Solver>::Reaction(
  const xt::xtensor<double,2> &x,
  const xt::xtensor<double,2> &b) const
{
  xt::xtensor<double,2> out = b;
  this->reaction(x, out);
  return out;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,1> MatrixPartitioned<Solver>::Reaction(
  const xt::xtensor<double,1> &x,
  const xt::xtensor<double,1> &b) const
{
  xt::xtensor<double,1> out = b;
  this->reaction(x, out);
  return out;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,1> MatrixPartitioned<Solver>::Reaction_p(
  const xt::xtensor<double,1> &x_u,
  const xt::xtensor<double,1> &x_p) const
{
  xt::xtensor<double,1> b_p = xt::empty<double>({m_nnp});
  this->reaction_p(x_u, x_p, b_p);
  return b_p;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitioned<Solver>::asDofs_u(
  const xt::xtensor<double,1> &dofval) const
{
  assert(dofval.size() == m_ndof);

  Eigen::VectorXd dofval_u(m_nnu,1);

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnu ; ++d)
    dofval_u(d) = dofval(m_iiu(d));

  return dofval_u;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitioned<Solver>::asDofs_u(
  const xt::xtensor<double,2> &nodevec) const
{
  assert(nodevec.shape() ==\
    std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

  Eigen::VectorXd dofval_u(m_nnu,1);

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      if (m_part(m,i) < m_nnu)
        dofval_u(m_part(m,i)) = nodevec(m,i);

  return dofval_u;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitioned<Solver>::asDofs_p(
  const xt::xtensor<double,1> &dofval) const
{
  assert(dofval.size() == m_ndof);

  Eigen::VectorXd dofval_p(m_nnp,1);

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnp ; ++d)
    dofval_p(d) = dofval(m_iip(d));

  return dofval_p;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitioned<Solver>::asDofs_p(
  const xt::xtensor<double,2> &nodevec) const
{
  assert(nodevec.shape() ==\
    std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

  Eigen::VectorXd dofval_p(m_nnp,1);

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      if (m_part(m,i) >= m_nnu)
        dofval_p(m_part(m,i)-m_nnu) = nodevec(m,i);

  return dofval_p;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

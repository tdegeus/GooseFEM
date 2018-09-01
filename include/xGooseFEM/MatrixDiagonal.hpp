/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_MATRIXDIAGONAL_CPP
#define XGOOSEFEM_MATRIXDIAGONAL_CPP

// -------------------------------------------------------------------------------------------------

#include "MatrixDiagonal.h"

// =========================================== GooseFEM ============================================

namespace xGooseFEM {

// ------------------------------------------ constructor ------------------------------------------

inline MatrixDiagonal::MatrixDiagonal(const xt::xtensor<size_t,2> &conn,
  const xt::xtensor<size_t,2> &dofs) : m_conn(conn), m_dofs(dofs)
{
  // mesh dimensions
  m_nelem = m_conn.shape()[0];
  m_nne   = m_conn.shape()[1];
  m_nnode = m_dofs.shape()[0];
  m_ndim  = m_dofs.shape()[1];

  // list with prescribed and unknown DOFs
  m_iip   = xt::empty<size_t>({0});
  m_iiu   = xt::unique(dofs);

  // dimensions of the system
  m_ndof  = xt::amax(m_dofs)[0] + 1;
  m_nnp   = m_iip.size();
  m_nnu   = m_iiu.size();

  // allocate matrix and its inverse
  m_data  = xt::empty<double>({m_ndof});
  m_inv   = xt::empty<double>({m_ndof});

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode );
  assert( m_ndof <= m_nnode * m_ndim );
}

// ------------------------------------------ constructor ------------------------------------------

inline MatrixDiagonal::MatrixDiagonal(const xt::xtensor<size_t,2> &conn,
  const xt::xtensor<size_t,2> &dofs, const xt::xtensor<size_t,1> &iip) :
  m_conn(conn), m_dofs(dofs), m_iip(iip)
{
  // mesh dimensions
  m_nelem = m_conn.shape()[0];
  m_nne   = m_conn.shape()[1];
  m_nnode = m_dofs.shape()[0];
  m_ndim  = m_dofs.shape()[1];

  // list with unknown DOFs
  m_iiu   = xt::setdiff1d(dofs, iip);

  // dimensions of the system
  m_ndof  = xt::amax(m_dofs)[0] + 1;
  m_nnp   = m_iip.size();
  m_nnu   = m_iiu.size();

  // allocate matrix and its inverse
  m_data  = xt::empty<double>({m_ndof});
  m_inv   = xt::empty<double>({m_ndof});

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode );
  assert( xt::amax(m_iip)[0] <= xt::amax(m_dofs)[0] );
  assert( m_ndof <= m_nnode * m_ndim );
}

// ---------------------------------------- index operator -----------------------------------------

inline double& MatrixDiagonal::operator[](size_t i)
{
  m_change = true;

  return m_data[i];
}

// ---------------------------------------- index operator -----------------------------------------

inline const double& MatrixDiagonal::operator[](size_t i) const
{
  return m_data[i];
}

// ---------------------------------------- index operator -----------------------------------------

inline double& MatrixDiagonal::operator()(size_t i)
{
  m_change = true;

  return m_data[i];
}

// ---------------------------------------- index operator -----------------------------------------

inline const double& MatrixDiagonal::operator()(size_t i) const
{
  return m_data[i];
}

// ---------------------------------------- index operator -----------------------------------------

inline double& MatrixDiagonal::operator()(size_t i, size_t j)
{
  assert( i == j );

  UNUSED(j);

  m_change = true;

  return m_data[i];
}

// ---------------------------------------- index operator -----------------------------------------

inline const double& MatrixDiagonal::operator()(size_t i, size_t j) const
{
  assert( i == j );

  UNUSED(j);

  return m_data[i];
}

// -------------------------------------- number of elements ---------------------------------------

inline size_t MatrixDiagonal::nelem() const
{
  return m_nelem;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t MatrixDiagonal::nne() const
{
  return m_nne;
}

// ---------------------------------------- number of nodes ----------------------------------------

inline size_t MatrixDiagonal::nnode() const
{
  return m_nnode;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t MatrixDiagonal::ndim() const
{
  return m_ndim;
}

// ---------------------------------------- number of DOFs -----------------------------------------

inline size_t MatrixDiagonal::ndof() const
{
  return m_ndof;
}

// ------------------------------------ number of unknown DOFs -------------------------------------

inline size_t MatrixDiagonal::nnu() const
{
  return m_nnu;
}

// ----------------------------------- number of prescribed DOFs -----------------------------------

inline size_t MatrixDiagonal::nnp() const
{
  return m_nnp;
}

// ------------------------------------------ return DOFs ------------------------------------------

inline xt::xtensor<size_t,2> MatrixDiagonal::dofs() const
{
  return m_dofs;
}

// -------------------------------------- return unknown DOFs --------------------------------------

inline xt::xtensor<size_t,1> MatrixDiagonal::iiu() const
{
  return m_iiu;
}

// ------------------------------------ return prescribed DOFs -------------------------------------

inline xt::xtensor<size_t,1> MatrixDiagonal::iip() const
{
  return m_iip;
}

// --------------------------------------- b_i = A_ij * x_j ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::dot(const xt::xtensor<double,1> &x) const
{
  // check input
  assert( x.size() == m_ndof );

  // compute product
  return m_data * x;
}

// --------------------------------------- b_i = A_ij * x_j ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::dot_u(
  const xt::xtensor<double,1> &x_u, const xt::xtensor<double,1> &x_p) const
{
  // suppress warning
  UNUSED(x_p);

  // check input
  assert( x_u.size() == m_nnu );
  assert( x_p.size() == m_nnp );

  // allocate output
  xt::xtensor<double,1> b_u = xt::empty<double>({m_nnu});

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    b_u(i) = m_data(m_iiu(i)) * x_u(i);

  return b_u;
}

// --------------------------------------- b_i = A_ij * x_j ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::dot_p(
  const xt::xtensor<double,1> &x_u, const xt::xtensor<double,1> &x_p) const
{
  // suppress warning
  UNUSED(x_u);

  // check input
  assert( x_u.size() == m_nnu );
  assert( x_p.size() == m_nnp );

  // allocate output
  xt::xtensor<double,1> b_p = xt::empty<double>({m_nnp});

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    b_p(i) = m_data(m_iip(i)) * x_p(i);

  return b_p;
}

// ----------------------------- assemble matrices stored per element ------------------------------

inline void MatrixDiagonal::assemble(const xt::xtensor<double,3> &elemmat)
{
  // check input
  assert( elemmat.shape()[0] == m_nelem      );
  assert( elemmat.shape()[1] == m_nne*m_ndim );
  assert( elemmat.shape()[2] == m_nne*m_ndim );
  assert( Element::isDiagonal(elemmat) );

  // zero-initialize matrix
  m_data.fill(0.0);

  // assemble
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        m_data(m_dofs(m_conn(e,m),i)) += elemmat(e,m*m_ndim+i,m*m_ndim+i);

  // signal change
  m_change = true;
}

// -------------------------------------------------------------------------------------------------

inline void MatrixDiagonal::factorize()
{
  // skip for unchanged "m_data"
  if ( ! m_change ) return;

  // invert if needed
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_ndof ; ++i )
    m_inv(i) = 1. / m_data(i);

  // reset signal
  m_change = false;
}

// --------------------------------------- solve: A * x = b ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::solve(const xt::xtensor<double,1> &b)
{
  // check input
  assert( m_nnp    == 0      );
  assert( b.size() == m_ndof );

  // compute inverse
  this->factorize();

  // solve
  xt::xtensor<double,1> x = m_inv * b;

  return x;
}

// --------------------------------------- solve: A * x = b ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::solve(
  const xt::xtensor<double,1> &b, const xt::xtensor<double,1> &x_p)
{
  // check input
  assert( x_p.size() == m_nnp  );
  assert( b  .size() == m_ndof );

  // compute inverse
  this->factorize();

  // solve
  xt::xtensor<double,1> x = m_inv * b;

  // set prescribed DOFs
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    x(m_iip(i)) = x_p(i);

  return x;
}

// --------------------------------------- solve: A * x = b ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::solve_u(
  const xt::xtensor<double,1> &b_u, const xt::xtensor<double,1> &x_p)
{
  // suppress warning
  UNUSED(x_p);

  // check input
  assert( x_p.size() == m_nnp );
  assert( b_u.size() == m_nnu );

  // compute inverse
  this->factorize();

  // allocate output
  xt::xtensor<double,1> x_u = xt::empty<double>({m_nnu});

  // solve
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    x_u(i) = m_inv(m_iiu(i)) * b_u(i);

  return x_u;
}

// ----------------------------------- return as diagonal matrix -----------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::asDiagonal() const
{
  return m_data;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

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
  // extract mesh dimensions
  m_nelem = m_conn.shape()[0];
  m_nne   = m_conn.shape()[1];
  m_nnode = m_dofs.shape()[0];
  m_ndim  = m_dofs.shape()[1];
  m_ndof  = xt::amax(m_dofs)[0] + 1;
  m_nnp   = 0;
  m_nnu   = m_ndof;
  m_iiu   = xt::arange<size_t>(m_ndof);
  m_iip   = xt::empty<size_t>({0});

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode );
  assert( m_ndof <= m_nnode * m_ndim );

  // allocate matrix and its inverse
  m_data = xt::empty<double>({m_ndof});
  m_inv  = xt::empty<double>({m_ndof});
}

// ------------------------------------------ constructor ------------------------------------------

inline MatrixDiagonal::MatrixDiagonal(const xt::xtensor<size_t,2> &conn,
  const xt::xtensor<size_t,2> &dofs, const xt::xtensor<size_t,1> &iip) :
  m_conn(conn), m_dofs(dofs), m_iip(iip)
{
  // extract mesh dimensions
  m_nelem = m_conn.shape()[0];
  m_nne   = m_conn.shape()[1];
  m_nnode = m_dofs.shape()[0];
  m_ndim  = m_dofs.shape()[1];
  m_ndof  = xt::amax(m_dofs)[0] + 1;
  m_nnp   = m_iip.size();
  m_nnu   = m_ndof - m_nnp;

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode );
  assert( m_ndof <= m_nnode * m_ndim );

  // reorder DOFs such that they can be used for partitioning; renumber such that
  //   "iiu" -> beginning
  //   "iip" -> end
  // (otherwise preserving the order)
  // this array can be used to assemble to/from partitioned arrays
  m_part = Mesh::reorder(m_dofs, m_iip, "end");

  // extract unknown DOFs
  // - allocate
  m_iiu = xt::empty<size_t>({m_nnu});
  // - set
  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) < m_nnu )
        m_iiu(m_part(n,i)) = m_dofs(n,i);

  // allocate matrix and its inverse
  m_data = xt::empty<double>({m_ndof});
  m_inv  = xt::empty<double>({m_ndof});
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

inline double& MatrixDiagonal::operator()(size_t a)
{
  m_change = true;

  return m_data[a];
}

// ---------------------------------------- index operator -----------------------------------------

inline const double& MatrixDiagonal::operator()(size_t a) const
{
  return m_data[a];
}

// ---------------------------------------- index operator -----------------------------------------

inline double& MatrixDiagonal::operator()(size_t a, size_t b)
{
  assert( a == b );

  UNUSED(b);

  m_change = true;

  return m_data[a];
}

// ---------------------------------------- index operator -----------------------------------------

inline const double& MatrixDiagonal::operator()(size_t a, size_t b) const
{
  assert( a == b );

  UNUSED(b);

  return m_data[a];
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

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::dot(const xt::xtensor<double,1> &b) const
{
  // check input
  assert( static_cast<size_t>(b.size()) == m_ndof );

  // compute product
  return m_data * b;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::dot_u(const xt::xtensor<double,1> &b) const
{
  // check input
  assert( static_cast<size_t>(b.size()) == m_ndof );

  // allocate output
  xt::xtensor<double,1> c_u = xt::empty<double>({m_nnu});

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    c_u(i) = m_data(m_iiu(i)) * b(m_iiu(i));

  return c_u;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::dot_u(const xt::xtensor<double,1> &b_u, const xt::xtensor<double,1> &b_p) const
{
  // suppress warning
  UNUSED(b_p);

  // check input
  assert( static_cast<size_t>(b_u.size()) == m_nnu );
  assert( static_cast<size_t>(b_p.size()) == m_nnp );

  // allocate output
  xt::xtensor<double,1> c_u = xt::empty<double>({m_nnu});

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    c_u(i) = m_data(m_iiu(i)) * b_u(i);

  return c_u;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::dot_p(const xt::xtensor<double,1> &b) const
{
  // check input
  assert( static_cast<size_t>(b.size()) == m_ndof );

  // allocate output
  xt::xtensor<double,1> c_p = xt::empty<double>({m_nnp});

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    c_p(i) = m_data(m_iip(i)) * b(m_iip(i));

  return c_p;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::dot_p(const xt::xtensor<double,1> &b_u, const xt::xtensor<double,1> &b_p) const
{
  // suppress warning
  UNUSED(b_u);

  // check input
  assert( static_cast<size_t>(b_u.size()) == m_nnu );
  assert( static_cast<size_t>(b_p.size()) == m_nnp );

  // allocate output
  xt::xtensor<double,1> c_p = xt::empty<double>({m_nnp});

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    c_p(i) = m_data(m_iip(i)) * b_p(i);

  return c_p;
}

// ------------------------ check structure of matrices stored per element -------------------------

inline void MatrixDiagonal::check_diagonal(const xt::xtensor<double,3> &elemmat) const
{
  // check input
  assert( elemmat.shape()[0] == m_nelem      );
  assert( elemmat.shape()[1] == m_nne*m_ndim );
  assert( elemmat.shape()[2] == m_nne*m_ndim );

  // get numerical precision
  double eps = std::numeric_limits<double>::epsilon();

  // loop over all entries
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        for ( size_t n = 0 ; n < m_nne ; ++n )
          for ( size_t j = 0 ; j < m_ndim ; ++j )
            if ( m*m_ndim+i != n*m_ndim+j )
              if ( std::abs(elemmat(e,m*m_ndim+i,n*m_ndim+j)) > eps )
                throw std::runtime_error("Element matrices are not diagonal");
}

// ----------------------------- assemble matrices stored per element ------------------------------

inline void MatrixDiagonal::assemble(const xt::xtensor<double,3> &elemmat)
{
  // check input
  assert( elemmat.shape()[0] == m_nelem      );
  assert( elemmat.shape()[1] == m_nne*m_ndim );
  assert( elemmat.shape()[2] == m_nne*m_ndim );

  // zero-initialize matrix
  m_data *= 0.0;

  // assemble
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        m_data(m_dofs(m_conn(e,m),i)) += elemmat(e,m*m_ndim+i,m*m_ndim+i);

  // signal change
  m_change = true;
}

// ------------------------------------- solve: Mat * u = rhs --------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::solve(const xt::xtensor<double,1> &rhs)
{
  // check input
  assert( m_nnp      == 0      );
  assert( rhs.size() == m_ndof );

  // invert if needed
  if ( m_change ) {
    for ( size_t i = 0 ; i < m_ndof ; ++i ) m_inv(i) = 1. / m_data(i);
  }

  // reset signal
  m_change = false;

  // solve
  xt::xtensor<double,1> u = m_inv * rhs;

  return u;
}

// ------------------------------------- solve: Mat * u = rhs --------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::solve(const xt::xtensor<double,1> &rhs, const xt::xtensor<double,1> &u_p)
{
  // check input
  assert( u_p.size() == m_nnp  );
  assert( rhs.size() == m_ndof );

  // invert if needed
  if ( m_change ) {
    for ( size_t i = 0 ; i < m_ndof ; ++i ) m_inv(i) = 1. / m_data(i);
  }

  // reset signal
  m_change = false;

  // solve
  xt::xtensor<double,1> u = m_inv * rhs;

  // set prescribed DOFs
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    u(m_iip(i)) = u_p(i);

  return u;
}

// ------------------------------------- solve: Mat * u = rhs --------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::solve_u(const xt::xtensor<double,1> &rhs_u)
{
  // check input
  assert( rhs_u.size() == m_nnu );

  // invert if needed
  if ( m_change ) {
    for ( size_t i = 0 ; i < m_ndof ; ++i ) m_inv(i) = 1. / m_data(i);
  }

  // reset signal
  m_change = false;

  // allocate output
  xt::xtensor<double,1> u_u = xt::empty<double>({m_nnu});

  // solve
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    u_u(i) = m_inv(m_iiu(i)) * rhs_u(i);

  return u_u;
}

// ------------------------------------- solve: Mat * u = rhs --------------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::solve_u(const xt::xtensor<double,1> &rhs_u, const xt::xtensor<double,1> &u_p)
{
  // suppress warning
  UNUSED(u_p);

  // check input
  assert( u_p  .size() == m_nnp );
  assert( rhs_u.size() == m_nnu );

  // invert if needed
  if ( m_change ) {
    for ( size_t i = 0 ; i < m_ndof ; ++i ) m_inv(i) = 1. / m_data(i);
  }

  // reset signal
  m_change = false;

  // allocate output
  xt::xtensor<double,1> u_u = xt::empty<double>({m_nnu});

  // solve
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    u_u(i) = m_inv(m_iiu(i)) * rhs_u(i);

  return u_u;
}

// ----------------------------------- return as diagonal matrix -----------------------------------

inline xt::xtensor<double,1> MatrixDiagonal::asDiagonal() const
{
  return m_data;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline xt::xtensor<double,1> operator* (const MatrixDiagonal &A, const xt::xtensor<double,1> &b)
{
  return A.dot(b);
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

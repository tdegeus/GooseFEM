/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIXDIAGONAL_CPP
#define GOOSEFEM_MATRIXDIAGONAL_CPP

// -------------------------------------------------------------------------------------------------

#include "MatrixDiagonal.h"

// =========================================== GooseFEM ============================================

namespace GooseFEM {

// ------------------------------------------ constructor ------------------------------------------

inline MatrixDiagonal::MatrixDiagonal(const MatS &conn, const MatS &dofs, const ColS &iip) :
m_conn(conn), m_dofs(dofs), m_iip(iip)
{
  // extract mesh dimensions
  m_nelem = static_cast<size_t>(m_conn.rows());
  m_nne   = static_cast<size_t>(m_conn.cols());
  m_nnode = static_cast<size_t>(m_dofs.rows());
  m_ndim  = static_cast<size_t>(m_dofs.cols());
  m_ndof  = static_cast<size_t>(m_dofs.maxCoeff() + 1);
  m_nnp   = static_cast<size_t>(m_iip .size());
  m_nnu   = m_ndof - m_nnp;

  // check consistency
  assert( m_conn.maxCoeff() + 1 == m_nnode );
  assert( m_ndof <= m_nnode * m_ndim );

  // reorder DOFs such that they can be used for partitioning; renumber such that
  //   "iiu" -> beginning
  //   "iip" -> end
  // (otherwise preserving the order)
  // this array can be used to assemble to/from partitioned arrays
  m_part = Mesh::reorder(m_dofs, m_iip, "end");

  // extract unknown DOFs
  // - allocate
  m_iiu.conservativeResize(m_nnu);
  // - set
  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) < m_nnu )
        m_iiu(m_part(n,i)) = m_dofs(n,i);

  // allocate matrix and its inverse
  m_data.conservativeResize(m_ndof);
  m_inv .conservativeResize(m_ndof);
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

inline ColS MatrixDiagonal::iiu() const
{
  return m_iiu;
}

// ------------------------------------ return prescribed DOFs -------------------------------------

inline ColS MatrixDiagonal::iip() const
{
  return m_iip;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline ColD MatrixDiagonal::dot(const ColD &b) const
{
  // check input
  assert( static_cast<size_t>(b.size()) == m_ndof );

  // compute product
  return m_data.cwiseProduct(b);
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline ColD MatrixDiagonal::dot_u(const ColD &b) const
{
  // check input
  assert( static_cast<size_t>(b.size()) == m_ndof );

  // allocate output
  ColD c_u(m_nnu);

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    c_u(i) = m_data(m_iiu(i)) * b(m_iiu(i));

  return c_u;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline ColD MatrixDiagonal::dot_u(const ColD &b_u, const ColD &b_p) const
{
  // suppress warning
  UNUSED(b_p);

  // check input
  assert( static_cast<size_t>(b_u.size()) == m_nnu );
  assert( static_cast<size_t>(b_p.size()) == m_nnp );

  // allocate output
  ColD c_u(m_nnu);

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    c_u(i) = m_data(m_iiu(i)) * b_u(i);

  return c_u;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline ColD MatrixDiagonal::dot_p(const ColD &b) const
{
  // check input
  assert( static_cast<size_t>(b.size()) == m_ndof );

  // allocate output
  ColD c_p(m_nnp);

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    c_p(i) = m_data(m_iip(i)) * b(m_iip(i));

  return c_p;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline ColD MatrixDiagonal::dot_p(const ColD &b_u, const ColD &b_p) const
{
  // suppress warning
  UNUSED(b_u);

  // check input
  assert( static_cast<size_t>(b_u.size()) == m_nnu );
  assert( static_cast<size_t>(b_p.size()) == m_nnp );

  // allocate output
  ColD c_p(m_nnp);

  // compute product
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    c_p(i) = m_data(m_iip(i)) * b_p(i);

  return c_p;
}

// ------------------------ check structure of matrices stored per element -------------------------

inline void MatrixDiagonal::check_diagonal(const ArrD &elemmat) const
{
  // check input
  assert( elemmat.ndim()   == 3            );
  assert( elemmat.shape(0) == m_nelem      );
  assert( elemmat.shape(1) == m_nne*m_ndim );
  assert( elemmat.shape(2) == m_nne*m_ndim );

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

inline void MatrixDiagonal::assemble(const ArrD &elemmat)
{
  // check input
  assert( elemmat.ndim()   == 3            );
  assert( elemmat.shape(0) == m_nelem      );
  assert( elemmat.shape(1) == m_nne*m_ndim );
  assert( elemmat.shape(2) == m_nne*m_ndim );

  // zero-initialize matrix
  m_data.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize matrix
    ColD t_mat = ColD::Zero(m_ndof);

    // assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          t_mat(m_dofs(m_conn(e,m),i)) += elemmat(e,m*m_ndim+i,m*m_ndim+i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      m_data += t_mat;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  // signal change
  m_change = true;
}

// ------------------------------------- solve: Mat * u = rhs --------------------------------------

inline ColD MatrixDiagonal::solve(const ColD &rhs, const ColD &u_p)
{
  // suppress warning
  UNUSED(u_p);

  // check input
  assert( static_cast<size_t>(u_p.size()) == m_nnp  );
  assert( static_cast<size_t>(rhs.size()) == m_ndof );

  // invert if needed
  if ( m_change ) m_inv = m_data.cwiseInverse();

  // reset signal
  m_change = false;

  // solve
  ColD u = m_inv.cwiseProduct(rhs);

  // set prescribed DOFs
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    u(m_iip(i)) = u_p(i);

  return u;
}

// ------------------------------------- solve: Mat * u = rhs --------------------------------------

inline ColD MatrixDiagonal::solve_u(const ColD &rhs_u, const ColD &u_p)
{
  // suppress warning
  UNUSED(u_p);

  // check input
  assert( static_cast<size_t>(u_p  .size()) == m_nnp );
  assert( static_cast<size_t>(rhs_u.size()) == m_nnu );

  // invert if needed
  if ( m_change ) m_inv = m_data.cwiseInverse();

  // reset signal
  m_change = false;

  // allocate output
  ColD u_u(m_nnu);

  // solve
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    u_u(i) = m_inv(m_iiu(i)) * rhs_u(i);

  return u_u;
}

// ----------------------------------- return as diagonal matrix -----------------------------------

inline ColD MatrixDiagonal::asDiagonal() const
{
  return m_data;
}

// --------------------------------------- c_i = A_ij * b_j ----------------------------------------

inline ColD operator* (const MatrixDiagonal &A, const ColD &b)
{
  return A.dot(b);
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

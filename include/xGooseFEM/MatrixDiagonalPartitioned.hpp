/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_MATRIXDIAGONALPARTITIONED_CPP
#define XGOOSEFEM_MATRIXDIAGONALPARTITIONED_CPP

// -------------------------------------------------------------------------------------------------

#include "MatrixDiagonalPartitioned.h"

// =================================================================================================

namespace xGooseFEM {

// -------------------------------------------------------------------------------------------------

inline MatrixDiagonalPartitioned::MatrixDiagonalPartitioned(const xt::xtensor<size_t,2> &conn,
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

  // DOFs per node, such that iiu = arange(nnu), iip = nnu + arange(nnp)
  m_part  = Mesh::reorder(m_dofs, m_iip, "end");

  // allocate matrix and its inverse
  m_data_uu = xt::empty<double>({m_nnu});
  m_data_pp = xt::empty<double>({m_nnp});
  m_inv_uu  = xt::empty<double>({m_nnu});

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode             );
  assert( xt::amax(m_iip)[0]      <= xt::amax(m_dofs)[0] );
  assert( m_ndof                  <= m_nnode * m_ndim    );
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixDiagonalPartitioned::nelem() const
{
  return m_nelem;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixDiagonalPartitioned::nne() const
{
  return m_nne;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixDiagonalPartitioned::nnode() const
{
  return m_nnode;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixDiagonalPartitioned::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixDiagonalPartitioned::ndof() const
{
  return m_ndof;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixDiagonalPartitioned::nnu() const
{
  return m_nnu;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixDiagonalPartitioned::nnp() const
{
  return m_nnp;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> MatrixDiagonalPartitioned::dofs() const
{
  return m_dofs;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> MatrixDiagonalPartitioned::iiu() const
{
  return m_iiu;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> MatrixDiagonalPartitioned::iip() const
{
  return m_iip;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> MatrixDiagonalPartitioned::dot(const xt::xtensor<double,1> &x) const
{
  // check input
  assert( x.size() == m_ndof );

  // allocate output
  xt::xtensor<double,1> b = xt::empty<double>({m_ndof});

  // compute product
  // -
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    b(m_iiu(i)) = m_data_uu(i) * x(m_iiu(i));
  // -
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    b(m_iip(i)) = m_data_pp(i) * x(m_iip(i));

  // return output
  return b;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> MatrixDiagonalPartitioned::dot_u(
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
    b_u(i) = m_data_uu(i) * x_u(i);

  return b_u;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> MatrixDiagonalPartitioned::dot_p(
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
    b_p(i) = m_data_pp(i) * x_p(i);

  return b_p;
}

// -------------------------------------------------------------------------------------------------

inline void MatrixDiagonalPartitioned::assemble(const xt::xtensor<double,3> &elemmat)
{
  // check input
  assert( elemmat.shape()[0] == m_nelem      );
  assert( elemmat.shape()[1] == m_nne*m_ndim );
  assert( elemmat.shape()[2] == m_nne*m_ndim );
  assert( Element::isDiagonal(elemmat) );

  // zero-initialize matrix
  m_data_uu.fill(0.0);
  m_data_pp.fill(0.0);

  // assemble
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    for ( size_t m = 0 ; m < m_nne ; ++m )
    {
      for ( size_t i = 0 ; i < m_ndim ; ++i )
      {
        size_t di = m_part(m_conn(e,m),i);

        if ( di < m_nnu )
          m_data_uu(di      ) += elemmat(e,m*m_ndim+i,m*m_ndim+i);
        else
          m_data_pp(di-m_nnu) += elemmat(e,m*m_ndim+i,m*m_ndim+i);
      }
    }
  }

  // signal change
  m_change = true;
}

// -------------------------------------------------------------------------------------------------

inline void MatrixDiagonalPartitioned::factorize()
{
  // skip for unchanged "m_data"
  if ( ! m_change ) return;

  // invert
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    m_inv_uu(i) = 1. / m_data_uu(i);

  // reset signal
  m_change = false;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> MatrixDiagonalPartitioned::solve(
  const xt::xtensor<double,1> &b_u, const xt::xtensor<double,1> &x_p)
{
  // suppress warning
  UNUSED(x_p);

  // check input
  assert( b_u.shape()[0] == m_nnu );
  assert( x_p.shape()[0] == m_nnp );

  // factorise (if needed)
  this->factorize();

  // allocate output
  xt::xtensor<double,1> x_u = xt::empty<double>({m_nnu});

  // solve
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    x_u(i) = m_inv_uu(i) * b_u(i);

  // return output
  return x_u;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> MatrixDiagonalPartitioned::asDiagonal() const
{
  // allocate output
  xt::xtensor<double,1> out = xt::zeros<double>({m_ndof});

  // assemble output
  // -
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i )
    out(m_iiu(i)) = m_data_uu(i);
  // -
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i )
    out(m_iip(i)) = m_data_pp(i);

  // return output
  return out;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

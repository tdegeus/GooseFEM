/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIXPARTITIONED_CPP
#define GOOSEFEM_MATRIXPARTITIONED_CPP

// -------------------------------------------------------------------------------------------------

#include "MatrixPartitioned.h"

#include <Eigen/SparseCholesky>

// =================================================================================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

inline MatrixPartitioned::MatrixPartitioned(const xt::xtensor<size_t,2> &conn,
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

  // allocate triplet list
  m_trip_uu.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_trip_up.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_trip_pu.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_trip_pp.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);

  // allocate sparse matrices
  m_data_uu.resize(m_nnu,m_nnu);
  m_data_up.resize(m_nnu,m_nnp);
  m_data_pu.resize(m_nnp,m_nnu);
  m_data_pp.resize(m_nnp,m_nnp);

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode             );
  assert( xt::amax(m_iip)[0]      <= xt::amax(m_dofs)[0] );
  assert( m_ndof                  <= m_nnode * m_ndim    );
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixPartitioned::nelem() const
{
  return m_nelem;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixPartitioned::nne() const
{
  return m_nne;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixPartitioned::nnode() const
{
  return m_nnode;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixPartitioned::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixPartitioned::ndof() const
{
  return m_ndof;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixPartitioned::nnu() const
{
  return m_nnu;
}

// -------------------------------------------------------------------------------------------------

inline size_t MatrixPartitioned::nnp() const
{
  return m_nnp;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> MatrixPartitioned::dofs() const
{
  return m_dofs;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> MatrixPartitioned::iiu() const
{
  return m_iiu;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> MatrixPartitioned::iip() const
{
  return m_iip;
}

// -------------------------------------------------------------------------------------------------

inline void MatrixPartitioned::assemble(const xt::xtensor<double,3> &elemmat)
{
  // check input
  assert( elemmat.shape()[0] == m_nelem      );
  assert( elemmat.shape()[1] == m_nne*m_ndim );
  assert( elemmat.shape()[2] == m_nne*m_ndim );

  // initialize triplet list
  m_trip_uu.clear();
  m_trip_up.clear();
  m_trip_pu.clear();
  m_trip_pp.clear();

  // assemble to triplet list
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    for ( size_t m = 0 ; m < m_nne ; ++m )
    {
      for ( size_t i = 0 ; i < m_ndim ; ++i )
      {
        size_t di = m_part(m_conn(e,m),i);

        for ( size_t n = 0 ; n < m_nne ; ++n )
        {
          for ( size_t j = 0 ; j < m_ndim ; ++j )
          {
            size_t dj = m_part(m_conn(e,n),j);

            if      ( di < m_nnu and dj < m_nnu )
              m_trip_uu.push_back(TripD(di      ,dj      ,elemmat(e,m*m_ndim+i,n*m_ndim+j)));
            else if ( di < m_nnu )
              m_trip_up.push_back(TripD(di      ,dj-m_nnu,elemmat(e,m*m_ndim+i,n*m_ndim+j)));
            else if ( dj < m_nnu )
              m_trip_pu.push_back(TripD(di-m_nnu,dj      ,elemmat(e,m*m_ndim+i,n*m_ndim+j)));
            else
              m_trip_pp.push_back(TripD(di-m_nnu,dj-m_nnu,elemmat(e,m*m_ndim+i,n*m_ndim+j)));
          }
        }
      }
    }
  }

  // convert to sparse matrix
  m_data_uu.setFromTriplets(m_trip_uu.begin(), m_trip_uu.end());
  m_data_up.setFromTriplets(m_trip_up.begin(), m_trip_up.end());
  m_data_pu.setFromTriplets(m_trip_pu.begin(), m_trip_pu.end());
  m_data_pp.setFromTriplets(m_trip_pp.begin(), m_trip_pp.end());

  // signal change
  m_change = true;
}

// -------------------------------------------------------------------------------------------------

inline void MatrixPartitioned::factorize()
{
  // skip for unchanged "m_data"
  if ( ! m_change ) return;

  // factorise
  m_solver.compute(m_data_uu);

  // reset signal
  m_change = false;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> MatrixPartitioned::solve(
  const xt::xtensor<double,1> &b_u, const xt::xtensor<double,1> &x_p)
{
  // check input
  assert( b_u.shape()[0] == m_nnu );
  assert( x_p.shape()[0] == m_nnp );

  // factorise (if needed)
  this->factorize();

  // copy of input as Eigen objects
  // - allocate
  Eigen::VectorXd B_u(m_nnu,1);
  Eigen::VectorXd X_p(m_nnp,1);
  // - copy
  std::copy(b_u.begin(), b_u.end(), B_u.data());
  std::copy(x_p.begin(), x_p.end(), X_p.data());

  // solve
  Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - m_data_up * X_p));

  // copy of output from Eigen object
  // - allocate
  xt::xtensor<double,1> x_u = xt::empty<double>({m_nnu});
  // - copy
  std::copy(X_u.data(), X_u.data()+m_nnu, x_u.begin());

  // return output
  return x_u;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

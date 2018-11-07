/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_VECTORPARTITIONED_CPP
#define XGOOSEFEM_VECTORPARTITIONED_CPP

// -------------------------------------------------------------------------------------------------

#include "VectorPartitioned.h"

// =================================================================================================

namespace xGooseFEM {

// -------------------------------------------------------------------------------------------------

inline VectorPartitioned::VectorPartitioned(const xt::xtensor<size_t,2> &conn,
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

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode             );
  assert( xt::amax(m_iip)[0]      <= xt::amax(m_dofs)[0] );
  assert( m_ndof                  <= m_nnode * m_ndim    );
}

// -------------------------------------------------------------------------------------------------

inline size_t VectorPartitioned::nelem() const
{
  return m_nelem;
}

// -------------------------------------------------------------------------------------------------

inline size_t VectorPartitioned::nne() const
{
  return m_nne;
}

// -------------------------------------------------------------------------------------------------

inline size_t VectorPartitioned::nnode() const
{
  return m_nnode;
}

// -------------------------------------------------------------------------------------------------

inline size_t VectorPartitioned::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

inline size_t VectorPartitioned::ndof() const
{
  return m_ndof;
}

// -------------------------------------------------------------------------------------------------

inline size_t VectorPartitioned::nnu() const
{
  return m_nnu;
}

// -------------------------------------------------------------------------------------------------

inline size_t VectorPartitioned::nnp() const
{
  return m_nnp;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> VectorPartitioned::dofs() const
{
  return m_dofs;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> VectorPartitioned::iiu() const
{
  return m_iiu;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> VectorPartitioned::iip() const
{
  return m_iip;
}

// ------------------- assemble: dofval[iiu] = dofval_u; dofval[iip] = dofval_p --------------------

inline void VectorPartitioned::asDofs(const xt::xtensor<double,1> &dofval_u,
  const xt::xtensor<double,1> &dofval_p, xt::xtensor<double,1> &dofval) const
{
  assert( dofval_u.size() == m_nnu  );
  assert( dofval_p.size() == m_nnp  );
  assert( dofval.size()   == m_ndof );

  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i ) dofval(m_iiu(i)) = dofval_u(i);

  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i ) dofval(m_iip(i)) = dofval_p(i);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::asDofs(const xt::xtensor<double,1> &dofval_u,
  const xt::xtensor<double,1> &dofval_p) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->asDofs(dofval_u, dofval_p, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asDofs(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,1> &dofval) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( dofval.size()      == m_ndof  );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) = nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asDofs_u(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,1> &dofval_u) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( dofval_u.size()    == m_nnu   );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) < m_nnu )
        dofval_u(m_part(n,i)) = nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asDofs_p(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,1> &dofval_p) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( dofval_p.size()    == m_nnp   );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) >= m_nnu )
        dofval_p(m_part(n,i)-m_nnu) = nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asDofs(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,1> &dofval) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( dofval.size()      == m_ndof  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        dofval(m_dofs(m_conn(e,m),i)) = elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asDofs_u(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,1> &dofval_u) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( dofval_u.size()    == m_nnu   );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) < m_nnu )
          dofval_u(m_part(m_conn(e,m),i)) = elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asDofs_p(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,1> &dofval_p) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( dofval_p.size()    == m_nnp   );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) >= m_nnu )
          dofval_p(m_part(m_conn(e,m),i)-m_nnu) = elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asNode(const xt::xtensor<double,1> &dofval,
  xt::xtensor<double,2> &nodevec) const
{
  assert( dofval.size()      == m_ndof  );
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      nodevec(n,i) = dofval(m_dofs(n,i));
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asNode(const xt::xtensor<double,1> &dofval_u,
  const xt::xtensor<double,1> &dofval_p, xt::xtensor<double,2> &nodevec) const
{
  assert( dofval_u.size()    == m_nnu   );
  assert( dofval_p.size()    == m_nnp   );
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n ) {
    for ( size_t i = 0 ; i < m_ndim ; ++i ) {
      if ( m_part(n,i) < m_nnu ) nodevec(n,i) = dofval_u(m_part(n,i)      );
      else                       nodevec(n,i) = dofval_p(m_part(n,i)-m_nnu);
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asNode(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,2> &nodevec) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        nodevec(m_conn(e,m),i) = elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asElement(const xt::xtensor<double,1> &dofval,
  xt::xtensor<double,3> &elemvec) const
{
  assert( dofval.size()      == m_ndof  );
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        elemvec(e,m,i) = dofval(m_dofs(m_conn(e,m),i));
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asElement(const xt::xtensor<double,1> &dofval_u,
  const xt::xtensor<double,1> &dofval_p, xt::xtensor<double,3> &elemvec) const
{
  assert( dofval_u.size()    == m_nnu   );
  assert( dofval_p.size()    == m_nnp   );
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e ) {
    for ( size_t m = 0 ; m < m_nne ; ++m ) {
      for ( size_t i = 0 ; i < m_ndim ; ++i ) {
        if ( m_part(m_conn(e,m),i)<m_nnu ) elemvec(e,m,i) = dofval_u(m_part(m_conn(e,m),i)      );
        else                               elemvec(e,m,i) = dofval_p(m_part(m_conn(e,m),i)-m_nnu);
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::asElement(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,3> &elemvec) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        elemvec(e,m,i) = nodevec(m_conn(e,m),i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::assembleDofs(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,1> &dofval) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( dofval.size()      == m_ndof  );

  dofval.fill(0.0);

  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) += nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::assembleDofs_u(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,1> &dofval_u) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( dofval_u.size()    == m_nnu   );

  dofval_u.fill(0.0);

  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) < m_nnu )
        dofval_u(m_part(n,i)) += nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::assembleDofs_p(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,1> &dofval_p) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( dofval_p.size()    == m_nnp   );

  dofval_p.fill(0.0);

  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) >= m_nnu )
        dofval_p(m_part(n,i)-m_nnu) += nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::assembleDofs(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,1> &dofval) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( dofval.size()      == m_ndof  );

  dofval.fill(0.0);

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
          dofval(m_dofs(m_conn(e,m),i)) += elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::assembleDofs_u(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,1> &dofval_u) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( dofval_u.size()    == m_nnu   );

  dofval_u.fill(0.0);

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) < m_nnu )
          dofval_u(m_dofs(m_conn(e,m),i)) += elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::assembleDofs_p(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,1> &dofval_p) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( dofval_p.size()    == m_nnp   );

  dofval_p.fill(0.0);

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) >= m_nnu )
          dofval_p(m_dofs(m_conn(e,m),i)-m_nnu) += elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline void VectorPartitioned::assembleNode(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,2> &nodevec) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );

  nodevec.fill(0.0);

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        nodevec(m_conn(e,m),i) += elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::asDofs(const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->asDofs(nodevec, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::asDofs_u(const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,1> dofval_u = xt::empty<double>({m_nnu});

  this->asDofs_u(nodevec, dofval_u);

  return dofval_u;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::asDofs_p(const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,1> dofval_p = xt::empty<double>({m_nnp});

  this->asDofs_p(nodevec, dofval_p);

  return dofval_p;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::asDofs(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->asDofs(elemvec, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::asDofs_u(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,1> dofval_u = xt::empty<double>({m_nnu});

  this->asDofs_u(elemvec, dofval_u);

  return dofval_u;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::asDofs_p(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,1> dofval_p = xt::empty<double>({m_nnp});

  this->asDofs_p(elemvec, dofval_p);

  return dofval_p;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> VectorPartitioned::asNode(const xt::xtensor<double,1> &dofval) const
{
  xt::xtensor<double,2> nodevec = xt::empty<double>({m_nnode, m_ndim});

  this->asNode(dofval, nodevec);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> VectorPartitioned::asNode(const xt::xtensor<double,1> &dofval_u,
  const xt::xtensor<double,1> &dofval_p) const
{
  xt::xtensor<double,2> nodevec = xt::empty<double>({m_nnode, m_ndim});

  this->asNode(dofval_u, dofval_p, nodevec);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> VectorPartitioned::asNode(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,2> nodevec = xt::empty<double>({m_nnode, m_ndim});

  this->asNode(elemvec, nodevec);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> VectorPartitioned::asElement(const xt::xtensor<double,1> &dofval) const
{
  xt::xtensor<double,3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});

  this->asElement(dofval, elemvec);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> VectorPartitioned::asElement(const xt::xtensor<double,1> &dofval_u,
  const xt::xtensor<double,1> &dofval_p) const
{
  xt::xtensor<double,3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});

  this->asElement(dofval_u, dofval_p, elemvec);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> VectorPartitioned::asElement(const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});

  this->asElement(nodevec, elemvec);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::assembleDofs(
  const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->assembleDofs(nodevec, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::assembleDofs_u(
  const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,1> dofval_u = xt::empty<double>({m_nnu});

  this->assembleDofs_u(nodevec, dofval_u);

  return dofval_u;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::assembleDofs_p(
  const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,1> dofval_p = xt::empty<double>({m_nnp});

  this->assembleDofs_p(nodevec, dofval_p);

  return dofval_p;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::assembleDofs(
  const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->assembleDofs(elemvec, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::assembleDofs_u(
  const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,1> dofval_u = xt::empty<double>({m_nnu});

  this->assembleDofs_u(elemvec, dofval_u);

  return dofval_u;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> VectorPartitioned::assembleDofs_p(
  const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,1> dofval_p = xt::empty<double>({m_nnp});

  this->assembleDofs_p(elemvec, dofval_p);

  return dofval_p;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> VectorPartitioned::assembleNode(
  const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,2> nodevec = xt::empty<double>({m_nnode, m_ndim});

  this->assembleNode(elemvec, nodevec);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

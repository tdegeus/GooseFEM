/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_VECTOR_HPP
#define GOOSEFEM_VECTOR_HPP

// -------------------------------------------------------------------------------------------------

#include "Vector.h"

// =================================================================================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

inline Vector::Vector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs) :
  m_conn(conn), m_dofs(dofs)
{
  // mesh dimensions
  m_nelem = m_conn.shape()[0];
  m_nne   = m_conn.shape()[1];
  m_nnode = m_dofs.shape()[0];
  m_ndim  = m_dofs.shape()[1];

  // dimensions of the system
  m_ndof  = xt::amax(m_dofs)[0] + 1;

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode          );
  assert( m_ndof                  <= m_nnode * m_ndim );
}

// -------------------------------------------------------------------------------------------------

inline size_t Vector::nelem() const { return m_nelem; }

inline size_t Vector::nne() const { return m_nne; }

inline size_t Vector::nnode() const { return m_nnode; }

inline size_t Vector::ndim() const { return m_ndim; }

inline size_t Vector::ndof() const { return m_ndof; }

inline xt::xtensor<size_t,2> Vector::dofs() const { return m_dofs; }

// -------------------------------------------------------------------------------------------------

inline void Vector::copy(const xt::xtensor<double,2> &nodevec_src,
    xt::xtensor<double,2> &nodevec_dest) const
{
  assert( nodevec_src .shape()[0] == m_nnode );
  assert( nodevec_src .shape()[1] == m_ndim  );
  assert( nodevec_dest.shape()[0] == m_nnode );
  assert( nodevec_dest.shape()[1] == m_ndim  );

  xt::noalias(nodevec_dest) = nodevec_src;
}

// -------------------------------------------------------------------------------------------------

inline void Vector::asDofs(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,1> &dofval) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( dofval.size()      == m_ndof  );

  #pragma omp parallel for
  for ( size_t m = 0 ; m < m_nnode ; ++m )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(m,i)) = nodevec(m,i);
}

// -------------------------------------------------------------------------------------------------

inline void Vector::asDofs(const xt::xtensor<double,3> &elemvec,
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

inline void Vector::asNode(const xt::xtensor<double,1> &dofval,
  xt::xtensor<double,2> &nodevec) const
{
  assert( dofval.size()      == m_ndof  );
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );

  #pragma omp parallel for
  for ( size_t m = 0 ; m < m_nnode ; ++m )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      nodevec(m,i) = dofval(m_dofs(m,i));
}

// -------------------------------------------------------------------------------------------------

inline void Vector::asNode(const xt::xtensor<double,3> &elemvec,
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

inline void Vector::asElement(const xt::xtensor<double,1> &dofval,
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

inline void Vector::asElement(const xt::xtensor<double,2> &nodevec,
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

inline void Vector::assembleDofs(const xt::xtensor<double,2> &nodevec,
  xt::xtensor<double,1> &dofval) const
{
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );
  assert( dofval.size()      == m_ndof  );

  dofval.fill(0.0);

  for ( size_t m = 0 ; m < m_nnode ; ++m )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(m,i)) += nodevec(m,i);
}

// -------------------------------------------------------------------------------------------------

inline void Vector::assembleDofs(const xt::xtensor<double,3> &elemvec,
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

inline void Vector::assembleNode(const xt::xtensor<double,3> &elemvec,
  xt::xtensor<double,2> &nodevec) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( nodevec.shape()[0] == m_nnode );
  assert( nodevec.shape()[1] == m_ndim  );

  // assemble to DOFs
  xt::xtensor<double,1> dofval = this->assembleDofs(elemvec);

  // read from DOFs
  this->asNode(dofval, nodevec);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Vector::asDofs(const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->asDofs(nodevec, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Vector::asDofs(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->asDofs(elemvec, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Vector::asNode(const xt::xtensor<double,1> &dofval) const
{
  xt::xtensor<double,2> nodevec = xt::empty<double>({m_nnode, m_ndim});

  this->asNode(dofval, nodevec);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Vector::asNode(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,2> nodevec = xt::empty<double>({m_nnode, m_ndim});

  this->asNode(elemvec, nodevec);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Vector::asElement(const xt::xtensor<double,1> &dofval) const
{
  xt::xtensor<double,3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});

  this->asElement(dofval, elemvec);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Vector::asElement(const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});

  this->asElement(nodevec, elemvec);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Vector::assembleDofs(
  const xt::xtensor<double,2> &nodevec) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->assembleDofs(nodevec, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Vector::assembleDofs(
  const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,1> dofval = xt::empty<double>({m_ndof});

  this->assembleDofs(elemvec, dofval);

  return dofval;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Vector::assembleNode(
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

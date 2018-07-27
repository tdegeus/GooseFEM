/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_VECTOR_CPP
#define GOOSEFEM_VECTOR_CPP

// -------------------------------------------------------------------------------------------------

#include "Vector.h"

// =========================================== GooseFEM ============================================

namespace GooseFEM {

// ------------------------------------------ constructor ------------------------------------------

inline Vector::Vector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs) :
m_conn(conn), m_dofs(dofs)
{
  // extract mesh dimensions
  m_nelem = static_cast<size_t>(m_conn.shape()[0]);
  m_nne   = static_cast<size_t>(m_conn.shape()[1]);
  m_nnode = static_cast<size_t>(m_dofs.shape()[0]);
  m_ndim  = static_cast<size_t>(m_dofs.shape()[1]);
  m_ndof  = static_cast<size_t>(xt::amax(m_dofs)[0] + 1);
  m_nnp   = 0;
  m_nnu   = m_ndof;
  m_iiu   = xt::arange<size_t>(m_ndof);
  m_iip   = xt::empty<size_t>({0});

  // check consistency
  assert( xt::amax(m_conn)[0] + 1 == m_nnode );
  assert( m_ndof <= m_nnode * m_ndim );
}

// ------------------------------------------ constructor ------------------------------------------

inline Vector::Vector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs, const xt::xtensor<size_t,1> &iip) :
m_conn(conn), m_dofs(dofs), m_iip(iip)
{
  // extract mesh dimensions
  m_nelem = static_cast<size_t>(m_conn.shape()[0]);
  m_nne   = static_cast<size_t>(m_conn.shape()[1]);
  m_nnode = static_cast<size_t>(m_dofs.shape()[0]);
  m_ndim  = static_cast<size_t>(m_dofs.shape()[1]);
  m_ndof  = static_cast<size_t>(xt::amax(m_dofs)[0] + 1);
  m_nnp   = static_cast<size_t>(m_iip .size());
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
}

// -------------------------------------- number of elements ---------------------------------------

inline size_t Vector::nelem() const
{
  return m_nelem;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t Vector::nne() const
{
  return m_nne;
}

// ---------------------------------------- number of nodes ----------------------------------------

inline size_t Vector::nnode() const
{
  return m_nnode;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t Vector::ndim() const
{
  return m_ndim;
}

// ---------------------------------------- number of DOFs -----------------------------------------

inline size_t Vector::ndof() const
{
  return m_ndof;
}

// ------------------------------------ number of unknown DOFs -------------------------------------

inline size_t Vector::nnu() const
{
  return m_nnu;
}

// ----------------------------------- number of prescribed DOFs -----------------------------------

inline size_t Vector::nnp() const
{
  return m_nnp;
}

// ------------------------------------------ return DOFs ------------------------------------------

inline xt::xtensor<size_t,2> Vector::dofs() const
{
  return m_dofs;
}

// -------------------------------------- return unknown DOFs --------------------------------------

inline xt::xtensor<size_t,1> Vector::iiu() const
{
  return m_iiu;
}

// ------------------------------------ return prescribed DOFs -------------------------------------

inline xt::xtensor<size_t,1> Vector::iip() const
{
  return m_iip;
}

// --------------------------------------- dofval -> dofval ----------------------------------------

inline xt::xtensor<double,1> Vector::asDofs(const xt::xtensor<double,1> &dofval_u, const xt::xtensor<double,1> &dofval_p) const
{
  // check input
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_ndof});

  // apply conversion
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i ) dofval(m_iiu(i)) = dofval_u(i);
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i ) dofval(m_iip(i)) = dofval_p(i);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::asDofs(const xt::xtensor<double,2> &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.shape()[0]) == m_nnode );
  assert( static_cast<size_t>(nodevec.shape()[1]) == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_ndof});

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) = nodevec(n,i);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::asDofs_u(const xt::xtensor<double,2> &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.shape()[0]) == m_nnode );
  assert( static_cast<size_t>(nodevec.shape()[1]) == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_nnu});

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) < m_nnu )
        dofval(m_part(n,i)) = nodevec(n,i);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::asDofs_p(const xt::xtensor<double,2> &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.shape()[0]) == m_nnode );
  assert( static_cast<size_t>(nodevec.shape()[1]) == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_nnp});

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) >= m_nnu )
        dofval(m_part(n,i)-m_nnu) = nodevec(n,i);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::asDofs(const xt::xtensor<double,3> &elemvec) const
{
  // check input
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_ndof});

  // apply conversion
  #pragma omp for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        dofval(m_dofs(m_conn(e,m),i)) = elemvec(e,m,i);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::asDofs_u(const xt::xtensor<double,3> &elemvec) const
{
  // check input
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_nnu});

  // apply conversion
  #pragma omp for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) < m_nnu )
          dofval(m_part(m_conn(e,m),i)) = elemvec(e,m,i);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::asDofs_p(const xt::xtensor<double,3> &elemvec) const
{
  // check input
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_nnp});

  // apply conversion
  #pragma omp for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) >= m_nnu )
          dofval(m_part(m_conn(e,m),i)-m_nnu) = elemvec(e,m,i);

  return dofval;
}

// --------------------------------------- dofval -> nodevec ---------------------------------------

inline xt::xtensor<double,2> Vector::asNode(const xt::xtensor<double,1> &dofval) const
{
  // check input
  assert( static_cast<size_t>(dofval.size()) == m_ndof );

  // zero-initialize output
  xt::xtensor<double,2> nodevec = xt::zeros<double>({m_nnode, m_ndim});

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      nodevec(n,i) = dofval(m_dofs(n,i));

  return nodevec;
}

// --------------------------------------- dofval -> nodevec ---------------------------------------

inline xt::xtensor<double,2> Vector::asNode(const xt::xtensor<double,1> &dofval_u, const xt::xtensor<double,1> &dofval_p) const
{
  // check input
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp );

  // zero-initialize output
  xt::xtensor<double,2> nodevec = xt::zeros<double>({m_nnode, m_ndim});

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n ) {
    for ( size_t i = 0 ; i < m_ndim ; ++i ) {
      if ( m_part(n,i) < m_nnu ) nodevec(n,i) = dofval_u(m_part(n,i)      );
      else                       nodevec(n,i) = dofval_p(m_part(n,i)-m_nnu);
    }
  }

  return nodevec;
}

// --------------------------------------- elemvec -> nodevec ---------------------------------------

inline xt::xtensor<double,2> Vector::asNode(const xt::xtensor<double,3> &elemvec) const
{
  // check input
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,2> nodevec = xt::zeros<double>({m_nnode, m_ndim});

  // apply conversion
  #pragma omp for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        nodevec(m_conn(e,m),i) = elemvec(e,m,i);

  return nodevec;
}

// --------------------------------------- dofval -> elemvec ---------------------------------------

inline xt::xtensor<double,3> Vector::asElement(const xt::xtensor<double,1> &dofval) const
{
  // check input
  assert( static_cast<size_t>(dofval.size()) == m_ndof );

  // zero-initialize output: nodal vectors stored per element
  xt::xtensor<double,3> elemvec = xt::zeros<double>({m_nelem, m_nne, m_ndim});

  // read from nodal vectors
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        elemvec(e,m,i) = dofval(m_dofs(m_conn(e,m),i));

  return elemvec;
}

// --------------------------------------- dofval -> elemvec ---------------------------------------

inline xt::xtensor<double,3> Vector::asElement(const xt::xtensor<double,1> &dofval_u, const xt::xtensor<double,1> &dofval_p) const
{
  // check input
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp );

  // zero-initialize output: nodal vectors stored per element
  xt::xtensor<double,3> elemvec = xt::zeros<double>({m_nelem, m_nne, m_ndim});

  // read from nodal vectors
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e ) {
    for ( size_t m = 0 ; m < m_nne ; ++m ) {
      for ( size_t i = 0 ; i < m_ndim ; ++i ) {
        if ( m_part(m_conn(e,m),i)<m_nnu ) elemvec(e,m,i) = dofval_u(m_part(m_conn(e,m),i)      );
        else                               elemvec(e,m,i) = dofval_p(m_part(m_conn(e,m),i)-m_nnu);
      }
    }
  }

  return elemvec;
}

// -------------------------------------- nodevec -> elemvec ---------------------------------------

inline xt::xtensor<double,3> Vector::asElement(const xt::xtensor<double,2> &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.shape()[0]) == m_nnode );
  assert( static_cast<size_t>(nodevec.shape()[1]) == m_ndim  );

  // zero-initialize output: nodal vectors stored per element
  xt::xtensor<double,3> elemvec = xt::zeros<double>({m_nelem, m_nne, m_ndim});

  // read from nodal vectors
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        elemvec(e,m,i) = nodevec(m_conn(e,m),i);

  return elemvec;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::assembleDofs(const xt::xtensor<double,2> &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.shape()[0]) == m_nnode );
  assert( static_cast<size_t>(nodevec.shape()[1]) == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_ndof});

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output
    xt::xtensor<double,1> t_dofval = xt::zeros<double>({m_ndof});

    // assemble
    #pragma omp for
    for ( size_t n = 0 ; n < m_nnode ; ++n )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        t_dofval(m_dofs(n,i)) += nodevec(n,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::assembleDofs_u(const xt::xtensor<double,2> &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.shape()[0]) == m_nnode );
  assert( static_cast<size_t>(nodevec.shape()[1]) == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_nnu});

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output
    xt::xtensor<double,1> t_dofval = xt::zeros<double>({m_nnu});

    // assemble
    #pragma omp for
    for ( size_t n = 0 ; n < m_nnode ; ++n )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(n,i) < m_nnu )
          t_dofval(m_part(n,i)) += nodevec(n,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::assembleDofs_p(const xt::xtensor<double,2> &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.shape()[0]) == m_nnode );
  assert( static_cast<size_t>(nodevec.shape()[1]) == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_nnp});

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output
    xt::xtensor<double,1> t_dofval = xt::zeros<double>({m_nnp});

    // assemble
    #pragma omp for
    for ( size_t n = 0 ; n < m_nnode ; ++n )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(n,i) >= m_nnu )
          t_dofval(m_part(n,i)-m_nnu) += nodevec(n,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::assembleDofs(const xt::xtensor<double,3> &elemvec) const
{
  // check input
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_ndof});

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output
    xt::xtensor<double,1> t_dofval = xt::zeros<double>({m_ndof});

    // assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          t_dofval(m_dofs(m_conn(e,m),i)) += elemvec(e,m,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::assembleDofs_u(const xt::xtensor<double,3> &elemvec) const
{
  // check input
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_nnu});

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output
    xt::xtensor<double,1> t_dofval = xt::zeros<double>({m_nnu});

    // assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          if ( m_part(m_conn(e,m),i) < m_nnu )
            t_dofval(m_dofs(m_conn(e,m),i)) += elemvec(e,m,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline xt::xtensor<double,1> Vector::assembleDofs_p(const xt::xtensor<double,3> &elemvec) const
{
  // check input
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,1> dofval = xt::zeros<double>({m_nnp});

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output
    xt::xtensor<double,1> t_dofval = xt::zeros<double>({m_nnp});

    // assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          if ( m_part(m_conn(e,m),i) >= m_nnu )
            t_dofval(m_dofs(m_conn(e,m),i)-m_nnu) += elemvec(e,m,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  return dofval;
}

// --------------------------------------- elemvec -> nodevec ---------------------------------------

inline xt::xtensor<double,2> Vector::assembleNode(const xt::xtensor<double,3> &elemvec) const
{
  // check input
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output
  xt::xtensor<double,2> nodevec = xt::zeros<double>({m_nnode, m_ndim});

  // start threads (all variables declared in this scope are local to each thread)
  #pragma omp parallel
  {
    // zero-initialize output
    xt::xtensor<double,2> t_nodevec = xt::zeros<double>({m_nnode, m_ndim});

    // assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          t_nodevec(m_conn(e,m),i) += elemvec(e,m,i);

    // reduce: combine result obtained on the different threads
    #pragma omp critical
      nodevec += t_nodevec;
  }

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

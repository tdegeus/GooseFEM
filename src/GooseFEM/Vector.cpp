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

inline Vector::Vector(const MatS &conn, const MatS &dofs, const ColS &iip) :
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
  // TODO: make more complete: it is also assumed that DOFs and iiu/iip have no missing numbers
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

// -------------------------------------- return unknown DOFs --------------------------------------

inline ColS Vector::iiu() const
{
  return m_iiu;
}

// ------------------------------------ return prescribed DOFs -------------------------------------

inline ColS Vector::iip() const
{
  return m_iip;
}

// --------------------------------------- dofval -> dofval ----------------------------------------

inline ColD Vector::asDofs(const ColD &dofval_u, const ColD &dofval_p) const
{
  // check input
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp );

  // allocate output
  ColD dofval(m_ndof);

  // apply conversion
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i ) dofval(m_iiu(i)) = dofval_u(i);
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i ) dofval(m_iip(i)) = dofval_p(i);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline ColD Vector::asDofs(const MatD &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  // allocate output
  ColD dofval(m_ndof);

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) = nodevec(n,i);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline ColD Vector::asDofs_u(const MatD &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  // allocate output
  ColD dofval(m_nnu);

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) < m_nnu )
        dofval(m_part(n,i)) = nodevec(n,i);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline ColD Vector::asDofs_p(const MatD &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  // allocate output
  ColD dofval(m_nnp);

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) >= m_nnu )
        dofval(m_part(n,i)-m_nnu) = nodevec(n,i);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline ColD Vector::asDofs(const ArrD &elemvec) const
{
  // check input
  assert( elemvec.ndim()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // allocate output
  ColD dofval(m_ndof);

  // apply conversion
  #pragma omp for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        dofval(m_dofs(m_conn(e,m),i)) = elemvec(e,m,i);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline ColD Vector::asDofs_u(const ArrD &elemvec) const
{
  // check input
  assert( elemvec.ndim()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // allocate output
  ColD dofval(m_nnu);

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

inline ColD Vector::asDofs_p(const ArrD &elemvec) const
{
  // check input
  assert( elemvec.ndim()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // allocate output
  ColD dofval(m_nnp);

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

inline MatD Vector::asNode(const ColD &dofval) const
{
  // check input
  assert( static_cast<size_t>(dofval.size()) == m_ndof );

  // allocate output
  MatD nodevec(m_nnode, m_ndim);

  // apply conversion
  #pragma omp for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      nodevec(n,i) = dofval(m_dofs(n,i));

  return nodevec;
}

// --------------------------------------- dofval -> nodevec ---------------------------------------

inline MatD Vector::asNode(const ColD &dofval_u, const ColD &dofval_p) const
{
  // check input
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp );

  // allocate output
  MatD nodevec(m_nnode, m_ndim);

  // apply conversion
  #pragma omp for
  {
    for ( size_t n = 0 ; n < m_nnode ; ++n ) {
      for ( size_t i = 0 ; i < m_ndim ; ++i ) {
        if ( m_part(n,i) < m_nnu ) nodevec(n,i) = dofval_u(m_part(n,i)      );
        else                       nodevec(n,i) = dofval_p(m_part(n,i)-m_nnu);
      }
    }
  }

  return nodevec;
}

// --------------------------------------- elemvec -> nodevec ---------------------------------------

inline MatD Vector::asNode(const ArrD &elemvec) const
{
  // check input
  assert( elemvec.ndim()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // allocate output
  MatD nodevec(m_nnode, m_ndim);

  // apply conversion
  #pragma omp for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        nodevec(m_conn(e,m),i) = elemvec(e,m,i);

  return nodevec;
}

// --------------------------------------- dofval -> elemvec ---------------------------------------

inline ArrD Vector::asElement(const ColD &dofval) const
{
  // check input
  assert( static_cast<size_t>(dofval.size()) == m_ndof );

  // allocate output: nodal vectors stored per element
  ArrD elemvec({m_nelem, m_nne, m_ndim});

  // read from nodal vectors
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        elemvec(e,m,i) = dofval(m_dofs(m_conn(e,m),i));

  return elemvec;
}

// --------------------------------------- dofval -> elemvec ---------------------------------------

inline ArrD Vector::asElement(const ColD &dofval_u, const ColD &dofval_p) const
{
  // check input
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp );

  // allocate output: nodal vectors stored per element
  ArrD elemvec({m_nelem, m_nne, m_ndim});

  // read from nodal vectors
  #pragma omp parallel for
  {
    for ( size_t e = 0 ; e < m_nelem ; ++e ) {
      for ( size_t m = 0 ; m < m_nne ; ++m ) {
        for ( size_t i = 0 ; i < m_ndim ; ++i ) {
          if ( m_part(m_conn(e,m),i)<m_nnu ) elemvec(e,m,i) = dofval_u(m_part(m_conn(e,m),i)      );
          else                               elemvec(e,m,i) = dofval_p(m_part(m_conn(e,m),i)-m_nnu);
        }
      }
    }
  }

  return elemvec;
}

// -------------------------------------- nodevec -> elemvec ---------------------------------------

inline ArrD Vector::asElement(const MatD &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  // allocate output: nodal vectors stored per element
  ArrD elemvec({m_nelem, m_nne, m_ndim});

  // read from nodal vectors
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        elemvec(e,m,i) = nodevec(m_conn(e,m),i);

  return elemvec;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline ColD Vector::assembleDofs(const MatD &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  // allocate output
  ColD dofval(m_ndof);

  // zero-initialize output
  dofval.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; allocate output
    ColD t_dofval(m_ndof);

    // - per thread; zero-initialize output
    t_dofval.setZero();

    // - per thread; assemble
    #pragma omp for
    for ( size_t n = 0 ; n < m_nnode ; ++n )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        t_dofval(m_dofs(n,i)) += nodevec(n,i);

    // - reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline ColD Vector::assembleDofs_u(const MatD &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  // allocate output
  ColD dofval(m_nnu);

  // zero-initialize output
  dofval.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; allocate output
    ColD t_dofval(m_nnu);

    // - per thread; zero-initialize output
    t_dofval.setZero();

    // - per thread; assemble
    #pragma omp for
    for ( size_t n = 0 ; n < m_nnode ; ++n )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(n,i) < m_nnu )
          t_dofval(m_part(n,i)) += nodevec(n,i);

    // - reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline ColD Vector::assembleDofs_p(const MatD &nodevec) const
{
  // check input
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  // allocate output
  ColD dofval(m_nnp);

  // zero-initialize output
  dofval.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; allocate output
    ColD t_dofval(m_nnp);

    // - per thread; zero-initialize output
    t_dofval.setZero();

    // - per thread; assemble
    #pragma omp for
    for ( size_t n = 0 ; n < m_nnode ; ++n )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(n,i) >= m_nnu )
          t_dofval(m_part(n,i)-m_nnu) += nodevec(n,i);

    // - reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline ColD Vector::assembleDofs(const ArrD &elemvec) const
{
  // check input
  assert( elemvec.ndim()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // allocate output
  ColD dofval(m_ndof);

  // zero-initialize output
  dofval.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; allocate output
    ColD t_dofval(m_ndof);

    // - per thread; zero-initialize output
    t_dofval.setZero();

    // - per thread; assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          t_dofval(m_dofs(m_conn(e,m),i)) += elemvec(e,m,i);

    // - reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline ColD Vector::assembleDofs_u(const ArrD &elemvec) const
{
  // check input
  assert( elemvec.ndim()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // allocate output
  ColD dofval(m_nnu);

  // zero-initialize output
  dofval.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; allocate output
    ColD t_dofval(m_nnu);

    // - per thread; zero-initialize output
    t_dofval.setZero();

    // - per thread; assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          if ( m_part(m_conn(e,m),i) < m_nnu )
            t_dofval(m_dofs(m_conn(e,m),i)) += elemvec(e,m,i);

    // - reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline ColD Vector::assembleDofs_p(const ArrD &elemvec) const
{
  // check input
  assert( elemvec.ndim()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // allocate output
  ColD dofval(m_nnp);

  // zero-initialize output
  dofval.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; allocate output
    ColD t_dofval(m_nnp);

    // - per thread; zero-initialize output
    t_dofval.setZero();

    // - per thread; assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          if ( m_part(m_conn(e,m),i) >= m_nnu )
            t_dofval(m_dofs(m_conn(e,m),i)-m_nnu) += elemvec(e,m,i);

    // - reduce: combine result obtained on the different threads
    #pragma omp critical
      dofval += t_dofval;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return dofval;
}

// --------------------------------------- elemvec -> nodevec ---------------------------------------

inline MatD Vector::assembleNode(const ArrD &elemvec) const
{
  // check input
  assert( elemvec.ndim()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // allocate output
  MatD nodevec(m_nnode, m_ndim);

  // zero-initialize output
  nodevec.setZero();

  // temporarily disable parallelization by Eigen
  Eigen::setNbThreads(1);

  // start threads
  #pragma omp parallel
  {
    // - per thread; allocate output
    MatD t_nodevec(m_nnode, m_ndim);

    // - per thread; zero-initialize output
    t_nodevec.setZero();

    // - per thread; assemble
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          t_nodevec(m_conn(e,m),i) += elemvec(e,m,i);

    // - reduce: combine result obtained on the different threads
    #pragma omp critical
      nodevec += t_nodevec;
  }

  // reset automatic parallelization by Eigen
  Eigen::setNbThreads(0);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

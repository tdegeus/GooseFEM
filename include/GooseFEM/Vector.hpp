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

inline Vector::Vector(const MatS &conn, const MatS &dofs) :
  m_conn(conn), m_dofs(dofs)
{
  m_nelem = static_cast<size_t>(m_conn.rows());
  m_nne   = static_cast<size_t>(m_conn.cols());
  m_nnode = static_cast<size_t>(m_dofs.rows());
  m_ndim  = static_cast<size_t>(m_dofs.cols());
  m_ndof  = static_cast<size_t>(m_dofs.maxCoeff() + 1);
  m_nnp   = 0;
  m_nnu   = m_ndof - m_nnp;
  m_iiu.conservativeResize(m_ndof);
  m_iip.conservativeResize(0);

  for ( size_t i = 0 ; i < m_ndof ; ++i ) m_iiu(i) = i;

  // check consistency
  assert( m_conn.maxCoeff() + 1 == m_nnode );
  assert( m_ndof <= m_nnode * m_ndim );
}

// ------------------------------------------ constructor ------------------------------------------

inline Vector::Vector(const MatS &conn, const MatS &dofs, const ColS &iip) :
  m_conn(conn), m_dofs(dofs), m_iip(iip)
{
  m_nelem = static_cast<size_t>(m_conn.rows());
  m_nne   = static_cast<size_t>(m_conn.cols());
  m_nnode = static_cast<size_t>(m_dofs.rows());
  m_ndim  = static_cast<size_t>(m_dofs.cols());
  m_ndof  = static_cast<size_t>(m_dofs.maxCoeff() + 1);
  m_nnp   = static_cast<size_t>(m_iip .size());
  m_nnu   = m_ndof - m_nnp;
  m_iiu.conservativeResize(m_nnu);

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

inline MatS Vector::dofs() const
{
  return m_dofs;
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

inline void Vector::asDofs(const ColD &dofval_u, const ColD &dofval_p, ColD &dofval) const
{
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu  );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp  );
  assert( static_cast<size_t>(dofval.size())   == m_ndof );

  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnu ; ++i ) dofval(m_iiu(i)) = dofval_u(i);
  #pragma omp parallel for
  for ( size_t i = 0 ; i < m_nnp ; ++i ) dofval(m_iip(i)) = dofval_p(i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::asDofs(const ColD &dofval_u, const ColD &dofval_p) const
{
  ColD dofval(m_ndof);

  this->asDofs(dofval_u, dofval_p, dofval);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline void Vector::asDofs(const MatD &nodevec, ColD &dofval) const
{
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );
  assert( static_cast<size_t>(dofval.size())  == m_ndof  );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) = nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::asDofs(const MatD &nodevec) const
{
  ColD dofval(m_ndof);

  this->asDofs(nodevec, dofval);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline void Vector::asDofs_u(const MatD &nodevec, ColD &dofval) const
{
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );
  assert( static_cast<size_t>(dofval.size())  == m_nnu   );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) < m_nnu )
        dofval(m_part(n,i)) = nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::asDofs_u(const MatD &nodevec) const
{
  ColD dofval(m_nnu);

  this->asDofs_u(nodevec, dofval);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline void Vector::asDofs_p(const MatD &nodevec, ColD &dofval) const
{
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );
  assert( static_cast<size_t>(dofval.size())  == m_nnp   );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) >= m_nnu )
        dofval(m_part(n,i)-m_nnu) = nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::asDofs_p(const MatD &nodevec) const
{
  ColD dofval(m_nnp);

  this->asDofs_p(nodevec, dofval);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline void Vector::asDofs(const ArrD &elemvec, ColD &dofval) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );
  assert( static_cast<size_t>(dofval.size()) == m_ndof );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        dofval(m_dofs(m_conn(e,m),i)) = elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::asDofs(const ArrD &elemvec) const
{
  ColD dofval(m_ndof);

  this->asDofs(elemvec, dofval);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline void Vector::asDofs_u(const ArrD &elemvec, ColD &dofval) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );
  assert( static_cast<size_t>(dofval.size()) == m_nnu );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) < m_nnu )
          dofval(m_part(m_conn(e,m),i)) = elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::asDofs_u(const ArrD &elemvec) const
{
  ColD dofval(m_nnu);

  this->asDofs_u(elemvec, dofval);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline void Vector::asDofs_p(const ArrD &elemvec, ColD &dofval) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );
  assert( static_cast<size_t>(dofval.size()) == m_nnp );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) >= m_nnu )
          dofval(m_part(m_conn(e,m),i)-m_nnu) = elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::asDofs_p(const ArrD &elemvec) const
{
  ColD dofval(m_nnp);

  this->asDofs_p(elemvec, dofval);

  return dofval;
}

// --------------------------------------- dofval -> nodevec ---------------------------------------

inline void Vector::asNode(const ColD &dofval, MatD &nodevec) const
{
  assert( static_cast<size_t>(dofval.size())  == m_ndof  );
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      nodevec(n,i) = dofval(m_dofs(n,i));
}

// -------------------------------------------------------------------------------------------------

inline MatD Vector::asNode(const ColD &dofval) const
{
  MatD nodevec(m_nnode, m_ndim);

  this->asNode(dofval, nodevec);

  return nodevec;
}

// --------------------------------------- dofval -> nodevec ---------------------------------------

inline void Vector::asNode(const ColD &dofval_u, const ColD &dofval_p, MatD &nodevec) const
{
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu  );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp  );
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  #pragma omp parallel for
  for ( size_t n = 0 ; n < m_nnode ; ++n ) {
    for ( size_t i = 0 ; i < m_ndim ; ++i ) {
      if ( m_part(n,i) < m_nnu ) nodevec(n,i) = dofval_u(m_part(n,i)      );
      else                       nodevec(n,i) = dofval_p(m_part(n,i)-m_nnu);
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline MatD Vector::asNode(const ColD &dofval_u, const ColD &dofval_p) const
{
  MatD nodevec(m_nnode, m_ndim);

  this->asNode(dofval_u, dofval_p, nodevec);

  return nodevec;
}

// --------------------------------------- elemvec -> nodevec ---------------------------------------

inline void Vector::asNode(const ArrD &elemvec, MatD &nodevec) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );
  assert( static_cast<size_t>(nodevec.rows())  == m_nnode );
  assert( static_cast<size_t>(nodevec.cols())  == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        nodevec(m_conn(e,m),i) = elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline MatD Vector::asNode(const ArrD &elemvec) const
{
  MatD nodevec(m_nnode, m_ndim);

  this->asNode(elemvec, nodevec);

  return nodevec;
}

// --------------------------------------- dofval -> elemvec ---------------------------------------

inline void Vector::asElement(const ColD &dofval, ArrD &elemvec) const
{
  assert( static_cast<size_t>(dofval.size()) == m_ndof );
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        elemvec(e,m,i) = dofval(m_dofs(m_conn(e,m),i));
}

// -------------------------------------------------------------------------------------------------

inline ArrD Vector::asElement(const ColD &dofval) const
{
  ArrD elemvec({m_nelem, m_nne, m_ndim});

  this->asElement(dofval, elemvec);

  return elemvec;
}

// --------------------------------------- dofval -> elemvec ---------------------------------------

inline void Vector::asElement(const ColD &dofval_u, const ColD &dofval_p, ArrD &elemvec) const
{
  assert( static_cast<size_t>(dofval_u.size()) == m_nnu );
  assert( static_cast<size_t>(dofval_p.size()) == m_nnp );
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

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

inline ArrD Vector::asElement(const ColD &dofval_u, const ColD &dofval_p) const
{
  ArrD elemvec({m_nelem, m_nne, m_ndim});

  this->asElement(dofval_u, dofval_p, elemvec);

  return elemvec;
}

// -------------------------------------- nodevec -> elemvec ---------------------------------------

inline void Vector::asElement(const MatD &nodevec, ArrD &elemvec) const
{
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        elemvec(e,m,i) = nodevec(m_conn(e,m),i);
}

// -------------------------------------------------------------------------------------------------

inline ArrD Vector::asElement(const MatD &nodevec) const
{
  ArrD elemvec({m_nelem, m_nne, m_ndim});

  this->asElement(nodevec, elemvec);

  return elemvec;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline void Vector::assembleDofs(const MatD &nodevec, ColD &dofval) const
{
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );
  assert( static_cast<size_t>(dofval.size())  == m_ndof  );

  dofval.setZero();

  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      dofval(m_dofs(n,i)) += nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::assembleDofs(const MatD &nodevec) const
{
  ColD dofval(m_ndof);

  this->assembleDofs(nodevec, dofval);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline void Vector::assembleDofs_u(const MatD &nodevec, ColD &dofval) const
{
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );
  assert( static_cast<size_t>(dofval.size())  == m_nnu   );

  dofval.setZero();

  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) < m_nnu )
        dofval(m_part(n,i)) += nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::assembleDofs_u(const MatD &nodevec) const
{
  ColD dofval(m_nnu);

  this->assembleDofs_u(nodevec, dofval);

  return dofval;
}

// --------------------------------------- nodevec -> dofval ---------------------------------------

inline void Vector::assembleDofs_p(const MatD &nodevec, ColD &dofval) const
{
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );
  assert( static_cast<size_t>(dofval.size())  == m_nnp   );

  dofval.setZero();

  for ( size_t n = 0 ; n < m_nnode ; ++n )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      if ( m_part(n,i) >= m_nnu )
        dofval(m_part(n,i)-m_nnu) += nodevec(n,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::assembleDofs_p(const MatD &nodevec) const
{
  ColD dofval(m_nnp);

  this->assembleDofs_p(nodevec, dofval);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline void Vector::assembleDofs(const ArrD &elemvec, ColD &dofval) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );
  assert( static_cast<size_t>(dofval.size()) == m_ndof );

  dofval.setZero();

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        dofval(m_dofs(m_conn(e,m),i)) += elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::assembleDofs(const ArrD &elemvec) const
{
  ColD dofval(m_ndof);

  this->assembleDofs(elemvec, dofval);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline void Vector::assembleDofs_u(const ArrD &elemvec, ColD &dofval) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );
  assert( static_cast<size_t>(dofval.size()) == m_nnu );

  dofval.setZero();

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) < m_nnu )
          dofval(m_dofs(m_conn(e,m),i)) += elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::assembleDofs_u(const ArrD &elemvec) const
{
  ColD dofval(m_nnu);

  this->assembleDofs_u(elemvec, dofval);

  return dofval;
}

// --------------------------------------- elemvec -> dofval ---------------------------------------

inline void Vector::assembleDofs_p(const ArrD &elemvec, ColD &dofval) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );
  assert( static_cast<size_t>(dofval.size()) == m_nnp );

  dofval.setZero();

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        if ( m_part(m_conn(e,m),i) >= m_nnu )
          dofval(m_dofs(m_conn(e,m),i)-m_nnu) += elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline ColD Vector::assembleDofs_p(const ArrD &elemvec) const
{
  ColD dofval(m_nnp);

  this->assembleDofs_p(elemvec, dofval);

  return dofval;
}

// --------------------------------------- elemvec -> nodevec ---------------------------------------

inline void Vector::assembleNode(const ArrD &elemvec, MatD &nodevec) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );
  assert( static_cast<size_t>(nodevec.rows()) == m_nnode );
  assert( static_cast<size_t>(nodevec.cols()) == m_ndim  );

  nodevec.setZero();

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t m = 0 ; m < m_nne ; ++m )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        nodevec(m_conn(e,m),i) += elemvec(e,m,i);
}

// -------------------------------------------------------------------------------------------------

inline MatD Vector::assembleNode(const ArrD &elemvec) const
{
  MatD nodevec(m_nnode, m_ndim);

  this->assembleNode(elemvec, nodevec);

  return nodevec;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

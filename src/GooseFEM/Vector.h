/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_VECTOR_H
#define GOOSEFEM_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// =========================================== GooseFEM ============================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

class Vector
{
private:
  MatS   m_conn;
  MatS   m_dofs;
  MatS   m_part;
  ColS   m_iip;
  ColS   m_iiu;

  size_t m_nelem;
  size_t m_nne;
  size_t m_nnode;
  size_t m_ndim;
  size_t m_ndof;
  size_t m_nnu;
  size_t m_nnp;

public:
  // constructor
  Vector(const MatS &conn, const MatS &dofs, const ColS &iip=ColS());

  // return dimensions
  size_t nelem();
  size_t nne();
  size_t nnode();
  size_t ndim();
  size_t ndof();
  size_t nnu();
  size_t nnp();

  // return DOF lists
  ColS iiu(); // unknown    DOFs
  ColS iip(); // prescribed DOFs

  // convert vectors; the following notation is used (unique overload from matrix-type):
  //
  // "nodevec"  -  nodal vectors                     -  MatD  -  [nnode, ndim]
  // "elemvec"  -  nodal vectors stored per element  -  ArrD  -  [nelem, nne, ndim]
  // "dofval"   -  DOF values                        -  ColD  -  [ndof]
  // "dofvalU"  -  DOF values (Unknown)              -  ColD  -  [nnu]
  // "dofvalP"  -  DOF values (Prescribed)           -  ColD  -  [nnp]
  //
  ColD asDofs   (const ColD &dofvalU, const ColD &dofvalP); // "dofval"   ->  "dofval"
  ColD asDofs   (const MatD &nodevec                     ); // "nodevec"  ->  "dofval"
  ColD asDofs   (const ArrD &elemvec                     ); // "elemvec"  ->  "dofval"
  ColD asDofs_u (const MatD &nodevec                     ); // "nodevec"  ->  "dofvalU"
  ColD asDofs_u (const ArrD &elemvec                     ); // "elemvec"  ->  "dofvalU"
  ColD asDofs_p (const MatD &nodevec                     ); // "nodevec"  ->  "dofvalP"
  ColD asDofs_p (const ArrD &elemvec                     ); // "elemvec"  ->  "dofvalP"
  MatD asNode   (const ColD &dofval                      ); // "dofval"   ->  "nodevec"
  MatD asNode   (const ColD &dofvalU, const ColD &dofvalP); // "dofval"   ->  "nodevec"
  MatD asNode   (const ArrD &elemvec                     ); // "elemvec"  ->  "nodevec"
  ArrD asElement(const ColD &dofval                      ); // "dofval"   ->  "elemvec"
  ArrD asElement(const ColD &dofvalU, const ColD &dofvalP); // "dofval"   ->  "elemvec"
  ArrD asElement(const MatD &nodevec                     ); // "nodevec"  ->  "elemvec"

  // assemble vectors (see notation and overload above)
  ColD assembleDofs  (const MatD &nodevec); // "nodevec"  ->  "dofval"
  ColD assembleDofs  (const ArrD &elemvec); // "elemvec"  ->  "dofval"
  ColD assembleDofs_u(const MatD &nodevec); // "nodevec"  ->  "dofval"
  ColD assembleDofs_u(const ArrD &elemvec); // "elemvec"  ->  "dofval"
  ColD assembleDofs_p(const MatD &nodevec); // "nodevec"  ->  "dofval"
  ColD assembleDofs_p(const ArrD &elemvec); // "elemvec"  ->  "dofval"
  MatD assembleNode  (const ArrD &elemvec); // "elemvec"  ->  "nodevec"

};

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

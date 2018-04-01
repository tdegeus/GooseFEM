/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIXDIAGONAL_H
#define GOOSEFEM_MATRIXDIAGONAL_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// =========================================== GooseFEM ============================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

class DiagonalMatrix
{
private:

public:
  // constructor
  DiagonalMatrix(const MatS &conn, const MatS &dofs, const ColS &iip=ColS());

  ColD solve(const ColD &rhs);
  ColD solve(const ColD &rhs, const ColD &up);

  ColD asDiagonal   ();
  ColD asDiagonal_uu();
  ColD asDiagonal_pp();

  SpMatD asSparse   ();
  SpMatD asSparse_uu();
  SpMatD asSparse_up();
  SpMatD asSparse_pu();
  SpMatD asSparse_pp();

  MatD asDense   ();
  MatD asDense_uu();
  MatD asDense_up();
  MatD asDense_pu();
  MatD asDense_pp();

};

// -------------------------------------------------------------------------------------------------

inline ColD operator* (const DiagonalMatrix &A, const ColD &b);

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

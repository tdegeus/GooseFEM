/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_VECTOR_H
#define XGOOSEFEM_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "GooseFEM.h"

// =========================================== GooseFEM ============================================

namespace xGooseFEM {

// -------------------------------------------------------------------------------------------------

/*
  "nodevec"   -  nodal vectors                            -  [nnode, ndim]
  "elemvec"   -  nodal vectors stored per element         -  [nelem, nne, ndim]
  "dofval"    -  DOF values                               -  [ndof]
  "dofval_u"  -  DOF values (Unknown)    "== dofval[iiu]" -  [nnu]
  "dofval_p"  -  DOF values (Prescribed) "== dofval[iiu]" -  [nnp]
*/

class Vector
{
public:

  // constructor

  Vector() = default;

  Vector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs);

  Vector(const xt::xtensor<size_t,2> &conn, const xt::xtensor<size_t,2> &dofs,
    const xt::xtensor<size_t,1> &iip);

  // dimensions

  size_t nelem() const; // number of elements
  size_t nne()   const; // number of nodes per element
  size_t nnode() const; // number of nodes
  size_t ndim()  const; // number of dimensions
  size_t ndof()  const; // number of DOFs
  size_t nnu()   const; // number of unknown DOFs
  size_t nnp()   const; // number of prescribed DOFs

  // DOF lists

  xt::xtensor<size_t,2> dofs() const; // DOFs
  xt::xtensor<size_t,1> iiu()  const; // unknown    DOFs
  xt::xtensor<size_t,1> iip()  const; // prescribed DOFs

  // convert to "dofval" (overwrite entries that occur more than once) -- (auto allocation below)

  void asDofs  (const xt::xtensor<double,1> &dofval_u, const xt::xtensor<double,1> &dofval_p,
    xt::xtensor<double,1> &dofval) const;

  void asDofs  (const xt::xtensor<double,2> &nodevec,
    xt::xtensor<double,1> &dofval) const;

  void asDofs  (const xt::xtensor<double,3> &elemvec,
    xt::xtensor<double,1> &dofval) const;

  void asDofs_u(const xt::xtensor<double,2> &nodevec,
    xt::xtensor<double,1> &dofval_u) const;

  void asDofs_u(const xt::xtensor<double,3> &elemvec,
    xt::xtensor<double,1> &dofval_u) const;

  void asDofs_p(const xt::xtensor<double,2> &nodevec,
    xt::xtensor<double,1> &dofval_p) const;

  void asDofs_p(const xt::xtensor<double,3> &elemvec,
    xt::xtensor<double,1> &dofval_p) const;

  // convert to "nodevec" (overwrite entries that occur more than once) -- (auto allocation below)

  void asNode(const xt::xtensor<double,1> &dofval,
    xt::xtensor<double,2> &nodevec) const;

  void asNode(const xt::xtensor<double,1> &dofval_u, const xt::xtensor<double,1> &dofval_p,
    xt::xtensor<double,2> &nodevec) const;

  void asNode(const xt::xtensor<double,3> &elemvec,
    xt::xtensor<double,2> &nodevec) const;

  // convert to "elemvec" (overwrite entries that occur more than once) -- (auto allocation below)

  void asElement(const xt::xtensor<double,1> &dofval,
    xt::xtensor<double,3> &elemvec) const;

  void asElement(const xt::xtensor<double,1> &dofval_u, const xt::xtensor<double,1> &dofval_p,
    xt::xtensor<double,3> &elemvec) const;

  void asElement(const xt::xtensor<double,2> &nodevec,
    xt::xtensor<double,3> &elemvec) const;

  // convert to "dofval" (overwrite entries that occur more than once)

  xt::xtensor<double,1> asDofs  (const xt::xtensor<double,1> &dofval_u,
    const xt::xtensor<double,1> &dofval_p) const;

  xt::xtensor<double,1> asDofs  (const xt::xtensor<double,2> &nodevec) const;

  xt::xtensor<double,1> asDofs  (const xt::xtensor<double,3> &elemvec) const;

  xt::xtensor<double,1> asDofs_u(const xt::xtensor<double,2> &nodevec) const;

  xt::xtensor<double,1> asDofs_u(const xt::xtensor<double,3> &elemvec) const;

  xt::xtensor<double,1> asDofs_p(const xt::xtensor<double,2> &nodevec) const;

  xt::xtensor<double,1> asDofs_p(const xt::xtensor<double,3> &elemvec) const;

  // convert to "nodevec" (overwrite entries that occur more than once)

  xt::xtensor<double,2> asNode(const xt::xtensor<double,1> &dofval_u,
    const xt::xtensor<double,1> &dofval_p) const;

  xt::xtensor<double,2> asNode(const xt::xtensor<double,1> &dofval) const;

  xt::xtensor<double,2> asNode(const xt::xtensor<double,3> &elemvec) const;

  // convert to "elemvec" (overwrite entries that occur more than once)

  xt::xtensor<double,3> asElement(const xt::xtensor<double,1> &dofval_u,
    const xt::xtensor<double,1> &dofval_p) const;

  xt::xtensor<double,3> asElement(const xt::xtensor<double,1> &dofval) const;

  xt::xtensor<double,3> asElement(const xt::xtensor<double,2> &nodevec) const;

  // assemble "dofval" (adds entries that occur more that once) -- (auto allocation below)

  void assembleDofs  (const xt::xtensor<double,2> &nodevec, xt::xtensor<double,1> &dofval  ) const;
  void assembleDofs  (const xt::xtensor<double,3> &elemvec, xt::xtensor<double,1> &dofval  ) const;
  void assembleDofs_u(const xt::xtensor<double,2> &nodevec, xt::xtensor<double,1> &dofval_u) const;
  void assembleDofs_u(const xt::xtensor<double,3> &elemvec, xt::xtensor<double,1> &dofval_u) const;
  void assembleDofs_p(const xt::xtensor<double,2> &nodevec, xt::xtensor<double,1> &dofval_p) const;
  void assembleDofs_p(const xt::xtensor<double,3> &elemvec, xt::xtensor<double,1> &dofval_p) const;

  // assemble "nodevec" (adds entries that occur more that once) -- (auto allocation below)

  void assembleNode(const xt::xtensor<double,3> &elemvec, xt::xtensor<double,2> &nodevec) const;

  // assemble "dofval" (adds entries that occur more that once)

  xt::xtensor<double,1> assembleDofs  (const xt::xtensor<double,2> &nodevec) const;
  xt::xtensor<double,1> assembleDofs  (const xt::xtensor<double,3> &elemvec) const;
  xt::xtensor<double,1> assembleDofs_u(const xt::xtensor<double,2> &nodevec) const;
  xt::xtensor<double,1> assembleDofs_u(const xt::xtensor<double,3> &elemvec) const;
  xt::xtensor<double,1> assembleDofs_p(const xt::xtensor<double,2> &nodevec) const;
  xt::xtensor<double,1> assembleDofs_p(const xt::xtensor<double,3> &elemvec) const;

  // assemble "nodevec" (adds entries that occur more that once)

  xt::xtensor<double,2> assembleNode(const xt::xtensor<double,3> &elemvec) const;

private:

  // bookkeeping
  xt::xtensor<size_t,2> m_conn; // connectivity                    [nelem, nne ]
  xt::xtensor<size_t,2> m_dofs; // DOF-numbers per node            [nnode, ndim]
  xt::xtensor<size_t,1> m_iiu;  // DOF-numbers that are unknown    [nnu]
  xt::xtensor<size_t,1> m_iip;  // DOF-numbers that are prescribed [nnp]

  // DOFs per node, such that iiu = arange(nnu), iip = nnu + arange(nnp)
  xt::xtensor<size_t,2> m_part;

  // dimensions
  size_t m_nelem; // number of elements
  size_t m_nne;   // number of nodes per element
  size_t m_nnode; // number of nodes
  size_t m_ndim;  // number of dimensions
  size_t m_ndof;  // number of DOFs
  size_t m_nnu;   // number of unknown DOFs
  size_t m_nnp;   // number of prescribed DOFs

};

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_VECTORPARTITIONED_H
#define GOOSEFEM_VECTORPARTITIONED_H

#include "config.h"
#include "Vector.h"

namespace GooseFEM {

/*
  "nodevec"   -  nodal vectors                            -  [nnode, ndim]
  "elemvec"   -  nodal vectors stored per element         -  [nelem, nne, ndim]
  "dofval"    -  DOF values                               -  [ndof]
  "dofval_u"  -  DOF values (Unknown)    "== dofval[iiu]" -  [nnu]
  "dofval_p"  -  DOF values (Prescribed) "== dofval[iiu]" -  [nnp]
*/

class VectorPartitioned : public Vector {
public:

    // making overloads from parent visible
    using Vector::asDofs;
    using Vector::asNode;
    using Vector::asElement;
    using Vector::AsDofs;
    using Vector::AsNode;
    using Vector::AsElement;

    // Constructor
    VectorPartitioned() = default;

    VectorPartitioned(
        const xt::xtensor<size_t, 2>& conn,
        const xt::xtensor<size_t, 2>& dofs,
        const xt::xtensor<size_t, 1>& iip);

    // Dimensions
    size_t nnu() const;   // number of unknown DOFs
    size_t nnp() const;   // number of prescribed DOFs

    // DOF lists
    xt::xtensor<size_t, 1> iiu() const;  // unknown    DOFs
    xt::xtensor<size_t, 1> iip() const;  // prescribed DOFs

    // Copy (part of) nodevec/dofval to another nodevec/dofval
    void copy_u(
        const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const; // "iiu" updated

    void copy_p(
        const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const; // "iip"  updated

    // Convert to "dofval" (overwrite entries that occur more than once)

    void asDofs(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p,
        xt::xtensor<double, 1>& dofval) const;

    void asDofs_u(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 1>& dofval_u) const;
    void asDofs_u(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval_u) const;
    void asDofs_u(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval_u) const;
    void asDofs_p(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 1>& dofval_p) const;
    void asDofs_p(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval_p) const;
    void asDofs_p(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval_p) const;

    // Convert to "nodevec" (overwrite entries that occur more than once) -- (auto allocation below)
    void asNode(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p,
        xt::xtensor<double, 2>& nodevec) const;

    // Convert to "elemvec" (overwrite entries that occur more than once) -- (auto allocation below)
    void asElement(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p,
        xt::xtensor<double, 3>& elemvec) const;

    // Assemble "dofval" (adds entries that occur more that once) -- (auto allocation below)
    void assembleDofs_u(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval_u) const;
    void assembleDofs_u(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval_u) const;
    void assembleDofs_p(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval_p) const;
    void assembleDofs_p(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval_p) const;

    // Auto-allocation of the functions above
    xt::xtensor<double, 1> AsDofs(
        const xt::xtensor<double, 1>& dofval_u, const xt::xtensor<double, 1>& dofval_p) const;

    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 1>& dofval) const;
    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 1>& dofval) const;
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 3>& elemvec) const;

    xt::xtensor<double, 2> AsNode(
        const xt::xtensor<double, 1>& dofval_u, const xt::xtensor<double, 1>& dofval_p) const;

    xt::xtensor<double, 3> AsElement(
        const xt::xtensor<double, 1>& dofval_u, const xt::xtensor<double, 1>& dofval_p) const;

    xt::xtensor<double, 1> AssembleDofs_u(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs_u(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 1> AssembleDofs_p(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs_p(const xt::xtensor<double, 3>& elemvec) const;

    xt::xtensor<double, 2> Copy_u(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

    /**
    Copy prescribed DOFs from an "src" to "dest".

    \param nodevec_src The input, from which to copy the prescribed DOFs.
    \param nodevec_dest The destination, to which to copy the prescribed DOFs.
    \returns The result after copying.
    */
    xt::xtensor<double, 2> Copy_p(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

protected:
    // Bookkeeping
    xt::xtensor<size_t, 1> m_iiu;  // DOF-numbers that are unknown    [nnu]
    xt::xtensor<size_t, 1> m_iip;  // DOF-numbers that are prescribed [nnp]

    // DOFs per node, such that iiu = arange(nnu), iip = nnu + arange(nnp)
    xt::xtensor<size_t, 2> m_part;

    // Dimensions
    size_t m_nnu;   // number of unknown DOFs
    size_t m_nnp;   // number of prescribed DOFs
};

} // namespace GooseFEM

#include "VectorPartitioned.hpp"

#endif

/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_VECTORPARTITIONED_H
#define GOOSEFEM_VECTORPARTITIONED_H

#include "config.h"

namespace GooseFEM {

/*
  "nodevec"   -  nodal vectors                            -  [nnode, ndim]
  "elemvec"   -  nodal vectors stored per element         -  [nelem, nne, ndim]
  "dofval"    -  DOF values                               -  [ndof]
  "dofval_u"  -  DOF values (Unknown)    "== dofval[iiu]" -  [nnu]
  "dofval_p"  -  DOF values (Prescribed) "== dofval[iiu]" -  [nnp]
*/

class VectorPartitioned {
public:
    // Constructor
    VectorPartitioned() = default;

    VectorPartitioned(
        const xt::xtensor<size_t, 2>& conn,
        const xt::xtensor<size_t, 2>& dofs,
        const xt::xtensor<size_t, 1>& iip);

    // Dimensions
    size_t nelem() const; // number of elements
    size_t nne() const;   // number of nodes per element
    size_t nnode() const; // number of nodes
    size_t ndim() const;  // number of dimensions
    size_t ndof() const;  // number of DOFs
    size_t nnu() const;   // number of unknown DOFs
    size_t nnp() const;   // number of prescribed DOFs

    // DOF lists
    xt::xtensor<size_t, 2> dofs() const; // DOFs
    xt::xtensor<size_t, 1> iiu() const;  // unknown    DOFs
    xt::xtensor<size_t, 1> iip() const;  // prescribed DOFs

    // Copy (part of) nodevec/dofval to another nodevec/dofval
    void copy(
        const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const;

    void copy_u(
        const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const; // "iiu" updated

    void copy_p(
        const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const; // "iip"  updated

    // Convert to "dofval" (overwrite entries that occur more than once)
    void asDofs(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p,
        xt::xtensor<double, 1>& dofval) const;

    void asDofs(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval) const;
    void asDofs(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval) const;
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

    void asNode(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 2>& nodevec) const;

    void asNode(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 2>& nodevec) const;

    // Convert to "elemvec" (overwrite entries that occur more than once) -- (auto allocation below)
    void asElement(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p,
        xt::xtensor<double, 3>& elemvec) const;

    void asElement(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 3>& elemvec) const;
    void asElement(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 3>& elemvec) const;

    // Assemble "dofval" (adds entries that occur more that once) -- (auto allocation below)
    void assembleDofs(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval) const;
    void assembleDofs(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval) const;
    void assembleDofs_u(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval_u) const;
    void assembleDofs_u(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval_u) const;
    void assembleDofs_p(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval_p) const;
    void assembleDofs_p(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval_p) const;

    // Assemble "nodevec" (adds entries that occur more that once) -- (auto allocation below)
    void assembleNode(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 2>& nodevec) const;

    // Auto-allocation of the functions above
    xt::xtensor<double, 1> AsDofs(
        const xt::xtensor<double, 1>& dofval_u, const xt::xtensor<double, 1>& dofval_p) const;

    xt::xtensor<double, 1> AsDofs(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AsDofs(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 1>& dofval) const;
    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 1>& dofval) const;
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 3>& elemvec) const;

    xt::xtensor<double, 2> AsNode(
        const xt::xtensor<double, 1>& dofval_u, const xt::xtensor<double, 1>& dofval_p) const;

    xt::xtensor<double, 2> AsNode(const xt::xtensor<double, 1>& dofval) const;
    xt::xtensor<double, 2> AsNode(const xt::xtensor<double, 3>& elemvec) const;

    xt::xtensor<double, 3> AsElement(
        const xt::xtensor<double, 1>& dofval_u, const xt::xtensor<double, 1>& dofval_p) const;

    xt::xtensor<double, 3> AsElement(const xt::xtensor<double, 1>& dofval) const;
    xt::xtensor<double, 3> AsElement(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 1> AssembleDofs_u(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs_u(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 1> AssembleDofs_p(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs_p(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 2> AssembleNode(const xt::xtensor<double, 3>& elemvec) const;

    xt::xtensor<double, 2> Copy(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

    xt::xtensor<double, 2> Copy_u(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

    xt::xtensor<double, 2> Copy_p(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

    // Get zero-allocated dofval, nodevec, elemvec
    xt::xtensor<double, 1> AllocateDofval() const;
    xt::xtensor<double, 2> AllocateNodevec() const;
    xt::xtensor<double, 3> AllocateElemvec() const;
    xt::xtensor<double, 1> AllocateDofval(double val) const;
    xt::xtensor<double, 2> AllocateNodevec(double val) const;
    xt::xtensor<double, 3> AllocateElemvec(double val) const;

private:
    // Bookkeeping
    xt::xtensor<size_t, 2> m_conn; // connectivity                    [nelem, nne ]
    xt::xtensor<size_t, 2> m_dofs; // DOF-numbers per node            [nnode, ndim]
    xt::xtensor<size_t, 1> m_iiu;  // DOF-numbers that are unknown    [nnu]
    xt::xtensor<size_t, 1> m_iip;  // DOF-numbers that are prescribed [nnp]

    // DOFs per node, such that iiu = arange(nnu), iip = nnu + arange(nnp)
    xt::xtensor<size_t, 2> m_part;

    // Dimensions
    size_t m_nelem; // number of elements
    size_t m_nne;   // number of nodes per element
    size_t m_nnode; // number of nodes
    size_t m_ndim;  // number of dimensions
    size_t m_ndof;  // number of DOFs
    size_t m_nnu;   // number of unknown DOFs
    size_t m_nnp;   // number of prescribed DOFs
};

} // namespace GooseFEM

#include "VectorPartitioned.hpp"

#endif

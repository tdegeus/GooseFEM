/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_VECTOR_H
#define GOOSEFEM_VECTOR_H

#include "config.h"

namespace GooseFEM {

/*
  "nodevec"  -  nodal vectors                     -  [nnode, ndim]
  "elemvec"  -  nodal vectors stored per element  -  [nelem, nne, ndim]
  "dofval"   -  DOF values                        -  [ndof]
*/

class Vector {
public:
    // Constructor
    Vector() = default;
    Vector(const xt::xtensor<size_t, 2>& conn, const xt::xtensor<size_t, 2>& dofs);

    // Dimensions
    size_t nelem() const; // number of elements
    size_t nne() const;   // number of nodes per element
    size_t nnode() const; // number of nodes
    size_t ndim() const;  // number of dimensions
    size_t ndof() const;  // number of DOFs

    // DOF lists
    xt::xtensor<size_t, 2> dofs() const; // DOFs

    // Copy nodevec to another nodevec
    void copy(const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const;

    // Convert to "dofval" (overwrite entries that occur more than once) -- (auto allocation below)
    void asDofs(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval) const;
    void asDofs(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval) const;

    // Convert to "nodevec" (overwrite entries that occur more than once) -- (auto allocation below)
    void asNode(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 2>& nodevec) const;
    void asNode(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 2>& nodevec) const;

    // Convert to "elemvec" (overwrite entries that occur more than once) -- (auto allocation below)
    void asElement(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 3>& elemvec) const;
    void asElement(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 3>& elemvec) const;

    // Assemble "dofval" (adds entries that occur more that once) -- (auto allocation below)
    void assembleDofs(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval) const;
    void assembleDofs(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval) const;

    // Assemble "nodevec" (adds entries that occur more that once) -- (auto allocation below)
    void assembleNode(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 2>& nodevec) const;

    // Auto-allocation of the functions above
    xt::xtensor<double, 1> AsDofs(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AsDofs(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 2> AsNode(const xt::xtensor<double, 1>& dofval) const;
    xt::xtensor<double, 2> AsNode(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 3> AsElement(const xt::xtensor<double, 1>& dofval) const;
    xt::xtensor<double, 3> AsElement(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 2> AssembleNode(const xt::xtensor<double, 3>& elemvec) const;

private:
    // Bookkeeping
    xt::xtensor<size_t, 2> m_conn; // connectivity         [nelem, nne ]
    xt::xtensor<size_t, 2> m_dofs; // DOF-numbers per node [nnode, ndim]

    // Dimensions
    size_t m_nelem; // number of elements
    size_t m_nne;   // number of nodes per element
    size_t m_nnode; // number of nodes
    size_t m_ndim;  // number of dimensions
    size_t m_ndof;  // number of DOFs
};

} // namespace GooseFEM

#include "Vector.hpp"

#endif

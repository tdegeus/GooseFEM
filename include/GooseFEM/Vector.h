/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_VECTOR_H
#define GOOSEFEM_VECTOR_H

#include "config.h"

namespace GooseFEM {

/**
Class to switch between:

-   "dofval": DOF values, shape: ``[ndof]``.
-   "nodevec": nodal vectors, shape ``[nnode, ndim]``.
-   "elemvec": nodal vectors stored per element, shape: ``[nelem, nne, ndim]``.
*/
class Vector {
public:

    Vector() = default;

    /**
    Constructor.

    \param conn Connectivity, shape ``[nelem, nne]``.
    \param dofs DOFs per node, shape ``[nnode, ndim]``.
    */
    Vector(const xt::xtensor<size_t, 2>& conn, const xt::xtensor<size_t, 2>& dofs);

    /**
    Get number of elements.

    \return unsigned int
    */
    size_t nelem() const;

    /**
    Get number of nodes per element.

    \return unsigned int
    */
    size_t nne() const;

    /**
    Get number of nodes.

    \return unsigned int
    */
    size_t nnode() const;

    /**
    Get number of dimensions.

    \return unsigned int
    */
    size_t ndim() const;

    /**
    Get number of DOFs.

    \return unsigned int
    */
    size_t ndof() const;

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

    /**
    Convert to ``elemvec`` (overwrite entries that occur more than once).
    This function writes to the fully allocated last argument.
    To use auto-allocation use AsElement().

    \param dofval input [ndof()].
    \param elemvec output [nelem(), nne(), ndim()].
    */
    void asElement(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 3>& elemvec) const;

    /**
    Convert to ``elemvec`` (overwrite entries that occur more than once).
    This function writes to the fully allocated last argument.
    To use auto-allocation use AsElement().

    \param nodevec input [nnode(), ndim()].
    \param elemvec output [nelem(), nne(), ndim()].
    */
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

    /**
    Convert to ``elemvec`` (overwrite entries that occur more than once).
    See asElement() to avoid allocation of return data.

    \param dofval [ndof()].
    \return ``elemvec`` [nelem(), nne(), ndim()].
    */
    xt::xtensor<double, 3> AsElement(const xt::xtensor<double, 1>& dofval) const;

    /**
    Convert to ``elemvec`` (overwrite entries that occur more than once).
    See asElement() to avoid allocation of return data.

    \param nodevec [nnode(), ndim()].
    \return ``elemvec`` [nelem(), nne(), ndim()].
    */
    xt::xtensor<double, 3> AsElement(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs(const xt::xtensor<double, 2>& nodevec) const;
    xt::xtensor<double, 1> AssembleDofs(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 2> AssembleNode(const xt::xtensor<double, 3>& elemvec) const;

    xt::xtensor<double, 2> Copy(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

    // Get shape of dofval, nodevec, elemvec
    std::array<size_t, 1> ShapeDofval() const;
    std::array<size_t, 2> ShapeNodevec() const;
    std::array<size_t, 3> ShapeElemvec() const;
    std::array<size_t, 3> ShapeElemmat() const;

    // Get zero-allocated dofval, nodevec, elemvec
    xt::xtensor<double, 1> AllocateDofval() const;
    xt::xtensor<double, 2> AllocateNodevec() const;
    xt::xtensor<double, 3> AllocateElemvec() const;
    xt::xtensor<double, 3> AllocateElemmat() const;
    xt::xtensor<double, 1> AllocateDofval(double val) const;
    xt::xtensor<double, 2> AllocateNodevec(double val) const;
    xt::xtensor<double, 3> AllocateElemvec(double val) const;
    xt::xtensor<double, 3> AllocateElemmat(double val) const;

protected:
    xt::xtensor<size_t, 2> m_conn; ///< Connectivity ``[nelem, nne]``
    xt::xtensor<size_t, 2> m_dofs; ///< DOF-numbers per node ``[nnode, ndim]``
    size_t m_nelem; ///< Number of elements
    size_t m_nne;   ///< Number of nodes per element
    size_t m_nnode; ///< Number of nodes
    size_t m_ndim;  ///< Number of dimensions
    size_t m_ndof;  ///< Number of DOFs
};

} // namespace GooseFEM

#include "Vector.hpp"

#endif

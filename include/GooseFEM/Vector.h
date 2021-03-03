/**
Methods to switch between storage types based on a mesh.

\file Vector.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VECTOR_H
#define GOOSEFEM_VECTOR_H

#include "config.h"

namespace GooseFEM {

/**
Class to switch between:

-   "dofval": DOF values [#ndof].
-   "nodevec": nodal vectors [#nnode, #ndim].
-   "elemvec": nodal vectors stored per element [#nelem, #nne, #ndim].
*/
class Vector {
public:

    Vector() = default;

    /**
    Constructor.

    \param conn Connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    */
    Vector(const xt::xtensor<size_t, 2>& conn, const xt::xtensor<size_t, 2>& dofs);

    /**
    Number of elements.

    \return unsigned int
    */
    size_t nelem() const;

    /**
    Number of nodes per element.

    \return unsigned int
    */
    size_t nne() const;

    /**
    Number of nodes.

    \return unsigned int
    */
    size_t nnode() const;

    /**
    Number of dimensions.

    \return unsigned int
    */
    size_t ndim() const;

    /**
    Number of DOFs.

    \return unsigned int
    */
    size_t ndof() const;

    /**
    DOFs per node.

    \return [#nnode, #ndim]
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    Copy "nodevec" to another "nodevec".

    \param nodevec_src input [#nnode, #ndim]
    \param nodevec_dest input [#nnode, #ndim]
    \return nodevec output [#nnode, #ndim]
    */
    xt::xtensor<double, 2> Copy(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

    /**
    Copy "nodevec" to another "nodevec".

    \param nodevec_src input [#nnode, #ndim]
    \param nodevec_dest output [#nnode, #ndim]
    */
    void copy(const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const;

    /**
    Convert "nodevec" to "dofval" (overwrite entries that occur more than once).

    \param nodevec input [#nnode, #ndim]
    \return dofval output [#ndof]
    */
    xt::xtensor<double, 1> AsDofs(const xt::xtensor<double, 2>& nodevec) const;

    /**
    Convert "nodevec" to "dofval" (overwrite entries that occur more than once).

    \param nodevec input [#nnode, #ndim]
    \param dofval output [#ndof]
    */
    void asDofs(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval) const;

    /**
    Convert "elemvec" to "dofval" (overwrite entries that occur more than once).

    \param elemvec input [#nelem, #nne, #ndim]
    \return dofval output [#ndof]
    */
    xt::xtensor<double, 1> AsDofs(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Convert "elemvec" to "dofval" (overwrite entries that occur more than once).

    \param elemvec input [#nelem, #nne, #ndim]
    \param dofval output [#ndof]
    */
    void asDofs(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval) const;

    /**
    Convert "dofval" to "nodevec" (overwrite entries that occur more than once).

    \param dofval input [#ndof]
    \return nodevec output [#nnode, #ndim]
    */
    xt::xtensor<double, 2> AsNode(const xt::xtensor<double, 1>& dofval) const;

    /**
    Convert "dofval" to "nodevec" (overwrite entries that occur more than once).

    \param dofval input [#ndof]
    \param nodevec input [#nnode, #ndim]
    */
    void asNode(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 2>& nodevec) const;

    /**
    Convert "elemvec" to "nodevec" (overwrite entries that occur more than once).

    \param elemvec input [#nelem, #nne, #ndim]
    \return nodevec output [#nnode, #ndim]
    */
    xt::xtensor<double, 2> AsNode(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Convert "elemvec" to "nodevec" (overwrite entries that occur more than once).

    \param elemvec input [#nelem, #nne, #ndim]
    \param nodevec [#nnode, #ndim]
    */
    void asNode(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 2>& nodevec) const;

    /**
    Convert "dofval" to "elemvec" (overwrite entries that occur more than once).

    \param dofval input [#ndof].
    \return elemvec output [#nelem, #nne, #ndim].
    */
    xt::xtensor<double, 3> AsElement(const xt::xtensor<double, 1>& dofval) const;

    /**
    Convert "dofval" to "elemvec" (overwrite entries that occur more than once).

    \param dofval input [#ndof].
    \param elemvec output [#nelem, #nne, #ndim].
    */
    void asElement(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 3>& elemvec) const;

    /**
    Convert "nodevec" to "elemvec" (overwrite entries that occur more than once).

    \param nodevec input [#nnode, #ndim].
    \return elemvec output [#nelem, #nne, #ndim].
    */
    xt::xtensor<double, 3> AsElement(const xt::xtensor<double, 2>& nodevec) const;

    /**
    Convert "nodevec" to "elemvec" (overwrite entries that occur more than once).

    \param nodevec input [#nnode, #ndim].
    \param elemvec output [#nelem, #nne, #ndim].
    */
    void asElement(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 3>& elemvec) const;

    /**
    Assemble "nodevec" to "dofval" (adds entries that occur more that once).

    \param nodevec input [#nnode, #ndim]
    \return dofval output [#ndof]
    */
    xt::xtensor<double, 1> AssembleDofs(const xt::xtensor<double, 2>& nodevec) const;

    /**
    Assemble "nodevec" to "dofval" (adds entries that occur more that once).

    \param nodevec input [#nnode, #ndim]
    \param dofval output [#ndof]
    */
    void assembleDofs(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval) const;

    /**
    Assemble "elemvec" to "dofval" (adds entries that occur more that once).

    \param elemvec input [#nelem, #nne, #ndim]
    \return dofval output [#ndof]
    */
    xt::xtensor<double, 1> AssembleDofs(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Assemble "elemvec" to "dofval" (adds entries that occur more that once).

    \param elemvec input [#nelem, #nne, #ndim]
    \param dofval output [#ndof]
    */
    void assembleDofs(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval) const;

    /**
    Assemble "elemvec" to "nodevec" (adds entries that occur more that once.

    \param elemvec input [#nelem, #nne, #ndim]
    \return nodevec output [#nnode, #ndim]
    */
    xt::xtensor<double, 2> AssembleNode(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Assemble "elemvec" to "nodevec" (adds entries that occur more that once.

    \param elemvec input [#nelem, #nne, #ndim]
    \param nodevec output [#nnode, #ndim]
    */
    void assembleNode(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 2>& nodevec) const;

    /**
    Shape of "dofval".

    \return [#ndof]
    */
    std::array<size_t, 1> ShapeDofval() const;

    /**
    Shape of "nodevec".

    \return [#nnode, #ndim]
    */
    std::array<size_t, 2> ShapeNodevec() const;

    /**
    Shape of "elemvec".

    \return [#nelem, #nne, #ndim]
    */
    std::array<size_t, 3> ShapeElemvec() const;

    /**
    Shape of "elemmat".

    \return [#nelem, #nne * #ndim, #nne * #ndim]
    */
    std::array<size_t, 3> ShapeElemmat() const;

    /**
    Allocated "dofval".

    \return [#ndof]
    */
    xt::xtensor<double, 1> AllocateDofval() const;

    /**
    Allocated "dofval".

    \param val Value with which to allocate.
    \return [#ndof]
    */
    xt::xtensor<double, 1> AllocateDofval(double val) const;

    /**
    Allocated "nodevec".

    \return [#nnode, #ndim]
    */
    xt::xtensor<double, 2> AllocateNodevec() const;

    /**
    Allocated "nodevec".

    \param val Value with which to allocate.
    \return [#nnode, #ndim]
    */
    xt::xtensor<double, 2> AllocateNodevec(double val) const;

    /**
    Allocated "elemvec".

    \return [#nelem, #nne, #ndim]
    */
    xt::xtensor<double, 3> AllocateElemvec() const;

    /**
    Allocated "elemvec".

    \param val Value with which to allocate.
    \return [#nelem, #nne, #ndim]
    */
    xt::xtensor<double, 3> AllocateElemvec(double val) const;

    /**
    Allocated "elemmat".

    \return [#nelem, #nne * #ndim, #nne * #ndim]
    */
    xt::xtensor<double, 3> AllocateElemmat() const;

    /**
    Allocated "elemmat".

    \param val Value with which to allocate.
    \return [#nelem, #nne * #ndim, #nne * #ndim]
    */
    xt::xtensor<double, 3> AllocateElemmat(double val) const;

protected:
    xt::xtensor<size_t, 2> m_conn; ///< See conn()
    xt::xtensor<size_t, 2> m_dofs; ///< See dofs()
    size_t m_nelem; ///< See #nelem
    size_t m_nne;   ///< See #nne
    size_t m_nnode; ///< See #nnode
    size_t m_ndim;  ///< See #ndim
    size_t m_ndof;  ///< See #ndof
};

} // namespace GooseFEM

#include "Vector.hpp"

#endif

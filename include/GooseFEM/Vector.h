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
Class to switch between storage types. In particular:

-   "dofval": DOF values [#ndof].
-   "nodevec": nodal vectors [#nnode, #ndim].
-   "elemvec": nodal vectors stored per element [#nelem, #nne, #ndim].
*/
class Vector {
public:
    Vector() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    */
    template <class S, class T>
    Vector(const S& conn, const T& dofs);

    /**
    \return  Number of elements.
    */
    size_t nelem() const;

    /**
    \return  Number of nodes per element.
    */
    size_t nne() const;

    /**
    \return  Number of nodes.
    */
    size_t nnode() const;

    /**
    \return  Number of dimensions.
    */
    size_t ndim() const;

    /**
    \return  Number of DOFs.
    */
    size_t ndof() const;

    /**
    \return Connectivity (nodes per element) [#nelem, #nne].
    */
    xt::xtensor<size_t, 2> conn() const;

    /**
    \return DOFs per node [#nnode, #ndim]
    */
    xt::xtensor<size_t, 2> dofs() const;

    /**
    Copy "nodevec" to another "nodevec".

    \param nodevec_src input [#nnode, #ndim]
    \param nodevec_dest input [#nnode, #ndim]
    \return nodevec output [#nnode, #ndim]
    */
    template <class T>
    T Copy(const T& nodevec_src, const T& nodevec_dest) const;

    /**
    Copy "nodevec" to another "nodevec".

    \param nodevec_src input [#nnode, #ndim]
    \param nodevec_dest output [#nnode, #ndim]
    */
    template <class T>
    void copy(const T& nodevec_src, T& nodevec_dest) const;

    /**
    Convert "nodevec" or "elemvec" to "dofval" (overwrite entries that occur more than once).

    \param arg nodevec [#nnode, #ndim] or elemvec [#nelem, #nne, #ndim]
    \return dofval [#ndof]
    */
    template <class T>
    xt::xtensor<double, 1> AsDofs(const T& arg) const;

    /**
    Convert "nodevec" or "elemvec" to "dofval" (overwrite entries that occur more than once).

    \param arg nodevec [#nnode, #ndim] or elemvec [#nelem, #nne, #ndim]
    \param dofval (output) [#ndof]
    */
    template <class T, class R>
    void asDofs(const T& arg, R& dofval) const;

    /**
    Convert "dofval" or "elemvec" to "nodevec" (overwrite entries that occur more than once).

    \param arg dofval [#ndof] or elemvec [#nelem, #nne, #ndim]
    \return nodevec output [#nnode, #ndim]
    */
    template <class T>
    xt::xtensor<double, 2> AsNode(const T& arg) const;

    /**
    Convert "dofval" or "elemvec" to "nodevec" (overwrite entries that occur more than once).

    \param arg dofval [#ndof] or elemvec [#nelem, #nne, #ndim]
    \param nodevec output [#nnode, #ndim]
    */
    template <class T, class R>
    void asNode(const T& arg, R& nodevec) const;

    /**
    Convert "dofval" or "nodevec" to "elemvec" (overwrite entries that occur more than once).

    \param arg dofval [#ndof] or nodevec [#nnode, #ndim].
    \return elemvec output [#nelem, #nne, #ndim].
    */
    template <class T>
    xt::xtensor<double, 3> AsElement(const T& arg) const;

    /**
    Convert "dofval" or "nodevec" to "elemvec" (overwrite entries that occur more than once).

    \param arg dofval [#ndof] or nodevec [#nnode, #ndim].
    \param elemvec output [#nelem, #nne, #ndim].
    */
    template <class T, class R>
    void asElement(const T& arg, R& elemvec) const;

    /**
    Assemble "nodevec" or "elemvec" to "dofval" (adds entries that occur more that once).

    \param arg nodevec [#nnode, #ndim] or elemvec [#nelem, #nne, #ndim]
    \return dofval output [#ndof]
    */
    template <class T>
    xt::xtensor<double, 1> AssembleDofs(const T& arg) const;

    /**
    Assemble "nodevec" or "elemvec" to "dofval" (adds entries that occur more that once).

    \param arg nodevec [#nnode, #ndim] or elemvec [#nelem, #nne, #ndim]
    \param dofval output [#ndof]
    */
    template <class T, class R>
    void assembleDofs(const T& arg, R& dofval) const;

    /**
    Assemble "elemvec" to "nodevec" (adds entries that occur more that once.

    \param arg elemvec [#nelem, #nne, #ndim]
    \return nodevec output [#nnode, #ndim]
    */
    template <class T>
    xt::xtensor<double, 2> AssembleNode(const T& arg) const;

    /**
    Assemble "elemvec" to "nodevec" (adds entries that occur more that once.

    \param arg elemvec [#nelem, #nne, #ndim]
    \param nodevec output [#nnode, #ndim]
    */
    template <class T, class R>
    void assembleNode(const T& arg, R& nodevec) const;

    /**
    Shape of "dofval".

    \return [#ndof]
    */
    std::array<size_t, 1> shape_dofval() const;

    /**
    Shape of "nodevec".

    \return [#nnode, #ndim]
    */
    std::array<size_t, 2> shape_nodevec() const;

    /**
    Shape of "elemvec".

    \return [#nelem, #nne, #ndim]
    */
    std::array<size_t, 3> shape_elemvec() const;

    /**
    Shape of "elemmat".

    \return [#nelem, #nne * #ndim, #nne * #ndim]
    */
    std::array<size_t, 3> shape_elemmat() const;

    /**
    Allocated "dofval".

    \return [#ndof]
    */
    xt::xtensor<double, 1> allocate_dofval() const;

    /**
    Allocated and initialised "dofval".

    \param val value to which to initialise.
    \return [#ndof]
    */
    xt::xtensor<double, 1> allocate_dofval(double val) const;

    /**
    Allocated "nodevec".

    \return [#nnode, #ndim]
    */
    xt::xtensor<double, 2> allocate_nodevec() const;

    /**
    Allocated and initialised "nodevec".

    \param val value to which to initialise.
    \return [#nnode, #ndim]
    */
    xt::xtensor<double, 2> allocate_nodevec(double val) const;

    /**
    Allocated "elemvec".

    \return [#nelem, #nne, #ndim]
    */
    xt::xtensor<double, 3> allocate_elemvec() const;

    /**
    Allocated and initialised "elemvec".

    \param val value to which to initialise.
    \return [#nelem, #nne, #ndim]
    */
    xt::xtensor<double, 3> allocate_elemvec(double val) const;

    /**
    Allocated "elemmat".

    \return [#nelem, #nne * #ndim, #nne * #ndim]
    */
    xt::xtensor<double, 3> allocate_elemmat() const;

    /**
    Allocated and initialised "elemmat".

    \param val value to which to initialise.
    \return [#nelem, #nne * #ndim, #nne * #ndim]
    */
    xt::xtensor<double, 3> allocate_elemmat(double val) const;

private:
    /** Distribution to relevant implementation of \copydoc asDofs(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void asDofs_impl(const T& arg, R& dofval) const;

    /** Distribution to relevant implementation of \copydoc asDofs(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 2, int> = 0>
    void asDofs_impl(const T& arg, R& dofval) const;

    /** Distribution to relevant implementation of \copydoc asDofs(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 3, int> = 0>
    void asDofs_impl(const T& arg, R& dofval) const;

    /** Distribution to relevant implementation of \copydoc asNode(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void asNode_impl(const T& arg, R& nodevec) const;

    /** Distribution to relevant implementation of \copydoc asNode(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 1, int> = 0>
    void asNode_impl(const T& arg, R& nodevec) const;

    /** Distribution to relevant implementation of \copydoc asNode(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 3, int> = 0>
    void asNode_impl(const T& arg, R& nodevec) const;

    /** Distribution to relevant implementation of \copydoc asElement(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void asElement_impl(const T& arg, R& elemvec) const;

    /** Distribution to relevant implementation of \copydoc asElement(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 1, int> = 0>
    void asElement_impl(const T& arg, R& elemvec) const;

    /** Distribution to relevant implementation of \copydoc asElement(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 2, int> = 0>
    void asElement_impl(const T& arg, R& elemvec) const;

    /** Distribution to relevant implementation of \copydoc assembleDofs(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void assembleDofs_impl(const T& arg, R& dofval) const;

    /** Distribution to relevant implementation of \copydoc assembleDofs(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 2, int> = 0>
    void assembleDofs_impl(const T& arg, R& dofval) const;

    /** Distribution to relevant implementation of \copydoc assembleDofs(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 3, int> = 0>
    void assembleDofs_impl(const T& arg, R& dofval) const;

    /** Distribution to relevant implementation of \copydoc assembleNode(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void assembleNode_impl(const T& arg, R& nodevec) const;

    /** Distribution to relevant implementation of \copydoc assembleNode(const T&, R&) const */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 3, int> = 0>
    void assembleNode_impl(const T& arg, R& nodevec) const;

    /** Implementation for 'nodevec' input of \copydoc asDofs(const T&, R&) const */
    template <class T, class R>
    void asDofs_impl_nodevec(const T& arg, R& dofval) const;

    /** Implementation for 'elemvec' input of \copydoc asDofs(const T&, R&) const */
    template <class T, class R>
    void asDofs_impl_elemvec(const T& arg, R& dofval) const;

    /** Implementation for 'dofval' input of \copydoc asNode(const T&, R&) const */
    template <class T, class R>
    void asNode_impl_dofval(const T& arg, R& nodevec) const;

    /** Implementation for 'elemvec' input of \copydoc asNode(const T&, R&) const */
    template <class T, class R>
    void asNode_impl_elemvec(const T& arg, R& nodevec) const;

    /** Implementation for 'dofval' input of \copydoc asElement(const T&, R&) const */
    template <class T, class R>
    void asElement_impl_dofval(const T& arg, R& elemvec) const;

    /** Implementation for 'nodevec' input of \copydoc asElement(const T&, R&) const */
    template <class T, class R>
    void asElement_impl_nodevec(const T& arg, R& elemvec) const;

    /** Implementation for 'nodevec' input of \copydoc assembleDofs(const T&, R&) const */
    template <class T, class R>
    void assembleDofs_impl_nodevec(const T& arg, R& dofval) const;

    /** Implementation for 'elemvec' input of \copydoc assembleDofs(const T&, R&) const */
    template <class T, class R>
    void assembleDofs_impl_elemvec(const T& arg, R& dofval) const;

    /** Implementation for 'elemvec' input of \copydoc assembleNode(const T&, R&) const */
    template <class T, class R>
    void assembleNode_impl_elemvec(const T& arg, R& nodevec) const;

protected:
    xt::xtensor<size_t, 2> m_conn; ///< See conn()
    xt::xtensor<size_t, 2> m_dofs; ///< See dofs()
    size_t m_nelem; ///< See #nelem
    size_t m_nne; ///< See #nne
    size_t m_nnode; ///< See #nnode
    size_t m_ndim; ///< See #ndim
    size_t m_ndof; ///< See #ndof
};

} // namespace GooseFEM

#include "Vector.hpp"

#endif

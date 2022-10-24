/**
 * Methods to switch between storage types based on a mesh.
 *
 * @file Vector.h
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#ifndef GOOSEFEM_VECTOR_H
#define GOOSEFEM_VECTOR_H

#include "config.h"

namespace GooseFEM {

/**
 * Class to switch between storage types. In particular:
 *
 * -   "dofval": DOF values [#ndof].
 * -   "nodevec": nodal vectors [#nnode, #ndim].
 * -   "elemvec": nodal vectors stored per element [#nelem, #nne, #ndim].
 */
class Vector {
public:
    Vector() = default;

    /**
     * Constructor.
     *
     * @param conn connectivity [#nelem, #nne].
     * @param dofs DOFs per node [#nnode, #ndim].
     */
    template <class S, class T>
    Vector(const S& conn, const T& dofs) : m_conn(conn), m_dofs(dofs)
    {
        GOOSEFEM_ASSERT(conn.dimension() == 2);
        GOOSEFEM_ASSERT(dofs.dimension() == 2);

        m_nelem = m_conn.shape(0);
        m_nne = m_conn.shape(1);
        m_nnode = m_dofs.shape(0);
        m_ndim = m_dofs.shape(1);
        m_ndof = xt::amax(m_dofs)() + 1;

        GOOSEFEM_ASSERT(xt::amax(m_conn)() + 1 <= m_nnode);
        GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
    }

    /**
     * @return  Number of elements.
     */
    size_t nelem() const
    {
        return m_nelem;
    }

    /**
     * @return  Number of nodes per element.
     */
    size_t nne() const
    {
        return m_nne;
    }

    /**
     * @return  Number of nodes.
     */
    size_t nnode() const
    {
        return m_nnode;
    }

    /**
     * @return  Number of dimensions.
     */
    size_t ndim() const
    {
        return m_ndim;
    }

    /**
     * @return  Number of DOFs.
     */
    size_t ndof() const
    {
        return m_ndof;
    }

    /**
     * @return Connectivity (nodes per element) [#nelem, #nne].
     */
    const array_type::tensor<size_t, 2>& conn() const
    {
        return m_conn;
    }

    /**
     * @return DOFs per node [#nnode, #ndim]
     */
    const array_type::tensor<size_t, 2>& dofs() const
    {
        return m_dofs;
    }

    /**
     * Copy "nodevec" to another "nodevec".
     *
     * @param nodevec_src input [#nnode, #ndim]
     * @param nodevec_dest input [#nnode, #ndim]
     * @return nodevec output [#nnode, #ndim]
     */
    template <class T>
    T Copy(const T& nodevec_src, const T& nodevec_dest) const
    {
        T ret = T::from_shape(nodevec_dest.shape());
        this->copy(nodevec_src, ret);
        return ret;
    }

    /**
     * Copy "nodevec" to another "nodevec".
     *
     * @param nodevec_src input [#nnode, #ndim]
     * @param nodevec_dest output [#nnode, #ndim]
     */
    template <class T>
    void copy(const T& nodevec_src, T& nodevec_dest) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(nodevec_src, this->shape_nodevec()));
        GOOSEFEM_ASSERT(xt::has_shape(nodevec_dest, this->shape_nodevec()));

        xt::noalias(nodevec_dest) = nodevec_src;
    }

    /**
     * Convert "nodevec" or "elemvec" to "dofval" (overwrite entries that occur more than once).
     *
     * @param arg nodevec [#nnode, #ndim] or elemvec [#nelem, #nne, #ndim]
     * @return dofval [#ndof]
     */
    template <class T>
    array_type::tensor<double, 1> AsDofs(const T& arg) const
    {
        array_type::tensor<double, 1> ret = xt::empty<double>(this->shape_dofval());
        this->asDofs_impl(arg, ret);
        return ret;
    }

    /**
     * Convert "nodevec" or "elemvec" to "dofval" (overwrite entries that occur more than once).
     *
     * @param arg nodevec [#nnode, #ndim] or elemvec [#nelem, #nne, #ndim]
     * @param ret dofval (output) [#ndof]
     */
    template <class T, class R>
    void asDofs(const T& arg, R& ret) const
    {
        this->asDofs_impl(arg, ret);
    }

    /**
     * Convert "dofval" or "elemvec" to "nodevec" (overwrite entries that occur more than once).
     *
     * @param arg dofval [#ndof] or elemvec [#nelem, #nne, #ndim]
     * @return nodevec output [#nnode, #ndim]
     */
    template <class T>
    array_type::tensor<double, 2> AsNode(const T& arg) const
    {
        array_type::tensor<double, 2> ret = xt::empty<double>(this->shape_nodevec());
        this->asNode_impl(arg, ret);
        return ret;
    }

    /**
     * Convert "dofval" or "elemvec" to "nodevec" (overwrite entries that occur more than once).
     *
     * @param arg dofval [#ndof] or elemvec [#nelem, #nne, #ndim]
     * @param ret nodevec, output [#nnode, #ndim]
     */
    template <class T, class R>
    void asNode(const T& arg, R& ret) const
    {
        this->asNode_impl(arg, ret);
    }

    /**
     * Convert "dofval" or "nodevec" to "elemvec" (overwrite entries that occur more than once).
     *
     * @param arg dofval [#ndof] or nodevec [#nnode, #ndim].
     * @return elemvec output [#nelem, #nne, #ndim].
     */
    template <class T>
    array_type::tensor<double, 3> AsElement(const T& arg) const
    {
        array_type::tensor<double, 3> ret = xt::empty<double>(this->shape_elemvec());
        this->asElement_impl(arg, ret);
        return ret;
    }

    /**
     * Convert "dofval" or "nodevec" to "elemvec" (overwrite entries that occur more than once).
     *
     * @param arg dofval [#ndof] or nodevec [#nnode, #ndim].
     * @param ret elemvec, output [#nelem, #nne, #ndim].
     */
    template <class T, class R>
    void asElement(const T& arg, R& ret) const
    {
        this->asElement_impl(arg, ret);
    }

    /**
     * Assemble "nodevec" or "elemvec" to "dofval" (adds entries that occur more that once).
     *
     * @param arg nodevec [#nnode, #ndim] or elemvec [#nelem, #nne, #ndim]
     * @return dofval output [#ndof]
     */
    template <class T>
    array_type::tensor<double, 1> AssembleDofs(const T& arg) const
    {
        array_type::tensor<double, 1> ret = xt::empty<double>(this->shape_dofval());
        this->assembleDofs_impl(arg, ret);
        return ret;
    }

    /**
     * Assemble "nodevec" or "elemvec" to "dofval" (adds entries that occur more that once).
     *
     * @param arg nodevec [#nnode, #ndim] or elemvec [#nelem, #nne, #ndim]
     * @param ret dofval, output [#ndof]
     */
    template <class T, class R>
    void assembleDofs(const T& arg, R& ret) const
    {
        this->assembleDofs_impl(arg, ret);
    }

    /**
     * Assemble "elemvec" to "nodevec" (adds entries that occur more that once.
     *
     * @param arg elemvec [#nelem, #nne, #ndim]
     * @return nodevec output [#nnode, #ndim]
     */
    template <class T>
    array_type::tensor<double, 2> AssembleNode(const T& arg) const
    {
        array_type::tensor<double, 2> ret = xt::empty<double>(this->shape_nodevec());
        this->assembleNode_impl(arg, ret);
        return ret;
    }

    /**
     * Assemble "elemvec" to "nodevec" (adds entries that occur more that once.
     *
     * @param arg elemvec [#nelem, #nne, #ndim]
     * @param ret nodevec, output [#nnode, #ndim]
     */
    template <class T, class R>
    void assembleNode(const T& arg, R& ret) const
    {
        this->assembleNode_impl(arg, ret);
    }

    /**
     * Shape of "dofval".
     *
     * @return [#ndof]
     */
    std::array<size_t, 1> shape_dofval() const
    {
        return std::array<size_t, 1>{m_ndof};
    }

    /**
     * Shape of "nodevec".
     *
     * @return [#nnode, #ndim]
     */
    std::array<size_t, 2> shape_nodevec() const
    {
        return std::array<size_t, 2>{m_nnode, m_ndim};
    }

    /**
     * Shape of "elemvec".
     *
     * @return [#nelem, #nne, #ndim]
     */
    std::array<size_t, 3> shape_elemvec() const
    {
        return std::array<size_t, 3>{m_nelem, m_nne, m_ndim};
    }

    /**
     * Shape of "elemmat".
     *
     * @return [#nelem, #nne * #ndim, #nne * #ndim]
     */
    std::array<size_t, 3> shape_elemmat() const
    {
        return std::array<size_t, 3>{m_nelem, m_nne * m_ndim, m_nne * m_ndim};
    }

    /**
     * Allocated "dofval".
     *
     * @return [#ndof]
     */
    array_type::tensor<double, 1> allocate_dofval() const
    {
        array_type::tensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
        return dofval;
    }

    /**
     * Allocated and initialised "dofval".
     *
     * @param val value to which to initialise.
     * @return [#ndof]
     */
    array_type::tensor<double, 1> allocate_dofval(double val) const
    {
        array_type::tensor<double, 1> dofval = xt::empty<double>(this->shape_dofval());
        dofval.fill(val);
        return dofval;
    }

    /**
     * Allocated "nodevec".
     *
     * @return [#nnode, #ndim]
     */
    array_type::tensor<double, 2> allocate_nodevec() const
    {
        array_type::tensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
        return nodevec;
    }

    /**
     * Allocated and initialised "nodevec".
     *
     * @param val value to which to initialise.
     * @return [#nnode, #ndim]
     */
    array_type::tensor<double, 2> allocate_nodevec(double val) const
    {
        array_type::tensor<double, 2> nodevec = xt::empty<double>(this->shape_nodevec());
        nodevec.fill(val);
        return nodevec;
    }

    /**
     * Allocated "elemvec".
     *
     * @return [#nelem, #nne, #ndim]
     */
    array_type::tensor<double, 3> allocate_elemvec() const
    {
        array_type::tensor<double, 3> elemvec = xt::empty<double>(this->shape_elemvec());
        return elemvec;
    }

    /**
     * Allocated and initialised "elemvec".
     *
     * @param val value to which to initialise.
     * @return [#nelem, #nne, #ndim]
     */
    array_type::tensor<double, 3> allocate_elemvec(double val) const
    {
        array_type::tensor<double, 3> elemvec = xt::empty<double>(this->shape_elemvec());
        elemvec.fill(val);
        return elemvec;
    }

    /**
     * Allocated "elemmat".
     *
     * @return [#nelem, #nne * #ndim, #nne * #ndim]
     */
    array_type::tensor<double, 3> allocate_elemmat() const
    {
        array_type::tensor<double, 3> elemmat = xt::empty<double>(this->shape_elemmat());
        return elemmat;
    }

    /**
     * Allocated and initialised "elemmat".
     *
     * @param val value to which to initialise.
     * @return [#nelem, #nne * #ndim, #nne * #ndim]
     */
    array_type::tensor<double, 3> allocate_elemmat(double val) const
    {
        array_type::tensor<double, 3> elemmat = xt::empty<double>(this->shape_elemmat());
        elemmat.fill(val);
        return elemmat;
    }

private:
    /**
     * Distribution to relevant implementation of \copydoc asDofs(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void asDofs_impl(const T& arg, R& ret) const
    {
        if (arg.dimension() == 2) {
            this->asDofs_impl_nodevec(arg, ret);
        }
        else if (arg.dimension() == 3) {
            this->asDofs_impl_elemvec(arg, ret);
        }
        else {
            throw std::runtime_error("Vector::asDofs unknown dimension for conversion");
        }
    }

    /**
     * Distribution to relevant implementation of \copydoc asDofs(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 2, int> = 0>
    void asDofs_impl(const T& arg, R& ret) const
    {
        this->asDofs_impl_nodevec(arg, ret);
    }

    /**
     * Distribution to relevant implementation of \copydoc asDofs(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 3, int> = 0>
    void asDofs_impl(const T& arg, R& ret) const
    {
        this->asDofs_impl_elemvec(arg, ret);
    }

    /**
     * Distribution to relevant implementation of \copydoc asNode(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void asNode_impl(const T& arg, R& ret) const
    {
        if (arg.dimension() == 1) {
            this->asNode_impl_dofval(arg, ret);
        }
        else if (arg.dimension() == 3) {
            this->asNode_impl_elemvec(arg, ret);
        }
        else {
            throw std::runtime_error("Vector::asNode unknown dimension for conversion");
        }
    }

    /**
     * Distribution to relevant implementation of \copydoc asNode(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 1, int> = 0>
    void asNode_impl(const T& arg, R& ret) const
    {
        this->asNode_impl_dofval(arg, ret);
    }

    /**
     * Distribution to relevant implementation of \copydoc asNode(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 3, int> = 0>
    void asNode_impl(const T& arg, R& ret) const
    {
        this->asNode_impl_elemvec(arg, ret);
    }

    /**
     * Distribution to relevant implementation of \copydoc asElement(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void asElement_impl(const T& arg, R& ret) const
    {
        if (arg.dimension() == 1) {
            this->asElement_impl_dofval(arg, ret);
        }
        else if (arg.dimension() == 2) {
            this->asElement_impl_nodevec(arg, ret);
        }
        else {
            throw std::runtime_error("Vector::asElement unknown dimension for conversion");
        }
    }

    /**
     * Distribution to relevant implementation of \copydoc asElement(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 1, int> = 0>
    void asElement_impl(const T& arg, R& ret) const
    {
        this->asElement_impl_dofval(arg, ret);
    }

    /**
     * Distribution to relevant implementation of \copydoc asElement(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 2, int> = 0>
    void asElement_impl(const T& arg, R& ret) const
    {
        this->asElement_impl_nodevec(arg, ret);
    }

    /**
     * Distribution to relevant implementation of \copydoc assembleDofs(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void assembleDofs_impl(const T& arg, R& ret) const
    {
        if (arg.dimension() == 2) {
            this->assembleDofs_impl_nodevec(arg, ret);
        }
        else if (arg.dimension() == 3) {
            this->assembleDofs_impl_elemvec(arg, ret);
        }
        else {
            throw std::runtime_error("Vector::assembleDofs unknown dimension for conversion");
        }
    }

    /**
     * Distribution to relevant implementation of \copydoc assembleDofs(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 2, int> = 0>
    void assembleDofs_impl(const T& arg, R& ret) const
    {
        this->assembleDofs_impl_nodevec(arg, ret);
    }

    /**
     * Distribution to relevant implementation of \copydoc assembleDofs(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 3, int> = 0>
    void assembleDofs_impl(const T& arg, R& ret) const
    {
        this->assembleDofs_impl_elemvec(arg, ret);
    }

    /**
     * Distribution to relevant implementation of \copydoc assembleNode(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<!xt::has_fixed_rank_t<T>::value, int> = 0>
    void assembleNode_impl(const T& arg, R& ret) const
    {
        if (arg.dimension() == 3) {
            this->assembleNode_impl_elemvec(arg, ret);
        }
        else {
            throw std::runtime_error("Vector::assembleNode unknown dimension for conversion");
        }
    }

    /**
     * Distribution to relevant implementation of \copydoc assembleNode(const T&, R&) const
     */
    template <class T, class R, typename std::enable_if_t<xt::get_rank<T>::value == 3, int> = 0>
    void assembleNode_impl(const T& arg, R& ret) const
    {
        this->assembleNode_impl_elemvec(arg, ret);
    }

    /**
     * Implementation for 'nodevec' input of \copydoc asDofs(const T&, R&) const
     */
    template <class T, class R>
    void asDofs_impl_nodevec(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 1 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_nodevec()));
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_dofval()));

        ret.fill(0.0);

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                ret(m_dofs(m, i)) = arg(m, i);
            }
        }
    }

    /**
     * Implementation for 'elemvec' input of \copydoc asDofs(const T&, R&) const
     */
    template <class T, class R>
    void asDofs_impl_elemvec(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 1 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_dofval()));

        ret.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    ret(m_dofs(m_conn(e, m), i)) = arg(e, m, i);
                }
            }
        }
    }

    /**
     * Implementation for 'dofval' input of \copydoc asNode(const T&, R&) const
     */
    template <class T, class R>
    void asNode_impl_dofval(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 2 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_dofval()));
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_nodevec()));

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                ret(m, i) = arg(m_dofs(m, i));
            }
        }
    }

    /**
     * Implementation for 'elemvec' input of \copydoc asNode(const T&, R&) const
     */
    template <class T, class R>
    void asNode_impl_elemvec(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 2 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_nodevec()));

        ret.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    ret(m_conn(e, m), i) = arg(e, m, i);
                }
            }
        }
    }

    /**
     * Implementation for 'dofval' input of \copydoc asElement(const T&, R&) const
     */
    template <class T, class R>
    void asElement_impl_dofval(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 3 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(arg.size() == m_ndof);
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_elemvec()));

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    ret(e, m, i) = arg(m_dofs(m_conn(e, m), i));
                }
            }
        }
    }

    /**
     * Implementation for 'nodevec' input of \copydoc asElement(const T&, R&) const
     */
    template <class T, class R>
    void asElement_impl_nodevec(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 3 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_nodevec()));
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_elemvec()));

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    ret(e, m, i) = arg(m_conn(e, m), i);
                }
            }
        }
    }

    /**
     * Implementation for 'nodevec' input of \copydoc assembleDofs(const T&, R&) const
     */
    template <class T, class R>
    void assembleDofs_impl_nodevec(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 1 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_nodevec()));
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_dofval()));

        ret.fill(0.0);

        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                ret(m_dofs(m, i)) += arg(m, i);
            }
        }
    }

    /**
     * Implementation for 'elemvec' input of \copydoc assembleDofs(const T&, R&) const
     */
    template <class T, class R>
    void assembleDofs_impl_elemvec(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 1 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_dofval()));

        ret.fill(0.0);

        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    ret(m_dofs(m_conn(e, m), i)) += arg(e, m, i);
                }
            }
        }
    }

    /**
     * Implementation for 'elemvec' input of \copydoc assembleNode(const T&, R&) const
     */
    template <class T, class R>
    void assembleNode_impl_elemvec(const T& arg, R& ret) const
    {
        static_assert(
            xt::get_rank<R>::value == 2 || !xt::has_fixed_rank_t<R>::value, "Unknown rank 'ret'");
        GOOSEFEM_ASSERT(xt::has_shape(arg, this->shape_elemvec()));
        GOOSEFEM_ASSERT(xt::has_shape(ret, this->shape_nodevec()));

        array_type::tensor<double, 1> dofval = this->AssembleDofs(arg);
        this->asNode(dofval, ret);
    }

protected:
    array_type::tensor<size_t, 2> m_conn; ///< See conn()
    array_type::tensor<size_t, 2> m_dofs; ///< See dofs()
    size_t m_nelem; ///< See #nelem
    size_t m_nne; ///< See #nne
    size_t m_nnode; ///< See #nnode
    size_t m_ndim; ///< See #ndim
    size_t m_ndof; ///< See #ndof
};

} // namespace GooseFEM

#endif

/**
 * Methods to switch between storage types based on a mesh and DOFs that are partitioned in:
 * -   unknown DOFs
 * -   prescribed DOFs
 *
 * @file VectorPartitioned.h
 * @copyright Copyright 2017. Tom de Geus. All rights reserved.
 * @license This project is released under the GNU Public License (GPLv3).
 */

#ifndef GOOSEFEM_VECTORPARTITIONED_H
#define GOOSEFEM_VECTORPARTITIONED_H

#include "Mesh.h"
#include "Vector.h"
#include "assertions.h"
#include "config.h"

namespace GooseFEM {

/**
 * Class to switch between storage types,
 * based on a mesh and DOFs that are partitioned in:
 *
 * -   unknown DOFs (iiu()), indicated with "u".
 * -   prescribed DOFs (iip()), indicated with "p".
 *
 * To this end some internal re-ordering of the DOFs has to be done, as follows:
 *
 *     iiu() -> arange(nnu())
 *     iip() -> nnu() + arange(nnp())
 *
 * which is relevant only if you interact using partitioned DOF-lists ("dofval_u" or "dofval_p").
 *
 * The "dofval", "nodevec", and "elemvec" are all stored in the 'normal' order.
 *
 * For reference:
 *
 * -   "dofval": DOF values [#ndof].
 * -   "dofval_u": unknown DOF values, `== dofval[iiu()]`, [#nnu].
 * -   "dofval_p": prescribed DOF values, `== dofval[iip()]`, [#nnp].
 * -   "nodevec": nodal vectors [#nnode, #ndim].
 * -   "elemvec": nodal vectors stored per element [#nelem, #nne, #ndim].
 */
class VectorPartitioned : public Vector {
protected:
    array_type::tensor<size_t, 1> m_iiu; ///< See iiu()
    array_type::tensor<size_t, 1> m_iip; ///< See iip()
    size_t m_nnu; ///< See #nnu
    size_t m_nnp; ///< See #nnp

    /**
     * Renumbered DOFs per node, such that
     *
     *     iiu = arange(nnu)
     *     iip = nnu + arange(nnp)
     *
     * making is much simpler to slice.
     */
    array_type::tensor<size_t, 2> m_part;

public:
    VectorPartitioned() = default;

    /**
     * Constructor.
     *
     * @param conn connectivity [#nelem, #nne].
     * @param dofs DOFs per node [#nnode, #ndim].
     * @param iip prescribed DOFs [#nnp].
     */
    VectorPartitioned(
        const array_type::tensor<size_t, 2>& conn,
        const array_type::tensor<size_t, 2>& dofs,
        const array_type::tensor<size_t, 1>& iip)
        : Vector(conn, dofs), m_iip(iip)
    {
        GOOSEFEM_ASSERT(is_unique(iip));

        m_iiu = xt::setdiff1d(m_dofs, m_iip);
        m_nnp = m_iip.size();
        m_nnu = m_iiu.size();
        m_part = Mesh::Reorder({m_iiu, m_iip}).apply(m_dofs);

        GOOSEFEM_ASSERT(xt::amax(m_iip)() <= xt::amax(m_dofs)());
    }

    /**
     * @return Number of unknown DOFs.
     */
    size_t nnu() const
    {
        return m_nnu;
    }

    /**
     * @return Number of prescribed DOFs.
     */
    size_t nnp() const
    {
        return m_nnp;
    }

    /**
     * @return Unknown DOFs [#nnu].
     */
    const array_type::tensor<size_t, 1>& iiu() const
    {
        return m_iiu;
    }

    /**
     * @return Prescribed DOFs [#nnp].
     */
    const array_type::tensor<size_t, 1>& iip() const
    {
        return m_iip;
    }

    /**
     * Per DOF (see Vector::dofs()) list if unknown ("u").
     *
     * @return Boolean "nodevec".
     */
    array_type::tensor<bool, 2> dofs_is_u() const
    {
        array_type::tensor<bool, 2> ret = xt::zeros<bool>(this->shape_nodevec());

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) < m_nnu) {
                    ret(m, i) = true;
                }
            }
        }

        return ret;
    }

    /**
     * Per DOF (see Vector::dofs()) list if prescribed ("p").
     *
     * @return Boolean "nodevec".
     */
    array_type::tensor<bool, 2> dofs_is_p() const
    {
        array_type::tensor<bool, 2> ret = xt::zeros<bool>(this->shape_nodevec());

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) >= m_nnu) {
                    ret(m, i) = true;
                }
            }
        }

        return ret;
    }

    /**
     * Copy unknown DOFs from "nodevec" to another "nodvec":
     *
     *     nodevec_dest[vector.dofs_is_u()] = nodevec_src
     *
     * the other DOFs are taken from `nodevec_dest`:
     *
     *     nodevec_dest[vector.dofs_is_p()] = nodevec_dest
     *
     * @param nodevec_src input [#nnode, #ndim]
     * @param nodevec_dest input [#nnode, #ndim]
     * @return nodevec output [#nnode, #ndim]
     */
    array_type::tensor<double, 2> Copy_u(
        const array_type::tensor<double, 2>& nodevec_src,
        const array_type::tensor<double, 2>& nodevec_dest) const
    {
        array_type::tensor<double, 2> ret = nodevec_dest;
        this->copy_u(nodevec_src, ret);
        return ret;
    }

    /**
     * Copy unknown DOFs from "nodevec" to another "nodvec":
     *
     *     nodevec_dest[vector.dofs_is_u()] = nodevec_src
     *
     * the other DOFs are taken from `nodevec_dest`:
     *
     *     nodevec_dest[vector.dofs_is_p()] = nodevec_dest
     *
     * @param nodevec_src input [#nnode, #ndim]
     * @param nodevec_dest input/output [#nnode, #ndim]
     */
    void copy_u(
        const array_type::tensor<double, 2>& nodevec_src,
        array_type::tensor<double, 2>& nodevec_dest) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(nodevec_src, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(nodevec_dest, {m_nnode, m_ndim}));

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) < m_nnu) {
                    nodevec_dest(m, i) = nodevec_src(m, i);
                }
            }
        }
    }

    /**
     * Copy prescribed DOFs from "nodevec" to another "nodvec":
     *
     *     nodevec_dest[vector.dofs_is_p()] = nodevec_src
     *
     * the other DOFs are taken from `nodevec_dest`:
     *
     *     nodevec_dest[vector.dofs_is_u()] = nodevec_dest
     *
     * @param nodevec_src input [#nnode, #ndim]
     * @param nodevec_dest input [#nnode, #ndim]
     * @return nodevec output [#nnode, #ndim]
     */
    array_type::tensor<double, 2> Copy_p(
        const array_type::tensor<double, 2>& nodevec_src,
        const array_type::tensor<double, 2>& nodevec_dest) const
    {
        array_type::tensor<double, 2> ret = nodevec_dest;
        this->copy_p(nodevec_src, ret);
        return ret;
    }

    /**
     * Copy prescribed DOFs from "nodevec" to another "nodvec":
     *
     *     nodevec_dest[vector.dofs_is_p()] = nodevec_src
     *
     * the other DOFs are taken from `nodevec_dest`:
     *
     *     nodevec_dest[vector.dofs_is_u()] = nodevec_dest
     *
     * @param nodevec_src input [#nnode, #ndim]
     * @param nodevec_dest input/output [#nnode, #ndim]
     */
    void copy_p(
        const array_type::tensor<double, 2>& nodevec_src,
        array_type::tensor<double, 2>& nodevec_dest) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(nodevec_src, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(xt::has_shape(nodevec_dest, {m_nnode, m_ndim}));

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) >= m_nnu) {
                    nodevec_dest(m, i) = nodevec_src(m, i);
                }
            }
        }
    }

    /**
     * Combine unknown and prescribed "dofval" into a single "dofval" list.
     *
     * @param dofval_u input [#nnu]
     * @param dofval_p input [#nnp]
     * @return dofval output [#ndof]
     */
    array_type::tensor<double, 1> DofsFromParitioned(
        const array_type::tensor<double, 1>& dofval_u,
        const array_type::tensor<double, 1>& dofval_p) const
    {
        array_type::tensor<double, 1> dofval = xt::empty<double>({m_ndof});
        this->dofsFromParitioned(dofval_u, dofval_p, dofval);
        return dofval;
    }

    /**
     * Combine unknown and prescribed "dofval" into a single "dofval" list.
     *
     * @param dofval_u input [#nnu]
     * @param dofval_p input [#nnp]
     * @param dofval output [#ndof]
     */
    void dofsFromParitioned(
        const array_type::tensor<double, 1>& dofval_u,
        const array_type::tensor<double, 1>& dofval_p,
        array_type::tensor<double, 1>& dofval) const
    {
        GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);
        GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);
        GOOSEFEM_ASSERT(dofval.size() == m_ndof);

        dofval.fill(0.0);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            dofval(m_iiu(d)) = dofval_u(d);
        }

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            dofval(m_iip(d)) = dofval_p(d);
        }
    }

    /**
     * Combine unknown and prescribed "dofval" into a single "dofval" list
     * and directly convert to "nodeval" without a temporary
     * (overwrite entries that occur more than once).
     *
     * @param dofval_u input [#nnu]
     * @param dofval_p input [#nnp]
     * @return nodevec output [#nnode, #ndim]
     */
    array_type::tensor<double, 2> NodeFromPartitioned(
        const array_type::tensor<double, 1>& dofval_u,
        const array_type::tensor<double, 1>& dofval_p) const
    {
        array_type::tensor<double, 2> nodevec = xt::empty<double>({m_nnode, m_ndim});
        this->nodeFromPartitioned(dofval_u, dofval_p, nodevec);
        return nodevec;
    }

    /**
     * Combine unknown and prescribed "dofval" into a single "dofval" list
     * and directly convert to "nodeval" without a temporary
     * (overwrite entries that occur more than once).
     *
     * @param dofval_u input [#nnu]
     * @param dofval_p input [#nnp]
     * @param nodevec output [#nnode, #ndim]
     */
    void nodeFromPartitioned(
        const array_type::tensor<double, 1>& dofval_u,
        const array_type::tensor<double, 1>& dofval_p,
        array_type::tensor<double, 2>& nodevec) const
    {
        GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);
        GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);
        GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) < m_nnu) {
                    nodevec(m, i) = dofval_u(m_part(m, i));
                }
                else {
                    nodevec(m, i) = dofval_p(m_part(m, i) - m_nnu);
                }
            }
        }
    }

    /**
     * Combine unknown and prescribed "dofval" into a single "dofval" list
     * and directly convert to "elemvec" without a temporary
     * (overwrite entries that occur more than once).
     *
     * @param dofval_u input [#nnu]
     * @param dofval_p input [#nnp]
     * @return elemvec output [#nelem, #nne, #ndim]
     */
    array_type::tensor<double, 3> ElementFromPartitioned(
        const array_type::tensor<double, 1>& dofval_u,
        const array_type::tensor<double, 1>& dofval_p) const
    {
        array_type::tensor<double, 3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});
        this->elementFromPartitioned(dofval_u, dofval_p, elemvec);
        return elemvec;
    }

    /**
     * Combine unknown and prescribed "dofval" into a single "dofval" list
     * and directly convert to "elemvec" without a temporary
     * (overwrite entries that occur more than once).
     *
     * @param dofval_u input [#nnu]
     * @param dofval_p input [#nnp]
     * @param elemvec output [#nelem, #nne, #ndim]
     */
    void elementFromPartitioned(
        const array_type::tensor<double, 1>& dofval_u,
        const array_type::tensor<double, 1>& dofval_p,
        array_type::tensor<double, 3>& elemvec) const
    {
        GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);
        GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    if (m_part(m_conn(e, m), i) < m_nnu) {
                        elemvec(e, m, i) = dofval_u(m_part(m_conn(e, m), i));
                    }
                    else {
                        elemvec(e, m, i) = dofval_p(m_part(m_conn(e, m), i) - m_nnu);
                    }
                }
            }
        }
    }

    /**
     * Extract the unknown "dofval":
     *
     *     dofval[iiu()]
     *
     * @param dofval input [#ndof]
     * @return dofval_u input [#nnu]
     */
    array_type::tensor<double, 1> AsDofs_u(const array_type::tensor<double, 1>& dofval) const
    {
        array_type::tensor<double, 1> dofval_u = xt::empty<double>({m_nnu});
        this->asDofs_u(dofval, dofval_u);
        return dofval_u;
    }

    /**
     * Extract the unknown "dofval":
     *
     *     dofval[iiu()]
     *
     * @param dofval input [#ndof]
     * @param dofval_u input [#nnu]
     */
    void asDofs_u(
        const array_type::tensor<double, 1>& dofval,
        array_type::tensor<double, 1>& dofval_u) const
    {
        GOOSEFEM_ASSERT(dofval.size() == m_ndof);
        GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnu; ++d) {
            dofval_u(d) = dofval(m_iiu(d));
        }
    }

    /**
     * Convert "nodevec" to "dofval" (overwrite entries that occur more than once)
     * and extract the unknown "dofval" without a temporary.
     *
     * @param nodevec input [#nnode, #ndim]
     * @return dofval_u input [#nnu]
     */
    array_type::tensor<double, 1> AsDofs_u(const array_type::tensor<double, 2>& nodevec) const
    {
        array_type::tensor<double, 1> dofval_u = xt::empty<double>({m_nnu});
        this->asDofs_u(nodevec, dofval_u);
        return dofval_u;
    }

    /**
     * Convert "nodevec" to "dofval" (overwrite entries that occur more than once)
     * and extract the unknown "dofval" without a temporary.
     *
     * @param nodevec input [#nnode, #ndim]
     * @param dofval_u input [#nnu]
     */
    void asDofs_u(
        const array_type::tensor<double, 2>& nodevec,
        array_type::tensor<double, 1>& dofval_u) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);

        dofval_u.fill(0.0);

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) < m_nnu) {
                    dofval_u(m_part(m, i)) = nodevec(m, i);
                }
            }
        }
    }

    /**
     * Convert "elemvec" to "dofval" (overwrite entries that occur more than once)
     * and extract the unknown "dofval" without a temporary.
     *
     * @param elemvec input [#nelem, #nne, #ndim]
     * @return dofval_u input [#nnu]
     */
    array_type::tensor<double, 1> AsDofs_u(const array_type::tensor<double, 3>& elemvec) const
    {
        array_type::tensor<double, 1> dofval_u = xt::empty<double>({m_nnu});
        this->asDofs_u(elemvec, dofval_u);
        return dofval_u;
    }

    /**
     * Convert "elemvec" to "dofval" (overwrite entries that occur more than once)
     * and extract the unknown "dofval" without a temporary.
     *
     * @param elemvec input [#nelem, #nne, #ndim]
     * @param dofval_u input [#nnu]
     */
    void asDofs_u(
        const array_type::tensor<double, 3>& elemvec,
        array_type::tensor<double, 1>& dofval_u) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));
        GOOSEFEM_ASSERT(dofval_u.size() == m_nnu);

        dofval_u.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    if (m_part(m_conn(e, m), i) < m_nnu) {
                        dofval_u(m_part(m_conn(e, m), i)) = elemvec(e, m, i);
                    }
                }
            }
        }
    }

    /**
     * Extract the prescribed "dofval":
     *
     *     dofval[iip()]
     *
     * @param dofval input [#ndof]
     * @return dofval_p input [#nnp]
     */
    array_type::tensor<double, 1> AsDofs_p(const array_type::tensor<double, 1>& dofval) const
    {
        array_type::tensor<double, 1> dofval_p = xt::empty<double>({m_nnp});
        this->asDofs_p(dofval, dofval_p);
        return dofval_p;
    }

    /**
     * Extract the prescribed "dofval":
     *
     *     dofval[iip()]
     *
     * @param dofval input [#ndof]
     * @param dofval_p input [#nnp]
     */
    void asDofs_p(
        const array_type::tensor<double, 1>& dofval,
        array_type::tensor<double, 1>& dofval_p) const
    {
        GOOSEFEM_ASSERT(dofval.size() == m_ndof);
        GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);

#pragma omp parallel for
        for (size_t d = 0; d < m_nnp; ++d) {
            dofval_p(d) = dofval(m_iip(d));
        }
    }

    /**
     * Convert "nodevec" to "dofval" (overwrite entries that occur more than once)
     * and extract the prescribed "dofval" without a temporary.
     *
     * @param nodevec input [#nnode, #ndim]
     * @return dofval_p input [#nnp]
     */
    array_type::tensor<double, 1> AsDofs_p(const array_type::tensor<double, 2>& nodevec) const
    {
        array_type::tensor<double, 1> dofval_p = xt::empty<double>({m_nnp});
        this->asDofs_p(nodevec, dofval_p);
        return dofval_p;
    }

    /**
     * Convert "nodevec" to "dofval" (overwrite entries that occur more than once)
     * and extract the prescribed "dofval" without a temporary.
     *
     * @param nodevec input [#nnode, #ndim]
     * @param dofval_p input [#nnp]
     */
    void asDofs_p(
        const array_type::tensor<double, 2>& nodevec,
        array_type::tensor<double, 1>& dofval_p) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(nodevec, {m_nnode, m_ndim}));
        GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);

        dofval_p.fill(0.0);

#pragma omp parallel for
        for (size_t m = 0; m < m_nnode; ++m) {
            for (size_t i = 0; i < m_ndim; ++i) {
                if (m_part(m, i) >= m_nnu) {
                    dofval_p(m_part(m, i) - m_nnu) = nodevec(m, i);
                }
            }
        }
    }

    /**
     * Convert "elemvec" to "dofval" (overwrite entries that occur more than once)
     * and extract the prescribed "dofval" without a temporary.
     *
     * @param elemvec input [#nelem, #nne, #ndim]
     * @return dofval_p input [#nnp]
     */
    array_type::tensor<double, 1> AsDofs_p(const array_type::tensor<double, 3>& elemvec) const
    {
        array_type::tensor<double, 1> dofval_p = xt::empty<double>({m_nnp});
        this->asDofs_p(elemvec, dofval_p);
        return dofval_p;
    }

    /**
     * Convert "elemvec" to "dofval" (overwrite entries that occur more than once)
     * and extract the prescribed "dofval" without a temporary.
     *
     * @param elemvec input [#nelem, #nne, #ndim]
     * @param dofval_p input [#nnp]
     */
    void asDofs_p(
        const array_type::tensor<double, 3>& elemvec,
        array_type::tensor<double, 1>& dofval_p) const
    {
        GOOSEFEM_ASSERT(xt::has_shape(elemvec, {m_nelem, m_nne, m_ndim}));
        GOOSEFEM_ASSERT(dofval_p.size() == m_nnp);

        dofval_p.fill(0.0);

#pragma omp parallel for
        for (size_t e = 0; e < m_nelem; ++e) {
            for (size_t m = 0; m < m_nne; ++m) {
                for (size_t i = 0; i < m_ndim; ++i) {
                    if (m_part(m_conn(e, m), i) >= m_nnu) {
                        dofval_p(m_part(m_conn(e, m), i) - m_nnu) = elemvec(e, m, i);
                    }
                }
            }
        }
    }
};

} // namespace GooseFEM

#endif

/**
Methods to switch between storage types based on a mesh and DOFs that are partitioned in:
-   unknown DOFs
-   prescribed DOFs

\file VectorPartitioned.h
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_VECTORPARTITIONED_H
#define GOOSEFEM_VECTORPARTITIONED_H

#include "Vector.h"
#include "config.h"

namespace GooseFEM {

/**
 *  Class to switch between storage types,
 *  based on a mesh and DOFs that are partitioned in:
 *
 *  -   unknown DOFs (iiu()), indicated with "u".
 *  -   prescribed DOFs (iip()), indicated with "p".
 *
 *  To this end some internal re-ordering of the DOFs has to be done, as follows:
 *
 *      iiu() -> arange(nnu())
 *      iip() -> nnu() + arange(nnp())
 *
 *  which is relevant only if you interact using partitioned DOF-lists ("dofval_u" or "dofval_p").
 *
 *  The "dofval", "nodevec", and "elemvec" are all stored in the 'normal' order.
 *
 *  For reference:
 *
 *  -   "dofval": DOF values [#ndof].
 *  -   "dofval_u": unknown DOF values, `== dofval[iiu()]`, [#nnu].
 *  -   "dofval_p": prescribed DOF values, `== dofval[iip()]`, [#nnp].
 *  -   "nodevec": nodal vectors [#nnode, #ndim].
 *  -   "elemvec": nodal vectors stored per element [#nelem, #nne, #ndim].
 *
 */
class VectorPartitioned : public Vector {
public:
    VectorPartitioned() = default;

    /**
    Constructor.

    \param conn connectivity [#nelem, #nne].
    \param dofs DOFs per node [#nnode, #ndim].
    \param iip prescribed DOFs [#nnp].
    */
    VectorPartitioned(
        const xt::xtensor<size_t, 2>& conn,
        const xt::xtensor<size_t, 2>& dofs,
        const xt::xtensor<size_t, 1>& iip);

    /**
    \return Number of unknown DOFs.
    */
    size_t nnu() const;

    /**
    \return Number of prescribed DOFs.
    */
    size_t nnp() const;

    /**
    \return Unknown DOFs [#nnu].
    */
    xt::xtensor<size_t, 1> iiu() const;

    /**
    \return Prescribed DOFs [#nnp].
    */
    xt::xtensor<size_t, 1> iip() const;

    /**
    Per DOF (see Vector::dofs()) list if unknown ("u").

    \return Boolean "nodevec".
    */
    xt::xtensor<bool, 2> dofs_is_u() const;

    /**
    Per DOF (see Vector::dofs()) list if prescribed ("p").

    \return Boolean "nodevec".
    */
    xt::xtensor<bool, 2> dofs_is_p() const;

    /**
    Copy unknown DOFs from "nodevec" to another "nodvec":

        nodevec_dest[vector.dofs_is_u()] = nodevec_src

    the other DOFs are taken from ``nodevec_dest``:

        nodevec_dest[vector.dofs_is_p()] = nodevec_dest

    \param nodevec_src input [#nnode, #ndim]
    \param nodevec_dest input [#nnode, #ndim]
    \return nodevec output [#nnode, #ndim]
    */
    xt::xtensor<double, 2> Copy_u(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

    /**
    Copy unknown DOFs from "nodevec" to another "nodvec":

        nodevec_dest[vector.dofs_is_u()] = nodevec_src

    the other DOFs are taken from ``nodevec_dest``:

        nodevec_dest[vector.dofs_is_p()] = nodevec_dest

    \param nodevec_src input [#nnode, #ndim]
    \param nodevec_dest input/output [#nnode, #ndim]
    */
    void
    copy_u(const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const;

    /**
    Copy prescribed DOFs from "nodevec" to another "nodvec":

        nodevec_dest[vector.dofs_is_p()] = nodevec_src

    the other DOFs are taken from ``nodevec_dest``:

        nodevec_dest[vector.dofs_is_u()] = nodevec_dest

    \param nodevec_src input [#nnode, #ndim]
    \param nodevec_dest input [#nnode, #ndim]
    \return nodevec output [#nnode, #ndim]
    */
    xt::xtensor<double, 2> Copy_p(
        const xt::xtensor<double, 2>& nodevec_src,
        const xt::xtensor<double, 2>& nodevec_dest) const;

    /**
    Copy prescribed DOFs from "nodevec" to another "nodvec":

        nodevec_dest[vector.dofs_is_p()] = nodevec_src

    the other DOFs are taken from ``nodevec_dest``:

        nodevec_dest[vector.dofs_is_u()] = nodevec_dest

    \param nodevec_src input [#nnode, #ndim]
    \param nodevec_dest input/output [#nnode, #ndim]
    */
    void
    copy_p(const xt::xtensor<double, 2>& nodevec_src, xt::xtensor<double, 2>& nodevec_dest) const;

    /**
    Combine unknown and prescribed "dofval" into a single "dofval" list.

    \param dofval_u input [#nnu]
    \param dofval_p input [#nnp]
    \return dofval output [#ndof]
    */
    xt::xtensor<double, 1> DofsFromParitioned(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p) const;

    /**
    Combine unknown and prescribed "dofval" into a single "dofval" list.

    \param dofval_u input [#nnu]
    \param dofval_p input [#nnp]
    \param dofval output [#ndof]
    */
    void dofsFromParitioned(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p,
        xt::xtensor<double, 1>& dofval) const;

    /**
    Combine unknown and prescribed "dofval" into a single "dofval" list
    and directly convert to "nodeval" without a temporary
    (overwrite entries that occur more than once).

    \param dofval_u input [#nnu]
    \param dofval_p input [#nnp]
    \return nodevec output [#nnode, #ndim]
    */
    xt::xtensor<double, 2> NodeFromPartitioned(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p) const;

    /**
    Combine unknown and prescribed "dofval" into a single "dofval" list
    and directly convert to "nodeval" without a temporary
    (overwrite entries that occur more than once).

    \param dofval_u input [#nnu]
    \param dofval_p input [#nnp]
    \param nodevec output [#nnode, #ndim]
    */
    void nodeFromPartitioned(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p,
        xt::xtensor<double, 2>& nodevec) const;

    /**
    Combine unknown and prescribed "dofval" into a single "dofval" list
    and directly convert to "elemvec" without a temporary
    (overwrite entries that occur more than once).

    \param dofval_u input [#nnu]
    \param dofval_p input [#nnp]
    \return elemvec output [#nelem, #nne, #ndim]
    */
    xt::xtensor<double, 3> ElementFromPartitioned(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p) const;

    /**
    Combine unknown and prescribed "dofval" into a single "dofval" list
    and directly convert to "elemvec" without a temporary
    (overwrite entries that occur more than once).

    \param dofval_u input [#nnu]
    \param dofval_p input [#nnp]
    \param elemvec output [#nelem, #nne, #ndim]
    */
    void elementFromPartitioned(
        const xt::xtensor<double, 1>& dofval_u,
        const xt::xtensor<double, 1>& dofval_p,
        xt::xtensor<double, 3>& elemvec) const;

    /**
    Extract the unknown "dofval":

        dofval[iiu()]

    \param dofval input [#ndof]
    \return dofval_u input [#nnu]
    */
    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 1>& dofval) const;

    /**
    Extract the unknown "dofval":

        dofval[iiu()]

    \param dofval input [#ndof]
    \param dofval_u input [#nnu]
    */
    void asDofs_u(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 1>& dofval_u) const;

    /**
    Convert "nodevec" to "dofval" (overwrite entries that occur more than once)
    and extract the unknown "dofval" without a temporary.

    \param nodevec input [#nnode, #ndim]
    \return dofval_u input [#nnu]
    */
    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 2>& nodevec) const;

    /**
    Convert "nodevec" to "dofval" (overwrite entries that occur more than once)
    and extract the unknown "dofval" without a temporary.

    \param nodevec input [#nnode, #ndim]
    \param dofval_u input [#nnu]
    */
    void asDofs_u(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval_u) const;

    /**
    Convert "elemvec" to "dofval" (overwrite entries that occur more than once)
    and extract the unknown "dofval" without a temporary.

    \param elemvec input [#nelem, #nne, #ndim]
    \return dofval_u input [#nnu]
    */
    xt::xtensor<double, 1> AsDofs_u(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Convert "elemvec" to "dofval" (overwrite entries that occur more than once)
    and extract the unknown "dofval" without a temporary.

    \param elemvec input [#nelem, #nne, #ndim]
    \param dofval_u input [#nnu]
    */
    void asDofs_u(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval_u) const;

    /**
    Extract the prescribed "dofval":

        dofval[iip()]

    \param dofval input [#ndof]
    \return dofval_p input [#nnp]
    */
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 1>& dofval) const;

    /**
    Extract the prescribed "dofval":

        dofval[iip()]

    \param dofval input [#ndof]
    \param dofval_p input [#nnp]
    */
    void asDofs_p(const xt::xtensor<double, 1>& dofval, xt::xtensor<double, 1>& dofval_p) const;

    /**
    Convert "nodevec" to "dofval" (overwrite entries that occur more than once)
    and extract the prescribed "dofval" without a temporary.

    \param nodevec input [#nnode, #ndim]
    \return dofval_p input [#nnp]
    */
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 2>& nodevec) const;

    /**
    Convert "nodevec" to "dofval" (overwrite entries that occur more than once)
    and extract the prescribed "dofval" without a temporary.

    \param nodevec input [#nnode, #ndim]
    \param dofval_p input [#nnp]
    */
    void asDofs_p(const xt::xtensor<double, 2>& nodevec, xt::xtensor<double, 1>& dofval_p) const;

    /**
    Convert "elemvec" to "dofval" (overwrite entries that occur more than once)
    and extract the prescribed "dofval" without a temporary.

    \param elemvec input [#nelem, #nne, #ndim]
    \return dofval_p input [#nnp]
    */
    xt::xtensor<double, 1> AsDofs_p(const xt::xtensor<double, 3>& elemvec) const;

    /**
    Convert "elemvec" to "dofval" (overwrite entries that occur more than once)
    and extract the prescribed "dofval" without a temporary.

    \param elemvec input [#nelem, #nne, #ndim]
    \param dofval_p input [#nnp]
    */
    void asDofs_p(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 1>& dofval_p) const;

protected:
    xt::xtensor<size_t, 1> m_iiu; ///< See iiu()
    xt::xtensor<size_t, 1> m_iip; ///< See iip()
    size_t m_nnu; ///< See #nnu
    size_t m_nnp; ///< See #nnp

    /**
    Renumbered DOFs per node, such that

        iiu = arange(nnu)
        iip = nnu + arange(nnp)

    making is much simpler to slice.
    */
    xt::xtensor<size_t, 2> m_part;
};

} // namespace GooseFEM

#include "VectorPartitioned.hpp"

#endif

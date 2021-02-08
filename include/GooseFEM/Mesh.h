/**
\file Mesh.h
Generic mesh operations.

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM
*/

#ifndef GOOSEFEM_MESH_H
#define GOOSEFEM_MESH_H

#include "config.h"

namespace GooseFEM {
namespace Mesh {


/**
Enumerator for element-types
*/
enum class ElementType {
    Quad4, /*!< Quadrilateral: 4-noded element in 2-d */
    Hex8, /*!< Hexahedron: 8-noded element in 3-d */
    Tri3 /*!< Triangle: 3-noded element in 2-d */
};

/**
Extract the element type based on the connectivity.

\param coor Nodal coordinates.
\param conn Connectivity.
\return ElementType
*/
inline ElementType defaultElementType(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn)
{
    if (coor.shape(1) == 2ul && conn.shape(1) == 3ul) {
        return ElementType::Tri3;
    }
    if (coor.shape(1) == 2ul && conn.shape(1) == 4ul) {
        return ElementType::Quad4;
    }
    if (coor.shape(1) == 3ul && conn.shape(1) == 8ul) {
        return ElementType::Hex8;
    }

    throw std::runtime_error("Element-type not implemented");
}

namespace detail {

    template <class T, class R>
    inline T renum(const T& arg, const R& mapping)
    {
        T ret = T::from_shape(arg.shape());

        auto jt = ret.begin();

        for (auto it = arg.begin(); it != arg.end(); ++it, ++jt) {
            *jt = mapping(*it);
        }

        return ret;
    }

} // namespace detail

/**
Find overlapping nodes. The output has the following structure:

    [[nodes_from_mesh_a],
     [nodes_from_mesh_b]]

\param coor_a Nodal coordinates of mesh "a".
\param coor_b Nodal coordinates of mesh "b".
\param rtol Relative tolerance for position match.
\param atol Absolute tolerance for position match.
\return Overlapping nodes.
*/
inline xt::xtensor<size_t, 2> overlapping(
    const xt::xtensor<double, 2>& coor_a,
    const xt::xtensor<double, 2>& coor_b,
    double rtol = 1e-5,
    double atol = 1e-8)
{
    GOOSEFEM_ASSERT(coor_a.shape(1) == coor_b.shape(1));

    std::vector<size_t> ret_a;
    std::vector<size_t> ret_b;

    for (size_t i = 0; i < coor_a.shape(0); ++i) {

        auto idx = xt::flatten_indices(xt::argwhere(xt::prod(xt::isclose(
            coor_b, xt::view(coor_a, i, xt::all()), rtol, atol), 1)));

        for (auto& j : idx) {
            ret_a.push_back(i);
            ret_b.push_back(j);
        }
    }

    xt::xtensor<size_t, 2> ret = xt::empty<size_t>({size_t(2), ret_a.size()});
    for (size_t i = 0; i < ret_a.size(); ++i) {
        ret(0, i) = ret_a[i];
        ret(1, i) = ret_b[i];
    }

    return ret;
}

/**
Stitch two mesh objects, specifying overlapping nodes by hand.
*/
class ManualStitch {
private:
    xt::xtensor<double, 2> m_coor;
    xt::xtensor<size_t, 2> m_conn;
    xt::xtensor<size_t, 1> m_map_b;
    size_t m_nnd_a;
    size_t m_nel_a;
    size_t m_nel_b;

public:
    ManualStitch() = default;

    /**
    \param coor_a Nodal coordinates of mesh "a".
    \param conn_a Connectivity of mesh "a".
    \param overlapping_nodes_a Node-numbers of mesh "a" that overlap with mesh "b".
    \param coor_b Nodal coordinates of mesh "b".
    \param conn_b Connectivity of mesh "b".
    \param overlapping_nodes_b Node-numbers of mesh "b" that overlap with mesh "a".
    \param check_position If ``true`` the nodes are checked for position overlap.
    \param rtol Relative tolerance for check on position overlap.
    \param atol Absolute tolerance for check on position overlap.
    */
    ManualStitch(
        const xt::xtensor<double, 2>& coor_a,
        const xt::xtensor<size_t, 2>& conn_a,
        const xt::xtensor<size_t, 1>& overlapping_nodes_a,
        const xt::xtensor<double, 2>& coor_b,
        const xt::xtensor<size_t, 2>& conn_b,
        const xt::xtensor<size_t, 1>& overlapping_nodes_b,
        bool check_position = true,
        double rtol = 1e-5,
        double atol = 1e-8)
    {
        UNUSED(rtol);
        UNUSED(atol);

        GOOSEFEM_ASSERT(xt::has_shape(overlapping_nodes_a, overlapping_nodes_b.shape()));
        GOOSEFEM_ASSERT(coor_a.shape(1) == coor_b.shape(1));
        GOOSEFEM_ASSERT(conn_a.shape(1) == conn_b.shape(1));

        if (check_position) {
            GOOSEFEM_ASSERT(xt::allclose(
                xt::view(coor_a, xt::keep(overlapping_nodes_a), xt::all()),
                xt::view(coor_b, xt::keep(overlapping_nodes_b), xt::all()),
                rtol,
                atol));
        }

        size_t nnda = coor_a.shape(0);
        size_t nndb = coor_b.shape(0);
        size_t ndim = coor_a.shape(1);
        size_t nelim = overlapping_nodes_a.size();

        size_t nela = conn_a.shape(0);
        size_t nelb = conn_b.shape(0);
        size_t nne = conn_a.shape(1);

        m_nel_a = nela;
        m_nel_b = nelb;
        m_nnd_a = nnda;

        xt::xtensor<size_t, 1> keep_b = xt::setdiff1d(xt::arange<size_t>(nndb), overlapping_nodes_b);

        m_map_b = xt::empty<size_t>({nndb});
        xt::view(m_map_b, xt::keep(overlapping_nodes_b)) = overlapping_nodes_a;
        xt::view(m_map_b, xt::keep(keep_b)) = xt::arange<size_t>(keep_b.size()) + nnda;

        m_conn = xt::empty<size_t>({nela + nelb, nne});
        xt::view(m_conn, xt::range(0, nela), xt::all()) = conn_a;
        xt::view(m_conn, xt::range(nela, nela + nelb), xt::all()) = detail::renum(conn_b, m_map_b);

        m_coor = xt::empty<size_t>({nnda + nndb - nelim, ndim});
        xt::view(m_coor, xt::range(0, nnda), xt::all()) = coor_a;
        xt::view(m_coor, xt::range(nnda, nnda + nndb - nelim), xt::all()) =
            xt::view(coor_b, xt::keep(keep_b), xt::all());
    }

    /**
    \return Nodal coordinates of stitched mesh.
    */
    xt::xtensor<double, 2> coor() const
    {
        return m_coor;
    }

    /**
    \return Connectivity of stitched mesh.
    */
    xt::xtensor<size_t, 2> conn() const
    {
        return m_conn;
    }

    /**
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return Node-map for a given mesh.
    */
    xt::xtensor<size_t, 1> nodemap(size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index <= 1);

        if (mesh_index == 0) {
            return xt::arange<size_t>(m_nnd_a);
        }

        return m_map_b;
    }

    /**
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return Element-map for a given mesh.
    */
    xt::xtensor<size_t, 1> elemmap(size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index <= 1);

        if (mesh_index == 0) {
            return xt::arange<size_t>(m_nel_a);
        }

        return xt::arange<size_t>(m_nel_b) + m_nel_a;
    }

    /**
    Convert set of node-numbers for an original mesh to the stitched mesh.

    \param set Set of node-numbers.
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return Set of node-number for the stitched mesh.
    */
    xt::xtensor<size_t, 1> nodeset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index <= 1);

        if (mesh_index == 0) {
            GOOSEFEM_ASSERT(xt::amax(set)() < m_nnd_a);
            return set;
        }

        GOOSEFEM_ASSERT(xt::amax(set)() < m_map_b.size());
        return detail::renum(set, m_map_b);
    }

    /**
    Convert set of element-numbers for an original mesh to the stitched mesh.

    \param set Set of element-numbers.
    \param mesh_index Index of the mesh ("a" = 1, "b" = 1).
    \return Set of element-number for the stitched mesh.
    */
    xt::xtensor<size_t, 1> elemset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index <= 1);

        if (mesh_index == 0) {
            GOOSEFEM_ASSERT(xt::amax(set)() < m_nel_a);
            return set;
        }

        GOOSEFEM_ASSERT(xt::amax(set)() < m_nel_b);
        return set + m_nel_a;
    }
};

/**
Stitch mesh objects, automatically searching for overlapping nodes.
*/
class Stitch {
private:
    xt::xtensor<double, 2> m_coor;
    xt::xtensor<size_t, 2> m_conn;
    std::vector<xt::xtensor<size_t, 1>> m_map;
    std::vector<size_t> m_nel;
    std::vector<size_t> m_el_offset;
    double m_rtol;
    double m_atol;

public:
    /**
    \param rtol Relative tolerance for position match.
    \param atol Absolute tolerance for position match.
    */
    Stitch(double rtol = 1e-5, double atol = 1e-8)
    {
        m_rtol = rtol;
        m_atol = atol;
    }

    /**
    Add mesh to be stitched.

    \param coor Nodal coordinates.
    \param conn Connectivity.
    */
    void push_back(const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn)
    {
        if (m_map.size() == 0) {
            m_coor = coor;
            m_conn = conn;
            m_map.push_back(xt::eval(xt::arange<size_t>(coor.shape(0))));
            m_nel.push_back(conn.shape(0));
            m_el_offset.push_back(0);
            return;
        }

        auto overlap = overlapping(m_coor, coor, m_rtol, m_atol);
        size_t index = m_map.size();

        ManualStitch stich(
            m_coor, m_conn, xt::view(overlap, 0, xt::all()),
            coor, conn, xt::view(overlap, 1, xt::all()),
            false);

        m_coor = stich.coor();
        m_conn = stich.conn();
        m_map.push_back(stich.nodemap(1));
        m_nel.push_back(conn.shape(0));
        m_el_offset.push_back(m_el_offset[index - 1] + m_nel[index - 1]);
    }

    /**
    \return Nodal coordinates.
    */
    xt::xtensor<double, 2> coor() const
    {
        return m_coor;
    }

    /**
    \return Connectivity.
    */
    xt::xtensor<size_t, 2> conn() const
    {
        return m_conn;
    }

    /**
    \param mesh_index Index of the mesh.
    \return Node-map for a given mesh.
    */
    xt::xtensor<size_t, 1> nodemap(size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index < m_map.size());
        return m_map[mesh_index];
    }

    /**
    \param mesh_index Index of the mesh.
    \return Element-map for a given mesh.
    */
    xt::xtensor<size_t, 1> elemmap(size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index < m_map.size());
        return xt::arange<size_t>(m_nel[mesh_index]) + m_el_offset[mesh_index];
    }

    /**
    Convert set of node-numbers for an original mesh to the stitched mesh.

    \param set Set of node-numbers.
    \param mesh_index Index of the mesh.
    \return Set of node-number for the stitched mesh.
    */
    xt::xtensor<size_t, 1> nodeset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index < m_map.size());
        GOOSEFEM_ASSERT(xt::amax(set)() < m_map[mesh_index].size());
        return detail::renum(set, m_map[mesh_index]);
    }

    /**
    Convert set of element-numbers for an original mesh to the stitched mesh.

    \param set Set of element-numbers.
    \param mesh_index Index of the mesh.
    \return Set of element-number for the stitched mesh.
    */
    xt::xtensor<size_t, 1> elemset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const
    {
        GOOSEFEM_ASSERT(mesh_index < m_map.size());
        GOOSEFEM_ASSERT(xt::amax(set)() < m_nel[mesh_index]);
        return set + m_el_offset[mesh_index];
    }

    /**
    Combine set of node-numbers for an original to the final mesh (removes duplicates).

    \param set List of node-sets per mesh.
    \return Combined node-set on the stitched mesh.
    */
    xt::xtensor<size_t, 1> nodeset(const std::vector<xt::xtensor<size_t, 1>>& set) const
    {
        GOOSEFEM_ASSERT(set.size() == m_map.size());

        size_t n = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            n += set[i].size();
        }

        xt::xtensor<size_t, 1> ret = xt::empty<size_t>({n});

        n = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            xt::view(ret, xt::range(n, n + set[i].size())) = this->nodeset(set[i], i);
            n += set[i].size();
        }

        return xt::unique(ret);
    }

    /**
    Combine set of element-numbers for an original to the final mesh.

    \param set List of element-sets per mesh.
    \return Combined element-set on the stitched mesh.
    */
    xt::xtensor<size_t, 1> elemset(const std::vector<xt::xtensor<size_t, 1>>& set) const
    {
        GOOSEFEM_ASSERT(set.size() == m_map.size());

        size_t n = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            n += set[i].size();
        }

        xt::xtensor<size_t, 1> ret = xt::empty<size_t>({n});

        n = 0;

        for (size_t i = 0; i < set.size(); ++i) {
            xt::view(ret, xt::range(n, n + set[i].size())) = this->elemset(set[i], i);
            n += set[i].size();
        }

        return ret;
    }
};

/**
\rst
Renumber indices to lowest possible index. For example:

.. math::

    \begin{bmatrix}
        0 & 1 \\
        5 & 4
    \end{bmatrix}

is renumbered to

.. math::

    \begin{bmatrix}
        0 & 1 \\
        3 & 2
    \end{bmatrix}

Or, in pseudo-code, the result of this function is that:

.. code-block:: python

    dofs = renumber(dofs)

    sort(unique(dofs[:])) == range(max(dofs+1))

.. tip::

    One can use the wrapper function :cpp:func:`GooseFEM::Mesh::renumber`.
    This class gives more advanced features.
\endrst
*/
class Renumber {
private:
    xt::xtensor<size_t, 1> m_renum;

public:
    Renumber() = default;

    /**
    \param dofs DOF-numbers.
    */
    template <class T>
    Renumber(const T& dofs)
    {
        size_t n = xt::amax(dofs)() + 1;
        size_t i = 0;

        xt::xtensor<size_t, 1> unique = xt::unique(dofs);

        m_renum = xt::empty<size_t>({n});

        for (auto& j : unique) {
            m_renum(j) = i;
            ++i;
        }
    }

    /**
    Get renumbered DOFs (same as ``Renumber::apply(dofs)``).

    \param dofs List of (DOF-)numbers.
    \return Renumbered list of (DOF-)numbers.
    */
    [[deprecated]]
    xt::xtensor<size_t, 2> get(const xt::xtensor<size_t, 2>& dofs) const
    {
        return this->apply(dofs);
    }

    /**
    Apply renumbering to other set.

    \param list List of (DOF-)numbers.
    \return Renumbered list of (DOF-)numbers.
    */
    template <class T>
    T apply(const T& list) const
    {
        return detail::renum(list, m_renum);
    }

    /**
    Get the list needed to renumber, e.g.:

        dofs_renumbered(i, j) = index(dofs(i, j))

    \return Renumber-index.
    */
    xt::xtensor<size_t, 1> index() const
    {
        return m_renum;
    }
};

/**
Renumber to lowest possible index (see GooseFEM::Mesh::Renumber).

\param dofs DOF-numbers.
\return Renumbered DOF-numbers.
*/
inline xt::xtensor<size_t, 2> renumber(const xt::xtensor<size_t, 2>& dofs)
{
    return Renumber(dofs).get(dofs);
}

/**
Reorder to lowest possible index, in specific order.

For example for ``Reorder({iiu, iip})`` after reordering:

    iiu = xt::range<size_t>(nnu);
    iip = xt::range<size_t>(nnp) + nnu;
*/
class Reorder {
private:
    xt::xtensor<size_t, 1> m_renum;

public:
    Reorder() = default;

    /**
    \param args List of (DOF-)numbers.
    */
    Reorder(const std::initializer_list<xt::xtensor<size_t, 1>> args)
    {
        size_t n = 0;
        size_t i = 0;

        for (auto& arg : args) {
            if (arg.size() == 0) {
                continue;
            }
            n = std::max(n, xt::amax(arg)() + 1);
        }

        #ifdef GOOSEFEM_ENABLE_ASSERT
        for (auto& arg : args) {
            GOOSEFEM_ASSERT(xt::unique(arg) == xt::sort(arg));
        }
        #endif

        m_renum = xt::empty<size_t>({n});

        for (auto& arg : args) {
            for (auto& j : arg) {
                m_renum(j) = i;
                ++i;
            }
        }
    }

    /**
    Get reordered DOFs (same as ``Reorder::apply(dofs)``).

    \param dofs List of (DOF-)numbers.
    \return Reordered list of (DOF-)numbers.
    */
    [[deprecated]]
    xt::xtensor<size_t, 2> get(const xt::xtensor<size_t, 2>& dofs) const
    {
        return this->apply(dofs);
    }

    /**
    Apply reordering to other set.

    \param list List of (DOF-)numbers.
    \return Reordered list of (DOF-)numbers.
    */
    template <class T>
    T apply(const T& list) const
    {
        T ret = T::from_shape(list.shape());

        auto jt = ret.begin();

        for (auto it = list.begin(); it != list.end(); ++it, ++jt) {
            *jt = m_renum(*it);
        }

        return ret;
    }

    /**
    Get the list needed to reorder, e.g.:

        dofs_reordered(i, j) = index(dofs(i, j))

    \return Reorder-index.
    */
    xt::xtensor<size_t, 1> index() const
    {
        return m_renum;
    }
};

/**
\rst
List with DOF-numbers in sequential order.
The output is a sequential list of DOF-numbers for each vector-component of each node.
For example for 3 nodes in 2 dimensions the output is

.. math::

    \begin{bmatrix}
        0 & 1 \\
        2 & 3 \\
        4 & 5
    \end{bmatrix}
\endrst

\param nnode Number of nodes.
\param ndim Number of dimensions.
\return DOF-numbers.
*/
inline xt::xtensor<size_t, 2> dofs(size_t nnode, size_t ndim)
{
    return xt::reshape_view(xt::arange<size_t>(nnode * ndim), {nnode, ndim});
}

/**
Number of elements connected to each node.

\param conn Connectivity.
\return Coordination per node.
*/
inline xt::xtensor<size_t, 1> coordination(const xt::xtensor<size_t, 2>& conn)
{
    size_t nnode = xt::amax(conn)() + 1;

    xt::xtensor<size_t, 1> N = xt::zeros<size_t>({nnode});

    for (auto it = conn.begin(); it != conn.end(); ++it) {
        N(*it) += 1;
    }

    return N;
}

/**
Elements connected to each node.

\param conn Connectivity.
\param sorted If ``true`` the output is sorted.
\return Elements per node.
*/
inline std::vector<std::vector<size_t>> elem2node(
    const xt::xtensor<size_t, 2>& conn,
    bool sorted = true)
{
    auto N = coordination(conn);
    auto nnode = N.size();

    std::vector<std::vector<size_t>> ret(nnode);
    for (size_t i = 0; i < nnode; ++i) {
        ret[i].reserve(N(i));
    }

    for (size_t e = 0; e < conn.shape(0); ++e) {
        for (size_t m = 0; m < conn.shape(1); ++m) {
            ret[conn(e, m)].push_back(e);
        }
    }

    if (sorted) {
        for (auto& row : ret) {
            std::sort(row.begin(), row.end());
        }
    }

    return ret;
}

/**
Return size of each element edge.

\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Edge-sizes per element.
*/
inline xt::xtensor<double, 2> edgesize(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn,
    ElementType type)
{
    GOOSEFEM_ASSERT(xt::amax(conn)() < coor.shape(0));

    if (type == ElementType::Quad4) {
        GOOSEFEM_ASSERT(coor.shape(1) == 2ul);
        GOOSEFEM_ASSERT(conn.shape(1) == 4ul);
        xt::xtensor<size_t, 1> n0 = xt::view(conn, xt::all(), 0);
        xt::xtensor<size_t, 1> n1 = xt::view(conn, xt::all(), 1);
        xt::xtensor<size_t, 1> n2 = xt::view(conn, xt::all(), 2);
        xt::xtensor<size_t, 1> n3 = xt::view(conn, xt::all(), 3);
        xt::xtensor<double, 1> x0 = xt::view(coor, xt::keep(n0), 0);
        xt::xtensor<double, 1> x1 = xt::view(coor, xt::keep(n1), 0);
        xt::xtensor<double, 1> x2 = xt::view(coor, xt::keep(n2), 0);
        xt::xtensor<double, 1> x3 = xt::view(coor, xt::keep(n3), 0);
        xt::xtensor<double, 1> y0 = xt::view(coor, xt::keep(n0), 1);
        xt::xtensor<double, 1> y1 = xt::view(coor, xt::keep(n1), 1);
        xt::xtensor<double, 1> y2 = xt::view(coor, xt::keep(n2), 1);
        xt::xtensor<double, 1> y3 = xt::view(coor, xt::keep(n3), 1);
        xt::xtensor<double, 2> ret = xt::empty<double>(conn.shape());
        xt::view(ret, xt::all(), 0) = xt::sqrt(xt::pow(x1 - x0, 2.0) + xt::pow(y1 - y0, 2.0));
        xt::view(ret, xt::all(), 1) = xt::sqrt(xt::pow(x2 - x1, 2.0) + xt::pow(y2 - y1, 2.0));
        xt::view(ret, xt::all(), 2) = xt::sqrt(xt::pow(x3 - x2, 2.0) + xt::pow(y3 - y2, 2.0));
        xt::view(ret, xt::all(), 3) = xt::sqrt(xt::pow(x0 - x3, 2.0) + xt::pow(y0 - y3, 2.0));
        return ret;
    }

    throw std::runtime_error("Element-type not implemented");
}

/**
Return size of each element edge.
The element-type is automatically determined, see GooseFEM::Mesh::defaultElementType.

\param coor Nodal coordinates.
\param conn Connectivity.
\return Edge-sizes per element.
*/
inline xt::xtensor<double, 2> edgesize(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn)
{
    return edgesize(coor, conn, defaultElementType(coor, conn));
}

/**
Coordinates of the center of each element.

\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Center of each element.
*/
inline xt::xtensor<double, 2> centers(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn,
    ElementType type)
{
    GOOSEFEM_ASSERT(xt::amax(conn)() < coor.shape(0));
    xt::xtensor<double, 2> ret = xt::zeros<double>({conn.shape(0), coor.shape(1)});

    if (type == ElementType::Quad4) {
        GOOSEFEM_ASSERT(coor.shape(1) == 2);
        GOOSEFEM_ASSERT(conn.shape(1) == 4);
        for (size_t i = 0; i < 4; ++i) {
            auto n = xt::view(conn, xt::all(), i);
            ret += xt::view(coor, xt::keep(n), xt::all());
        }
        ret /= 4.0;
        return ret;
    }

    throw std::runtime_error("Element-type not implemented");
}

/**
Coordinates of the center of each element.
The element-type is automatically determined, see GooseFEM::Mesh::defaultElementType.

\param coor Nodal coordinates.
\param conn Connectivity.
\return Center of each element.
*/
inline xt::xtensor<double, 2> centers(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn)
{
    return centers(coor, conn, defaultElementType(coor, conn));
}

/**
Convert an element-map to a node-map.

\param elem_map Element-map such that ``new_elvar = elvar[elem_map]``.
\param coor Nodal coordinates.
\param conn Connectivity.
\param type ElementType.
\return Node-map such that ``new_nodevar = nodevar[node_map]``
*/
inline xt::xtensor<size_t, 1> elemmap2nodemap(
    const xt::xtensor<size_t, 1>& elem_map,
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn,
    ElementType type)
{
    GOOSEFEM_ASSERT(xt::amax(conn)() < coor.shape(0));
    GOOSEFEM_ASSERT(elem_map.size() == conn.shape(0));
    size_t N = coor.shape(0);

    xt::xtensor<size_t, 1> ret = N * xt::ones<size_t>({N});

    if (type == ElementType::Quad4) {
        GOOSEFEM_ASSERT(coor.shape(1) == 2);
        GOOSEFEM_ASSERT(conn.shape(1) == 4);

        for (size_t i = 0; i < 4; ++i) {
            xt::xtensor<size_t, 1> t = N * xt::ones<size_t>({N});
            auto old_nd = xt::view(conn, xt::all(), i);
            auto new_nd = xt::view(conn, xt::keep(elem_map), i);
            xt::view(t, xt::keep(old_nd)) = new_nd;
            ret = xt::where(xt::equal(ret, N), t, ret);
        }

        return ret;
    }

    throw std::runtime_error("Element-type not implemented");
}

/**
Convert an element-map to a node-map.
The element-type is automatically determined, see GooseFEM::Mesh::defaultElementType.

\param elem_map Element-map such that ``new_elvar = elvar[elem_map]``.
\param coor Nodal coordinates.
\param conn Connectivity.
\return Node-map such that ``new_nodevar = nodevar[node_map]``
*/
inline xt::xtensor<size_t, 1> elemmap2nodemap(
    const xt::xtensor<size_t, 1>& elem_map,
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn)
{
    return elemmap2nodemap(elem_map, coor, conn, defaultElementType(coor, conn));
}

} // namespace Mesh
} // namespace GooseFEM

#endif

/**
Implementation of Mesh.h

\file Mesh.hpp
\copyright Copyright 2017. Tom de Geus. All rights reserved.
\license This project is released under the GNU Public License (GPLv3).
*/

#ifndef GOOSEFEM_MESH_HPP
#define GOOSEFEM_MESH_HPP

#include "Mesh.h"

namespace GooseFEM {
namespace Mesh {

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

inline xt::xtensor<size_t, 2> overlapping(
    const xt::xtensor<double, 2>& coor_a,
    const xt::xtensor<double, 2>& coor_b,
    double rtol,
    double atol)
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

inline ManualStitch::ManualStitch(
    const xt::xtensor<double, 2>& coor_a,
    const xt::xtensor<size_t, 2>& conn_a,
    const xt::xtensor<size_t, 1>& overlapping_nodes_a,
    const xt::xtensor<double, 2>& coor_b,
    const xt::xtensor<size_t, 2>& conn_b,
    const xt::xtensor<size_t, 1>& overlapping_nodes_b,
    bool check_position,
    double rtol,
    double atol)
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

inline xt::xtensor<double, 2> ManualStitch::coor() const
{
    return m_coor;
}

inline xt::xtensor<size_t, 2> ManualStitch::conn() const
{
    return m_conn;
}

inline size_t ManualStitch::nmesh() const
{
    return 2;
}

inline size_t ManualStitch::nelem() const
{
    return m_conn.shape(0);
}

inline size_t ManualStitch::nnode() const
{
    return m_coor.shape(0);
}

inline size_t ManualStitch::nne() const
{
    return m_conn.shape(1);
}

inline size_t ManualStitch::ndim() const
{
    return m_coor.shape(1);
}

inline xt::xtensor<size_t, 2> ManualStitch::dofs() const
{
    size_t nnode = this->nnode();
    size_t ndim = this->ndim();
    return xt::reshape_view(xt::arange<size_t>(nnode * ndim), {nnode, ndim});
}

inline xt::xtensor<size_t, 1> ManualStitch::nodemap(size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index <= 1);

    if (mesh_index == 0) {
        return xt::arange<size_t>(m_nnd_a);
    }

    return m_map_b;
}

inline xt::xtensor<size_t, 1> ManualStitch::elemmap(size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index <= 1);

    if (mesh_index == 0) {
        return xt::arange<size_t>(m_nel_a);
    }

    return xt::arange<size_t>(m_nel_b) + m_nel_a;
}

inline xt::xtensor<size_t, 1>
ManualStitch::nodeset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index <= 1);

    if (mesh_index == 0) {
        GOOSEFEM_ASSERT(xt::amax(set)() < m_nnd_a);
        return set;
    }

    GOOSEFEM_ASSERT(xt::amax(set)() < m_map_b.size());
    return detail::renum(set, m_map_b);
}

inline xt::xtensor<size_t, 1>
ManualStitch::elemset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index <= 1);

    if (mesh_index == 0) {
        GOOSEFEM_ASSERT(xt::amax(set)() < m_nel_a);
        return set;
    }

    GOOSEFEM_ASSERT(xt::amax(set)() < m_nel_b);
    return set + m_nel_a;
}

inline Stitch::Stitch(double rtol, double atol)
{
    m_rtol = rtol;
    m_atol = atol;
}

inline void Stitch::push_back(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn)
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

inline size_t Stitch::nmesh() const
{
    return m_map.size();
}

inline xt::xtensor<double, 2> Stitch::coor() const
{
    return m_coor;
}

inline xt::xtensor<size_t, 2> Stitch::conn() const
{
    return m_conn;
}

inline size_t Stitch::nelem() const
{
    return m_conn.shape(0);
}

inline size_t Stitch::nnode() const
{
    return m_coor.shape(0);
}

inline size_t Stitch::nne() const
{
    return m_conn.shape(1);
}

inline size_t Stitch::ndim() const
{
    return m_coor.shape(1);
}

inline xt::xtensor<size_t, 2> Stitch::dofs() const
{
    size_t nnode = this->nnode();
    size_t ndim = this->ndim();
    return xt::reshape_view(xt::arange<size_t>(nnode * ndim), {nnode, ndim});
}

inline xt::xtensor<size_t, 1> Stitch::nodemap(size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index < m_map.size());
    return m_map[mesh_index];
}

inline xt::xtensor<size_t, 1> Stitch::elemmap(size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index < m_map.size());
    return xt::arange<size_t>(m_nel[mesh_index]) + m_el_offset[mesh_index];
}

inline xt::xtensor<size_t, 1> Stitch::nodeset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index < m_map.size());
    GOOSEFEM_ASSERT(xt::amax(set)() < m_map[mesh_index].size());
    return detail::renum(set, m_map[mesh_index]);
}

inline xt::xtensor<size_t, 1> Stitch::elemset(const xt::xtensor<size_t, 1>& set, size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index < m_map.size());
    GOOSEFEM_ASSERT(xt::amax(set)() < m_nel[mesh_index]);
    return set + m_el_offset[mesh_index];
}

inline xt::xtensor<size_t, 1> Stitch::nodeset(const std::vector<xt::xtensor<size_t, 1>>& set) const
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

inline xt::xtensor<size_t, 1> Stitch::elemset(const std::vector<xt::xtensor<size_t, 1>>& set) const
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

template <class T>
inline Renumber::Renumber(const T& dofs)
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

inline xt::xtensor<size_t, 2> Renumber::get(const xt::xtensor<size_t, 2>& dofs) const
{
    GOOSEFEM_WARNING("Renumber::get is deprecated, use Renumber::apply");
    return this->apply(dofs);
}

template <class T>
inline T Renumber::apply(const T& list) const
{
    return detail::renum(list, m_renum);
}

inline xt::xtensor<size_t, 1> Renumber::index() const
{
    return m_renum;
}

inline xt::xtensor<size_t, 2> renumber(const xt::xtensor<size_t, 2>& dofs)
{
    return Renumber(dofs).apply(dofs);
}

inline Reorder::Reorder(const std::initializer_list<xt::xtensor<size_t, 1>> args)
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

inline xt::xtensor<size_t, 2> Reorder::get(const xt::xtensor<size_t, 2>& dofs) const
{
    GOOSEFEM_WARNING("Reorder::get is deprecated, use Reorder::apply");
    return this->apply(dofs);
}

template <class T>
inline T Reorder::apply(const T& list) const
{
    T ret = T::from_shape(list.shape());

    auto jt = ret.begin();

    for (auto it = list.begin(); it != list.end(); ++it, ++jt) {
        *jt = m_renum(*it);
    }

    return ret;
}

inline xt::xtensor<size_t, 1> Reorder::index() const
{
    return m_renum;
}

inline xt::xtensor<size_t, 2> dofs(size_t nnode, size_t ndim)
{
    return xt::reshape_view(xt::arange<size_t>(nnode * ndim), {nnode, ndim});
}

inline xt::xtensor<size_t, 1> coordination(const xt::xtensor<size_t, 2>& conn)
{
    size_t nnode = xt::amax(conn)() + 1;

    xt::xtensor<size_t, 1> N = xt::zeros<size_t>({nnode});

    for (auto it = conn.begin(); it != conn.end(); ++it) {
        N(*it) += 1;
    }

    return N;
}

inline std::vector<std::vector<size_t>> elem2node(const xt::xtensor<size_t, 2>& conn, bool sorted)
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

inline xt::xtensor<double, 2> edgesize(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn)
{
    return edgesize(coor, conn, defaultElementType(coor, conn));
}

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

inline xt::xtensor<double, 2> centers(
    const xt::xtensor<double, 2>& coor,
    const xt::xtensor<size_t, 2>& conn)
{
    return centers(coor, conn, defaultElementType(coor, conn));
}

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

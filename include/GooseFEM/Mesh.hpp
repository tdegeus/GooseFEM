/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

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

inline Stitch::Stitch(
        const xt::xtensor<double, 2>& coor_a,
        const xt::xtensor<size_t, 2>& conn_a,
        const xt::xtensor<size_t, 1>& overlapping_nodes_a,
        const xt::xtensor<double, 2>& coor_b,
        const xt::xtensor<size_t, 2>& conn_b,
        const xt::xtensor<size_t, 1>& overlapping_nodes_b,
        bool check_position)
{
    GOOSEFEM_ASSERT(xt::has_shape(overlapping_nodes_a, overlapping_nodes_b.shape()));
    GOOSEFEM_ASSERT(coor_a.shape(1) == coor_b.shape(1));
    GOOSEFEM_ASSERT(conn_a.shape(1) == conn_b.shape(1));

    if (check_position) {
        GOOSEFEM_ASSERT(xt::allclose(
            xt::view(coor_a, xt::keep(overlapping_nodes_a), xt::all()),
            xt::view(coor_b, xt::keep(overlapping_nodes_b), xt::all())));
    }

    size_t nnda = coor_a.shape(0);
    size_t nndb = coor_b.shape(0);
    size_t ndim = coor_a.shape(1);
    size_t nelim = overlapping_nodes_a.size();

    size_t nela = conn_a.shape(0);
    size_t nelb = conn_b.shape(0);
    size_t nne = conn_a.shape(1);

    m_nel_a = nela;

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

inline xt::xtensor<double, 2> Stitch::coor() const
{
    return m_coor;
}

inline xt::xtensor<size_t, 2> Stitch::conn() const
{
    return m_conn;
}

inline xt::xtensor<size_t, 1> Stitch::nodeset(const xt::xtensor<size_t, 1>& set, size_t mesh) const
{
    GOOSEFEM_ASSERT(mesh <= 1);

    if (mesh == 0) {
        return set;
    }

    return detail::renum(set, m_map_b);
}

inline xt::xtensor<size_t, 1> Stitch::elementset(const xt::xtensor<size_t, 1>& set, size_t mesh) const
{
    GOOSEFEM_ASSERT(mesh <= 1);

    if (mesh == 0) {
        return set;
    }

    return set + m_nel_a;
}

inline Renumber::Renumber(const xt::xarray<size_t>& dofs)
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

// ret(i,j) = renum(list(i,j))
template <class T>
T Renumber::apply(const T& list) const
{
    return detail::renum(list, m_renum);
}

inline xt::xtensor<size_t, 2> Renumber::get(const xt::xtensor<size_t, 2>& dofs) const
{
    return this->apply(dofs);
}

inline xt::xtensor<size_t, 1> Renumber::index() const
{
    return m_renum;
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
    return this->apply(dofs);
}

inline xt::xtensor<size_t, 1> Reorder::index() const
{
    return m_renum;
}

// apply renumbering, e.g. for a matrix:
//
//   ret(i,j) = renum(list(i,j))

template <class T>
T Reorder::apply(const T& list) const
{
    T ret = T::from_shape(list.shape());

    auto jt = ret.begin();

    for (auto it = list.begin(); it != list.end(); ++it, ++jt) {
        *jt = m_renum(*it);
    }

    return ret;
}

inline xt::xtensor<size_t, 2> renumber(const xt::xtensor<size_t, 2>& dofs)
{
    return Renumber(dofs).get(dofs);
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
    const xt::xtensor<double, 2>& coor, const xt::xtensor<size_t, 2>& conn, ElementType type)
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

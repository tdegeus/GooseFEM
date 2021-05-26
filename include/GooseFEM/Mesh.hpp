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

template <class S, class T>
inline ElementType defaultElementType(const S& coor, const T& conn)
{
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);

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

template <class D>
inline auto RegularBase<D>::nelem() const
{
    return derived_cast().m_nelem;
}

template <class D>
inline auto RegularBase<D>::nnode() const
{
    return derived_cast().m_nnode;
}

template <class D>
inline auto RegularBase<D>::nne() const
{
    return derived_cast().m_nne;
}

template <class D>
inline auto RegularBase<D>::ndim() const
{
    return derived_cast().m_ndim;
}

template <class D>
inline auto RegularBase<D>::nelx() const
{
    return derived_cast().nelx_impl();
}

template <class D>
inline auto RegularBase<D>::nely() const
{
    return derived_cast().nely_impl();
}

template <class D>
inline auto RegularBase<D>::h() const
{
    return derived_cast().m_h;
}

template <class D>
inline auto RegularBase<D>::getElementType() const
{
    return derived_cast().getElementType_impl();
}

template <class D>
inline auto RegularBase<D>::coor() const
{
    return derived_cast().coor_impl();
}

template <class D>
inline auto RegularBase<D>::conn() const
{
    return derived_cast().conn_impl();
}

template <class D>
inline auto RegularBase<D>::dofs() const
{
    return GooseFEM::Mesh::dofs(this->nnode(), this->ndim());
}

template <class D>
inline auto RegularBase<D>::dofsPeriodic() const
{
    xt::xtensor<size_t, 2> ret = this->dofs();
    xt::xtensor<size_t, 2> nodePer = this->nodesPeriodic();
    xt::xtensor<size_t, 1> independent = xt::view(nodePer, xt::all(), 0);
    xt::xtensor<size_t, 1> dependent = xt::view(nodePer, xt::all(), 1);

    for (size_t j = 0; j < this->ndim(); ++j) {
        xt::view(ret, xt::keep(dependent), j) = xt::view(ret, xt::keep(independent), j);
    }

    return GooseFEM::Mesh::renumber(ret);
}

template <class D>
inline auto RegularBase<D>::derived_cast() -> derived_type&
{
    return *static_cast<derived_type*>(this);
}

template <class D>
inline auto RegularBase<D>::derived_cast() const -> const derived_type&
{
    return *static_cast<const derived_type*>(this);
}

template <class D>
inline auto RegularBase<D>::nodesPeriodic() const
{
    return derived_cast().nodesPeriodic_impl();
}

template <class D>
inline auto RegularBase<D>::nodesOrigin() const
{
    return derived_cast().nodesOrigin_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesBottomEdge() const
{
    return derived_cast().nodesBottomEdge_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesTopEdge() const
{
    return derived_cast().nodesTopEdge_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesLeftEdge() const
{
    return derived_cast().nodesLeftEdge_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesRightEdge() const
{
    return derived_cast().nodesRightEdge_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesBottomOpenEdge() const
{
    return derived_cast().nodesBottomOpenEdge_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesTopOpenEdge() const
{
    return derived_cast().nodesTopOpenEdge_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesLeftOpenEdge() const
{
    return derived_cast().nodesLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesRightOpenEdge() const
{
    return derived_cast().nodesRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesBottomLeftCorner() const
{
    return derived_cast().nodesBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesTopLeftCorner() const
{
    return derived_cast().nodesTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesBottomRightCorner() const
{
    return derived_cast().nodesBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesTopRightCorner() const
{
    return derived_cast().nodesTopRightCorner_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesLeftBottomCorner() const
{
    return derived_cast().nodesBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesLeftTopCorner() const
{
    return derived_cast().nodesTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesRightBottomCorner() const
{
    return derived_cast().nodesBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase2d<D>::nodesRightTopCorner() const
{
    return derived_cast().nodesTopRightCorner_impl();
}

template <class D>
inline xt::xtensor<size_t, 2> RegularBase2d<D>::nodesPeriodic_impl() const
{
    xt::xtensor<size_t, 1> bot = derived_cast().nodesBottomOpenEdge_impl();
    xt::xtensor<size_t, 1> top = derived_cast().nodesTopOpenEdge_impl();
    xt::xtensor<size_t, 1> lft = derived_cast().nodesLeftOpenEdge_impl();
    xt::xtensor<size_t, 1> rgt = derived_cast().nodesRightOpenEdge_impl();
    std::array<size_t, 2> shape = {bot.size() + lft.size() + size_t(3), size_t(2)};
    auto ret = xt::xtensor<size_t, 2>::from_shape(shape);

    ret(0, 0) = derived_cast().nodesBottomLeftCorner_impl();
    ret(0, 1) = derived_cast().nodesBottomRightCorner_impl();

    ret(1, 0) = derived_cast().nodesBottomLeftCorner_impl();
    ret(1, 1) = derived_cast().nodesTopRightCorner_impl();

    ret(2, 0) = derived_cast().nodesBottomLeftCorner_impl();
    ret(2, 1) = derived_cast().nodesTopLeftCorner_impl();

    size_t i = 3;

    xt::view(ret, xt::range(i, i + bot.size()), 0) = bot;
    xt::view(ret, xt::range(i, i + bot.size()), 1) = top;

    i += bot.size();

    xt::view(ret, xt::range(i, i + lft.size()), 0) = lft;
    xt::view(ret, xt::range(i, i + lft.size()), 1) = rgt;

    return ret;
}

template <class D>
inline auto RegularBase2d<D>::nodesOrigin_impl() const
{
    return derived_cast().nodesBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase2d<D>::derived_cast() -> derived_type&
{
    return *static_cast<derived_type*>(this);
}

template <class D>
inline auto RegularBase2d<D>::derived_cast() const -> const derived_type&
{
    return *static_cast<const derived_type*>(this);
}

template <class D>
inline auto RegularBase3d<D>::nelz() const
{
    return derived_cast().nelz_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottom() const
{
    return derived_cast().nodesBottom_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTop() const
{
    return derived_cast().nodesTop_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeft() const
{
    return derived_cast().nodesLeft_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRight() const
{
    return derived_cast().nodesRight_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFront() const
{
    return derived_cast().nodesFront_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBack() const
{
    return derived_cast().nodesBack_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontBottomEdge() const
{
    return derived_cast().nodesFrontBottomEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontTopEdge() const
{
    return derived_cast().nodesFrontTopEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontLeftEdge() const
{
    return derived_cast().nodesFrontLeftEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontRightEdge() const
{
    return derived_cast().nodesFrontRightEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackBottomEdge() const
{
    return derived_cast().nodesBackBottomEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackTopEdge() const
{
    return derived_cast().nodesBackTopEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackLeftEdge() const
{
    return derived_cast().nodesBackLeftEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackRightEdge() const
{
    return derived_cast().nodesBackRightEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomLeftEdge() const
{
    return derived_cast().nodesBottomLeftEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomRightEdge() const
{
    return derived_cast().nodesBottomRightEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopLeftEdge() const
{
    return derived_cast().nodesTopLeftEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopRightEdge() const
{
    return derived_cast().nodesTopRightEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomFrontEdge() const
{
    return derived_cast().nodesFrontBottomEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomBackEdge() const
{
    return derived_cast().nodesBackBottomEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopFrontEdge() const
{
    return derived_cast().nodesFrontTopEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopBackEdge() const
{
    return derived_cast().nodesBackTopEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftBottomEdge() const
{
    return derived_cast().nodesBottomLeftEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftFrontEdge() const
{
    return derived_cast().nodesFrontLeftEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftBackEdge() const
{
    return derived_cast().nodesBackLeftEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftTopEdge() const
{
    return derived_cast().nodesTopLeftEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightBottomEdge() const
{
    return derived_cast().nodesBottomRightEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightTopEdge() const
{
    return derived_cast().nodesTopRightEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightFrontEdge() const
{
    return derived_cast().nodesFrontRightEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightBackEdge() const
{
    return derived_cast().nodesBackRightEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontFace() const
{
    return derived_cast().nodesFrontFace_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackFace() const
{
    return derived_cast().nodesBackFace_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftFace() const
{
    return derived_cast().nodesLeftFace_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightFace() const
{
    return derived_cast().nodesRightFace_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomFace() const
{
    return derived_cast().nodesBottomFace_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopFace() const
{
    return derived_cast().nodesTopFace_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontBottomOpenEdge() const
{
    return derived_cast().nodesFrontBottomOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontTopOpenEdge() const
{
    return derived_cast().nodesFrontTopOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontLeftOpenEdge() const
{
    return derived_cast().nodesFrontLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontRightOpenEdge() const
{
    return derived_cast().nodesFrontRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackBottomOpenEdge() const
{
    return derived_cast().nodesBackBottomOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackTopOpenEdge() const
{
    return derived_cast().nodesBackTopOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackLeftOpenEdge() const
{
    return derived_cast().nodesBackLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackRightOpenEdge() const
{
    return derived_cast().nodesBackRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomLeftOpenEdge() const
{
    return derived_cast().nodesBottomLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomRightOpenEdge() const
{
    return derived_cast().nodesBottomRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopLeftOpenEdge() const
{
    return derived_cast().nodesTopLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopRightOpenEdge() const
{
    return derived_cast().nodesTopRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomFrontOpenEdge() const
{
    return derived_cast().nodesFrontBottomOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomBackOpenEdge() const
{
    return derived_cast().nodesBackBottomOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopFrontOpenEdge() const
{
    return derived_cast().nodesFrontTopOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopBackOpenEdge() const
{
    return derived_cast().nodesBackTopOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftBottomOpenEdge() const
{
    return derived_cast().nodesBottomLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftFrontOpenEdge() const
{
    return derived_cast().nodesFrontLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftBackOpenEdge() const
{
    return derived_cast().nodesBackLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftTopOpenEdge() const
{
    return derived_cast().nodesTopLeftOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightBottomOpenEdge() const
{
    return derived_cast().nodesBottomRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightTopOpenEdge() const
{
    return derived_cast().nodesTopRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightFrontOpenEdge() const
{
    return derived_cast().nodesFrontRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightBackOpenEdge() const
{
    return derived_cast().nodesBackRightOpenEdge_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontBottomLeftCorner() const
{
    return derived_cast().nodesFrontBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontBottomRightCorner() const
{
    return derived_cast().nodesFrontBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontTopLeftCorner() const
{
    return derived_cast().nodesFrontTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontTopRightCorner() const
{
    return derived_cast().nodesFrontTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackBottomLeftCorner() const
{
    return derived_cast().nodesBackBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackBottomRightCorner() const
{
    return derived_cast().nodesBackBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackTopLeftCorner() const
{
    return derived_cast().nodesBackTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackTopRightCorner() const
{
    return derived_cast().nodesBackTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontLeftBottomCorner() const
{
    return derived_cast().nodesFrontBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomFrontLeftCorner() const
{
    return derived_cast().nodesFrontBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomLeftFrontCorner() const
{
    return derived_cast().nodesFrontBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftFrontBottomCorner() const
{
    return derived_cast().nodesFrontBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftBottomFrontCorner() const
{
    return derived_cast().nodesFrontBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontRightBottomCorner() const
{
    return derived_cast().nodesFrontBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomFrontRightCorner() const
{
    return derived_cast().nodesFrontBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomRightFrontCorner() const
{
    return derived_cast().nodesFrontBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightFrontBottomCorner() const
{
    return derived_cast().nodesFrontBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightBottomFrontCorner() const
{
    return derived_cast().nodesFrontBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontLeftTopCorner() const
{
    return derived_cast().nodesFrontTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopFrontLeftCorner() const
{
    return derived_cast().nodesFrontTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopLeftFrontCorner() const
{
    return derived_cast().nodesFrontTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftFrontTopCorner() const
{
    return derived_cast().nodesFrontTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftTopFrontCorner() const
{
    return derived_cast().nodesFrontTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesFrontRightTopCorner() const
{
    return derived_cast().nodesFrontTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopFrontRightCorner() const
{
    return derived_cast().nodesFrontTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopRightFrontCorner() const
{
    return derived_cast().nodesFrontTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightFrontTopCorner() const
{
    return derived_cast().nodesFrontTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightTopFrontCorner() const
{
    return derived_cast().nodesFrontTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackLeftBottomCorner() const
{
    return derived_cast().nodesBackBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomBackLeftCorner() const
{
    return derived_cast().nodesBackBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomLeftBackCorner() const
{
    return derived_cast().nodesBackBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftBackBottomCorner() const
{
    return derived_cast().nodesBackBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftBottomBackCorner() const
{
    return derived_cast().nodesBackBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackRightBottomCorner() const
{
    return derived_cast().nodesBackBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomBackRightCorner() const
{
    return derived_cast().nodesBackBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBottomRightBackCorner() const
{
    return derived_cast().nodesBackBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightBackBottomCorner() const
{
    return derived_cast().nodesBackBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightBottomBackCorner() const
{
    return derived_cast().nodesBackBottomRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackLeftTopCorner() const
{
    return derived_cast().nodesBackTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopBackLeftCorner() const
{
    return derived_cast().nodesBackTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopLeftBackCorner() const
{
    return derived_cast().nodesBackTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftBackTopCorner() const
{
    return derived_cast().nodesBackTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesLeftTopBackCorner() const
{
    return derived_cast().nodesBackTopLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesBackRightTopCorner() const
{
    return derived_cast().nodesBackTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopBackRightCorner() const
{
    return derived_cast().nodesBackTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesTopRightBackCorner() const
{
    return derived_cast().nodesBackTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightBackTopCorner() const
{
    return derived_cast().nodesBackTopRightCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::nodesRightTopBackCorner() const
{
    return derived_cast().nodesBackTopRightCorner_impl();
}

template <class D>
inline xt::xtensor<size_t, 2> RegularBase3d<D>::nodesPeriodic_impl() const
{
    xt::xtensor<size_t, 1> fro = derived_cast().nodesFrontFace_impl();
    xt::xtensor<size_t, 1> bck = derived_cast().nodesBackFace_impl();
    xt::xtensor<size_t, 1> lft = derived_cast().nodesLeftFace_impl();
    xt::xtensor<size_t, 1> rgt = derived_cast().nodesRightFace_impl();
    xt::xtensor<size_t, 1> bot = derived_cast().nodesBottomFace_impl();
    xt::xtensor<size_t, 1> top = derived_cast().nodesTopFace_impl();

    xt::xtensor<size_t, 1> froBot = derived_cast().nodesFrontBottomOpenEdge_impl();
    xt::xtensor<size_t, 1> froTop = derived_cast().nodesFrontTopOpenEdge_impl();
    xt::xtensor<size_t, 1> froLft = derived_cast().nodesFrontLeftOpenEdge_impl();
    xt::xtensor<size_t, 1> froRgt = derived_cast().nodesFrontRightOpenEdge_impl();
    xt::xtensor<size_t, 1> bckBot = derived_cast().nodesBackBottomOpenEdge_impl();
    xt::xtensor<size_t, 1> bckTop = derived_cast().nodesBackTopOpenEdge_impl();
    xt::xtensor<size_t, 1> bckLft = derived_cast().nodesBackLeftOpenEdge_impl();
    xt::xtensor<size_t, 1> bckRgt = derived_cast().nodesBackRightOpenEdge_impl();
    xt::xtensor<size_t, 1> botLft = derived_cast().nodesBottomLeftOpenEdge_impl();
    xt::xtensor<size_t, 1> botRgt = derived_cast().nodesBottomRightOpenEdge_impl();
    xt::xtensor<size_t, 1> topLft = derived_cast().nodesTopLeftOpenEdge_impl();
    xt::xtensor<size_t, 1> topRgt = derived_cast().nodesTopRightOpenEdge_impl();

    size_t tface = fro.size() + lft.size() + bot.size();
    size_t tedge = 3 * froBot.size() + 3 * froLft.size() + 3 * botLft.size();
    size_t tnode = 7;
    xt::xtensor<size_t, 2> ret = xt::empty<size_t>({tface + tedge + tnode, std::size_t(2)});

    size_t i = 0;

    ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
    ret(i, 1) = derived_cast().nodesFrontBottomRightCorner_impl();
    ++i;
    ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
    ret(i, 1) = derived_cast().nodesBackBottomRightCorner_impl();
    ++i;
    ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
    ret(i, 1) = derived_cast().nodesBackBottomLeftCorner_impl();
    ++i;
    ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
    ret(i, 1) = derived_cast().nodesFrontTopLeftCorner_impl();
    ++i;
    ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
    ret(i, 1) = derived_cast().nodesFrontTopRightCorner_impl();
    ++i;
    ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
    ret(i, 1) = derived_cast().nodesBackTopRightCorner_impl();
    ++i;
    ret(i, 0) = derived_cast().nodesFrontBottomLeftCorner_impl();
    ret(i, 1) = derived_cast().nodesBackTopLeftCorner_impl();
    ++i;

    for (size_t j = 0; j < froBot.size(); ++j) {
        ret(i, 0) = froBot(j);
        ret(i, 1) = bckBot(j);
        ++i;
    }
    for (size_t j = 0; j < froBot.size(); ++j) {
        ret(i, 0) = froBot(j);
        ret(i, 1) = bckTop(j);
        ++i;
    }
    for (size_t j = 0; j < froBot.size(); ++j) {
        ret(i, 0) = froBot(j);
        ret(i, 1) = froTop(j);
        ++i;
    }
    for (size_t j = 0; j < botLft.size(); ++j) {
        ret(i, 0) = botLft(j);
        ret(i, 1) = botRgt(j);
        ++i;
    }
    for (size_t j = 0; j < botLft.size(); ++j) {
        ret(i, 0) = botLft(j);
        ret(i, 1) = topRgt(j);
        ++i;
    }
    for (size_t j = 0; j < botLft.size(); ++j) {
        ret(i, 0) = botLft(j);
        ret(i, 1) = topLft(j);
        ++i;
    }
    for (size_t j = 0; j < froLft.size(); ++j) {
        ret(i, 0) = froLft(j);
        ret(i, 1) = froRgt(j);
        ++i;
    }
    for (size_t j = 0; j < froLft.size(); ++j) {
        ret(i, 0) = froLft(j);
        ret(i, 1) = bckRgt(j);
        ++i;
    }
    for (size_t j = 0; j < froLft.size(); ++j) {
        ret(i, 0) = froLft(j);
        ret(i, 1) = bckLft(j);
        ++i;
    }

    for (size_t j = 0; j < fro.size(); ++j) {
        ret(i, 0) = fro(j);
        ret(i, 1) = bck(j);
        ++i;
    }
    for (size_t j = 0; j < lft.size(); ++j) {
        ret(i, 0) = lft(j);
        ret(i, 1) = rgt(j);
        ++i;
    }
    for (size_t j = 0; j < bot.size(); ++j) {
        ret(i, 0) = bot(j);
        ret(i, 1) = top(j);
        ++i;
    }

    return ret;
}

template <class D>
inline auto RegularBase3d<D>::nodesOrigin_impl() const
{
    return derived_cast().nodesFrontBottomLeftCorner_impl();
}

template <class D>
inline auto RegularBase3d<D>::derived_cast() -> derived_type&
{
    return *static_cast<derived_type*>(this);
}

template <class D>
inline auto RegularBase3d<D>::derived_cast() const -> const derived_type&
{
    return *static_cast<const derived_type*>(this);
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

template <class S, class T>
inline xt::xtensor<size_t, 2> overlapping(
    const S& coor_a,
    const T& coor_b,
    double rtol,
    double atol)
{
    GOOSEFEM_ASSERT(coor_a.dimension() == 2);
    GOOSEFEM_ASSERT(coor_b.dimension() == 2);
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

template <class CA, class EA, class NA, class CB, class EB, class NB>
inline ManualStitch::ManualStitch(
    const CA& coor_a,
    const EA& conn_a,
    const NA& overlapping_nodes_a,
    const CB& coor_b,
    const EB& conn_b,
    const NB& overlapping_nodes_b,
    bool check_position,
    double rtol,
    double atol)
{
    UNUSED(rtol);
    UNUSED(atol);

    GOOSEFEM_ASSERT(coor_a.dimension() == 2);
    GOOSEFEM_ASSERT(conn_a.dimension() == 2);
    GOOSEFEM_ASSERT(overlapping_nodes_a.dimension() == 1);
    GOOSEFEM_ASSERT(coor_b.dimension() == 2);
    GOOSEFEM_ASSERT(conn_b.dimension() == 2);
    GOOSEFEM_ASSERT(overlapping_nodes_b.dimension() == 1);
    GOOSEFEM_ASSERT(xt::has_shape(overlapping_nodes_a, overlapping_nodes_b.shape()));
    GOOSEFEM_ASSERT(coor_a.shape(1) == coor_b.shape(1));
    GOOSEFEM_ASSERT(conn_a.shape(1) == conn_b.shape(1));

    if (check_position) {
        GOOSEFEM_CHECK(xt::allclose(
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

inline std::vector<xt::xtensor<size_t, 1>> ManualStitch::nodemap() const
{
    std::vector<xt::xtensor<size_t, 1>> ret(this->nmesh());
    for (size_t i = 0; i < this->nmesh(); ++i) {
        ret[i] = this->nodemap(i);
    }
    return ret;
}

inline std::vector<xt::xtensor<size_t, 1>> ManualStitch::elemmap() const
{
    std::vector<xt::xtensor<size_t, 1>> ret(this->nmesh());
    for (size_t i = 0; i < this->nmesh(); ++i) {
        ret[i] = this->elemmap(i);
    }
    return ret;
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

template <class T>
inline T ManualStitch::nodeset(const T& set, size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index <= 1);

    if (mesh_index == 0) {
        GOOSEFEM_ASSERT(xt::amax(set)() < m_nnd_a);
        return set;
    }

    GOOSEFEM_ASSERT(xt::amax(set)() < m_map_b.size());
    return detail::renum(set, m_map_b);
}

template <class T>
inline T ManualStitch::elemset(const T& set, size_t mesh_index) const
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

template <class C, class E>
inline void Stitch::push_back(const C& coor, const E& conn)
{
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);

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

    ManualStitch stitch(
        m_coor, m_conn, xt::eval(xt::view(overlap, 0, xt::all())),
        coor, conn, xt::eval(xt::view(overlap, 1, xt::all())),
        false);

    m_coor = stitch.coor();
    m_conn = stitch.conn();
    m_map.push_back(stitch.nodemap(1));
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

inline std::vector<xt::xtensor<size_t, 1>> Stitch::nodemap() const
{
    std::vector<xt::xtensor<size_t, 1>> ret(this->nmesh());
    for (size_t i = 0; i < this->nmesh(); ++i) {
        ret[i] = this->nodemap(i);
    }
    return ret;
}

inline std::vector<xt::xtensor<size_t, 1>> Stitch::elemmap() const
{
    std::vector<xt::xtensor<size_t, 1>> ret(this->nmesh());
    for (size_t i = 0; i < this->nmesh(); ++i) {
        ret[i] = this->elemmap(i);
    }
    return ret;
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

template <class T>
inline T Stitch::nodeset(const T& set, size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index < m_map.size());
    GOOSEFEM_ASSERT(xt::amax(set)() < m_map[mesh_index].size());
    return detail::renum(set, m_map[mesh_index]);
}

template <class T>
inline T Stitch::elemset(const T& set, size_t mesh_index) const
{
    GOOSEFEM_ASSERT(mesh_index < m_map.size());
    GOOSEFEM_ASSERT(xt::amax(set)() < m_nel[mesh_index]);
    return set + m_el_offset[mesh_index];
}

template <class T>
inline T Stitch::nodeset(const std::vector<T>& set) const
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

template <class T>
inline T Stitch::elemset(const std::vector<T>& set) const
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
inline T Stitch::nodeset(std::initializer_list<T> set) const
{
    return this->nodeset(std::vector<T>(set));
}

template <class T>
inline T Stitch::elemset(std::initializer_list<T> set) const
{
    return this->elemset(std::vector<T>(set));
}

inline Vstack::Vstack(bool check_overlap, double rtol, double atol)
{
    m_check_overlap = check_overlap;
    m_rtol = rtol;
    m_atol = atol;
}

template <class C, class E, class N>
inline void Vstack::push_back(const C& coor, const E& conn, const N& nodes_bot, const N& nodes_top)
{
    if (m_map.size() == 0) {
        m_coor = coor;
        m_conn = conn;
        m_map.push_back(xt::eval(xt::arange<size_t>(coor.shape(0))));
        m_nel.push_back(conn.shape(0));
        m_el_offset.push_back(0);
        m_nodes_bot.push_back(nodes_bot);
        m_nodes_top.push_back(nodes_top);
        return;
    }

    GOOSEFEM_ASSERT(nodes_bot.size() == m_nodes_top.back().size());

    size_t index = m_map.size();

    double shift = xt::amax(xt::view(m_coor, xt::all(), 1))();
    auto x = coor;
    xt::view(x, xt::all(), 1) += shift;

    ManualStitch stitch(
        m_coor, m_conn, m_nodes_top.back(),
        x, conn, nodes_bot,
        m_check_overlap, m_rtol, m_atol);

    m_nodes_bot.push_back(stitch.nodeset(nodes_bot, 1));
    m_nodes_top.push_back(stitch.nodeset(nodes_top, 1));

    m_coor = stitch.coor();
    m_conn = stitch.conn();
    m_map.push_back(stitch.nodemap(1));
    m_nel.push_back(conn.shape(0));
    m_el_offset.push_back(m_el_offset[index - 1] + m_nel[index - 1]);
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

template <class T>
inline T Renumber::apply(const T& list) const
{
    return detail::renum(list, m_renum);
}

inline xt::xtensor<size_t, 1> Renumber::index() const
{
    return m_renum;
}

template <class T>
inline T renumber(const T& dofs)
{
    return Renumber(dofs).apply(dofs);
}

template <class T>
inline Reorder::Reorder(const std::initializer_list<T> args)
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

template <class E>
inline xt::xtensor<size_t, 1> coordination(const E& conn)
{
    GOOSEFEM_ASSERT(conn.dimension() == 2);

    size_t nnode = xt::amax(conn)() + 1;

    xt::xtensor<size_t, 1> N = xt::zeros<size_t>({nnode});

    for (auto it = conn.begin(); it != conn.end(); ++it) {
        N(*it) += 1;
    }

    return N;
}

template <class E>
inline std::vector<std::vector<size_t>> elem2node(const E& conn, bool sorted)
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

template <class C, class E>
inline xt::xtensor<double, 2> edgesize(const C& coor, const E& conn, ElementType type)
{
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);
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

template <class C, class E>
inline xt::xtensor<double, 2> edgesize(const C& coor, const E& conn)
{
    return edgesize(coor, conn, defaultElementType(coor, conn));
}

template <class C, class E>
inline xt::xtensor<double, 2> centers(const C& coor, const E& conn, ElementType type)
{
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);
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

template <class C, class E>
inline xt::xtensor<double, 2> centers(const C& coor, const E& conn)
{
    return centers(coor, conn, defaultElementType(coor, conn));
}

template <class T, class C, class E>
inline xt::xtensor<size_t, 1> elemmap2nodemap(
    const T& elem_map,
    const C& coor,
    const E& conn,
    ElementType type)
{
    GOOSEFEM_ASSERT(elem_map.dimension() == 1);
    GOOSEFEM_ASSERT(coor.dimension() == 2);
    GOOSEFEM_ASSERT(conn.dimension() == 2);
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

template <class T, class C, class E>
inline xt::xtensor<size_t, 1> elemmap2nodemap(const T& elem_map,const C& coor, const E& conn)
{
    return elemmap2nodemap(elem_map, coor, conn, defaultElementType(coor, conn));
}

} // namespace Mesh
} // namespace GooseFEM

#endif

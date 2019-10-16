/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHHEX8_HPP
#define GOOSEFEM_MESHHEX8_HPP

// -------------------------------------------------------------------------------------------------

#include "MeshHex8.h"

// =================================================================================================

namespace GooseFEM {
namespace Mesh {
namespace Hex8 {

// -------------------------------------------------------------------------------------------------

inline Regular::Regular(size_t nelx, size_t nely, size_t nelz, double h):
m_h(h), m_nelx(nelx), m_nely(nely), m_nelz(nelz)
{
  GOOSEFEM_ASSERT(m_nelx >= 1ul);
  GOOSEFEM_ASSERT(m_nely >= 1ul);
  GOOSEFEM_ASSERT(m_nelz >= 1ul);

  m_nnode = (m_nelx+1) * (m_nely+1) * (m_nelz+1);
  m_nelem =  m_nelx    *  m_nely    *  m_nelz   ;
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nelem() const
{
  return m_nelem;
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nnode() const
{
  return m_nnode;
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nne() const
{
  return m_nne;
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

inline ElementType Regular::getElementType() const
{
  return ElementType::Hex8;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Regular::coor() const
{
  xt::xtensor<double,2> out = xt::empty<double>({m_nnode, m_ndim});

  xt::xtensor<double,1> x = xt::linspace<double>(0.0, m_h*static_cast<double>(m_nelx), m_nelx+1);
  xt::xtensor<double,1> y = xt::linspace<double>(0.0, m_h*static_cast<double>(m_nely), m_nely+1);
  xt::xtensor<double,1> z = xt::linspace<double>(0.0, m_h*static_cast<double>(m_nelz), m_nelz+1);

  size_t inode = 0;

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz ) {
    for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy ) {
      for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix ) {
        out(inode,0) = x(ix);
        out(inode,1) = y(iy);
        out(inode,2) = z(iz);
        ++inode;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> Regular::conn() const
{
  xt::xtensor<size_t,2> out = xt::empty<size_t>({m_nelem,m_nne});

  size_t ielem = 0;

  for ( size_t iz = 0 ; iz < m_nelz ; ++iz ) {
    for ( size_t iy = 0 ; iy < m_nely ; ++iy ) {
      for ( size_t ix = 0 ; ix < m_nelx ; ++ix ) {
        out(ielem,0) = (iy  )*(m_nelx+1) + (ix  ) + (iz  )*(m_nely+1)*(m_nelx+1);
        out(ielem,1) = (iy  )*(m_nelx+1) + (ix+1) + (iz  )*(m_nely+1)*(m_nelx+1);
        out(ielem,3) = (iy+1)*(m_nelx+1) + (ix  ) + (iz  )*(m_nely+1)*(m_nelx+1);
        out(ielem,2) = (iy+1)*(m_nelx+1) + (ix+1) + (iz  )*(m_nely+1)*(m_nelx+1);
        out(ielem,4) = (iy  )*(m_nelx+1) + (ix  ) + (iz+1)*(m_nely+1)*(m_nelx+1);
        out(ielem,5) = (iy  )*(m_nelx+1) + (ix+1) + (iz+1)*(m_nely+1)*(m_nelx+1);
        out(ielem,7) = (iy+1)*(m_nelx+1) + (ix  ) + (iz+1)*(m_nely+1)*(m_nelx+1);
        out(ielem,6) = (iy+1)*(m_nelx+1) + (ix+1) + (iz+1)*(m_nely+1)*(m_nelx+1);
        ++ielem;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesFront() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx+1)*(m_nely+1)});

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
      out(iy*(m_nelx+1)+ix) = iy*(m_nelx+1) + ix;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBack() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx+1)*(m_nely+1)});

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
      out(iy*(m_nelx+1)+ix) = iy*(m_nelx+1) + ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesLeft() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nely+1)*(m_nelz+1)});

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
      out(iz*(m_nely+1)+iy) = iy*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesRight() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nely+1)*(m_nelz+1)});

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
      out(iz*(m_nely+1)+iy) = iy*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBottom() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx+1)*(m_nelz+1)});

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
      out(iz*(m_nelx+1)+ix) = ix + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesTop() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx+1)*(m_nelz+1)});

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
      out(iz*(m_nelx+1)+ix) = ix + m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontFace() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx-1)*(m_nely-1)});

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
      out((iy-1)*(m_nelx-1)+(ix-1)) = iy*(m_nelx+1) + ix;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBackFace() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx-1)*(m_nely-1)});

  for ( size_t iy = 1 ; iy < m_nely ; ++iy ) {
    for ( size_t ix = 1 ; ix < m_nelx ; ++ix ) {
      out((iy-1)*(m_nelx-1)+(ix-1)) = iy*(m_nelx+1) + ix + m_nelz*(m_nely+1)*(m_nelx+1);
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesLeftFace() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nely-1)*(m_nelz-1)});

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    for ( size_t iy = 1 ; iy < m_nely ; ++iy )
      out((iz-1)*(m_nely-1)+(iy-1)) = iy*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesRightFace() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nely-1)*(m_nelz-1)});

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    for ( size_t iy = 1 ; iy < m_nely ; ++iy )
      out((iz-1)*(m_nely-1)+(iy-1)) = iy*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomFace() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx-1)*(m_nelz-1)});

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
      out((iz-1)*(m_nelx-1)+(ix-1)) = ix + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesTopFace() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx-1)*(m_nelz-1)});

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
      out((iz-1)*(m_nelx-1)+(ix-1)) = ix + m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontBottomEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx+1});

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = ix;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontTopEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx+1});

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = ix + m_nely*(m_nelx+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontLeftEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely+1});

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontRightEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely+1});

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1) + m_nelx;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBackBottomEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx+1});

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBackTopEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx+1});

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = m_nely*(m_nelx+1) + ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBackLeftEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely+1});

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1) + m_nelz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBackRightEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely+1});

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1) + m_nelz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomLeftEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz+1});

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    out(iz) = iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomRightEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz+1});

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    out(iz) = iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesTopLeftEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz+1});

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    out(iz) = m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesTopRightEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz+1});

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    out(iz) = m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomFrontEdge() const { return nodesFrontBottomEdge(); }
inline xt::xtensor<size_t,1> Regular::nodesBottomBackEdge()  const { return nodesBackBottomEdge();  }
inline xt::xtensor<size_t,1> Regular::nodesTopFrontEdge()    const { return nodesFrontTopEdge();    }
inline xt::xtensor<size_t,1> Regular::nodesTopBackEdge()     const { return nodesBackTopEdge();     }
inline xt::xtensor<size_t,1> Regular::nodesLeftBottomEdge()  const { return nodesBottomLeftEdge();  }
inline xt::xtensor<size_t,1> Regular::nodesLeftFrontEdge()   const { return nodesFrontLeftEdge();   }
inline xt::xtensor<size_t,1> Regular::nodesLeftBackEdge()    const { return nodesBackLeftEdge();    }
inline xt::xtensor<size_t,1> Regular::nodesLeftTopEdge()     const { return nodesTopLeftEdge();     }
inline xt::xtensor<size_t,1> Regular::nodesRightBottomEdge() const { return nodesBottomRightEdge(); }
inline xt::xtensor<size_t,1> Regular::nodesRightTopEdge()    const { return nodesTopRightEdge();    }
inline xt::xtensor<size_t,1> Regular::nodesRightFrontEdge()  const { return nodesFrontRightEdge();  }
inline xt::xtensor<size_t,1> Regular::nodesRightBackEdge()   const { return nodesBackRightEdge();   }

// ------------------- node-numbers along the front-bottom edge, without corners -------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontBottomOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx-1});

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = ix;

  return out;
}

// -------------------- node-numbers along the front-top edge, without corners ---------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontTopOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx-1});

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = ix + m_nely*(m_nelx+1);

  return out;
}

// -------------------- node-numbers along the front-left edge, without corners --------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontLeftOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely-1});

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1);

  return out;
}

// ------------------- node-numbers along the front-right edge, without corners --------------------

inline xt::xtensor<size_t,1> Regular::nodesFrontRightOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely-1});

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1) + m_nelx;

  return out;
}

// ------------------- node-numbers along the back-bottom edge, without corners --------------------

inline xt::xtensor<size_t,1> Regular::nodesBackBottomOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx-1});

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBackTopOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx-1});

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = m_nely*(m_nelx+1) + ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// -------------------- node-numbers along the back-left edge, without corners ---------------------

inline xt::xtensor<size_t,1> Regular::nodesBackLeftOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely-1});

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1) + m_nelz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------- node-numbers along the back-right edge, without corners --------------------

inline xt::xtensor<size_t,1> Regular::nodesBackRightOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely-1});

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1) + m_nelz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// ------------------- node-numbers along the bottom-left edge, without corners --------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomLeftOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz-1});

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    out(iz-1) = iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// ------------------- node-numbers along the bottom-right edge, without corners -------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomRightOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz-1});

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    out(iz-1) = iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesTopLeftOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz-1});

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    out(iz-1) = m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------- node-numbers along the top-right edge, without corners ---------------------

inline xt::xtensor<size_t,1> Regular::nodesTopRightOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz-1});

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    out(iz-1) = m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomFrontOpenEdge() const { return nodesFrontBottomOpenEdge(); }
inline xt::xtensor<size_t,1> Regular::nodesBottomBackOpenEdge()  const { return nodesBackBottomOpenEdge();  }
inline xt::xtensor<size_t,1> Regular::nodesTopFrontOpenEdge()    const { return nodesFrontTopOpenEdge();    }
inline xt::xtensor<size_t,1> Regular::nodesTopBackOpenEdge()     const { return nodesBackTopOpenEdge();     }
inline xt::xtensor<size_t,1> Regular::nodesLeftBottomOpenEdge()  const { return nodesBottomLeftOpenEdge();  }
inline xt::xtensor<size_t,1> Regular::nodesLeftFrontOpenEdge()   const { return nodesFrontLeftOpenEdge();   }
inline xt::xtensor<size_t,1> Regular::nodesLeftBackOpenEdge()    const { return nodesBackLeftOpenEdge();    }
inline xt::xtensor<size_t,1> Regular::nodesLeftTopOpenEdge()     const { return nodesTopLeftOpenEdge();     }
inline xt::xtensor<size_t,1> Regular::nodesRightBottomOpenEdge() const { return nodesBottomRightOpenEdge(); }
inline xt::xtensor<size_t,1> Regular::nodesRightTopOpenEdge()    const { return nodesTopRightOpenEdge();    }
inline xt::xtensor<size_t,1> Regular::nodesRightFrontOpenEdge()  const { return nodesFrontRightOpenEdge();  }
inline xt::xtensor<size_t,1> Regular::nodesRightBackOpenEdge()   const { return nodesBackRightOpenEdge();   }

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesFrontBottomLeftCorner() const
{
  return 0;
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesFrontBottomRightCorner() const
{
  return m_nelx;
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesFrontTopLeftCorner() const
{
  return m_nely*(m_nelx+1);
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesFrontTopRightCorner() const
{
  return m_nely*(m_nelx+1) + m_nelx;
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesBackBottomLeftCorner() const
{
  return m_nelz*(m_nely+1)*(m_nelx+1);
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesBackBottomRightCorner() const
{
  return m_nelx + m_nelz*(m_nely+1)*(m_nelx+1);
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesBackTopLeftCorner() const
{
  return m_nely*(m_nelx+1) + m_nelz*(m_nely+1)*(m_nelx+1);
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesBackTopRightCorner() const
{
  return m_nely*(m_nelx+1) + m_nelx + m_nelz*(m_nely+1)*(m_nelx+1);
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesFrontLeftBottomCorner() const { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesBottomFrontLeftCorner() const { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesBottomLeftFrontCorner() const { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesLeftFrontBottomCorner() const { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesLeftBottomFrontCorner() const { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesFrontRightBottomCorner()const { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesBottomFrontRightCorner()const { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesBottomRightFrontCorner()const { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesRightFrontBottomCorner()const { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesRightBottomFrontCorner()const { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesFrontLeftTopCorner()    const { return nodesFrontTopLeftCorner();     }
inline size_t Regular::nodesTopFrontLeftCorner()    const { return nodesFrontTopLeftCorner();     }
inline size_t Regular::nodesTopLeftFrontCorner()    const { return nodesFrontTopLeftCorner();     }
inline size_t Regular::nodesLeftFrontTopCorner()    const { return nodesFrontTopLeftCorner();     }
inline size_t Regular::nodesLeftTopFrontCorner()    const { return nodesFrontTopLeftCorner();     }
inline size_t Regular::nodesFrontRightTopCorner()   const { return nodesFrontTopRightCorner();    }
inline size_t Regular::nodesTopFrontRightCorner()   const { return nodesFrontTopRightCorner();    }
inline size_t Regular::nodesTopRightFrontCorner()   const { return nodesFrontTopRightCorner();    }
inline size_t Regular::nodesRightFrontTopCorner()   const { return nodesFrontTopRightCorner();    }
inline size_t Regular::nodesRightTopFrontCorner()   const { return nodesFrontTopRightCorner();    }
inline size_t Regular::nodesBackLeftBottomCorner()  const { return nodesBackBottomLeftCorner();   }
inline size_t Regular::nodesBottomBackLeftCorner()  const { return nodesBackBottomLeftCorner();   }
inline size_t Regular::nodesBottomLeftBackCorner()  const { return nodesBackBottomLeftCorner();   }
inline size_t Regular::nodesLeftBackBottomCorner()  const { return nodesBackBottomLeftCorner();   }
inline size_t Regular::nodesLeftBottomBackCorner()  const { return nodesBackBottomLeftCorner();   }
inline size_t Regular::nodesBackRightBottomCorner() const { return nodesBackBottomRightCorner();  }
inline size_t Regular::nodesBottomBackRightCorner() const { return nodesBackBottomRightCorner();  }
inline size_t Regular::nodesBottomRightBackCorner() const { return nodesBackBottomRightCorner();  }
inline size_t Regular::nodesRightBackBottomCorner() const { return nodesBackBottomRightCorner();  }
inline size_t Regular::nodesRightBottomBackCorner() const { return nodesBackBottomRightCorner();  }
inline size_t Regular::nodesBackLeftTopCorner()     const { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesTopBackLeftCorner()     const { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesTopLeftBackCorner()     const { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesLeftBackTopCorner()     const { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesLeftTopBackCorner()     const { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesBackRightTopCorner()    const { return nodesBackTopRightCorner();     }
inline size_t Regular::nodesTopBackRightCorner()    const { return nodesBackTopRightCorner();     }
inline size_t Regular::nodesTopRightBackCorner()    const { return nodesBackTopRightCorner();     }
inline size_t Regular::nodesRightBackTopCorner()    const { return nodesBackTopRightCorner();     }
inline size_t Regular::nodesRightTopBackCorner()    const { return nodesBackTopRightCorner();     }

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> Regular::nodesPeriodic() const
{
  // faces
  xt::xtensor<size_t,1> fro = nodesFrontFace();
  xt::xtensor<size_t,1> bck = nodesBackFace();
  xt::xtensor<size_t,1> lft = nodesLeftFace();
  xt::xtensor<size_t,1> rgt = nodesRightFace();
  xt::xtensor<size_t,1> bot = nodesBottomFace();
  xt::xtensor<size_t,1> top = nodesTopFace();

  // edges
  xt::xtensor<size_t,1> froBot = nodesFrontBottomOpenEdge();
  xt::xtensor<size_t,1> froTop = nodesFrontTopOpenEdge();
  xt::xtensor<size_t,1> froLft = nodesFrontLeftOpenEdge();
  xt::xtensor<size_t,1> froRgt = nodesFrontRightOpenEdge();
  xt::xtensor<size_t,1> bckBot = nodesBackBottomOpenEdge();
  xt::xtensor<size_t,1> bckTop = nodesBackTopOpenEdge();
  xt::xtensor<size_t,1> bckLft = nodesBackLeftOpenEdge();
  xt::xtensor<size_t,1> bckRgt = nodesBackRightOpenEdge();
  xt::xtensor<size_t,1> botLft = nodesBottomLeftOpenEdge();
  xt::xtensor<size_t,1> botRgt = nodesBottomRightOpenEdge();
  xt::xtensor<size_t,1> topLft = nodesTopLeftOpenEdge();
  xt::xtensor<size_t,1> topRgt = nodesTopRightOpenEdge();

  // allocate nodal ties
  // - number of tying per category
  size_t tface = fro.size() + lft.size() + bot.size();
  size_t tedge = 3*froBot.size() + 3*froLft.size() + 3*botLft.size();
  size_t tnode = 7;
  // - allocate
  xt::xtensor<size_t,2> out = xt::empty<size_t>({tface+tedge+tnode, std::size_t(2)});

  // counter
  size_t i = 0;

  // tie all corners to one corner
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesFrontBottomRightCorner(); ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesBackBottomRightCorner();  ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesBackBottomLeftCorner();   ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesFrontTopLeftCorner();     ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesFrontTopRightCorner();    ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesBackTopRightCorner();     ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesBackTopLeftCorner();      ++i;

  // tie all corresponding edges to each other (exclude corners)
  for ( size_t j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = bckBot(j); ++i; }
  for ( size_t j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = bckTop(j); ++i; }
  for ( size_t j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = froTop(j); ++i; }
  for ( size_t j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = botRgt(j); ++i; }
  for ( size_t j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = topRgt(j); ++i; }
  for ( size_t j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = topLft(j); ++i; }
  for ( size_t j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = froRgt(j); ++i; }
  for ( size_t j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = bckRgt(j); ++i; }
  for ( size_t j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = bckLft(j); ++i; }

  // tie faces to each-other
  for ( size_t j = 0 ; j<fro.size()    ; ++j ){ out(i,0) = fro(j);    out(i,1) = bck(j);    ++i; }
  for ( size_t j = 0 ; j<lft.size()    ; ++j ){ out(i,0) = lft(j);    out(i,1) = rgt(j);    ++i; }
  for ( size_t j = 0 ; j<bot.size()    ; ++j ){ out(i,0) = bot(j);    out(i,1) = top(j);    ++i; }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline size_t Regular::nodesOrigin() const
{
  return nodesFrontBottomLeftCorner();
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> Regular::dofs() const
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> Regular::dofsPeriodic() const
{
  // DOF-numbers for each component of each node (sequential)
  xt::xtensor<size_t,2> out = GooseFEM::Mesh::dofs(m_nnode,m_ndim);

  // periodic node-pairs
  xt::xtensor<size_t,2> nodePer = nodesPeriodic();

  // eliminate 'dependent' DOFs; renumber "out" to be sequential for the remaining DOFs
  for (size_t i = 0; i < nodePer.shape(0); ++i)
    for (size_t j = 0; j < m_ndim; ++j)
      out(nodePer(i,1),j) = out(nodePer(i,0),j);

  // renumber "out" to be sequential
  return GooseFEM::Mesh::renumber(out);
}

// -------------------------------------------------------------------------------------------------

inline FineLayer::FineLayer(size_t nelx, size_t nely, size_t nelz, double h, size_t nfine):
m_h(h)
{
  // basic assumptions
  GOOSEFEM_ASSERT(nelx >= 1ul);
  GOOSEFEM_ASSERT(nely >= 1ul);
  GOOSEFEM_ASSERT(nelz >= 1ul);

  // store basic info
  m_Lx = m_h * static_cast<double>(nelx);
  m_Lz = m_h * static_cast<double>(nelz);

  // compute element size in y-direction (use symmetry, compute upper half)
  // -------------------------------------------------------------------------------------------------

  // temporary variables
  size_t nmin, ntot;
  xt::xtensor<size_t,1> nhx    =      xt::ones<size_t>({nely});
  xt::xtensor<size_t,1> nhy    =      xt::ones<size_t>({nely});
  xt::xtensor<size_t,1> nhz    =      xt::ones<size_t>({nely});
  xt::xtensor<int   ,1> refine = -1 * xt::ones<int>   ({nely});

  // minimum height in y-direction (half of the height because of symmetry)
  if ( nely  % 2 == 0 ) nmin  =  nely    /2;
  else                  nmin  = (nely +1)/2;

  // minimum number of fine layers in y-direction (minimum 1, middle layer part of this half)
  if ( nfine % 2 == 0 ) nfine =  nfine   /2 + 1;
  else                  nfine = (nfine+1)/2;
  if ( nfine < 1      ) nfine = 1;
  if ( nfine > nmin   ) nfine = nmin;

  // loop over element layers in y-direction, try to coarsen using these rules:
  // (1) element size in y-direction <= distance to origin in y-direction
  // (2) element size in x-(z-)direction should fit the total number of elements in x-(z-)direction
  // (3) a certain number of layers have the minimum size "1" (are fine)
  for ( size_t iy = nfine ; ; )
  {
    // initialize current size in y-direction
    if ( iy == nfine ) ntot = nfine;
    // check to stop
    if ( iy >= nely or ntot >= nmin ) { nely = iy; break; }

    // rules (1,2) satisfied: coarsen in x-direction (and z-direction)
    if ( 3*nhy(iy) <= ntot and nelx%(3*nhx(iy)) == 0 and ntot+nhy(iy) < nmin )
    {
      // - process refinement in x-direction
      refine(iy)  = 0;
      nhy   (iy) *= 2;
      auto vnhy = xt::view(nhy, xt::range(iy+1, _));
      auto vnhx = xt::view(nhx, xt::range(iy  , _));
      vnhy *= 3;
      vnhx *= 3;

      // - rule (2) satisfied: coarsen next element layer in z-direction
      if ( iy+1 < nely and ntot+2*nhy(iy) < nmin )
      {
        if ( nelz%(3*nhz(iy+1)) == 0 )
        {
          // - update the number of elements in y-direction
          ntot += nhy(iy);
          // - proceed to next element layer in y-direction
          ++iy;
          // - process refinement in z-direction
          refine(iy) = 2;
          nhy   (iy) = nhy(iy-1);
          auto vnhz = xt::view(nhz, xt::range(iy, _));
          vnhz *= 3;
        }
      }
    }

    // rules (1,2) satisfied: coarse in z-direction
    else if ( 3*nhy(iy) <= ntot and nelz%(3*nhz(iy)) == 0 and ntot+nhy(iy) < nmin )
    {
      // - process refinement in z-direction
      refine(iy)  = 2;
      nhy   (iy) *= 2;
      auto vnhy = xt::view(nhy, xt::range(iy+1, _));
      auto vnhz = xt::view(nhz, xt::range(iy  , _));
      vnhy *= 3;
      vnhz *= 3;
    }

    // update the number of elements in y-direction
    ntot += nhy(iy);
    // proceed to next element layer in y-direction
    ++iy;
    // check to stop
    if ( iy >= nely or ntot >= nmin ) { nely = iy; break; }
  }

  // symmetrize, compute full information
  // -------------------------------------------------------------------------------------------------

  // allocate mesh constructor parameters
  m_nhx       = xt::empty<size_t>({nely*2-1});
  m_nhy       = xt::empty<size_t>({nely*2-1});
  m_nhz       = xt::empty<size_t>({nely*2-1});
  m_refine    = xt::empty<int>   ({nely*2-1});
  m_nelx      = xt::empty<size_t>({nely*2-1});
  m_nelz      = xt::empty<size_t>({nely*2-1});
  m_nnd       = xt::empty<size_t>({nely*2  });
  m_startElem = xt::empty<size_t>({nely*2-1});
  m_startNode = xt::empty<size_t>({nely*2  });

  // fill
  // - lower half
  for ( size_t iy = 0 ; iy < nely ; ++iy )
  {
    m_nhx   (iy) = nhx   (nely-iy-1);
    m_nhy   (iy) = nhy   (nely-iy-1);
    m_nhz   (iy) = nhz   (nely-iy-1);
    m_refine(iy) = refine(nely-iy-1);
  }
  // - upper half
  for ( size_t iy = 0 ; iy < nely-1 ; ++iy )
  {
    m_nhx   (iy+nely) = nhx   (iy+1);
    m_nhy   (iy+nely) = nhy   (iy+1);
    m_nhz   (iy+nely) = nhz   (iy+1);
    m_refine(iy+nely) = refine(iy+1);
  }

  // update size
  nely = m_nhx.size();

  // compute the number of elements per element layer in y-direction
  for ( size_t iy = 0 ; iy < nely ; ++iy )
  {
    m_nelx(iy) = nelx / m_nhx(iy);
    m_nelz(iy) = nelz / m_nhz(iy);
  }

  // compute the number of nodes per node layer in y-direction
  // - bottom half
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    m_nnd(iy) = (m_nelx(iy)+1) * (m_nelz(iy)+1);
  // - top half
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    m_nnd(iy+1) = (m_nelx(iy)+1) * (m_nelz(iy)+1);

  // compute mesh dimensions
  // -------------------------------------------------------------------------------------------------

  // initialize
  m_nnode        = 0;
  m_nelem        = 0;
  m_startNode(0) = 0;

  // loop over element layers (bottom -> middle, elements become finer)
  for ( size_t i = 0 ; i < (nely-1)/2 ; ++i )
  {
    // - store the first element of the layer
    m_startElem(i) = m_nelem;
    // - add the nodes of this layer
    if      ( m_refine(i) == 0 ) { m_nnode += (3*m_nelx(i)+1) * (  m_nelz(i)+1); }
    else if ( m_refine(i) == 2 ) { m_nnode += (  m_nelx(i)+1) * (3*m_nelz(i)+1); }
    else                         { m_nnode += (  m_nelx(i)+1) * (  m_nelz(i)+1); }
    // - add the elements of this layer
    if      ( m_refine(i) == 0 ) { m_nelem += (4*m_nelx(i)  ) * (  m_nelz(i)  ); }
    else if ( m_refine(i) == 2 ) { m_nelem += (  m_nelx(i)  ) * (4*m_nelz(i)  ); }
    else                         { m_nelem += (  m_nelx(i)  ) * (  m_nelz(i)  ); }
    // - store the starting node of the next layer
    m_startNode(i+1) = m_nnode;
  }

  // loop over element layers (middle -> top, elements become coarser)
  for ( size_t i = (nely-1)/2 ; i < nely ; ++i )
  {
    // - store the first element of the layer
    m_startElem(i) = m_nelem;
    // - add the nodes of this layer
    if      ( m_refine(i) == 0 ) { m_nnode += (5*m_nelx(i)+1) * (  m_nelz(i)+1); }
    else if ( m_refine(i) == 2 ) { m_nnode += (  m_nelx(i)+1) * (5*m_nelz(i)+1); }
    else                         { m_nnode += (  m_nelx(i)+1) * (  m_nelz(i)+1); }
    // - add the elements of this layer
    if      ( m_refine(i) == 0 ) { m_nelem += (4*m_nelx(i)  ) * (  m_nelz(i)  ); }
    else if ( m_refine(i) == 2 ) { m_nelem += (  m_nelx(i)  ) * (4*m_nelz(i)  ); }
    else                         { m_nelem += (  m_nelx(i)  ) * (  m_nelz(i)  ); }
    // - store the starting node of the next layer
    m_startNode(i+1) = m_nnode;
  }
  // - add the top row of nodes
  m_nnode += (m_nelx(nely-1)+1) * (m_nelz(nely-1)+1);
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nelem() const
{
  return m_nelem;
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nnode() const
{
  return m_nnode;
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nne() const
{
  return m_nne;
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::shape(size_t i) const
{
  GOOSEFEM_ASSERT(i <= 2ul);

  if      ( i == 0 ) return xt::amax(m_nelx)[0];
  else if ( i == 2 ) return xt::amax(m_nelz)[0];
  else               return xt::sum (m_nhy )[0];
}

// -------------------------------------------------------------------------------------------------

inline ElementType FineLayer::getElementType() const
{
  return ElementType::Hex8;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> FineLayer::coor() const
{
  // allocate output
  xt::xtensor<double,2> out = xt::empty<double>({m_nnode, m_ndim});

  // current node, number of element layers
  size_t inode = 0;
  size_t nely  = static_cast<size_t>(m_nhy.size());

  // y-position of each main node layer (i.e. excluding node layers for refinement/coarsening)
  // - allocate
  xt::xtensor<double,1> y = xt::empty<double>({nely+1});
  // - initialize
  y(0) = 0.0;
  // - compute
  for ( size_t iy = 1 ; iy < nely+1 ; ++iy )
    y(iy) = y(iy-1) + m_nhy(iy-1) * m_h;

  // loop over element layers (bottom -> middle) : add bottom layer (+ refinement layer) of nodes
  // -------------------------------------------------------------------------------------------------

  for ( size_t iy = 0 ; ; ++iy )
  {
    // get positions along the x- and z-axis
    xt::xtensor<double,1> x = xt::linspace<double>(0.0, m_Lx, m_nelx(iy)+1);
    xt::xtensor<double,1> z = xt::linspace<double>(0.0, m_Lz, m_nelz(iy)+1);

    // add nodes of the bottom layer of this element
    for ( size_t iz = 0 ; iz < m_nelz(iy)+1 ; ++iz ) {
      for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
        out(inode,0) = x(ix);
        out(inode,1) = y(iy);
        out(inode,2) = z(iz);
        ++inode;
      }
    }

    // stop at middle layer
    if ( iy == (nely-1)/2 )
      break;

    // add extra nodes of the intermediate layer, for refinement in x-direction
    if ( m_refine(iy) == 0 )
    {
      // - get position offset in x- and y-direction
      double dx = m_h * static_cast<double>(m_nhx(iy)/3);
      double dy = m_h * static_cast<double>(m_nhy(iy)/2);
      // - add nodes of the intermediate layer
      for ( size_t iz = 0 ; iz < m_nelz(iy)+1 ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
          for ( size_t j = 0 ; j < 2 ; ++j ) {
            out(inode,0) = x(ix) + dx * static_cast<double>(j+1);
            out(inode,1) = y(iy) + dy;
            out(inode,2) = z(iz);
            ++inode;
          }
        }
      }
    }

    // add extra nodes of the intermediate layer, for refinement in z-direction
    else if ( m_refine(iy) == 2 )
    {
      // - get position offset in y- and z-direction
      double dz = m_h * static_cast<double>(m_nhz(iy)/3);
      double dy = m_h * static_cast<double>(m_nhy(iy)/2);
      // - add nodes of the intermediate layer
      for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz ) {
        for ( size_t j = 0 ; j < 2 ; ++j ) {
          for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
            out(inode,0) = x(ix);
            out(inode,1) = y(iy) + dy;
            out(inode,2) = z(iz) + dz * static_cast<double>(j+1);
            ++inode;
          }
        }
      }
    }
  }

  // loop over element layers (middle -> top) : add (refinement layer +) top layer of nodes
  // -------------------------------------------------------------------------------------------------

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    // get positions along the x- and z-axis
    xt::xtensor<double,1> x = xt::linspace<double>(0.0, m_Lx, m_nelx(iy)+1);
    xt::xtensor<double,1> z = xt::linspace<double>(0.0, m_Lz, m_nelz(iy)+1);

    // add extra nodes of the intermediate layer, for refinement in x-direction
    if ( m_refine(iy) == 0 )
    {
      // - get position offset in x- and y-direction
      double dx = m_h * static_cast<double>(m_nhx(iy)/3);
      double dy = m_h * static_cast<double>(m_nhy(iy)/2);
      // - add nodes of the intermediate layer
      for ( size_t iz = 0 ; iz < m_nelz(iy)+1 ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
          for ( size_t j = 0 ; j < 2 ; ++j ) {
            out(inode,0) = x(ix) + dx * static_cast<double>(j+1);
            out(inode,1) = y(iy) + dy;
            out(inode,2) = z(iz);
            ++inode;
          }
        }
      }
    }

    // add extra nodes of the intermediate layer, for refinement in z-direction
    else if ( m_refine(iy) == 2 )
    {
      // - get position offset in y- and z-direction
      double dz = m_h * static_cast<double>(m_nhz(iy)/3);
      double dy = m_h * static_cast<double>(m_nhy(iy)/2);
      // - add nodes of the intermediate layer
      for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz ) {
        for ( size_t j = 0 ; j < 2 ; ++j ) {
          for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
            out(inode,0) = x(ix);
            out(inode,1) = y(iy) + dy;
            out(inode,2) = z(iz) + dz * static_cast<double>(j+1);
            ++inode;
          }
        }
      }
    }

    // add nodes of the top layer of this element
    for ( size_t iz = 0 ; iz < m_nelz(iy)+1 ; ++iz ) {
      for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
        out(inode,0) = x(ix  );
        out(inode,1) = y(iy+1);
        out(inode,2) = z(iz  );
        ++inode;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> FineLayer::conn() const
{
  // allocate output
  xt::xtensor<size_t,2> out = xt::empty<size_t>({m_nelem, m_nne});

  // current element, number of element layers, starting nodes of each node layer
  size_t ielem = 0;
  size_t nely  = static_cast<size_t>(m_nhy.size());
  size_t bot,mid,top;

  // loop over all element layers
  for ( size_t iy = 0 ; iy < nely ; ++iy )
  {
    // - get: starting nodes of bottom(, middle) and top layer
    bot = m_startNode(iy  );
    mid = m_startNode(iy  ) + m_nnd(iy);
    top = m_startNode(iy+1);

    // - define connectivity: no coarsening/refinement
    if ( m_refine(iy) == -1 )
    {
      for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
          out(ielem,0) = bot + (ix  ) + (iz  ) * (m_nelx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (iz  ) * (m_nelx(iy)+1);
          out(ielem,2) = top + (ix+1) + (iz  ) * (m_nelx(iy)+1);
          out(ielem,3) = top + (ix  ) + (iz  ) * (m_nelx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (iz+1) * (m_nelx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (iz+1) * (m_nelx(iy)+1);
          out(ielem,6) = top + (ix+1) + (iz+1) * (m_nelx(iy)+1);
          out(ielem,7) = top + (ix  ) + (iz+1) * (m_nelx(iy)+1);
          ielem++;
        }
      }
    }

    // - define connectivity: refinement along the x-direction (below the middle layer)
    else if ( m_refine(iy) == 0 and iy <= (nely-1)/2 )
    {
      for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
          // -- bottom element
          out(ielem,0) = bot + (  ix  ) + (iz  ) * (  m_nelx(iy)+1);
          out(ielem,1) = bot + (  ix+1) + (iz  ) * (  m_nelx(iy)+1);
          out(ielem,2) = mid + (2*ix+1) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,3) = mid + (2*ix  ) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,4) = bot + (  ix  ) + (iz+1) * (  m_nelx(iy)+1);
          out(ielem,5) = bot + (  ix+1) + (iz+1) * (  m_nelx(iy)+1);
          out(ielem,6) = mid + (2*ix+1) + (iz+1) * (2*m_nelx(iy)  );
          out(ielem,7) = mid + (2*ix  ) + (iz+1) * (2*m_nelx(iy)  );
          ielem++;
          // -- top-right element
          out(ielem,0) = bot + (  ix+1) + (iz  ) * (  m_nelx(iy)+1);
          out(ielem,1) = top + (3*ix+3) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,2) = top + (3*ix+2) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,3) = mid + (2*ix+1) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,4) = bot + (  ix+1) + (iz+1) * (  m_nelx(iy)+1);
          out(ielem,5) = top + (3*ix+3) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,6) = top + (3*ix+2) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,7) = mid + (2*ix+1) + (iz+1) * (2*m_nelx(iy)  );
          ielem++;
          // -- top-center element
          out(ielem,0) = mid + (2*ix  ) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,1) = mid + (2*ix+1) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,2) = top + (3*ix+2) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,3) = top + (3*ix+1) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,4) = mid + (2*ix  ) + (iz+1) * (2*m_nelx(iy)  );
          out(ielem,5) = mid + (2*ix+1) + (iz+1) * (2*m_nelx(iy)  );
          out(ielem,6) = top + (3*ix+2) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,7) = top + (3*ix+1) + (iz+1) * (3*m_nelx(iy)+1);
          ielem++;
          // -- top-left element
          out(ielem,0) = bot + (  ix  ) + (iz  ) * (  m_nelx(iy)+1);
          out(ielem,1) = mid + (2*ix  ) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,2) = top + (3*ix+1) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,3) = top + (3*ix  ) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,4) = bot + (  ix  ) + (iz+1) * (  m_nelx(iy)+1);
          out(ielem,5) = mid + (2*ix  ) + (iz+1) * (2*m_nelx(iy)  );
          out(ielem,6) = top + (3*ix+1) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,7) = top + (3*ix  ) + (iz+1) * (3*m_nelx(iy)+1);
          ielem++;
        }
      }
    }

    // - define connectivity: coarsening along the x-direction (above the middle layer)
    else if ( m_refine(iy) == 0 and iy > (nely-1)/2 )
    {
      for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
          // -- lower-left element
          out(ielem,0) = bot + (3*ix  ) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,1) = bot + (3*ix+1) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,2) = mid + (2*ix  ) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,3) = top + (  ix  ) + (iz  ) * (  m_nelx(iy)+1);
          out(ielem,4) = bot + (3*ix  ) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,5) = bot + (3*ix+1) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,6) = mid + (2*ix  ) + (iz+1) * (2*m_nelx(iy)  );
          out(ielem,7) = top + (  ix  ) + (iz+1) * (  m_nelx(iy)+1);
          ielem++;
          // -- lower-center element
          out(ielem,0) = bot + (3*ix+1) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,1) = bot + (3*ix+2) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,2) = mid + (2*ix+1) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,3) = mid + (2*ix  ) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,4) = bot + (3*ix+1) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,5) = bot + (3*ix+2) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,6) = mid + (2*ix+1) + (iz+1) * (2*m_nelx(iy)  );
          out(ielem,7) = mid + (2*ix  ) + (iz+1) * (2*m_nelx(iy)  );
          ielem++;
          // -- lower-right element
          out(ielem,0) = bot + (3*ix+2) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,1) = bot + (3*ix+3) + (iz  ) * (3*m_nelx(iy)+1);
          out(ielem,2) = top + (  ix+1) + (iz  ) * (  m_nelx(iy)+1);
          out(ielem,3) = mid + (2*ix+1) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,4) = bot + (3*ix+2) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,5) = bot + (3*ix+3) + (iz+1) * (3*m_nelx(iy)+1);
          out(ielem,6) = top + (  ix+1) + (iz+1) * (  m_nelx(iy)+1);
          out(ielem,7) = mid + (2*ix+1) + (iz+1) * (2*m_nelx(iy)  );
          ielem++;
          // -- upper element
          out(ielem,0) = mid + (2*ix  ) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,1) = mid + (2*ix+1) + (iz  ) * (2*m_nelx(iy)  );
          out(ielem,2) = top + (  ix+1) + (iz  ) * (  m_nelx(iy)+1);
          out(ielem,3) = top + (  ix  ) + (iz  ) * (  m_nelx(iy)+1);
          out(ielem,4) = mid + (2*ix  ) + (iz+1) * (2*m_nelx(iy)  );
          out(ielem,5) = mid + (2*ix+1) + (iz+1) * (2*m_nelx(iy)  );
          out(ielem,6) = top + (  ix+1) + (iz+1) * (  m_nelx(iy)+1);
          out(ielem,7) = top + (  ix  ) + (iz+1) * (  m_nelx(iy)+1);
          ielem++;
        }
      }
    }

    // - define connectivity: refinement along the z-direction (below the middle layer)
    else if ( m_refine(iy) == 2 and iy <= (nely-1)/2 )
    {
      for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
          // -- bottom element
          out(ielem,0) = bot + (ix  ) +    iz    * (m_nelx(iy)+1);
          out(ielem,1) = bot + (ix  ) + (  iz+1) * (m_nelx(iy)+1);
          out(ielem,2) = bot + (ix+1) + (  iz+1) * (m_nelx(iy)+1);
          out(ielem,3) = bot + (ix+1) +    iz    * (m_nelx(iy)+1);
          out(ielem,4) = mid + (ix  ) +  2*iz    * (m_nelx(iy)+1);
          out(ielem,5) = mid + (ix  ) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,6) = mid + (ix+1) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,7) = mid + (ix+1) +  2*iz    * (m_nelx(iy)+1);
          ielem++;
          // -- top-back element
          out(ielem,0) = mid + (ix  ) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,1) = mid + (ix+1) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,2) = top + (ix+1) + (3*iz+2) * (m_nelx(iy)+1);
          out(ielem,3) = top + (ix  ) + (3*iz+2) * (m_nelx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (  iz+1) * (m_nelx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (  iz+1) * (m_nelx(iy)+1);
          out(ielem,6) = top + (ix+1) + (3*iz+3) * (m_nelx(iy)+1);
          out(ielem,7) = top + (ix  ) + (3*iz+3) * (m_nelx(iy)+1);
          ielem++;
          // -- top-center element
          out(ielem,0) = mid + (ix  ) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,1) = mid + (ix+1) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,2) = top + (ix+1) + (3*iz+1) * (m_nelx(iy)+1);
          out(ielem,3) = top + (ix  ) + (3*iz+1) * (m_nelx(iy)+1);
          out(ielem,4) = mid + (ix  ) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,5) = mid + (ix+1) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,6) = top + (ix+1) + (3*iz+2) * (m_nelx(iy)+1);
          out(ielem,7) = top + (ix  ) + (3*iz+2) * (m_nelx(iy)+1);
          ielem++;
          // -- top-front element
          out(ielem,0) = bot + (ix  ) + (  iz  ) * (m_nelx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (  iz  ) * (m_nelx(iy)+1);
          out(ielem,2) = top + (ix+1) + (3*iz  ) * (m_nelx(iy)+1);
          out(ielem,3) = top + (ix  ) + (3*iz  ) * (m_nelx(iy)+1);
          out(ielem,4) = mid + (ix  ) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,5) = mid + (ix+1) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,6) = top + (ix+1) + (3*iz+1) * (m_nelx(iy)+1);
          out(ielem,7) = top + (ix  ) + (3*iz+1) * (m_nelx(iy)+1);
          ielem++;
        }
      }
    }

    // - define connectivity: coarsening along the z-direction (above the middle layer)
    else if ( m_refine(iy) == 2 and iy > (nely-1)/2 )
    {
      for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
          // -- bottom-front element
          out(ielem,0) = bot + (ix  ) + (3*iz  ) * (m_nelx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (3*iz  ) * (m_nelx(iy)+1);
          out(ielem,2) = top + (ix+1) + (  iz  ) * (m_nelx(iy)+1);
          out(ielem,3) = top + (ix  ) + (  iz  ) * (m_nelx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (3*iz+1) * (m_nelx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (3*iz+1) * (m_nelx(iy)+1);
          out(ielem,6) = mid + (ix+1) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,7) = mid + (ix  ) + (2*iz  ) * (m_nelx(iy)+1);
          ielem++;
          // -- bottom-center element
          out(ielem,0) = bot + (ix  ) + (3*iz+1) * (m_nelx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (3*iz+1) * (m_nelx(iy)+1);
          out(ielem,2) = mid + (ix+1) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,3) = mid + (ix  ) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (3*iz+2) * (m_nelx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (3*iz+2) * (m_nelx(iy)+1);
          out(ielem,6) = mid + (ix+1) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,7) = mid + (ix  ) + (2*iz+1) * (m_nelx(iy)+1);
          ielem++;
          // -- bottom-back element
          out(ielem,0) = bot + (ix  ) + (3*iz+2) * (m_nelx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (3*iz+2) * (m_nelx(iy)+1);
          out(ielem,2) = mid + (ix+1) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,3) = mid + (ix  ) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (3*iz+3) * (m_nelx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (3*iz+3) * (m_nelx(iy)+1);
          out(ielem,6) = top + (ix+1) + (  iz+1) * (m_nelx(iy)+1);
          out(ielem,7) = top + (ix  ) + (  iz+1) * (m_nelx(iy)+1);
          ielem++;
          // -- top element
          out(ielem,0) = mid + (ix  ) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,1) = mid + (ix+1) + (2*iz  ) * (m_nelx(iy)+1);
          out(ielem,2) = top + (ix+1) + (  iz  ) * (m_nelx(iy)+1);
          out(ielem,3) = top + (ix  ) + (  iz  ) * (m_nelx(iy)+1);
          out(ielem,4) = mid + (ix  ) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,5) = mid + (ix+1) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,6) = top + (ix+1) + (  iz+1) * (m_nelx(iy)+1);
          out(ielem,7) = top + (ix  ) + (  iz+1) * (m_nelx(iy)+1);
          ielem++;
        }
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::elementsMiddleLayer() const
{
  // number of element layers in y-direction, the index of the middle layer
  size_t nely = static_cast<size_t>(m_nhy.size());
  size_t iy   = (nely-1)/2;

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(iy)*m_nelz(iy)});

  for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix )
    for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz )
      out(ix+iz*m_nelx(iy)) = m_startElem(iy) + ix + iz*m_nelx(iy);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFront() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // number of boundary nodes
  // - initialize
  size_t n = 0;
  // - bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
  {
    if ( m_refine(iy) == 0 ) n += m_nelx(iy) * 3 + 1;
    else                     n += m_nelx(iy)     + 1;
  }
  // - top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    if ( m_refine(iy) == 0 ) n += m_nelx(iy) * 3 + 1;
    else                     n += m_nelx(iy)     + 1;
  }

  // allocate node-list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({n});

  // initialize counter: current index in the node-list "out"
  size_t j = 0;

  // bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
  {
    // -- bottom node layer
    for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
      out(j) = m_startNode(iy) + ix;
      ++j;
    }
    // -- refinement layer
    if ( m_refine(iy) == 0 ) {
      for ( size_t ix = 0 ; ix < 2*m_nelx(iy) ; ++ix ) {
        out(j) = m_startNode(iy) + ix + m_nnd(iy);
        ++j;
      }
    }
  }

  // top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    // -- refinement layer
    if ( m_refine(iy) == 0 ) {
      for ( size_t ix = 0 ; ix < 2*m_nelx(iy) ; ++ix ) {
        out(j) = m_startNode(iy) + ix + m_nnd(iy);
        ++j;
      }
    }
    // -- top node layer
    for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
      out(j) = m_startNode(iy+1) + ix;
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBack() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // number of boundary nodes
  // - initialize
  size_t n = 0;
  // - bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
  {
    if ( m_refine(iy) == 0 ) n += m_nelx(iy) * 3 + 1;
    else                     n += m_nelx(iy)     + 1;
  }
  // - top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    if ( m_refine(iy) == 0 ) n += m_nelx(iy) * 3 + 1;
    else                     n += m_nelx(iy)     + 1;
  }

  // allocate node-list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({n});

  // initialize counter: current index in the node-list "out"
  size_t j = 0;

  // bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
  {
    // -- bottom node layer
    for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
      out(j) = m_startNode(iy) + ix + (m_nelx(iy)+1)*m_nelz(iy);
      ++j;
    }
    // -- refinement layer
    if ( m_refine(iy) == 0 ) {
      for ( size_t ix = 0 ; ix < 2*m_nelx(iy) ; ++ix ) {
        out(j) = m_startNode(iy) + ix + m_nnd(iy) + 2*m_nelx(iy)*m_nelz(iy);
        ++j;
      }
    }
  }

  // top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    // -- refinement layer
    if ( m_refine(iy) == 0 ) {
      for ( size_t ix = 0 ; ix < 2*m_nelx(iy) ; ++ix ) {
        out(j) = m_startNode(iy) + ix + m_nnd(iy) + 2*m_nelx(iy)*m_nelz(iy);
        ++j;
      }
    }
    // -- top node layer
    for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
      out(j) = m_startNode(iy+1) + ix + (m_nelx(iy)+1)*m_nelz(iy);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesLeft() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // number of boundary nodes
  // - initialize
  size_t n = 0;
  // - bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
  {
    if ( m_refine(iy) == 2 ) n += m_nelz(iy) * 3 + 1;
    else                     n += m_nelz(iy)     + 1;
  }
  // - top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    if ( m_refine(iy) == 2 ) n += m_nelz(iy) * 3 + 1;
    else                     n += m_nelz(iy)     + 1;
  }

  // allocate node-list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({n});

  // initialize counter: current index in the node-list "out"
  size_t j = 0;

  // bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
  {
    // -- bottom node layer
    for ( size_t iz = 0 ; iz < m_nelz(iy)+1 ; ++iz ) {
      out(j) = m_startNode(iy) + iz * (m_nelx(iy)+1);
      ++j;
    }
    // -- refinement layer
    if ( m_refine(iy) == 2 ) {
      for ( size_t iz = 0 ; iz < 2*m_nelz(iy) ; ++iz ) {
        out(j) = m_startNode(iy) + iz * (m_nelx(iy)+1) + m_nnd(iy);
        ++j;
      }
    }
  }

  // top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    // -- refinement layer
    if ( m_refine(iy) == 2 ) {
      for ( size_t iz = 0 ; iz < 2*m_nelz(iy) ; ++iz ) {
        out(j) = m_startNode(iy) + iz * (m_nelx(iy)+1) + m_nnd(iy);
        ++j;
      }
    }
    // -- top node layer
    for ( size_t iz = 0 ; iz < m_nelz(iy)+1 ; ++iz ) {
      out(j) = m_startNode(iy+1) + iz * (m_nelx(iy)+1);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesRight() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // number of boundary nodes
  // - initialize
  size_t n = 0;
  // - bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
  {
    if ( m_refine(iy) == 2 ) n += m_nelz(iy) * 3 + 1;
    else                     n += m_nelz(iy)     + 1;
  }
  // - top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    if ( m_refine(iy) == 2 ) n += m_nelz(iy) * 3 + 1;
    else                     n += m_nelz(iy)     + 1;
  }

  // allocate node-list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({n});

  // initialize counter: current index in the node-list "out"
  size_t j = 0;

  // bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
  {
    // -- bottom node layer
    for ( size_t iz = 0 ; iz < m_nelz(iy)+1 ; ++iz ) {
      out(j) = m_startNode(iy) + iz * (m_nelx(iy)+1) + m_nelx(iy);
      ++j;
    }
    // -- refinement layer
    if ( m_refine(iy) == 2 ) {
      for ( size_t iz = 0 ; iz < 2*m_nelz(iy) ; ++iz ) {
        out(j) = m_startNode(iy) + m_nnd(iy) + iz * (m_nelx(iy)+1) + m_nelx(iy);
        ++j;
      }
    }
  }

  // top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    // -- refinement layer
    if ( m_refine(iy) == 2 ) {
      for ( size_t iz = 0 ; iz < 2*m_nelz(iy) ; ++iz ) {
        out(j) = m_startNode(iy) + m_nnd(iy) + iz * (m_nelx(iy)+1) + m_nelx(iy);
        ++j;
      }
    }
    // -- top node layer
    for ( size_t iz = 0 ; iz < m_nelz(iy)+1 ; ++iz ) {
      out(j) = m_startNode(iy+1) + iz * (m_nelx(iy)+1) + m_nelx(iy);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottom() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // allocate node list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nnd(nely)});

  // counter
  size_t j = 0;

  // fill node list
  for ( size_t ix = 0 ; ix < m_nelx(0)+1 ; ++ix ) {
    for ( size_t iz = 0 ; iz < m_nelz(0)+1 ; ++iz ) {
      out(j) = m_startNode(0) + ix + iz * (m_nelx(0)+1);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesTop() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // allocate node list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nnd(nely)});

  // counter
  size_t j = 0;

  // fill node list
  for ( size_t ix = 0 ; ix < m_nelx(nely-1)+1 ; ++ix ) {
    for ( size_t iz = 0 ; iz < m_nelz(nely-1)+1 ; ++iz ) {
      out(j) = m_startNode(nely) + ix + iz * (m_nelx(nely-1)+1);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontFace() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // number of boundary nodes
  // - initialize
  size_t n = 0;
  // - bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
  {
    if ( m_refine(iy) == 0 ) n += m_nelx(iy) * 3 - 1;
    else                     n += m_nelx(iy)     - 1;
  }
  // - top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
  {
    if ( m_refine(iy) == 0 ) n += m_nelx(iy) * 3 - 1;
    else                     n += m_nelx(iy)     - 1;
  }

  // allocate node-list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({n});

  // initialize counter: current index in the node-list "out"
  size_t j = 0;

  // bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
  {
    // -- bottom node layer
    for ( size_t ix = 1 ; ix < m_nelx(iy) ; ++ix ) {
      out(j) = m_startNode(iy) + ix;
      ++j;
    }
    // -- refinement layer
    if ( m_refine(iy) == 0 ) {
      for ( size_t ix = 0 ; ix < 2*m_nelx(iy) ; ++ix ) {
        out(j) = m_startNode(iy) + ix + m_nnd(iy);
        ++j;
      }
    }
  }

  // top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
  {
    // -- refinement layer
    if ( m_refine(iy) == 0 ) {
      for ( size_t ix = 0 ; ix < 2*m_nelx(iy) ; ++ix ) {
        out(j) = m_startNode(iy) + ix + m_nnd(iy);
        ++j;
      }
    }
    // -- top node layer
    for ( size_t ix = 1 ; ix < m_nelx(iy) ; ++ix ) {
      out(j) = m_startNode(iy+1) + ix;
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackFace() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // number of boundary nodes
  // - initialize
  size_t n = 0;
  // - bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
  {
    if ( m_refine(iy) == 0 ) n += m_nelx(iy) * 3 - 1;
    else                     n += m_nelx(iy)     - 1;
  }
  // - top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
  {
    if ( m_refine(iy) == 0 ) n += m_nelx(iy) * 3 - 1;
    else                     n += m_nelx(iy)     - 1;
  }

  // allocate node-list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({n});

  // initialize counter: current index in the node-list "out"
  size_t j = 0;

  // bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
  {
    // -- bottom node layer
    for ( size_t ix = 1 ; ix < m_nelx(iy) ; ++ix ) {
      out(j) = m_startNode(iy) + ix + (m_nelx(iy)+1)*m_nelz(iy);
      ++j;
    }
    // -- refinement layer
    if ( m_refine(iy) == 0 ) {
      for ( size_t ix = 0 ; ix < 2*m_nelx(iy) ; ++ix ) {
        out(j) = m_startNode(iy) + ix + m_nnd(iy) + 2*m_nelx(iy)*m_nelz(iy);
        ++j;
      }
    }
  }

  // top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
  {
    // -- refinement layer
    if ( m_refine(iy) == 0 ) {
      for ( size_t ix = 0 ; ix < 2*m_nelx(iy) ; ++ix ) {
        out(j) = m_startNode(iy) + ix + m_nnd(iy) + 2*m_nelx(iy)*m_nelz(iy);
        ++j;
      }
    }
    // -- top node layer
    for ( size_t ix = 1 ; ix < m_nelx(iy) ; ++ix ) {
      out(j) = m_startNode(iy+1) + ix + (m_nelx(iy)+1)*m_nelz(iy);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesLeftFace() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // number of boundary nodes
  // - initialize
  size_t n = 0;
  // - bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
  {
    if ( m_refine(iy) == 2 ) n += m_nelz(iy) * 3 - 1;
    else                     n += m_nelz(iy)     - 1;
  }
  // - top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
  {
    if ( m_refine(iy) == 2 ) n += m_nelz(iy) * 3 - 1;
    else                     n += m_nelz(iy)     - 1;
  }

  // allocate node-list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({n});

  // initialize counter: current index in the node-list "out"
  size_t j = 0;

  // bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
  {
    // -- bottom node layer
    for ( size_t iz = 1 ; iz < m_nelz(iy) ; ++iz ) {
      out(j) = m_startNode(iy) + iz * (m_nelx(iy)+1);
      ++j;
    }
    // -- refinement layer
    if ( m_refine(iy) == 2 ) {
      for ( size_t iz = 0 ; iz < 2*m_nelz(iy) ; ++iz ) {
        out(j) = m_startNode(iy) + iz * (m_nelx(iy)+1) + m_nnd(iy);
        ++j;
      }
    }
  }

  // top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
  {
    // -- refinement layer
    if ( m_refine(iy) == 2 ) {
      for ( size_t iz = 0 ; iz < 2*m_nelz(iy) ; ++iz ) {
        out(j) = m_startNode(iy) + iz * (m_nelx(iy)+1) + m_nnd(iy);
        ++j;
      }
    }
    // -- top node layer
    for ( size_t iz = 1 ; iz < m_nelz(iy) ; ++iz ) {
      out(j) = m_startNode(iy+1) + iz * (m_nelx(iy)+1);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesRightFace() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // number of boundary nodes
  // - initialize
  size_t n = 0;
  // - bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
  {
    if ( m_refine(iy) == 2 ) n += m_nelz(iy) * 3 - 1;
    else                     n += m_nelz(iy)     - 1;
  }
  // - top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
  {
    if ( m_refine(iy) == 2 ) n += m_nelz(iy) * 3 - 1;
    else                     n += m_nelz(iy)     - 1;
  }

  // allocate node-list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({n});

  // initialize counter: current index in the node-list "out"
  size_t j = 0;

  // bottom half: bottom node layer (+ middle node layer)
  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
  {
    // -- bottom node layer
    for ( size_t iz = 1 ; iz < m_nelz(iy) ; ++iz ) {
      out(j) = m_startNode(iy) + iz * (m_nelx(iy)+1) + m_nelx(iy);
      ++j;
    }
    // -- refinement layer
    if ( m_refine(iy) == 2 ) {
      for ( size_t iz = 0 ; iz < 2*m_nelz(iy) ; ++iz ) {
        out(j) = m_startNode(iy) + m_nnd(iy) + iz * (m_nelx(iy)+1) + m_nelx(iy);
        ++j;
      }
    }
  }

  // top half: (middle node layer +) top node layer
  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
  {
    // -- refinement layer
    if ( m_refine(iy) == 2 ) {
      for ( size_t iz = 0 ; iz < 2*m_nelz(iy) ; ++iz ) {
        out(j) = m_startNode(iy) + m_nnd(iy) + iz * (m_nelx(iy)+1) + m_nelx(iy);
        ++j;
      }
    }
    // -- top node layer
    for ( size_t iz = 1 ; iz < m_nelz(iy) ; ++iz ) {
      out(j) = m_startNode(iy+1) + iz * (m_nelx(iy)+1) + m_nelx(iy);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomFace() const
{
  // allocate node list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx(0)-1)*(m_nelz(0)-1)});

  // counter
  size_t j = 0;

  // fill node list
  for ( size_t ix = 1 ; ix < m_nelx(0) ; ++ix ) {
    for ( size_t iz = 1 ; iz < m_nelz(0) ; ++iz ) {
      out(j) = m_startNode(0) + ix + iz * (m_nelx(0)+1);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesTopFace() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // allocate node list
  xt::xtensor<size_t,1> out = xt::empty<size_t>({(m_nelx(nely-1)-1)*(m_nelz(nely-1)-1)});

  // counter
  size_t j = 0;

  // fill node list
  for ( size_t ix = 1 ; ix < m_nelx(nely-1) ; ++ix ) {
    for ( size_t iz = 1 ; iz < m_nelz(nely-1) ; ++iz ) {
      out(j) = m_startNode(nely) + ix + iz * (m_nelx(nely-1)+1);
      ++j;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontBottomEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(0)+1});

  for ( size_t ix = 0 ; ix < m_nelx(0)+1 ; ++ix )
    out(ix) = m_startNode(0) + ix;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontTopEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(nely-1)+1});

  for ( size_t ix = 0 ; ix < m_nelx(nely-1)+1 ; ++ix )
    out(ix) = m_startNode(nely) + ix;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontLeftEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely+1});

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontRightEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely+1});

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy) + m_nelx(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1) + m_nelx(iy);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackBottomEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(0)+1});

  for ( size_t ix = 0 ; ix < m_nelx(0)+1 ; ++ix )
    out(ix) = m_startNode(0) + ix + (m_nelx(0)+1)*(m_nelz(0));

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackTopEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(nely-1)+1});

  for ( size_t ix = 0 ; ix < m_nelx(nely-1)+1 ; ++ix )
    out(ix) = m_startNode(nely) + ix + (m_nelx(nely-1)+1)*(m_nelz(nely-1));

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackLeftEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely+1});

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1) + (m_nelx(iy)+1)*(m_nelz(iy));

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackRightEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely+1});

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy) + m_nelx(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1) + m_nelx(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomLeftEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz(0)+1});

  for ( size_t iz = 0 ; iz < m_nelz(0)+1 ; ++iz )
    out(iz) = m_startNode(0) + iz * (m_nelx(0)+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomRightEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz(0)+1});

  for ( size_t iz = 0 ; iz < m_nelz(0)+1 ; ++iz )
    out(iz) = m_startNode(0) + m_nelx(0) + iz * (m_nelx(0)+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesTopLeftEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz(nely-1)+1});

  for ( size_t iz = 0 ; iz < m_nelz(nely-1)+1 ; ++iz )
    out(iz) = m_startNode(nely) + iz * (m_nelx(nely-1)+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesTopRightEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz(nely-1)+1});

  for ( size_t iz = 0 ; iz < m_nelz(nely-1)+1 ; ++iz )
    out(iz) = m_startNode(nely) + m_nelx(nely-1) + iz * (m_nelx(nely-1)+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomFrontEdge() const { return nodesFrontBottomEdge(); }
inline xt::xtensor<size_t,1> FineLayer::nodesBottomBackEdge()  const { return nodesBackBottomEdge();  }
inline xt::xtensor<size_t,1> FineLayer::nodesTopFrontEdge()    const { return nodesFrontTopEdge();    }
inline xt::xtensor<size_t,1> FineLayer::nodesTopBackEdge()     const { return nodesBackTopEdge();     }
inline xt::xtensor<size_t,1> FineLayer::nodesLeftBottomEdge()  const { return nodesBottomLeftEdge();  }
inline xt::xtensor<size_t,1> FineLayer::nodesLeftFrontEdge()   const { return nodesFrontLeftEdge();   }
inline xt::xtensor<size_t,1> FineLayer::nodesLeftBackEdge()    const { return nodesBackLeftEdge();    }
inline xt::xtensor<size_t,1> FineLayer::nodesLeftTopEdge()     const { return nodesTopLeftEdge();     }
inline xt::xtensor<size_t,1> FineLayer::nodesRightBottomEdge() const { return nodesBottomRightEdge(); }
inline xt::xtensor<size_t,1> FineLayer::nodesRightTopEdge()    const { return nodesTopRightEdge();    }
inline xt::xtensor<size_t,1> FineLayer::nodesRightFrontEdge()  const { return nodesFrontRightEdge();  }
inline xt::xtensor<size_t,1> FineLayer::nodesRightBackEdge()   const { return nodesBackRightEdge();   }

// ------------------- node-numbers along the front-bottom edge, without corners -------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontBottomOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(0)-1});

  for ( size_t ix = 1 ; ix < m_nelx(0) ; ++ix )
    out(ix-1) = m_startNode(0) + ix;

  return out;
}

// -------------------- node-numbers along the front-top edge, without corners ---------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontTopOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(nely-1)-1});

  for ( size_t ix = 1 ; ix < m_nelx(nely-1) ; ++ix )
    out(ix-1) = m_startNode(nely) + ix;

  return out;
}

// -------------------- node-numbers along the front-left edge, without corners --------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontLeftOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely-1});

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1);

  return out;
}

// ------------------- node-numbers along the front-right edge, without corners --------------------

inline xt::xtensor<size_t,1> FineLayer::nodesFrontRightOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely-1});

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy) + m_nelx(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1) + m_nelx(iy);

  return out;
}

// ------------------- node-numbers along the back-bottom edge, without corners --------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackBottomOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(0)-1});

  for ( size_t ix = 1 ; ix < m_nelx(0) ; ++ix )
    out(ix-1) = m_startNode(0) + ix + (m_nelx(0)+1)*(m_nelz(0));

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackTopOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(nely-1)-1});

  for ( size_t ix = 1 ; ix < m_nelx(nely-1) ; ++ix )
    out(ix-1) = m_startNode(nely) + ix + (m_nelx(nely-1)+1)*(m_nelz(nely-1));

  return out;
}

// -------------------- node-numbers along the back-left edge, without corners ---------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackLeftOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely-1});

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1) + (m_nelx(iy)+1)*(m_nelz(iy));

  return out;
}

// -------------------- node-numbers along the back-right edge, without corners --------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBackRightOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely-1});

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy) + m_nelx(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1) + m_nelx(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  return out;
}

// ------------------- node-numbers along the bottom-left edge, without corners --------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomLeftOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz(0)-1});

  for ( size_t iz = 1 ; iz < m_nelz(0) ; ++iz )
    out(iz-1) = m_startNode(0) + iz * (m_nelx(0)+1);

  return out;
}

// ------------------- node-numbers along the bottom-right edge, without corners -------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomRightOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz(0)-1});

  for ( size_t iz = 1 ; iz < m_nelz(0) ; ++iz )
    out(iz-1) = m_startNode(0) + m_nelx(0) + iz * (m_nelx(0)+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesTopLeftOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz(nely-1)-1});

  for ( size_t iz = 1 ; iz < m_nelz(nely-1) ; ++iz )
    out(iz-1) = m_startNode(nely) + iz * (m_nelx(nely-1)+1);

  return out;
}

// -------------------- node-numbers along the top-right edge, without corners ---------------------

inline xt::xtensor<size_t,1> FineLayer::nodesTopRightOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelz(nely-1)-1});

  for ( size_t iz = 1 ; iz < m_nelz(nely-1) ; ++iz )
    out(iz-1) = m_startNode(nely) + m_nelx(nely-1) + iz * (m_nelx(nely-1)+1);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomFrontOpenEdge() const { return nodesFrontBottomOpenEdge(); }
inline xt::xtensor<size_t,1> FineLayer::nodesBottomBackOpenEdge() const  { return nodesBackBottomOpenEdge();  }
inline xt::xtensor<size_t,1> FineLayer::nodesTopFrontOpenEdge() const    { return nodesFrontTopOpenEdge();    }
inline xt::xtensor<size_t,1> FineLayer::nodesTopBackOpenEdge() const     { return nodesBackTopOpenEdge();     }
inline xt::xtensor<size_t,1> FineLayer::nodesLeftBottomOpenEdge() const  { return nodesBottomLeftOpenEdge();  }
inline xt::xtensor<size_t,1> FineLayer::nodesLeftFrontOpenEdge() const   { return nodesFrontLeftOpenEdge();   }
inline xt::xtensor<size_t,1> FineLayer::nodesLeftBackOpenEdge() const    { return nodesBackLeftOpenEdge();    }
inline xt::xtensor<size_t,1> FineLayer::nodesLeftTopOpenEdge() const     { return nodesTopLeftOpenEdge();     }
inline xt::xtensor<size_t,1> FineLayer::nodesRightBottomOpenEdge() const { return nodesBottomRightOpenEdge(); }
inline xt::xtensor<size_t,1> FineLayer::nodesRightTopOpenEdge() const    { return nodesTopRightOpenEdge();    }
inline xt::xtensor<size_t,1> FineLayer::nodesRightFrontOpenEdge() const  { return nodesFrontRightOpenEdge();  }
inline xt::xtensor<size_t,1> FineLayer::nodesRightBackOpenEdge() const   { return nodesBackRightOpenEdge();   }

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesFrontBottomLeftCorner() const
{
  return m_startNode(0);
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesFrontBottomRightCorner() const
{
  return m_startNode(0) + m_nelx(0);
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesFrontTopLeftCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely);
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesFrontTopRightCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely) + m_nelx(nely-1);
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesBackBottomLeftCorner() const
{
  return m_startNode(0) + (m_nelx(0)+1)*(m_nelz(0));
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesBackBottomRightCorner() const
{
  return m_startNode(0) + m_nelx(0) + (m_nelx(0)+1)*(m_nelz(0));
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesBackTopLeftCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely) + (m_nelx(nely-1)+1)*(m_nelz(nely-1));
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesBackTopRightCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely) + m_nelx(nely-1) + (m_nelx(nely-1)+1)*(m_nelz(nely-1));
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesFrontLeftBottomCorner() const   { return nodesFrontBottomLeftCorner();  }
inline size_t FineLayer::nodesBottomFrontLeftCorner() const   { return nodesFrontBottomLeftCorner();  }
inline size_t FineLayer::nodesBottomLeftFrontCorner() const   { return nodesFrontBottomLeftCorner();  }
inline size_t FineLayer::nodesLeftFrontBottomCorner() const   { return nodesFrontBottomLeftCorner();  }
inline size_t FineLayer::nodesLeftBottomFrontCorner() const   { return nodesFrontBottomLeftCorner();  }
inline size_t FineLayer::nodesFrontRightBottomCorner() const  { return nodesFrontBottomRightCorner(); }
inline size_t FineLayer::nodesBottomFrontRightCorner() const  { return nodesFrontBottomRightCorner(); }
inline size_t FineLayer::nodesBottomRightFrontCorner() const  { return nodesFrontBottomRightCorner(); }
inline size_t FineLayer::nodesRightFrontBottomCorner() const  { return nodesFrontBottomRightCorner(); }
inline size_t FineLayer::nodesRightBottomFrontCorner() const  { return nodesFrontBottomRightCorner(); }
inline size_t FineLayer::nodesFrontLeftTopCorner() const      { return nodesFrontTopLeftCorner();     }
inline size_t FineLayer::nodesTopFrontLeftCorner() const      { return nodesFrontTopLeftCorner();     }
inline size_t FineLayer::nodesTopLeftFrontCorner() const      { return nodesFrontTopLeftCorner();     }
inline size_t FineLayer::nodesLeftFrontTopCorner() const      { return nodesFrontTopLeftCorner();     }
inline size_t FineLayer::nodesLeftTopFrontCorner() const      { return nodesFrontTopLeftCorner();     }
inline size_t FineLayer::nodesFrontRightTopCorner() const     { return nodesFrontTopRightCorner();    }
inline size_t FineLayer::nodesTopFrontRightCorner() const     { return nodesFrontTopRightCorner();    }
inline size_t FineLayer::nodesTopRightFrontCorner() const     { return nodesFrontTopRightCorner();    }
inline size_t FineLayer::nodesRightFrontTopCorner() const     { return nodesFrontTopRightCorner();    }
inline size_t FineLayer::nodesRightTopFrontCorner() const     { return nodesFrontTopRightCorner();    }
inline size_t FineLayer::nodesBackLeftBottomCorner() const    { return nodesBackBottomLeftCorner();   }
inline size_t FineLayer::nodesBottomBackLeftCorner() const    { return nodesBackBottomLeftCorner();   }
inline size_t FineLayer::nodesBottomLeftBackCorner() const    { return nodesBackBottomLeftCorner();   }
inline size_t FineLayer::nodesLeftBackBottomCorner() const    { return nodesBackBottomLeftCorner();   }
inline size_t FineLayer::nodesLeftBottomBackCorner() const    { return nodesBackBottomLeftCorner();   }
inline size_t FineLayer::nodesBackRightBottomCorner() const   { return nodesBackBottomRightCorner();  }
inline size_t FineLayer::nodesBottomBackRightCorner() const   { return nodesBackBottomRightCorner();  }
inline size_t FineLayer::nodesBottomRightBackCorner() const   { return nodesBackBottomRightCorner();  }
inline size_t FineLayer::nodesRightBackBottomCorner() const   { return nodesBackBottomRightCorner();  }
inline size_t FineLayer::nodesRightBottomBackCorner() const   { return nodesBackBottomRightCorner();  }
inline size_t FineLayer::nodesBackLeftTopCorner() const       { return nodesBackTopLeftCorner();      }
inline size_t FineLayer::nodesTopBackLeftCorner() const       { return nodesBackTopLeftCorner();      }
inline size_t FineLayer::nodesTopLeftBackCorner() const       { return nodesBackTopLeftCorner();      }
inline size_t FineLayer::nodesLeftBackTopCorner() const       { return nodesBackTopLeftCorner();      }
inline size_t FineLayer::nodesLeftTopBackCorner() const       { return nodesBackTopLeftCorner();      }
inline size_t FineLayer::nodesBackRightTopCorner() const      { return nodesBackTopRightCorner();     }
inline size_t FineLayer::nodesTopBackRightCorner() const      { return nodesBackTopRightCorner();     }
inline size_t FineLayer::nodesTopRightBackCorner() const      { return nodesBackTopRightCorner();     }
inline size_t FineLayer::nodesRightBackTopCorner() const      { return nodesBackTopRightCorner();     }
inline size_t FineLayer::nodesRightTopBackCorner() const      { return nodesBackTopRightCorner();     }

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> FineLayer::nodesPeriodic() const
{
  // faces
  xt::xtensor<size_t,1> fro = nodesFrontFace();
  xt::xtensor<size_t,1> bck = nodesBackFace();
  xt::xtensor<size_t,1> lft = nodesLeftFace();
  xt::xtensor<size_t,1> rgt = nodesRightFace();
  xt::xtensor<size_t,1> bot = nodesBottomFace();
  xt::xtensor<size_t,1> top = nodesTopFace();

  // edges
  xt::xtensor<size_t,1> froBot = nodesFrontBottomOpenEdge();
  xt::xtensor<size_t,1> froTop = nodesFrontTopOpenEdge();
  xt::xtensor<size_t,1> froLft = nodesFrontLeftOpenEdge();
  xt::xtensor<size_t,1> froRgt = nodesFrontRightOpenEdge();
  xt::xtensor<size_t,1> bckBot = nodesBackBottomOpenEdge();
  xt::xtensor<size_t,1> bckTop = nodesBackTopOpenEdge();
  xt::xtensor<size_t,1> bckLft = nodesBackLeftOpenEdge();
  xt::xtensor<size_t,1> bckRgt = nodesBackRightOpenEdge();
  xt::xtensor<size_t,1> botLft = nodesBottomLeftOpenEdge();
  xt::xtensor<size_t,1> botRgt = nodesBottomRightOpenEdge();
  xt::xtensor<size_t,1> topLft = nodesTopLeftOpenEdge();
  xt::xtensor<size_t,1> topRgt = nodesTopRightOpenEdge();

  // allocate nodal ties
  // - number of tying per category
  size_t tface = fro.size() + lft.size() + bot.size();
  size_t tedge = 3*froBot.size() + 3*froLft.size() + 3*botLft.size();
  size_t tnode = 7;
  // - allocate
  xt::xtensor<size_t,2> out = xt::empty<size_t>({tface+tedge+tnode, std::size_t(2)});

  // counter
  size_t i = 0;

  // tie all corners to one corner
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesFrontBottomRightCorner(); ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesBackBottomRightCorner();  ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesBackBottomLeftCorner();   ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesFrontTopLeftCorner();     ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesFrontTopRightCorner();    ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesBackTopRightCorner();     ++i;
  out(i,0) = nodesFrontBottomLeftCorner(); out(i,1) = nodesBackTopLeftCorner();      ++i;

  // tie all corresponding edges to each other (exclude corners)
  for ( size_t j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = bckBot(j); ++i; }
  for ( size_t j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = bckTop(j); ++i; }
  for ( size_t j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = froTop(j); ++i; }
  for ( size_t j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = botRgt(j); ++i; }
  for ( size_t j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = topRgt(j); ++i; }
  for ( size_t j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = topLft(j); ++i; }
  for ( size_t j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = froRgt(j); ++i; }
  for ( size_t j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = bckRgt(j); ++i; }
  for ( size_t j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = bckLft(j); ++i; }

  // tie faces to each-other
  for ( size_t j = 0 ; j<fro.size()    ; ++j ){ out(i,0) = fro(j);    out(i,1) = bck(j);    ++i; }
  for ( size_t j = 0 ; j<lft.size()    ; ++j ){ out(i,0) = lft(j);    out(i,1) = rgt(j);    ++i; }
  for ( size_t j = 0 ; j<bot.size()    ; ++j ){ out(i,0) = bot(j);    out(i,1) = top(j);    ++i; }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::nodesOrigin() const
{
  return nodesFrontBottomLeftCorner();
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> FineLayer::dofs() const
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<size_t,2> FineLayer::dofsPeriodic() const
{
  // DOF-numbers for each component of each node (sequential)
  xt::xtensor<size_t,2> out = GooseFEM::Mesh::dofs(m_nnode,m_ndim);

  // periodic node-pairs
  xt::xtensor<size_t,2>   nodePer = nodesPeriodic();

  // eliminate 'dependent' DOFs; renumber "out" to be sequential for the remaining DOFs
  for (size_t i = 0; i < nodePer.shape(0); ++i)
    for (size_t j = 0; j < m_ndim; ++j)
      out(nodePer(i,1),j) = out(nodePer(i,0),j);

  // renumber "out" to be sequential
  return GooseFEM::Mesh::renumber(out);
}

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHHEX8_CPP
#define GOOSEFEM_MESHHEX8_CPP

// -------------------------------------------------------------------------------------------------

#include "MeshHex8.h"

// ===================================== GooseFEM::Mesh::Hex8 ======================================

namespace GooseFEM {
namespace Mesh {
namespace Hex8 {

// ------------------------------------------ constructor ------------------------------------------

inline Regular::Regular(size_t nelx, size_t nely, size_t nelz, double h):
m_h(h), m_nelx(nelx), m_nely(nely), m_nelz(nelz)
{
  assert( m_nelx >= 1 );
  assert( m_nely >= 1 );
  assert( m_nelz >= 1 );

  m_nnode = (m_nelx+1) * (m_nely+1) * (m_nelz+1);
  m_nelem =  m_nelx    *  m_nely    *  m_nelz   ;
}

// -------------------------------------- number of elements ---------------------------------------

inline size_t Regular::nelem() const
{
  return m_nelem;
}

// ---------------------------------------- number of nodes ----------------------------------------

inline size_t Regular::nnode() const
{
  return m_nnode;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t Regular::nne() const
{
  return m_nne;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t Regular::ndim() const
{
  return m_ndim;
}

// --------------------------------- coordinates (nodal positions) ---------------------------------

inline MatD Regular::coor() const
{
  MatD out(m_nnode,m_ndim);

  ColD x = ColD::LinSpaced(m_nelx+1, 0.0, m_h*static_cast<double>(m_nelx));
  ColD y = ColD::LinSpaced(m_nely+1, 0.0, m_h*static_cast<double>(m_nely));
  ColD z = ColD::LinSpaced(m_nelz+1, 0.0, m_h*static_cast<double>(m_nelz));

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

// ---------------------------- connectivity (node-numbers per element) ----------------------------

inline MatS Regular::conn() const
{
  MatS out(m_nelem,m_nne);

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

// ------------------------------ node-numbers along the front plane -------------------------------

inline ColS Regular::nodesFront() const
{
  ColS out((m_nelx+1)*(m_nely+1));

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
      out(iy*(m_nelx+1)+ix) = iy*(m_nelx+1) + ix;

  return out;
}

// ------------------------------- node-numbers along the back plane -------------------------------

inline ColS Regular::nodesBack() const
{
  ColS out((m_nelx+1)*(m_nely+1));

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
      out(iy*(m_nelx+1)+ix) = iy*(m_nelx+1) + ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// ------------------------------- node-numbers along the left plane -------------------------------

inline ColS Regular::nodesLeft() const
{
  ColS out((m_nely+1)*(m_nelz+1));

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
      out(iz*(m_nely+1)+iy) = iy*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// ------------------------------ node-numbers along the right plane -------------------------------

inline ColS Regular::nodesRight() const
{
  ColS out((m_nely+1)*(m_nelz+1));

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
      out(iz*(m_nely+1)+iy) = iy*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// ------------------------------ node-numbers along the bottom plane ------------------------------

inline ColS Regular::nodesBottom() const
{
  ColS out((m_nelx+1)*(m_nelz+1));

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
      out(iz*(m_nelx+1)+ix) = ix + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// ------------------------------- node-numbers along the top plane --------------------------------

inline ColS Regular::nodesTop() const
{
  ColS out((m_nelx+1)*(m_nelz+1));

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
      out(iz*(m_nelx+1)+ix) = ix + m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// ------------------------------- node-numbers along the front face -------------------------------

inline ColS Regular::nodesFrontFace() const
{
  ColS out((m_nelx-1)*(m_nely-1));

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
      out((iy-1)*(m_nelx-1)+(ix-1)) = iy*(m_nelx+1) + ix;

  return out;
}

// ------------------------------- node-numbers along the back face --------------------------------

inline ColS Regular::nodesBackFace() const
{
  ColS out((m_nelx-1)*(m_nely-1));

  for ( size_t iy = 1 ; iy < m_nely ; ++iy ) {
    for ( size_t ix = 1 ; ix < m_nelx ; ++ix ) {
      out((iy-1)*(m_nelx-1)+(ix-1)) = iy*(m_nelx+1) + ix + m_nelz*(m_nely+1)*(m_nelx+1);
    }
  }

  return out;
}

// ------------------------------- node-numbers along the left face --------------------------------

inline ColS Regular::nodesLeftFace() const
{
  ColS out((m_nely-1)*(m_nelz-1));

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    for ( size_t iy = 1 ; iy < m_nely ; ++iy )
      out((iz-1)*(m_nely-1)+(iy-1)) = iy*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// ------------------------------- node-numbers along the right face -------------------------------

inline ColS Regular::nodesRightFace() const
{
  ColS out((m_nely-1)*(m_nelz-1));

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    for ( size_t iy = 1 ; iy < m_nely ; ++iy )
      out((iz-1)*(m_nely-1)+(iy-1)) = iy*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// ------------------------------ node-numbers along the bottom face -------------------------------

inline ColS Regular::nodesBottomFace() const
{
  ColS out((m_nelx-1)*(m_nelz-1));

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
      out((iz-1)*(m_nelx-1)+(ix-1)) = ix + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------------------- node-numbers along the top face --------------------------------

inline ColS Regular::nodesTopFace() const
{
  ColS out((m_nelx-1)*(m_nelz-1));

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
      out((iz-1)*(m_nelx-1)+(ix-1)) = ix + m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// --------------------------- node-numbers along the front-bottom edge ----------------------------

inline ColS Regular::nodesFrontBottomEdge() const
{
  ColS out(m_nelx+1);

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = ix;

  return out;
}

// ----------------------------- node-numbers along the front-top edge -----------------------------

inline ColS Regular::nodesFrontTopEdge() const
{
  ColS out(m_nelx+1);

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = ix + m_nely*(m_nelx+1);

  return out;
}

// ---------------------------- node-numbers along the front-left edge -----------------------------

inline ColS Regular::nodesFrontLeftEdge() const
{
  ColS out(m_nely+1);

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1);

  return out;
}

// ---------------------------- node-numbers along the front-right edge ----------------------------

inline ColS Regular::nodesFrontRightEdge() const
{
  ColS out(m_nely+1);

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1) + m_nelx;

  return out;
}

// ---------------------------- node-numbers along the back-bottom edge ----------------------------

inline ColS Regular::nodesBackBottomEdge() const
{
  ColS out(m_nelx+1);

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// ----------------------------- node-numbers along the back-top edge ------------------------------

inline ColS Regular::nodesBackTopEdge() const
{
  ColS out(m_nelx+1);

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = m_nely*(m_nelx+1) + ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// ----------------------------- node-numbers along the back-left edge -----------------------------

inline ColS Regular::nodesBackLeftEdge() const
{
  ColS out(m_nely+1);

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1) + m_nelz*(m_nelx+1)*(m_nely+1);

  return out;
}

// ---------------------------- node-numbers along the back-right edge -----------------------------

inline ColS Regular::nodesBackRightEdge() const
{
  ColS out(m_nely+1);

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1) + m_nelz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// ---------------------------- node-numbers along the bottom-left edge ----------------------------

inline ColS Regular::nodesBottomLeftEdge() const
{
  ColS out(m_nelz+1);

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    out(iz) = iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// --------------------------- node-numbers along the bottom-right edge ----------------------------

inline ColS Regular::nodesBottomRightEdge() const
{
  ColS out(m_nelz+1);

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    out(iz) = iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// ----------------------------- node-numbers along the top-left edge ------------------------------

inline ColS Regular::nodesTopLeftEdge() const
{
  ColS out(m_nelz+1);

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    out(iz) = m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// ----------------------------- node-numbers along the top-right edge -----------------------------

inline ColS Regular::nodesTopRightEdge() const
{
  ColS out(m_nelz+1);

  for ( size_t iz = 0 ; iz < m_nelz+1 ; ++iz )
    out(iz) = m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------- aliases --------------------------------------------

inline ColS Regular::nodesBottomFrontEdge() const { return nodesFrontBottomEdge(); }
inline ColS Regular::nodesBottomBackEdge()  const { return nodesBackBottomEdge();  }
inline ColS Regular::nodesTopFrontEdge()    const { return nodesFrontTopEdge();    }
inline ColS Regular::nodesTopBackEdge()     const { return nodesBackTopEdge();     }
inline ColS Regular::nodesLeftBottomEdge()  const { return nodesBottomLeftEdge();  }
inline ColS Regular::nodesLeftFrontEdge()   const { return nodesFrontLeftEdge();   }
inline ColS Regular::nodesLeftBackEdge()    const { return nodesBackLeftEdge();    }
inline ColS Regular::nodesLeftTopEdge()     const { return nodesTopLeftEdge();     }
inline ColS Regular::nodesRightBottomEdge() const { return nodesBottomRightEdge(); }
inline ColS Regular::nodesRightTopEdge()    const { return nodesTopRightEdge();    }
inline ColS Regular::nodesRightFrontEdge()  const { return nodesFrontRightEdge();  }
inline ColS Regular::nodesRightBackEdge()   const { return nodesBackRightEdge();   }

// ------------------- node-numbers along the front-bottom edge, without corners -------------------

inline ColS Regular::nodesFrontBottomOpenEdge() const
{
  ColS out(m_nelx-1);

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = ix;

  return out;
}

// -------------------- node-numbers along the front-top edge, without corners ---------------------

inline ColS Regular::nodesFrontTopOpenEdge() const
{
  ColS out(m_nelx-1);

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = ix + m_nely*(m_nelx+1);

  return out;
}

// -------------------- node-numbers along the front-left edge, without corners --------------------

inline ColS Regular::nodesFrontLeftOpenEdge() const
{
  ColS out(m_nely-1);

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1);

  return out;
}

// ------------------- node-numbers along the front-right edge, without corners --------------------

inline ColS Regular::nodesFrontRightOpenEdge() const
{
  ColS out(m_nely-1);

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1) + m_nelx;

  return out;
}

// ------------------- node-numbers along the back-bottom edge, without corners --------------------

inline ColS Regular::nodesBackBottomOpenEdge() const
{
  ColS out(m_nelx-1);

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// --------------------- node-numbers along the back-top edge, without corners ---------------------

inline ColS Regular::nodesBackTopOpenEdge() const
{
  ColS out(m_nelx-1);

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = m_nely*(m_nelx+1) + ix + m_nelz*(m_nely+1)*(m_nelx+1);

  return out;
}

// -------------------- node-numbers along the back-left edge, without corners ---------------------

inline ColS Regular::nodesBackLeftOpenEdge() const
{
  ColS out(m_nely-1);

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1) + m_nelz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------- node-numbers along the back-right edge, without corners --------------------

inline ColS Regular::nodesBackRightOpenEdge() const
{
  ColS out(m_nely-1);

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1) + m_nelz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// ------------------- node-numbers along the bottom-left edge, without corners --------------------

inline ColS Regular::nodesBottomLeftOpenEdge() const
{
  ColS out(m_nelz-1);

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    out(iz-1) = iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// ------------------- node-numbers along the bottom-right edge, without corners -------------------

inline ColS Regular::nodesBottomRightOpenEdge() const
{
  ColS out(m_nelz-1);

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    out(iz-1) = iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// --------------------- node-numbers along the top-left edge, without corners ---------------------

inline ColS Regular::nodesTopLeftOpenEdge() const
{
  ColS out(m_nelz-1);

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    out(iz-1) = m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1);

  return out;
}

// -------------------- node-numbers along the top-right edge, without corners ---------------------

inline ColS Regular::nodesTopRightOpenEdge() const
{
  ColS out(m_nelz-1);

  for ( size_t iz = 1 ; iz < m_nelz ; ++iz )
    out(iz-1) = m_nely*(m_nelx+1) + iz*(m_nelx+1)*(m_nely+1) + m_nelx;

  return out;
}

// -------------------------------------------- aliases --------------------------------------------

inline ColS Regular::nodesBottomFrontOpenEdge() const { return nodesFrontBottomOpenEdge(); }
inline ColS Regular::nodesBottomBackOpenEdge()  const { return nodesBackBottomOpenEdge();  }
inline ColS Regular::nodesTopFrontOpenEdge()    const { return nodesFrontTopOpenEdge();    }
inline ColS Regular::nodesTopBackOpenEdge()     const { return nodesBackTopOpenEdge();     }
inline ColS Regular::nodesLeftBottomOpenEdge()  const { return nodesBottomLeftOpenEdge();  }
inline ColS Regular::nodesLeftFrontOpenEdge()   const { return nodesFrontLeftOpenEdge();   }
inline ColS Regular::nodesLeftBackOpenEdge()    const { return nodesBackLeftOpenEdge();    }
inline ColS Regular::nodesLeftTopOpenEdge()     const { return nodesTopLeftOpenEdge();     }
inline ColS Regular::nodesRightBottomOpenEdge() const { return nodesBottomRightOpenEdge(); }
inline ColS Regular::nodesRightTopOpenEdge()    const { return nodesTopRightOpenEdge();    }
inline ColS Regular::nodesRightFrontOpenEdge()  const { return nodesFrontRightOpenEdge();  }
inline ColS Regular::nodesRightBackOpenEdge()   const { return nodesBackRightOpenEdge();   }

// -------------------------- node-number of the front-bottom-left corner --------------------------

inline size_t Regular::nodesFrontBottomLeftCorner() const
{
  return 0;
}

// ------------------------- node-number of the front-bottom-right corner --------------------------

inline size_t Regular::nodesFrontBottomRightCorner() const
{
  return m_nelx;
}

// --------------------------- node-number of the front-top-left corner ----------------------------

inline size_t Regular::nodesFrontTopLeftCorner() const
{
  return m_nely*(m_nelx+1);
}

// --------------------------- node-number of the front-top-right corner ---------------------------

inline size_t Regular::nodesFrontTopRightCorner() const
{
  return m_nely*(m_nelx+1) + m_nelx;
}

// -------------------------- node-number of the back-bottom-left corner ---------------------------

inline size_t Regular::nodesBackBottomLeftCorner() const
{
  return m_nelz*(m_nely+1)*(m_nelx+1);
}

// -------------------------- node-number of the back-bottom-right corner --------------------------

inline size_t Regular::nodesBackBottomRightCorner() const
{
  return m_nelx + m_nelz*(m_nely+1)*(m_nelx+1);
}

// ---------------------------- node-number of the back-top-left corner ----------------------------

inline size_t Regular::nodesBackTopLeftCorner() const
{
  return m_nely*(m_nelx+1) + m_nelz*(m_nely+1)*(m_nelx+1);
}

// --------------------------- node-number of the back-top-right corner ----------------------------

inline size_t Regular::nodesBackTopRightCorner() const
{
  return m_nely*(m_nelx+1) + m_nelx + m_nelz*(m_nely+1)*(m_nelx+1);
}

// -------------------------------------------- aliases --------------------------------------------

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

// ------------------------------ node-numbers of periodic node-pairs ------------------------------

inline MatS Regular::nodesPeriodic() const
{
  // faces
  ColS fro = nodesFrontFace();
  ColS bck = nodesBackFace();
  ColS lft = nodesLeftFace();
  ColS rgt = nodesRightFace();
  ColS bot = nodesBottomFace();
  ColS top = nodesTopFace();

  // edges
  ColS froBot = nodesFrontBottomOpenEdge();
  ColS froTop = nodesFrontTopOpenEdge();
  ColS froLft = nodesFrontLeftOpenEdge();
  ColS froRgt = nodesFrontRightOpenEdge();
  ColS bckBot = nodesBackBottomOpenEdge();
  ColS bckTop = nodesBackTopOpenEdge();
  ColS bckLft = nodesBackLeftOpenEdge();
  ColS bckRgt = nodesBackRightOpenEdge();
  ColS botLft = nodesBottomLeftOpenEdge();
  ColS botRgt = nodesBottomRightOpenEdge();
  ColS topLft = nodesTopLeftOpenEdge();
  ColS topRgt = nodesTopRightOpenEdge();

  // allocate nodal ties
  // - number of tying per category
  size_t tface = fro.size() + lft.size() + bot.size();
  size_t tedge = 3*froBot.size() + 3*froLft.size() + 3*botLft.size();
  size_t tnode = 7;
  // - allocate
  MatS out(tface+tedge+tnode, 2);

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
  for ( auto j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = bckBot(j); ++i; }
  for ( auto j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = bckTop(j); ++i; }
  for ( auto j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = froTop(j); ++i; }
  for ( auto j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = botRgt(j); ++i; }
  for ( auto j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = topRgt(j); ++i; }
  for ( auto j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = topLft(j); ++i; }
  for ( auto j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = froRgt(j); ++i; }
  for ( auto j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = bckRgt(j); ++i; }
  for ( auto j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = bckLft(j); ++i; }

  // tie faces to each-other
  for ( auto j = 0 ; j<fro.size()    ; ++j ){ out(i,0) = fro(j);    out(i,1) = bck(j);    ++i; }
  for ( auto j = 0 ; j<lft.size()    ; ++j ){ out(i,0) = lft(j);    out(i,1) = rgt(j);    ++i; }
  for ( auto j = 0 ; j<bot.size()    ; ++j ){ out(i,0) = bot(j);    out(i,1) = top(j);    ++i; }

  return out;
}

// ------------------------------ node-number that lies in the origin ------------------------------

inline size_t Regular::nodesOrigin() const
{
  return nodesFrontBottomLeftCorner();
}

// ------------------------- DOF numbers per node (sequentially numbered) --------------------------

inline MatS Regular::dofs() const
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// ------------------------ DOP-numbers with periodic dependencies removed -------------------------

inline MatS Regular::dofsPeriodic() const
{
  // DOF-numbers for each component of each node (sequential)
  MatS out = GooseFEM::Mesh::dofs(m_nnode,m_ndim);

  // periodic node-pairs
  MatS   nodePer = nodesPeriodic();
  size_t nper    = static_cast<size_t>(nodePer.rows());

  // eliminate 'dependent' DOFs; renumber "out" to be sequential for the remaining DOFs
  for ( size_t i = 0 ; i < nper ; ++i )
    for ( size_t j = 0 ; j < m_ndim ; ++j )
      out(nodePer(i,1),j) = out(nodePer(i,0),j);

  // renumber "out" to be sequential
  return GooseFEM::Mesh::renumber(out);
}

// ------------------------------------------ constructor ------------------------------------------

inline FineLayer::FineLayer(size_t nelx, size_t nely, size_t nelz, double h, size_t nfine):
m_h(h), m_nelx(nelx), m_nelz(nelz)
{
  // basic assumptions
  assert( nelx >= 1 );
  assert( nely >= 1 );
  assert( nelz >= 1 );

  // store basic info
  m_Lx = m_h * static_cast<double>(nelx);
  m_Lz = m_h * static_cast<double>(nelz);

  // compute element size in y-direction (use symmetry, compute upper half)
  // ----------------------------------------------------------------------

  // temporary variables
  size_t nmin, ntot;
  ColS nhx(nely), nhy(nely), nhz(nely);
  ColI refine(nely);

  // minimum height in y-direction (half of the height because of symmetry)
  if ( nely  % 2 == 0 ) nmin  =  nely    /2;
  else                  nmin  = (nely +1)/2;

  // minimum number of fine layers in y-direction (minimum 1, middle layer part of this half)
  if ( nfine % 2 == 0 ) nfine =  nfine   /2 + 1;
  else                  nfine = (nfine+1)/2;
  if ( nfine < 1      ) nfine = 1;
  if ( nfine > nmin   ) nfine = nmin;

  // initialize to state with only fine elements
  nhx   .setOnes();
  nhy   .setOnes();
  nhz   .setOnes();
  refine.setConstant(-1);

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
      refine     (iy            )  = 0;
      nhy        (iy            ) *= 2;
      nhy.segment(iy+1,nely-iy-1) *= 3;
      nhx.segment(iy  ,nely-iy  ) *= 3;

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
          refine     (iy        )  = 2;
          nhy        (iy        )  = nhy(iy-1);
          nhz.segment(iy,nely-iy) *= 3;
        }
      }
    }

    // rules (1,2) satisfied: coarse in z-direction
    else if ( 3*nhy(iy) <= ntot and nelz%(3*nhz(iy)) == 0 and ntot+nhy(iy) < nmin )
    {
      // - process refinement in z-direction
      refine     (iy            )  = 2;
      nhy        (iy            ) *= 2;
      nhy.segment(iy+1,nely-iy-1) *= 3;
      nhz.segment(iy  ,nely-iy  ) *= 3;
    }

    // update the number of elements in y-direction
    ntot += nhy(iy);
    // proceed to next element layer in y-direction
    ++iy;
    // check to stop
    if ( iy >= nely or ntot >= nmin ) { nely = iy; break; }
  }

  // symmetrize, compute full information
  // ------------------------------------

  // allocate mesh constructor parameters
  m_nhx      .conservativeResize(nely*2-1);
  m_nhy      .conservativeResize(nely*2-1);
  m_nhz      .conservativeResize(nely*2-1);
  m_refine   .conservativeResize(nely*2-1);
  m_nelx     .conservativeResize(nely*2-1);
  m_nelz     .conservativeResize(nely*2-1);
  m_nnd      .conservativeResize(nely*2  );
  m_startElem.conservativeResize(nely*2-1);
  m_startNode.conservativeResize(nely*2  );

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
  // -----------------------

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

// -------------------------------------- number of elements ---------------------------------------

inline size_t FineLayer::nelem() const
{
  return m_nelem;
}

// ---------------------------------------- number of nodes ----------------------------------------

inline size_t FineLayer::nnode() const
{
  return m_nnode;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t FineLayer::nne() const
{
  return m_nne;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t FineLayer::ndim() const
{
  return m_ndim;
}

// ---------------------------- actual number of nodes in one direction ----------------------------

inline size_t FineLayer::shape(size_t i) const
{
  assert( i >= 0 and i <= 2 );

  if      ( i == 0 ) return m_nelx.maxCoeff();
  else if ( i == 2 ) return m_nelz.maxCoeff();
  else               return m_nhy .sum();

}

// --------------------------------- coordinates (nodal positions) ---------------------------------

inline MatD FineLayer::coor() const
{
  // allocate output
  MatD out(m_nnode, m_ndim);

  // current node, number of element layers
  size_t inode = 0;
  size_t nely  = static_cast<size_t>(m_nhy.size());

  // y-position of each main node layer (i.e. excluding node layers for refinement/coarsening)
  // - allocate
  ColD y(nely+1);
  // - initialize
  y(0) = 0.0;
  // - compute
  for ( size_t iy = 1 ; iy < nely+1 ; ++iy )
    y(iy) = y(iy-1) + m_nhy(iy-1) * m_h;

  // loop over element layers (bottom -> middle) : add bottom layer (+ refinement layer) of nodes
  // --------------------------------------------------------------------------------------------

  for ( size_t iy = 0 ; ; ++iy )
  {
    // get positions along the x- and z-axis
    ColD x = ColD::LinSpaced(m_nelx(iy)+1, 0.0, m_Lx);
    ColD z = ColD::LinSpaced(m_nelz(iy)+1, 0.0, m_Lz);

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
  // --------------------------------------------------------------------------------------

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    // get positions along the x- and z-axis
    ColD x = ColD::LinSpaced(m_nelx(iy)+1, 0.0, m_Lx);
    ColD z = ColD::LinSpaced(m_nelz(iy)+1, 0.0, m_Lz);

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

// ---------------------------- connectivity (node-numbers per element) ----------------------------

inline MatS FineLayer::conn() const
{
  // allocate output
  MatS out(m_nelem, m_nne);

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
          out(ielem,1) = bot + (ix+1) +    iz    * (m_nelx(iy)+1);
          out(ielem,2) = bot + (ix+1) + (  iz+1) * (m_nelx(iy)+1);
          out(ielem,3) = bot + (ix  ) + (  iz+1) * (m_nelx(iy)+1);
          out(ielem,4) = mid + (ix  ) +  2*iz    * (m_nelx(iy)+1);
          out(ielem,5) = mid + (ix+1) +  2*iz    * (m_nelx(iy)+1);
          out(ielem,6) = mid + (ix+1) + (2*iz+1) * (m_nelx(iy)+1);
          out(ielem,7) = mid + (ix  ) + (2*iz+1) * (m_nelx(iy)+1);
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

// ------------------------------ element numbers of the middle layer ------------------------------

inline ColS FineLayer::elementsMiddleLayer() const
{
  // number of element layers in y-direction, the index of the middle layer
  size_t nely = static_cast<size_t>(m_nhy.size());
  size_t iy   = (nely-1)/2;

  ColS out(m_nelx(iy)*m_nelz(iy));

  for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix )
    for ( size_t iz = 0 ; iz < m_nelz(iy) ; ++iz )
      out(ix+iz*m_nelx(iy)) = m_startElem(iy) + ix + iz*m_nelx(iy);

  return out;
}

// ------------------------------ node-numbers along the front plane -------------------------------

inline ColS FineLayer::nodesFront() const
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
  ColS out(n);

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

// ------------------------------- node-numbers along the back plane -------------------------------

inline ColS FineLayer::nodesBack() const
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
  ColS out(n);

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

// ------------------------------- node-numbers along the left plane -------------------------------

inline ColS FineLayer::nodesLeft() const
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
  ColS out(n);

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

// ------------------------------ node-numbers along the right plane -------------------------------

inline ColS FineLayer::nodesRight() const
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
  ColS out(n);

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

// ------------------------------ node-numbers along the bottom plane ------------------------------

inline ColS FineLayer::nodesBottom() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // allocate node list
  ColS out(m_nnd(nely));

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

// ------------------------------- node-numbers along the top plane --------------------------------

inline ColS FineLayer::nodesTop() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // allocate node list
  ColS out(m_nnd(nely));

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

// ------------------------------- node-numbers along the front face -------------------------------

inline ColS FineLayer::nodesFrontFace() const
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
  ColS out(n);

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

// ------------------------------- node-numbers along the back face --------------------------------

inline ColS FineLayer::nodesBackFace() const
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
  ColS out(n);

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

// ------------------------------- node-numbers along the left face --------------------------------

inline ColS FineLayer::nodesLeftFace() const
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
  ColS out(n);

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

// ------------------------------- node-numbers along the right face -------------------------------

inline ColS FineLayer::nodesRightFace() const
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
  ColS out(n);

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

// ------------------------------ node-numbers along the bottom face -------------------------------

inline ColS FineLayer::nodesBottomFace() const
{
  // allocate node list
  ColS out((m_nelx(0)-1)*(m_nelz(0)-1));

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

// -------------------------------- node-numbers along the top face --------------------------------

inline ColS FineLayer::nodesTopFace() const
{
  // number of element layers in y-direction
  size_t nely = static_cast<size_t>(m_nhy.size());

  // allocate node list
  ColS out((m_nelx(nely-1)-1)*(m_nelz(nely-1)-1));

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

// --------------------------- node-numbers along the front-bottom edge ----------------------------

inline ColS FineLayer::nodesFrontBottomEdge() const
{
  ColS out(m_nelx(0)+1);

  for ( size_t ix = 0 ; ix < m_nelx(0)+1 ; ++ix )
    out(ix) = m_startNode(0) + ix;

  return out;
}

// ----------------------------- node-numbers along the front-top edge -----------------------------

inline ColS FineLayer::nodesFrontTopEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(m_nelx(nely-1)+1);

  for ( size_t ix = 0 ; ix < m_nelx(nely-1)+1 ; ++ix )
    out(ix) = m_startNode(nely) + ix;

  return out;
}

// ---------------------------- node-numbers along the front-left edge -----------------------------

inline ColS FineLayer::nodesFrontLeftEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(nely+1);

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1);

  return out;
}

// ---------------------------- node-numbers along the front-right edge ----------------------------

inline ColS FineLayer::nodesFrontRightEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(nely+1);

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy) + m_nelx(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1) + m_nelx(iy);

  return out;
}

// ---------------------------- node-numbers along the back-bottom edge ----------------------------

inline ColS FineLayer::nodesBackBottomEdge() const
{
  ColS out(m_nelx(0)+1);

  for ( size_t ix = 0 ; ix < m_nelx(0)+1 ; ++ix )
    out(ix) = m_startNode(0) + ix + (m_nelx(0)+1)*(m_nelz(0));

  return out;
}

// ----------------------------- node-numbers along the back-top edge ------------------------------

inline ColS FineLayer::nodesBackTopEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(m_nelx(nely-1)+1);

  for ( size_t ix = 0 ; ix < m_nelx(nely-1)+1 ; ++ix )
    out(ix) = m_startNode(nely) + ix + (m_nelx(nely-1)+1)*(m_nelz(nely-1));

  return out;
}

// ----------------------------- node-numbers along the back-left edge -----------------------------

inline ColS FineLayer::nodesBackLeftEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(nely+1);

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1) + (m_nelx(iy)+1)*(m_nelz(iy));

  return out;
}

// ---------------------------- node-numbers along the back-right edge -----------------------------

inline ColS FineLayer::nodesBackRightEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(nely+1);

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy) + m_nelx(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1) + m_nelx(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  return out;
}

// ---------------------------- node-numbers along the bottom-left edge ----------------------------

inline ColS FineLayer::nodesBottomLeftEdge() const
{
  ColS out(m_nelz(0)+1);

  for ( size_t iz = 0 ; iz < m_nelz(0)+1 ; ++iz )
    out(iz) = m_startNode(0) + iz * (m_nelx(0)+1);

  return out;
}

// --------------------------- node-numbers along the bottom-right edge ----------------------------

inline ColS FineLayer::nodesBottomRightEdge() const
{
  ColS out(m_nelz(0)+1);

  for ( size_t iz = 0 ; iz < m_nelz(0)+1 ; ++iz )
    out(iz) = m_startNode(0) + m_nelx(0) + iz * (m_nelx(0)+1);

  return out;
}

// ----------------------------- node-numbers along the top-left edge ------------------------------

inline ColS FineLayer::nodesTopLeftEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(m_nelz(nely-1)+1);

  for ( size_t iz = 0 ; iz < m_nelz(nely-1)+1 ; ++iz )
    out(iz) = m_startNode(nely) + iz * (m_nelx(nely-1)+1);

  return out;
}

// ----------------------------- node-numbers along the top-right edge -----------------------------

inline ColS FineLayer::nodesTopRightEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(m_nelz(nely-1)+1);

  for ( size_t iz = 0 ; iz < m_nelz(nely-1)+1 ; ++iz )
    out(iz) = m_startNode(nely) + m_nelx(nely-1) + iz * (m_nelx(nely-1)+1);

  return out;
}

// -------------------------------------------- aliases --------------------------------------------

inline ColS FineLayer::nodesBottomFrontEdge() const { return nodesFrontBottomEdge(); }
inline ColS FineLayer::nodesBottomBackEdge()  const { return nodesBackBottomEdge();  }
inline ColS FineLayer::nodesTopFrontEdge()    const { return nodesFrontTopEdge();    }
inline ColS FineLayer::nodesTopBackEdge()     const { return nodesBackTopEdge();     }
inline ColS FineLayer::nodesLeftBottomEdge()  const { return nodesBottomLeftEdge();  }
inline ColS FineLayer::nodesLeftFrontEdge()   const { return nodesFrontLeftEdge();   }
inline ColS FineLayer::nodesLeftBackEdge()    const { return nodesBackLeftEdge();    }
inline ColS FineLayer::nodesLeftTopEdge()     const { return nodesTopLeftEdge();     }
inline ColS FineLayer::nodesRightBottomEdge() const { return nodesBottomRightEdge(); }
inline ColS FineLayer::nodesRightTopEdge()    const { return nodesTopRightEdge();    }
inline ColS FineLayer::nodesRightFrontEdge()  const { return nodesFrontRightEdge();  }
inline ColS FineLayer::nodesRightBackEdge()   const { return nodesBackRightEdge();   }

// ------------------- node-numbers along the front-bottom edge, without corners -------------------

inline ColS FineLayer::nodesFrontBottomOpenEdge() const
{
  ColS out(m_nelx(0)-1);

  for ( size_t ix = 1 ; ix < m_nelx(0) ; ++ix )
    out(ix-1) = m_startNode(0) + ix;

  return out;
}

// -------------------- node-numbers along the front-top edge, without corners ---------------------

inline ColS FineLayer::nodesFrontTopOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(m_nelx(nely-1)-1);

  for ( size_t ix = 1 ; ix < m_nelx(nely-1) ; ++ix )
    out(ix-1) = m_startNode(nely) + ix;

  return out;
}

// -------------------- node-numbers along the front-left edge, without corners --------------------

inline ColS FineLayer::nodesFrontLeftOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(nely-1);

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1);

  return out;
}

// ------------------- node-numbers along the front-right edge, without corners --------------------

inline ColS FineLayer::nodesFrontRightOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(nely-1);

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy) + m_nelx(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1) + m_nelx(iy);

  return out;
}

// ------------------- node-numbers along the back-bottom edge, without corners --------------------

inline ColS FineLayer::nodesBackBottomOpenEdge() const
{
  ColS out(m_nelx(0)-1);

  for ( size_t ix = 1 ; ix < m_nelx(0) ; ++ix )
    out(ix-1) = m_startNode(0) + ix + (m_nelx(0)+1)*(m_nelz(0));

  return out;
}

// --------------------- node-numbers along the back-top edge, without corners ---------------------

inline ColS FineLayer::nodesBackTopOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(m_nelx(nely-1)-1);

  for ( size_t ix = 1 ; ix < m_nelx(nely-1) ; ++ix )
    out(ix-1) = m_startNode(nely) + ix + (m_nelx(nely-1)+1)*(m_nelz(nely-1));

  return out;
}

// -------------------- node-numbers along the back-left edge, without corners ---------------------

inline ColS FineLayer::nodesBackLeftOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(nely-1);

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1) + (m_nelx(iy)+1)*(m_nelz(iy));

  return out;
}

// -------------------- node-numbers along the back-right edge, without corners --------------------

inline ColS FineLayer::nodesBackRightOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(nely-1);

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy) + m_nelx(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1) + m_nelx(iy) + (m_nelx(iy)+1)*(m_nelz(iy));

  return out;
}

// ------------------- node-numbers along the bottom-left edge, without corners --------------------

inline ColS FineLayer::nodesBottomLeftOpenEdge() const
{
  ColS out(m_nelz(0)-1);

  for ( size_t iz = 1 ; iz < m_nelz(0) ; ++iz )
    out(iz-1) = m_startNode(0) + iz * (m_nelx(0)+1);

  return out;
}

// ------------------- node-numbers along the bottom-right edge, without corners -------------------

inline ColS FineLayer::nodesBottomRightOpenEdge() const
{
  ColS out(m_nelz(0)-1);

  for ( size_t iz = 1 ; iz < m_nelz(0) ; ++iz )
    out(iz-1) = m_startNode(0) + m_nelx(0) + iz * (m_nelx(0)+1);

  return out;
}

// --------------------- node-numbers along the top-left edge, without corners ---------------------

inline ColS FineLayer::nodesTopLeftOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(m_nelz(nely-1)-1);

  for ( size_t iz = 1 ; iz < m_nelz(nely-1) ; ++iz )
    out(iz-1) = m_startNode(nely) + iz * (m_nelx(nely-1)+1);

  return out;
}

// -------------------- node-numbers along the top-right edge, without corners ---------------------

inline ColS FineLayer::nodesTopRightOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  ColS out(m_nelz(nely-1)-1);

  for ( size_t iz = 1 ; iz < m_nelz(nely-1) ; ++iz )
    out(iz-1) = m_startNode(nely) + m_nelx(nely-1) + iz * (m_nelx(nely-1)+1);

  return out;
}

// -------------------------------------------- aliases --------------------------------------------

inline ColS FineLayer::nodesBottomFrontOpenEdge() const { return nodesFrontBottomOpenEdge(); }
inline ColS FineLayer::nodesBottomBackOpenEdge() const  { return nodesBackBottomOpenEdge();  }
inline ColS FineLayer::nodesTopFrontOpenEdge() const    { return nodesFrontTopOpenEdge();    }
inline ColS FineLayer::nodesTopBackOpenEdge() const     { return nodesBackTopOpenEdge();     }
inline ColS FineLayer::nodesLeftBottomOpenEdge() const  { return nodesBottomLeftOpenEdge();  }
inline ColS FineLayer::nodesLeftFrontOpenEdge() const   { return nodesFrontLeftOpenEdge();   }
inline ColS FineLayer::nodesLeftBackOpenEdge() const    { return nodesBackLeftOpenEdge();    }
inline ColS FineLayer::nodesLeftTopOpenEdge() const     { return nodesTopLeftOpenEdge();     }
inline ColS FineLayer::nodesRightBottomOpenEdge() const { return nodesBottomRightOpenEdge(); }
inline ColS FineLayer::nodesRightTopOpenEdge() const    { return nodesTopRightOpenEdge();    }
inline ColS FineLayer::nodesRightFrontOpenEdge() const  { return nodesFrontRightOpenEdge();  }
inline ColS FineLayer::nodesRightBackOpenEdge() const   { return nodesBackRightOpenEdge();   }

// -------------------------- node-number of the front-bottom-left corner --------------------------

inline size_t FineLayer::nodesFrontBottomLeftCorner() const
{
  return m_startNode(0);
}

// ------------------------- node-number of the front-bottom-right corner --------------------------

inline size_t FineLayer::nodesFrontBottomRightCorner() const
{
  return m_startNode(0) + m_nelx(0);
}

// --------------------------- node-number of the front-top-left corner ----------------------------

inline size_t FineLayer::nodesFrontTopLeftCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely);
}

// --------------------------- node-number of the front-top-right corner ---------------------------

inline size_t FineLayer::nodesFrontTopRightCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely) + m_nelx(nely-1);
}

// -------------------------- node-number of the back-bottom-left corner ---------------------------

inline size_t FineLayer::nodesBackBottomLeftCorner() const
{
  return m_startNode(0) + (m_nelx(0)+1)*(m_nelz(0));
}

// -------------------------- node-number of the back-bottom-right corner --------------------------

inline size_t FineLayer::nodesBackBottomRightCorner() const
{
  return m_startNode(0) + m_nelx(0) + (m_nelx(0)+1)*(m_nelz(0));
}

// ---------------------------- node-number of the back-top-left corner ----------------------------

inline size_t FineLayer::nodesBackTopLeftCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely) + (m_nelx(nely-1)+1)*(m_nelz(nely-1));
}

// --------------------------- node-number of the back-top-right corner ----------------------------

inline size_t FineLayer::nodesBackTopRightCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely) + m_nelx(nely-1) + (m_nelx(nely-1)+1)*(m_nelz(nely-1));
}

// -------------------------------------------- aliases --------------------------------------------

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

// ------------------------------ node-numbers of periodic node-pairs ------------------------------

inline MatS FineLayer::nodesPeriodic() const
{
  // faces
  ColS fro = nodesFrontFace();
  ColS bck = nodesBackFace();
  ColS lft = nodesLeftFace();
  ColS rgt = nodesRightFace();
  ColS bot = nodesBottomFace();
  ColS top = nodesTopFace();

  // edges
  ColS froBot = nodesFrontBottomOpenEdge();
  ColS froTop = nodesFrontTopOpenEdge();
  ColS froLft = nodesFrontLeftOpenEdge();
  ColS froRgt = nodesFrontRightOpenEdge();
  ColS bckBot = nodesBackBottomOpenEdge();
  ColS bckTop = nodesBackTopOpenEdge();
  ColS bckLft = nodesBackLeftOpenEdge();
  ColS bckRgt = nodesBackRightOpenEdge();
  ColS botLft = nodesBottomLeftOpenEdge();
  ColS botRgt = nodesBottomRightOpenEdge();
  ColS topLft = nodesTopLeftOpenEdge();
  ColS topRgt = nodesTopRightOpenEdge();

  // allocate nodal ties
  // - number of tying per category
  size_t tface = fro.size() + lft.size() + bot.size();
  size_t tedge = 3*froBot.size() + 3*froLft.size() + 3*botLft.size();
  size_t tnode = 7;
  // - allocate
  MatS out(tface+tedge+tnode, 2);

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
  for ( auto j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = bckBot(j); ++i; }
  for ( auto j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = bckTop(j); ++i; }
  for ( auto j = 0 ; j<froBot.size() ; ++j ){ out(i,0) = froBot(j); out(i,1) = froTop(j); ++i; }
  for ( auto j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = botRgt(j); ++i; }
  for ( auto j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = topRgt(j); ++i; }
  for ( auto j = 0 ; j<botLft.size() ; ++j ){ out(i,0) = botLft(j); out(i,1) = topLft(j); ++i; }
  for ( auto j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = froRgt(j); ++i; }
  for ( auto j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = bckRgt(j); ++i; }
  for ( auto j = 0 ; j<froLft.size() ; ++j ){ out(i,0) = froLft(j); out(i,1) = bckLft(j); ++i; }

  // tie faces to each-other
  for ( auto j = 0 ; j<fro.size()    ; ++j ){ out(i,0) = fro(j);    out(i,1) = bck(j);    ++i; }
  for ( auto j = 0 ; j<lft.size()    ; ++j ){ out(i,0) = lft(j);    out(i,1) = rgt(j);    ++i; }
  for ( auto j = 0 ; j<bot.size()    ; ++j ){ out(i,0) = bot(j);    out(i,1) = top(j);    ++i; }

  return out;
}

// ------------------------------ node-number that lies in the origin ------------------------------

inline size_t FineLayer::nodesOrigin() const
{
  return nodesFrontBottomLeftCorner();
}

// ------------------------- DOF numbers per node (sequentially numbered) --------------------------

inline MatS FineLayer::dofs() const
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// ------------------------ DOP-numbers with periodic dependencies removed -------------------------

inline MatS FineLayer::dofsPeriodic() const
{
  // DOF-numbers for each component of each node (sequential)
  MatS out = GooseFEM::Mesh::dofs(m_nnode,m_ndim);

  // periodic node-pairs
  MatS   nodePer = nodesPeriodic();
  size_t nper    = static_cast<size_t>(nodePer.rows());

  // eliminate 'dependent' DOFs; renumber "out" to be sequential for the remaining DOFs
  for ( size_t i = 0 ; i < nper ; ++i )
    for ( size_t j = 0 ; j < m_ndim ; ++j )
      out(nodePer(i,1),j) = out(nodePer(i,0),j);

  // renumber "out" to be sequential
  return GooseFEM::Mesh::renumber(out);
}

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

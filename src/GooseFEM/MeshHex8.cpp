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

// ===================================== CLASS - REGULAR MESH ======================================

// ------------------------------------------ constructor ------------------------------------------

inline Regular::Regular(size_t nx, size_t ny, size_t nz, double h):
m_nx(nx), m_ny(ny), m_nz(nz), m_h(h)
{
  assert( m_nx >= 1 );
  assert( m_ny >= 1 );
  assert( m_nz >= 1 );

  m_nnode = (m_nx+1) * (m_ny+1) * (m_nz+1);
  m_nelem =  m_nx    *  m_ny    *  m_nz   ;
}

// -------------------------------------- number of elements ---------------------------------------

inline size_t Regular::nelem()
{
  return m_nelem;
}

// ---------------------------------------- number of nodes ----------------------------------------

inline size_t Regular::nnode()
{
  return m_nnode;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t Regular::nne()
{
  return m_nne;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t Regular::ndim()
{
  return m_ndim;
}

// --------------------------------- coordinates (nodal positions) ---------------------------------

inline MatD Regular::coor()
{
  MatD coor(m_nnode,m_ndim);

  ColD x = ColD::LinSpaced(m_nx+1, 0.0, m_h*static_cast<double>(m_nx));
  ColD y = ColD::LinSpaced(m_ny+1, 0.0, m_h*static_cast<double>(m_ny));
  ColD z = ColD::LinSpaced(m_nz+1, 0.0, m_h*static_cast<double>(m_nz));

  size_t inode = 0;

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz ) {
    for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy ) {
      for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix ) {
        coor(inode,0) = x(ix);
        coor(inode,1) = y(iy);
        coor(inode,2) = z(iz);
        ++inode;
      }
    }
  }

  return coor;
}

// ---------------------------- connectivity (node-numbers per element) ----------------------------

inline MatS Regular::conn()
{
  MatS conn(m_nelem,m_nne);

  size_t ielem = 0;

  for ( size_t iz = 0 ; iz < m_nz ; ++iz ) {
    for ( size_t iy = 0 ; iy < m_ny ; ++iy ) {
      for ( size_t ix = 0 ; ix < m_nx ; ++ix ) {
        conn(ielem,0) = (iz+0)*(m_ny+1)*(m_nx+1) + (iy+0)*(m_nx+1) + (ix+0);
        conn(ielem,1) = (iz+0)*(m_ny+1)*(m_nx+1) + (iy+0)*(m_nx+1) + (ix+1);
        conn(ielem,3) = (iz+0)*(m_ny+1)*(m_nx+1) + (iy+1)*(m_nx+1) + (ix+0);
        conn(ielem,2) = (iz+0)*(m_ny+1)*(m_nx+1) + (iy+1)*(m_nx+1) + (ix+1);
        conn(ielem,4) = (iz+1)*(m_ny+1)*(m_nx+1) + (iy+0)*(m_nx+1) + (ix+0);
        conn(ielem,5) = (iz+1)*(m_ny+1)*(m_nx+1) + (iy+0)*(m_nx+1) + (ix+1);
        conn(ielem,7) = (iz+1)*(m_ny+1)*(m_nx+1) + (iy+1)*(m_nx+1) + (ix+0);
        conn(ielem,6) = (iz+1)*(m_ny+1)*(m_nx+1) + (iy+1)*(m_nx+1) + (ix+1);
        ++ielem;
      }
    }
  }

  return conn;
}

// ------------------------------ node-numbers along the bottom plane ------------------------------

inline ColS Regular::nodesBottom()
{
  ColS nodes((m_nx+1)*(m_ny+1));

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
      nodes(iy*(m_nx+1)+ix) = iy*(m_nx+1) + ix;

  return nodes;
}

// ------------------------------- node-numbers along the top plane --------------------------------

inline ColS Regular::nodesTop()
{
  ColS nodes((m_nx+1)*(m_ny+1));

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
      nodes(iy*(m_nx+1)+ix) = iy*(m_nx+1) + ix + m_nz*(m_ny+1)*(m_nx+1);

  return nodes;
}

// ------------------------------- node-numbers along the left plane -------------------------------

inline ColS Regular::nodesLeft()
{
  ColS nodes((m_ny+1)*(m_nz+1));

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
      nodes(iz*(m_ny+1)+iy) = iy*(m_nx+1) + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ------------------------------ node-numbers along the right plane -------------------------------

inline ColS Regular::nodesRight()
{
  ColS nodes((m_ny+1)*(m_nz+1));

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
      nodes(iz*(m_ny+1)+iy) = iy*(m_nx+1) + iz*(m_nx+1)*(m_ny+1) + m_nx;

  return nodes;
}

// ------------------------------ node-numbers along the front plane -------------------------------

inline ColS Regular::nodesFront()
{
  ColS nodes((m_nx+1)*(m_nz+1));

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
      nodes(iz*(m_nx+1)+ix) = ix + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ------------------------------- node-numbers along the back plane -------------------------------

inline ColS Regular::nodesBack()
{
  ColS nodes((m_nx+1)*(m_nz+1));

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
      nodes(iz*(m_nx+1)+ix) = ix + m_ny*(m_nx+1) + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ------------------------------ node-numbers along the bottom face -------------------------------

inline ColS Regular::nodesBottomFace()
{
  ColS nodes((m_nx-1)*(m_ny-1));

  for ( size_t iy = 1 ; iy < m_ny ; ++iy )
    for ( size_t ix = 1 ; ix < m_nx ; ++ix )
      nodes((iy-1)*(m_nx-1)+(ix-1)) = iy*(m_nx+1) + ix;

  return nodes;
}

// -------------------------------- node-numbers along the top face --------------------------------

inline ColS Regular::nodesTopFace()
{
  ColS nodes((m_nx-1)*(m_ny-1));

  for ( size_t iy = 1 ; iy < m_ny ; ++iy ) {
    for ( size_t ix = 1 ; ix < m_nx ; ++ix ) {
      nodes((iy-1)*(m_nx-1)+(ix-1)) = iy*(m_nx+1) + ix + m_nz*(m_ny+1)*(m_nx+1);
    }
  }

  return nodes;
}

// ------------------------------- node-numbers along the left face --------------------------------

inline ColS Regular::nodesLeftFace()
{
  ColS nodes((m_ny-1)*(m_nz-1));

  for ( size_t iz = 1 ; iz < m_nz ; ++iz )
    for ( size_t iy = 1 ; iy < m_ny ; ++iy )
      nodes((iz-1)*(m_ny-1)+(iy-1)) = iy*(m_nx+1) + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ------------------------------- node-numbers along the right face -------------------------------

inline ColS Regular::nodesRightFace()
{
  ColS nodes((m_ny-1)*(m_nz-1));

  for ( size_t iz = 1 ; iz < m_nz ; ++iz )
    for ( size_t iy = 1 ; iy < m_ny ; ++iy )
      nodes((iz-1)*(m_ny-1)+(iy-1)) = iy*(m_nx+1) + iz*(m_nx+1)*(m_ny+1) + m_nx;

  return nodes;
}

// ------------------------------- node-numbers along the front face -------------------------------

inline ColS Regular::nodesFrontFace()
{
  ColS nodes((m_nx-1)*(m_nz-1));

  for ( size_t iz = 1 ; iz < m_nz ; ++iz )
    for ( size_t ix = 1 ; ix < m_nx ; ++ix )
      nodes((iz-1)*(m_nx-1)+(ix-1)) = ix + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ------------------------------- node-numbers along the back face --------------------------------

inline ColS Regular::nodesBackFace()
{
  ColS nodes((m_nx-1)*(m_nz-1));

  for ( size_t iz = 1 ; iz < m_nz ; ++iz )
    for ( size_t ix = 1 ; ix < m_nx ; ++ix )
      nodes((iz-1)*(m_nx-1)+(ix-1)) = ix + m_ny*(m_nx+1) + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// --------------------------- node-numbers along the bottom-front edge ----------------------------

inline ColS Regular::nodesBottomFrontEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = ix;

  return nodes;
}

// ---------------------------- node-numbers along the bottom-back edge ----------------------------

inline ColS Regular::nodesBottomBackEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = m_ny*(m_nx+1) + ix;

  return nodes;
}

// ---------------------------- node-numbers along the bottom-left edge ----------------------------

inline ColS Regular::nodesBottomLeftEdge()
{
  ColS nodes(m_ny+1);

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    nodes(iy) = iy*(m_nx+1);

  return nodes;
}

// --------------------------- node-numbers along the bottom-right edge ----------------------------

inline ColS Regular::nodesBottomRightEdge()
{
  ColS nodes(m_ny+1);

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    nodes(iy) = iy*(m_nx+1) + m_nx;

  return nodes;
}

// ----------------------------- node-numbers along the top-front edge -----------------------------

inline ColS Regular::nodesTopFrontEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = ix + m_nz*(m_ny+1)*(m_nx+1);

  return nodes;
}

// ----------------------------- node-numbers along the top-back edge ------------------------------

inline ColS Regular::nodesTopBackEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = m_ny*(m_nx+1) + ix + m_nz*(m_ny+1)*(m_nx+1);

  return nodes;
}

// ----------------------------- node-numbers along the top-left edge ------------------------------

inline ColS Regular::nodesTopLeftEdge()
{
  ColS nodes(m_ny+1);

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    nodes(iy) = iy*(m_nx+1) + m_nz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ----------------------------- node-numbers along the top-right edge -----------------------------

inline ColS Regular::nodesTopRightEdge()
{
  ColS nodes(m_ny+1);

    for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
      nodes(iy) = iy*(m_nx+1) + m_nz*(m_nx+1)*(m_ny+1) + m_nx;

  return nodes;
}

// ---------------------------- node-numbers along the front-left edge -----------------------------

inline ColS Regular::nodesFrontLeftEdge()
{
  ColS nodes(m_nz+1);

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    nodes(iz) = iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ---------------------------- node-numbers along the front-right edge ----------------------------

inline ColS Regular::nodesFrontRightEdge()
{
  ColS nodes(m_nz+1);

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    nodes(iz) = iz*(m_nx+1)*(m_ny+1) + m_nx;

  return nodes;
}

// -------------------------- node-numbers along the node back-left edge ---------------------------

inline ColS Regular::nodesBackLeftEdge()
{
  ColS nodes(m_nz+1);

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    nodes(iz) = m_ny*(m_nx+1) + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ---------------------------- node-numbers along the back-right edge -----------------------------

inline ColS Regular::nodesBackRightEdge()
{
  ColS nodes(m_nz+1);

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    nodes(iz) = m_ny*(m_nx+1) + iz*(m_nx+1)*(m_ny+1) + m_nx;

  return nodes;
}

// -------------------------------------------- aliases --------------------------------------------

inline ColS Regular::nodesFrontBottomEdge() { return nodesBottomFrontEdge(); }
inline ColS Regular::nodesFrontTopEdge()    { return nodesTopFrontEdge();    }
inline ColS Regular::nodesBackBottomEdge()  { return nodesBottomBackEdge();  }
inline ColS Regular::nodesBackTopEdge()     { return nodesTopBackEdge();     }
inline ColS Regular::nodesLeftFrontEdge()   { return nodesFrontLeftEdge();   }
inline ColS Regular::nodesLeftBottomEdge()  { return nodesBottomLeftEdge();  }
inline ColS Regular::nodesLeftTopEdge()     { return nodesTopLeftEdge();     }
inline ColS Regular::nodesLeftBackEdge()    { return nodesBackLeftEdge();    }
inline ColS Regular::nodesRightFrontEdge()  { return nodesFrontRightEdge();  }
inline ColS Regular::nodesRightBackEdge()   { return nodesBackRightEdge();   }
inline ColS Regular::nodesRightBottomEdge() { return nodesBottomRightEdge(); }
inline ColS Regular::nodesRightTopEdge()    { return nodesTopRightEdge();    }

// -------------------------- node-number of the bottom-front-left corner --------------------------

inline size_t Regular::nodesBottomFrontLeftCorner()
{
  return 0;
}

// -------------------------- node-number of the bottom-front-right corner --------------------------

inline size_t Regular::nodesBottomFrontRightCorner()
{
  return m_nx;
}

// -------------------------- node-number of the bottom-back-left corner ---------------------------

inline size_t Regular::nodesBottomBackLeftCorner()
{
  return m_ny*(m_nx+1);
}

// -------------------------- node-number of the bottom-back-right corner --------------------------

inline size_t Regular::nodesBottomBackRightCorner()
{
  return m_ny*(m_nx+1) + m_nx;
}

// -------------------------- node-number of the top-front-left corner --------------------------

inline size_t Regular::nodesTopFrontLeftCorner()
{
  return m_nz*(m_ny+1)*(m_nx+1);
}

// -------------------------- node-number of the top-front-right corner --------------------------

inline size_t Regular::nodesTopFrontRightCorner()
{
  return m_nx + m_nz*(m_ny+1)*(m_nx+1);
}

// -------------------------- node-number of the top-back-left corner ---------------------------

inline size_t Regular::nodesTopBackLeftCorner()
{
  return m_ny*(m_nx+1) + m_nz*(m_ny+1)*(m_nx+1);
}

// -------------------------- node-number of the top-back-right corner --------------------------

inline size_t Regular::nodesTopBackRightCorner()
{
  return m_ny*(m_nx+1) + m_nx + m_nz*(m_ny+1)*(m_nx+1);
}

// -------------------------------------------- aliases --------------------------------------------

inline size_t Regular::nodesBottomLeftFrontCorner()  { return nodesBottomFrontLeftCorner();  }
inline size_t Regular::nodesFrontBottomLeftCorner()  { return nodesBottomFrontLeftCorner();  }
inline size_t Regular::nodesFrontLeftBottomCorner()  { return nodesBottomFrontLeftCorner();  }
inline size_t Regular::nodesLeftBottomFrontCorner()  { return nodesBottomFrontLeftCorner();  }
inline size_t Regular::nodesLeftFrontBottomCorner()  { return nodesBottomFrontLeftCorner();  }
inline size_t Regular::nodesBottomRightFrontCorner() { return nodesBottomFrontRightCorner(); }
inline size_t Regular::nodesFrontBottomRightCorner() { return nodesBottomFrontRightCorner(); }
inline size_t Regular::nodesFrontRightBottomCorner() { return nodesBottomFrontRightCorner(); }
inline size_t Regular::nodesRightBottomFrontCorner() { return nodesBottomFrontRightCorner(); }
inline size_t Regular::nodesRightFrontBottomCorner() { return nodesBottomFrontRightCorner(); }
inline size_t Regular::nodesBottomLeftBackCorner()   { return nodesBottomBackLeftCorner();   }
inline size_t Regular::nodesBackBottomLeftCorner()   { return nodesBottomBackLeftCorner();   }
inline size_t Regular::nodesBackLeftBottomCorner()   { return nodesBottomBackLeftCorner();   }
inline size_t Regular::nodesLeftBottomBackCorner()   { return nodesBottomBackLeftCorner();   }
inline size_t Regular::nodesLeftBackBottomCorner()   { return nodesBottomBackLeftCorner();   }
inline size_t Regular::nodesBottomRightBackCorner()  { return nodesBottomBackRightCorner();  }
inline size_t Regular::nodesBackBottomRightCorner()  { return nodesBottomBackRightCorner();  }
inline size_t Regular::nodesBackRightBottomCorner()  { return nodesBottomBackRightCorner();  }
inline size_t Regular::nodesRightBottomBackCorner()  { return nodesBottomBackRightCorner();  }
inline size_t Regular::nodesRightBackBottomCorner()  { return nodesBottomBackRightCorner();  }
inline size_t Regular::nodesTopLeftFrontCorner()     { return nodesTopFrontLeftCorner();     }
inline size_t Regular::nodesFrontTopLeftCorner()     { return nodesTopFrontLeftCorner();     }
inline size_t Regular::nodesFrontLeftTopCorner()     { return nodesTopFrontLeftCorner();     }
inline size_t Regular::nodesLeftTopFrontCorner()     { return nodesTopFrontLeftCorner();     }
inline size_t Regular::nodesLeftFrontTopCorner()     { return nodesTopFrontLeftCorner();     }
inline size_t Regular::nodesTopRightFrontCorner()    { return nodesTopFrontRightCorner();    }
inline size_t Regular::nodesFrontTopRightCorner()    { return nodesTopFrontRightCorner();    }
inline size_t Regular::nodesFrontRightTopCorner()    { return nodesTopFrontRightCorner();    }
inline size_t Regular::nodesRightTopFrontCorner()    { return nodesTopFrontRightCorner();    }
inline size_t Regular::nodesRightFrontTopCorner()    { return nodesTopFrontRightCorner();    }
inline size_t Regular::nodesTopLeftBackCorner()      { return nodesTopBackLeftCorner();      }
inline size_t Regular::nodesBackTopLeftCorner()      { return nodesTopBackLeftCorner();      }
inline size_t Regular::nodesBackLeftTopCorner()      { return nodesTopBackLeftCorner();      }
inline size_t Regular::nodesLeftTopBackCorner()      { return nodesTopBackLeftCorner();      }
inline size_t Regular::nodesLeftBackTopCorner()      { return nodesTopBackLeftCorner();      }
inline size_t Regular::nodesTopRightBackCorner()     { return nodesTopBackRightCorner();     }
inline size_t Regular::nodesBackTopRightCorner()     { return nodesTopBackRightCorner();     }
inline size_t Regular::nodesBackRightTopCorner()     { return nodesTopBackRightCorner();     }
inline size_t Regular::nodesRightTopBackCorner()     { return nodesTopBackRightCorner();     }
inline size_t Regular::nodesRightBackTopCorner()     { return nodesTopBackRightCorner();     }

// ------------------------------ node-numbers of periodic node-pairs ------------------------------

inline MatS Regular::nodesPeriodic()
{
  // faces
  ColS bot = nodesBottomFace();
  ColS top = nodesTopFace();
  ColS lft = nodesLeftFace();
  ColS rgt = nodesRightFace();
  ColS fro = nodesFrontFace();
  ColS bck = nodesBackFace();

  // edges
  ColS botFro = nodesBottomFrontEdge();
  ColS botBck = nodesBottomBackEdge();
  ColS botLft = nodesBottomLeftEdge();
  ColS botRgt = nodesBottomRightEdge();
  ColS topFro = nodesTopFrontEdge();
  ColS topBck = nodesTopBackEdge();
  ColS topLft = nodesTopLeftEdge();
  ColS topRgt = nodesTopRightEdge();
  ColS froLft = nodesFrontLeftEdge();
  ColS froRgt = nodesFrontRightEdge();
  ColS bckLft = nodesBackLeftEdge();
  ColS bckRgt = nodesBackRightEdge();

  // allocate nodal ties
  // - number of tying per category
  size_t tface = bot.size() + lft.size() + fro.size();
  size_t tedge = 3*(botLft.size()-2) + 3*(botFro.size()-2) + 3*(froLft.size()-2);
  size_t tnode = 7;
  // - allocate
  MatS nodes(tface+tedge+tnode, 2);

  // counter
  size_t i = 0;

  // tie all corners to one corner
  nodes(i,0) = nodesBottomFrontLeftCorner(); nodes(i,1) = nodesBottomFrontRightCorner(); ++i;
  nodes(i,0) = nodesBottomFrontLeftCorner(); nodes(i,1) = nodesBottomBackRightCorner (); ++i;
  nodes(i,0) = nodesBottomFrontLeftCorner(); nodes(i,1) = nodesBottomBackLeftCorner  (); ++i;
  nodes(i,0) = nodesBottomFrontLeftCorner(); nodes(i,1) = nodesTopFrontLeftCorner    (); ++i;
  nodes(i,0) = nodesBottomFrontLeftCorner(); nodes(i,1) = nodesTopFrontRightCorner   (); ++i;
  nodes(i,0) = nodesBottomFrontLeftCorner(); nodes(i,1) = nodesTopBackRightCorner    (); ++i;
  nodes(i,0) = nodesBottomFrontLeftCorner(); nodes(i,1) = nodesTopBackLeftCorner     (); ++i;

  // tie all corresponding edges to each other (exclude corners)
  for ( auto j=1; j<botFro.size()-1; ++j ) { nodes(i,0) = botFro(j); nodes(i,1) = botBck(j); ++i; }
  for ( auto j=1; j<botFro.size()-1; ++j ) { nodes(i,0) = botFro(j); nodes(i,1) = topBck(j); ++i; }
  for ( auto j=1; j<botFro.size()-1; ++j ) { nodes(i,0) = botFro(j); nodes(i,1) = topFro(j); ++i; }
  for ( auto j=1; j<botLft.size()-1; ++j ) { nodes(i,0) = botLft(j); nodes(i,1) = botRgt(j); ++i; }
  for ( auto j=1; j<botLft.size()-1; ++j ) { nodes(i,0) = botLft(j); nodes(i,1) = topRgt(j); ++i; }
  for ( auto j=1; j<botLft.size()-1; ++j ) { nodes(i,0) = botLft(j); nodes(i,1) = topLft(j); ++i; }
  for ( auto j=1; j<froLft.size()-1; ++j ) { nodes(i,0) = froLft(j); nodes(i,1) = froRgt(j); ++i; }
  for ( auto j=1; j<froLft.size()-1; ++j ) { nodes(i,0) = froLft(j); nodes(i,1) = bckRgt(j); ++i; }
  for ( auto j=1; j<froLft.size()-1; ++j ) { nodes(i,0) = froLft(j); nodes(i,1) = bckLft(j); ++i; }

  // tie faces to each-other
  for ( auto j = 0 ; j < bot.size() ; ++j ) { nodes(i,0) = bot(j); nodes(i,1) = top(j); ++i; }
  for ( auto j = 0 ; j < lft.size() ; ++j ) { nodes(i,0) = lft(j); nodes(i,1) = rgt(j); ++i; }
  for ( auto j = 0 ; j < fro.size() ; ++j ) { nodes(i,0) = fro(j); nodes(i,1) = bck(j); ++i; }

  return nodes;
}

// ------------------------------ node-number that lies in the origin ------------------------------

inline size_t Regular::nodesOrigin()
{
  return nodesBottomFrontLeftCorner();
}

// ------------------------- DOF numbers per node (sequentially numbered) --------------------------

inline MatS Regular::dofs()
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// ------------------------ DOP-numbers with periodic dependencies removed -------------------------

inline MatS Regular::dofsPeriodic()
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

// =================================================================================================

}}} // namespace ...

// =================================================================================================

#endif

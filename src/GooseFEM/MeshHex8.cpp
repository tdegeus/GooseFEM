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

// ------------------------------ node-numbers along the front plane ------------------------------

inline ColS Regular::nodesFront()
{
  ColS nodes((m_nx+1)*(m_ny+1));

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
      nodes(iy*(m_nx+1)+ix) = iy*(m_nx+1) + ix;

  return nodes;
}

// ------------------------------- node-numbers along the back plane --------------------------------

inline ColS Regular::nodesBack()
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

// ------------------------------ node-numbers along the bottom plane -------------------------------

inline ColS Regular::nodesBottom()
{
  ColS nodes((m_nx+1)*(m_nz+1));

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
      nodes(iz*(m_nx+1)+ix) = ix + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ------------------------------- node-numbers along the top plane -------------------------------

inline ColS Regular::nodesTop()
{
  ColS nodes((m_nx+1)*(m_nz+1));

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
      nodes(iz*(m_nx+1)+ix) = ix + m_ny*(m_nx+1) + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ------------------------------ node-numbers along the front face -------------------------------

inline ColS Regular::nodesFrontFace()
{
  ColS nodes((m_nx-1)*(m_ny-1));

  for ( size_t iy = 1 ; iy < m_ny ; ++iy )
    for ( size_t ix = 1 ; ix < m_nx ; ++ix )
      nodes((iy-1)*(m_nx-1)+(ix-1)) = iy*(m_nx+1) + ix;

  return nodes;
}

// -------------------------------- node-numbers along the back face --------------------------------

inline ColS Regular::nodesBackFace()
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

// ------------------------------- node-numbers along the bottom face -------------------------------

inline ColS Regular::nodesBottomFace()
{
  ColS nodes((m_nx-1)*(m_nz-1));

  for ( size_t iz = 1 ; iz < m_nz ; ++iz )
    for ( size_t ix = 1 ; ix < m_nx ; ++ix )
      nodes((iz-1)*(m_nx-1)+(ix-1)) = ix + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ------------------------------- node-numbers along the top face --------------------------------

inline ColS Regular::nodesTopFace()
{
  ColS nodes((m_nx-1)*(m_nz-1));

  for ( size_t iz = 1 ; iz < m_nz ; ++iz )
    for ( size_t ix = 1 ; ix < m_nx ; ++ix )
      nodes((iz-1)*(m_nx-1)+(ix-1)) = ix + m_ny*(m_nx+1) + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// --------------------------- node-numbers along the front-bottom edge ----------------------------

inline ColS Regular::nodesFrontBottomEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = ix;

  return nodes;
}

// ---------------------------- node-numbers along the front-top edge ----------------------------

inline ColS Regular::nodesFrontTopEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = m_ny*(m_nx+1) + ix;

  return nodes;
}

// ---------------------------- node-numbers along the front-left edge ----------------------------

inline ColS Regular::nodesFrontLeftEdge()
{
  ColS nodes(m_ny+1);

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    nodes(iy) = iy*(m_nx+1);

  return nodes;
}

// --------------------------- node-numbers along the front-right edge ----------------------------

inline ColS Regular::nodesFrontRightEdge()
{
  ColS nodes(m_ny+1);

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    nodes(iy) = iy*(m_nx+1) + m_nx;

  return nodes;
}

// ----------------------------- node-numbers along the back-bottom edge -----------------------------

inline ColS Regular::nodesBackBottomEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = ix + m_nz*(m_ny+1)*(m_nx+1);

  return nodes;
}

// ----------------------------- node-numbers along the back-top edge ------------------------------

inline ColS Regular::nodesBackTopEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = m_ny*(m_nx+1) + ix + m_nz*(m_ny+1)*(m_nx+1);

  return nodes;
}

// ----------------------------- node-numbers along the back-left edge ------------------------------

inline ColS Regular::nodesBackLeftEdge()
{
  ColS nodes(m_ny+1);

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    nodes(iy) = iy*(m_nx+1) + m_nz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ----------------------------- node-numbers along the back-right edge -----------------------------

inline ColS Regular::nodesBackRightEdge()
{
  ColS nodes(m_ny+1);

    for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
      nodes(iy) = iy*(m_nx+1) + m_nz*(m_nx+1)*(m_ny+1) + m_nx;

  return nodes;
}

// ---------------------------- node-numbers along the bottom-left edge -----------------------------

inline ColS Regular::nodesBottomLeftEdge()
{
  ColS nodes(m_nz+1);

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    nodes(iz) = iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ---------------------------- node-numbers along the bottom-right edge ----------------------------

inline ColS Regular::nodesBottomRightEdge()
{
  ColS nodes(m_nz+1);

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    nodes(iz) = iz*(m_nx+1)*(m_ny+1) + m_nx;

  return nodes;
}

// -------------------------- node-numbers along the node top-left edge ---------------------------

inline ColS Regular::nodesTopLeftEdge()
{
  ColS nodes(m_nz+1);

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    nodes(iz) = m_ny*(m_nx+1) + iz*(m_nx+1)*(m_ny+1);

  return nodes;
}

// ---------------------------- node-numbers along the top-right edge -----------------------------

inline ColS Regular::nodesTopRightEdge()
{
  ColS nodes(m_nz+1);

  for ( size_t iz = 0 ; iz < m_nz+1 ; ++iz )
    nodes(iz) = m_ny*(m_nx+1) + iz*(m_nx+1)*(m_ny+1) + m_nx;

  return nodes;
}

// -------------------------------------------- aliases --------------------------------------------

inline ColS Regular::nodesBottomFrontEdge() { return nodesFrontBottomEdge(); }
inline ColS Regular::nodesBottomBackEdge()    { return nodesBackBottomEdge();    }
inline ColS Regular::nodesTopFrontEdge()  { return nodesFrontTopEdge();  }
inline ColS Regular::nodesTopBackEdge()     { return nodesBackTopEdge();     }
inline ColS Regular::nodesLeftBottomEdge()   { return nodesBottomLeftEdge();   }
inline ColS Regular::nodesLeftFrontEdge()  { return nodesFrontLeftEdge();  }
inline ColS Regular::nodesLeftBackEdge()     { return nodesBackLeftEdge();     }
inline ColS Regular::nodesLeftTopEdge()    { return nodesTopLeftEdge();    }
inline ColS Regular::nodesRightBottomEdge()  { return nodesBottomRightEdge();  }
inline ColS Regular::nodesRightTopEdge()   { return nodesTopRightEdge();   }
inline ColS Regular::nodesRightFrontEdge() { return nodesFrontRightEdge(); }
inline ColS Regular::nodesRightBackEdge()    { return nodesBackRightEdge();    }

// -------------------------- node-number of the front-bottom-left corner --------------------------

inline size_t Regular::nodesFrontBottomLeftCorner()
{
  return 0;
}

// -------------------------- node-number of the front-bottom-right corner --------------------------

inline size_t Regular::nodesFrontBottomRightCorner()
{
  return m_nx;
}

// -------------------------- node-number of the front-top-left corner ---------------------------

inline size_t Regular::nodesFrontTopLeftCorner()
{
  return m_ny*(m_nx+1);
}

// -------------------------- node-number of the front-top-right corner --------------------------

inline size_t Regular::nodesFrontTopRightCorner()
{
  return m_ny*(m_nx+1) + m_nx;
}

// -------------------------- node-number of the back-bottom-left corner --------------------------

inline size_t Regular::nodesBackBottomLeftCorner()
{
  return m_nz*(m_ny+1)*(m_nx+1);
}

// -------------------------- node-number of the back-bottom-right corner --------------------------

inline size_t Regular::nodesBackBottomRightCorner()
{
  return m_nx + m_nz*(m_ny+1)*(m_nx+1);
}

// -------------------------- node-number of the back-top-left corner ---------------------------

inline size_t Regular::nodesBackTopLeftCorner()
{
  return m_ny*(m_nx+1) + m_nz*(m_ny+1)*(m_nx+1);
}

// -------------------------- node-number of the back-top-right corner --------------------------

inline size_t Regular::nodesBackTopRightCorner()
{
  return m_ny*(m_nx+1) + m_nx + m_nz*(m_ny+1)*(m_nx+1);
}

// -------------------------------------------- aliases --------------------------------------------

inline size_t Regular::nodesFrontLeftBottomCorner()  { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesBottomFrontLeftCorner()  { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesBottomLeftFrontCorner()  { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesLeftFrontBottomCorner()  { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesLeftBottomFrontCorner()  { return nodesFrontBottomLeftCorner();  }
inline size_t Regular::nodesFrontRightBottomCorner() { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesBottomFrontRightCorner() { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesBottomRightFrontCorner() { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesRightFrontBottomCorner() { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesRightBottomFrontCorner() { return nodesFrontBottomRightCorner(); }
inline size_t Regular::nodesFrontLeftTopCorner()   { return nodesFrontTopLeftCorner();   }
inline size_t Regular::nodesTopFrontLeftCorner()   { return nodesFrontTopLeftCorner();   }
inline size_t Regular::nodesTopLeftFrontCorner()   { return nodesFrontTopLeftCorner();   }
inline size_t Regular::nodesLeftFrontTopCorner()   { return nodesFrontTopLeftCorner();   }
inline size_t Regular::nodesLeftTopFrontCorner()   { return nodesFrontTopLeftCorner();   }
inline size_t Regular::nodesFrontRightTopCorner()  { return nodesFrontTopRightCorner();  }
inline size_t Regular::nodesTopFrontRightCorner()  { return nodesFrontTopRightCorner();  }
inline size_t Regular::nodesTopRightFrontCorner()  { return nodesFrontTopRightCorner();  }
inline size_t Regular::nodesRightFrontTopCorner()  { return nodesFrontTopRightCorner();  }
inline size_t Regular::nodesRightTopFrontCorner()  { return nodesFrontTopRightCorner();  }
inline size_t Regular::nodesBackLeftBottomCorner()     { return nodesBackBottomLeftCorner();     }
inline size_t Regular::nodesBottomBackLeftCorner()     { return nodesBackBottomLeftCorner();     }
inline size_t Regular::nodesBottomLeftBackCorner()     { return nodesBackBottomLeftCorner();     }
inline size_t Regular::nodesLeftBackBottomCorner()     { return nodesBackBottomLeftCorner();     }
inline size_t Regular::nodesLeftBottomBackCorner()     { return nodesBackBottomLeftCorner();     }
inline size_t Regular::nodesBackRightBottomCorner()    { return nodesBackBottomRightCorner();    }
inline size_t Regular::nodesBottomBackRightCorner()    { return nodesBackBottomRightCorner();    }
inline size_t Regular::nodesBottomRightBackCorner()    { return nodesBackBottomRightCorner();    }
inline size_t Regular::nodesRightBackBottomCorner()    { return nodesBackBottomRightCorner();    }
inline size_t Regular::nodesRightBottomBackCorner()    { return nodesBackBottomRightCorner();    }
inline size_t Regular::nodesBackLeftTopCorner()      { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesTopBackLeftCorner()      { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesTopLeftBackCorner()      { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesLeftBackTopCorner()      { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesLeftTopBackCorner()      { return nodesBackTopLeftCorner();      }
inline size_t Regular::nodesBackRightTopCorner()     { return nodesBackTopRightCorner();     }
inline size_t Regular::nodesTopBackRightCorner()     { return nodesBackTopRightCorner();     }
inline size_t Regular::nodesTopRightBackCorner()     { return nodesBackTopRightCorner();     }
inline size_t Regular::nodesRightBackTopCorner()     { return nodesBackTopRightCorner();     }
inline size_t Regular::nodesRightTopBackCorner()     { return nodesBackTopRightCorner();     }

// ------------------------------ node-numbers of periodic node-pairs ------------------------------

inline MatS Regular::nodesPeriodic()
{
  // faces
  ColS fro = nodesFrontFace();
  ColS bck = nodesBackFace();
  ColS lft = nodesLeftFace();
  ColS rgt = nodesRightFace();
  ColS bot = nodesBottomFace();
  ColS top = nodesTopFace();

  // edges
  ColS froBot = nodesFrontBottomEdge();
  ColS froTop = nodesFrontTopEdge();
  ColS froLft = nodesFrontLeftEdge();
  ColS froRgt = nodesFrontRightEdge();
  ColS bckBot = nodesBackBottomEdge();
  ColS bckTop = nodesBackTopEdge();
  ColS bckLft = nodesBackLeftEdge();
  ColS bckRgt = nodesBackRightEdge();
  ColS botLft = nodesBottomLeftEdge();
  ColS botRgt = nodesBottomRightEdge();
  ColS topLft = nodesTopLeftEdge();
  ColS topRgt = nodesTopRightEdge();

  // allocate nodal ties
  // - number of tying per category
  size_t tface = fro.size() + lft.size() + bot.size();
  size_t tedge = 3*(froBot.size()-2) + 3*(froLft.size()-2) + 3*(botLft.size()-2);
  size_t tnode = 7;
  // - allocate
  MatS nodes(tface+tedge+tnode, 2);

  // counter
  size_t i = 0;

  // tie all corners to one corner
  nodes(i,0) = nodesFrontBottomLeftCorner(); nodes(i,1) = nodesFrontBottomRightCorner(); ++i;
  nodes(i,0) = nodesFrontBottomLeftCorner(); nodes(i,1) = nodesBackBottomRightCorner();  ++i;
  nodes(i,0) = nodesFrontBottomLeftCorner(); nodes(i,1) = nodesBackBottomLeftCorner();   ++i;
  nodes(i,0) = nodesFrontBottomLeftCorner(); nodes(i,1) = nodesFrontTopLeftCorner();     ++i;
  nodes(i,0) = nodesFrontBottomLeftCorner(); nodes(i,1) = nodesFrontTopRightCorner();    ++i;
  nodes(i,0) = nodesFrontBottomLeftCorner(); nodes(i,1) = nodesBackTopRightCorner();     ++i;
  nodes(i,0) = nodesFrontBottomLeftCorner(); nodes(i,1) = nodesBackTopLeftCorner();      ++i;

  // tie all corresponding edges to each other (exclude corners)
  for ( auto j=1; j<froBot.size()-1; ++j ) { nodes(i,0) = froBot(j); nodes(i,1) = bckBot(j); ++i; }
  for ( auto j=1; j<froBot.size()-1; ++j ) { nodes(i,0) = froBot(j); nodes(i,1) = bckTop(j); ++i; }
  for ( auto j=1; j<froBot.size()-1; ++j ) { nodes(i,0) = froBot(j); nodes(i,1) = froTop(j); ++i; }
  for ( auto j=1; j<botLft.size()-1; ++j ) { nodes(i,0) = botLft(j); nodes(i,1) = botRgt(j); ++i; }
  for ( auto j=1; j<botLft.size()-1; ++j ) { nodes(i,0) = botLft(j); nodes(i,1) = topRgt(j); ++i; }
  for ( auto j=1; j<botLft.size()-1; ++j ) { nodes(i,0) = botLft(j); nodes(i,1) = topLft(j); ++i; }
  for ( auto j=1; j<froLft.size()-1; ++j ) { nodes(i,0) = froLft(j); nodes(i,1) = froRgt(j); ++i; }
  for ( auto j=1; j<froLft.size()-1; ++j ) { nodes(i,0) = froLft(j); nodes(i,1) = bckRgt(j); ++i; }
  for ( auto j=1; j<froLft.size()-1; ++j ) { nodes(i,0) = froLft(j); nodes(i,1) = bckLft(j); ++i; }

  // tie faces to each-other
  for ( auto j = 0 ; j < fro.size() ; ++j ) { nodes(i,0) = fro(j); nodes(i,1) = bck(j); ++i; }
  for ( auto j = 0 ; j < lft.size() ; ++j ) { nodes(i,0) = lft(j); nodes(i,1) = rgt(j); ++i; }
  for ( auto j = 0 ; j < bot.size() ; ++j ) { nodes(i,0) = bot(j); nodes(i,1) = top(j); ++i; }

  return nodes;
}

// ------------------------------ node-number that lies in the origin ------------------------------

inline size_t Regular::nodesOrigin()
{
  return nodesFrontBottomLeftCorner();
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

// ==================================== CLASS - FINELAYER MESH =====================================

// ------------------------------------------ constructor ------------------------------------------

inline FineLayer::FineLayer(size_t nx, size_t ny, size_t nz, double h, size_t nfine, size_t nskip):
m_h(h), m_nx(nx), m_nz(nz)
{
  assert( nx >= 1 );
  assert( ny >= 1 );
  assert( nz >= 1 );

  // TODO: fake data
  size_t N = 9;

  // allocate counters
  m_nx       .conservativeResize(N  );
  m_nz       .conservativeResize(N  );
  m_nhx      .conservativeResize(N  );
  m_nhy      .conservativeResize(N  );
  m_nhz      .conservativeResize(N  );
  m_refine   .conservativeResize(N  );
  m_startElem.conservativeResize(N  );
  m_startNode.conservativeResize(N+1);

  // TODO: fake data
  m_Lx = m_h * 6;
  m_Lz = m_h * 6;

  // TODO: fake data
  m_nhx(0) = 3;  m_nhz(0) = 3;  m_nhy(0) = 3;  m_refine(0) = -1; m_nx(0) = 2; m_nz(0) = 2;
  m_nhx(1) = 3;  m_nhz(1) = 3;  m_nhy(1) = 2;  m_refine(1) =  2; m_nx(1) = 2; m_nz(1) = 2;
  m_nhx(2) = 3;  m_nhz(2) = 1;  m_nhy(2) = 2;  m_refine(2) =  0; m_nx(2) = 2; m_nz(2) = 6;
  m_nhx(3) = 1;  m_nhz(3) = 1;  m_nhy(3) = 1;  m_refine(3) = -1; m_nx(3) = 6; m_nz(3) = 6;
  m_nhx(4) = 1;  m_nhz(4) = 1;  m_nhy(4) = 1;  m_refine(4) = -1; m_nx(4) = 6; m_nz(4) = 6;
  m_nhx(5) = 1;  m_nhz(5) = 1;  m_nhy(5) = 1;  m_refine(5) = -1; m_nx(5) = 6; m_nz(5) = 6;
  m_nhx(6) = 3;  m_nhz(6) = 1;  m_nhy(6) = 2;  m_refine(6) =  0; m_nx(6) = 2; m_nz(6) = 6;
  m_nhx(7) = 3;  m_nhz(7) = 3;  m_nhy(7) = 2;  m_refine(7) =  2; m_nx(7) = 2; m_nz(7) = 2;
  m_nhx(8) = 3;  m_nhz(8) = 3;  m_nhy(8) = 3;  m_refine(8) = -1; m_nx(8) = 2; m_nz(8) = 2;

  // compute mesh dimensions
  // -----------------------

  // initialize
  m_nnode        = 0;
  m_nelem        = 0;
  m_startNode(0) = 0;

  // loop from bottom to middle : elements become finer
  for ( size_t i = 0 ; i < (N-1)/2 ; ++i )
  {
    // - store the first element of the layer
    m_startElem(i) = m_nelem;
    // - add the nodes of this layer
    if      ( m_refine(i) == 0 ) { m_nnode += (3*m_nx(i)+1) * (  m_nz(i)+1); }
    else if ( m_refine(i) == 2 ) { m_nnode += (  m_nx(i)+1) * (3*m_nz(i)+1); }
    else                         { m_nnode += (  m_nx(i)+1) * (  m_nz(i)+1); }
    // - add the elements of this layer
    if      ( m_refine(i) == 0 ) { m_nelem += 4*m_nx(i) *   m_nz(i); }
    else if ( m_refine(i) == 2 ) { m_nelem +=   m_nx(i) * 4*m_nz(i); }
    else                         { m_nelem +=   m_nx(i) *   m_nz(i); }
    // - store the starting node of the next layer
    m_startNode(i+1) = m_nnode;
  }

  // loop from middle to top : elements become coarser
  for ( size_t i = (N-1)/2 ; i < N ; ++i )
  {
    // - store the first element of the layer
    m_startElem(i) = m_nelem;
    // - add the nodes of this layer
    if      ( m_refine(i) == 0 ) { m_nnode += (5*m_nx(i)+1) * (  m_nz(i)+1); }
    else if ( m_refine(i) == 2 ) { m_nnode += (  m_nx(i)+1) * (5*m_nz(i)+1); }
    else                         { m_nnode += (  m_nx(i)+1) * (  m_nz(i)+1); }
    // - add the elements of this layer
    if      ( m_refine(i) == 0 ) { m_nelem += 4*m_nx(i) *   m_nz(i); }
    else if ( m_refine(i) == 2 ) { m_nelem +=   m_nx(i) * 4*m_nz(i); }
    else                         { m_nelem +=   m_nx(i) *   m_nz(i); }
    // - store the starting node of the next layer
    m_startNode(i+1) = m_nnode;
  }
  // - add the top row of nodes
  m_nnode += (m_nx(N-1)+1) * (m_nz(N-1)+1);
}

// -------------------------------------- number of elements ---------------------------------------

inline size_t FineLayer::nelem()
{
  return m_nelem;
}

// ---------------------------------------- number of nodes ----------------------------------------

inline size_t FineLayer::nnode()
{
  return m_nnode;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t FineLayer::nne()
{
  return m_nne;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t FineLayer::ndim()
{
  return m_ndim;
}

// --------------------------------- coordinates (nodal positions) ---------------------------------

inline MatD FineLayer::coor()
{
  // allocate output
  MatD out(m_nnode, m_ndim);

  // current node, number of element layers
  size_t inode = 0;
  size_t N     = static_cast<size_t>(m_nhy.size());

  // y-position of each main node layer (i.e. excluding node layers for refinement/coarsening)
  // - allocate
  ColD y(N+1);
  // - initialize
  y(0) = 0.0;
  // - compute
  for ( size_t iy = 1 ; iy < N+1 ; ++iy )
    y(iy) = y(iy-1) + m_nhy(iy-1) * m_h;

  // loop over element layers (bottom -> middle) : add bottom layer (+ refinement layer) of nodes
  // --------------------------------------------------------------------------------------------

  for ( size_t iy = 0 ; ; ++iy )
  {
    // get positions along the x- and z-axis
    ColD x = ColD::LinSpaced(m_nx(iy)+1, 0.0, m_Lx);
    ColD z = ColD::LinSpaced(m_nz(iy)+1, 0.0, m_Lz);

    // add nodes of the bottom layer of this element
    for ( size_t iz = 0 ; iz < m_nz(iy)+1 ; ++iz ) {
      for ( size_t ix = 0 ; ix < m_nx(iy)+1 ; ++ix ) {
        out(inode,0) = x(ix);
        out(inode,1) = y(iy);
        out(inode,2) = z(iz);
        ++inode;
      }
    }

    // stop at middle layer
    if ( iy == (N-1)/2 )
      break;

    // add extra nodes of the intermediate layer, for refinement in x-direction
    if ( m_refine(iy) == 0 )
    {
      // - get position offset in x- and y-direction
      double dx = m_h * static_cast<double>(m_nhx(iy)/3);
      double dy = m_h * static_cast<double>(m_nhy(iy)/2);
      // - add nodes of the intermediate layer
      for ( size_t iz = 0 ; iz < m_nz(iy)+1 ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nx(iy) ; ++ix ) {
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
      for ( size_t iz = 0 ; iz < m_nz(iy) ; ++iz ) {
        for ( size_t j = 0 ; j < 2 ; ++j ) {
          for ( size_t ix = 0 ; ix < m_nx(iy)+1 ; ++ix ) {
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

  for ( size_t iy = (N-1)/2 ; iy < N ; ++iy )
  {
    // get positions along the x- and z-axis
    ColD x = ColD::LinSpaced(m_nx(iy)+1, 0.0, m_Lx);
    ColD z = ColD::LinSpaced(m_nz(iy)+1, 0.0, m_Lz);

    // add extra nodes of the intermediate layer, for refinement in x-direction
    if ( m_refine(iy) == 0 )
    {
      // - get position offset in x- and y-direction
      double dx = m_h * static_cast<double>(m_nhx(iy)/3);
      double dy = m_h * static_cast<double>(m_nhy(iy)/2);
      // - add nodes of the intermediate layer
      for ( size_t iz = 0 ; iz < m_nz(iy)+1 ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nx(iy) ; ++ix ) {
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
      for ( size_t iz = 0 ; iz < m_nz(iy) ; ++iz ) {
        for ( size_t j = 0 ; j < 2 ; ++j ) {
          for ( size_t ix = 0 ; ix < m_nx(iy)+1 ; ++ix ) {
            out(inode,0) = x(ix);
            out(inode,1) = y(iy) + dy;
            out(inode,2) = z(iz) + dz * static_cast<double>(j+1);
            ++inode;
          }
        }
      }
    }

    // add nodes of the top layer of this element
    for ( size_t iz = 0 ; iz < m_nz(iy)+1 ; ++iz ) {
      for ( size_t ix = 0 ; ix < m_nx(iy)+1 ; ++ix ) {
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

inline MatS FineLayer::conn()
{
  // allocate output
  MatS out(m_nelem, m_nne);

  // current element, number of element layers, starting nodes of each node layer
  size_t ielem = 0;
  size_t N     = static_cast<size_t>(m_nhy.size());
  size_t bot,mid,top;

  // loop over all element layers
  for ( size_t iy = 0 ; iy < N ; ++iy )
  {
    // - get: starting nodes of bottom and top layer
    bot = m_startNode(iy  );
    top = m_startNode(iy+1);
    // - get: starting nodes of the middle layer (if present)
    if ( iy <= (N-1)/2 ) mid = m_startNode(iy) + (m_nx(iy  )+1) * (m_nz(iy  )+1);
    else                 mid = m_startNode(iy) + (m_nx(iy-1)+1) * (m_nz(iy-1)+1);

    // - define connectivity: no coarsening/refinement
    if ( m_refine(iy) == -1 )
    {
      for ( size_t iz = 0 ; iz < m_nz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nx(iy) ; ++ix ) {
          out(ielem,0) = bot + (ix  ) + (iz  ) * (m_nx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (iz  ) * (m_nx(iy)+1);
          out(ielem,2) = top + (ix+1) + (iz  ) * (m_nx(iy)+1);
          out(ielem,3) = top + (ix  ) + (iz  ) * (m_nx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (iz+1) * (m_nx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (iz+1) * (m_nx(iy)+1);
          out(ielem,6) = top + (ix+1) + (iz+1) * (m_nx(iy)+1);
          out(ielem,7) = top + (ix  ) + (iz+1) * (m_nx(iy)+1);
          ielem++;
        }
      }
    }

    // - define connectivity: refinement along the x-direction (below the middle layer)
    else if ( m_refine(iy) == 0 and iy <= (N-1)/2 )
    {
      for ( size_t iz = 0 ; iz < m_nz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nx(iy) ; ++ix ) {
          // -- bottom element
          out(ielem,0) = bot + (  ix  ) + (iz  ) * (  m_nx(iy)+1);
          out(ielem,1) = bot + (  ix+1) + (iz  ) * (  m_nx(iy)+1);
          out(ielem,2) = mid + (2*ix+1) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,3) = mid + (2*ix  ) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,4) = bot + (  ix  ) + (iz+1) * (  m_nx(iy)+1);
          out(ielem,5) = bot + (  ix+1) + (iz+1) * (  m_nx(iy)+1);
          out(ielem,6) = mid + (2*ix+1) + (iz+1) * (2*m_nx(iy)  );
          out(ielem,7) = mid + (2*ix  ) + (iz+1) * (2*m_nx(iy)  );
          ielem++;
          // -- top-right element
          out(ielem,0) = bot + (  ix+1) + (iz  ) * (  m_nx(iy)+1);
          out(ielem,1) = top + (3*ix+3) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,2) = top + (3*ix+2) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,3) = mid + (2*ix+1) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,4) = bot + (  ix+1) + (iz+1) * (  m_nx(iy)+1);
          out(ielem,5) = top + (3*ix+3) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,6) = top + (3*ix+2) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,7) = mid + (2*ix+1) + (iz+1) * (2*m_nx(iy)  );
          ielem++;
          // -- top-center element
          out(ielem,0) = mid + (2*ix  ) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,1) = mid + (2*ix+1) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,2) = top + (3*ix+2) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,3) = top + (3*ix+1) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,4) = mid + (2*ix  ) + (iz+1) * (2*m_nx(iy)  );
          out(ielem,5) = mid + (2*ix+1) + (iz+1) * (2*m_nx(iy)  );
          out(ielem,6) = top + (3*ix+2) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,7) = top + (3*ix+1) + (iz+1) * (3*m_nx(iy)+1);
          ielem++;
          // -- top-left element
          out(ielem,0) = bot + (  ix  ) + (iz  ) * (  m_nx(iy)+1);
          out(ielem,1) = mid + (2*ix  ) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,2) = top + (3*ix+1) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,3) = top + (3*ix  ) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,4) = bot + (  ix  ) + (iz+1) * (  m_nx(iy)+1);
          out(ielem,5) = mid + (2*ix  ) + (iz+1) * (2*m_nx(iy)  );
          out(ielem,6) = top + (3*ix+1) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,7) = top + (3*ix  ) + (iz+1) * (3*m_nx(iy)+1);
          ielem++;
        }
      }
    }

    // - define connectivity: coarsening along the x-direction (above the middle layer)
    else if ( m_refine(iy) == 0 and iy > (N-1)/2 )
    {
      for ( size_t iz = 0 ; iz < m_nz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nx(iy) ; ++ix ) {
          // -- lower-left element
          out(ielem,0) = bot + (3*ix  ) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,1) = bot + (3*ix+1) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,2) = mid + (2*ix  ) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,3) = top + (  ix  ) + (iz  ) * (  m_nx(iy)+1);
          out(ielem,4) = bot + (3*ix  ) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,5) = bot + (3*ix+1) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,6) = mid + (2*ix  ) + (iz+1) * (2*m_nx(iy)  );
          out(ielem,7) = top + (  ix  ) + (iz+1) * (  m_nx(iy)+1);
          ielem++;
          // -- lower-center element
          out(ielem,0) = bot + (3*ix+1) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,1) = bot + (3*ix+2) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,2) = mid + (2*ix+1) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,3) = mid + (2*ix  ) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,4) = bot + (3*ix+1) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,5) = bot + (3*ix+2) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,6) = mid + (2*ix+1) + (iz+1) * (2*m_nx(iy)  );
          out(ielem,7) = mid + (2*ix  ) + (iz+1) * (2*m_nx(iy)  );
          ielem++;
          // -- lower-right element
          out(ielem,0) = bot + (3*ix+2) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,1) = bot + (3*ix+3) + (iz  ) * (3*m_nx(iy)+1);
          out(ielem,2) = top + (  ix+1) + (iz  ) * (  m_nx(iy)+1);
          out(ielem,3) = mid + (2*ix+1) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,4) = bot + (3*ix+2) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,5) = bot + (3*ix+3) + (iz+1) * (3*m_nx(iy)+1);
          out(ielem,6) = top + (  ix+1) + (iz+1) * (  m_nx(iy)+1);
          out(ielem,7) = mid + (2*ix+1) + (iz+1) * (2*m_nx(iy)  );
          ielem++;
          // -- upper element
          out(ielem,0) = mid + (2*ix  ) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,1) = mid + (2*ix+1) + (iz  ) * (2*m_nx(iy)  );
          out(ielem,2) = top + (  ix+1) + (iz  ) * (  m_nx(iy)+1);
          out(ielem,3) = top + (  ix  ) + (iz  ) * (  m_nx(iy)+1);
          out(ielem,4) = mid + (2*ix  ) + (iz+1) * (2*m_nx(iy)  );
          out(ielem,5) = mid + (2*ix+1) + (iz+1) * (2*m_nx(iy)  );
          out(ielem,6) = top + (  ix+1) + (iz+1) * (  m_nx(iy)+1);
          out(ielem,7) = top + (  ix  ) + (iz+1) * (  m_nx(iy)+1);
          ielem++;
        }
      }
    }

    // - define connectivity: refinement along the z-direction (below the middle layer)
    else if ( m_refine(iy) == 2 and iy <= (N-1)/2 )
    {
      for ( size_t iz = 0 ; iz < m_nz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nx(iy) ; ++ix ) {
          // -- bottom element
          out(ielem,0) = bot + (ix  ) +    iz    * (m_nx(iy)+1);
          out(ielem,1) = bot + (ix+1) +    iz    * (m_nx(iy)+1);
          out(ielem,2) = bot + (ix+1) + (  iz+1) * (m_nx(iy)+1);
          out(ielem,3) = bot + (ix  ) + (  iz+1) * (m_nx(iy)+1);
          out(ielem,4) = mid + (ix  ) +  2*iz    * (m_nx(iy)+1);
          out(ielem,5) = mid + (ix+1) +  2*iz    * (m_nx(iy)+1);
          out(ielem,6) = mid + (ix+1) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,7) = mid + (ix  ) + (2*iz+1) * (m_nx(iy)+1);
          ielem++;
          // -- top-back element
          out(ielem,0) = mid + (ix  ) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,1) = mid + (ix+1) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,2) = top + (ix+1) + (3*iz+2) * (m_nx(iy)+1);
          out(ielem,3) = top + (ix  ) + (3*iz+2) * (m_nx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (  iz+1) * (m_nx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (  iz+1) * (m_nx(iy)+1);
          out(ielem,6) = top + (ix+1) + (3*iz+3) * (m_nx(iy)+1);
          out(ielem,7) = top + (ix  ) + (3*iz+3) * (m_nx(iy)+1);
          ielem++;
          // -- top-center element
          out(ielem,0) = mid + (ix  ) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,1) = mid + (ix+1) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,2) = top + (ix+1) + (3*iz+1) * (m_nx(iy)+1);
          out(ielem,3) = top + (ix  ) + (3*iz+1) * (m_nx(iy)+1);
          out(ielem,4) = mid + (ix  ) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,5) = mid + (ix+1) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,6) = top + (ix+1) + (3*iz+2) * (m_nx(iy)+1);
          out(ielem,7) = top + (ix  ) + (3*iz+2) * (m_nx(iy)+1);
          ielem++;
          // -- top-front element
          out(ielem,0) = bot + (ix  ) + (  iz  ) * (m_nx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (  iz  ) * (m_nx(iy)+1);
          out(ielem,2) = top + (ix+1) + (3*iz  ) * (m_nx(iy)+1);
          out(ielem,3) = top + (ix  ) + (3*iz  ) * (m_nx(iy)+1);
          out(ielem,4) = mid + (ix  ) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,5) = mid + (ix+1) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,6) = top + (ix+1) + (3*iz+1) * (m_nx(iy)+1);
          out(ielem,7) = top + (ix  ) + (3*iz+1) * (m_nx(iy)+1);
          ielem++;
        }
      }
    }

    // - define connectivity: coarsening along the z-direction (above the middle layer)
    else if ( m_refine(iy) == 2 and iy > (N-1)/2 )
    {
      for ( size_t iz = 0 ; iz < m_nz(iy) ; ++iz ) {
        for ( size_t ix = 0 ; ix < m_nx(iy) ; ++ix ) {
          // -- bottom-front element
          out(ielem,0) = bot + (ix  ) + (3*iz  ) * (m_nx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (3*iz  ) * (m_nx(iy)+1);
          out(ielem,2) = top + (ix+1) + (  iz  ) * (m_nx(iy)+1);
          out(ielem,3) = top + (ix  ) + (  iz  ) * (m_nx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (3*iz+1) * (m_nx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (3*iz+1) * (m_nx(iy)+1);
          out(ielem,6) = mid + (ix+1) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,7) = mid + (ix  ) + (2*iz  ) * (m_nx(iy)+1);
          ielem++;
          // -- bottom-center element
          out(ielem,0) = bot + (ix  ) + (3*iz+1) * (m_nx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (3*iz+1) * (m_nx(iy)+1);
          out(ielem,2) = mid + (ix+1) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,3) = mid + (ix  ) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (3*iz+2) * (m_nx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (3*iz+2) * (m_nx(iy)+1);
          out(ielem,6) = mid + (ix+1) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,7) = mid + (ix  ) + (2*iz+1) * (m_nx(iy)+1);
          ielem++;
          // -- bottom-back element
          out(ielem,0) = bot + (ix  ) + (3*iz+2) * (m_nx(iy)+1);
          out(ielem,1) = bot + (ix+1) + (3*iz+2) * (m_nx(iy)+1);
          out(ielem,2) = mid + (ix+1) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,3) = mid + (ix  ) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,4) = bot + (ix  ) + (3*iz+3) * (m_nx(iy)+1);
          out(ielem,5) = bot + (ix+1) + (3*iz+3) * (m_nx(iy)+1);
          out(ielem,6) = top + (ix+1) + (  iz+1) * (m_nx(iy)+1);
          out(ielem,7) = top + (ix  ) + (  iz+1) * (m_nx(iy)+1);
          ielem++;
          // -- top element
          out(ielem,0) = mid + (ix  ) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,1) = mid + (ix+1) + (2*iz  ) * (m_nx(iy)+1);
          out(ielem,2) = top + (ix+1) + (  iz  ) * (m_nx(iy)+1);
          out(ielem,3) = top + (ix  ) + (  iz  ) * (m_nx(iy)+1);
          out(ielem,4) = mid + (ix  ) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,5) = mid + (ix+1) + (2*iz+1) * (m_nx(iy)+1);
          out(ielem,6) = top + (ix+1) + (  iz+1) * (m_nx(iy)+1);
          out(ielem,7) = top + (ix  ) + (  iz+1) * (m_nx(iy)+1);
          ielem++;
        }
      }
    }
  }

  return out;
}

// =================================================================================================

}}} // namespace ...

// =================================================================================================

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHQUAD4_H
#define GOOSEFEM_MESHQUAD4_H

#include "Macros.h"
#include "Mesh.h"

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {

// ========================================== REGULAR MESH =========================================

class Regular
{
private:
  size_t m_nx;     // number of 'pixel' horizontal direction
  size_t m_ny;     // number of 'pixel' vertical direction
  double m_Lx;     // length in horizontal direction
  double m_Ly;     // length in vertical direction
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=4;  // number of nodes-per-element
  size_t m_ndim=2; // number of dimensions

public:
  // Regular            (const Regular &) = default;
  // Regular& operator= (const Regular &) = default;
  // Regular(){};
  // mesh with "nx" pixels in horizontal direction (total length "Lx"), idem in vertical direction
  Regular(size_t nx, size_t ny, double Lx=1., double Ly=1.);

  size_t nelem() { return m_nelem;}; // number of elements
  size_t nnode() { return m_nnode;}; // number of nodes
  size_t nne  () { return m_nne;  }; // number of nodes-per-element
  size_t ndim () { return m_ndim; }; // number of dimensions
  MatD   coor         ();            // nodal positions [ nnode , ndim ]
  MatS   conn         ();            // connectivity    [ nelem , nne  ]
  ColS   nodesBottom  ();            // nodes along the bottom edge
  ColS   nodesTop     ();            // nodes along the top    edge
  ColS   nodesLeft    ();            // nodes along the left   edge
  ColS   nodesRight   ();            // nodes along the right  edge
  MatS   nodesPeriodic();            // periodic node pairs [ : , 2 ]: ( independent , dependent )
  size_t nodesRef     ();            // lower-left node, to be used as reference for periodicity
  MatS   dofs         ();            // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic ();            // DOF-numbers for each component of each node (sequential)
};

// ====================== MESH WITH A FINE LAYER THAT EXPONENTIALLY COARSENS =======================


// ========================================== SOURCE CODE ==========================================

Regular::Regular(size_t nx, size_t ny, double Lx, double Ly): m_nx(nx), m_ny(ny), m_Lx(Lx), m_Ly(Ly)
{
  assert( m_nx >= 1 );
  assert( m_ny >= 1 );

  m_nnode = (m_nx+1) * (m_ny+1);
  m_nelem =  m_nx    *  m_ny   ;
}

// -------------------------------------------------------------------------------------------------

MatD Regular::coor()
{
  MatD coor( m_nnode , m_ndim );

  ColD x = ColD::LinSpaced( m_nx+1 , 0.0 , m_Lx );
  ColD y = ColD::LinSpaced( m_ny+1 , 0.0 , m_Ly );

  size_t inode = 0;

  for ( size_t row = 0 ; row < m_ny+1 ; ++row ) {
    for ( size_t col = 0 ; col < m_nx+1 ; ++col ) {
      coor(inode,0) = x(col);
      coor(inode,1) = y(row);
      ++inode;
    }
  }

  return coor;
}

// -------------------------------------------------------------------------------------------------

MatS Regular::conn()
{
  MatS conn( m_nelem , m_nne );

  size_t ielem = 0;

  for ( size_t row = 0 ; row < m_ny ; ++row ) {
    for ( size_t col = 0 ; col < m_nx ; ++col ) {
      conn(ielem,0) = (row+0)*(m_nx+1)+col+0;
      conn(ielem,1) = (row+0)*(m_nx+1)+col+1;
      conn(ielem,3) = (row+1)*(m_nx+1)+col+0;
      conn(ielem,2) = (row+1)*(m_nx+1)+col+1;
      ++ielem;
    }
  }

  return conn;
}

// -------------------------------------------------------------------------------------------------

ColS Regular::nodesBottom()
{
  ColS nodes(m_nx+1);

  for ( size_t col = 0 ; col < m_nx+1 ; ++col ) nodes(col) = col;

  return nodes;
}

// -------------------------------------------------------------------------------------------------

ColS Regular::nodesTop()
{
  ColS nodes(m_nx+1);

  for ( size_t col = 0 ; col < m_nx+1 ; ++col ) nodes(col) = col+m_ny*(m_nx+1);

  return nodes;
}

// -------------------------------------------------------------------------------------------------

ColS Regular::nodesLeft()
{
  ColS nodes(m_ny+1);

  for ( size_t row = 0 ; row < m_ny+1 ; ++row ) nodes(row) = row*(m_nx+1);

  return nodes;
}

// -------------------------------------------------------------------------------------------------

ColS Regular::nodesRight()
{
  ColS nodes(m_ny+1);

  for ( size_t row = 0 ; row < m_ny+1 ; ++row ) nodes(row) = row*(m_nx+1)+m_nx;

  return nodes;
}

// -------------------------------------------------------------------------------------------------

MatS Regular::nodesPeriodic()
{
  ColS bot = nodesBottom();
  ColS top = nodesTop   ();
  ColS rgt = nodesRight ();
  ColS lft = nodesLeft  ();

  MatS nodes( bot.size()-2 + lft.size()-2 + 3 , 2 );

  size_t i = 0;

  nodes(i,0) = bot(0); nodes(i,1) = bot(bot.size()-1); ++i;
  nodes(i,0) = bot(0); nodes(i,1) = top(top.size()-1); ++i;
  nodes(i,0) = bot(0); nodes(i,1) = top(0           ); ++i;

  for ( auto j = 1 ; j < bot.size()-1 ; ++j ) { nodes(i,0) = bot(j); nodes(i,1) = top(j); ++i; }
  for ( auto j = 1 ; j < lft.size()-1 ; ++j ) { nodes(i,0) = lft(j); nodes(i,1) = rgt(j); ++i; }

  return nodes;
}

// -------------------------------------------------------------------------------------------------

size_t Regular::nodesRef()
{
  return 0;
}

// -------------------------------------------------------------------------------------------------

MatS Regular::dofs()
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// -------------------------------------------------------------------------------------------------

MatS Regular::dofsPeriodic()
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

} // namespace Quad4
} // namespace Mesh
} // namespace GooseFEM

#endif

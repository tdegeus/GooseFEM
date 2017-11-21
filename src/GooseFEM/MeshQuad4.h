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
  size_t m_nx;     // number of 'pixels' horizontal direction (length == "m_nx * m_h")
  size_t m_ny;     // number of 'pixels' vertical direction   (length == "m_ny * m_h")
  double m_h;      // size of the element edge (equal in both directions)
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=4;  // number of nodes-per-element
  size_t m_ndim=2; // number of dimensions

public:
  // mesh with "nx" pixels in horizontal direction, "ny" in vertical direction and "h" the edge size
  Regular(size_t nx, size_t ny, double h=1.);

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
  size_t nodeOrigin   ();            // bottom-left node, to be used as reference for periodicity
  MatS   dofs         ();            // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic ();            // DOF-numbers for each component of each node (sequential)
};

// ====================== MESH WITH A FINE LAYER THAT EXPONENTIALLY COARSENS =======================

class FineLayer
{
private:
  double m_h;      // base size of the element edge (equal in both directions)
  size_t m_nx;     // number of elements in vertical direction
  ColS   m_nh;     // element size in vertical direction (number of time "h")
  ColS   m_startNode;  // start node of each row
  ColS   m_startElem;  // start element of each row
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=4;  // number of nodes-per-element
  size_t m_ndim=2; // number of dimensions

public:
  // mesh with "nx" pixels in horizontal direction, "ny" in vertical direction and "h" the edge size
  FineLayer(size_t nx, size_t ny, double h=1., size_t nfine=0, size_t nskip=0);

  size_t nelem() { return m_nelem;}; // number of elements
  size_t nnode() { return m_nnode;}; // number of nodes
  size_t nne  () { return m_nne;  }; // number of nodes-per-element
  size_t ndim () { return m_ndim; }; // number of dimensions
  size_t shape(size_t i);            // actual shape in horizontal and vertical direction
  MatD   coor         ();            // nodal positions [ nnode , ndim ]
  MatS   conn         ();            // connectivity    [ nelem , nne  ]
  ColS   elementsFine ();            // elements in the middle, fine, layer
  ColS   nodesBottom  ();            // nodes along the bottom edge
  ColS   nodesTop     ();            // nodes along the top    edge
  ColS   nodesLeft    ();            // nodes along the left   edge
  ColS   nodesRight   ();            // nodes along the right  edge
  MatS   nodesPeriodic();            // periodic node pairs [ : , 2 ]: ( independent , dependent )
  size_t nodeOrigin   ();            // bottom-left node, to be used as reference for periodicity
  MatS   dofs         ();            // DOF-numbers for each component of each node (sequential)
  MatS   dofsPeriodic ();            // DOF-numbers for each component of each node (sequential)
};


// ===================================== SOURCE CODE : REGULAR =====================================

inline Regular::Regular(size_t nx, size_t ny, double h): m_nx(nx), m_ny(ny), m_h(h)
{
  assert( m_nx >= 1 );
  assert( m_ny >= 1 );

  m_nnode = (m_nx+1) * (m_ny+1);
  m_nelem =  m_nx    *  m_ny   ;
}

// -------------------------------------------------------------------------------------------------

inline MatD Regular::coor()
{
  MatD coor( m_nnode , m_ndim );

  ColD x = ColD::LinSpaced( m_nx+1 , 0.0 , m_h * static_cast<double>(m_nx) );
  ColD y = ColD::LinSpaced( m_ny+1 , 0.0 , m_h * static_cast<double>(m_ny) );

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

inline MatS Regular::conn()
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

inline ColS Regular::nodesBottom()
{
  ColS nodes(m_nx+1);

  for ( size_t col = 0 ; col < m_nx+1 ; ++col ) nodes(col) = col;

  return nodes;
}

// -------------------------------------------------------------------------------------------------

inline ColS Regular::nodesTop()
{
  ColS nodes(m_nx+1);

  for ( size_t col = 0 ; col < m_nx+1 ; ++col ) nodes(col) = col+m_ny*(m_nx+1);

  return nodes;
}

// -------------------------------------------------------------------------------------------------

inline ColS Regular::nodesLeft()
{
  ColS nodes(m_ny+1);

  for ( size_t row = 0 ; row < m_ny+1 ; ++row ) nodes(row) = row*(m_nx+1);

  return nodes;
}

// -------------------------------------------------------------------------------------------------

inline ColS Regular::nodesRight()
{
  ColS nodes(m_ny+1);

  for ( size_t row = 0 ; row < m_ny+1 ; ++row ) nodes(row) = row*(m_nx+1)+m_nx;

  return nodes;
}

// -------------------------------------------------------------------------------------------------

inline MatS Regular::nodesPeriodic()
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

inline size_t Regular::nodeOrigin()
{
  return 0;
}

// -------------------------------------------------------------------------------------------------

inline MatS Regular::dofs()
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// -------------------------------------------------------------------------------------------------

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

// ==================================== SOURCE CODE : FINELAYER ====================================

inline FineLayer::FineLayer(size_t nx, size_t ny, double h, size_t nfine, size_t nskip)
{
  assert( nx >= 1 );
  assert( ny >= 1 );

  m_nx = nx;
  m_h  = h;

  // compute the element size : based on symmetric half
  // --------------------------------------------------

  // convert height to half of the height
  if ( ny % 2 == 0 ) ny =  ny   /2;
  else               ny = (ny+1)/2;

  // check the number of fine layers from the center
  assert( nfine <= ny );

  // define arrays to determine to coarsen
  ColS n   (ny+1); // size of the element (number of times "h")
  ColS nsum(ny+1); // current size of the mesh in y-direction (in number of times "nsum")

  // initialize all elements of size "1"
  size_t next = 1;
  n.setOnes();

  // size of the mesh in y-direction after "i" element
  // - initialize
  nsum(0) = 1;
  // - compute based on the initially equi-sized elements
  for ( size_t i = 1 ; i < ny+1 ; ++i )
    nsum(i) = nsum(i-1) + n(i-1);

  // loop in vertical direction and check to coarsen; rules:
  // - the size of the element cannot be greater than the distance
  // - the size of the element should match the horizontal direction
  // - skip a certain number of elements
  size_t min = 1;
  if ( nfine > min ) min = nfine;
  for ( size_t i = min ; i < ny+1 ; ++i )
  {
    if ( 3*n(i-1) <= nsum(i-1) and m_nx % ( 3*next ) == 0 )
    {
      n   (i) = n   (i-1) * 2;
      nsum(i) = nsum(i-1) + n(i);
      next   *= 3;
    }
    else
    {
      n   (i) = next;
      nsum(i) = nsum(i-1) + n(i);
    }
  }

  // skip the last "nskip" coarsening steps
  // - counter : number of steps skipped
  size_t iskip = 0;
  // - loop to skip
  for ( size_t i = ny+1 ; i-- > 1 ; )
  {
    // -- check -> quit if the number of steps is reached
    if ( iskip >= nskip )
      break;
    // -- if coarsening step found -> remove
    if ( n(i) % 2 == 0 )
    {
      for ( size_t j = i ; j < ny+1 ; ++j ) n(j) = n(i-1);
      iskip++;
    }
  }
  // - update the height after "i" elements
  for ( size_t i = 1 ; i < ny+1 ; ++i )
    nsum(i) = nsum(i-1) + n(i);

  // truncate such that the height does not exceed "ny"
  size_t N;
  for ( N = 0 ; N < ny+1 ; ++N )
    if ( nsum(N) > ny+1 )
      break;

  // truncate
  n.conservativeResize(N);

  // compute mesh dimensions : based on the entire mesh
  // --------------------------------------------------

  // symmetrize
  // - get size
  N = static_cast<size_t>(n.size());
  // - allocate
  m_nh.conservativeResize(2*N-1);
  // - store symmetrized
  for ( size_t i = 0 ; i < N ; i++ ) m_nh(  i  ) = n(N-1-i);
  for ( size_t i = 1 ; i < N ; i++ ) m_nh(N+i-1) = n(    i);

  // allocate counters
  m_startElem.conservativeResize(m_nh.size()  );
  m_startNode.conservativeResize(m_nh.size()+1);

  // compute the dimensions of the mesh
  // - initialize
  m_nnode        = 0;
  m_nelem        = 0;
  N              = static_cast<size_t>(m_nh.size());
  m_startNode(0) = 0;

  // loop from bottom to middle : elements become finer
  for ( size_t i = 0 ; i < (N-1)/2 ; ++i )
  {
    // - store the first element of the row
    m_startElem(i) = m_nelem;
    // - add the nodes of this row
    if ( m_nh(i)%2 == 0 ) { m_nnode += 3 * m_nx/(3*m_nh(i)/2) + 1; }
    else                  { m_nnode +=     m_nx/   m_nh(i)    + 1; }
    // - add the elements of this row
    if ( m_nh(i)%2 == 0 ) { m_nelem += 4 * m_nx/(3*m_nh(i)/2); }
    else                  { m_nelem +=     m_nx/   m_nh(i)   ; }
    // - store the starting node of the next row
    m_startNode(i+1) = m_nnode;
  }

  // loop from middle to top : elements become coarser
  for ( size_t i = (N-1)/2 ; i < N ; ++i )
  {
    // - store the first element of the row
    m_startElem(i) = m_nelem;
    // - add the nodes of this row
    if ( m_nh(i)%2 == 0 ) { m_nnode += 1 + m_nx/(m_nh(i)/2) + 2 * m_nx/(3*m_nh(i)/2); }
    else                  { m_nnode += 1 + m_nx/ m_nh(i);                             }
    // - add the elements of this row
    if ( m_nh(i)%2 == 0 ) { m_nelem += 4 * m_nx/(3*m_nh(i)/2); }
    else                  { m_nelem +=     m_nx/   m_nh(i)   ; }
    // - store the starting node of the next row
    m_startNode(i+1) = m_nnode;
  }

  // add the top row of nodes to the number of rows (starting node already stored)
  if ( m_nh(N-1)%2 == 0 ) { m_nnode += 1 + m_nx/(3*m_nh(N-1)/2); }
  else                    { m_nnode += 1 + m_nx/   m_nh(N-1);    }
}

// -------------------------------------------------------------------------------------------------

inline size_t FineLayer::shape(size_t i)
{
  // check the index
  assert( i >= 0 and i < 2 );

  // horizontal direction
  if ( i == 0 )
    return m_nx;

  // vertical direction
  // - zero-initialize
  size_t n = 0;
  // - cumulative sum of the size per row
  for ( size_t i = 0 ; i < static_cast<size_t>(m_nh.size()) ; ++i )
    n += m_nh(i);

  return n;
}

// -------------------------------------------------------------------------------------------------

inline MatD FineLayer::coor()
{
  // allocate output
  MatD out( m_nnode , m_ndim );

  // initialize position in horizontal direction (of the finest nodes)
  ColD x = ColD::LinSpaced( m_nx+1 , 0.0 , m_h * static_cast<double>(m_nx) );

  // zero-initialize current node and height: loop from top to bottom; number of rows
  double h     = 0;
  size_t inode = 0;
  size_t N     = static_cast<size_t>(m_nh.size());

  // lower half
  for ( size_t i = 0 ; i < (N-1)/2 ; ++i )
  {
    // - refinement row
    if ( m_nh(i)%2 == 0 )
    {
      // -- bottom row: coarse nodes
      for ( size_t j = 0 ; j < m_nx+1 ; j += 3*m_nh(i)/2 )
      {
        out(inode,0) = x(j); out(inode,1) = h*m_h; inode++;
      }
      h += static_cast<double>(m_nh(i)/2);
      // -- middle row: part the fine nodes
      for ( size_t j = 0 ; j < m_nx   ; j += 3*m_nh(i)/2 )
      {
        out(inode,0) = x(j+m_nh(i)/2); out(inode,1) = h*m_h; inode++;
        out(inode,0) = x(j+m_nh(i)  ); out(inode,1) = h*m_h; inode++;
      }
      h += static_cast<double>(m_nh(i)/2);
    }
    // - normal row
    else
    {
      for ( size_t j = 0 ; j < m_nx+1 ; j += m_nh(i) )
      {
        out(inode,0) = x(j); out(inode,1) = h*m_h; inode++;
      }
      h += static_cast<double>(m_nh(i));
    }
  }

  // top half
  for ( size_t i = (N-1)/2 ; i < N ; ++i )
  {
    // - coarsening row
    if ( m_nh(i)%2 == 0 )
    {
      // -- bottom row: fine nodes
      for ( size_t j = 0 ; j < m_nx+1 ; j += m_nh(i)/2 )
      {
        out(inode,0) = x(j); out(inode,1) = h*m_h; inode++;
      }
      h += static_cast<double>(m_nh(i)/2);
      // -- middle row: part the fine nodes
      for ( size_t j = 0 ; j < m_nx   ; j += 3*m_nh(i)/2 )
      {
        out(inode,0) = x(j+m_nh(i)/2); out(inode,1) = h*m_h; inode++;
        out(inode,0) = x(j+m_nh(i)  ); out(inode,1) = h*m_h; inode++;
      }
      h += static_cast<double>(m_nh(i)/2);
    }
    // - normal row
    else
    {
      for ( size_t j = 0 ; j < m_nx+1 ; j += m_nh(i) )
      {
        out(inode,0) = x(j); out(inode,1) = h*m_h; inode++;
      }
      h += static_cast<double>(m_nh(i));
    }
  }

  // top row
  if ( m_nh(N-1)%2 == 0 )
  {
    for ( size_t j = 0 ; j < m_nx+1 ; j += 3*m_nh(N-1)/2 )
    {
      out(inode,0) = x(j); out(inode,1) = h*m_h; inode++;
    }
  }
  else
  {
    for ( size_t j = 0 ; j < m_nx+1 ; j += m_nh(N-1) )
    {
      out(inode,0) = x(j); out(inode,1) = h*m_h; inode++;
    }
  }

  return out;
}

// // -------------------------------------------------------------------------------------------------

inline MatS FineLayer::conn()
{
  // allocate output
  MatS out( m_nelem , m_nne );

  // current element, number of rows, beginning nodes
  size_t ielem = 0;
  size_t N     = static_cast<size_t>(m_nh.size());
  size_t bot,mid,top;

  // bottom half
  for ( size_t i = 0 ; i < (N-1)/2 ; ++i )
  {
    // - starting nodes of bottom and top layer (and middle layer, only relevant for even shape)
    if ( true           ) bot = m_startNode(i );
    if ( m_nh(i)%2 == 0 ) mid = m_startNode(i )+m_nx/(3*m_nh(i)/2)+1;
    if ( true           ) top = m_startNode(i+1);

    // - even sized shape: refinement layer
    if ( m_nh(i)%2 == 0 )
    {
      for ( size_t j = 0 ; j < m_nx/(3*m_nh(i)/2) ; ++j )
      {
        // -- lower
        out(ielem,0) = j    +bot;
        out(ielem,1) = j  +1+bot;
        out(ielem,2) = j*2+1+mid;
        out(ielem,3) = j*2  +mid;
        ielem++;
        // -- upper right
        out(ielem,0) = j  +1+bot;
        out(ielem,1) = j*3+3+top;
        out(ielem,2) = j*3+2+top;
        out(ielem,3) = j*2+1+mid;
        ielem++;
        // -- upper middle
        out(ielem,0) = j*2  +mid;
        out(ielem,1) = j*2+1+mid;
        out(ielem,2) = j*3+2+top;
        out(ielem,3) = j*3+1+top;
        ielem++;
        // -- upper left
        out(ielem,0) = j    +bot;
        out(ielem,1) = j*2  +mid;
        out(ielem,2) = j*3+1+top;
        out(ielem,3) = j*3  +top;
        ielem++;
      }
    }
    // - odd sized shape: normal layer
    else
    {
      for ( size_t j = 0 ; j < m_nx/m_nh(i) ; ++j )
      {
        out(ielem,0) = j  +bot;
        out(ielem,1) = j+1+bot;
        out(ielem,2) = j+1+top;
        out(ielem,3) = j  +top;
        ielem++;
      }
    }
  }

  // top half
  for ( size_t i = (N-1)/2 ; i < N ; ++i )
  {
    // - starting nodes of bottom and top layer (and middle layer, only relevant for even shape)
    if ( true           ) bot = m_startNode(i );
    if ( m_nh(i)%2 == 0 ) mid = m_startNode(i )+m_nx/(m_nh(i)/2)+1;
    if ( true           ) top = m_startNode(i+1);

    // - even sized shape: refinement layer
    if ( m_nh(i)%2 == 0 )
    {
      for ( size_t j = 0 ; j < m_nx/(3*m_nh(i)/2) ; ++j )
      {
        // -- lower left
        out(ielem,0) = j*3  +bot;
        out(ielem,1) = j*3+1+bot;
        out(ielem,2) = j*2  +mid;
        out(ielem,3) = j    +top;
        ielem++;
        // -- lower middle
        out(ielem,0) = j*3+1+bot;
        out(ielem,1) = j*3+2+bot;
        out(ielem,2) = j*2+1+mid;
        out(ielem,3) = j*2  +mid;
        ielem++;
        // -- lower right
        out(ielem,0) = j*3+2+bot;
        out(ielem,1) = j*3+3+bot;
        out(ielem,2) = j  +1+top;
        out(ielem,3) = j*2+1+mid;
        ielem++;
        // -- upper
        out(ielem,0) = j    +top;
        out(ielem,1) = j*2  +mid;
        out(ielem,2) = j*2+1+mid;
        out(ielem,3) = j  +1+top;
        ielem++;
      }
    }
    // - odd sized shape: normal layer
    else
    {
      for ( size_t j = 0 ; j < m_nx/m_nh(i) ; ++j )
      {
        out(ielem,0) = j  +bot;
        out(ielem,1) = j+1+bot;
        out(ielem,2) = j+1+top;
        out(ielem,3) = j  +top;
        ielem++;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline ColS FineLayer::elementsFine()
{
  ColS out(m_nx);

  size_t j = m_startElem((m_startElem.size()-1)/2);

  for ( size_t i = 0 ; i < m_nx ; ++i ) out(i) = i+j;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline ColS FineLayer::nodesBottom()
{
  ColS nodes;

  if ( m_nh(0)%2 == 0 ) nodes.conservativeResize(m_nx/(3*m_nh(0)/2)+1);
  else                  nodes.conservativeResize(m_nx/   m_nh(0)   +1);

  for ( size_t i = 0 ; i < static_cast<size_t>(nodes.size()) ; ++i ) nodes(i) = i;

  return nodes;
}

// -------------------------------------------------------------------------------------------------

inline ColS FineLayer::nodesTop()
{
  ColS nodes;

  if ( m_nh(0)%2 == 0 ) nodes.conservativeResize(m_nx/(3*m_nh(0)/2)+1);
  else                  nodes.conservativeResize(m_nx/   m_nh(0)   +1);

  size_t j = m_startNode(m_startNode.size()-1);

  for ( size_t i = 0 ; i < static_cast<size_t>(nodes.size()) ; ++i ) nodes(i) = i+j;

  return nodes;
}

// -------------------------------------------------------------------------------------------------

inline ColS FineLayer::nodesLeft()
{
  return m_startNode;
}

// -------------------------------------------------------------------------------------------------

inline ColS FineLayer::nodesRight()
{
  size_t N = static_cast<size_t>(m_nh.size());
  ColS   nodes(N+1);

  // bottom half
  for ( size_t i = 0 ; i < (N-1)/2 ; ++i )
  {
    if ( m_nh(i)%2 == 0 ) nodes(i) = m_startNode(i) + m_nx/(3*m_nh(i)/2);
    else                  nodes(i) = m_startNode(i) + m_nx/   m_nh(i)   ;
  }

  // top half
  for ( size_t i = (N-1)/2 ; i < N ; ++i )
  {
    if ( m_nh(i)%2 == 0 ) nodes(i) = m_startNode(i) + m_nx/(m_nh(i)/2);
    else                  nodes(i) = m_startNode(i) + m_nx/ m_nh(i)   ;
  }

  // last node
  nodes(N) = m_nnode-1;

  return nodes;
}

// -------------------------------------------------------------------------------------------------

inline MatS FineLayer::nodesPeriodic()
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

inline size_t FineLayer::nodeOrigin()
{
  return 0;
}

// -------------------------------------------------------------------------------------------------

inline MatS FineLayer::dofs()
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// -------------------------------------------------------------------------------------------------

inline MatS FineLayer::dofsPeriodic()
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

}}} // namespace GooseFEM::Mesh::Quad4

#endif

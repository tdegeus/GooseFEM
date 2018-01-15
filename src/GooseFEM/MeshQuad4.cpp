/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHQUAD4_CPP
#define GOOSEFEM_MESHQUAD4_CPP

// -------------------------------------------------------------------------------------------------

#include "MeshQuad4.h"

// ===================================== GooseFEM::Mesh::Quad4 =====================================

namespace GooseFEM {
namespace Mesh {
namespace Quad4 {

// ===================================== CLASS - REGULAR MESH ======================================

// ------------------------------------------ constructor ------------------------------------------

inline Regular::Regular(size_t nx, size_t ny, double h):
m_nx(nx), m_ny(ny), m_h(h)
{
  assert( m_nx >= 1 );
  assert( m_ny >= 1 );

  m_nnode = (m_nx+1) * (m_ny+1);
  m_nelem =  m_nx    *  m_ny   ;
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

  size_t inode = 0;

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy ) {
    for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix ) {
      coor(inode,0) = x(ix);
      coor(inode,1) = y(iy);
      ++inode;
    }
  }

  return coor;
}

// ---------------------------- connectivity (node-numbers per element) ----------------------------

inline MatS Regular::conn()
{
  MatS conn(m_nelem,m_nne);

  size_t ielem = 0;

  for ( size_t iy = 0 ; iy < m_ny ; ++iy ) {
    for ( size_t ix = 0 ; ix < m_nx ; ++ix ) {
      conn(ielem,0) = (iy+0)*(m_nx+1) + (ix+0);
      conn(ielem,1) = (iy+0)*(m_nx+1) + (ix+1);
      conn(ielem,3) = (iy+1)*(m_nx+1) + (ix+0);
      conn(ielem,2) = (iy+1)*(m_nx+1) + (ix+1);
      ++ielem;
    }
  }

  return conn;
}

// ------------------------------ node-numbers along the bottom edge -------------------------------

inline ColS Regular::nodesBottomEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = ix;

  return nodes;
}

// -------------------------------- node-numbers along the top edge --------------------------------

inline ColS Regular::nodesTopEdge()
{
  ColS nodes(m_nx+1);

  for ( size_t ix = 0 ; ix < m_nx+1 ; ++ix )
    nodes(ix) = ix + m_ny*(m_nx+1);

  return nodes;
}

// ----------------------------- node-numbers along the node left edge -----------------------------

inline ColS Regular::nodesLeftEdge()
{
  ColS nodes(m_ny+1);

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    nodes(iy) = iy*(m_nx+1);

  return nodes;
}

// ------------------------------- node-numbers along the right edge -------------------------------

inline ColS Regular::nodesRightEdge()
{
  ColS nodes(m_ny+1);

  for ( size_t iy = 0 ; iy < m_ny+1 ; ++iy )
    nodes(iy) = iy*(m_nx+1) + m_nx;

  return nodes;
}

// ----------------------------- node-number of the bottom-left corner -----------------------------

inline size_t Regular::nodesBottomLeftCorner()
{
  return 0;
}

// ---------------------------- node-number of the bottom-right corner -----------------------------

inline size_t Regular::nodesBottomRightCorner()
{
  return m_nx;
}

// ------------------------------ node-number of the top-left corner -------------------------------

inline size_t Regular::nodesTopLeftCorner()
{
  return m_ny*(m_nx+1);
}

// ------------------------------ node-number of the top-right corner ------------------------------

inline size_t Regular::nodesTopRightCorner()
{
  return (m_ny+1)*(m_nx+1) - 1;
}

// ----------------------------- node-number of the corners (aliases) ------------------------------

inline size_t Regular::nodesLeftBottomCorner()  { return nodesBottomLeftCorner();  }
inline size_t Regular::nodesLeftTopCorner()     { return nodesTopLeftCorner();     }
inline size_t Regular::nodesRightBottomCorner() { return nodesBottomRightCorner(); }
inline size_t Regular::nodesRightTopCorner()    { return nodesTopRightCorner();    }

// ------------------------------ node-numbers of periodic node-pairs ------------------------------

inline MatS Regular::nodesPeriodic()
{
  // edges
  ColS bot = nodesBottomEdge();
  ColS top = nodesTopEdge();
  ColS rgt = nodesRightEdge();
  ColS lft = nodesLeftEdge();

  // allocate nodal ties
  MatS nodes(bot.size()-2 + lft.size()-2 + 3, 2);

  // counter
  size_t i = 0;

  // tie all corners to one corner
  nodes(i,0) = nodesBottomLeftCorner(); nodes(i,1) = nodesBottomRightCorner(); ++i;
  nodes(i,0) = nodesBottomLeftCorner(); nodes(i,1) = nodesTopRightCorner();    ++i;
  nodes(i,0) = nodesBottomLeftCorner(); nodes(i,1) = nodesTopLeftCorner();     ++i;

  // tie edges to each other (exclude corners)
  for ( auto j = 1 ; j < bot.size()-1 ; ++j ) { nodes(i,0) = bot(j); nodes(i,1) = top(j); ++i; }
  for ( auto j = 1 ; j < lft.size()-1 ; ++j ) { nodes(i,0) = lft(j); nodes(i,1) = rgt(j); ++i; }

  return nodes;
}

// ------------------------------ node-number that lies in the origin ------------------------------

inline size_t Regular::nodesOrigin()
{
  return nodesBottomLeftCorner();
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

inline FineLayer::FineLayer(size_t nx, size_t ny, double h, size_t nfine, size_t nskip):
m_h(h), m_nx(nx)
{
  assert( nx >= 1 );
  assert( ny >= 1 );

  // local counters
  size_t N;

  // compute the element size : based on symmetric half
  // --------------------------------------------------

  // convert height to half of the height
  if ( ny % 2 == 0 ) ny =  ny   /2;
  else               ny = (ny+1)/2;

  // convert to half the number of fine layer (minimum 1)
  if ( nfine % 2 == 0 ) nfine =  nfine   /2 + 1;
  else                  nfine = (nfine+1)/2;
  if ( nfine < 1      ) nfine = 1;

  // check the number of fine layers from the center
  assert( nfine <= ny );

  // define arrays to determine to coarsen
  ColS n   (ny+1); // size of the element (number of times "h"): even sizes -> coarsen
  ColS nsum(ny+1); // current 'size' of the mesh in y-direction

  // initialize all elements of size "1"
  n.setOnes();

  // size of the mesh in y-direction after "i" elements
  // - initialize
  nsum(0) = 1;
  // - compute based on the current value of "n"
  for ( size_t i = 1 ; i < ny+1 ; ++i )
    nsum(i) = nsum(i-1) + n(i-1);

  // loop in vertical direction and check to coarsen; rules:
  // * the size of the element cannot be greater than the distance
  // * the size of the element should match the horizontal direction
  // * skip a certain number of layers that have the minimum size "1" (are fine)
  // - initialize next element size
  size_t next = 1;
  // - loop to coarsen
  for ( size_t i = nfine ; i < ny+1 ; ++i )
  {
    // -- set element size, based on whether to coarsen or not
    if ( 3*n(i-1) <= nsum(i-1) and m_nx % ( 3*next ) == 0 ) { n(i) = next * 2; next *= 3; }
    else                                                    { n(i) = next;                }
    // -- compute the total number of elements
    nsum(i) = nsum(i-1) + n(i);
  }

  // truncate to size "N" such that the height does not exceed "ny"
  for ( N = 0 ; N < ny+1 ; ++N )
    if ( nsum(N) > ny+1 )
      break;

  // fiddle a little bit to approximate the size as best as possible
  if ( N > 1 and N < ny )
  {
    // - set all sequential steps to last value, make sure not to use a coarsening step
    if ( n(N-1) % 2 == 0 ) for ( size_t i = N ; i < ny+1 ; ++i ) n(i) = n(i-1)/2*3;
    else                   for ( size_t i = N ; i < ny+1 ; ++i ) n(i) = n(i-1);
    // - fix total size
    for ( size_t i = N ; i < ny+1 ; ++i ) nsum(i) = nsum(i-1) + n(i);
    // - truncate such that the height does not exceed "ny"
    for ( N = 0 ; N < ny ; ++N )
      if ( nsum(N) > ny+1 )
        break;
    // - check to extend one
    if ( nsum(N) - ny < ny - nsum(N-1) )
      ++N;
  }

  // skip the last "nskip" coarsening steps
  if ( nskip > 0 )
  {
    // - counter : number of steps skipped
    size_t iskip = 0;
    // - loop to skip
    for ( size_t i = N ; i-- > 1 ; )
    {
      // -- check -> quit if the number of steps is reached
      if ( iskip >= nskip ) break;
      // -- check -> quit if the finest layer has been reached
      if ( n(i) == 1 ) break;
      // -- if coarsening step found -> remove
      if ( n(i) % 2 == 0 )
      {
        size_t k = n(i) / 2;
        for ( size_t j = i ; j < ny+1 ; ++j ) n(j) = k;
        iskip++;
      }
    }
    // - update the height after "i" elements
    for ( size_t i = 1 ; i < ny+1 ; ++i )
      nsum(i) = nsum(i-1) + n(i);
    // - truncate such that the height does not exceed "ny"
    for ( N = 0 ; N < ny+1 ; ++N )
      if ( nsum(N) > ny+1 )
        break;
  }

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

// ---------------------------- actual number of nodes in one direction ----------------------------

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

// --------------------------------- coordinates (nodal positions) ---------------------------------

inline MatD FineLayer::coor()
{
  // allocate output
  MatD out(m_nnode, m_ndim);

  // initialize position in horizontal direction (of the finest nodes)
  ColD x = ColD::LinSpaced(m_nx+1, 0.0, m_h*static_cast<double>(m_nx));

  // zero-initialize current node and height: loop from bottom to top; number of rows
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

// ---------------------------- connectivity (node-numbers per element) ----------------------------

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

// ------------------------------ element numbers of the middle layer ------------------------------

inline ColS FineLayer::elementsMiddleLayer()
{
  ColS out(m_nx);

  size_t j = m_startElem((m_startElem.size()-1)/2);

  for ( size_t i = 0 ; i < m_nx ; ++i ) out(i) = i+j;

  return out;
}

// ------------------------------ node-numbers along the bottom edge -------------------------------

inline ColS FineLayer::nodesBottomEdge()
{
  ColS nodes;

  if ( m_nh(0)%2 == 0 ) nodes.conservativeResize(m_nx/(3*m_nh(0)/2)+1);
  else                  nodes.conservativeResize(m_nx/   m_nh(0)   +1);

  for ( size_t i = 0 ; i < static_cast<size_t>(nodes.size()) ; ++i ) nodes(i) = i;

  return nodes;
}

// -------------------------------- node-numbers along the top edge --------------------------------

inline ColS FineLayer::nodesTopEdge()
{
  ColS nodes;

  if ( m_nh(0)%2 == 0 ) nodes.conservativeResize(m_nx/(3*m_nh(0)/2)+1);
  else                  nodes.conservativeResize(m_nx/   m_nh(0)   +1);

  size_t j = m_startNode(m_startNode.size()-1);

  for ( size_t i = 0 ; i < static_cast<size_t>(nodes.size()) ; ++i ) nodes(i) = i+j;

  return nodes;
}

// ----------------------------- node-numbers along the node left edge -----------------------------

inline ColS FineLayer::nodesLeftEdge()
{
  return m_startNode;
}

// ------------------------------- node-numbers along the right edge -------------------------------

inline ColS FineLayer::nodesRightEdge()
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

// ----------------------------- node-number of the bottom-left corner -----------------------------

inline size_t FineLayer::nodesBottomLeftCorner()
{
  ColS nodes = nodesBottomEdge();

  return nodes(0);
}

// ---------------------------- node-number of the bottom-right corner -----------------------------

inline size_t FineLayer::nodesBottomRightCorner()
{
  ColS nodes = nodesBottomEdge();

  return nodes(nodes.size()-1);
}

// ------------------------------ node-number of the top-left corner -------------------------------

inline size_t FineLayer::nodesTopLeftCorner()
{
  ColS nodes = nodesTopEdge();

  return nodes(0);
}

// ------------------------------ node-number of the top-right corner ------------------------------

inline size_t FineLayer::nodesTopRightCorner()
{
  ColS nodes = nodesTopEdge();

  return nodes(nodes.size()-1);
}

// ----------------------------- node-number of the corners (aliases) ------------------------------

inline size_t FineLayer::nodesLeftBottomCorner()  { return nodesBottomLeftCorner();  }
inline size_t FineLayer::nodesLeftTopCorner()     { return nodesTopLeftCorner();     }
inline size_t FineLayer::nodesRightBottomCorner() { return nodesBottomRightCorner(); }
inline size_t FineLayer::nodesRightTopCorner()    { return nodesTopRightCorner();    }

// ------------------------------ node-numbers of periodic node-pairs ------------------------------

inline MatS FineLayer::nodesPeriodic()
{
  // edges
  ColS bot = nodesBottomEdge();
  ColS top = nodesTopEdge();
  ColS rgt = nodesRightEdge();
  ColS lft = nodesLeftEdge();

  // allocate nodal ties
  MatS nodes(bot.size()-2 + lft.size()-2 + 3, 2);

  // counter
  size_t i = 0;

  // tie all corners to one corner
  nodes(i,0) = nodesBottomLeftCorner(); nodes(i,1) = nodesBottomRightCorner(); ++i;
  nodes(i,0) = nodesBottomLeftCorner(); nodes(i,1) = nodesTopRightCorner();    ++i;
  nodes(i,0) = nodesBottomLeftCorner(); nodes(i,1) = nodesTopLeftCorner();     ++i;

  // tie edges to each other (exclude corners)
  for ( auto j = 1 ; j < bot.size()-1 ; ++j ) { nodes(i,0) = bot(j); nodes(i,1) = top(j); ++i; }
  for ( auto j = 1 ; j < lft.size()-1 ; ++j ) { nodes(i,0) = lft(j); nodes(i,1) = rgt(j); ++i; }

  return nodes;
}

// ------------------------------ node-number that lies in the origin ------------------------------

inline size_t FineLayer::nodesOrigin()
{
  return nodesBottomLeftCorner();
}

// ------------------------- DOF numbers per node (sequentially numbered) --------------------------

inline MatS FineLayer::dofs()
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// ------------------------ DOP-numbers with periodic dependencies removed -------------------------

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

}}} // namespace ...

// =================================================================================================

#endif



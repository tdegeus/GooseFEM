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

// ------------------------------------------ constructor ------------------------------------------

inline Regular::Regular(size_t nelx, size_t nely, double h):
m_h(h), m_nelx(nelx), m_nely(nely)
{
  assert( m_nelx >= 1 );
  assert( m_nely >= 1 );

  m_nnode = (m_nelx+1) * (m_nely+1);
  m_nelem =  m_nelx    *  m_nely   ;
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

// ------------------------ number of nodes, after eliminating periodicity -------------------------

inline size_t Regular::nnodePeriodic() const
{
  return (m_nelx+1) * (m_nely+1) - (m_nely+1) - (m_nelx);
}

// --------------------------------- coordinates (nodal positions) ---------------------------------

inline xt::xtensor<double,2> Regular::coor() const
{
  xt::xtensor<double,2> out = xt::empty<double>({m_nnode, m_ndim});

  xt::xtensor<double,1> x = xt::linspace<double>(0.0, m_h*static_cast<double>(m_nelx), m_nelx+1);
  xt::xtensor<double,1> y = xt::linspace<double>(0.0, m_h*static_cast<double>(m_nely), m_nely+1);

  size_t inode = 0;

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy ) {
    for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix ) {
      out(inode,0) = x(ix);
      out(inode,1) = y(iy);
      ++inode;
    }
  }

  return out;
}

// ---------------------------- connectivity (node-numbers per element) ----------------------------

inline xt::xtensor<size_t,2> Regular::conn() const
{
  xt::xtensor<size_t,2> out = xt::empty<size_t>({m_nelem,m_nne});

  size_t ielem = 0;

  for ( size_t iy = 0 ; iy < m_nely ; ++iy ) {
    for ( size_t ix = 0 ; ix < m_nelx ; ++ix ) {
      out(ielem,0) = (iy  )*(m_nelx+1) + (ix  );
      out(ielem,1) = (iy  )*(m_nelx+1) + (ix+1);
      out(ielem,3) = (iy+1)*(m_nelx+1) + (ix  );
      out(ielem,2) = (iy+1)*(m_nelx+1) + (ix+1);
      ++ielem;
    }
  }

  return out;
}

// ------------------------------ node-numbers along the bottom edge -------------------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx+1});

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = ix;

  return out;
}

// -------------------------------- node-numbers along the top edge --------------------------------

inline xt::xtensor<size_t,1> Regular::nodesTopEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx+1});

  for ( size_t ix = 0 ; ix < m_nelx+1 ; ++ix )
    out(ix) = ix + m_nely*(m_nelx+1);

  return out;
}

// ------------------------------- node-numbers along the left edge --------------------------------

inline xt::xtensor<size_t,1> Regular::nodesLeftEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely+1});

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1);

  return out;
}

// ------------------------------- node-numbers along the right edge -------------------------------

inline xt::xtensor<size_t,1> Regular::nodesRightEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely+1});

  for ( size_t iy = 0 ; iy < m_nely+1 ; ++iy )
    out(iy) = iy*(m_nelx+1) + m_nelx;

  return out;
}

// ---------------------- node-numbers along the bottom edge, without corners ----------------------

inline xt::xtensor<size_t,1> Regular::nodesBottomOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx-1});

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = ix;

  return out;
}

// ----------------------- node-numbers along the top edge, without corners ------------------------

inline xt::xtensor<size_t,1> Regular::nodesTopOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx-1});

  for ( size_t ix = 1 ; ix < m_nelx ; ++ix )
    out(ix-1) = ix + m_nely*(m_nelx+1);

  return out;
}

// ----------------------- node-numbers along the left edge, without corners -----------------------

inline xt::xtensor<size_t,1> Regular::nodesLeftOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely-1});

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1);

  return out;
}

// ---------------------- node-numbers along the right edge, without corners -----------------------

inline xt::xtensor<size_t,1> Regular::nodesRightOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nely-1});

  for ( size_t iy = 1 ; iy < m_nely ; ++iy )
    out(iy-1) = iy*(m_nelx+1) + m_nelx;

  return out;
}

// ----------------------------- node-number of the bottom-left corner -----------------------------

inline size_t Regular::nodesBottomLeftCorner() const
{
  return 0;
}

// ---------------------------- node-number of the bottom-right corner -----------------------------

inline size_t Regular::nodesBottomRightCorner() const
{
  return m_nelx;
}

// ------------------------------ node-number of the top-left corner -------------------------------

inline size_t Regular::nodesTopLeftCorner() const
{
  return m_nely*(m_nelx+1);
}

// ------------------------------ node-number of the top-right corner ------------------------------

inline size_t Regular::nodesTopRightCorner() const
{
  return m_nely*(m_nelx+1) + m_nelx;
}

// ----------------------------- node-number of the corners (aliases) ------------------------------

inline size_t Regular::nodesLeftBottomCorner() const  { return nodesBottomLeftCorner();  }
inline size_t Regular::nodesLeftTopCorner() const     { return nodesTopLeftCorner();     }
inline size_t Regular::nodesRightBottomCorner() const { return nodesBottomRightCorner(); }
inline size_t Regular::nodesRightTopCorner() const    { return nodesTopRightCorner();    }

// ------------------------------ node-numbers of periodic node-pairs ------------------------------

inline xt::xtensor<size_t,2> Regular::nodesPeriodic() const
{
  // edges (without corners)
  xt::xtensor<size_t,1> bot = nodesBottomOpenEdge();
  xt::xtensor<size_t,1> top = nodesTopOpenEdge();
  xt::xtensor<size_t,1> lft = nodesLeftOpenEdge();
  xt::xtensor<size_t,1> rgt = nodesRightOpenEdge();

  // allocate nodal ties
  // - number of tying per category
  size_t tedge = bot.size() + lft.size();
  size_t tnode = 3;
  // - allocate
  xt::xtensor<size_t,2> out = xt::empty<size_t>({tedge+tnode, std::size_t(2)});

  // counter
  size_t i = 0;

  // tie all corners to one corner
  out(i,0) = nodesBottomLeftCorner(); out(i,1) = nodesBottomRightCorner(); ++i;
  out(i,0) = nodesBottomLeftCorner(); out(i,1) = nodesTopRightCorner();    ++i;
  out(i,0) = nodesBottomLeftCorner(); out(i,1) = nodesTopLeftCorner();     ++i;

  // tie all corresponding edges to each other
  for ( size_t j = 0 ; j<bot.size() ; ++j ){ out(i,0) = bot(j); out(i,1) = top(j); ++i; }
  for ( size_t j = 0 ; j<lft.size() ; ++j ){ out(i,0) = lft(j); out(i,1) = rgt(j); ++i; }

  return out;
}

// ------------------------------ node-number that lies in the origin ------------------------------

inline size_t Regular::nodesOrigin() const
{
  return nodesBottomLeftCorner();
}

// ------------------------- DOF numbers per node (sequentially numbered) --------------------------

inline xt::xtensor<size_t,2> Regular::dofs() const
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// ------------------------ DOP-numbers with periodic dependencies removed -------------------------

inline xt::xtensor<size_t,2> Regular::dofsPeriodic() const
{
  // DOF-numbers for each component of each node (sequential)
  xt::xtensor<size_t,2> out = GooseFEM::Mesh::dofs(m_nnode,m_ndim);

  // periodic node-pairs
  xt::xtensor<size_t,2> nodePer = nodesPeriodic();

  // eliminate 'dependent' DOFs; renumber "out" to be sequential for the remaining DOFs
  for ( size_t i = 0 ; i < nodePer.shape()[0] ; ++i )
    for ( size_t j = 0 ; j < m_ndim ; ++j )
      out(nodePer(i,1),j) = out(nodePer(i,0),j);

  // renumber "out" to be sequential
  return GooseFEM::Mesh::renumber(out);
}

// ------------------------------------------ constructor ------------------------------------------

inline FineLayer::FineLayer(size_t nelx, size_t nely, double h, size_t nfine):
m_h(h)
{
  // basic assumptions
  assert( nelx >= 1 );
  assert( nely >= 1 );

  // store basic info
  m_Lx = m_h * static_cast<double>(nelx);

  // compute element size in y-direction (use symmetry, compute upper half)
  // ----------------------------------------------------------------------

  // temporary variables
  size_t nmin, ntot;
  xt::xtensor<size_t,1> nhx    =      xt::ones<size_t>({nely});
  xt::xtensor<size_t,1> nhy    =      xt::ones<size_t>({nely});
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
  // (2) element size in x-direction should fit the total number of elements in x-direction
  // (3) a certain number of layers have the minimum size "1" (are fine)
  for ( size_t iy = nfine ; ; )
  {
    // initialize current size in y-direction
    if ( iy == nfine ) ntot = nfine;
    // check to stop
    if ( iy >= nely or ntot >= nmin ) { nely = iy; break; }

    // rules (1,2) satisfied: coarsen in x-direction
    if ( 3*nhy(iy) <= ntot and nelx%(3*nhx(iy)) == 0 and ntot+nhy(iy) < nmin )
    {
      refine(iy)  = 0;
      nhy   (iy) *= 2;
      auto vnhy = xt::view(nhy, xt::range(iy+1, _));
      auto vnhx = xt::view(nhx, xt::range(iy  , _));
      vnhy *= 3;
      vnhx *= 3;
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
  m_nhx       = xt::empty<size_t>({nely*2-1});
  m_nhy       = xt::empty<size_t>({nely*2-1});
  m_refine    = xt::empty<int>   ({nely*2-1});
  m_nelx      = xt::empty<size_t>({nely*2-1});
  m_nnd       = xt::empty<size_t>({nely*2  });
  m_startElem = xt::empty<size_t>({nely*2-1});
  m_startNode = xt::empty<size_t>({nely*2  });

  // fill
  // - lower half
  for ( size_t iy = 0 ; iy < nely ; ++iy )
  {
    m_nhx   (iy) = nhx   (nely-iy-1);
    m_nhy   (iy) = nhy   (nely-iy-1);
    m_refine(iy) = refine(nely-iy-1);
  }
  // - upper half
  for ( size_t iy = 0 ; iy < nely-1 ; ++iy )
  {
    m_nhx   (iy+nely) = nhx   (iy+1);
    m_nhy   (iy+nely) = nhy   (iy+1);
    m_refine(iy+nely) = refine(iy+1);
  }

  // update size
  nely = m_nhx.size();

  // compute the number of elements per element layer in y-direction
  for ( size_t iy = 0 ; iy < nely ; ++iy )
    m_nelx(iy) = nelx / m_nhx(iy);

  // compute the number of nodes per node layer in y-direction
  for ( size_t iy = 0          ; iy < (nely+1)/2 ; ++iy ) m_nnd(iy  ) = m_nelx(iy)+1;
  for ( size_t iy = (nely-1)/2 ; iy < nely       ; ++iy ) m_nnd(iy+1) = m_nelx(iy)+1;

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
    if ( m_refine(i) == 0 ) { m_nnode += (3*m_nelx(i)+1); }
    else                    { m_nnode += (  m_nelx(i)+1); }
    // - add the elements of this layer
    if ( m_refine(i) == 0 ) { m_nelem += (4*m_nelx(i)  ); }
    else                    { m_nelem += (  m_nelx(i)  ); }
    // - store the starting node of the next layer
    m_startNode(i+1) = m_nnode;
  }

  // loop over element layers (middle -> top, elements become coarser)
  for ( size_t i = (nely-1)/2 ; i < nely ; ++i )
  {
    // - store the first element of the layer
    m_startElem(i) = m_nelem;
    // - add the nodes of this layer
    if ( m_refine(i) == 0 ) { m_nnode += (5*m_nelx(i)+1); }
    else                    { m_nnode += (  m_nelx(i)+1); }
    // - add the elements of this layer
    if ( m_refine(i) == 0 ) { m_nelem += (4*m_nelx(i)  ); }
    else                    { m_nelem += (  m_nelx(i)  ); }
    // - store the starting node of the next layer
    m_startNode(i+1) = m_nnode;
  }
  // - add the top row of nodes
  m_nnode += m_nelx(nely-1)+1;
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
  assert( i >= 0 and i <= 1 );

  if ( i == 0 ) return xt::amax(m_nelx)[0];
  else          return xt::sum (m_nhy )[0];

}

// --------------------------------- coordinates (nodal positions) ---------------------------------

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
  // --------------------------------------------------------------------------------------------

  for ( size_t iy = 0 ; ; ++iy )
  {
    // get positions along the x- and z-axis
    xt::xtensor<double,1> x = xt::linspace<double>(0.0, m_Lx, m_nelx(iy)+1);

    // add nodes of the bottom layer of this element
    for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
      out(inode,0) = x(ix);
      out(inode,1) = y(iy);
      ++inode;
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
      for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
        for ( size_t j = 0 ; j < 2 ; ++j ) {
          out(inode,0) = x(ix) + dx * static_cast<double>(j+1);
          out(inode,1) = y(iy) + dy;
          ++inode;
        }
      }
    }
  }

  // loop over element layers (middle -> top) : add (refinement layer +) top layer of nodes
  // --------------------------------------------------------------------------------------

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
  {
    // get positions along the x- and z-axis
    xt::xtensor<double,1> x = xt::linspace<double>(0.0, m_Lx, m_nelx(iy)+1);

    // add extra nodes of the intermediate layer, for refinement in x-direction
    if ( m_refine(iy) == 0 )
    {
      // - get position offset in x- and y-direction
      double dx = m_h * static_cast<double>(m_nhx(iy)/3);
      double dy = m_h * static_cast<double>(m_nhy(iy)/2);
      // - add nodes of the intermediate layer
      for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
        for ( size_t j = 0 ; j < 2 ; ++j ) {
          out(inode,0) = x(ix) + dx * static_cast<double>(j+1);
          out(inode,1) = y(iy) + dy;
          ++inode;
        }
      }
    }

    // add nodes of the top layer of this element
    for ( size_t ix = 0 ; ix < m_nelx(iy)+1 ; ++ix ) {
      out(inode,0) = x(ix  );
      out(inode,1) = y(iy+1);
      ++inode;
    }
  }

  return out;
}

// ---------------------------- connectivity (node-numbers per element) ----------------------------

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
      for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
        out(ielem,0) = bot + (ix  );
        out(ielem,1) = bot + (ix+1);
        out(ielem,2) = top + (ix+1);
        out(ielem,3) = top + (ix  );
        ielem++;
      }
    }

    // - define connectivity: refinement along the x-direction (below the middle layer)
    else if ( m_refine(iy) == 0 and iy <= (nely-1)/2 )
    {
      for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
        // -- bottom element
        out(ielem,0) = bot + (  ix  );
        out(ielem,1) = bot + (  ix+1);
        out(ielem,2) = mid + (2*ix+1);
        out(ielem,3) = mid + (2*ix  );
        ielem++;
        // -- top-right element
        out(ielem,0) = bot + (  ix+1);
        out(ielem,1) = top + (3*ix+3);
        out(ielem,2) = top + (3*ix+2);
        out(ielem,3) = mid + (2*ix+1);
        ielem++;
        // -- top-center element
        out(ielem,0) = mid + (2*ix  );
        out(ielem,1) = mid + (2*ix+1);
        out(ielem,2) = top + (3*ix+2);
        out(ielem,3) = top + (3*ix+1);
        ielem++;
        // -- top-left element
        out(ielem,0) = bot + (  ix  );
        out(ielem,1) = mid + (2*ix  );
        out(ielem,2) = top + (3*ix+1);
        out(ielem,3) = top + (3*ix  );
        ielem++;
      }
    }

    // - define connectivity: coarsening along the x-direction (above the middle layer)
    else if ( m_refine(iy) == 0 and iy > (nely-1)/2 )
    {
      for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix ) {
        // -- lower-left element
        out(ielem,0) = bot + (3*ix  );
        out(ielem,1) = bot + (3*ix+1);
        out(ielem,2) = mid + (2*ix  );
        out(ielem,3) = top + (  ix  );
        ielem++;
        // -- lower-center element
        out(ielem,0) = bot + (3*ix+1);
        out(ielem,1) = bot + (3*ix+2);
        out(ielem,2) = mid + (2*ix+1);
        out(ielem,3) = mid + (2*ix  );
        ielem++;
        // -- lower-right element
        out(ielem,0) = bot + (3*ix+2);
        out(ielem,1) = bot + (3*ix+3);
        out(ielem,2) = top + (  ix+1);
        out(ielem,3) = mid + (2*ix+1);
        ielem++;
        // -- upper element
        out(ielem,0) = mid + (2*ix  );
        out(ielem,1) = mid + (2*ix+1);
        out(ielem,2) = top + (  ix+1);
        out(ielem,3) = top + (  ix  );
        ielem++;
      }
    }
  }

  return out;
}

// ------------------------------ element numbers of the middle layer ------------------------------

inline xt::xtensor<size_t,1> FineLayer::elementsMiddleLayer() const
{
  // number of element layers in y-direction, the index of the middle layer
  size_t nely = static_cast<size_t>(m_nhy.size());
  size_t iy   = (nely-1)/2;

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(iy)});

  for ( size_t ix = 0 ; ix < m_nelx(iy) ; ++ix )
    out(ix) = m_startElem(iy) + ix;

  return out;
}

// ------------------------------ node-numbers along the bottom edge -------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(0)+1});

  for ( size_t ix = 0 ; ix < m_nelx(0)+1 ; ++ix )
    out(ix) = m_startNode(0) + ix;

  return out;
}

// -------------------------------- node-numbers along the top edge --------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesTopEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(nely-1)+1});

  for ( size_t ix = 0 ; ix < m_nelx(nely-1)+1 ; ++ix )
    out(ix) = m_startNode(nely) + ix;

  return out;
}

// ------------------------------- node-numbers along the left edge --------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesLeftEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely+1});

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1);

  return out;
}

// ------------------------------- node-numbers along the right edge -------------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesRightEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely+1});

  for ( size_t iy = 0 ; iy < (nely+1)/2 ; ++iy )
    out(iy) = m_startNode(iy) + m_nelx(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely ; ++iy )
    out(iy+1) = m_startNode(iy+1) + m_nelx(iy);

  return out;
}

// ---------------------- node-numbers along the bottom edge, without corners ----------------------

inline xt::xtensor<size_t,1> FineLayer::nodesBottomOpenEdge() const
{
  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(0)-1});

  for ( size_t ix = 1 ; ix < m_nelx(0) ; ++ix )
    out(ix-1) = m_startNode(0) + ix;

  return out;
}

// ----------------------- node-numbers along the top edge, without corners ------------------------

inline xt::xtensor<size_t,1> FineLayer::nodesTopOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({m_nelx(nely-1)-1});

  for ( size_t ix = 1 ; ix < m_nelx(nely-1) ; ++ix )
    out(ix-1) = m_startNode(nely) + ix;

  return out;
}

// ----------------------- node-numbers along the left edge, without corners -----------------------

inline xt::xtensor<size_t,1> FineLayer::nodesLeftOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely-1});

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1);

  return out;
}

// ---------------------- node-numbers along the right edge, without corners -----------------------

inline xt::xtensor<size_t,1> FineLayer::nodesRightOpenEdge() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  xt::xtensor<size_t,1> out = xt::empty<size_t>({nely-1});

  for ( size_t iy = 1 ; iy < (nely+1)/2 ; ++iy )
    out(iy-1) = m_startNode(iy) + m_nelx(iy);

  for ( size_t iy = (nely-1)/2 ; iy < nely-1 ; ++iy )
    out(iy) = m_startNode(iy+1) + m_nelx(iy);

  return out;
}

// ----------------------------- node-number of the bottom-left corner -----------------------------

inline size_t FineLayer::nodesBottomLeftCorner() const
{
  return m_startNode(0);
}

// ---------------------------- node-number of the bottom-right corner -----------------------------

inline size_t FineLayer::nodesBottomRightCorner() const
{
  return m_startNode(0) + m_nelx(0);
}

// ------------------------------ node-number of the top-left corner -------------------------------

inline size_t FineLayer::nodesTopLeftCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely);
}

// ------------------------------ node-number of the top-right corner ------------------------------

inline size_t FineLayer::nodesTopRightCorner() const
{
  size_t nely = static_cast<size_t>(m_nhy.size());

  return m_startNode(nely) + m_nelx(nely-1);
}

// -------------------------------------------- aliases --------------------------------------------

inline size_t FineLayer::nodesLeftBottomCorner() const  { return nodesBottomLeftCorner();  }
inline size_t FineLayer::nodesRightBottomCorner() const { return nodesBottomRightCorner(); }
inline size_t FineLayer::nodesLeftTopCorner() const     { return nodesTopLeftCorner();     }
inline size_t FineLayer::nodesRightTopCorner() const    { return nodesTopRightCorner();    }

// ------------------------------ node-numbers of periodic node-pairs ------------------------------

inline xt::xtensor<size_t,2> FineLayer::nodesPeriodic() const
{
  // edges (without corners)
  xt::xtensor<size_t,1> bot = nodesBottomOpenEdge();
  xt::xtensor<size_t,1> top = nodesTopOpenEdge();
  xt::xtensor<size_t,1> lft = nodesLeftOpenEdge();
  xt::xtensor<size_t,1> rgt = nodesRightOpenEdge();

  // allocate nodal ties
  // - number of tying per category
  size_t tedge = bot.size() + lft.size();
  size_t tnode = 3;
  // - allocate
  xt::xtensor<size_t,2> out = xt::empty<size_t>({tedge+tnode, std::size_t(2)});

  // counter
  size_t i = 0;

  // tie all corners to one corner
  out(i,0) = nodesBottomLeftCorner(); out(i,1) = nodesBottomRightCorner(); ++i;
  out(i,0) = nodesBottomLeftCorner(); out(i,1) = nodesTopRightCorner();    ++i;
  out(i,0) = nodesBottomLeftCorner(); out(i,1) = nodesTopLeftCorner();     ++i;

  // tie all corresponding edges to each other
  for ( size_t j = 0 ; j<bot.size() ; ++j ){ out(i,0) = bot(j); out(i,1) = top(j); ++i; }
  for ( size_t j = 0 ; j<lft.size() ; ++j ){ out(i,0) = lft(j); out(i,1) = rgt(j); ++i; }

  return out;
}

// ------------------------------ node-number that lies in the origin ------------------------------

inline size_t FineLayer::nodesOrigin() const
{
  return nodesBottomLeftCorner();
}

// ------------------------- DOF numbers per node (sequentially numbered) --------------------------

inline xt::xtensor<size_t,2> FineLayer::dofs() const
{
  return GooseFEM::Mesh::dofs(m_nnode,m_ndim);
}

// ------------------------ DOP-numbers with periodic dependencies removed -------------------------

inline xt::xtensor<size_t,2> FineLayer::dofsPeriodic() const
{
  // DOF-numbers for each component of each node (sequential)
  xt::xtensor<size_t,2> out = GooseFEM::Mesh::dofs(m_nnode,m_ndim);

  // periodic node-pairs
  xt::xtensor<size_t,2>   nodePer = nodesPeriodic();

  // eliminate 'dependent' DOFs; renumber "out" to be sequential for the remaining DOFs
  for ( size_t i = 0 ; i < nodePer.shape()[0] ; ++i )
    for ( size_t j = 0 ; j < m_ndim ; ++j )
      out(nodePer(i,1),j) = out(nodePer(i,0),j);

  // renumber "out" to be sequential
  return GooseFEM::Mesh::renumber(out);
}

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

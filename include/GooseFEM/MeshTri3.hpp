/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHTRI3_CPP
#define GOOSEFEM_MESHTRI3_CPP

// -------------------------------------------------------------------------------------------------

#include "MeshTri3.h"

// ===================================== GooseFEM::Mesh::Tri3 ======================================

namespace GooseFEM {
namespace Mesh {
namespace Tri3 {

// ------------------------------------------ constructor ------------------------------------------

inline Regular::Regular(size_t nelx, size_t nely, double h):
m_h(h), m_nelx(nelx), m_nely(nely)
{
  assert( m_nelx >= 1 );
  assert( m_nely >= 1 );

  m_nnode = (m_nelx+1) * (m_nely+1);
  m_nelem =  m_nelx    *  m_nely * 2;
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
      out(ielem,2) = (iy+1)*(m_nelx+1) + (ix  );
      ++ielem;
      out(ielem,0) = (iy  )*(m_nelx+1) + (ix+1);
      out(ielem,1) = (iy+1)*(m_nelx+1) + (ix+1);
      out(ielem,2) = (iy+1)*(m_nelx+1) + (ix  );
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

// ------------------------------ get the orientation of each element ------------------------------

inline xt::xtensor<int,1> getOrientation(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn)
{
  assert( conn.shape()[1] == 3 );
  assert( coor.shape()[1] == 2 );

  double k;
  size_t nelem = conn.shape()[0];

  xt::xtensor<int,1> out = xt::empty<int>({nelem});

  for ( size_t ielem = 0 ; ielem < nelem ; ++ielem )
  {
    auto v1 = xt::view(coor, conn(ielem,0), xt::all()) - xt::view(coor, conn(ielem,1), xt::all());
    auto v2 = xt::view(coor, conn(ielem,2), xt::all()) - xt::view(coor, conn(ielem,1), xt::all());

    k = v1(0) * v2(1) - v2(0) * v1(1);

    if ( k < 0 ) out(ielem) = -1;
    else         out(ielem) = +1;
  }

  return out;
}

// ------------------------------ set the orientation of each element ------------------------------

inline xt::xtensor<size_t,2> setOrientation(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn, int orientation)
{
  assert( conn.shape()[1] == 3 );
  assert( coor.shape()[1] == 2 );
  assert( orientation == -1 || orientation == +1 );

  xt::xtensor<int,1> val = getOrientation(coor, conn);

  return setOrientation(coor, conn, val, orientation);
}

// -------------------- set the orientation of each element to a certain value ---------------------

inline xt::xtensor<size_t,2> setOrientation(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn, const xt::xtensor<int,1> &val, int orientation)
{
  assert( conn.shape()[1] == 3 );
  assert( coor.shape()[1] == 2 );
  assert( conn.shape()[0] == val.size() );
  assert( orientation == -1 || orientation == +1 );

  // avoid compiler warning
  UNUSED(coor);

  size_t nelem = conn.shape()[0];
  xt::xtensor<size_t,2> out = conn;

  for ( size_t ielem = 0 ; ielem < nelem ; ++ielem )
    if ( ( orientation == -1 and val(ielem) > 0 ) or ( orientation == +1 and val(ielem) < 0 ) )
      std::swap( out(ielem,2) , out(ielem,1) );

  return out;
}

// ------------------------------------- re-triangulation API --------------------------------------

inline xt::xtensor<size_t,2> retriangulate(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn, int orientation)
{
  // get the orientation of all elements
  xt::xtensor<int,1> dir = getOrientation(coor, conn);
  // check the orientation
  bool eq = static_cast<size_t>(std::abs(xt::sum(dir)[0])) == conn.shape()[0];

  // new connectivity
  xt::xtensor<size_t,2> out;

  // perform re-triangulation
  // - use "TriUpdate"
  if ( eq )
  {
    Private::TriUpdate obj(coor,conn);
    obj.eval();
    out = obj.conn();
  }
  // - using TriCompute
  else
  {
    throw std::runtime_error("Work-in-progress, has to be re-triangulated using 'TriCompute'");
  }

  return setOrientation(coor,out,orientation);
}

// ================================= GooseFEM::Mesh::Tri3::Private =================================

namespace Private {

// ------------------------------------------ constructor ------------------------------------------

inline TriUpdate::TriUpdate(const xt::xtensor<double,2> &coor, const xt::xtensor<size_t,2> &conn): m_conn(conn), m_coor(coor)
{
  assert( conn.shape()[1] == 3 );
  assert( coor.shape()[1] == 2 );

  // store shapes
  m_nnode = coor.shape()[0];
  m_ndim  = coor.shape()[1];
  m_nelem = conn.shape()[0];
  m_nne   = conn.shape()[1];

  // set default to out-of-bounds, to make clear that nothing happened yet
  m_elem = m_nelem * xt::ones<size_t>({2});
  m_node = m_nnode * xt::ones<size_t>({4});

  edge();
}

// -------------------------- compute neighbors per edge of all elements ---------------------------

inline void TriUpdate::edge()
{
  // signal that nothing has been set
  m_edge = m_nelem * xt::ones<size_t>({m_nelem , m_nne});

  std::vector<size_t> idx = {0,1,2}; // lists to convert connectivity -> edges
  std::vector<size_t> jdx = {1,2,0}; // lists to convert connectivity -> edges

  std::vector<Edge> edge;
  edge.reserve(m_nelem*idx.size());

  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t i = 0 ; i < m_nne ; ++i )
      edge.push_back( Edge( m_conn(e,idx[i]), m_conn(e,jdx[i]) , e , i , true ) );

  std::sort( edge.begin() , edge.end() , Edge_sort );

  for ( size_t i = 0 ; i < edge.size()-1 ; ++i )
  {
    if ( edge[i].n1 == edge[i+1].n1 and edge[i].n2 == edge[i+1].n2 )
    {
      m_edge( edge[i  ].elem , edge[i  ].edge ) = edge[i+1].elem;
      m_edge( edge[i+1].elem , edge[i+1].edge ) = edge[i  ].elem;
    }
  }
}

// ---------------------------- update edges around renumbered elements ----------------------------

inline void TriUpdate::chedge(size_t edge, size_t old_elem, size_t new_elem)
{
  size_t m;
  size_t neigh = m_edge(old_elem , edge);

  if ( neigh >= m_nelem ) return;

  for ( m = 0 ; m < m_nne ; ++m )
    if ( m_edge( neigh , m ) == old_elem )
      break;

  m_edge( neigh , m ) = new_elem;
}

// --------------------------------- re-triangulate the full mesh ----------------------------------

inline bool TriUpdate::eval()
{
  bool change = false;

  while ( increment() ) { change = true; }

  return change;
}

// ----------------------------------- one re-triangulation step -----------------------------------

inline bool TriUpdate::increment()
{
  size_t ielem,jelem,iedge,jedge;
  double phi1,phi2;

  xt::xtensor_fixed<size_t,xt::xshape<4>> c = xt::empty<size_t>({4});
  xt::xtensor_fixed<size_t,xt::xshape<4>> n = xt::empty<size_t>({4});

  // loop over all elements
  for ( ielem = 0 ; ielem < m_nelem ; ++ielem )
  {
    // loop over all edges
    for ( iedge = 0 ; iedge < m_nne ; ++iedge )
    {
      // only proceed if the edge is shared with another element
      if ( m_edge(ielem,iedge) >= m_nelem )
        continue;

      // read "jelem"
      jelem = m_edge(ielem,iedge);

      // find the edge involved for "jelem"
      for ( jedge=0; jedge<m_nne; ++jedge ) if ( m_edge(jelem,jedge) == ielem ) break;

      // convert to four static nodes
      // - read first three from "ielem"
      if      ( iedge==0 ) { c(0)=m_conn(ielem,2); c(1)=m_conn(ielem,0); c(2)=m_conn(ielem,1); }
      else if ( iedge==1 ) { c(0)=m_conn(ielem,0); c(1)=m_conn(ielem,1); c(2)=m_conn(ielem,2); }
      else if ( iedge==2 ) { c(0)=m_conn(ielem,1); c(1)=m_conn(ielem,2); c(2)=m_conn(ielem,0); }
      // - read last from "jelem"
      if      ( jedge==0 ) { c(3)=m_conn(jelem,2); }
      else if ( jedge==1 ) { c(3)=m_conn(jelem,0); }
      else if ( jedge==2 ) { c(3)=m_conn(jelem,1); }

      // construct edge vectors
      auto a1 = xt::view(m_coor, c(1), xt::all()) - xt::view(m_coor, c(0), xt::all());
      auto b1 = xt::view(m_coor, c(2), xt::all()) - xt::view(m_coor, c(0), xt::all());
      auto a2 = xt::view(m_coor, c(1), xt::all()) - xt::view(m_coor, c(3), xt::all());
      auto b2 = xt::view(m_coor, c(2), xt::all()) - xt::view(m_coor, c(3), xt::all());

      // compute angles of the relevant corners
      phi1 = std::acos((a1(0)*b1(0)+a1(1)*b1(1))/(std::pow((a1(0)*a1(0)+a1(1)*a1(1)),0.5)*std::pow((b1(0)*b1(0)+b1(1)*b1(1)),0.5)));
      phi2 = std::acos((a2(0)*b2(0)+a2(1)*b2(1))/(std::pow((a2(0)*a2(0)+a2(1)*a2(1)),0.5)*std::pow((b2(0)*b2(0)+b2(1)*b2(1)),0.5)));

      // update mesh if needed
      if ( phi1+phi2 > M_PI )
      {
        // update connectivity
        m_conn(ielem,0) = c(0); m_conn(ielem,1) = c(3); m_conn(ielem,2) = c(2);
        m_conn(jelem,0) = c(0); m_conn(jelem,1) = c(1); m_conn(jelem,2) = c(3);

        // change list with neighbors for the elements around (only two neighbors change)
        if      ( iedge==0 ) { chedge(2,ielem,jelem); }
        else if ( iedge==1 ) { chedge(0,ielem,jelem); }
        else if ( iedge==2 ) { chedge(1,ielem,jelem); }

        if      ( jedge==0 ) { chedge(2,jelem,ielem); }
        else if ( jedge==1 ) { chedge(0,jelem,ielem); }
        else if ( jedge==2 ) { chedge(1,jelem,ielem); }

        // convert to four static nodes
        if      ( iedge==0 ) { n(0)=m_edge(ielem,2); n(3)=m_edge(ielem,1); }
        else if ( iedge==1 ) { n(0)=m_edge(ielem,0); n(3)=m_edge(ielem,2); }
        else if ( iedge==2 ) { n(0)=m_edge(ielem,1); n(3)=m_edge(ielem,0); }

        if      ( jedge==0 ) { n(1)=m_edge(jelem,1); n(2)=m_edge(jelem,2); }
        else if ( jedge==1 ) { n(1)=m_edge(jelem,2); n(2)=m_edge(jelem,0); }
        else if ( jedge==2 ) { n(1)=m_edge(jelem,0); n(2)=m_edge(jelem,1); }

        // store the neighbors for the changed elements
        m_edge(ielem,0) = jelem; m_edge(jelem,0) = n(0) ;
        m_edge(ielem,1) = n(2) ; m_edge(jelem,1) = n(1) ;
        m_edge(ielem,2) = n(3) ; m_edge(jelem,2) = ielem;

        // store information for transfer algorithm
        m_node    = c;
        m_elem(0) = ielem;
        m_elem(1) = jelem;

        return true;
      }
    }
  }

  return false;
}

// ------------------------------------------ constructor ------------------------------------------

inline Edge::Edge(size_t i, size_t j, size_t el, size_t ed, bool sort):
n1(i), n2(j), elem(el), edge(ed)
{
  if ( sort && n1>n2 )
    std::swap(n1,n2);
}

// --------------------------------------- compare two edges ---------------------------------------

inline bool Edge_cmp(Edge a, Edge b)
{
  if ( a.n1 == b.n1 and a.n2 == b.n2 )
    return true;

  return false;
}

// ----------------------- sort edges by comparing the first and second node -----------------------

inline bool Edge_sort(Edge a, Edge b)
{
  if ( a.n1 < b.n1 or a.n2 < b.n2 )
    return true;

  return false;
}

// -------------------------------------------------------------------------------------------------

} // namespace Private

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

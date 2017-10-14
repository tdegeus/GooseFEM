/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESHTRI3_H
#define GOOSEFEM_MESHTRI3_H

#include "Macros.h"

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Mesh {
namespace Tri3 {

// ========================================== REGULAR MESH =========================================

class Regular
{
private:
  size_t m_nrow;   // number of 'pixel' rows
  size_t m_ncol;   // number of 'pixel' columns
  size_t m_nelem;  // number of elements
  size_t m_nnode;  // number of nodes
  size_t m_nne=3;  // number of nodes-per-element
  size_t m_ndim=2; // number of dimensions

public:
  Regular            (const Regular &) = default;
  Regular& operator= (const Regular &) = default;
  Regular(){};
  Regular(size_t nrow, size_t ncol); // create mesh with [nrow,ncol] 'pixels' (nrow*ncol*2 elements)

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
};

// ========================================= MESH ANALYSIS =========================================

// read / set the orientation (-1 / +1) of all triangles
ColI getOrientation ( const MatD &coor, const MatS &conn                                      );
MatS setOrientation ( const MatD &coor, const MatS &conn,                  int orientation=-1 );
MatS setOrientation ( const MatD &coor, const MatS &conn, const ColI &val, int orientation=-1 );

// ======================================= RE-TRIANGULATION ========================================

// simple interface to compute the full re-triangulation; it uses, depending on the input mesh:
// (1) the minimal evasive "TriUpdate"
// (2) the more rigorous "TriCompute"
MatS retriangulate ( const MatD &coor, const MatS &conn, int orientation=-1 );

// -------------------------------------------------------------------------------------------------

// minimal evasive re-triangulation which only flips edges of the existing connectivity
class TriUpdate
{
private:
  MatS   m_edge;  // the element that neighbors along each edge (m_nelem: no neighbor)
  MatS   m_conn;  // connectivity (updated)
  MatD   m_coor;  // nodal positions (does not change)
  size_t m_nelem; // #elements
  size_t m_nnode; // #nodes
  size_t m_nne;   // #nodes-per-element
  size_t m_ndim;  // #dimensions
  ColS   m_elem;  // the two elements involved in the last element change (see below)
  ColS   m_node;  // the four nodes   involved in the last element change (see below)
  // old: m_elem(0) = [ m_node(0) , m_node(1) , m_node(2) ]
  //      m_elem(1) = [ m_node(1) , m_node(3) , m_node(2) ]
  // new: m_elem(0) = [ m_node(0) , m_node(3) , m_node(2) ]
  //      m_elem(1) = [ m_node(0) , m_node(1) , m_node(3) ]

  // compute neighbors per edge of all elements
  void edge();

  // update edges around renumbered elements
  void chedge(size_t edge, size_t old_elem, size_t new_elem);

public:
  TriUpdate(){};
  TriUpdate(const MatD &coor, const MatS &conn);

  bool eval     ();                     // re-triangulate the full mesh (returns "true" if changed)
  bool increment();                     // one re-triangulation step    (returns "true" if changed)
  MatS conn     () { return m_conn; };  // return (new) connectivity
  MatS ch_elem  () { return m_elem; };  // return element involved in last element change
  MatS ch_node  () { return m_node; };  // return nodes   involved in last element change
};

// ================================== REGULAR MESH - SOURCE CODE ===================================

Regular::Regular(size_t nrow, size_t ncol): m_nrow(nrow), m_ncol(ncol)
{
  assert( m_nrow >= 1 );
  assert( m_ncol >= 1 );

  m_nnode = (m_nrow+1) * (m_ncol+1)   ;
  m_nelem =  m_nrow    *  m_ncol   * 2;
}

// -------------------------------------------------------------------------------------------------

MatD Regular::coor()
{
  MatD coor( m_nnode , m_ndim );

  ColD x = ColD::LinSpaced( m_ncol+1 , 0.0 , 1.0 );
  ColD y = ColD::LinSpaced( m_nrow+1 , 0.0 , 1.0 );

  size_t inode = 0;

  for ( size_t row = 0 ; row < m_nrow+1 ; ++row ) {
    for ( size_t col = 0 ; col < m_ncol+1 ; ++col ) {
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

  for ( size_t row = 0 ; row < m_nrow ; ++row ) {
    for ( size_t col = 0 ; col < m_ncol ; ++col ) {
      conn(ielem,0) = (row+0)*(m_ncol+1)+col+0;
      conn(ielem,1) = (row+0)*(m_ncol+1)+col+1;
      conn(ielem,2) = (row+1)*(m_ncol+1)+col+0;
      ++ielem;
      conn(ielem,0) = (row+0)*(m_ncol+1)+col+1;
      conn(ielem,1) = (row+1)*(m_ncol+1)+col+1;
      conn(ielem,2) = (row+1)*(m_ncol+1)+col+0;
      ++ielem;
    }
  }

  return conn;
}

// -------------------------------------------------------------------------------------------------

ColS Regular::nodesBottom()
{
  ColS nodes(m_ncol+1);

  for ( size_t col = 0 ; col < m_ncol+1 ; ++col ) nodes(col) = col;

  return nodes;
}

// -------------------------------------------------------------------------------------------------

ColS Regular::nodesTop()
{
  ColS nodes(m_ncol+1);

  for ( size_t col = 0 ; col < m_ncol+1 ; ++col ) nodes(col) = col+m_nrow*(m_ncol+1);

  return nodes;
}

// -------------------------------------------------------------------------------------------------

ColS Regular::nodesLeft()
{
  ColS nodes(m_nrow+1);

  for ( size_t row = 0 ; row < m_nrow+1 ; ++row ) nodes(row) = row*(m_ncol+1);

  return nodes;
}

// -------------------------------------------------------------------------------------------------

ColS Regular::nodesRight()
{
  ColS nodes(m_nrow+1);

  for ( size_t row = 0 ; row < m_nrow+1 ; ++row ) nodes(row) = row*(m_ncol+1)+m_ncol;

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

// ================================== MESH ANALYSIS - SOURCE CODE ==================================

ColI getOrientation ( const MatD &coor, const MatS &conn )
{
  assert( conn.cols() == 3 );
  assert( coor.cols() == 2 );

  Eigen::Vector2d v1,v2;
  double k;
  size_t nelem = static_cast<size_t>( conn.rows() );

  ColI out( nelem );

  for ( size_t ielem = 0 ; ielem < nelem ; ++ielem )
  {
    v1 = coor.row( conn(ielem,0) ) - coor.row( conn(ielem,1) );
    v2 = coor.row( conn(ielem,2) ) - coor.row( conn(ielem,1) );

    k  = v1(0) * v2(1) - v2(0) * v1(1);

    if ( k < 0 ) out( ielem ) = -1;
    else         out( ielem ) = +1;
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

MatS setOrientation ( const MatD &coor, const MatS &conn, int orientation )
{
  assert( conn.cols() == 3 );
  assert( coor.cols() == 2 );
  assert( orientation == -1 || orientation == +1 );

  ColI val = getOrientation( coor, conn );

  return setOrientation( coor, conn, val, orientation );
}

// -------------------------------------------------------------------------------------------------

MatS setOrientation ( const MatD &coor, const MatS &conn, const ColI &val, int orientation )
{
  assert( conn.cols() == 3 );
  assert( coor.cols() == 2 );
  assert( conn.rows() == val.size() );
  assert( orientation == -1 || orientation == +1 );

  size_t nelem = static_cast<size_t>( conn.rows() );
  MatS   out   = conn;

  for ( size_t ielem = 0 ; ielem < nelem ; ++ielem )
    if ( ( orientation == -1 and val(ielem) > 0 ) or ( orientation == +1 and val(ielem) < 0 ) )
      std::swap( out(ielem,2) , out(ielem,1) );

  return out;
}

// ================================ RE-TRIANGULATION - SOURCE CODE =================================

MatS retriangulate ( const MatD &coor, const MatS &conn, int orientation )
{
  // get the orientation of all elements
  ColI dir = getOrientation( coor, conn );
  // check the orientation
  bool eq  = std::abs( dir.sum() ) == conn.rows();

  // new connectivity
  MatS out;

  // perform re-triangulation
  // - use "TriUpdate"
  if ( eq )
  {
    TriUpdate obj(coor,conn);
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

// =================================================================================================

// support class to allow the storage of a list of edges
class Edge {
public:
  size_t n1  ; // node 1 (edge from node 1-2)
  size_t n2  ; // node 2 (edge from node 1-2)
  size_t elem; // element to which the edge belong
  size_t edge; // edge index within the element (e.g. edge==1 -> n1=conn(0,elem), n2=conn(1,elem))

  Edge(){}

  Edge(size_t i, size_t j, size_t el, size_t ed, bool sort=false): n1(i), n2(j), elem(el), edge(ed)
  {
    if ( sort && n1>n2 ) { std::swap(n1,n2); }
  }
};

// -------------------------------------------------------------------------------------------------

bool Edge_cmp( Edge a , Edge b )
{
  if ( a.n1 == b.n1 and a.n2 == b.n2 )
    return true;

  return false;
}

// -------------------------------------------------------------------------------------------------

bool Edge_sort( Edge a , Edge b )
{
  if ( a.n1 < b.n1 or a.n2 < b.n2 )
    return true;

  return false;
}

// =================================================================================================

TriUpdate::TriUpdate(const MatD &coor, const MatS &conn): m_conn(conn), m_coor(coor)
{
  assert( conn.cols() == 3 );
	assert( coor.cols() == 2 );

	// store shapes
  m_nnode = static_cast<size_t>( coor.rows() );
  m_ndim  = static_cast<size_t>( coor.cols() );
  m_nelem = static_cast<size_t>( conn.rows() );
  m_nne   = static_cast<size_t>( conn.cols() );

  // resize internal arrays
  m_elem.resize(2);
  m_node.resize(4);

  // set default to out-of-bounds, to make clear that nothing happened yet
  m_elem.setConstant(m_nelem);
  m_node.setConstant(m_nnode);

  edge();
}

// -------------------------------------------------------------------------------------------------

void TriUpdate::edge()
{
  m_edge.resize( m_nelem , m_nne );
  m_edge.setConstant( m_nelem ); // signal that nothing has been set

  std::vector<size_t> idx = {0,1,2}; // lists to convert connectivity -> edges
  std::vector<size_t> jdx = {1,2,0}; // lists to convert connectivity -> edges

  std::vector<Edge> edge;
  edge.reserve( m_nelem*idx.size() );

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

// -------------------------------------------------------------------------------------------------

void TriUpdate::chedge(size_t edge, size_t old_elem, size_t new_elem)
{
  size_t m;
  size_t neigh = m_edge( old_elem , edge );

  if ( neigh >= m_nelem ) return;

  for ( m = 0 ; m < m_nne ; ++m )
    if ( m_edge( neigh , m ) == old_elem )
      break;

  m_edge( neigh , m ) = new_elem;
}

// -------------------------------------------------------------------------------------------------

bool TriUpdate::increment()
{
  size_t ielem,jelem,iedge,jedge;
  double phi1,phi2;

  ColS c(4);
  ColS n(4);

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
      Eigen::Vector2d a1 = m_coor.row( c(1) ) - m_coor.row( c(0) );
      Eigen::Vector2d b1 = m_coor.row( c(2) ) - m_coor.row( c(0) );
      Eigen::Vector2d a2 = m_coor.row( c(1) ) - m_coor.row( c(3) );
      Eigen::Vector2d b2 = m_coor.row( c(2) ) - m_coor.row( c(3) );

      // compute angles of the relevant corners
      phi1 = std::acos(a1.dot(b1)/(std::pow(a1.dot(a1),0.5)*std::pow(b1.dot(b1),0.5)));
      phi2 = std::acos(a2.dot(b2)/(std::pow(a2.dot(a2),0.5)*std::pow(b2.dot(b2),0.5)));

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

// -------------------------------------------------------------------------------------------------

bool TriUpdate::eval()
{
  bool change = false;

  while ( increment() ) { change = true; }

  return change;
}

// =================================================================================================

} // namespace Tri3
} // namespace Mesh
} // namespace GooseFEM

#endif

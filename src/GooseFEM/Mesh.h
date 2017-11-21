/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESH_H
#define GOOSEFEM_MESH_H

#include "Macros.h"

// -------------------------------------------------------------------------------------------------

namespace GooseFEM {
namespace Mesh {

// ====================================== EXTRACT INFORMATION ======================================

// elements connected to each node:
// out[ : , 0   ] = number of elements connected to each node
// out[ j , i+1 ] = "i"th element connected to node "j"
inline MatS elem2node ( const MatS &conn );

// list with DOF-numbers in sequential order
inline MatS dofs ( size_t nnode , size_t ndim );

// renumber list of indices
// - renumber to lowest possible numbers (e.g. [0,3,4,2] -> [0,2,3,1])
inline MatS renumber ( const MatS &in );
// - renumber to begin [ idx , ... ] or end [ ... , idx ] (e.g. [0,1,3,2] -> [3,0,2,1]; with idx=0)
inline MatS renumber ( const MatS &in , const ColS &idx , std::string location="end" );

// ========================================== SOURCE CODE ==========================================

inline MatS elem2node ( const MatS &conn )
{
  size_t nnode = conn.maxCoeff() + 1;

  ColS N( nnode );

  N *= 0;

  for ( auto e = 0 ; e < conn.rows() ; ++e )
    for ( auto m = 0 ; m < conn.cols() ; ++m )
      N( conn(e,m) ) += 1;

  MatS out( nnode , N.maxCoeff()+1 );
  out *= 0;

  for ( auto e = 0 ; e < conn.rows() ; ++e ) {
    for ( auto m = 0 ; m < conn.cols() ; ++m ) {
      auto nd = conn(e,m);
      out( nd, 0         ) += 1;
      out( nd, out(nd,0) )  = e;
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline MatS dofs ( size_t nnode , size_t ndim )
{
  ColS dofs_vec = ColS::LinSpaced( nnode*ndim , 0 , nnode*ndim );

  Eigen::Map<MatS> dofs( dofs_vec.data() , nnode , ndim );

  return dofs;
}

// -------------------------------------------------------------------------------------------------

inline MatS renumber ( const MatS &in )
{
  // size parameters
  size_t N = in.size();          // number of items in "in" (and in "out")
  size_t M = in.maxCoeff() + 1;  // size of renumber-list

  // allocate output
  MatS out ( in.rows() , in.cols() );

  // allocate renumber-list
  ColS renum ( M );

  // set all items as being 'excluded' ( == 0 )
  renum.setZero();

  // selectively set items that are 'included' ( == 1 )
  for ( size_t i = 0 ; i < N ; ++i ) renum ( in(i) ) = 1;

  // convert to indices, step 1: correct such that first index will start at 0
  renum[0] -= 1;
  // convert to indices, step 2: cumulative sum
  for ( size_t i = 1 ; i < M ; ++i ) renum ( i ) += renum ( i-1 );

  // renumber DOFs
  for ( size_t i = 0 ; i < N ; ++i ) out ( i ) = renum ( in(i) );

  return out;
}

// -------------------------------------------------------------------------------------------------

inline MatS renumber ( const MatS &in , const ColS &idx , std::string loc )
{
  // copy input to make modifications
  ColS tmp = idx;
  MatS out = in;
  ColS dat = ColS::LinSpaced(in.maxCoeff()+1, 0, in.maxCoeff()+1);

  // store size
  size_t ndx = static_cast<size_t>( idx.size() );
  size_t N   = static_cast<size_t>( out.size() );

  // sort input
  std::sort ( tmp.data(), tmp.data()+tmp.size() );

  // find "jdx = setdiff ( dat , idx )"
  // - allocate
  ColS jdx ( dat.size() - ndx );
  // - compute
  std::set_difference (
    dat.data() , dat.data()+dat.size() ,
    tmp.data() , tmp.data()+tmp.size() ,
    jdx.data()
  );

  // store size
  size_t mdx = static_cast<size_t>( jdx.size() );

  // renumber-list
  // - allocate
  ColS renum ( dat.size() );
  // - fill
  if ( loc == "end" )
  {
    // -- order [ jdx , idx ]
    for ( size_t i = 0 ; i < mdx ; ++i ) renum ( jdx(i) ) = i;
    for ( size_t i = 0 ; i < ndx ; ++i ) renum ( idx(i) ) = i + mdx;
  }
  else if ( loc == "begin" )
  {
    // -- order [ jdx, idx ]
    for ( size_t i = 0 ; i < ndx ; ++i ) renum ( idx(i) ) = i;
    for ( size_t i = 0 ; i < mdx ; ++i ) renum ( jdx(i) ) = i  + ndx;
  }

  // renumber
  for ( size_t i = 0 ; i < N ; ++i ) out ( i ) = renum ( out (i) );

  return out;
}

// -------------------------------------------------------------------------------------------------

} // namespace Mesh
} // namespace GooseFEM

#endif

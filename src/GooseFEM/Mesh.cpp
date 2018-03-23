/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MESH_CPP
#define GOOSEFEM_MESH_CPP

// -------------------------------------------------------------------------------------------------

#include "Mesh.h"

// ======================================== GooseFEM::Mesh =========================================

namespace GooseFEM {
namespace Mesh {

// ========================================= LIST OF DOFS ==========================================

inline MatS dofs(size_t nnode, size_t ndim)
{
  ColS dofs_vec = ColS::LinSpaced(nnode*ndim, 0, nnode*ndim);

  Eigen::Map<MatS> dofs(dofs_vec.data(), nnode, ndim);

  return dofs;
}

// ================================== RENUMBER TO LOWEST POSSIBLE ==================================

// --------------------------------------------- core ----------------------------------------------

template<class InputIterator, class OutputIterator>
inline void renumber(
  const InputIterator first, const InputIterator last, const OutputIterator result
)
{
  // size of input and the the renumber list
  auto N = last - first;
  auto M = (*std::max_element(first, last)) + 1;

  // allocate presence list and list of new indices
  std::vector<size_t> inList(M, 0); // initialize all items as 'excluded' (== 0)
  std::vector<size_t> index (M);

  // selectively set items that are 'included' (== 1)
  for ( auto it = first; it != last; ++it ) inList[*it] = 1;

  // new indices: cumulative sum of presence list
  std::partial_sum(inList.begin(), inList.end(), index.begin());
  // make index start at 0
  for ( auto &i: index ) i--;

  // apply renumbering
  for ( auto i = 0; i < N; ++i ) *(result+i) = index[*(first+i)];
};

// ------------------------------------------- interface -------------------------------------------

inline MatS renumber(const MatS &in)
{
  MatS out(in.rows(), in.cols());

  renumber(in.data(), in.data()+in.size(), out.data());

  return out;
}

// ============================================ REORDER ============================================

// --------------------------------------------- core ----------------------------------------------

template<class InputIterator, class OutputIterator, class IndexIterator>
inline void reorder(
  const InputIterator first, const InputIterator last, const OutputIterator result,
  const IndexIterator first_index, const IndexIterator last_index, std::string location
)
{
  // size of input and the the renumber list
  auto N = last - first;
  auto M = (*std::max_element(first, last)) + 1;

  // allocate presence list and list of new indices
  std::vector<size_t> inList(M, 0); // initialize all items as 'excluded' (== 0)
  std::vector<size_t> index (M);

  // selectively set items that are 'included' (== 1)
  for ( auto it = first; it != last; ++it ) inList[*it] = 1;
  // check that all indices are actually included
  #ifndef NDEBUG
  for ( auto it = first_index; it != last_index; ++it ) assert( inList[*it] == 1 );
  #endif
  // remove indices whose new index will be fixed
  for ( auto it = first_index; it != last_index; ++it ) inList[*it] = 0;

  // new indices: cumulative sum of presence list
  std::partial_sum(inList.begin(), inList.end(), index.begin());
  // make index start at 0
  for ( auto &i: index ) i--;
  // optionally move to end
  if ( location == "begin" ) for ( auto &i: index ) i += (last_index - first_index);

  // manually set indices
  // - allocate offset
  size_t offset;
  // - set offset
  if ( location == "begin" ) offset = 0;
  else                       offset = N - (last_index - first_index);
  // - apply
  for ( auto it = first_index; it != last_index; ++it ) { index[*it] = offset; ++offset; }

  // apply renumbering
  for ( auto i = 0; i < N; ++i ) *(result+i) = index[*(first+i)];
};

// ------------------------------------------- interface -------------------------------------------

inline MatS reorder(const MatS &in, const ColS &idx, std::string loc)
{
  MatS out(in.rows(), in.cols());

  reorder(in.data(), in.data()+in.size(), out.data(), idx.data(), idx.data()+idx.size(), loc);

  return out;
}

// ============================ LIST OF ELEMENTS CONNECTED TO EACH NODE ============================

inline SpMatS elem2node(const MatS &conn)
{
  // get number of nodes
  size_t nnode = conn.maxCoeff() + 1;

  // number of elements connected to each node
  // - allocate
  ColS N(nnode), idx(nnode);
  // - initialize
  N  .setZero();
  idx.setZero();
  // - fill from connectivity
  for ( auto it = conn.data(); it != conn.data()+conn.size(); ++it ) N(*it) += 1;

  // triplet list, with elements per node
  // - type
  typedef Eigen::Triplet<size_t> T;
  // - allocate
  std::vector<T> triplets;
  // - predict size
  triplets.reserve(N.sum());
  // - fill
  for ( auto e = 0 ; e < conn.rows() ; ++e ) {
    for ( auto m = 0 ; m < conn.cols() ; ++m ) {
      auto nd = conn(e,m);
      triplets.push_back(T(nd, idx(nd), e));
      idx(nd)++;
    }
  }

  // spare matrix
  // - allocate
  SpMatS mat(nnode, N.maxCoeff());
  // - fill
  mat.setFromTriplets(triplets.begin(), triplets.end());

  return mat;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENTHEX8_CPP
#define GOOSEFEM_ELEMENTHEX8_CPP

// -------------------------------------------------------------------------------------------------

#include "ElementHex8.h"

// ==================================== GooseFEM::Element::Hex8 ====================================

namespace GooseFEM {
namespace Element {
namespace Hex8 {

// ======================================== tensor algebra =========================================

inline double det(const T2 &A)
{
  return ( A[0] * A[4] * A[8] +
           A[1] * A[5] * A[6] +
           A[2] * A[3] * A[7] ) -
         ( A[2] * A[4] * A[6] +
           A[1] * A[3] * A[8] +
           A[0] * A[5] * A[7] );
}

// -------------------------------------------------------------------------------------------------

inline T2 inv(const T2 &A)
{
  // compute determinant
  double D = det(A);

  // allocate result
  T2 C = xt::empty<double>({ndim,ndim});

  // compute inverse
  C[0] = (A[4]*A[8]-A[5]*A[7]) / D;
  C[1] = (A[2]*A[7]-A[1]*A[8]) / D;
  C[2] = (A[1]*A[5]-A[2]*A[4]) / D;
  C[3] = (A[5]*A[6]-A[3]*A[8]) / D;
  C[4] = (A[0]*A[8]-A[2]*A[6]) / D;
  C[5] = (A[2]*A[3]-A[0]*A[5]) / D;
  C[6] = (A[3]*A[7]-A[4]*A[6]) / D;
  C[7] = (A[1]*A[6]-A[0]*A[7]) / D;
  C[8] = (A[0]*A[4]-A[1]*A[3]) / D;

  return C;
}

// ================================ GooseFEM::Element::Hex8::Gauss =================================

namespace Gauss {

// --------------------------------- number of integration points ----------------------------------

inline size_t nip()
{
  return 8;
}

// ----------------------- integration point coordinates (local coordinates) -----------------------

inline xt::xtensor<double,2> xi()
{
  static const size_t nip  = 8;
  static const size_t ndim = 3;

  xt::xtensor<double,2> xi = xt::empty<double>({nip,ndim});

  xi(0,0) = -1./std::sqrt(3.);    xi(0,1) = -1./std::sqrt(3.);    xi(0,2) = -1./std::sqrt(3.);
  xi(1,0) = +1./std::sqrt(3.);    xi(1,1) = -1./std::sqrt(3.);    xi(1,2) = -1./std::sqrt(3.);
  xi(2,0) = +1./std::sqrt(3.);    xi(2,1) = +1./std::sqrt(3.);    xi(2,2) = -1./std::sqrt(3.);
  xi(3,0) = -1./std::sqrt(3.);    xi(3,1) = +1./std::sqrt(3.);    xi(3,2) = -1./std::sqrt(3.);
  xi(4,0) = -1./std::sqrt(3.);    xi(4,1) = -1./std::sqrt(3.);    xi(4,2) = +1./std::sqrt(3.);
  xi(5,0) = +1./std::sqrt(3.);    xi(5,1) = -1./std::sqrt(3.);    xi(5,2) = +1./std::sqrt(3.);
  xi(6,0) = +1./std::sqrt(3.);    xi(6,1) = +1./std::sqrt(3.);    xi(6,2) = +1./std::sqrt(3.);
  xi(7,0) = -1./std::sqrt(3.);    xi(7,1) = +1./std::sqrt(3.);    xi(7,2) = +1./std::sqrt(3.);

  return xi;
}

// ----------------------------------- integration point weights -----------------------------------

inline xt::xtensor<double,1> w()
{
  static const size_t nip = 8;

  xt::xtensor<double,1> w = xt::empty<double>({nip});

  w(0) = 1.;
  w(1) = 1.;
  w(2) = 1.;
  w(3) = 1.;
  w(4) = 1.;
  w(5) = 1.;
  w(6) = 1.;
  w(7) = 1.;

  return w;
}

// -------------------------------------------------------------------------------------------------

}

// ================================ GooseFEM::Element::Hex8::Nodal ================================

namespace Nodal {

// --------------------------------- number of integration points ----------------------------------

inline size_t nip()
{
  return 8;
}

// ----------------------- integration point coordinates (local coordinates) -----------------------

inline xt::xtensor<double,2> xi()
{
  static const size_t nip  = 8;
  static const size_t ndim = 3;

  xt::xtensor<double,2> xi = xt::empty<double>({nip,ndim});

  xi(0,0) = -1.;    xi(0,1) = -1.;    xi(0,2) = -1.;
  xi(1,0) = +1.;    xi(1,1) = -1.;    xi(1,2) = -1.;
  xi(2,0) = +1.;    xi(2,1) = +1.;    xi(2,2) = -1.;
  xi(3,0) = -1.;    xi(3,1) = +1.;    xi(3,2) = -1.;
  xi(4,0) = -1.;    xi(4,1) = -1.;    xi(4,2) = +1.;
  xi(5,0) = +1.;    xi(5,1) = -1.;    xi(5,2) = +1.;
  xi(6,0) = +1.;    xi(6,1) = +1.;    xi(6,2) = +1.;
  xi(7,0) = -1.;    xi(7,1) = +1.;    xi(7,2) = +1.;

  return xi;
}

// ----------------------------------- integration point weights -----------------------------------

inline xt::xtensor<double,1> w()
{
  static const size_t nip = 8;

  xt::xtensor<double,1> w = xt::empty<double>({nip});

  w(0) = 1.;
  w(1) = 1.;
  w(2) = 1.;
  w(3) = 1.;
  w(4) = 1.;
  w(5) = 1.;
  w(6) = 1.;
  w(7) = 1.;

  return w;
}

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

// ------------------------------------------ constructor ------------------------------------------

inline Quadrature::Quadrature(const xt::xtensor<double,3> &x) : m_x(x)
{
  assert( m_x.shape()[1] == m_nne  );
  assert( m_x.shape()[2] == m_ndim );

  // set integration scheme
  m_xi = Gauss::xi();
  m_w  = Gauss::w ();

  // extract number of elements
  m_nelem = m_x.shape()[0];
  m_nip   = m_w.size();

  // allocate arrays
  // - shape functions
  m_N    = xt::empty<double>({m_nip,m_nne});
  // - shape function gradients in local coordinates
  m_dNxi = xt::empty<double>({m_nip,m_nne,m_ndim});
  // - shape function gradients in global coordinates
  m_dNx  = xt::empty<double>({m_nelem,m_nip,m_nne,m_ndim});
  // - integration point volume
  m_vol  = xt::empty<double>({m_nelem,m_nip});

  // shape functions
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    m_N(k,0) = .125 * (1.-m_xi(k,0)) * (1.-m_xi(k,1)) * (1.-m_xi(k,2));
    m_N(k,1) = .125 * (1.+m_xi(k,0)) * (1.-m_xi(k,1)) * (1.-m_xi(k,2));
    m_N(k,2) = .125 * (1.+m_xi(k,0)) * (1.+m_xi(k,1)) * (1.-m_xi(k,2));
    m_N(k,3) = .125 * (1.-m_xi(k,0)) * (1.+m_xi(k,1)) * (1.-m_xi(k,2));
    m_N(k,4) = .125 * (1.-m_xi(k,0)) * (1.-m_xi(k,1)) * (1.+m_xi(k,2));
    m_N(k,5) = .125 * (1.+m_xi(k,0)) * (1.-m_xi(k,1)) * (1.+m_xi(k,2));
    m_N(k,6) = .125 * (1.+m_xi(k,0)) * (1.+m_xi(k,1)) * (1.+m_xi(k,2));
    m_N(k,7) = .125 * (1.-m_xi(k,0)) * (1.+m_xi(k,1)) * (1.+m_xi(k,2));
  }

  // shape function gradients in local coordinates
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    // - dN / dxi_0
    m_dNxi(k,0,0) = -.125*(1.-m_xi(k,1))*(1.-m_xi(k,2));
    m_dNxi(k,1,0) = +.125*(1.-m_xi(k,1))*(1.-m_xi(k,2));
    m_dNxi(k,2,0) = +.125*(1.+m_xi(k,1))*(1.-m_xi(k,2));
    m_dNxi(k,3,0) = -.125*(1.+m_xi(k,1))*(1.-m_xi(k,2));
    m_dNxi(k,4,0) = -.125*(1.-m_xi(k,1))*(1.+m_xi(k,2));
    m_dNxi(k,5,0) = +.125*(1.-m_xi(k,1))*(1.+m_xi(k,2));
    m_dNxi(k,6,0) = +.125*(1.+m_xi(k,1))*(1.+m_xi(k,2));
    m_dNxi(k,7,0) = -.125*(1.+m_xi(k,1))*(1.+m_xi(k,2));
    // - dN / dxi_1
    m_dNxi(k,0,1) = -.125*(1.-m_xi(k,0))*(1.-m_xi(k,2));
    m_dNxi(k,1,1) = -.125*(1.+m_xi(k,0))*(1.-m_xi(k,2));
    m_dNxi(k,2,1) = +.125*(1.+m_xi(k,0))*(1.-m_xi(k,2));
    m_dNxi(k,3,1) = +.125*(1.-m_xi(k,0))*(1.-m_xi(k,2));
    m_dNxi(k,4,1) = -.125*(1.-m_xi(k,0))*(1.+m_xi(k,2));
    m_dNxi(k,5,1) = -.125*(1.+m_xi(k,0))*(1.+m_xi(k,2));
    m_dNxi(k,6,1) = +.125*(1.+m_xi(k,0))*(1.+m_xi(k,2));
    m_dNxi(k,7,1) = +.125*(1.-m_xi(k,0))*(1.+m_xi(k,2));
    // - dN / dxi_2
    m_dNxi(k,0,2) = -.125*(1.-m_xi(k,0))*(1.-m_xi(k,1));
    m_dNxi(k,1,2) = -.125*(1.+m_xi(k,0))*(1.-m_xi(k,1));
    m_dNxi(k,2,2) = -.125*(1.+m_xi(k,0))*(1.+m_xi(k,1));
    m_dNxi(k,3,2) = -.125*(1.-m_xi(k,0))*(1.+m_xi(k,1));
    m_dNxi(k,4,2) = +.125*(1.-m_xi(k,0))*(1.-m_xi(k,1));
    m_dNxi(k,5,2) = +.125*(1.+m_xi(k,0))*(1.-m_xi(k,1));
    m_dNxi(k,6,2) = +.125*(1.+m_xi(k,0))*(1.+m_xi(k,1));
    m_dNxi(k,7,2) = +.125*(1.-m_xi(k,0))*(1.+m_xi(k,1));
  }

  // compute the shape function gradients, based on "x"
  compute_dN();
}

// ------------------------------------------ constructor ------------------------------------------

inline Quadrature::Quadrature(const xt::xtensor<double,3> &x, const xt::xtensor<double,2> &xi,
  const xt::xtensor<double,1> &w) : m_x(x), m_w(w), m_xi(xi)
{
  assert( m_x.shape()[1] == m_nne  );
  assert( m_x.shape()[2] == m_ndim );

  // extract shape
  m_nelem = m_x.shape()[0];
  m_nip   = m_w.size();

  assert( m_xi.shape()[0] == m_nip  );
  assert( m_xi.shape()[1] == m_ndim );
  assert( m_w .size()     == m_nip  );

  // allocate arrays
  // - shape functions
  m_N    = xt::empty<double>({m_nip,m_nne});
  // - shape function gradients in local coordinates
  m_dNxi = xt::empty<double>({m_nip,m_nne,m_ndim});
  // - shape function gradients in global coordinates
  m_dNx  = xt::empty<double>({m_nelem,m_nip,m_nne,m_ndim});
  // - integration point volume
  m_vol  = xt::empty<double>({m_nelem,m_nip});

  // shape functions
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    m_N(k,0) = .125 * (1.-m_xi(k,0)) * (1.-m_xi(k,1)) * (1.-m_xi(k,2));
    m_N(k,1) = .125 * (1.+m_xi(k,0)) * (1.-m_xi(k,1)) * (1.-m_xi(k,2));
    m_N(k,2) = .125 * (1.+m_xi(k,0)) * (1.+m_xi(k,1)) * (1.-m_xi(k,2));
    m_N(k,3) = .125 * (1.-m_xi(k,0)) * (1.+m_xi(k,1)) * (1.-m_xi(k,2));
    m_N(k,4) = .125 * (1.-m_xi(k,0)) * (1.-m_xi(k,1)) * (1.+m_xi(k,2));
    m_N(k,5) = .125 * (1.+m_xi(k,0)) * (1.-m_xi(k,1)) * (1.+m_xi(k,2));
    m_N(k,6) = .125 * (1.+m_xi(k,0)) * (1.+m_xi(k,1)) * (1.+m_xi(k,2));
    m_N(k,7) = .125 * (1.-m_xi(k,0)) * (1.+m_xi(k,1)) * (1.+m_xi(k,2));
  }

  // shape function gradients in local coordinates
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    // - dN / dxi_0
    m_dNxi(k,0,0) = -.125*(1.-m_xi(k,1))*(1.-m_xi(k,2));
    m_dNxi(k,1,0) = +.125*(1.-m_xi(k,1))*(1.-m_xi(k,2));
    m_dNxi(k,2,0) = +.125*(1.+m_xi(k,1))*(1.-m_xi(k,2));
    m_dNxi(k,3,0) = -.125*(1.+m_xi(k,1))*(1.-m_xi(k,2));
    m_dNxi(k,4,0) = -.125*(1.-m_xi(k,1))*(1.+m_xi(k,2));
    m_dNxi(k,5,0) = +.125*(1.-m_xi(k,1))*(1.+m_xi(k,2));
    m_dNxi(k,6,0) = +.125*(1.+m_xi(k,1))*(1.+m_xi(k,2));
    m_dNxi(k,7,0) = -.125*(1.+m_xi(k,1))*(1.+m_xi(k,2));
    // - dN / dxi_1
    m_dNxi(k,0,1) = -.125*(1.-m_xi(k,0))*(1.-m_xi(k,2));
    m_dNxi(k,1,1) = -.125*(1.+m_xi(k,0))*(1.-m_xi(k,2));
    m_dNxi(k,2,1) = +.125*(1.+m_xi(k,0))*(1.-m_xi(k,2));
    m_dNxi(k,3,1) = +.125*(1.-m_xi(k,0))*(1.-m_xi(k,2));
    m_dNxi(k,4,1) = -.125*(1.-m_xi(k,0))*(1.+m_xi(k,2));
    m_dNxi(k,5,1) = -.125*(1.+m_xi(k,0))*(1.+m_xi(k,2));
    m_dNxi(k,6,1) = +.125*(1.+m_xi(k,0))*(1.+m_xi(k,2));
    m_dNxi(k,7,1) = +.125*(1.-m_xi(k,0))*(1.+m_xi(k,2));
    // - dN / dxi_2
    m_dNxi(k,0,2) = -.125*(1.-m_xi(k,0))*(1.-m_xi(k,1));
    m_dNxi(k,1,2) = -.125*(1.+m_xi(k,0))*(1.-m_xi(k,1));
    m_dNxi(k,2,2) = -.125*(1.+m_xi(k,0))*(1.+m_xi(k,1));
    m_dNxi(k,3,2) = -.125*(1.-m_xi(k,0))*(1.+m_xi(k,1));
    m_dNxi(k,4,2) = +.125*(1.-m_xi(k,0))*(1.-m_xi(k,1));
    m_dNxi(k,5,2) = +.125*(1.+m_xi(k,0))*(1.-m_xi(k,1));
    m_dNxi(k,6,2) = +.125*(1.+m_xi(k,0))*(1.+m_xi(k,1));
    m_dNxi(k,7,2) = +.125*(1.-m_xi(k,0))*(1.+m_xi(k,1));
  }

  // compute the shape function gradients, based on "x"
  compute_dN();
}

// --------------------------- integration volume (per tensor-component) ---------------------------

inline void Quadrature::dV(xt::xtensor<double,2> &qscalar) const
{
  assert( qscalar.shape()[0] == m_nelem );
  assert( qscalar.shape()[1] == m_nip   );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t k = 0 ; k < m_nip ; ++k )
      qscalar(e,k) = m_vol(e,k);
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::dV(xt::xtensor<double,4> &qtensor) const
{
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nne   );
  assert( qtensor.shape()[2] >= m_ndim  );
  assert( qtensor.shape()[3] >= m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t k = 0 ; k < m_nip ; ++k )
      for ( size_t i = 0 ; i < qtensor.shape()[2] ; ++i )
        for ( size_t j = 0 ; j < qtensor.shape()[3] ; ++j )
          qtensor(e,k,i,j) = m_vol(e,k);
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Quadrature::dV() const
{
  xt::xtensor<double,2> out = xt::empty<double>({m_nelem, m_nip});

  this->dV(out);

  return out;
}

// -------------------------------------- number of elements ---------------------------------------

inline size_t Quadrature::nelem() const
{
  return m_nelem;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t Quadrature::nne() const
{
  return m_nne;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t Quadrature::ndim() const
{
  return m_ndim;
}

// --------------------------------- number of integration points ----------------------------------

inline size_t Quadrature::nip() const
{
  return m_nip;
}

// --------------------------------------- update positions ----------------------------------------

inline void Quadrature::update_x(const xt::xtensor<double,3> &x)
{
  assert( x.shape()[0] == m_nelem    );
  assert( x.shape()[1] == m_nne      );
  assert( x.shape()[2] == m_ndim     );
  assert( x.size()     == m_x.size() );

  // update positions
  m_x = x;

  // update the shape function gradients for the new "x"
  compute_dN();
}

// ------------------------ shape function gradients in global coordinates -------------------------

inline void Quadrature::compute_dN()
{
  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias nodal positions
    auto x = xt::view(m_x, e, xt::all(), xt::all());

    // loop over integration points
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto dNxi = xt::view(m_dNxi,    k, xt::all(), xt::all());
      auto dNx  = xt::view(m_dNx , e, k, xt::all(), xt::all());

      // - Jacobian
      T2 J = xt::zeros<double>({m_ndim,m_ndim});
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          for ( size_t j = 0 ; j < m_ndim ; ++j )
            J(i,j) += dNxi(m,i) * x(m,j);

      // - determinant and inverse of the Jacobian
      double Jdet = det(J);
      T2     Jinv = inv(J);

      // - shape function gradients wrt global coordinates
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          for ( size_t j = 0 ; j < m_ndim ; ++j )
            dNx(m,i) += Jinv(i,j) * dNxi(m,j);

      // - copy to matrix: integration point volume
      m_vol(e,k) = m_w(k) * Jdet;
    }
  }
}

// ------------------- dyadic product "qtensor(i,j) = dNdx(m,i) * elemvec(m,j)" --------------------

inline void Quadrature::gradN_vector(
  const xt::xtensor<double,3> &elemvec, xt::xtensor<double,4> &qtensor) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nne   );
  assert( qtensor.shape()[2] >= m_ndim  );
  assert( qtensor.shape()[3] >= m_ndim  );

  // zero-initialize output: matrix of tensors
  qtensor *= 0.0;

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias element vector (e.g. nodal displacements)
    auto u = xt::view(elemvec, e, xt::all(), xt::all());

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto dNx   = xt::view(m_dNx  , e, k, xt::all()          , xt::all()          );
      auto gradu = xt::view(qtensor, e, k, xt::range(0,m_ndim), xt::range(0,m_ndim));

      // - evaluate dyadic product
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          for ( size_t j = 0 ; j < m_ndim ; ++j )
            gradu(i,j) += dNx(m,i) * u(m,j);
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Quadrature::gradN_vector(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,4> qtensor = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  this->gradN_vector(elemvec, qtensor);

  return qtensor;
}

// ---------------------------------- transpose of "gradN_vector" ----------------------------------

inline void Quadrature::gradN_vector_T(
  const xt::xtensor<double,3> &elemvec, xt::xtensor<double,4> &qtensor) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nne   );
  assert( qtensor.shape()[2] >= m_ndim  );
  assert( qtensor.shape()[3] >= m_ndim  );

  // zero-initialize output: matrix of tensors
  qtensor *= 0.0;

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias element vector (e.g. nodal displacements)
    auto u = xt::view(elemvec, e, xt::all(), xt::all());

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto dNx   = xt::view(m_dNx  , e, k, xt::all()          , xt::all()          );
      auto gradu = xt::view(qtensor, e, k, xt::range(0,m_ndim), xt::range(0,m_ndim));

      // - evaluate dyadic product
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          for ( size_t j = 0 ; j < m_ndim ; ++j )
            gradu(j,i) += dNx(m,i) * u(m,j);
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Quadrature::gradN_vector_T(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,4> qtensor = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  this->gradN_vector_T(elemvec, qtensor);

  return qtensor;
}

// ------------------------------- symmetric part of "gradN_vector" --------------------------------

inline void Quadrature::symGradN_vector(
  const xt::xtensor<double,3> &elemvec, xt::xtensor<double,4> &qtensor) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nne   );
  assert( qtensor.shape()[2] >= m_ndim  );
  assert( qtensor.shape()[3] >= m_ndim  );

  // zero-initialize output: matrix of tensors
  qtensor *= 0.0;

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias element vector (e.g. nodal displacements)
    auto u = xt::view(elemvec, e, xt::all(), xt::all());

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias shape function gradients (local coordinates)
      auto dNx = xt::view(m_dNx  , e, k, xt::all()          , xt::all()          );
      auto eps = xt::view(qtensor, e, k, xt::range(0,m_ndim), xt::range(0,m_ndim));

      // - evaluate dyadic product
      for ( size_t m = 0 ; m < m_nne ; ++m ) {
        for ( size_t i = 0 ; i < m_ndim ; ++i ) {
          for ( size_t j = 0 ; j < m_ndim ; ++j ) {
            eps(i,j) += dNx(m,i) * u(m,j) / 2.;
            eps(j,i) += dNx(m,i) * u(m,j) / 2.;
          }
        }
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Quadrature::symGradN_vector(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,4> qtensor = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  this->symGradN_vector(elemvec, qtensor);

  return qtensor;
}

// ------- scalar product "elemmat(m*ndim+i,n*ndim+i) = N(m) * qscalar * N(n)"; for all "i" --------

inline void Quadrature::int_N_scalar_NT_dV(
  const xt::xtensor<double,2> &qscalar, xt::xtensor<double,3> &elemmat) const
{
  assert( qscalar.shape()[0] == m_nelem      );
  assert( qscalar.shape()[1] == m_nip        );
  assert( elemmat.shape()[0] == m_nelem      );
  assert( elemmat.shape()[1] == m_nne*m_ndim );
  assert( elemmat.shape()[2] == m_nne*m_ndim );

  // zero-initialize: matrix of matrices
  elemmat *= 0.0;

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias (e.g. mass matrix)
    auto M = xt::view(elemmat, e, xt::all(), xt::all());

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias shape functions
      auto N = xt::view(m_N, k, xt::all());

      // - alias
      double vol = m_vol  (e,k);  // integration point volume
      double rho = qscalar(e,k);  // integration point scalar (e.g. density)

      // - evaluate scalar product, for all dimensions, and assemble
      //   M(m*ndim+i,n*ndim+i) += N(m) * scalar * N(n) * dV
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t n = 0 ; n < m_nne ; ++n )
          for ( size_t i = 0 ; i < m_ndim ; ++i )
            M(m*m_ndim+i, n*m_ndim+i) += N(m) * rho * N(n) * vol;
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Quadrature::int_N_scalar_NT_dV(const xt::xtensor<double,2> &qscalar) const
{
  xt::xtensor<double,3> elemmat = xt::empty<double>({m_nelem, m_nne*m_ndim, m_nne*m_ndim});

  this->int_N_scalar_NT_dV(qscalar, elemmat);

  return elemmat;
}

// ------------ integral of dot product "elemvec(m,j) += dNdx(m,i) * qtensor(i,j) * dV" ------------

inline void Quadrature::int_gradN_dot_tensor2_dV(const xt::xtensor<double,4> &qtensor,
  xt::xtensor<double,3> &elemvec) const
{
  assert( qtensor.shape()[0] == m_nelem ); // number of elements
  assert( qtensor.shape()[1] == m_nip   ); // number of integration points
  assert( qtensor.shape()[2] >= m_ndim  ); // number of dimensions
  assert( qtensor.shape()[3] >= m_ndim  ); // number of dimensions
  assert( elemvec.shape()[0] == m_nelem ); // number of elements
  assert( elemvec.shape()[1] == m_nne   ); // number of nodes per element
  assert( elemvec.shape()[2] == m_ndim  ); // number of dimensions

  // zero-initialize output: matrix of vectors
  elemvec *= 0.0;

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias (e.g. nodal force)
    auto f = xt::view(elemvec, e, xt::all(), xt::all());

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto   dNx = xt::view(m_dNx  , e, k, xt::all(), xt::all());
      auto   sig = xt::view(qtensor, e, k, xt::range(0,m_ndim), xt::range(0,m_ndim));
      double vol = m_vol(e,k);

      // - evaluate dot product, and assemble
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          for ( size_t j = 0 ; j < m_ndim ; ++j )
            f(m,j) += dNx(m,i) * sig(i,j) * vol;
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Quadrature::int_gradN_dot_tensor2_dV(const xt::xtensor<double,4> &qtensor) const
{
  xt::xtensor<double,3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});

  this->int_gradN_dot_tensor2_dV(qtensor, elemvec);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

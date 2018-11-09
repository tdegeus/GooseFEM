/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENTHEX8_CPP
#define GOOSEFEM_ELEMENTHEX8_CPP

// -------------------------------------------------------------------------------------------------

#include "ElementHex8.h"

// =================================================================================================

namespace GooseFEM {
namespace Element {
namespace Hex8 {

// =================================================================================================

inline double inv(const T2 &A, T2 &Ainv)
{
  // compute determinant
  double det = ( A(0,0) * A(1,1) * A(2,2) +
                 A(0,1) * A(1,2) * A(2,0) +
                 A(0,2) * A(1,0) * A(2,1) ) -
               ( A(0,2) * A(1,1) * A(2,0) +
                 A(0,1) * A(1,0) * A(2,2) +
                 A(0,0) * A(1,2) * A(2,1) );

  // compute inverse
  Ainv(0,0) = (A(1,1)*A(2,2)-A(1,2)*A(2,1)) / det;
  Ainv(0,1) = (A(0,2)*A(2,1)-A(0,1)*A(2,2)) / det;
  Ainv(0,2) = (A(0,1)*A(1,2)-A(0,2)*A(1,1)) / det;
  Ainv(1,0) = (A(1,2)*A(2,0)-A(1,0)*A(2,2)) / det;
  Ainv(1,1) = (A(0,0)*A(2,2)-A(0,2)*A(2,0)) / det;
  Ainv(1,2) = (A(0,2)*A(1,0)-A(0,0)*A(1,2)) / det;
  Ainv(2,0) = (A(1,0)*A(2,1)-A(1,1)*A(2,0)) / det;
  Ainv(2,1) = (A(0,1)*A(2,0)-A(0,0)*A(2,1)) / det;
  Ainv(2,2) = (A(0,0)*A(1,1)-A(0,1)*A(1,0)) / det;

  return det;
}

// =================================================================================================

namespace Gauss {

// -------------------------------------------------------------------------------------------------

inline size_t nip()
{
  return 8;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> xi()
{
  size_t nip  = 8;
  size_t ndim = 3;

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

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> w()
{
  size_t nip = 8;

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

namespace Nodal {

// -------------------------------------------------------------------------------------------------

inline size_t nip()
{
  return 8;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> xi()
{
  size_t nip  = 8;
  size_t ndim = 3;

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

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> w()
{
  size_t nip = 8;

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

inline Quadrature::Quadrature(const xt::xtensor<double,3> &x, const xt::xtensor<double,2> &xi,
  const xt::xtensor<double,1> &w) : m_x(x), m_w(w), m_xi(xi)
{
  // check input
  assert( m_x.shape()[1] == m_nne  );
  assert( m_x.shape()[2] == m_ndim );

  // extract number of elements and number of integration points
  m_nelem = m_x.shape()[0];
  m_nip   = m_w.size();

  // check input
  assert( m_xi.shape()[0] == m_nip  );
  assert( m_xi.shape()[1] == m_ndim );
  assert( m_w .size()     == m_nip  );

  // allocate arrays
  m_N    = xt::empty<double>({         m_nip, m_nne        });
  m_dNxi = xt::empty<double>({         m_nip, m_nne, m_ndim});
  m_dNx  = xt::empty<double>({m_nelem, m_nip, m_nne, m_ndim});
  m_vol  = xt::empty<double>({m_nelem, m_nip               });

  // shape functions
  for ( size_t q = 0 ; q < m_nip ; ++q )
  {
    m_N(q,0) = .125 * (1.-m_xi(q,0)) * (1.-m_xi(q,1)) * (1.-m_xi(q,2));
    m_N(q,1) = .125 * (1.+m_xi(q,0)) * (1.-m_xi(q,1)) * (1.-m_xi(q,2));
    m_N(q,2) = .125 * (1.+m_xi(q,0)) * (1.+m_xi(q,1)) * (1.-m_xi(q,2));
    m_N(q,3) = .125 * (1.-m_xi(q,0)) * (1.+m_xi(q,1)) * (1.-m_xi(q,2));
    m_N(q,4) = .125 * (1.-m_xi(q,0)) * (1.-m_xi(q,1)) * (1.+m_xi(q,2));
    m_N(q,5) = .125 * (1.+m_xi(q,0)) * (1.-m_xi(q,1)) * (1.+m_xi(q,2));
    m_N(q,6) = .125 * (1.+m_xi(q,0)) * (1.+m_xi(q,1)) * (1.+m_xi(q,2));
    m_N(q,7) = .125 * (1.-m_xi(q,0)) * (1.+m_xi(q,1)) * (1.+m_xi(q,2));
  }

  // shape function gradients in local coordinates
  for ( size_t q = 0 ; q < m_nip ; ++q )
  {
    // - dN / dxi_0
    m_dNxi(q,0,0) = -.125*(1.-m_xi(q,1))*(1.-m_xi(q,2));
    m_dNxi(q,1,0) = +.125*(1.-m_xi(q,1))*(1.-m_xi(q,2));
    m_dNxi(q,2,0) = +.125*(1.+m_xi(q,1))*(1.-m_xi(q,2));
    m_dNxi(q,3,0) = -.125*(1.+m_xi(q,1))*(1.-m_xi(q,2));
    m_dNxi(q,4,0) = -.125*(1.-m_xi(q,1))*(1.+m_xi(q,2));
    m_dNxi(q,5,0) = +.125*(1.-m_xi(q,1))*(1.+m_xi(q,2));
    m_dNxi(q,6,0) = +.125*(1.+m_xi(q,1))*(1.+m_xi(q,2));
    m_dNxi(q,7,0) = -.125*(1.+m_xi(q,1))*(1.+m_xi(q,2));
    // - dN / dxi_1
    m_dNxi(q,0,1) = -.125*(1.-m_xi(q,0))*(1.-m_xi(q,2));
    m_dNxi(q,1,1) = -.125*(1.+m_xi(q,0))*(1.-m_xi(q,2));
    m_dNxi(q,2,1) = +.125*(1.+m_xi(q,0))*(1.-m_xi(q,2));
    m_dNxi(q,3,1) = +.125*(1.-m_xi(q,0))*(1.-m_xi(q,2));
    m_dNxi(q,4,1) = -.125*(1.-m_xi(q,0))*(1.+m_xi(q,2));
    m_dNxi(q,5,1) = -.125*(1.+m_xi(q,0))*(1.+m_xi(q,2));
    m_dNxi(q,6,1) = +.125*(1.+m_xi(q,0))*(1.+m_xi(q,2));
    m_dNxi(q,7,1) = +.125*(1.-m_xi(q,0))*(1.+m_xi(q,2));
    // - dN / dxi_2
    m_dNxi(q,0,2) = -.125*(1.-m_xi(q,0))*(1.-m_xi(q,1));
    m_dNxi(q,1,2) = -.125*(1.+m_xi(q,0))*(1.-m_xi(q,1));
    m_dNxi(q,2,2) = -.125*(1.+m_xi(q,0))*(1.+m_xi(q,1));
    m_dNxi(q,3,2) = -.125*(1.-m_xi(q,0))*(1.+m_xi(q,1));
    m_dNxi(q,4,2) = +.125*(1.-m_xi(q,0))*(1.-m_xi(q,1));
    m_dNxi(q,5,2) = +.125*(1.+m_xi(q,0))*(1.-m_xi(q,1));
    m_dNxi(q,6,2) = +.125*(1.+m_xi(q,0))*(1.+m_xi(q,1));
    m_dNxi(q,7,2) = +.125*(1.-m_xi(q,0))*(1.+m_xi(q,1));
  }

  // compute the shape function gradients, based on "x"
  compute_dN();
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::dV(xt::xtensor<double,2> &qscalar) const
{
  assert( qscalar.shape()[0] == m_nelem );
  assert( qscalar.shape()[1] == m_nip   );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t q = 0 ; q < m_nip ; ++q )
      qscalar(e,q) = m_vol(e,q);
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::dV(xt::xtensor<double,4> &qtensor) const
{
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nne   );
  assert( qtensor.shape()[2] == m_ndim  );
  assert( qtensor.shape()[3] == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t q = 0 ; q < m_nip ; ++q )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        for ( size_t j = 0 ; j < m_ndim ; ++j )
          qtensor(e,q,i,j) = m_vol(e,q);
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::update_x(const xt::xtensor<double,3> &x)
{
  assert( x.shape()[0] == m_nelem    );
  assert( x.shape()[1] == m_nne      );
  assert( x.shape()[2] == m_ndim     );
  assert( x.size()     == m_x.size() );

  // update positions
  xt::noalias(m_x) = x;

  // update the shape function gradients for the new "x"
  compute_dN();
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::compute_dN()
{
  // loop over all elements (in parallel)
  #pragma omp parallel
  {
    // allocate local variables
    T2 J, Jinv;

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias nodal positions
      auto x = xt::adapt(&m_x(e,0,0), xt::xshape<m_nne,m_ndim>());

      // loop over integration points
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        // - alias
        auto dNxi = xt::adapt(&m_dNxi(  q,0,0), xt::xshape<m_nne,m_ndim>());
        auto dNx  = xt::adapt(&m_dNx (e,q,0,0), xt::xshape<m_nne,m_ndim>());

        // - zero-initialize
        J.fill(0.0);

        // - Jacobian
        for ( size_t m = 0 ; m < m_nne ; ++m )
          for ( size_t i = 0 ; i < m_ndim ; ++i )
            for ( size_t j = 0 ; j < m_ndim ; ++j )
              J(i,j) += dNxi(m,i) * x(m,j);

        // - determinant and inverse of the Jacobian
        double Jdet = inv(J, Jinv);

        // - shape function gradients wrt global coordinates (loops partly unrolled for efficiency)
        //   dNx(m,i) += Jinv(i,j) * dNxi(m,j);
        for ( size_t m = 0 ; m < m_nne ; ++m )
        {
          dNx(m,0) = Jinv(0,0) * dNxi(m,0) + Jinv(0,1) * dNxi(m,1) + Jinv(0,2) * dNxi(m,2);
          dNx(m,1) = Jinv(1,0) * dNxi(m,0) + Jinv(1,1) * dNxi(m,1) + Jinv(1,2) * dNxi(m,2);
          dNx(m,2) = Jinv(2,0) * dNxi(m,0) + Jinv(2,1) * dNxi(m,1) + Jinv(2,2) * dNxi(m,2);
        }

        // - integration point volume
        m_vol(e,q) = m_w(q) * Jdet;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::gradN_vector(
  const xt::xtensor<double,3> &elemvec, xt::xtensor<double,4> &qtensor) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nne   );
  assert( qtensor.shape()[2] == m_ndim  );
  assert( qtensor.shape()[3] == m_ndim  );

  // zero-initialize output: matrix of tensors
  qtensor.fill(0.0);

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias element vector (e.g. nodal displacements)
    auto u = xt::adapt(&elemvec(e,0,0), xt::xshape<m_nne,m_ndim>());

    // loop over all integration points in element "e"
    for ( size_t q = 0 ; q < m_nip ; ++q )
    {
      // - alias
      auto dNx   = xt::adapt(&m_dNx  (e,q,0,0), xt::xshape<m_nne ,m_ndim>());
      auto gradu = xt::adapt(&qtensor(e,q,0,0), xt::xshape<m_ndim,m_ndim>());

      // - evaluate dyadic product
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          for ( size_t j = 0 ; j < m_ndim ; ++j )
            gradu(i,j) += dNx(m,i) * u(m,j);
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::gradN_vector_T(
  const xt::xtensor<double,3> &elemvec, xt::xtensor<double,4> &qtensor) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nne   );
  assert( qtensor.shape()[2] == m_ndim  );
  assert( qtensor.shape()[3] == m_ndim  );

  // zero-initialize output: matrix of tensors
  qtensor.fill(0.0);

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias element vector (e.g. nodal displacements)
    auto u = xt::adapt(&elemvec(e,0,0), xt::xshape<m_nne,m_ndim>());

    // loop over all integration points in element "e"
    for ( size_t q = 0 ; q < m_nip ; ++q )
    {
      // - alias
      auto dNx   = xt::adapt(&m_dNx  (e,q,0,0), xt::xshape<m_nne ,m_ndim>());
      auto gradu = xt::adapt(&qtensor(e,q,0,0), xt::xshape<m_ndim,m_ndim>());

      // - evaluate transpose of dyadic product
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t i = 0 ; i < m_ndim ; ++i )
          for ( size_t j = 0 ; j < m_ndim ; ++j )
            gradu(j,i) += dNx(m,i) * u(m,j);
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::symGradN_vector(
  const xt::xtensor<double,3> &elemvec, xt::xtensor<double,4> &qtensor) const
{
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nne   );
  assert( qtensor.shape()[2] == m_ndim  );
  assert( qtensor.shape()[3] == m_ndim  );

  // zero-initialize output: matrix of tensors
  qtensor.fill(0.0);

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias element vector (e.g. nodal displacements)
    auto u = xt::adapt(&elemvec(e,0,0), xt::xshape<m_nne,m_ndim>());

    // loop over all integration points in element "e"
    for ( size_t q = 0 ; q < m_nip ; ++q )
    {
      // - alias
      auto dNx = xt::adapt(&m_dNx  (e,q,0,0), xt::xshape<m_nne ,m_ndim>());
      auto eps = xt::adapt(&qtensor(e,q,0,0), xt::xshape<m_ndim,m_ndim>());

      // - evaluate symmetrized dyadic product
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

inline void Quadrature::int_N_scalar_NT_dV(
  const xt::xtensor<double,2> &qscalar, xt::xtensor<double,3> &elemmat) const
{
  assert( qscalar.shape()[0] == m_nelem      );
  assert( qscalar.shape()[1] == m_nip        );
  assert( elemmat.shape()[0] == m_nelem      );
  assert( elemmat.shape()[1] == m_nne*m_ndim );
  assert( elemmat.shape()[2] == m_nne*m_ndim );

  // zero-initialize: matrix of matrices
  elemmat.fill(0.0);

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias (e.g. mass matrix)
    auto M = xt::adapt(&elemmat(e,0,0), xt::xshape<m_nne*m_ndim,m_nne*m_ndim>());

    // loop over all integration points in element "e"
    for ( size_t q = 0 ; q < m_nip ; ++q )
    {
      // - alias
      auto  N   = xt::adapt(&m_N(q,0), xt::xshape<m_nne>());
      auto& vol = m_vol  (e,q);
      auto& rho = qscalar(e,q);

      // - evaluate scalar product, for all dimensions, and assemble
      //   M(m*ndim+i,n*ndim+i) += N(m) * scalar * N(n) * dV
      for ( size_t m = 0 ; m < m_nne ; ++m ) {
        for ( size_t n = 0 ; n < m_nne ; ++n ) {
          M(m*m_ndim+0, n*m_ndim+0) += N(m) * rho * N(n) * vol;
          M(m*m_ndim+1, n*m_ndim+1) += N(m) * rho * N(n) * vol;
          M(m*m_ndim+2, n*m_ndim+2) += N(m) * rho * N(n) * vol;
        }
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::int_gradN_dot_tensor2_dV(const xt::xtensor<double,4> &qtensor,
  xt::xtensor<double,3> &elemvec) const
{
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nip   );
  assert( qtensor.shape()[2] == m_ndim  );
  assert( qtensor.shape()[3] == m_ndim  );
  assert( elemvec.shape()[0] == m_nelem );
  assert( elemvec.shape()[1] == m_nne   );
  assert( elemvec.shape()[2] == m_ndim  );

  // zero-initialize output: matrix of vectors
  elemvec.fill(0.0);

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias (e.g. nodal force)
    auto f = xt::adapt(&elemvec(e,0,0), xt::xshape<m_nne,m_ndim>());

    // loop over all integration points in element "e"
    for ( size_t q = 0 ; q < m_nip ; ++q )
    {
      // - alias
      auto  dNx = xt::adapt(&m_dNx  (e,q,0,0), xt::xshape<m_nne ,m_ndim>());
      auto  sig = xt::adapt(&qtensor(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
      auto& vol = m_vol(e,q);

      // - evaluate dot product, and assemble
      for ( size_t m = 0 ; m < m_nne ; ++m )
      {
        f(m,0) += ( dNx(m,0) * sig(0,0) + dNx(m,1) * sig(1,0) + dNx(m,2) * sig(2,0) ) * vol;
        f(m,1) += ( dNx(m,0) * sig(0,1) + dNx(m,1) * sig(1,1) + dNx(m,2) * sig(2,1) ) * vol;
        f(m,2) += ( dNx(m,0) * sig(0,2) + dNx(m,1) * sig(1,2) + dNx(m,2) * sig(2,2) ) * vol;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::int_gradN_dot_tensor4_dot_gradNT_dV(const xt::xtensor<double,6> &qtensor,
  xt::xtensor<double,3> &elemmat) const
{
  assert( qtensor.shape()[0] == m_nelem );
  assert( qtensor.shape()[1] == m_nip   );
  assert( qtensor.shape()[2] == m_ndim  );
  assert( qtensor.shape()[3] == m_ndim  );
  assert( qtensor.shape()[4] == m_ndim  );
  assert( qtensor.shape()[5] == m_ndim  );

  assert( elemmat.shape()[0] == m_nelem      );
  assert( elemmat.shape()[1] == m_nne*m_ndim );
  assert( elemmat.shape()[2] == m_nne*m_ndim );

  // zero-initialize output: matrix of vector
  elemmat.fill(0.0);

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias (e.g. nodal force)
    auto K = xt::adapt(&elemmat(e,0,0), xt::xshape<m_nne*m_ndim,m_nne*m_ndim>());

    // loop over all integration points in element "e"
    for ( size_t q = 0 ; q < m_nip ; ++q ){

      // - alias
      auto  dNx = xt::adapt(&m_dNx(e,q,0,0), xt::xshape<m_nne,m_ndim>());
      auto  C   = xt::adapt(&qtensor(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
      auto& vol = m_vol(e,q);

      // - evaluate dot product, and assemble
      for ( size_t m = 0 ; m < m_nne ; ++m )
        for ( size_t n = 0 ; n < m_nne ; ++n )
          for ( size_t i = 0 ; i < m_ndim ; ++i )
            for ( size_t j = 0 ; j < m_ndim ; ++j )
              for ( size_t k = 0 ; k < m_ndim ; ++k )
                for ( size_t l = 0 ; l < m_ndim ; ++l )
                  K(m*m_ndim+j, n*m_ndim+k) += dNx(m,i) * C(i,j,k,l) * dNx(n,l) * vol;
     }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Quadrature::dV() const
{
  xt::xtensor<double,2> out = xt::empty<double>({m_nelem, m_nip});

  this->dV(out);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Quadrature::dVtensor() const
{
  xt::xtensor<double,4> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  this->dV(out);

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Quadrature::gradN_vector(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,4> qtensor = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  this->gradN_vector(elemvec, qtensor);

  return qtensor;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Quadrature::gradN_vector_T(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,4> qtensor = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  this->gradN_vector_T(elemvec, qtensor);

  return qtensor;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Quadrature::symGradN_vector(const xt::xtensor<double,3> &elemvec) const
{
  xt::xtensor<double,4> qtensor = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  this->symGradN_vector(elemvec, qtensor);

  return qtensor;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Quadrature::int_N_scalar_NT_dV(
  const xt::xtensor<double,2> &qscalar) const
{
  xt::xtensor<double,3> elemmat = xt::empty<double>({m_nelem, m_nne*m_ndim, m_nne*m_ndim});

  this->int_N_scalar_NT_dV(qscalar, elemmat);

  return elemmat;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Quadrature::int_gradN_dot_tensor2_dV(
  const xt::xtensor<double,4> &qtensor) const
{
  xt::xtensor<double,3> elemvec = xt::empty<double>({m_nelem, m_nne, m_ndim});

  this->int_gradN_dot_tensor2_dV(qtensor, elemvec);

  return elemvec;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Quadrature::int_gradN_dot_tensor4_dot_gradNT_dV(
  const xt::xtensor<double,6> &qtensor) const
 {
   xt::xtensor<double,3> elemmat = xt::empty<double>({m_nelem, m_ndim*m_nne, m_ndim*m_nne});

   this->int_gradN_dot_tensor4_dot_gradNT_dV(qtensor, elemmat);

   return elemmat;
 }

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

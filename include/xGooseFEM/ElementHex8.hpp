/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef XGOOSEFEM_ELEMENTHEX8_CPP
#define XGOOSEFEM_ELEMENTHEX8_CPP

// -------------------------------------------------------------------------------------------------

#include "ElementHex8.h"

// ==================================== GooseFEM::Element::Hex8 ====================================

namespace xGooseFEM {
namespace Element {
namespace Hex8 {

// ======================================== tensor algebra =========================================

inline double inv(const T2 &A, T2 &Ainv)
{
  // compute determinant
  double det = ( A[0] * A[4] * A[8] +
                 A[1] * A[5] * A[6] +
                 A[2] * A[3] * A[7] ) -
               ( A[2] * A[4] * A[6] +
                 A[1] * A[3] * A[8] +
                 A[0] * A[5] * A[7] );

  // compute inverse
  Ainv[0] = (A[4]*A[8]-A[5]*A[7]) / det;
  Ainv[1] = (A[2]*A[7]-A[1]*A[8]) / det;
  Ainv[2] = (A[1]*A[5]-A[2]*A[4]) / det;
  Ainv[3] = (A[5]*A[6]-A[3]*A[8]) / det;
  Ainv[4] = (A[0]*A[8]-A[2]*A[6]) / det;
  Ainv[5] = (A[2]*A[3]-A[0]*A[5]) / det;
  Ainv[6] = (A[3]*A[7]-A[4]*A[6]) / det;
  Ainv[7] = (A[1]*A[6]-A[0]*A[7]) / det;
  Ainv[8] = (A[0]*A[4]-A[1]*A[3]) / det;

  return det;
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

// ----------------------------------- integration point weights -----------------------------------

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

// ----------------------------------- integration point weights -----------------------------------

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

// ------------------------------------------ constructor ------------------------------------------

inline Quadrature::Quadrature(const xt::xtensor<double,3> &x) : m_x(x)
{
  assert( m_x.shape()[1] == m_nne  );
  assert( m_x.shape()[2] == m_ndim );

  // integration scheme
  m_nip = Gauss::nip();
  m_xi  = Gauss::xi();
  m_w   = Gauss::w();

  // extract number of elements
  m_nelem = m_x.shape()[0];

  // allocate arrays
  m_N    = xt::empty<double>({         m_nip, m_nne        });
  m_dNxi = xt::empty<double>({         m_nip, m_nne, m_ndim});
  m_dNx  = xt::empty<double>({m_nelem, m_nip, m_nne, m_ndim});
  m_vol  = xt::empty<double>({m_nelem, m_nip               });

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

  // extract number of elements and number of integration points
  m_nelem = m_x.shape()[0];
  m_nip   = m_w.size();

  assert( m_xi.shape()[0] == m_nip  );
  assert( m_xi.shape()[1] == m_ndim );
  assert( m_w .size()     == m_nip  );

  // allocate arrays
  m_N    = xt::empty<double>({         m_nip, m_nne        });
  m_dNxi = xt::empty<double>({         m_nip, m_nne, m_ndim});
  m_dNx  = xt::empty<double>({m_nelem, m_nip, m_nne, m_ndim});
  m_vol  = xt::empty<double>({m_nelem, m_nip               });

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
  assert( qtensor.shape()[2] == m_ndim  );
  assert( qtensor.shape()[3] == m_ndim  );

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t k = 0 ; k < m_nip ; ++k )
      for ( size_t i = 0 ; i < m_ndim ; ++i )
        for ( size_t j = 0 ; j < m_ndim ; ++j )
          qtensor(e,k,i,j) = m_vol(e,k);
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
  xt::noalias(m_x) = x;

  // update the shape function gradients for the new "x"
  compute_dN();
}

// ------------------------ shape function gradients in global coordinates -------------------------

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
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        auto dNxi = xt::adapt(&m_dNxi(  k,0,0), xt::xshape<m_nne,m_ndim>());
        auto dNx  = xt::adapt(&m_dNx (e,k,0,0), xt::xshape<m_nne,m_ndim>());

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
        m_vol(e,k) = m_w(k) * Jdet;
      }
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
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto dNx   = xt::adapt(&m_dNx  (e,k,0,0), xt::xshape<m_nne ,m_ndim>());
      auto gradu = xt::adapt(&qtensor(e,k,0,0), xt::xshape<m_ndim,m_ndim>());

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
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto dNx   = xt::adapt(&m_dNx  (e,k,0,0), xt::xshape<m_nne ,m_ndim>());
      auto gradu = xt::adapt(&qtensor(e,k,0,0), xt::xshape<m_ndim,m_ndim>());

      // - evaluate transpose of dyadic product
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
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto dNx = xt::adapt(&m_dNx  (e,k,0,0), xt::xshape<m_nne ,m_ndim>());
      auto eps = xt::adapt(&qtensor(e,k,0,0), xt::xshape<m_ndim,m_ndim>());

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
  elemmat.fill(0.0);

  // loop over all elements (in parallel)
  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // alias (e.g. mass matrix)
    auto M = xt::adapt(&elemmat(e,0,0), xt::xshape<m_nne*m_ndim,m_nne*m_ndim>());

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto  N   = xt::adapt(&m_N(k,0), xt::xshape<m_nne>());
      auto& vol = m_vol  (e,k);
      auto& rho = qscalar(e,k);

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
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      auto  dNx = xt::adapt(&m_dNx  (e,k,0,0), xt::xshape<m_nne ,m_ndim>());
      auto  sig = xt::adapt(&qtensor(e,k,0,0), xt::xshape<m_ndim,m_ndim>());
      auto& vol = m_vol(e,k);

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

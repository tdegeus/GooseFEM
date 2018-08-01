/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_ELEMENTQUAD4_CPP
#define GOOSEFEM_ELEMENTQUAD4_CPP

// -------------------------------------------------------------------------------------------------

#include "ElementQuad4.h"

// =================================== GooseFEM::Element::Quad4 ====================================

namespace GooseFEM {
namespace Element {
namespace Quad4 {

// ======================================== tensor algebra =========================================

inline double inv(const T2 &A, T2 &Ainv)
{
  // compute determinant
  double det = A[0] * A[3] - A[1] * A[2];

  // compute inverse
  Ainv[0] =       A[3] / det;
  Ainv[1] = -1. * A[1] / det;
  Ainv[2] = -1. * A[2] / det;
  Ainv[3] =       A[0] / det;

  return det;
}

// ================================ GooseFEM::Element::Quad4::Gauss ================================

namespace Gauss {

// --------------------------------- number of integration points ----------------------------------

inline size_t nip()
{
  return 4;
}

// ----------------------- integration point coordinates (local coordinates) -----------------------

inline ArrD xi()
{
  size_t nip  = 4;
  size_t ndim = 2;

  ArrD xi({nip,ndim});

  xi(0,0) = -1./std::sqrt(3.);    xi(0,1) = -1./std::sqrt(3.);
  xi(1,0) = +1./std::sqrt(3.);    xi(1,1) = -1./std::sqrt(3.);
  xi(2,0) = +1./std::sqrt(3.);    xi(2,1) = +1./std::sqrt(3.);
  xi(3,0) = -1./std::sqrt(3.);    xi(3,1) = +1./std::sqrt(3.);

  return xi;
}

// ----------------------------------- integration point weights -----------------------------------

inline ArrD w()
{
  size_t nip = 4;

  ArrD w({nip});

  w(0) = 1.;
  w(1) = 1.;
  w(2) = 1.;
  w(3) = 1.;

  return w;
}

// -------------------------------------------------------------------------------------------------

}

// ================================ GooseFEM::Element::Quad4::Nodal ================================

namespace Nodal {

// --------------------------------- number of integration points ----------------------------------

inline size_t nip()
{
  return 4;
}

// ----------------------- integration point coordinates (local coordinates) -----------------------

inline ArrD xi()
{
  size_t nip  = 4;
  size_t ndim = 2;

  ArrD xi({nip,ndim});

  xi(0,0) = -1.;    xi(0,1) = -1.;
  xi(1,0) = +1.;    xi(1,1) = -1.;
  xi(2,0) = +1.;    xi(2,1) = +1.;
  xi(3,0) = -1.;    xi(3,1) = +1.;

  return xi;
}

// ----------------------------------- integration point weights -----------------------------------

inline ArrD w()
{
  size_t nip = 4;

  ArrD w({nip});

  w(0) = 1.;
  w(1) = 1.;
  w(2) = 1.;
  w(3) = 1.;

  return w;
}

// -------------------------------------------------------------------------------------------------

}

// =================================================================================================

// ------------------------------------------ constructor ------------------------------------------

inline Quadrature::Quadrature(const ArrD &x)
: m_x(x)
{
  assert( m_x.rank()   == 3      );
  assert( m_x.shape(1) == m_nne  );
  assert( m_x.shape(2) == m_ndim );

  // extract shape
  m_nelem = m_x.shape(0);

  // integration scheme
  m_nip = Gauss::nip();
  m_xi  = Gauss::xi();
  m_w   = Gauss::w();

  // allocate arrays
  // - shape functions
  m_N.resize({m_nip,m_nne});
  // - shape function gradients in local coordinates
  m_dNxi.resize({m_nip,m_nne,m_ndim});
  // - shape function gradients in global coordinates
  m_dNx.resize({m_nelem,m_nip,m_nne,m_ndim});
  // - integration point volume
  m_vol.resize({m_nelem,m_nip});

  // shape functions
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    m_N(k,0) = .25 * (1.-m_xi(k,0)) * (1.-m_xi(k,1));
    m_N(k,1) = .25 * (1.+m_xi(k,0)) * (1.-m_xi(k,1));
    m_N(k,2) = .25 * (1.+m_xi(k,0)) * (1.+m_xi(k,1));
    m_N(k,3) = .25 * (1.-m_xi(k,0)) * (1.+m_xi(k,1));
  }

  // shape function gradients in local coordinates
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    // - dN / dxi_0
    m_dNxi(k,0,0) = -.25*(1.-m_xi(k,1));
    m_dNxi(k,1,0) = +.25*(1.-m_xi(k,1));
    m_dNxi(k,2,0) = +.25*(1.+m_xi(k,1));
    m_dNxi(k,3,0) = -.25*(1.+m_xi(k,1));
    // - dN / dxi_1
    m_dNxi(k,0,1) = -.25*(1.-m_xi(k,0));
    m_dNxi(k,1,1) = -.25*(1.+m_xi(k,0));
    m_dNxi(k,2,1) = +.25*(1.+m_xi(k,0));
    m_dNxi(k,3,1) = +.25*(1.-m_xi(k,0));
  }

  // compute the shape function gradients, based on "x"
  compute_dN();
}

// ------------------------------------------ constructor ------------------------------------------

inline Quadrature::Quadrature(const ArrD &x, const ArrD &xi, const ArrD &w)
: m_x(x), m_w(w), m_xi(xi)
{
  assert( m_x.rank()   == 3      );
  assert( m_x.shape(1) == m_nne  );
  assert( m_x.shape(2) == m_ndim );

  // extract shape
  m_nelem = m_x.shape(0);
  m_nip   = m_w.size();

  assert( m_xi.rank()   == 2      );
  assert( m_xi.shape(0) == m_nip  );
  assert( m_xi.shape(1) == m_ndim );
  assert( m_w .rank()   == 1      );
  assert( m_w .size()   == m_nip  );

  // allocate arrays
  // - shape functions
  m_N.resize({m_nip,m_nne});
  // - shape function gradients in local coordinates
  m_dNxi.resize({m_nip,m_nne,m_ndim});
  // - shape function gradients in global coordinates
  m_dNx.resize({m_nelem,m_nip,m_nne,m_ndim});
  // - integration point volume
  m_vol.resize({m_nelem,m_nip});

  // shape functions
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    m_N(k,0) = .25 * (1.-m_xi(k,0)) * (1.-m_xi(k,1));
    m_N(k,1) = .25 * (1.+m_xi(k,0)) * (1.-m_xi(k,1));
    m_N(k,2) = .25 * (1.+m_xi(k,0)) * (1.+m_xi(k,1));
    m_N(k,3) = .25 * (1.-m_xi(k,0)) * (1.+m_xi(k,1));
  }

  // shape function gradients in local coordinates
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    // - dN / dxi_0
    m_dNxi(k,0,0) = -.25*(1.-m_xi(k,1));
    m_dNxi(k,1,0) = +.25*(1.-m_xi(k,1));
    m_dNxi(k,2,0) = +.25*(1.+m_xi(k,1));
    m_dNxi(k,3,0) = -.25*(1.+m_xi(k,1));
    // - dN / dxi_1
    m_dNxi(k,0,1) = -.25*(1.-m_xi(k,0));
    m_dNxi(k,1,1) = -.25*(1.+m_xi(k,0));
    m_dNxi(k,2,1) = +.25*(1.+m_xi(k,0));
    m_dNxi(k,3,1) = +.25*(1.-m_xi(k,0));
  }

  // compute the shape function gradients, based on "x"
  compute_dN();
}

// --------------------------- integration volume (per tensor-component) ---------------------------

inline ArrD Quadrature::dV(size_t ncomp) const
{
  if ( ncomp == 0 ) return m_vol;

  ArrD out = ArrD::Zero({m_nelem, m_nip, ncomp});

  #pragma omp parallel for
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t k = 0 ; k < m_nip ; ++k )
      for ( size_t i = 0 ; i < ncomp ; ++i )
        out(e,k,i) = m_vol(e,k);

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

inline void Quadrature::update_x(const ArrD &x)
{
  assert( x.rank()   == 3          );
  assert( x.shape(0) == m_nelem    );
  assert( x.shape(1) == m_nne      );
  assert( x.shape(2) == m_ndim     );
  assert( x.size()   == m_x.size() );

  // update positions
  m_x.setCopy(x.begin(), x.end());

  // update the shape function gradients for the new "x"
  compute_dN();
}

// ------------------------ shape function gradients in global coordinates -------------------------

inline void Quadrature::compute_dN()
{
  #pragma omp parallel
  {
    // - allocate
    T2 J, Jinv;
    cppmat::tiny::matrix<double,m_nne,m_ndim> dNx;
    cppmat::view::matrix<double,m_nne,m_ndim> dNxi, x;

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias nodal positions
      x.setMap(&m_x(e));

      // loop over integration points
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        dNxi.setMap(&m_dNxi(k));

        // - Jacobian (loops unrolled for efficiency)
        //   J(i,j) += dNxi(m,i) * x(m,j);
        J(0,0) = dNxi(0,0)*x(0,0) + dNxi(1,0)*x(1,0) + dNxi(2,0)*x(2,0) + dNxi(3,0)*x(3,0);
        J(0,1) = dNxi(0,0)*x(0,1) + dNxi(1,0)*x(1,1) + dNxi(2,0)*x(2,1) + dNxi(3,0)*x(3,1);
        J(1,0) = dNxi(0,1)*x(0,0) + dNxi(1,1)*x(1,0) + dNxi(2,1)*x(2,0) + dNxi(3,1)*x(3,0);
        J(1,1) = dNxi(0,1)*x(0,1) + dNxi(1,1)*x(1,1) + dNxi(2,1)*x(2,1) + dNxi(3,1)*x(3,1);

        // - determinant and inverse of the Jacobian
        double Jdet = inv(J, Jinv);

        // - shape function gradients wrt global coordinates (loops partly unrolled for efficiency)
        //   dNx(m,i) += Jinv(i,j) * dNxi(m,j);
        for ( size_t m = 0 ; m < m_nne ; ++m )
        {
          dNx(m,0) = Jinv(0,0) * dNxi(m,0) + Jinv(0,1) * dNxi(m,1);
          dNx(m,1) = Jinv(1,0) * dNxi(m,0) + Jinv(1,1) * dNxi(m,1);
        }

        // - copy to matrix: shape function gradients (global coordinates)
        dNx.copyTo(m_dNx.item(e,k));

        // - copy to matrix: integration point volume
        m_vol(e,k) = m_w(k) * Jdet;
      }
    }
  }
}

// ------------------- dyadic product "qtensor(i,j) = dNdx(m,i) * elemvec(m,j)" --------------------

template<class T>
inline void Quadrature::gradN_vector(const ArrD &elemvec, ArrD &qtensor) const
{
  assert( elemvec.rank()   == 3         );
  assert( elemvec.shape(0) == m_nelem   );
  assert( elemvec.shape(1) == m_nne     );
  assert( elemvec.shape(2) == m_ndim    );
  assert( qtensor.rank()   == 3         );
  assert( qtensor.shape(0) == m_nelem   );
  assert( qtensor.shape(1) == m_nip     );
  assert( qtensor.shape(2) == T::Size() );

  // zero-initialize output: matrix of tensors
  qtensor *= 0.0;

  #pragma omp parallel
  {
    // intermediate quantities and local views
    T gradu;
    cppmat::view::matrix<double,m_nne,m_ndim> dNx, u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias element vector (e.g. nodal displacements)
      u.setMap(&elemvec(e));

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        dNx.setMap(&m_dNx(e,k));

        // - evaluate dyadic product (loops unrolled for efficiency)
        //   gradu(i,j) += dNx(m,i) * u(m,j)
        gradu(0,0) = dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        gradu(0,1) = dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1);
        gradu(1,0) = dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0);
        gradu(1,1) = dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);

        // - copy resulting integration point tensor
        std::copy(gradu.begin(), gradu.end(), qtensor.item(e,k));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class T>
inline ArrD Quadrature::gradN_vector(const ArrD &elemvec) const
{
  ArrD qtensor({m_nelem, m_nip, T::Size()});

  this->gradN_vector(elemvec, qtensor);

  return qtensor;
}

// ---------------------------------- transpose of "GradN_vector" ----------------------------------

template<class T>
inline void Quadrature::gradN_vector_T(const ArrD &elemvec, ArrD &qtensor) const
{
  assert( elemvec.rank()   == 3       );
  assert( elemvec.shape(0) == m_nelem );
  assert( elemvec.shape(1) == m_nne   );
  assert( elemvec.shape(2) == m_ndim  );

  // zero-initialize output: matrix of tensors
  qtensor *= 0.0;

  #pragma omp parallel
  {
    // intermediate quantities and local views
    T gradu;
    cppmat::view::matrix<double,m_nne,m_ndim> dNx, u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias element vector (e.g. nodal displacements)
      u.setMap(&elemvec(e));

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias shape function gradients (global coordinates)
        dNx.setMap(&m_dNx(e,k));

        // - evaluate dyadic product (loops unrolled for efficiency)
        //   gradu(j,i) += dNx(m,i) * u(m,j)
        gradu(0,0) = dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        gradu(1,0) = dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1);
        gradu(0,1) = dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0);
        gradu(1,1) = dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);

        // - copy resulting integration point tensor
        std::copy(gradu.begin(), gradu.end(), qtensor.item(e,k));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class T>
inline ArrD Quadrature::gradN_vector_T(const ArrD &elemvec) const
{
  ArrD qtensor({m_nelem, m_nip, T::Size()});

  this->gradN_vector_T(elemvec, qtensor);

  return qtensor;
}

// ------------------------------- symmetric part of "gradN_vector" --------------------------------

template<class T>
inline void Quadrature::symGradN_vector(const ArrD &elemvec, ArrD &qtensor) const
{
  assert( elemvec.rank()   == 3         );
  assert( elemvec.shape(0) == m_nelem   );
  assert( elemvec.shape(1) == m_nne     );
  assert( elemvec.shape(2) == m_ndim    );
  assert( qtensor.rank()   == 3         );
  assert( qtensor.shape(0) == m_nelem   );
  assert( qtensor.shape(1) == m_nip     );
  assert( qtensor.shape(2) == T::Size() );

  // zero-initialize output: matrix of tensors
  qtensor *= 0.0;

  #pragma omp parallel
  {
    // intermediate quantities and local views
    T eps;
    cppmat::view::matrix<double,m_nne,m_ndim> dNx, u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias element vector (e.g. nodal displacements)
      u.setMap(&elemvec(e));

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias shape function gradients (global coordinates)
        dNx.setMap(&m_dNx(e,k));

        // - evaluate symmetrized dyadic product (loops unrolled for efficiency)
        //   grad(i,j) += dNx(m,i) * u(m,j)
        //   eps (j,i)  = 0.5 * ( grad(i,j) + grad(j,i) )
        eps(0,0) =   dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        eps(1,1) =   dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);
        eps(0,1) = ( dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1) +
                     dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0) ) / 2.;
        eps(1,0) =   eps(0,1);

        // - copy resulting integration point tensor
        std::copy(eps.begin(), eps.end(), qtensor.item(e,k));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class T>
inline ArrD Quadrature::symGradN_vector(const ArrD &elemvec) const
{
  ArrD qtensor({m_nelem, m_nip, T::Size()});

  this->symGradN_vector(elemvec, qtensor);

  return qtensor;
}

// ------- scalar product "elemmat(m*ndim+i,n*ndim+i) = N(m) * qscalar * N(n)"; for all "i" --------

inline void Quadrature::int_N_scalar_NT_dV(const ArrD &qscalar, ArrD &elemmat) const
{
  assert( qscalar.rank()   == 2            );
  assert( qscalar.shape(0) == m_nelem      );
  assert( qscalar.shape(1) == m_nip        );
  assert( elemmat.rank()   == 3            );
  assert( elemmat.shape(0) == m_nelem      );
  assert( elemmat.shape(1) == m_nne*m_ndim );
  assert( elemmat.shape(2) == m_nne*m_ndim );

  // zero-initialize: matrix of matrices
  elemmat *= 0.0;

  #pragma omp parallel
  {
    // intermediate quantities and local views
    cppmat::tiny::matrix<double,m_nne*m_ndim,m_nne*m_ndim> M;
    cppmat::view::vector<double,m_nne> N;
    double rho, vol;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // zero-initialize (e.g. mass matrix)
      M.setZero();

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias shape functions
        N.setMap(&m_N(k));

        // - alias
        vol = m_vol  (e,k);  // integration point volume
        rho = qscalar(e,k);  // integration point scalar (e.g. density)

        // - evaluate scalar product, for all dimensions, and assemble
        //   M(m*ndim+i,n*ndim+i) += N(m) * scalar * N(n) * dV
        for ( size_t m = 0 ; m < m_nne ; ++m ) {
          for ( size_t n = 0 ; n < m_nne ; ++n ) {
            M(m*m_ndim+0, n*m_ndim+0) += N(m) * rho * N(n) * vol;
            M(m*m_ndim+1, n*m_ndim+1) += N(m) * rho * N(n) * vol;
          }
        }
      }

      // copy result to element matrix
      std::copy(M.begin(), M.end(), elemmat.item(e));
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::int_N_scalar_NT_dV(const ArrD &qscalar) const
{
  ArrD elemmat({m_nelem, m_nne*m_ndim, m_nne*m_ndim});

  this->int_N_scalar_NT_dV(qscalar, elemmat);

  return elemmat;
}

// ------------ integral of dot product "elemvec(m,j) += dNdx(m,i) * qtensor(i,j) * dV" ------------

template<class T>
inline void Quadrature::int_gradN_dot_tensor2_dV(const ArrD &qtensor, ArrD &elemvec) const
{
  assert( qtensor.rank()   == 3         );
  assert( qtensor.shape(0) == m_nelem   );
  assert( qtensor.shape(1) == m_nip     );
  assert( qtensor.shape(2) == T::Size() );
  assert( elemvec.rank()   == 3         );
  assert( elemvec.shape(0) == m_nelem   );
  assert( elemvec.shape(1) == m_nne     );
  assert( elemvec.shape(2) == m_ndim    );

  // zero-initialize output: matrix of vectors
  elemvec *= 0.0;

  #pragma omp parallel
  {
    // intermediate quantities and local views
    cppmat::view::matrix<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix<double,m_nne,m_ndim> f;
    double vol;
    T sig;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // zero-initialize (e.g. nodal force)
      f.setZero();

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        dNx.setMap (&m_dNx  (e,k)); // shape function gradients (global coordinates)
        sig.setCopy(&qtensor(e,k)); // integration point tensor (e.g. stress)
        vol = m_vol         (e,k);  // integration point volume

        // - evaluate dot product, and assemble (loops partly unrolled for efficiency)
        //   f(m,j) += dNx(m,i) * sig(i,j) * vol;
        for ( size_t m = 0 ; m < m_nne ; ++m )
        {
          f(m,0) += ( dNx(m,0) * sig(0,0) + dNx(m,1) * sig(1,0) ) * vol;
          f(m,1) += ( dNx(m,0) * sig(0,1) + dNx(m,1) * sig(1,1) ) * vol;
        }
      }

      // copy result to element vector
      std::copy(f.begin(), f.end(), elemvec.item(e));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class T>
inline ArrD Quadrature::int_gradN_dot_tensor2_dV(const ArrD &qtensor) const
{
  ArrD elemvec({m_nelem, m_nne, m_ndim});

  this->int_gradN_dot_tensor2_dV(qtensor, elemvec);

  return elemvec;
}

// ---------------------- wrappers with default storage (no template needed) -----------------------

inline void Quadrature::gradN_vector(const ArrD &elemvec, ArrD &qtensor) const
{
  return gradN_vector<cppmat::tiny::cartesian::tensor2<double,2>>(elemvec, qtensor);
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::gradN_vector(const ArrD &elemvec) const
{
  return gradN_vector<cppmat::tiny::cartesian::tensor2<double,2>>(elemvec);
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::gradN_vector_T(const ArrD &elemvec, ArrD &qtensor) const
{
  return gradN_vector_T<cppmat::tiny::cartesian::tensor2<double,2>>(elemvec, qtensor);
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::gradN_vector_T(const ArrD &elemvec) const
{
  return gradN_vector_T<cppmat::tiny::cartesian::tensor2<double,2>>(elemvec);
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::symGradN_vector(const ArrD &elemvec, ArrD &qtensor) const
{
  return symGradN_vector<cppmat::tiny::cartesian::tensor2s<double,2>>(elemvec, qtensor);
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::symGradN_vector(const ArrD &elemvec) const
{
  return symGradN_vector<cppmat::tiny::cartesian::tensor2s<double,2>>(elemvec);
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::int_gradN_dot_tensor2_dV(const ArrD &qtensor, ArrD &elemvec) const
{
  assert( qtensor.rank() == 3 ); // shape: [nelem, nip, #tensor-components]

  if ( qtensor.shape(2) == m_ndim*m_ndim )

    return int_gradN_dot_tensor2_dV<cppmat::tiny::cartesian::tensor2<double,2>>(qtensor, elemvec);

  else if ( qtensor.shape(2) == (m_ndim+1)*m_ndim/2 )

    return int_gradN_dot_tensor2_dV<cppmat::tiny::cartesian::tensor2s<double,2>>(qtensor, elemvec);

  else

    throw std::runtime_error("assert: qtensor.shape(2) == 4 or qtensor.shape(2) == 3");
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::int_gradN_dot_tensor2_dV(const ArrD &qtensor) const
{
  assert( qtensor.rank() == 3 ); // shape: [nelem, nip, #tensor-components]

  if ( qtensor.shape(2) == m_ndim*m_ndim )

    return int_gradN_dot_tensor2_dV<cppmat::tiny::cartesian::tensor2<double,2>>(qtensor);

  else if ( qtensor.shape(2) == (m_ndim+1)*m_ndim/2 )

    return int_gradN_dot_tensor2_dV<cppmat::tiny::cartesian::tensor2s<double,2>>(qtensor);

  else

    throw std::runtime_error("assert: qtensor.shape(2) == 4 or qtensor.shape(2) == 3");
}

// -------------------------------------------------------------------------------------------------

inline void Quadrature::int_gradN_dot_tensor2s_dV(const ArrD &qtensor, ArrD &elemvec) const
{
  return int_gradN_dot_tensor2_dV<cppmat::tiny::cartesian::tensor2s<double,2>>(qtensor, elemvec);
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::int_gradN_dot_tensor2s_dV(const ArrD &qtensor) const
{
  return int_gradN_dot_tensor2_dV<cppmat::tiny::cartesian::tensor2s<double,2>>(qtensor);
}

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

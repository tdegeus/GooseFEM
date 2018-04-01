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

// ================================ GooseFEM::Element::Quad4::Gauss ================================

namespace Gauss {

// --------------------------------- number of integration points ----------------------------------

size_t nip()
{
  return 4;
}

// ----------------------- integration point coordinates (local coordinates) -----------------------

ArrD coordinates()
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

ArrD weights()
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

size_t nip()
{
  return 4;
}

// ----------------------- integration point coordinates (local coordinates) -----------------------

ArrD coordinates()
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

ArrD weights()
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

inline Quadrature::Quadrature(const ArrD &x, const ArrD &xi, const ArrD &w)
: m_x(x), m_w(w), m_xi(xi)
{
  // check input
  assert( m_x.ndim()   == 3      ); // shape: [nelem, nne, ndim]
  assert( m_x.shape(1) == m_nne  ); // number of nodes per element
  assert( m_x.shape(2) == m_ndim ); // number of dimensions

  // extract number of elements
  m_nelem = m_x.shape(0);

  // integration scheme
  // - default
  if ( m_w.size() == 0 and m_xi.size() == 0 )
  {
    m_nip = Gauss::nip();
    m_xi  = Gauss::coordinates();
    m_w   = Gauss::weights();
  }
  // - input
  else if ( m_w.size() > 0 and m_xi.size() > 0 )
  {
    m_nip = m_w.size();
  }
  // - unknown
  else
  {
    throw std::runtime_error("Input integration point coordinates and weights");
  }

  // check input
  assert( m_xi.ndim()   == 2      ); // shape: [nip, ndim]
  assert( m_xi.shape(0) == m_nip  ); // number of integration points
  assert( m_xi.shape(1) == m_ndim ); // number of dimensions
  assert( m_w .ndim()   == 1      ); // shape: [nip]
  assert( m_w .size()   == m_nip  ); // number of integration points

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
    m_dNxi(k,0,0) = -.25*(1.-m_xi(k,1));    m_dNxi(k,0,1) = -.25*(1.-m_xi(k,0));
    m_dNxi(k,1,0) = +.25*(1.-m_xi(k,1));    m_dNxi(k,1,1) = -.25*(1.+m_xi(k,0));
    m_dNxi(k,2,0) = +.25*(1.+m_xi(k,1));    m_dNxi(k,2,1) = +.25*(1.+m_xi(k,0));
    m_dNxi(k,3,0) = -.25*(1.+m_xi(k,1));    m_dNxi(k,3,1) = +.25*(1.-m_xi(k,0));
  }

  // compute the shape function gradients, based on "x"
  compute_dN();
}

// -------------------------------------- number of elements ---------------------------------------

inline size_t Quadrature::nelem()
{
  return m_nelem;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t Quadrature::nne()
{
  return m_nne;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t Quadrature::ndim()
{
  return m_ndim;
}

// --------------------------------- number of integration points ----------------------------------

inline size_t Quadrature::nip()
{
  return m_nip;
}

// --------------------------------------- update positions ----------------------------------------

inline void Quadrature::update_x(const ArrD &x)
{
  // check input
  assert( x.ndim()   == 3          ); // shape: [nelem, nne, ndim]
  assert( x.shape(0) == m_nelem    ); // number of elements
  assert( x.shape(1) == m_nne      ); // number of nodes per element
  assert( x.shape(2) == m_ndim     ); // number of dimensions
  assert( x.size()   == m_x.size() ); // total number of components (redundant)

  // update positions
  std::copy(x.begin(), x.end(), m_x.begin());

  // update the shape function gradients for the new "x"
  compute_dN();
}

// ------------------------ shape function gradients in global coordinates -------------------------

inline void Quadrature::compute_dN()
{
  #pragma omp parallel
  {
    // intermediate quantities and local views
    double Jdet, w, *vol;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNxi, dNx, x;
    cppmat::cartesian2d::tensor2<double> J, Jinv;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias
      x.map(&m_x(e)); // nodal positions

      // loop over integration points
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        dNxi.map(&m_dNxi(  k)); // shape function gradients (local  coordinates)
        dNx .map(&m_dNx (e,k)); // shape function gradients (global coordinates)
        vol =    &m_vol (e,k);  // volume
        w   =     m_w   (  k);  // weight factor

        // - Jacobian (loops unrolled for efficiency)
        //   J(i,j) += dNxi(m,i) * xe(m,j)
        J(0,0) = dNxi(0,0)*x(0,0) + dNxi(1,0)*x(1,0) + dNxi(2,0)*x(2,0) + dNxi(3,0)*x(3,0);
        J(0,1) = dNxi(0,0)*x(0,1) + dNxi(1,0)*x(1,1) + dNxi(2,0)*x(2,1) + dNxi(3,0)*x(3,1);
        J(1,0) = dNxi(0,1)*x(0,0) + dNxi(1,1)*x(1,0) + dNxi(2,1)*x(2,0) + dNxi(3,1)*x(3,0);
        J(1,1) = dNxi(0,1)*x(0,1) + dNxi(1,1)*x(1,1) + dNxi(2,1)*x(2,1) + dNxi(3,1)*x(3,1);

        // - determinant and inverse of the Jacobian
        Jdet = J.det();
        Jinv = J.inv();

        // - integration point volume
        (*vol) = w * Jdet;

        // - shape function gradients wrt global coordinates (loops partly unrolled for efficiency)
        //   dNx(m,i) += Jinv(i,j) * dNxi(m,j)
        for ( size_t m = 0 ; m < m_nne ; ++m )
        {
          dNx(m,0) = Jinv(0,0) * dNxi(m,0) + Jinv(0,1) * dNxi(m,1);
          dNx(m,1) = Jinv(1,0) * dNxi(m,0) + Jinv(1,1) * dNxi(m,1);
        }
      }
    }
  } // #pragma omp parallel
}

// ------------------- dyadic product "qtensor(i,j) = dNdx(m,i) * elemvec(m,j)" --------------------

template<class T>
inline ArrD Quadrature::gradN_vector(const ArrD &elemvec)
{
  // check input
  assert( elemvec.ndim()   == 3       ); // shape: [nelem, nne, ndim]
  assert( elemvec.shape(0) == m_nelem ); // number of elements
  assert( elemvec.shape(1) == m_nne   ); // number of nodes per element
  assert( elemvec.shape(2) == m_ndim  ); // number of dimensions

  // temporary tensor, to deduce size of the output
  T dummytensor;
  size_t tdim = dummytensor.size();

  // allocate output: matrix of tensors
  ArrD qtensor({m_nelem, m_nip, tdim});

  // zero-initialize (if needed)
  if ( tdim != m_ndim*m_ndim ) qtensor.setZero();

  #pragma omp parallel
  {
    // intermediate quantities and local views
    T gradu;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<const double,m_nne,m_ndim> u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias
      u.map(&elemvec(e)); // element vector (e.g. nodal displacements)

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        dNx  .map(&m_dNx  (e,k)); // shape function gradients (global coordinates)
        gradu.map(&qtensor(e,k)); // integration point tensor (e.g. deformation gradient)

        // - evaluate dyadic product (loops unrolled for efficiency)
        //   gradu(i,j) += dNx(m,i) * ue(m,j)
        gradu(0,0) = dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        gradu(0,1) = dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1);
        gradu(1,0) = dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0);
        gradu(1,1) = dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);
      }
    }
  } // #pragma omp parallel

  return qtensor;
}

// ---------------------------------- transpose of "GradN_vector" ----------------------------------

template<class T>
inline ArrD Quadrature::gradN_vector_T(const ArrD &elemvec)
{
  // check input
  assert( elemvec.ndim()   == 3       ); // shape: [nelem, nne, ndim]
  assert( elemvec.shape(0) == m_nelem ); // number of elements
  assert( elemvec.shape(1) == m_nne   ); // number of nodes per element
  assert( elemvec.shape(2) == m_ndim  ); // number of dimensions

  // temporary tensor, to deduce size of the output
  T dummytensor;
  size_t tdim = dummytensor.size();

  // allocate output: matrix of tensors
  ArrD qtensor({m_nelem, m_nip, tdim});

  // zero-initialize (if needed)
  if ( tdim != m_ndim*m_ndim ) qtensor.setZero();

  #pragma omp parallel
  {
    // intermediate quantities and local views
    T gradu;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<const double,m_nne,m_ndim> u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias
      u.map(&elemvec(e)); // element vector (e.g. nodal displacements)

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        dNx  .map(&m_dNx  (e,k)); // shape function gradients (global coordinates)
        gradu.map(&qtensor(e,k)); // integration point tensor (e.g. deformation gradient)

        // - evaluate dyadic product (loops unrolled for efficiency)
        //   gradu(j,i) += dNx(m,i) * ue(m,j)
        gradu(0,0) = dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        gradu(1,0) = dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1);
        gradu(0,1) = dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0);
        gradu(1,1) = dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);
      }
    }
  } // #pragma omp parallel

  return qtensor;
}

// ------------------------------- symmetric part of "GradN_vector" --------------------------------

template<class T>
inline ArrD Quadrature::symGradN_vector(const ArrD &elemvec)
{
  // check input
  assert( elemvec.ndim()   == 3       ); // shape: [nelem, nne, ndim]
  assert( elemvec.shape(0) == m_nelem ); // number of elements
  assert( elemvec.shape(1) == m_nne   ); // number of nodes per element
  assert( elemvec.shape(2) == m_ndim  ); // number of dimensions

  // temporary tensor, to deduce size of the output
  T dummytensor;
  size_t tdim = dummytensor.size();

  // allocate output: matrix of tensors
  ArrD qtensor({m_nelem, m_nip, tdim});

  // zero-initialize (if needed)
  if ( !( tdim == (m_ndim+1)*m_ndim/2 or tdim == m_ndim*m_ndim ) ) qtensor.setZero();

  #pragma omp parallel
  {
    // intermediate quantities and local views
    T eps;
    cppmat::cartesian2d::tensor2<double> gradu;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<const double,m_nne,m_ndim> u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias
      u.map(&elemvec(e)); // element vector (e.g. nodal displacements)

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        dNx.map(&m_dNx  (e,k)); // shape function gradients (global coordinates)
        eps.map(&qtensor(e,k)); // integration point tensor (e.g. strain)

        // - evaluate dyadic product (loops unrolled for efficiency)
        //   gradu(i,j) += dNx(m,i) * ue(m,j)
        gradu(0,0) = dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        gradu(0,1) = dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1);
        gradu(1,0) = dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0);
        gradu(1,1) = dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);

        // - symmetrize (loops unrolled for efficiency)
        //   eps(i,j) = .5 * ( gradu(i,j) + gradu(j,i) )
        eps(0,0) =        gradu(0,0);
        eps(0,1) = .5 * ( gradu(0,1) + gradu(1,0) );
        eps(1,0) =        eps  (0,1);
        eps(1,1) =        gradu(1,1);
      }
    }
  } // #pragma omp parallel

  return qtensor;
}

// ------- scalar product "elemmat(m*ndim+i,n*ndim+i) = N(m) * qscalar * N(n)"; for all "i" --------

inline ArrD Quadrature::int_N_scalar_NT_dV(const ArrD &qscalar)
{
  // check input
  assert( qscalar.ndim()   == 2          ); // shape: [nelem, nip]
  assert( qscalar.shape(0) == m_nelem    ); // number of elements
  assert( qscalar.shape(1) == m_nip      ); // number of integration points

  // allocate output: matrix of matrices
  ArrD elemmat({m_nelem, m_nne*m_ndim, m_nne*m_ndim});

  // zero-initialize
  elemmat.setZero();

  #pragma omp parallel
  {
    // intermediate quantities and local views
    cppmat::tiny::matrix2<double,m_nne*m_ndim,m_nne*m_ndim> M;
    cppmat::tiny::vector<double,m_nne> N;
    double rho, vol;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias
      M.map(&elemmat(e)); // element matrix (e.g. mass matrix)

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        N.map(&m_N   (  k)); // shape functions
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
    }
  } // #pragma omp parallel

  return elemmat;
}

// ------------ integral of dot product "elemvec(m,j) += dNdx(m,i) * qtensor(i,j) * dV" ------------

template<class T>
inline ArrD Quadrature::int_gradN_dot_tensor2_dV(const ArrD &qtensor)
{
  #ifndef NDEBUG
  // dummy variable
  T tmp;

  // check input
  assert( qtensor.ndim()   == 3          ); // shape: [nelem, nip, #tensor-components]
  assert( qtensor.shape(0) == m_nelem    ); // number of elements
  assert( qtensor.shape(1) == m_nip      ); // number of integration points
  assert( qtensor.shape(2) == tmp.size() ); // tensor dimensions
  #endif

  // allocate output: matrix of vectors
  ArrD elemvec({m_nelem, m_nne, m_ndim});

  // zero-initialize
  elemvec.setZero();

  #pragma omp parallel
  {
    // intermediate quantities and local views
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> f;
    double vol;
    T sig;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      // alias
      f.map(&elemvec(e)); // element vector (e.g. nodal force)

      // loop over all integration points in element "e"
      for ( size_t k = 0 ; k < m_nip ; ++k )
      {
        // - alias
        dNx.map (&m_dNx  (e,k)); // shape function gradients (global coordinates)
        sig.copy(&qtensor(e,k)); // integration point tensor (e.g. stress)
        vol = m_vol      (e,k);  // integration point volume

        // - evaluate dot product, and assemble (loops partly unrolled for efficiency)
        //   f(m,j) += dNdx(m,i) * sig(i,j) * dV;
        for ( size_t m = 0 ; m < m_nne ; ++m )
        {
          f(m,0) += ( dNx(m,0) * sig(0,0) + dNx(m,1) * sig(1,0) ) * vol;
          f(m,1) += ( dNx(m,0) * sig(0,1) + dNx(m,1) * sig(1,1) ) * vol;
        }
      }
    }
  } // #pragma omp parallel

  return elemvec;
}

// -------------------------- element integral of tensor (volume average) --------------------------

template<class T>
inline T Quadrature::int_tensor2_dV(const ArrD &qtensor, size_t e)
{
  #ifndef NDEBUG
  // dummy variable
  T dummytensor;

  // check input
  assert( qtensor.ndim()   == 3                  ); // shape: [nelem, nip, #tensor-components]
  assert( qtensor.shape(0) == m_nelem            ); // number of elements
  assert( qtensor.shape(1) == m_nip              ); // number of integration points
  assert( qtensor.shape(2) == dummytensor.size() ); // tensor dimensions
  #endif

  // intermediate quantities and local views
  double vol, VOL;
  T sig, SIG;

  // zero-initialize
  SIG.setZero();
  VOL = 0.0;

  // loop over all integration points in element "e"
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    // - alias
    sig.copy(&qtensor(e,k)); // integration point tensor (e.g. stress)
    vol  =  m_vol    (e,k);  // integration point volume

    // - add to average
    SIG += vol * sig;
    VOL += vol;
  }

  // return volume average
  return SIG / VOL;
}

// ------------------------------ integral of tensor (volume average) ------------------------------

template<class T>
inline T Quadrature::int_tensor2_dV(const ArrD &qtensor)
{
  #ifndef NDEBUG
  // dummy variable
  T tmp;

  // check input
  assert( qtensor.ndim()   == 3          ); // shape: [nelem, nip, #tensor-components]
  assert( qtensor.shape(0) == m_nelem    ); // number of elements
  assert( qtensor.shape(1) == m_nip      ); // number of integration points
  assert( qtensor.shape(2) == tmp.size() ); // tensor dimensions
  #endif

  // intermediate quantities and local views
  double vol, VOL;
  T sig, SIG;

  // zero-initialize
  SIG.setZero();
  VOL = 0.0;

  // loop over all elements
  for ( size_t e = 0 ; e < m_nelem ; ++e )
  {
    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - alias
      sig.copy(&qtensor(e,k)); // integration point tensor (e.g. stress)
      vol  =  m_vol    (e,k);  // integration point volume

      // - add to average
      SIG += vol * sig;
      VOL += vol;
    }
  }

  // return volume average
  return SIG / VOL;
}

// ---------------------- wrappers with default storage (no template needed) -----------------------

inline ArrD Quadrature::gradN_vector(const ArrD &elemvec)
{
  return gradN_vector<cppmat::cartesian2d::tensor2<double>>(elemvec);
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::gradN_vector_T(const ArrD &elemvec)
{
  return gradN_vector_T<cppmat::cartesian2d::tensor2<double>>(elemvec);
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::symGradN_vector(const ArrD &elemvec)
{
  return symGradN_vector<cppmat::cartesian2d::tensor2s<double>>(elemvec);
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::int_gradN_dot_tensor2_dV(const ArrD &qtensor)
{
  assert( qtensor.ndim() == 3 ); // shape: [nelem, nip, #tensor-components]

  if ( qtensor.shape(2) == m_ndim*m_ndim )
    return int_gradN_dot_tensor2_dV<cppmat::cartesian2d::tensor2<double>>(qtensor);
  else if ( qtensor.shape(2) == (m_ndim+1)*m_ndim/2 )
    return int_gradN_dot_tensor2_dV<cppmat::cartesian2d::tensor2s<double>>(qtensor);
  else
    throw std::runtime_error("assert: qtensor.shape(2) == 4 or qtensor.shape(2) == 3");
}

// -------------------------------------------------------------------------------------------------

inline ArrD Quadrature::int_gradN_dot_tensor2s_dV(const ArrD &qtensor)
{
  return int_gradN_dot_tensor2_dV<cppmat::cartesian2d::tensor2s<double>>(qtensor);
}

// -------------------------------------------------------------------------------------------------

inline cppmat::cartesian2d::tensor2<double> Quadrature::int_tensor2_dV(
  const ArrD &qtensor, size_t e)
{
  return int_tensor2_dV<cppmat::cartesian2d::tensor2<double>>(qtensor,e);
}

// -------------------------------------------------------------------------------------------------

inline cppmat::cartesian2d::tensor2s<double> Quadrature::int_tensor2s_dV(
  const ArrD &qtensor, size_t e)
{
  return int_tensor2_dV<cppmat::cartesian2d::tensor2s<double>>(qtensor,e);
}

// -------------------------------------------------------------------------------------------------

inline cppmat::cartesian2d::tensor2<double> Quadrature::int_tensor2_dV(const ArrD &qtensor)
{
  return int_tensor2_dV<cppmat::cartesian2d::tensor2<double>>(qtensor);
}


// -------------------------------------------------------------------------------------------------

inline cppmat::cartesian2d::tensor2s<double> Quadrature::int_tensor2s_dV(const ArrD &qtensor)
{
  return int_tensor2_dV<cppmat::cartesian2d::tensor2s<double>>(qtensor);
}

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

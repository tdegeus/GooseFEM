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

// ------------------------------------------ constructor ------------------------------------------

inline Gauss::Gauss(const ArrD &x) : m_x(x)
{
  // check input
  assert( m_x.ndim()   == 3      ); // shape of the matrix [nelem, nne, ndim]
  assert( m_x.shape(1) == m_nne  ); // number of nodes per element
  assert( m_x.shape(2) == m_ndim ); // dimensions

  // extract mesh size
  m_nelem = m_x.shape(0);

  // integration point weights
  m_w.resize({m_nip});
  // integration point coordinates in local coordinates
  m_xi.resize({m_nip,m_ndim});
  // shape function gradients in local coordinates
  m_dNxi.resize({m_nip,m_nne,m_ndim});
  // shape function gradients in global coordinates
  m_dNx.resize({m_nelem,m_nip,m_nne,m_ndim});
  // volume of each integration point
  m_vol.resize({m_nelem,m_nip});

  // set integration point coordinates in local coordinates and weights
  m_xi(0,0) = -1./std::sqrt(3.); m_xi(0,1) = -1./std::sqrt(3.); m_w(0) = 1.;
  m_xi(1,0) = +1./std::sqrt(3.); m_xi(1,1) = -1./std::sqrt(3.); m_w(1) = 1.;
  m_xi(2,0) = +1./std::sqrt(3.); m_xi(2,1) = +1./std::sqrt(3.); m_w(2) = 1.;
  m_xi(3,0) = -1./std::sqrt(3.); m_xi(3,1) = +1./std::sqrt(3.); m_w(3) = 1.;

  // set shape function gradients in local coordinates
  for ( auto k = 0 ; k < m_nip ; ++k )
  {
    m_dNxi(k,0,0) = -.25*(1.-m_xi(k,1)); m_dNxi(k,0,1) = -.25*(1.-m_xi(k,0));
    m_dNxi(k,1,0) = +.25*(1.-m_xi(k,1)); m_dNxi(k,1,1) = -.25*(1.+m_xi(k,0));
    m_dNxi(k,2,0) = +.25*(1.+m_xi(k,1)); m_dNxi(k,2,1) = +.25*(1.+m_xi(k,0));
    m_dNxi(k,3,0) = -.25*(1.+m_xi(k,1)); m_dNxi(k,3,1) = +.25*(1.-m_xi(k,0));
  }

  // compute the shape function gradients
  compute_dN();
}

// -------------------------------------- number of elements ---------------------------------------

inline size_t Gauss::nelem()
{
  return m_nelem;
}

// ---------------------------------- number of nodes per element ----------------------------------

inline size_t Gauss::nne()
{
  return m_nne;
}

// ------------------------------------- number of dimensions --------------------------------------

inline size_t Gauss::ndim()
{
  return m_ndim;
}

// --------------------------------- number of integration points ----------------------------------

inline size_t Gauss::nip()
{
  return m_nip;
}

// --------------------------------------- update positions ----------------------------------------

inline void Gauss::update_x(const ArrD &x)
{
  // check input
  assert( x.ndim()   == 3       ); // shape of the matrix [nelem, nne, ndim]
  assert( x.shape(0) == m_nelem ); // number of elements
  assert( x.shape(1) == m_nne   ); // number of nodes per element
  assert( x.shape(2) == m_ndim  ); // dimensions

  // update positions
  std::copy(x.begin(), x.end(), m_x.begin());

  // update the shape function gradients
  compute_dN();
}

// ------------------------ shape function gradients in global coordinates -------------------------

inline void Gauss::compute_dN()
{
  #pragma omp parallel
  {
    // intermediate quantities
    cppmat::cartesian2d::tensor2<double> J, Jinv;
    double Jdet;
    // local views of the global arrays (speeds up indexing, and increases readability)
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNxi, dNx, x;
    cppmat::tiny::vector<double,m_nne> w, vol;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( auto e = 0 ; e < m_nelem ; ++e )
    {
      // pointers to element positions, element integration volume, and weight
      x  .map(&m_x  (e));
      vol.map(&m_vol(e));
      w  .map(&m_w  (0));

      // loop over Gauss points
      for ( auto k = 0 ; k < m_nip ; ++k )
      {
        // - pointer to the shape function gradients (local/global coordinates)
        dNxi.map(&m_dNxi(  k));
        dNx .map(&m_dNx (e,k));

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
        vol(k) = w(k) * Jdet;

        // - shape function gradients wrt global coordinates (loops partly unrolled for efficiency)
        //   dNx(m,i) += Jinv(i,j) * dNxi(m,j)
        for ( auto m = 0 ; m < m_nne ; ++m )
        {
          dNx(m,0) = Jinv(0,0) * dNxi(m,0) + Jinv(0,1) * dNxi(m,1);
          dNx(m,1) = Jinv(1,0) * dNxi(m,0) + Jinv(1,1) * dNxi(m,1);
        }
      }
    }
  } // #pragma omp parallel
}

// ----------------------- dyadic product "out(i,j) = dNdx(m,i) * inp(m,j)" ------------------------

template<class T>
inline ArrD Gauss::gradN_vector(const ArrD &inp)
{
  // check input
  assert( inp.ndim()   == 3       ); // shape of the matrix [nelem,nne,ndim]
  assert( inp.shape(0) == m_nelem ); // number of elements
  assert( inp.shape(1) == m_nne   ); // number of nodes per element
  assert( inp.shape(2) == m_ndim  ); // dimensions

  // temporary tensor, to deduce size of the output
  T tmp;

  // allocate output: matrix of tensors
  ArrD out({m_nelem, m_nip, tmp.size()});

  // zero-initialize
  if ( tmp.size() != m_ndim*m_ndim ) out.setZero();

  #pragma omp parallel
  {
    // local views of the global arrays (speeds up indexing, and increases readability)
    T gradu;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<const double,m_nne,m_ndim> u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( auto e = 0 ; e < m_nelem ; ++e )
    {
      // pointer to element vector
      // (e.g. nodal displacements of each element)
      u.map(&inp(e));

      // loop over all integration points in element "e"
      for ( auto k = 0 ; k < m_nip ; ++k )
      {
        // - pointer to the shape function gradients and integration point tensor
        //   (e.g. integration point deformation gradient)
        dNx  .map(&m_dNx(e,k));
        gradu.map(&out  (e,k));

        // - evaluate dyadic product (loops unrolled for efficiency)
        //   (e.g. integration point deformation gradient)
        //   gradu(i,j) += dNx(m,i) * ue(m,j)
        gradu(0,0) = dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        gradu(0,1) = dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1);
        gradu(1,0) = dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0);
        gradu(1,1) = dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);
      }
    }
  } // #pragma omp parallel

  return out;
}

// ----------------------- wrapper with default storage (no template needed) -----------------------

inline ArrD Gauss::gradN_vector(const ArrD &inp)
{
  return gradN_vector<cppmat::cartesian2d::tensor2<double>>(inp);
}

// ----------------------- dyadic product "out(j,i) = dNdx(m,i) * inp(m,j)" ------------------------

template<class T>
inline ArrD Gauss::vector_GradN(const ArrD &inp)
{
  // check input
  assert( inp.ndim()   == 3       ); // shape of the matrix [nelem,nne,ndim]
  assert( inp.shape(0) == m_nelem ); // number of elements
  assert( inp.shape(1) == m_nne   ); // number of nodes per element
  assert( inp.shape(2) == m_ndim  ); // dimensions

  // temporary tensor, to deduce size of the output
  T tmp;

  // allocate output: matrix of tensors
  ArrD out({m_nelem, m_nip, tmp.size()});

  // zero-initialize
  if ( tmp.size() != m_ndim*m_ndim ) out.setZero();

  #pragma omp parallel
  {
    // local views of the global arrays (speeds up indexing, and increases readability)
    T gradu;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<const double,m_nne,m_ndim> u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( auto e = 0 ; e < m_nelem ; ++e )
    {
      // pointer to element vector
      // (e.g. nodal displacements of each element)
      u.map(&inp(e));

      // loop over all integration points in element "e"
      for ( auto k = 0 ; k < m_nip ; ++k )
      {
        // - pointer to the shape function gradients and integration point tensor
        //   (e.g. integration point deformation gradient)
        dNx  .map(&m_dNx(e,k));
        gradu.map(&out  (e,k));

        // - evaluate dyadic product (loops unrolled for efficiency)
        //   (e.g. integration point deformation gradient)
        //   gradu(j,i) += dNx(m,i) * ue(m,j)
        gradu(0,0) = dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        gradu(1,0) = dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1);
        gradu(0,1) = dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0);
        gradu(1,1) = dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);
      }
    }
  } // #pragma omp parallel

  return out;
}

// ----------------------- wrapper with default storage (no template needed) -----------------------

inline ArrD Gauss::vector_GradN(const ArrD &inp)
{
  return vector_GradN<cppmat::cartesian2d::tensor2<double>>(inp);
}

// ------------------ symmetric dyadic product "out(i,j) = dNdx(m,i) * inp(m,j)" -------------------

template<class T>
inline ArrD Gauss::symGradN_vector(const ArrD &inp)
{
  // check input
  assert( inp.ndim()   == 3       ); // shape of the matrix [nelem,nne,ndim]
  assert( inp.shape(0) == m_nelem ); // number of elements
  assert( inp.shape(1) == m_nne   ); // number of nodes per element
  assert( inp.shape(2) == m_ndim  ); // dimensions

  // temporary tensor, to deduce size of the output
  T tmp;

  // allocate output: matrix of tensors
  ArrD out({m_nelem, m_nip, tmp.size()});

  // zero-initialize
  if ( !( tmp.size() == (m_ndim+1)*m_ndim/2 or tmp.size() == m_ndim ) ) out.setZero();

  #pragma omp parallel
  {
    // local views of the global arrays (speeds up indexing, and increases readability)
    T eps;
    cppmat::cartesian2d::tensor2<double> gradu;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<const double,m_nne,m_ndim> u;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( auto e = 0 ; e < m_nelem ; ++e )
    {
      // pointer to element vector
      // (e.g. nodal displacements of each element)
      u.map(&inp(e));

      // loop over all integration points in element "e"
      for ( auto k = 0 ; k < m_nip ; ++k )
      {
        // - pointer to the shape function gradients and integration point tensor
        //   (e.g. integration point strain)
        dNx.map(&m_dNx(e,k));
        eps.map(&out  (e,k));

        // - evaluate dyadic product (loops unrolled for efficiency)
        //   (e.g. integration point deformation gradient)
        //   gradu(i,j) += dNx(m,i) * ue(m,j)
        gradu(0,0) = dNx(0,0)*u(0,0) + dNx(1,0)*u(1,0) + dNx(2,0)*u(2,0) + dNx(3,0)*u(3,0);
        gradu(0,1) = dNx(0,0)*u(0,1) + dNx(1,0)*u(1,1) + dNx(2,0)*u(2,1) + dNx(3,0)*u(3,1);
        gradu(1,0) = dNx(0,1)*u(0,0) + dNx(1,1)*u(1,0) + dNx(2,1)*u(2,0) + dNx(3,1)*u(3,0);
        gradu(1,1) = dNx(0,1)*u(0,1) + dNx(1,1)*u(1,1) + dNx(2,1)*u(2,1) + dNx(3,1)*u(3,1);

        // - symmetrize, store symmetric (loops unrolled for efficiency)
        //   (e.g. strain)
        //   eps(i,j) = .5 * ( gradu(i,j) + gradu(j,i) )
        eps(0,0) =        gradu(0,0);
        eps(0,1) = .5 * ( gradu(0,1) + gradu(1,0) );
        eps(1,0) =        eps  (0,1);
        eps(1,1) =        gradu(1,1);
      }
    }
  } // #pragma omp parallel

  return out;
}

// ----------------------- wrapper with default storage (no template needed) -----------------------

inline ArrD Gauss::symGradN_vector(const ArrD &inp)
{
  return symGradN_vector<cppmat::cartesian2d::tensor2s<double>>(inp);
}

// ------------------------- dot product "out(m,j) = dNdx(m,i) * inp(i,j)" -------------------------

template<class T>
inline ArrD Gauss::gradN_dot_tensor2(const ArrD &inp)
{
  #ifndef NDEBUG
  // dummy variable
  T tmp;

  // check input
  assert( inp.ndim()   == 3          ); // shape of the matrix [nelem,nip,tensor-dim]
  assert( inp.shape(0) == m_nelem    ); // number of elements
  assert( inp.shape(1) == m_nip      ); // number of integration points
  assert( inp.shape(2) == tmp.size() ); // tensor dimensions
  #endif

  // allocate output: matrix of vectors
  ArrD out({m_nelem, m_nne, m_ndim});

  #pragma omp parallel
  {
    // local views of the global arrays (speeds up indexing, and increases readability)
    T sig;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> f;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( auto e = 0 ; e < m_nelem ; ++e )
    {
      // pointer to element vector
      // (e.g. nodal force of each element)
      f.map(&out(e));

      // zero-initialize
      f.setZero();

      // loop over all integration points in element "e"
      for ( auto k = 0 ; k < m_nip ; ++k )
      {
        // - pointer/copy to the shape function gradients and integration point tensor
        dNx.map (&m_dNx(e,k));
        sig.copy(&inp  (e,k));

        // - assemble to element vector
        //   (e.g. element force)
        //   f(m,j) += dNdx(m,i) * sig(i,j);
        for ( size_t m = 0 ; m < m_nne ; ++m )
        {
          f(m,0) += dNx(m,0) * sig(0,0) + dNx(m,1) * sig(1,0);
          f(m,1) += dNx(m,0) * sig(0,1) + dNx(m,1) * sig(1,1);
        }
      }
    }
  } // #pragma omp parallel

  return out;
}

// ----------------------- wrapper with default storage (no template needed) -----------------------

inline ArrD Gauss::gradN_dot_tensor2(const ArrD &inp)
{
  // check input
  assert( inp.ndim() == 3 ); // shape of the matrix [nelem,nip,tensor-dim]

  // general tensor
  if ( inp.shape(2) == m_ndim*m_ndim )
    return gradN_dot_tensor2<cppmat::cartesian2d::tensor2<double>>(inp);
  // symmetric tensor
  else if ( inp.shape(2) == (m_ndim+1)*m_ndim/2 )
    return gradN_dot_tensor2<cppmat::cartesian2d::tensor2s<double>>(inp);
  // unknown input
  else
    throw std::runtime_error("assert: inp.shape(2) == 4 or inp.shape(2) == 3");
}

// ----------------------- wrapper with default storage (no template needed) -----------------------

inline ArrD Gauss::gradN_dot_tensor2s(const ArrD &inp)
{
  return gradN_dot_tensor2<cppmat::cartesian2d::tensor2s<double>>(inp);
}

// ---------------------- dot product "out(m,j) = dNdx(m,i) * inp(i,j) * dV" -----------------------

template<class T>
inline ArrD Gauss::gradN_dot_tensor2_dV(const ArrD &inp)
{
  #ifndef NDEBUG
  // dummy variable
  T tmp;

  // check input
  assert( inp.ndim()   == 3          ); // shape of the matrix [nelem,nip,tensor-dim]
  assert( inp.shape(0) == m_nelem    ); // number of elements
  assert( inp.shape(1) == m_nip      ); // number of integration points
  assert( inp.shape(2) == tmp.size() ); // tensor dimensions
  #endif

  // allocate output: matrix of vectors
  ArrD out({m_nelem, m_nne, m_ndim});

  #pragma omp parallel
  {
    // local views of the global arrays (speeds up indexing, and increases readability)
    T sig;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> dNx;
    cppmat::tiny::matrix2<double,m_nne,m_ndim> f;
    cppmat::tiny::vector<double,m_nne> vol;

    // loop over all elements (in parallel)
    #pragma omp for
    for ( auto e = 0 ; e < m_nelem ; ++e )
    {
      // pointer to element vector, integration volume
      // (e.g. nodal force of each element)
      f  .map(&out  (e));
      vol.map(&m_vol(e));

      // zero-initialize
      f.setZero();

      // loop over all integration points in element "e"
      for ( auto k = 0 ; k < m_nip ; ++k )
      {
        // - pointer/copy to the shape function gradients and integration point tensor
        dNx.map (&m_dNx(e,k));
        sig.copy(&inp  (e,k));

        // - assemble to element vector
        //   (e.g. element force)
        //   f(m,j) += dNdx(m,i) * sig(i,j);
        for ( size_t m = 0 ; m < m_nne ; ++m )
        {
          f(m,0) += dNx(m,0) * sig(0,0) * vol(k) + dNx(m,1) * sig(1,0) * vol(k);
          f(m,1) += dNx(m,0) * sig(0,1) * vol(k) + dNx(m,1) * sig(1,1) * vol(k);
        }
      }
    }
  } // #pragma omp parallel

  return out;
}

// ----------------------- wrapper with default storage (no template needed) -----------------------

inline ArrD Gauss::gradN_dot_tensor2_dV(const ArrD &inp)
{
  // check input
  assert( inp.ndim() == 3 ); // shape of the matrix [nelem,nip,tensor-dim]

  // general tensor
  if ( inp.shape(2) == m_ndim*m_ndim )
    return gradN_dot_tensor2_dV<cppmat::cartesian2d::tensor2<double>>(inp);
  // symmetric tensor
  else if ( inp.shape(2) == (m_ndim+1)*m_ndim/2 )
    return gradN_dot_tensor2_dV<cppmat::cartesian2d::tensor2s<double>>(inp);
  // unknown input
  else
    throw std::runtime_error("assert: inp.shape(2) == 4 or inp.shape(2) == 3");
}

// ------------------------------------ element volume average -------------------------------------

template<class T>
inline T Gauss::average_tensor2(const ArrD &inp, size_t e)
{
  T sig, SIG;
  cppmat::tiny::vector<double,m_nne> vol;
  double VOL;

  // zero-initialize
  SIG.setZero();
  VOL = 0.0;

  // pointer to integration volume
  vol.map(&m_vol(e));

  // loop over all integration points in element "e"
  for ( size_t k = 0 ; k < m_nip ; ++k )
  {
    // - copy element tensor (copy is needed because of the input that is marker "const")
    //   (e.g. stress)
    sig.copy(&inp(e,k));

    // - add to average
    SIG += vol(k) * sig;
    VOL += vol(k);
  }

  // return volume average
  return SIG / VOL;
}

// ----------------------- wrapper with default storage (no template needed) -----------------------

inline cppmat::cartesian2d::tensor2<double> Gauss::average_tensor2(const ArrD &inp, size_t e)
{
  return average_tensor2<cppmat::cartesian2d::tensor2<double>>(inp,e);
}


// ----------------------- wrapper with default storage (no template needed) -----------------------

inline cppmat::cartesian2d::tensor2s<double> Gauss::average_tensor2s(const ArrD &inp, size_t e)
{
  return average_tensor2<cppmat::cartesian2d::tensor2s<double>>(inp,e);
}

// ---------------------------------------- volume average -----------------------------------------

template<class T>
inline T Gauss::average_tensor2(const ArrD &inp)
{
  T sig, SIG;
  cppmat::tiny::vector<double,m_nne> vol;
  double VOL;

  // zero-initialize
  SIG.setZero();
  VOL = 0.0;

  // loop over all elements
  for ( auto e = 0 ; e < m_nelem ; ++e )
  {
    // pointer to integration volume
    vol.map(&m_vol(e));

    // loop over all integration points in element "e"
    for ( size_t k = 0 ; k < m_nip ; ++k )
    {
      // - copy element tensor (copy is needed because of the input that is marker "const")
      //   (e.g. stress)
      sig.copy(&inp(e,k));

      // - add to average
      SIG += vol(k) * sig;
      VOL += vol(k);
    }
  }

  // return volume average
  return SIG / VOL;
}

// ----------------------- wrapper with default storage (no template needed) -----------------------

inline cppmat::cartesian2d::tensor2<double> Gauss::average_tensor2(const ArrD &inp)
{
  return average_tensor2<cppmat::cartesian2d::tensor2<double>>(inp);
}


// ----------------------- wrapper with default storage (no template needed) -----------------------

inline cppmat::cartesian2d::tensor2s<double> Gauss::average_tensor2s(const ArrD &inp)
{
  return average_tensor2<cppmat::cartesian2d::tensor2s<double>>(inp);
}

// -------------------------------------------------------------------------------------------------

}}} // namespace ...

// =================================================================================================

#endif

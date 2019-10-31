/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#ifndef GOOSEFEM_MATRIXPARTITIONEDTYINGS_HPP
#define GOOSEFEM_MATRIXPARTITIONEDTYINGS_HPP

// -------------------------------------------------------------------------------------------------

#include "MatrixPartitionedTyings.h"

// =================================================================================================

namespace GooseFEM {

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline MatrixPartitionedTyings<Solver>::MatrixPartitionedTyings(
  const xt::xtensor<size_t,2>& conn,
  const xt::xtensor<size_t,2>& dofs,
  const Eigen::SparseMatrix<double>& Cdu,
  const Eigen::SparseMatrix<double>& Cdp) :
  m_conn(conn),
  m_dofs(dofs),
  m_Cdu(Cdu),
  m_Cdp(Cdp)
{
  GOOSEFEM_ASSERT(Cdu.rows() == Cdp.rows());

  m_nnu = static_cast<size_t>(m_Cdu.cols());
  m_nnp = static_cast<size_t>(m_Cdp.cols());
  m_nnd = static_cast<size_t>(m_Cdp.rows());
  m_nni = m_nnu + m_nnp;
  m_ndof = m_nni + m_nnd;

  m_iiu = xt::arange<size_t>(m_nnu);
  m_iip = xt::arange<size_t>(m_nnu, m_nnu + m_nnp);
  m_iid = xt::arange<size_t>(m_nni, m_nni + m_nnd);

  m_nelem = m_conn.shape(0);
  m_nne = m_conn.shape(1);
  m_nnode = m_dofs.shape(0);
  m_ndim = m_dofs.shape(1);

  m_Cud = m_Cdu.transpose();
  m_Cpd = m_Cdp.transpose();

  m_Tuu.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tup.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tpu.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tpp.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tud.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tpd.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tdu.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tdp.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);
  m_Tdd.reserve(m_nelem*m_nne*m_ndim*m_nne*m_ndim);

  m_Auu.resize(m_nnu,m_nnu);
  m_Aup.resize(m_nnu,m_nnp);
  m_Apu.resize(m_nnp,m_nnu);
  m_App.resize(m_nnp,m_nnp);
  m_Aud.resize(m_nnu,m_nnd);
  m_Apd.resize(m_nnp,m_nnd);
  m_Adu.resize(m_nnd,m_nnu);
  m_Adp.resize(m_nnd,m_nnp);
  m_Add.resize(m_nnd,m_nnd);

  GOOSEFEM_ASSERT(m_ndof <= m_nnode * m_ndim);
  GOOSEFEM_ASSERT(m_ndof == xt::amax(m_dofs)[0] + 1);
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::nelem() const
{
  return m_nelem;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::nne() const
{
  return m_nne;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::nnode() const
{
  return m_nnode;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::ndof() const
{
  return m_ndof;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::nnu() const
{
  return m_nnu;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::nnp() const
{
  return m_nnp;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::nni() const
{
  return m_nni;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline size_t MatrixPartitionedTyings<Solver>::nnd() const
{
  return m_nnd;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<size_t,2> MatrixPartitionedTyings<Solver>::dofs() const
{
  return m_dofs;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<size_t,1> MatrixPartitionedTyings<Solver>::iiu() const
{
  return m_iiu;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<size_t,1> MatrixPartitionedTyings<Solver>::iip() const
{
  return m_iip;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<size_t,1> MatrixPartitionedTyings<Solver>::iii() const
{
  return xt::arange<size_t>(m_nni);
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<size_t,1> MatrixPartitionedTyings<Solver>::iid() const
{
  return m_iid;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitionedTyings<Solver>::factorize()
{
  if (!m_factor) return;

  m_ACuu = m_Auu + m_Aud * m_Cdu + m_Cud * m_Adu + m_Cud * m_Add * m_Cdu;
  m_ACup = m_Aup + m_Aud * m_Cdp + m_Cud * m_Adp + m_Cud * m_Add * m_Cdp;
  // m_ACpu = m_Apu + m_Apd * m_Cdu + m_Cpd * m_Adu + m_Cpd * m_Add * m_Cdu;
  // m_ACpp = m_App + m_Apd * m_Cdp + m_Cpd * m_Adp + m_Cpd * m_Add * m_Cdp;

  m_solver.compute(m_ACuu);

  m_factor = false;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitionedTyings<Solver>::assemble(const xt::xtensor<double,3>& elemmat)
{
  GOOSEFEM_ASSERT(elemmat.shape() ==\
    std::decay_t<decltype(elemmat)>::shape_type({m_nelem, m_nne*m_ndim, m_nne*m_ndim}));

  m_Tuu.clear();
  m_Tup.clear();
  m_Tpu.clear();
  m_Tpp.clear();
  m_Tud.clear();
  m_Tpd.clear();
  m_Tdu.clear();
  m_Tdp.clear();
  m_Tdd.clear();

  for (size_t e = 0 ; e < m_nelem ; ++e) {
    for (size_t m = 0 ; m < m_nne ; ++m) {
      for (size_t i = 0 ; i < m_ndim ; ++i) {

        size_t di = m_dofs(m_conn(e,m),i);

        for (size_t n = 0 ; n < m_nne ; ++n) {
          for (size_t j = 0 ; j < m_ndim ; ++j) {

            size_t dj = m_dofs(m_conn(e,n),j);

            if (di < m_nnu && dj < m_nnu) {
              m_Tuu.push_back(Eigen::Triplet<double>(
                di,
                dj,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
            else if (di < m_nnu && dj < m_nni) {
              m_Tup.push_back(Eigen::Triplet<double>(
                di,
                dj - m_nnu,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
            else if (di < m_nnu) {
              m_Tud.push_back(Eigen::Triplet<double>(
                di,
                dj - m_nni,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
            else if (di < m_nni && dj < m_nnu) {
              m_Tpu.push_back(Eigen::Triplet<double>(
                di - m_nnu,
                dj,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
            else if (di < m_nni && dj < m_nni) {
              m_Tpp.push_back(Eigen::Triplet<double>(
                di - m_nnu,
                dj - m_nnu,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
            else if (di < m_nni) {
              m_Tpd.push_back(Eigen::Triplet<double>(
                di - m_nnu,
                dj - m_nni,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
            else if (dj < m_nnu) {
              m_Tdu.push_back(Eigen::Triplet<double>(
                di - m_nni,
                dj,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
            else if (dj < m_nni) {
              m_Tdp.push_back(Eigen::Triplet<double>(
                di - m_nni,
                dj - m_nnu,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
            else {
              m_Tdd.push_back(Eigen::Triplet<double>(
                di - m_nni,
                dj - m_nni,
                elemmat(e, m * m_ndim + i, n * m_ndim + j)));
            }
          }
        }
      }
    }
  }

  m_Auu.setFromTriplets(m_Tuu.begin(), m_Tuu.end());
  m_Aup.setFromTriplets(m_Tup.begin(), m_Tup.end());
  m_Apu.setFromTriplets(m_Tpu.begin(), m_Tpu.end());
  m_App.setFromTriplets(m_Tpp.begin(), m_Tpp.end());
  m_Aud.setFromTriplets(m_Tud.begin(), m_Tud.end());
  m_Apd.setFromTriplets(m_Tpd.begin(), m_Tpd.end());
  m_Adu.setFromTriplets(m_Tdu.begin(), m_Tdu.end());
  m_Adp.setFromTriplets(m_Tdp.begin(), m_Tdp.end());
  m_Add.setFromTriplets(m_Tdd.begin(), m_Tdd.end());

  m_factor = true;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitionedTyings<Solver>::solve(
  const xt::xtensor<double,2>& b,
        xt::xtensor<double,2>& x)
{
  GOOSEFEM_ASSERT(b.shape() ==\
    std::decay_t<decltype(b)>::shape_type({m_nnode, m_ndim}));
  GOOSEFEM_ASSERT(x.shape() ==\
    std::decay_t<decltype(x)>::shape_type({m_nnode, m_ndim}));

  this->factorize();

  Eigen::VectorXd B_u = this->asDofs_u(b);
  Eigen::VectorXd B_d = this->asDofs_d(b);
  Eigen::VectorXd X_p = this->asDofs_p(x);

  B_u += m_Cud * B_d;

  Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - m_ACup * X_p));
  Eigen::VectorXd X_d = m_Cdu * X_u + m_Cdp * X_p;

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m) {
    for (size_t i = 0 ; i < m_ndim ; ++i) {
      if (m_dofs(m,i) < m_nnu)
        x(m,i) = X_u(m_dofs(m,i));
      else if (m_dofs(m,i) >= m_nni)
        x(m,i) = X_d(m_dofs(m,i)-m_nni);
    }
  }
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitionedTyings<Solver>::solve(
  const xt::xtensor<double,1>& b,
        xt::xtensor<double,1>& x)
{
  GOOSEFEM_ASSERT(b.size() == m_ndof);
  GOOSEFEM_ASSERT(x.size() == m_ndof);

  this->factorize();

  Eigen::VectorXd B_u = this->asDofs_u(b);
  Eigen::VectorXd B_d = this->asDofs_d(b);
  Eigen::VectorXd X_p = this->asDofs_p(x);

  Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - m_ACup * X_p));
  Eigen::VectorXd X_d = m_Cdu * X_u + m_Cdp * X_p;

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnu ; ++d)
    x(m_iiu(d)) = X_u(d);

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnd ; ++d)
    x(m_iid(d)) = X_d(d);
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline void MatrixPartitionedTyings<Solver>::solve_u(
  const xt::xtensor<double,1>& b_u,
  const xt::xtensor<double,1>& b_d,
  const xt::xtensor<double,1>& x_p,
        xt::xtensor<double,1>& x_u)
{
  GOOSEFEM_ASSERT(b_u.size() == m_nnu);
  GOOSEFEM_ASSERT(b_d.size() == m_nnd);
  GOOSEFEM_ASSERT(x_p.size() == m_nnp);
  GOOSEFEM_ASSERT(x_u.size() == m_nnu);

  this->factorize();

  Eigen::VectorXd B_u(m_nnu,1);
  Eigen::VectorXd B_d(m_nnd,1);
  Eigen::VectorXd X_p(m_nnp,1);

  std::copy(b_u.begin(), b_u.end(), B_u.data());
  std::copy(b_d.begin(), b_d.end(), B_d.data());
  std::copy(x_p.begin(), x_p.end(), X_p.data());

  Eigen::VectorXd X_u = m_solver.solve(Eigen::VectorXd(B_u - m_ACup * X_p));

  std::copy(X_u.data(), X_u.data()+m_nnu, x_u.begin());
}

// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,2> MatrixPartitionedTyings<Solver>::Solve(
  const xt::xtensor<double,2> &b,
  const xt::xtensor<double,2> &x)
{
  xt::xtensor<double,2> out = x;
  this->solve(b, out);
  return out;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,1> MatrixPartitionedTyings<Solver>::Solve(
  const xt::xtensor<double,1> &b,
  const xt::xtensor<double,1> &x)
{
  xt::xtensor<double,1> out = x;
  this->solve(b, out);
  return out;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline xt::xtensor<double,1> MatrixPartitionedTyings<Solver>::Solve_u(
  const xt::xtensor<double,1> &b_u,
  const xt::xtensor<double,1> &b_d,
  const xt::xtensor<double,1> &x_p)
{
  xt::xtensor<double,1> x_u = xt::empty<double>({m_nnu});
  this->solve_u(b_u, b_d, x_p, x_u);
  return x_u;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitionedTyings<Solver>::asDofs_u(
  const xt::xtensor<double,1>& dofval) const
{
  assert(dofval.size() == m_ndof);

  Eigen::VectorXd dofval_u(m_nnu,1);

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnu ; ++d)
    dofval_u(d) = dofval(m_iiu(d));

  return dofval_u;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitionedTyings<Solver>::asDofs_u(
  const xt::xtensor<double,2>& nodevec) const
{
  assert(nodevec.shape() ==\
    std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

  Eigen::VectorXd dofval_u(m_nnu,1);

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      if (m_dofs(m,i) < m_nnu)
        dofval_u(m_dofs(m,i)) = nodevec(m,i);

  return dofval_u;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitionedTyings<Solver>::asDofs_p(
  const xt::xtensor<double,1>& dofval) const
{
  assert(dofval.size() == m_ndof);

  Eigen::VectorXd dofval_p(m_nnp,1);

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnp ; ++d)
    dofval_p(d) = dofval(m_iip(d));

  return dofval_p;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitionedTyings<Solver>::asDofs_p(
  const xt::xtensor<double,2>& nodevec) const
{
  assert(nodevec.shape() ==\
    std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

  Eigen::VectorXd dofval_p(m_nnp,1);

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      if (m_dofs(m,i) >= m_nnu && m_dofs(m,i) < m_nni)
        dofval_p(m_dofs(m,i)-m_nnu) = nodevec(m,i);

  return dofval_p;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitionedTyings<Solver>::asDofs_d(
  const xt::xtensor<double,1>& dofval) const
{
  assert(dofval.size() == m_ndof);

  Eigen::VectorXd dofval_d(m_nnd,1);

  #pragma omp parallel for
  for (size_t d = 0 ; d < m_nnd ; ++d)
    dofval_d(d) = dofval(m_iip(d));

  return dofval_d;
}

// -------------------------------------------------------------------------------------------------

template <class Solver>
inline Eigen::VectorXd MatrixPartitionedTyings<Solver>::asDofs_d(
  const xt::xtensor<double,2>& nodevec) const
{
  assert(nodevec.shape() ==\
    std::decay_t<decltype(nodevec)>::shape_type({m_nnode, m_ndim}));

  Eigen::VectorXd dofval_d(m_nnd,1);

  #pragma omp parallel for
  for (size_t m = 0 ; m < m_nnode ; ++m)
    for (size_t i = 0 ; i < m_ndim ; ++i)
      if (m_dofs(m,i) >= m_nni)
        dofval_d(m_dofs(m,i)-m_nni) = nodevec(m,i);

  return dofval_d;
}

// -------------------------------------------------------------------------------------------------

} // namespace ...

// =================================================================================================

#endif

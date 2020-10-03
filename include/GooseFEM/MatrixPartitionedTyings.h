/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_MATRIXPARTITIONEDTYINGS_H
#define GOOSEFEM_MATRIXPARTITIONEDTYINGS_H

#include "config.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace GooseFEM {

// forward declaration
template <class> class MatrixPartitionedTyingsSolver;

class MatrixPartitionedTyings {
public:
    // Constructors
    MatrixPartitionedTyings() = default;

    MatrixPartitionedTyings(
        const xt::xtensor<size_t, 2>& conn,
        const xt::xtensor<size_t, 2>& dofs,
        const Eigen::SparseMatrix<double>& Cdu,
        const Eigen::SparseMatrix<double>& Cdp);

    // Dimensions
    size_t nelem() const; // number of elements
    size_t nne() const;   // number of nodes per element
    size_t nnode() const; // number of nodes
    size_t ndim() const;  // number of dimensions
    size_t ndof() const;  // number of DOFs
    size_t nnu() const;   // number of independent, unknown DOFs
    size_t nnp() const;   // number of independent, prescribed DOFs
    size_t nni() const;   // number of independent DOFs
    size_t nnd() const;   // number of dependent DOFs

    // DOF lists
    xt::xtensor<size_t, 2> dofs() const; // DOFs
    xt::xtensor<size_t, 1> iiu() const;  // independent, unknown DOFs
    xt::xtensor<size_t, 1> iip() const;  // independent, prescribed DOFs
    xt::xtensor<size_t, 1> iii() const;  // independent DOFs
    xt::xtensor<size_t, 1> iid() const;  // dependent DOFs

    // Assemble from matrices stored per element [nelem, nne*ndim, nne*ndim]
    void assemble(const xt::xtensor<double, 3>& elemmat);

private:
    // The matrix
    Eigen::SparseMatrix<double> m_Auu;
    Eigen::SparseMatrix<double> m_Aup;
    Eigen::SparseMatrix<double> m_Apu;
    Eigen::SparseMatrix<double> m_App;
    Eigen::SparseMatrix<double> m_Aud;
    Eigen::SparseMatrix<double> m_Apd;
    Eigen::SparseMatrix<double> m_Adu;
    Eigen::SparseMatrix<double> m_Adp;
    Eigen::SparseMatrix<double> m_Add;

    // The matrix for which the tyings have been applied
    Eigen::SparseMatrix<double> m_ACuu;
    Eigen::SparseMatrix<double> m_ACup;
    Eigen::SparseMatrix<double> m_ACpu;
    Eigen::SparseMatrix<double> m_ACpp;

    // Matrix entries
    std::vector<Eigen::Triplet<double>> m_Tuu;
    std::vector<Eigen::Triplet<double>> m_Tup;
    std::vector<Eigen::Triplet<double>> m_Tpu;
    std::vector<Eigen::Triplet<double>> m_Tpp;
    std::vector<Eigen::Triplet<double>> m_Tud;
    std::vector<Eigen::Triplet<double>> m_Tpd;
    std::vector<Eigen::Triplet<double>> m_Tdu;
    std::vector<Eigen::Triplet<double>> m_Tdp;
    std::vector<Eigen::Triplet<double>> m_Tdd;

    // Signal changes to data
    bool m_changed = true;

    // Bookkeeping
    xt::xtensor<size_t, 2> m_conn; // connectivity          [nelem, nne ]
    xt::xtensor<size_t, 2> m_dofs; // DOF-numbers per node  [nnode, ndim]
    xt::xtensor<size_t, 1> m_iiu;  // unknown     DOFs      [nnu]
    xt::xtensor<size_t, 1> m_iip;  // prescribed  DOFs      [nnp]
    xt::xtensor<size_t, 1> m_iid;  // dependent   DOFs      [nnd]

    // Dimensions
    size_t m_nelem; // number of elements
    size_t m_nne;   // number of nodes per element
    size_t m_nnode; // number of nodes
    size_t m_ndim;  // number of dimensions
    size_t m_ndof;  // number of DOFs
    size_t m_nnu;   // number of independent, unknown DOFs
    size_t m_nnp;   // number of independent, prescribed DOFs
    size_t m_nni;   // number of independent DOFs
    size_t m_nnd;   // number of dependent DOFs

    // Tyings
    Eigen::SparseMatrix<double> m_Cdu;
    Eigen::SparseMatrix<double> m_Cdp;
    Eigen::SparseMatrix<double> m_Cud;
    Eigen::SparseMatrix<double> m_Cpd;

    // grant access to solver class
    template <class> friend class MatrixPartitionedTyingsSolver;

    // Convert arrays (Eigen version of VectorPartitioned, which contains public functions)
    Eigen::VectorXd asDofs_u(const xt::xtensor<double, 1>& dofval) const;
    Eigen::VectorXd asDofs_u(const xt::xtensor<double, 2>& nodevec) const;
    Eigen::VectorXd asDofs_p(const xt::xtensor<double, 1>& dofval) const;
    Eigen::VectorXd asDofs_p(const xt::xtensor<double, 2>& nodevec) const;
    Eigen::VectorXd asDofs_d(const xt::xtensor<double, 1>& dofval) const;
    Eigen::VectorXd asDofs_d(const xt::xtensor<double, 2>& nodevec) const;
};

template <class Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>
class MatrixPartitionedTyingsSolver {
public:
    // Constructors
    MatrixPartitionedTyingsSolver() = default;

    // Solve:
    // A' = A_ii + K_id * C_di + C_di^T * K_di + C_di^T * K_dd * C_di
    // b' = b_i + C_di^T * b_d
    // x_u = A'_uu \ ( b'_u - A'_up * x_p )
    // x_i = [x_u, x_p]
    // x_d = C_di * x_i
    void solve(
        MatrixPartitionedTyings& matrix,
        const xt::xtensor<double, 2>& b,
        xt::xtensor<double, 2>& x); // updates x_u and x_d

    void solve(
        MatrixPartitionedTyings& matrix,
        const xt::xtensor<double, 1>& b,
        xt::xtensor<double, 1>& x); // updates x_u and x_d

    void solve_u(
        MatrixPartitionedTyings& matrix,
        const xt::xtensor<double, 1>& b_u,
        const xt::xtensor<double, 1>& b_d,
        const xt::xtensor<double, 1>& x_p,
        xt::xtensor<double, 1>& x_u);

    // Auto-allocation of the functions above
    xt::xtensor<double, 2> Solve(
        MatrixPartitionedTyings& matrix,
        const xt::xtensor<double, 2>& b,
        const xt::xtensor<double, 2>& x);

    xt::xtensor<double, 1> Solve(
        MatrixPartitionedTyings& matrix,
        const xt::xtensor<double, 1>& b,
        const xt::xtensor<double, 1>& x);

    xt::xtensor<double, 1> Solve_u(
        MatrixPartitionedTyings& matrix,
        const xt::xtensor<double, 1>& b_u,
        const xt::xtensor<double, 1>& b_d,
        const xt::xtensor<double, 1>& x_p);

private:
    Solver m_solver; // solver
    bool m_factor = false; // signal to force factorization
    void factorize(MatrixPartitionedTyings& matrix); // compute inverse (evaluated by "solve")
};

} // namespace GooseFEM

#include "MatrixPartitionedTyings.hpp"

#endif

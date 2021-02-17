/*

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

*/

#ifndef GOOSEFEM_ELEMENTQUAD4PLANAR_H
#define GOOSEFEM_ELEMENTQUAD4PLANAR_H

#include "config.h"

namespace GooseFEM {
namespace Element {
namespace Quad4 {

class QuadraturePlanar : public GooseFEM::Element::QuadratureBase<4, 2, 3> {
public:
    // Fixed dimensions:
    //    ndim = 2   -  number of dimensions
    //    nne  = 4   -  number of nodes per element
    //    tdim = 3   -  number of dimensions of tensors
    //
    // Naming convention:
    //    "elemmat"  -  matrices stored per element       -  [nelem, nne*ndim, nne*ndim]
    //    "elemvec"  -  nodal vectors stored per element  -  [nelem, nne, ndim]
    //    "qtensor"  -  integration point tensor          -  [nelem, nip, tdim, tdim]
    //    "qscalar"  -  integration point scalar          -  [nelem, nip]

    // Constructor: integration point coordinates and weights are optional (default: Gauss)
    QuadraturePlanar() = default;

    QuadraturePlanar(const xt::xtensor<double, 3>& x, double thick = 1.0);

    QuadraturePlanar(
        const xt::xtensor<double, 3>& x,
        const xt::xtensor<double, 2>& xi,
        const xt::xtensor<double, 1>& w,
        double thick = 1.0);

    // Update the nodal positions (shape of "x" should match the earlier definition)
    void update_x(const xt::xtensor<double, 3>& x);

    // Return shape function gradients
    xt::xtensor<double, 4> GradN() const;

    // Return integration volume
    xt::xtensor<double, 2> dV() const;

    // Dyadic product (and its transpose and symmetric part)
    // qtensor(i,j) += dNdx(m,i) * elemvec(m,j)
    void gradN_vector(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;
    void gradN_vector_T(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;
    void symGradN_vector(const xt::xtensor<double, 3>& elemvec, xt::xtensor<double, 4>& qtensor) const;

    // Integral of the scalar product
    // elemmat(m*ndim+i,n*ndim+i) += N(m) * qscalar * N(n) * dV
    void int_N_scalar_NT_dV(
        const xt::xtensor<double, 2>& qscalar, xt::xtensor<double, 3>& elemmat) const;

    // Integral of the dot product
    // elemvec(m,j) += dNdx(m,i) * qtensor(i,j) * dV
    void int_gradN_dot_tensor2_dV(
        const xt::xtensor<double, 4>& qtensor, xt::xtensor<double, 3>& elemvec) const;

    // Integral of the dot product
    // elemmat(m*2+j, n*2+k) += dNdx(m,i) * qtensor(i,j,k,l) * dNdx(n,l) * dV
    void int_gradN_dot_tensor4_dot_gradNT_dV(
        const xt::xtensor<double, 6>& qtensor, xt::xtensor<double, 3>& elemmat) const;

    // Auto-allocation of the functions above
    xt::xtensor<double, 4> GradN_vector(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 4> GradN_vector_T(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 4> SymGradN_vector(const xt::xtensor<double, 3>& elemvec) const;
    xt::xtensor<double, 3> Int_N_scalar_NT_dV(const xt::xtensor<double, 2>& qscalar) const;
    xt::xtensor<double, 3> Int_gradN_dot_tensor2_dV(const xt::xtensor<double, 4>& qtensor) const;
    xt::xtensor<double, 3> Int_gradN_dot_tensor4_dot_gradNT_dV(const xt::xtensor<double, 6>& qtensor) const;

private:
    // Compute "vol" and "dNdx" based on current "x"
    void compute_dN();

private:
    xt::xtensor<double, 3> m_x;    // nodal positions stored per element [nelem, nne, ndim]
    xt::xtensor<double, 1> m_w;    // weight of each integration point [nip]
    xt::xtensor<double, 2> m_xi;   // local coordinate of each integration point [nip, ndim]
    xt::xtensor<double, 2> m_N;    // shape functions [nip, nne]
    xt::xtensor<double, 3> m_dNxi; // shape function grad. wrt local  coor. [nip, nne, ndim]
    xt::xtensor<double, 4> m_dNx;  // shape function grad. wrt global coor. [nelem, nip, nne, ndim]
    xt::xtensor<double, 2> m_vol;  // integration point volume [nelem, nip]

    // Thickness
    double m_thick;
};

} // namespace Quad4
} // namespace Element
} // namespace GooseFEM

#include "ElementQuad4Planar.hpp"

#endif

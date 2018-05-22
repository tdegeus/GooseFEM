
import numpy as np

np.set_printoptions(linewidth=200)

# ==================================================================================================
# tensor products
# ==================================================================================================

def ddot43(A,B):

  nd = A.shape[0]

  C = np.zeros((nd,nd,nd))

  for i in range(nd):
    for j in range(nd):
      for k in range(nd):
        for l in range(nd):
          for m in range(nd):
            C[i,j,m] += A[i,j,k,l] * B[l,k,m]

  return C

# --------------------------------------------------------------------------------------------------

def ddot33(A,B):

  nd = A.shape[0]

  C = np.zeros((nd,nd))

  for i in range(nd):
    for j in range(nd):
      for k in range(nd):
        for l in range(nd):
          C[i,l] += A[i,j,k] * B[k,j,l]

  return C

# ==================================================================================================
# identity tensors
# ==================================================================================================

II   = np.zeros((3,3,3,3))
I4   = np.zeros((3,3,3,3))
I4rt = np.zeros((3,3,3,3))
Is   = np.zeros((3,3,3,3))

for i in range(3):
  for j in range(3):
    for k in range(3):
      for l in range(3):
        if ( i == j and k == l ):
          II[i,i,j,l] += 1.

for i in range(3):
  for j in range(3):
    for k in range(3):
      for l in range(3):
        if ( i == l and j == k ):
          I4[i,i,j,l] += 1.

for i in range(3):
  for j in range(3):
    for k in range(3):
      for l in range(3):
        if ( i == k and j == l ):
          I4rt[i,i,j,l] += 1.

I4s = ( I4  + I4rt  ) / 2.
I4d = ( I4s - II/3. )

# ==================================================================================================
# elasticity tensor
# ==================================================================================================

K = 1.
G = 1.

C4 = K * II + 2. * G * I4d

# ==================================================================================================
# mesh definition
# ==================================================================================================

nne  = 4
ndim = 2

xe = np.array([
  [0., 0.],
  [1., 0.],
  [1., 1.],
  [0., 1.],
])

Xi = np.array([
  [-1./np.sqrt(3.), -1./np.sqrt(3.)],
  [+1./np.sqrt(3.), -1./np.sqrt(3.)],
  [+1./np.sqrt(3.), +1./np.sqrt(3.)],
  [-1./np.sqrt(3.), +1./np.sqrt(3.)],
])

W = np.array([
  [1.],
  [1.],
  [1.],
  [1.],
])

Ke = np.zeros((nne*ndim, nne*ndim))

for w, xi in zip(W, Xi):

  N = np.array([
    [.25 * (1.-xi[0]) * (1.-xi[1])],
    [.25 * (1.+xi[0]) * (1.-xi[1])],
    [.25 * (1.+xi[0]) * (1.+xi[1])],
    [.25 * (1.-xi[0]) * (1.+xi[1])],
  ])

  dNdxi = np.array([
    [-.25 * (1.-xi[1]), -.25 * (1.-xi[0])],
    [+.25 * (1.-xi[1]), -.25 * (1.+xi[0])],
    [+.25 * (1.+xi[1]), +.25 * (1.+xi[0])],
    [-.25 * (1.+xi[1]), +.25 * (1.-xi[0])],
  ])

  J = np.zeros((ndim, ndim))

  for n in range(nne):
    for i in range(ndim):
      for j in range(ndim):
        J[i,j] += dNdxi[n,i] * xe[n,j]

  Jinv = np.linalg.inv(J)
  Jdet = np.linalg.det(J)

  dNdx = np.zeros((nne,ndim))

  for n in range(nne):
    for i in range(ndim):
      for j in range(ndim):
        dNdx[n,i] += Jinv[i,j] * dNdxi[n,j]

  xk = np.zeros((ndim))

  for n in range(nne):
    for i in range(ndim):
      xk[i] += N[n] * xe[n,i]

  rk = xk[0]

  for n in range(nne):

    Bn = np.zeros((3,3,3))

    Bn[0,0,0] = dNdx[n,0]        # B(n, r      r      r      )
    Bn[0,1,0] = -1./rk * N[n]    # B(n, r      \theta r      )
    Bn[0,1,1] = dNdx[n,0]        # B(n, r      \theta \theta )
    Bn[0,2,2] = dNdx[n,0]        # B(n, r      z      z      )

    Bn[1,0,0] = 0.               # B(n, \theta r      r      ) - axisymmetric
    Bn[1,1,0] = +1./rk * N[n]    # B(n, \theta \theta r      )
    Bn[1,1,1] = 0.               # B(n, \theta \theta \theta ) - axisymmetric
    Bn[1,2,2] = 0.               # B(n, \theta z      z      ) - axisymmetric

    Bn[2,0,0] = dNdx[n,1]        # B(n, z      r      r      )
    Bn[2,1,1] = dNdx[n,1]        # B(n, z      \theta \theta )
    Bn[2,2,2] = dNdx[n,1]        # B(n, z      z      z      )

    for m in range(nne):

      Bm = np.zeros((3,3,3))

      Bm[0,0,0] = dNdx[m,0]        # B(n, r      r      r      )
      Bm[0,1,0] = -1./rk * N[n]    # B(n, r      \theta r      )
      Bm[0,1,1] = dNdx[m,0]        # B(n, r      \theta \theta )
      Bm[0,2,2] = dNdx[m,0]        # B(n, r      z      z      )

      Bm[1,0,0] = 0.               # B(n, \theta r      r      ) - axisymmetric
      Bm[1,1,0] = +1./rk * N[n]    # B(n, \theta \theta r      )
      Bm[1,1,1] = 0.               # B(n, \theta \theta \theta ) - axisymmetric
      Bm[1,2,2] = 0.               # B(n, \theta z      z      ) - axisymmetric

      Bm[2,0,0] = dNdx[m,1]        # B(n, z      r      r      )
      Bm[2,1,1] = dNdx[m,1]        # B(n, z      \theta \theta )
      Bm[2,2,2] = dNdx[m,1]        # B(n, z      z      z      )

      K = ddot33(Bn,ddot43(C4,Bm))

      iin = np.array([ 2*n+0, 2*n+1 ])
      iim = np.array([ 2*m+0, 2*m+1 ])

      Ke[np.ix_(iin,iim)] += K[np.ix_([0,2], [0,2])]

  # apply the integration point volume
  Ke *= w * Jdet * rk

print(Ke)


  # Ke[]










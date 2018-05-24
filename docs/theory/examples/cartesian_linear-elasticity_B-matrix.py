
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

# --------------------------------------------------------------------------------------------------

def transpose3(A):

  nd = A.shape[0]

  C = np.zeros((nd,nd,nd))

  for i in range(nd):
    for j in range(nd):
      for k in range(nd):
        C[k,j,i] = A[i,j,k]

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
          II[i,j,k,l] = 1.

for i in range(3):
  for j in range(3):
    for k in range(3):
      for l in range(3):
        if ( i == l and j == k ):
          I4[i,j,k,l] = 1.

for i in range(3):
  for j in range(3):
    for k in range(3):
      for l in range(3):
        if ( i == k and j == l ):
          I4rt[i,j,k,l] = 1.

I4s = ( I4  + I4rt  ) / 2.
I4d = ( I4s - II/3. )

# ==================================================================================================
# elasticity tensor
# ==================================================================================================

K = 0.8333333333333333
G = 0.3846153846153846

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

  dNdxi = np.array([
    [-.25*(1.-xi[1]), -.25*(1.-xi[0])],
    [+.25*(1.-xi[1]), -.25*(1.+xi[0])],
    [+.25*(1.+xi[1]), +.25*(1.+xi[0])],
    [-.25*(1.+xi[1]), +.25*(1.-xi[0])],
  ])

  J = np.zeros((ndim, ndim))

  for m in range(nne):
    for i in range(ndim):
      for j in range(ndim):
        J[i,j] += dNdxi[m,i] * xe[m,j]

  Jinv = np.linalg.inv(J)
  Jdet = np.linalg.det(J)

  dNdx = np.zeros((nne,ndim))

  for m in range(nne):
    for i in range(ndim):
      for j in range(ndim):
        dNdx[m,i] += Jinv[i,j] * dNdxi[m,j]

  for m in range(nne):

    Bm = np.zeros((3,3,3))

    Bm[0,0,0] = dNdx[m,0]
    Bm[0,1,1] = dNdx[m,0]

    Bm[1,0,0] = dNdx[m,1]
    Bm[1,1,1] = dNdx[m,1]

    for n in range(nne):

      Bn = np.zeros((3,3,3))

      Bn[0,0,0] = dNdx[n,0]
      Bn[0,1,1] = dNdx[n,0]

      Bn[1,0,0] = dNdx[n,1]
      Bn[1,1,1] = dNdx[n,1]

      Kmn = ddot33(transpose3(Bm),ddot43(C4,Bn))

      iim = np.array([ m*ndim+0, m*ndim+1 ])
      iin = np.array([ n*ndim+0, n*ndim+1 ])

      Ke[np.ix_(iim,iin)] += Kmn[np.ix_([0,1], [0,1])] * w * Jdet

print(Ke)

# print(np.linalg.inv(Ke))


  # Ke[]










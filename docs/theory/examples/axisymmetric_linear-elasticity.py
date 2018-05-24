
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

# bulk and shear modulus
K = 0.8333333333333333
G = 0.3846153846153846

# elasticity tensor (3d)
C4 = K * II + 2. * G * I4d

# ==================================================================================================
# mesh definition
# ==================================================================================================

# number of elements
nz = 3
nr = 3

# mesh dimensions
nelem =  nz    *  nr     # number of elements
nnode = (nz+1) * (nr+1)  # number of nodes
nne   = 4                # number of nodes per element
nd    = 2                # number of dimensions
ndof  = nnode * nd       # total number of degrees of freedom

# coordinates and connectivity: zero-initialize
coor = np.zeros((nnode,nd ), dtype='float')
conn = np.zeros((nelem,nne), dtype='int'  )

# coordinates: set
# - grid of points
z,r = np.meshgrid(np.linspace(0,1,nz+1), np.linspace(0,1,nr+1))
# - store from grid of points
coor[:,0] = z.ravel()
coor[:,1] = r.ravel()

# connectivity: set
# - grid of node numbers
inode = np.arange(nnode).reshape(nr+1, nz+1)
# - store from grid of node numbers
conn[:,0] = inode[ :-1, :-1].ravel()
conn[:,1] = inode[ :-1,1:  ].ravel()
conn[:,2] = inode[1:  ,1:  ].ravel()
conn[:,3] = inode[1:  , :-1].ravel()

# DOF-numbers per node
dofs = np.arange(nnode*nd).reshape(nnode,nd)

# ==================================================================================================
# quadrature definition
# ==================================================================================================

# integration point coordinates (local element coordinates)
Xi = np.array([
  [-1./np.sqrt(3.), -1./np.sqrt(3.)],
  [+1./np.sqrt(3.), -1./np.sqrt(3.)],
  [+1./np.sqrt(3.), +1./np.sqrt(3.)],
  [-1./np.sqrt(3.), +1./np.sqrt(3.)],
])

# integration point weights
W = np.array([
  [1.],
  [1.],
  [1.],
  [1.],
])

# ==================================================================================================
# assemble stiffness matrix
# ==================================================================================================

# allocate matrix
K = np.zeros((ndof, ndof))

# loop over nodes
for e in conn:

  # - nodal coordinates
  xe = coor[e,:]

  # - allocate element stiffness matrix
  Ke = np.zeros((nne*nd, nne*nd))

  # - numerical quadrature
  for w, xi in zip(W, Xi):

    # -- shape functions
    N = np.array([
      [.25 * (1.-xi[0]) * (1.-xi[1])],
      [.25 * (1.+xi[0]) * (1.-xi[1])],
      [.25 * (1.+xi[0]) * (1.+xi[1])],
      [.25 * (1.-xi[0]) * (1.+xi[1])],
    ])

    # -- shape function gradients (w.r.t. the local element coordinates)
    dNdxi = np.array([
      [-.25*(1.-xi[1]), -.25*(1.-xi[0])],
      [+.25*(1.-xi[1]), -.25*(1.+xi[0])],
      [+.25*(1.+xi[1]), +.25*(1.+xi[0])],
      [-.25*(1.+xi[1]), +.25*(1.-xi[0])],
    ])

    # -- Jacobian
    J = np.zeros((nd, nd))

    for m in range(nne):
      for i in range(nd):
        for j in range(nd):
          J[i,j] += dNdxi[m,i] * xe[m,j]

    Jinv = np.linalg.inv(J)
    Jdet = np.linalg.det(J)

    # -- shape function gradients (w.r.t. the global coordinates)
    dNdx = np.zeros((nne,nd))

    for m in range(nne):
      for i in range(nd):
        for j in range(nd):
          dNdx[m,i] += Jinv[i,j] * dNdxi[m,j]

    # -- global coordinates of the integration point
    xk = np.zeros((nd))

    for n in range(nne):
      for i in range(nd):
        xk[i] += N[n] * xe[n,i]

    rk = xk[1]

    # -- assemble to element stiffness matrix
    for m in range(nne):

      Bm = np.zeros((3,3,3))

      Bm[0,0,0] = dNdx[m,1]        # B(m, r      r      r      )
      Bm[0,1,0] = -1./rk * N[m]    # B(m, r      \theta r      )
      Bm[0,1,1] = dNdx[m,1]        # B(m, r      \theta \theta )
      Bm[0,2,2] = dNdx[m,1]        # B(m, r      z      z      )

      Bm[1,0,0] = 0.               # B(m, \theta r      r      ) - axisymmetric
      Bm[1,1,0] = +1./rk * N[m]    # B(m, \theta \theta r      )
      Bm[1,1,1] = 0.               # B(m, \theta \theta \theta ) - axisymmetric
      Bm[1,2,2] = 0.               # B(m, \theta z      z      ) - axisymmetric

      Bm[2,0,0] = dNdx[m,0]        # B(m, z      r      r      )
      Bm[2,1,1] = dNdx[m,0]        # B(m, z      \theta \theta )
      Bm[2,2,2] = dNdx[m,0]        # B(m, z      z      z      )

      for n in range(nne):

        Bn = np.zeros((3,3,3))

        Bn[0,0,0] = dNdx[n,1]        # B(n, r      r      r      )
        Bn[0,1,0] = -1./rk * N[n]    # B(n, r      \theta r      )
        Bn[0,1,1] = dNdx[n,1]        # B(n, r      \theta \theta )
        Bn[0,2,2] = dNdx[n,1]        # B(n, r      z      z      )

        Bn[1,0,0] = 0.               # B(n, \theta r      r      ) - axisymmetric
        Bn[1,1,0] = +1./rk * N[n]    # B(n, \theta \theta r      )
        Bn[1,1,1] = 0.               # B(n, \theta \theta \theta ) - axisymmetric
        Bn[1,2,2] = 0.               # B(n, \theta z      z      ) - axisymmetric

        Bn[2,0,0] = dNdx[n,0]        # B(n, z      r      r      )
        Bn[2,1,1] = dNdx[n,0]        # B(n, z      \theta \theta )
        Bn[2,2,2] = dNdx[n,0]        # B(n, z      z      z      )

        Kmn = ddot33(transpose3(Bm),ddot43(C4,Bn))

        iim = np.array([ m*nd+0, m*nd+1 ])
        iin = np.array([ n*nd+0, n*nd+1 ])

        Ke[np.ix_(iim,iin)] += Kmn[np.ix_([2,0], [2,0])] * w * Jdet * rk * 2. * np.pi

  # - assemble to global stiffness matrix
  iie = dofs[e,:].ravel()

  K[np.ix_(iie,iie)] += Ke

# ==================================================================================================
# partition and solve
# ==================================================================================================

# zero-initialize displacements and forces
f = np.zeros((ndof))
u = np.zeros((ndof))

# fixed displacements
# - zero-initialize
iip  = np.empty((0),dtype='int'  )
up   = np.empty((0),dtype='float')
# - 'r = 0' : 'u_r = 0'
idof = dofs[inode[0,:], 1]
iip  = np.hstack(( iip, idof                      ))
up   = np.hstack(( up , 0.0 * np.ones(idof.shape) ))
# - 'z = 0' : 'u_z = 0'
idof = dofs[inode[:,0], 0]
iip  = np.hstack(( iip, idof                      ))
up   = np.hstack(( up , 0.0 * np.ones(idof.shape) ))
# - 'z = 1' : 'u_z = 0.1'
idof = dofs[inode[:,-1], 0]
iip  = np.hstack(( iip, idof                      ))
up   = np.hstack(( up , 0.1 * np.ones(idof.shape) ))

# free displacements
iiu  = np.setdiff1d(dofs.ravel(), iip)

# partition
# - stiffness matrix
Kuu  = K[np.ix_(iiu, iiu)]
Kup  = K[np.ix_(iiu, iip)]
Kpu  = K[np.ix_(iip, iiu)]
Kpp  = K[np.ix_(iip, iip)]
# - external force
fu   = f[iiu]

# solve
uu   = np.linalg.solve(Kuu, fu - Kup.dot(up))
fp   = Kpu.dot(uu) + Kpp.dot(up)

# assemble
u[iiu] = uu
u[iip] = up
f[iip] = fp

# convert to nodal displacement and forces
disp = u[dofs]
fext = f[dofs]

# ==================================================================================================
# plot
# ==================================================================================================

import matplotlib.pyplot as plt

fig, axes =  plt.subplots(ncols=2, figsize=(16,8))

axes[0].scatter(coor[:,0], coor[:,1])
axes[0].quiver (coor[:,0], coor[:,1], disp[:,0], disp[:,1])

axes[1].scatter((coor+disp)[:,0], (coor+disp)[:,1])
axes[1].quiver ((coor+disp)[:,0], (coor+disp)[:,1], fext[:,0], fext[:,1], scale=.1)

for axis in axes:
  axis.set_xlim([-0.4, 1.5])
  axis.set_ylim([-0.4, 1.5])

plt.show()

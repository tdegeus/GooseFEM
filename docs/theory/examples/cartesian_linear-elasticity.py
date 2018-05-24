
import numpy as np

np.set_printoptions(linewidth=200)

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
nx = 10
ny = 10

# mesh dimensions
nelem =  nx    *  ny     # number of elements
nnode = (nx+1) * (ny+1)  # number of nodes
nne   = 4                # number of nodes per element
nd    = 2                # number of dimensions
ndof  = nnode * nd       # total number of degrees of freedom

# coordinates and connectivity: zero-initialize
coor = np.zeros((nnode,nd ), dtype='float')
conn = np.zeros((nelem,nne), dtype='int'  )

# coordinates: set
# - grid of points
x,y = np.meshgrid(np.linspace(0,1,nx+1), np.linspace(0,1,ny+1))
# - store from grid of points
coor[:,0] = x.ravel()
coor[:,1] = y.ravel()

# connectivity: set
# - grid of node numbers
inode = np.arange(nnode).reshape(ny+1, nx+1)
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

    # -- assemble to element stiffness matrix
    for m in range(nne):
      for n in range(nne):
        for i in range(nd):
          for j in range(nd):
            for k in range(nd):
              for l in range(nd):
                Ke[m*nd+j, n*nd+k] += dNdx[m,i] * C4[i,j,k,l] * dNdx[n,l] * w * Jdet

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
# - 'y = 0' : 'u_y = 0'
idof = dofs[inode[0,:], 1]
iip  = np.hstack(( iip, idof                      ))
up   = np.hstack(( up , 0.0 * np.ones(idof.shape) ))
# - 'x = 0' : 'u_x = 0'
idof = dofs[inode[:,0], 0]
iip  = np.hstack(( iip, idof                      ))
up   = np.hstack(( up , 0.0 * np.ones(idof.shape) ))
# - 'x = 1' : 'u_x = 0.1'
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

print(fext[inode[:,-1],0].sum())

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

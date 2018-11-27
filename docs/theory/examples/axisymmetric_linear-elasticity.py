
import numpy as np

np.set_printoptions(linewidth=200)

# ==================================================================================================
# mesh definition
# ==================================================================================================

# number of elements
nz = 20
nr = 16

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

# number of integration points
nip = 4

# ==================================================================================================
# integration point tensors and operations
# ==================================================================================================

# tensors products
ddot22     = lambda A2,B2: np.einsum('...ij  ,...ji->...    ',A2,B2)
ddot42     = lambda A4,B2: np.einsum('...ijkl,...lk->...ij  ',A4,B2)
dyad22     = lambda A2,B2: np.einsum('...ij  ,...kl->...ijkl',A2,B2)
ddot43     = lambda A4,B3: np.einsum('...ijkl,...lkm->...ijm',A4,B3)
ddot33     = lambda A3,B3: np.einsum('...ijk ,...kjl->...il ',A3,B3)
dot31      = lambda A3,B1: np.einsum('...ijk ,...k  ->...ij ',A3,B1)
transpose3 = lambda A3   : np.einsum('...ijk        ->...kji',A3   )

# identity tensors (per integration point)
i    = np.eye(3)
I    = np.einsum('ij  ,...'         ,                  i   ,np.ones([nelem,nip]))
I4   = np.einsum('ijkl,...->...ijkl',np.einsum('il,jk',i,i),np.ones([nelem,nip]))
I4rt = np.einsum('ijkl,...->...ijkl',np.einsum('ik,jl',i,i),np.ones([nelem,nip]))
I4s  = (I4+I4rt)/2.
II   = dyad22(I,I)
I4d  = I4s-II/3.

# bulk and shear modulus
K = 0.8333333333333333
G = 0.3846153846153846

# elasticity tensor (per integration point)
C4 = K * II + 2. * G * I4d

# ==================================================================================================
# B-matrices at each integration point
# ==================================================================================================

# allocate matrix
B   = np.zeros((nelem,nip,nne,3,3,3))
vol = np.zeros((nelem,nip))

# loop over nodes
for e, nodes in enumerate(conn):

  # nodal coordinates
  xe = coor[nodes,:]

  # loop over integration points
  for q, (w, xi) in enumerate(zip(W, Xi)):

    # shape functions
    N = np.array([
      .25 * (1.-xi[0]) * (1.-xi[1]),
      .25 * (1.+xi[0]) * (1.-xi[1]),
      .25 * (1.+xi[0]) * (1.+xi[1]),
      .25 * (1.-xi[0]) * (1.+xi[1]),
    ])

    # shape function gradients (w.r.t. the local element coordinates)
    dNdxi = np.array([
      [-.25*(1.-xi[1]), -.25*(1.-xi[0])],
      [+.25*(1.-xi[1]), -.25*(1.+xi[0])],
      [+.25*(1.+xi[1]), +.25*(1.+xi[0])],
      [-.25*(1.+xi[1]), +.25*(1.-xi[0])],
    ])

    # Jacobian
    J    = np.einsum('mi,mj->ij', dNdxi, xe)
    Jinv = np.linalg.inv(J)
    Jdet = np.linalg.det(J)

    # shape function gradients (w.r.t. the global coordinates)
    dNdx = np.einsum('ij,mj->mi', Jinv, dNdxi)

    # global coordinates of the integration point
    xq = np.einsum('m,mi->i', N, xe)
    rq = xq[1]

    # compute B-matrix and integration-point volume and store for later use
    for m in range(nne):

      Bm = np.zeros((3,3,3))

      Bm[0,1,1] = -1./rq * N[m]    # B(m, r      \theta \theta )
      Bm[1,1,0] = +1./rq * N[m]    # B(m, \theta \theta r      )

      Bm[0,0,0] = dNdx[m,1]        # B(m, r      r      r      )
      Bm[0,1,1] = dNdx[m,1]        # B(m, r      \theta \theta )
      Bm[0,2,2] = dNdx[m,1]        # B(m, r      z      z      )

      Bm[2,0,0] = dNdx[m,0]        # B(m, z      r      r      )
      Bm[2,1,1] = dNdx[m,0]        # B(m, z      \theta \theta )
      Bm[2,2,2] = dNdx[m,0]        # B(m, z      z      z      )

      B[e,q,m,:,:,:] = Bm

      vol[e,q] = w * Jdet * rq * 2. * np.pi

# ==================================================================================================
# assemble stiffness matrix
# ==================================================================================================

# allocate matrix
K = np.zeros((ndof, ndof))

# loop over nodes
for e, nodes in enumerate(conn):

  # allocate element stiffness matrix
  Ke = np.zeros((nne*nd, nne*nd))

  # loop over integration points
  for q, (w, xi) in enumerate(zip(W, Xi)):

    # alias integration point values
    C4q = C4[e,q,:,:,:,:]
    dV  = vol[e,q]

    # assemble to element stiffness matrix
    for m in range(nne):

      Bm = B[e,q,m,:,:,:]

      for n in range(nne):

        Bn = B[e,q,n,:,:,:]

        Kmn = ddot33(transpose3(Bm),ddot43(C4q,Bn))

        iim = np.array([ m*nd+0, m*nd+1 ])
        iin = np.array([ n*nd+0, n*nd+1 ])

        Ke[np.ix_(iim,iin)] += Kmn[np.ix_([2,0], [2,0])] * dV

  # assemble to global stiffness matrix
  iie = dofs[nodes,:].ravel()

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
# compute strain and stress
# ==================================================================================================

# zero-initialize
# - strain tensor per integration point
eps = np.zeros((nelem,nip,3,3))
# - nodal displacement in 3-d
um = np.zeros((3))

# loop over nodes
for e, nodes in enumerate(conn):

  # nodal displacements
  ue = disp[nodes,:]

  # loop over integration points
  for q, (w, xi) in enumerate(zip(W, Xi)):

    # zero-initialize displacements gradient
    gradu = np.zeros((3,3))

    # compute displacement gradient
    for m in range(nne):

      # alias
      # - B-matrix
      Bm = B[e,q,m,:,:,:]
      # - nodal displacement in 3-d (theta-component always zero)
      um[np.ix_([2,0])] = ue[m,:]

      # update
      gradu += dot31(Bm, um)

    # compute strain tensor, and store per integration point
    eps[e,q] = .5 * ( gradu + gradu.T )

# compute stress tensor (per integration point)
sig = ddot42(C4,eps)

# ==================================================================================================
# plot
# ==================================================================================================

import matplotlib.pyplot as plt

fig, axes =  plt.subplots(ncols=2, figsize=(16,8))

axes[0].scatter(coor[:,0], coor[:,1])
axes[0].quiver (coor[:,0], coor[:,1], disp[:,0], disp[:,1])

axes[1].scatter((coor+disp)[:,0], (coor+disp)[:,1])
axes[1].quiver ((coor+disp)[:,0], (coor+disp)[:,1], fext[:,0], fext[:,1], scale=.4)

for axis in axes:
  axis.set_xlim([-0.4, 1.5])
  axis.set_ylim([-0.4, 1.5])

plt.show()

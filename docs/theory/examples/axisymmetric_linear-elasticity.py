# ==================================================================================================
#
# (c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM
#
# ==================================================================================================

import numpy as np

np.set_printoptions(linewidth=200)

# ==================================================================================================
# tensor products
# ==================================================================================================

ddot22     = lambda A2,B2: np.einsum('...ij  ,...ji->...    ',A2,B2)
ddot42     = lambda A4,B2: np.einsum('...ijkl,...lk->...ij  ',A4,B2)
dyad22     = lambda A2,B2: np.einsum('...ij  ,...kl->...ijkl',A2,B2)
ddot43     = lambda A4,B3: np.einsum('...ijkl,...lkm->...ijm',A4,B3)
ddot33     = lambda A3,B3: np.einsum('...ijk ,...kjl->...il ',A3,B3)
dot31      = lambda A3,B1: np.einsum('...ijk ,...k  ->...ij ',A3,B1)
transpose3 = lambda A3   : np.einsum('...ijk        ->...kji',A3   )

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
ndim  = 2                # number of dimensions
ndof  = nnode * ndim     # total number of degrees of freedom

# zero-initialise coordinates, displacements, and connectivity
coor = np.zeros((nnode,ndim), dtype='float')
disp = np.zeros((nnode,ndim), dtype='float')
conn = np.zeros((nelem,nne ), dtype='int'  )

# coordinates
# - grid of points
z,r = np.meshgrid(np.linspace(0,1,nz+1), np.linspace(0,1,nr+1))
# - store from grid of points
coor[:,0] = z.ravel()
coor[:,1] = r.ravel()

# connectivity
# - grid of node numbers
inode = np.arange(nnode).reshape(nr+1, nz+1)
# - store from grid of node numbers
conn[:,0] = inode[ :-1, :-1].ravel()
conn[:,1] = inode[ :-1,1:  ].ravel()
conn[:,2] = inode[1:  ,1:  ].ravel()
conn[:,3] = inode[1:  , :-1].ravel()
# - node sets
nodesLeft   = inode[ :, 0]
nodesRight  = inode[ :,-1]
nodesBottom = inode[ 0, :]
nodesTop    = inode[-1, :]

# DOF-numbers per node
dofs = np.arange(nnode*ndim).reshape(nnode,ndim)

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
# B-matrix at each integration point
# ==================================================================================================

# allocate
B   = np.empty((nelem,nip,nne,3,3,3))
vol = np.empty((nelem,nip))

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
    Je    = np.einsum('mi,mj->ij', dNdxi, xe)
    invJe = np.linalg.inv(Je)
    detJe = np.linalg.det(Je)

    # shape function gradients (w.r.t. the global coordinates)
    dNdxe = np.einsum('ij,mj->mi', invJe, dNdxi)

    # global coordinates of the integration point, extract the radius
    xq = np.einsum('m,mi->i', N, xe)
    rq = xq[1]

    # compute B-matrix and integration-point volume and store for later use
    for m in range(nne):

      Be = np.zeros((3,3,3))

      Be[0,0,0] = dNdxe[m,1]    # B(m, r, r, r)
      Be[0,2,2] = dNdxe[m,1]    # B(m, r, z, z)
      Be[1,1,0] = +1./rq * N[m] # B(m, t, t, r)
      Be[2,0,0] = dNdxe[m,0]    # B(m, z, r, r)
      Be[2,2,2] = dNdxe[m,0]    # B(m, z, z, z)

      B[e,q,m,:,:,:] = Be

      vol[e,q] = w * detJe * rq * 2. * np.pi

# ==================================================================================================
# stiffness tensor at each integration point (provides constitutive response and 'tangent')
# ==================================================================================================

# identity tensors (per integration point)
i    = np.eye(3)
I    = np.einsum('ij  ,...'         ,                  i   ,np.ones([nelem,nip]))
I4   = np.einsum('ijkl,...->...ijkl',np.einsum('il,jk',i,i),np.ones([nelem,nip]))
I4rt = np.einsum('ijkl,...->...ijkl',np.einsum('ik,jl',i,i),np.ones([nelem,nip]))
I4s  = (I4+I4rt)/2.
II   = dyad22(I,I)
I4d  = I4s-II/3.

# bulk and shear modulus (homogeneous)
K = 0.8333333333333333 * np.ones((nelem, nip))
G = 0.3846153846153846 * np.ones((nelem, nip))

# elasticity tensor (per integration point)
C4 = np.einsum('eq,eqijkl->eqijkl', K, II) + 2. * np.einsum('eq,eqijkl->eqijkl', G, I4d)

# ==================================================================================================
# strain from nodal displacements, stress from constitutive response
# ==================================================================================================

# allocate strain tensor per integration point
Eps = np.empty((nelem,nip,3,3))

# loop over nodes
for e, nodes in enumerate(conn):

  # nodal displacements in 3-d
  #   u_theta = 0
  #   (z,r) -> (r,theta,z)
  ue      = np.zeros((nne,3))
  ue[:,0] = disp[nodes,1]
  ue[:,2] = disp[nodes,0]

  # loop over integration points
  for q, (w, xi) in enumerate(zip(W, Xi)):

    # alias integration point values
    Be = B[e,q,:,:,:,:]

    # displacement gradient
    gradu = np.einsum('mijk,mk->ij', Be, ue)

    # compute strain tensor, and store per integration point
    Eps[e,q] = .5 * ( gradu + gradu.T )

# constitutive response: strain tensor -> stress tensor (per integration point)
Sig = ddot42(C4, Eps)

# ==================================================================================================
# internal force from stress
# ==================================================================================================

# allocate internal force
fint = np.zeros((ndof))

# loop over nodes
for e, nodes in enumerate(conn):

  # allocate element internal force
  fe = np.zeros((nne*ndim))

  # loop over integration points
  for q, (w, xi) in enumerate(zip(W, Xi)):

    # alias integration point values
    sig = Sig[e,q,:,:]
    Be  = B  [e,q,:,:,:,:]
    dV  = vol[e,q]

    # assemble to element internal force
    #   Bm = B[e,q,m,:,:,:]
    #   fm = ddot32(transpose3(Bm), sig) * dV
    #   fe[[m*ndim+0, m*ndim+1]] = [fm[m,2], fm[m,0]]
    fq  = np.einsum('mijk,ij->mk', Be, sig) * dV
    fe += fq[:,[2,0]].reshape(nne*ndim)

  # assemble to global stiffness matrix
  iie = dofs[nodes,:].ravel()
  fint[iie] += fe

# ==================================================================================================
# stiffness matrix from 'tangent'
# ==================================================================================================

# allocate stiffness matrix
K = np.zeros((ndof, ndof))

# loop over nodes
for e, nodes in enumerate(conn):

  # allocate element stiffness matrix
  Ke = np.zeros((nne*ndim, nne*ndim))

  # loop over integration points
  for q, (w, xi) in enumerate(zip(W, Xi)):

    # alias integration point values
    c4  = C4 [e,q,:,:,:,:]
    Be  = B  [e,q,:,:,:,:]
    dV  = vol[e,q]

    # assemble to element stiffness matrix
    #   Bm  = B[e,q,m,:,:,:]
    #   Bn  = B[e,q,n,:,:,:]
    #   Kmn = ddot33(transpose3(Bm),ddot43(C4,Bn)) * dV
    #   Ke[[m*ndim+0,m*ndim+1], [n*ndim+0,n*ndim+1]] += Kmn[[2,0], [2,0]]
    Kq  = np.einsum('mabc,abde,nedf->mcnf', Be, c4, Be) * dV
    Ke += Kq[np.ix_(np.arange(nne),[2,0],np.arange(nne),[2,0])].reshape(nne*ndim, nne*ndim)

  # assemble to global stiffness matrix
  iie = dofs[nodes,:].ravel()
  K[np.ix_(iie,iie)] += Ke

# ==================================================================================================
# partition and solve
# ==================================================================================================

# prescribed external force: zero on all free DOFs
# (other DOFs ignored in solution, the reaction forces on the DOFs are computed below)
fext = np.zeros((ndof))

# DOF-partitioning: ['u'nknown, 'p'rescribed]
# - prescribed
iip = np.hstack((
  dofs[nodesBottom,1],
  dofs[nodesLeft  ,0],
  dofs[nodesRight ,0],
))
# - unknown
iiu = np.setdiff1d(dofs.ravel(), iip)

# fixed displacements
up = np.hstack((
  0.0 * np.ones((len(nodesBottom))),
  0.0 * np.ones((len(nodesLeft  ))),
  0.1 * np.ones((len(nodesRight ))),
))

# residual force
r = fext - fint

# partition
# - stiffness matrix
Kuu = K[np.ix_(iiu, iiu)]
Kup = K[np.ix_(iiu, iip)]
Kpu = K[np.ix_(iip, iiu)]
Kpp = K[np.ix_(iip, iip)]
# - residual force
ru = r[iiu]

# solve for unknown displacement DOFs
uu = np.linalg.solve(Kuu, ru - Kup.dot(up))

# assemble
u      = np.empty((ndof))
u[iiu] = uu
u[iip] = up

# convert to nodal displacements
disp = u[dofs]

# ==================================================================================================
# strain from nodal displacements, stress from constitutive response
# ==================================================================================================

# allocate strain tensor per integration point
Eps = np.empty((nelem,nip,3,3))

# loop over nodes
for e, nodes in enumerate(conn):

  # nodal displacements in 3-d
  #   u_theta = 0
  #   (z,r) -> (r,theta,z)
  ue      = np.zeros((nne,3))
  ue[:,0] = disp[nodes,1]
  ue[:,2] = disp[nodes,0]

  # loop over integration points
  for q, (w, xi) in enumerate(zip(W, Xi)):

    # alias integration point values
    Be = B[e,q,:,:,:,:]

    # displacement gradient
    gradu = np.einsum('mijk,mk->ij', Be, ue)

    # compute strain tensor, and store per integration point
    Eps[e,q] = .5 * ( gradu + gradu.T )

# constitutive response: strain tensor -> stress tensor (per integration point)
Sig = ddot42(C4, Eps)

# ==================================================================================================
# internal force from stress, reaction (external) force from internal force on fixed nodes
# ==================================================================================================

# allocate internal force
fint = np.zeros((ndof))

# loop over nodes
for e, nodes in enumerate(conn):

  # allocate element internal force
  fe = np.zeros((nne*ndim))

  # loop over integration points
  for q, (w, xi) in enumerate(zip(W, Xi)):

    # alias integration point values
    sig = Sig[e,q,:,:]
    Be  = B  [e,q,:,:,:,:]
    dV  = vol[e,q]

    # assemble to element internal force
    #   Bm = B[e,q,m,:,:,:]
    #   fm = ddot32(transpose3(Bm), sig) * dV
    #   fe[[m*ndim+0, m*ndim+1]] = [fm[m,2], fm[m,0]]
    fq  = np.einsum('mijk,ij->mk', Be, sig) * dV
    fe += fq[:,[2,0]].reshape(nne*ndim)

  # assemble to global stiffness matrix
  iie = dofs[nodes,:].ravel()
  fint[iie] += fe

# reaction force
fext[iip] = fint[iip]

# ==================================================================================================
# plot
# ==================================================================================================

import matplotlib.pyplot as plt

fig, axes = plt.subplots(ncols=2, figsize=(16,8))

# reconstruct external force as nodal vectors
fext = fext[dofs]

# plot original nodal positions + displacement field as arrows
axes[0].scatter(coor[:,0], coor[:,1])
axes[0].quiver (coor[:,0], coor[:,1], disp[:,0], disp[:,1])

# plot new nodal positions + external (reaction) force field as arrows
axes[1].scatter((coor+disp)[:,0], (coor+disp)[:,1])
axes[1].quiver ((coor+disp)[:,0], (coor+disp)[:,1], fext[:,0], fext[:,1], scale=.4)

# fix axes limits
for axis in axes:
  axis.set_xlim([-0.4, 1.5])
  axis.set_ylim([-0.4, 1.5])

plt.show()

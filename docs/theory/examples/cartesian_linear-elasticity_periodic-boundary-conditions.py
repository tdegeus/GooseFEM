# ==================================================================================================
#
# (c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM
#
# ==================================================================================================
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(linewidth=200)

# ==================================================================================================
# tensor products
# ==================================================================================================


def ddot22(A2, B2):
    return np.einsum("...ij,...ji->...", A2, B2)


def ddot42(A4, B2):
    return np.einsum("...ijkl,...lk->...ij", A4, B2)


def dyad22(A2, B2):
    return np.einsum("...ij,...kl->...ijkl", A2, B2)


# ==================================================================================================
# mesh definition
# ==================================================================================================

# number of elements
nx = 20
ny = 20

# mesh dimensions
nelem = nx * ny  # number of elements
nnode = (nx + 1) * (ny + 1)  # number of nodes
nne = 4  # number of nodes per element
ndim = 2  # number of dimensions
ndof = nnode * ndim  # total number of degrees of freedom

# out-of-plane thickness
thick = 1.0

# zero-initialise coordinates, displacements, and connectivity
coor = np.zeros((nnode, ndim), dtype="float")
disp = np.zeros((nnode, ndim), dtype="float")
conn = np.zeros((nelem, nne), dtype="int")

# coordinates
# - grid of points
x, y = np.meshgrid(np.linspace(0, 1, nx + 1), np.linspace(0, 1, ny + 1))
# - store from grid of points
coor[:, 0] = x.ravel()
coor[:, 1] = y.ravel()

# connectivity
# - grid of node numbers
inode = np.arange(nnode).reshape(ny + 1, nx + 1)
# - store from grid of node numbers
conn[:, 0] = inode[:-1, :-1].ravel()
conn[:, 1] = inode[:-1, 1:].ravel()
conn[:, 2] = inode[1:, 1:].ravel()
conn[:, 3] = inode[1:, :-1].ravel()
# - node sets
nodesLeft = inode[:, 0]
nodesRight = inode[:, -1]
nodesBottom = inode[0, :]
nodesTop = inode[-1, :]

# DOF-numbers per node
dofs = np.arange(nnode * ndim).reshape(nnode, ndim)

# ==================================================================================================
# periodicity
# ==================================================================================================

# add virtual nodes
# - node set
nodesVirtual = nnode + np.arange(2)
# - update size
nnode += 2
ndof += 4
# - add coordinates (position is completely arbitrary)
coor = np.vstack((coor, np.array([[0.0, 1.1], [0.1, 1.1]])))
# - add DOF numbers
dofs = np.vstack((dofs, np.max(dofs) + 1 + np.arange(4).reshape(2, 2)))

# tyings (dependent, independent)
tyings = np.vstack(
    (
        np.hstack((nodesRight[1:-1].reshape(-1, 1), nodesLeft[1:-1].reshape(-1, 1))),
        np.hstack((nodesTop[1:-1].reshape(-1, 1), nodesBottom[1:-1].reshape(-1, 1))),
        np.hstack((nodesRight[0].reshape(-1, 1), nodesBottom[0].reshape(-1, 1))),
        np.hstack((nodesRight[-1].reshape(-1, 1), nodesBottom[0].reshape(-1, 1))),
        np.hstack((nodesLeft[-1].reshape(-1, 1), nodesBottom[0].reshape(-1, 1))),
    )
)

# DOF-sets
# - dependent DOFs
iid = dofs[tyings[:, 0], :].ravel()
# - prescribed DOfs
iip = np.hstack(
    (
        dofs[nodesBottom[0], :].ravel(),
        dofs[nodesVirtual, :].ravel(),
    )
)
# - independent DOFs
iii = np.setdiff1d(dofs.ravel(), iid)
# - unknown DOFs
iiu = np.setdiff1d(iii, iip)

# renumber
# - list
renum = np.zeros(dofs.size, dtype="int")
renum[iiu] = np.arange(iiu.size)
renum[iip] = np.arange(iip.size) + iiu.size
renum[iid] = np.arange(iid.size) + iiu.size + iip.size
# - apply
iiu = renum[iiu]
iip = renum[iip]
dofs = renum[dofs]
# - define auxiliary DOF sets
iii = np.hstack((iiu, iip))
iid = np.arange(iid.size) + iii.size

# nodal tyings
# - zero-initialize
Cdi = np.zeros((iid.size, iii.size))
# - tie
iie = 2 * np.arange(tyings.shape[0])
Cdi[iie, dofs[tyings[:, 1], 0]] = 1.0
Cdi[iie + 1, dofs[tyings[:, 1], 1]] = 1.0
Cdi[iie, dofs[nodesVirtual[0], 0]] = coor[tyings[:, 0], 0] - coor[tyings[:, 1], 0]  # (F-I)_xx
Cdi[iie, dofs[nodesVirtual[0], 1]] = coor[tyings[:, 0], 1] - coor[tyings[:, 1], 1]  # (F-I)_xy
Cdi[iie + 1, dofs[nodesVirtual[1], 0]] = coor[tyings[:, 0], 0] - coor[tyings[:, 1], 0]  # (F-I)_yx
Cdi[iie + 1, dofs[nodesVirtual[1], 1]] = coor[tyings[:, 0], 1] - coor[tyings[:, 1], 1]  # (F-I)_yy

# ==================================================================================================
# quadrature definition
# ==================================================================================================

# integration point coordinates (local element coordinates)
Xi = np.array(
    [
        [-1.0 / np.sqrt(3.0), -1.0 / np.sqrt(3.0)],
        [+1.0 / np.sqrt(3.0), -1.0 / np.sqrt(3.0)],
        [+1.0 / np.sqrt(3.0), +1.0 / np.sqrt(3.0)],
        [-1.0 / np.sqrt(3.0), +1.0 / np.sqrt(3.0)],
    ]
)

# integration point weights
W = np.array(
    [
        [1.0],
        [1.0],
        [1.0],
        [1.0],
    ]
)

# number of integration points
nip = 4

# ==================================================================================================
# gradient of the shape functions at each integration point
# ==================================================================================================

# allocate
dNx = np.empty((nelem, nip, nne, ndim))
vol = np.empty((nelem, nip))

# loop over nodes
for e, nodes in enumerate(conn):
    # nodal coordinates
    xe = coor[nodes, :]

    # loop over integration points
    for q, (w, xi) in enumerate(zip(W, Xi)):
        # shape function gradients (w.r.t. the local element coordinates)
        dNdxi = np.array(
            [
                [-0.25 * (1.0 - xi[1]), -0.25 * (1.0 - xi[0])],
                [+0.25 * (1.0 - xi[1]), -0.25 * (1.0 + xi[0])],
                [+0.25 * (1.0 + xi[1]), +0.25 * (1.0 + xi[0])],
                [-0.25 * (1.0 + xi[1]), +0.25 * (1.0 - xi[0])],
            ]
        )

        # Jacobian
        Je = np.einsum("mi,mj->ij", dNdxi, xe)
        invJe = np.linalg.inv(Je)
        detJe = np.linalg.det(Je)

        # shape function gradients (w.r.t. the global coordinates)
        dNdxe = np.einsum("ij,mj->mi", invJe, dNdxi)

        # store for later use
        dNx[e, q, :, :] = dNdxe
        vol[e, q] = w * detJe * thick

# ==================================================================================================
# stiffness tensor at each integration point (provides constitutive response and 'tangent')
# ==================================================================================================

# identity tensors (per integration point)
i = np.eye(3)
I2 = np.einsum("ij  ,...", i, np.ones([nelem, nip]))
I4 = np.einsum("ijkl,...->...ijkl", np.einsum("il,jk", i, i), np.ones([nelem, nip]))
I4rt = np.einsum("ijkl,...->...ijkl", np.einsum("ik,jl", i, i), np.ones([nelem, nip]))
I4s = (I4 + I4rt) / 2.0
II = dyad22(I2, I2)
I4d = I4s - II / 3.0

# bulk and shear modulus (homogeneous)
K = 0.8333333333333333 * np.ones((nelem, nip))
G = 0.3846153846153846 * np.ones((nelem, nip))

# elasticity tensor (per integration point)
C4 = np.einsum("eq,eqijkl->eqijkl", K, II) + 2.0 * np.einsum("eq,eqijkl->eqijkl", G, I4d)

# ==================================================================================================
# strain from nodal displacements, stress from constitutive response
# ==================================================================================================

# zero-initialise strain tensor per integration point
# (plain strain -> all z-components are not written and should remain zero at all times)
Eps = np.zeros((nelem, nip, 3, 3))

# loop over nodes
for e, nodes in enumerate(conn):
    # nodal displacements
    ue = disp[nodes, :]

    # loop over integration points
    for q, (w, xi) in enumerate(zip(W, Xi)):
        # alias integration point values
        dNdxe = dNx[e, q, :, :]

        # displacement gradient
        gradu = np.einsum("mi,mj->ij", dNdxe, ue)

        # compute strain tensor, and store per integration point
        # (use plain strain to convert 2-d to 3-d tensor)
        Eps[e, q, :2, :2] = 0.5 * (gradu + gradu.T)

# constitutive response: strain tensor -> stress tensor (per integration point)
Sig = ddot42(C4, Eps)

# ==================================================================================================
# internal force from stress
# ==================================================================================================

# allocate internal force
fint = np.zeros(ndof)

# loop over nodes
for e, nodes in enumerate(conn):
    # allocate element internal force
    fe = np.zeros(nne * ndim)

    # loop over integration points
    for q, (w, xi) in enumerate(zip(W, Xi)):
        # alias integration point values
        # (plane strain: stress in z-direction irrelevant)
        sig = Sig[e, q, :2, :2]
        dNdxe = dNx[e, q, :, :]
        dV = vol[e, q]

        # assemble to element internal force
        #   fe[m*nd+j] += dNdx[m,i] * sig[i,j] * dV
        fe += (np.einsum("mi,ij->mj", dNdxe, sig) * dV).reshape(nne * ndim)

    # assemble to global stiffness matrix
    iie = dofs[nodes, :].ravel()
    fint[iie] += fe

# ==================================================================================================
# stiffness matrix from 'tangent'
# ==================================================================================================

# allocate stiffness matrix
K = np.zeros((ndof, ndof))

# loop over nodes
for e, nodes in enumerate(conn):
    # allocate element stiffness matrix
    Ke = np.zeros((nne * ndim, nne * ndim))

    # loop over integration points
    for q, (w, xi) in enumerate(zip(W, Xi)):
        # alias integration point values
        c4 = C4[e, q, :2, :2, :2, :2]
        dNdxe = dNx[e, q, :, :]
        dV = vol[e, q]

        # assemble to element stiffness matrix
        #   Ke[m*nd+j, n*nd+k] += dNdx[m,i] * C4[i,j,k,l] * dNdx[n,l] * dV
        Ke += (np.einsum("mi,ijkl,nl->mjnk", dNdxe, c4, dNdxe) * dV).reshape(nne * ndim, nne * ndim)

    # assemble to global stiffness matrix
    iie = dofs[nodes, :].ravel()
    K[np.ix_(iie, iie)] += Ke

# ==================================================================================================
# partition and solve
# ==================================================================================================

# prescribed external force: zero on all free DOFs
# (other DOFs ignored in solution, the reaction forces on the DOFs are computed below)
fext = np.zeros(ndof)

# fixed displacements
up = np.array(
    [
        0.0,  # suppress rigid-body motion in x
        0.0,  # suppress rigid-body motion in x
        0.0,  # (F-I)_xx
        0.1,  # (F-I)_xy
        0.0,  # (F-I)_yx
        0.0,  # (F-I)_yy
    ]
)

# residual force
r = fext - fint

# partition in independent and dependent part
# - stiffness matrix
Kii = K[np.ix_(iii, iii)]
Kid = K[np.ix_(iii, iid)]
Kdi = K[np.ix_(iid, iii)]
Kdd = K[np.ix_(iid, iid)]
# - residual force
ri = r[iii]
rd = r[iid]

# condense system
Kii = Kii + np.dot(Kid, Cdi) + np.dot(Cdi.T, Kdi) + np.dot(Cdi.T, np.dot(Kdd, Cdi))
ri = ri + np.dot(Cdi.T, rd)

# partition
# - stiffness matrix
Kuu = Kii[np.ix_(iiu, iiu)]
Kup = Kii[np.ix_(iiu, iip)]
Kpu = Kii[np.ix_(iip, iiu)]
Kpp = Kii[np.ix_(iip, iip)]
# - residual force
ru = ri[iiu]

# solve for unknown displacement DOFs
uu = np.linalg.solve(Kuu, ru - Kup.dot(up))

# assemble
ui = np.empty(iii.size)
ui[iiu] = uu
ui[iip] = up

# reconstruct
# - dependent displacement DOFs
ud = np.dot(Cdi, ui)
# - assemble
u = np.empty(ndof)
u[iii] = ui
u[iid] = ud

# convert to nodal displacements
disp = u[dofs]

# ==================================================================================================
# strain from nodal displacements, stress from constitutive response
# ==================================================================================================

# zero-initialise strain tensor per integration point
# (plain strain -> all z-components are not written and should remain zero at all times)
Eps = np.zeros((nelem, nip, 3, 3))

# loop over nodes
for e, nodes in enumerate(conn):
    # nodal displacements
    ue = disp[nodes, :]

    # loop over integration points
    for q, (w, xi) in enumerate(zip(W, Xi)):
        # alias integration point values
        dNdxe = dNx[e, q, :, :]

        # displacement gradient
        gradu = np.einsum("mi,mj->ij", dNdxe, ue)

        # compute strain tensor, and store per integration point
        # (use plain strain to convert 2-d to 3-d tensor)
        Eps[e, q, :2, :2] = 0.5 * (gradu + gradu.T)

# constitutive response: strain tensor -> stress tensor (per integration point)
Sig = ddot42(C4, Eps)

# ==================================================================================================
# internal force from stress, reaction (external) force from internal force on fixed nodes
# ==================================================================================================

# allocate internal force
fint = np.zeros(ndof)

# loop over nodes
for e, nodes in enumerate(conn):
    # allocate element internal force
    fe = np.zeros(nne * ndim)

    # loop over integration points
    for q, (w, xi) in enumerate(zip(W, Xi)):
        # alias integration point values
        # (plane strain: stress in z-direction irrelevant)
        sig = Sig[e, q, :2, :2]
        dNdxe = dNx[e, q, :, :]
        dV = vol[e, q]

        # assemble to element internal force
        #   fe[m*nd+j] += dNdx[m,i] * sig[i,j] * dV
        fe += (np.einsum("mi,ij->mj", dNdxe, sig) * dV).reshape(nne * ndim)

    # assemble to global stiffness matrix
    iie = dofs[nodes, :].ravel()
    fint[iie] += fe

# residual force
# - all DOFs
r = fext - fint
# - partition [independent, dependent]
ri = r[iii]
rd = r[iid]
# - condense system: apply tying relations
ri = ri + np.dot(Cdi.T, rd)

# reaction force
fext[iip] = -ri[iip]

# ==================================================================================================
# plot
# ==================================================================================================

fig, axes = plt.subplots(ncols=2, figsize=(16, 8))

# reconstruct external force as nodal vectors (for all nodes)
fext = fext[dofs]

# plot original nodal positions + displacement field as arrows
axes[0].scatter(coor[:, 0], coor[:, 1])
axes[0].quiver(coor[:, 0], coor[:, 1], disp[:, 0], disp[:, 1])

# plot new nodal positions + external (reaction) force field as arrows
axes[1].scatter((coor + disp)[:, 0], (coor + disp)[:, 1])
axes[1].quiver((coor + disp)[:, 0], (coor + disp)[:, 1], fext[:, 0], fext[:, 1], scale=0.1)

# fix axes limits
for axis in axes:
    axis.set_xlim([-0.4, 1.5])
    axis.set_ylim([-0.4, 1.5])

plt.show()

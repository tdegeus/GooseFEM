import sys
import GooseFEM
import GMatElastic
import numpy as np

# mesh
# ----

# define mesh
mesh = GooseFEM.Mesh.Quad4.Regular(5, 5)

# mesh dimensions
nelem = mesh.nelem()
nne = mesh.nne()
ndim = mesh.ndim()

# mesh definitions
coor = mesh.coor()
conn = mesh.conn()
dofs = mesh.dofs()

# node sets
nodesLft = mesh.nodesLeftEdge()
nodesRgt = mesh.nodesRightEdge()
nodesTop = mesh.nodesTopEdge()
nodesBot = mesh.nodesBottomEdge()

# fixed displacements DOFs
# ------------------------

iip = np.concatenate((
    dofs[nodesRgt, 0],
    dofs[nodesTop, 1],
    dofs[nodesLft, 0],
    dofs[nodesBot, 1]
))

# simulation variables
# --------------------

# vector definition
vector = GooseFEM.VectorPartitioned(conn, dofs, iip)

# allocate system matrix
K = GooseFEM.MatrixPartitioned(conn, dofs, iip)
Solver = GooseFEM.MatrixPartitionedSolver()

# nodal quantities
disp = np.zeros(coor.shape)
fint = np.zeros(coor.shape)
fext = np.zeros(coor.shape)
fres = np.zeros(coor.shape)

# DOF values
fext_p = np.zeros(vector.nnp())
fres_u = np.zeros(vector.nnu())
fext_p = np.zeros(vector.nnp())

# element vectors
ue = np.empty((nelem, nne, ndim))
fe = np.empty((nelem, nne, ndim))
Ke = np.empty((nelem, nne * ndim, nne * ndim))

# element/material definition
# ---------------------------

# element definition
elem = GooseFEM.Element.Quad4.QuadraturePlanar(vector.AsElement(coor))
nip = elem.nip()

# material definition
mat = GMatElastic.Cartesian3d.Array2d([nelem, nip], 1.0, 1.0)

# integration point tensors
Eps = np.empty((nelem, nip, 3, 3))
Sig = np.empty((nelem, nip, 3, 3))
C = np.empty((nelem, nip, 3, 3, 3, 3))

# solve
# -----

# strain
ue = vector.AsElement(disp)
Eps = elem.SymGradN_vector(ue)

# stress & tangent
mat.setStrain(Eps)
Sig = mat.Stress()
C = mat.Tangent()

# internal force
fe = elem.Int_gradN_dot_tensor2_dV(Sig)
fint = vector.AssembleNode(fe)

# stiffness matrix
Ke = elem.Int_gradN_dot_tensor4_dot_gradNT_dV(C)
K.assemble(Ke)

# set fixed displacements
u_p = np.concatenate((
    +0.1 * np.ones(nodesRgt.size),
    -0.1 * np.ones(nodesTop.size),
    np.zeros(nodesLft.size),
    np.zeros(nodesBot.size)
))

# residual
fres = fext - fint

# partition
fres_u = vector.AsDofs_u(fres)

# solve
u_u = Solver.Solve_u(K, fres_u, u_p)

# assemble to nodal vector
disp = vector.AsNode(u_u, u_p)

# post-process
# ------------

# compute strain and stress
ue = vector.AsElement(disp)
Eps = elem.SymGradN_vector(ue)
mat.setStrain(Eps)
Sig = mat.Stress()

# internal force
fe = elem.Int_gradN_dot_tensor2_dV(Sig)
fint = vector.AssembleNode(fe)

# apply reaction force
fext_p = vector.AsDofs_p(fint)

# residual
fres = fext - fint

# partition
fres_u = vector.AsDofs_u(fres)

# print residual
print(np.sum(np.abs(fres_u)) / np.sum(np.abs(fext_p)))

# average stress per node
dV = elem.AsTensor(2, elem.dV())
Sig = np.average(Sig, weights=dV, axis=1)

# skip plot with "--no-plot" command line argument
# ------------------------------------------------

if len(sys.argv) == 2:
    if sys.argv[1] == "--no-plot":
        sys.exit(0)

# plot
# ----

import matplotlib.pyplot as plt
import GooseMPL as gplt

plt.style.use(['goose', 'goose-latex'])

# extract dimension

nelem = conn.shape[0]

# tensor products

def ddot22(A2, B2):
    return np.einsum('eij, eji -> e', A2, B2)

def ddot42(A4, B2):
    return np.einsum('eijkl, elk -> eij', A4, B2)

def dyad22(A2, B2):
    return np.einsum('eij, ekl -> eijkl', A2, B2)

# identity tensors

i = np.eye(3)
I = np.einsum('ij, e', i, np.ones([nelem]))
I4 = np.einsum('ijkl, e -> eijkl', np.einsum('il, jk', i, i), np.ones([nelem]))
I4rt = np.einsum('ijkl, e -> eijkl', np.einsum('ik,jl', i, i), np.ones([nelem]))
I4s = 0.5 * (I4 + I4rt)
II = dyad22(I, I)
I4d = I4s - II / 3.0

# compute equivalent stress

Sigd = ddot42(I4d, Sig)
sigeq = np.sqrt(3.0 / 2.0 * ddot22(Sigd, Sigd))

# plot

fig, ax = plt.subplots()
gplt.patch(coor=coor + disp, conn=conn, cindex=sigeq, cmap='jet', axis=ax, clim=(0, 0.1))
gplt.patch(coor=coor, conn=conn, linestyle='--', axis=ax)
plt.savefig('plot.pdf')
plt.close()

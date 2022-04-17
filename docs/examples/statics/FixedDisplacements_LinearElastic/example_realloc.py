import sys

import GMatElastic
import GooseFEM
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose", "goose-latex"])

# mesh
# ----

# define mesh
mesh = GooseFEM.Mesh.Quad4.Regular(5, 5)

# mesh dimensions
nelem = mesh.nelem()
nne = mesh.nne()
ndim = mesh.ndim()

# mesh definition, displacement, external forces
coor = mesh.coor()
conn = mesh.conn()
dofs = mesh.dofs()
disp = np.zeros_like(coor)
fext = np.zeros_like(coor)

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
    dofs[nodesBot, 1],
))

# simulation variables
# --------------------

# vector definition
vector = GooseFEM.VectorPartitioned(conn, dofs, iip)

# allocate system matrix
K = GooseFEM.MatrixPartitioned(conn, dofs, iip)
Solver = GooseFEM.MatrixPartitionedSolver()

# element/material definition
# ---------------------------

# element definition
elem = GooseFEM.Element.Quad4.QuadraturePlanar(vector.AsElement(coor))
nip = elem.nip()

# material definition
mat = GMatElastic.Cartesian3d.Array2d([nelem, nip], 1.0, 1.0)

# solve
# -----

# strain
Eps = elem.SymGradN_vector(vector.AsElement(disp))

# stress & tangent
mat.setStrain(Eps)
Sig = mat.Stress()
C = mat.Tangent()

# internal force
fint = vector.AssembleNode(elem.Int_gradN_dot_tensor2_dV(Sig))

# stiffness matrix
K.assemble(elem.Int_gradN_dot_tensor4_dot_gradNT_dV(C))

# set fixed displacements
disp[nodesRgt, 0] = +0.1
disp[nodesTop, 1] = -0.1
disp[nodesLft, 0] = 0.0
disp[nodesBot, 1] = 0.0

# residual
fres = fext - fint

# solve
disp = Solver.Solve(K, fres, disp)

# post-process
# ------------

# compute strain and stress
Eps = elem.SymGradN_vector(vector.AsElement(disp))
mat.setStrain(Eps)
Sig = mat.Stress()

# internal force
fint = vector.AssembleNode(elem.Int_gradN_dot_tensor2_dV(Sig))

# apply reaction force
fext = vector.Copy_p(fint, fext)

# residual
fres = fext - fint

# print residual
print(np.sum(np.abs(fres)) / np.sum(np.abs(fext)))

# average stress per element
dV = elem.AsTensor(2, elem.dV())
Sig = np.average(Sig, weights=dV, axis=1)

# skip plot with "--no-plot" command line argument
# ------------------------------------------------

if len(sys.argv) == 2:
    if sys.argv[1] == "--no-plot":
        sys.exit(0)

# plot
# ----

# extract dimension

nelem = conn.shape[0]

# tensor products

def ddot22(A2, B2):
    return np.einsum("eij, eji -> e", A2, B2)


def ddot42(A4, B2):
    return np.einsum("eijkl, elk -> eij", A4, B2)


def dyad22(A2, B2):
    return np.einsum("eij, ekl -> eijkl", A2, B2)


# identity tensors

i = np.eye(3)
I2 = np.einsum("ij, e", i, np.ones([nelem]))
I4 = np.einsum("ijkl, e -> eijkl", np.einsum("il, jk", i, i), np.ones([nelem]))
I4rt = np.einsum("ijkl, e -> eijkl", np.einsum("ik,jl", i, i), np.ones([nelem]))
I4s = 0.5 * (I4 + I4rt)
II = dyad22(I2, I2)
I4d = I4s - II / 3.0

# compute equivalent stress

Sigd = ddot42(I4d, Sig)
sigeq = np.sqrt(3.0 / 2.0 * ddot22(Sigd, Sigd))

# plot

fig, ax = plt.subplots()
gplt.patch(coor=coor + disp, conn=conn, cindex=sigeq, cmap="jet", axis=ax, clim=(0, 0.1))
gplt.patch(coor=coor, conn=conn, linestyle="--", axis=ax)
plt.savefig("plot.pdf")
plt.close()

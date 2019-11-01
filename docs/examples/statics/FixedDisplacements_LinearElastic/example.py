import GooseFEM
import GMatElastic
import numpy as np

# mesh
# ----

# define mesh
mesh = GooseFEM.Mesh.Quad4.Regular(5,5)

# mesh dimensions
nelem = mesh.nelem()
nne   = mesh.nne()
ndim  = mesh.ndim()

# mesh definitions
coor = mesh.coor()
conn = mesh.conn()
dofs = mesh.dofs()

# node sets
nodesLeft   = mesh.nodesLeftEdge()
nodesRight  = mesh.nodesRightEdge()
nodesTop    = mesh.nodesTopEdge()
nodesBottom = mesh.nodesBottomEdge()

# fixed displacements DOFs
# ------------------------

iip = np.concatenate((
  dofs[nodesRight , 0],
  dofs[nodesTop   , 1],
  dofs[nodesLeft  , 0],
  dofs[nodesBottom, 1]
))

# simulation variables
# --------------------

# vector definition
vector = GooseFEM.VectorPartitioned(conn, dofs, iip)

# allocate system matrix
K = GooseFEM.MatrixPartitioned(conn, dofs, iip)

# nodal quantities
disp = np.zeros(coor.shape)
fint = np.zeros(coor.shape)
fext = np.zeros(coor.shape)
fres = np.zeros(coor.shape)

# element vectors
ue = np.empty((nelem, nne, ndim))
fe = np.empty((nelem, nne, ndim))
Ke = np.empty((nelem, nne*ndim, nne*ndim))

# element/material definition
# ---------------------------

# element definition
elem = GooseFEM.Element.Quad4.QuadraturePlanar(vector.AsElement(coor))
nip = elem.nip()

# material definition
mat = GMatElastic.Cartesian3d.Matrix(nelem, nip, 1., 1.)

# integration point tensors
d = 3
Eps = np.empty((nelem, nip, d, d      ))
Sig = np.empty((nelem, nip, d, d      ))
C   = np.empty((nelem, nip, d, d, d, d))

# solve
# -----

# strain
ue = vector.AsElement(disp)
Eps = elem.SymGradN_vector(ue)

# stress & tangent
(Sig, C) = mat.Tangent(Eps)

# internal force
fe = elem.Int_gradN_dot_tensor2_dV(Sig)
fint = vector.AssembleNode(fe)

# stiffness matrix
Ke = elem.Int_gradN_dot_tensor4_dot_gradNT_dV(C)
K.assemble(Ke)

# set fixed displacements
disp[nodesRight , 0] = +0.1
disp[nodesTop   , 1] = -0.1
disp[nodesLeft  , 0] =  0.0
disp[nodesBottom, 1] =  0.0

# residual
fres = fext - fint

# solve
disp = K.Solve(fres, disp)

# post-process
# ------------

# compute strain and stress
ue = vector.AsElement(disp)
Eps = elem.SymGradN_vector(ue)
Sig = mat.Stress(Eps)

# internal force
fe = elem.Int_gradN_dot_tensor2_dV(Sig)
fint = vector.AssembleNode(fe)

# apply reaction force
fext = vector.Copy_p(fint, fext)

# residual
fres = fext - fint

# print residual
print(np.sum(np.abs(fres)) / np.sum(np.abs(fext)))

# average stress per node
dV = elem.DV(2)
Sig = np.average(Sig, weights=dV, axis=1)

# plot
# ----

import matplotlib.pyplot as plt
import GooseMPL as gplt

plt.style.use(['goose', 'goose-latex'])

# extract dimension
nelem = conn.shape[0]

# tensor products
ddot22 = lambda A2,B2: np.einsum('eij  ,eji->e    ',A2,B2)
ddot42 = lambda A4,B2: np.einsum('eijkl,elk->eij  ',A4,B2)
dyad22 = lambda A2,B2: np.einsum('eij  ,ekl->eijkl',A2,B2)

# identity tensor (single tensor)
i    = np.eye(3)

# identity tensors (grid)
I    = np.einsum('ij  ,e'       ,                  i   ,np.ones([nelem]))
I4   = np.einsum('ijkl,e->eijkl',np.einsum('il,jk',i,i),np.ones([nelem]))
I4rt = np.einsum('ijkl,e->eijkl',np.einsum('ik,jl',i,i),np.ones([nelem]))
I4s  = (I4+I4rt)/2.
II   = dyad22(I,I)
I4d  = I4s-II/3.

# compute equivalent stress
Sigd  = ddot42(I4d, Sig)
sigeq = np.sqrt(3./2.*ddot22(Sigd,Sigd))

# plot
fig, ax = plt.subplots()
gplt.patch(coor=coor+disp, conn=conn, cindex=sigeq, cmap='jet', axis=ax, clim=(0,0.1))
gplt.patch(coor=coor, conn=conn, linestyle='--', axis=ax)
plt.savefig('plot.pdf')
plt.close()

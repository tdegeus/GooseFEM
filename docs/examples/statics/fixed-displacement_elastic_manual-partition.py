import argparse
import sys

import GMatElastic
import GooseFEM
import numpy as np

# mesh
# ----

# define mesh
mesh = GooseFEM.Mesh.Quad4.Regular(5, 5)

# mesh dimensions
nelem = mesh.nelem
nne = mesh.nne
ndim = mesh.ndim

# mesh definition, displacement, external forces
coor = mesh.coor()
conn = mesh.conn()
dofs = mesh.dofs()
disp = np.zeros_like(coor)
fext = np.zeros_like(coor)

# list of prescribed DOFs
iip = np.concatenate(
    (
        dofs[mesh.nodesRightEdge(), 0],
        dofs[mesh.nodesTopEdge(), 1],
        dofs[mesh.nodesLeftEdge(), 0],
        dofs[mesh.nodesBottomEdge(), 1],
    )
)

# simulation variables
# --------------------

# vector definition
vector = GooseFEM.VectorPartitioned(conn, dofs, iip)

# allocate system matrix
K = GooseFEM.MatrixPartitioned(conn, dofs, iip)
Solver = GooseFEM.MatrixPartitionedSolver()

# element definition
elem = GooseFEM.Element.Quad4.QuadraturePlanar(vector.AsElement(coor))
nip = elem.nip

# material definition
# -------------------

mat = GMatElastic.Cartesian3d.Elastic2d(K=np.ones([nelem, nip]), G=np.ones([nelem, nip]))

# solve
# -----

# strain
ue = vector.AsElement(disp)
elem.symGradN_vector(ue, mat.Eps)
mat.refresh()

# internal force
fe = elem.Int_gradN_dot_tensor2_dV(mat.Sig)
fint = vector.AssembleNode(fe)

# stiffness matrix
Ke = elem.Int_gradN_dot_tensor4_dot_gradNT_dV(mat.C)
K.assemble(Ke)

# set fixed displacements
u_p = np.concatenate(
    (
        +0.2 * np.ones(mesh.nodesRightEdge().size),
        -0.2 * np.ones(mesh.nodesTopEdge().size),
        np.zeros(mesh.nodesLeftEdge().size),
        np.zeros(mesh.nodesBottomEdge().size),
    )
)

# residual
fres = fext - fint

# partition
fres_u = vector.AsDofs_u(fres)

# solve
u_u = Solver.Solve_u(K, fres_u, u_p)

# assemble to nodal vector
u = np.empty(vector.shape_dofval())
u[vector.iiu] = u_u
u[vector.iip] = u_p
disp = vector.AsNode(u)

# post-process
# ------------

# strain
vector.asElement(disp, ue)
elem.symGradN_vector(ue, mat.Eps)
mat.refresh()

# internal force
elem.int_gradN_dot_tensor2_dV(mat.Sig, fe)
vector.assembleNode(fe, fint)

# apply reaction force
vector.copy_p(fint, fext)

# residual
fres = fext - fint

# print residual
assert np.isclose(np.sum(np.abs(fres)) / np.sum(np.abs(fext)), 0)

# plot
# ----

parser = argparse.ArgumentParser()
parser.add_argument("--plot", action="store_true", help="Plot result")
parser.add_argument("--save", type=str, help="Save plot (plot not shown)")
args = parser.parse_args(sys.argv[1:])

if args.plot:

    import GooseMPL as gplt
    import matplotlib.pyplot as plt

    plt.style.use(["goose", "goose-latex"])

    # average equivalent stress per element
    dV = elem.AsTensor(2, elem.dV)
    Sigav = np.average(mat.Sig, weights=dV, axis=1)
    sigeq = GMatElastic.Cartesian3d.Sigeq(Sigav)

    # plot
    fig, ax = plt.subplots()
    gplt.patch(coor=coor + disp, conn=conn, cindex=sigeq, cmap="jet", axis=ax, clim=(0, 0.1))
    gplt.patch(coor=coor, conn=conn, linestyle="--", axis=ax)

    # optional save
    if args.save is not None:
        fig.savefig(args.save)
    else:
        plt.show()

    plt.close(fig)

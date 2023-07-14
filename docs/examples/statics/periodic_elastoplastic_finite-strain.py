import argparse
import sys

import GMatElastoPlasticFiniteStrainSimo
import GMatTensor
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

# mesh definitions
coor = mesh.coor
conn = mesh.conn
dofs = mesh.dofs

# periodicity and fixed displacements DOFs
# ----------------------------------------

# add control nodes
control = GooseFEM.Tyings.Control(coor, dofs)
coor = control.coor
dofs = control.dofs
control_dofs = control.controlDofs
control_nodes = control.controlNodes

# extract fixed DOFs:
# - all control nodes: to prescribe the deformation gradient
# - one node of the mesh: to remove rigid body modes
iip = np.concatenate((control_dofs.ravel(), dofs[mesh.nodesOrigin(), :].ravel))

# get DOF-tyings, reorganise system
tyings = GooseFEM.Tyings.Periodic(coor, dofs, control_dofs, mesh.nodesPeriodic, iip)
dofs = tyings.dofs

# simulation variables
# --------------------

# vector definition
vector = GooseFEM.VectorPartitionedTyings(conn, dofs, tyings.Cdu(), tyings.Cdp(), tyings.Cdi())

# element definition
elem0 = GooseFEM.Element.Quad4.QuadraturePlanar(vector.AsElement(coor))
elem = GooseFEM.Element.Quad4.QuadraturePlanar(vector.AsElement(coor))
nip = elem.nip

# nodal quantities
disp = np.zeros_like(coor)  # nodal displacement
du = np.zeros_like(coor)  # iterative displacement update
fint = np.zeros_like(coor)  # internal force
fext = np.zeros_like(coor)  # external force

# element vectors / matrix
ue = np.empty([nelem, nne, ndim])
fe = np.empty([nelem, nne, ndim])
Ke = np.empty([nelem, nne * ndim, nne * ndim])

# DOF values
Fext = np.zeros([tyings.nni])
Fint = np.zeros([tyings.nni])

# material definition
# -------------------

kappa = np.ones([nelem, nip])
mu = np.ones([nelem, nip])
tauy0 = 0.05 * np.ones([nelem, nip])
H = np.ones([nelem, nip])

elmat = mesh.elementgrid()
ehard = elmat[:2, :2].ravel()
tauy0[ehard, :] = 100

mat = GMatElastoPlasticFiniteStrainSimo.Cartesian3d.LinearHardening2d(kappa, mu, tauy0, H)

# solve
# -----

# allocate system matrix
K = GooseFEM.MatrixPartitionedTyings(conn, dofs, tyings.Cdu(), tyings.Cdp())
Solver = GooseFEM.MatrixPartitionedTyingsSolver()

# array of unit tensor
I2 = GMatTensor.Cartesian3d.Array2d(mat.shape).I2

# some shear strain history
dgamma = 0.001 * np.ones([101])
dgamma[0] = 0

# # loop over increments
for inc in range(dgamma.size):
    # update history
    mat.increment()

    for iiter in range(100):
        # deformation gradient
        vector.asElement(disp, ue)
        elem0.symGradN_vector(ue, mat.F)
        mat.F += I2
        mat.refresh()

        # internal force
        elem.int_gradN_dot_tensor2_dV(mat.Sig, fe)
        vector.assembleNode(fe, fint)

        # stiffness matrix
        elem.int_gradN_dot_tensor4_dot_gradNT_dV(mat.C, Ke)
        K.assemble(Ke)

        # residual
        fres = fext - fint

        # check for convergence (skip the zeroth iteration, as the residual still vanishes)
        if iiter > 0:
            # - internal/external force as DOFs (account for periodicity)
            vector.asDofs_i(fext, Fext)
            vector.asDofs_i(fint, Fint)
            # - extract reaction force
            vector.copy_p(Fint, Fext)
            # - norm of the residual and the reaction force
            nfres = np.sum(np.abs(Fext - Fint))
            nfext = np.sum(np.abs(Fext))
            # - relative residual, for convergence check
            if nfext:
                res = nfres / nfext
            else:
                res = nfres
            # - check for convergence
            if res < 1e-5:
                break
            # - safe-guard from infinite loop
            if iiter > 20:
                raise OSError("Maximal number of iterations exceeded")

        # initialise displacement update
        du.fill(0.0)

        # set fixed displacements
        if iiter == 0:
            du[control_nodes[0], 1] = dgamma[inc]

        # solve
        Solver.solve(K, fres, du)

        # add displacement update
        disp += du

        # update shape functions
        elem.update_x(vector.AsElement(coor + disp))

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

    # strain
    vector.asElement(disp, ue)
    elem0.symGradN_vector(ue, mat.F)
    mat.F += I2
    mat.refresh()

    # average equivalent stress per element
    dV = elem.AsTensor(2, elem.dV)
    Sigav = np.average(mat.Sig, weights=dV, axis=1)
    sigeq = GMatElastoPlasticFiniteStrainSimo.Cartesian3d.Sigeq(Sigav)

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

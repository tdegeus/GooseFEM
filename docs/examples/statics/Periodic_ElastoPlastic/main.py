import GMatElastoPlastic.Cartesian3d as GMat
import GooseFEM
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose", "goose-autolayout"])

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
elmat = mesh.elementgrid()

# periodicity and fixed displacements DOFs
# ----------------------------------------

# add control nodes
control = GooseFEM.Tyings.Control(coor, dofs)
coor = control.coor()
dofs = control.dofs()
control_dofs = control.controlDofs()
control_nodes = control.controlNodes()

# extract fixed DOFs:
# - all control nodes: to prescribe the deformation gradient
# - one node of the mesh: to remove rigid body modes
iip = np.concatenate((control_dofs.ravel(), dofs[mesh.nodesOrigin(), :].ravel()))


# get DOF-tyings, reorganise system
tyings = GooseFEM.Tyings.Periodic(coor, dofs, control_dofs, mesh.nodesPeriodic(), iip)
dofs = tyings.dofs()

# simulation variables
# --------------------

# vector definition:
# provides methods to switch between dofval/nodeval/elemvec, or to manipulate a part of them
vector = GooseFEM.VectorPartitionedTyings(conn, dofs, tyings.Cdu(), tyings.Cdp(), tyings.Cdi())

# nodal quantities
disp = np.zeros_like(coor)  # nodal displacement
du = np.zeros_like(coor)  # iterative displacement update
fint = np.zeros_like(coor)  # internal force
fext = np.zeros_like(coor)  # external force
fres = np.zeros_like(coor)  # residual force

# element vectors / matrix
ue = np.empty([nelem, nne, ndim])
fe = np.empty([nelem, nne, ndim])
Ke = np.empty([nelem, nne * ndim, nne * ndim])

# DOF values
Fext = np.zeros([tyings.nni()])
Fint = np.zeros([tyings.nni()])

# element/material definition
# ---------------------------

# FEM quadrature
elem = GooseFEM.Element.Quad4.QuadraturePlanar(vector.AsElement(coor))
nip = elem.nip()
dV = elem.AsTensor(2, elem.dV())

# material model
# even though the problem is 2-d, the material model is 3-d, plane strain is implicitly assumed
mat = GMat.Array2d([nelem, nip])
tdim = 3

# some artificial material definition
ehard = elmat[:2, :2]
Ihard = np.zeros([nelem, nip], dtype=bool)
Ihard[ehard, :] = True
Isoft = ~Ihard

mat.setLinearHardening(Isoft, 1.0, 1.0, 0.05, 0.05)
mat.setElastic(Ihard, 1.0, 1.0)

# solve
# -----

# allocate tensors
Eps = np.empty([nelem, nip, tdim, tdim])
Sig = np.empty([nelem, nip, tdim, tdim])
C = np.empty([nelem, nip, tdim, tdim, tdim, tdim])

# allocate system matrix
K = GooseFEM.MatrixPartitionedTyings(conn, dofs, tyings.Cdu(), tyings.Cdp())
Solver = GooseFEM.MatrixPartitionedTyingsSolver()

# # allocate internal variables
# double res

# some shear strain history
dgamma = 0.001 * np.ones([101])
dgamma[0] = 0.0

# # loop over increments
for inc in range(dgamma.size):

    # update history
    mat.increment()

    for iiter in range(100):

        # strain
        vector.asElement(disp, ue)
        elem.symGradN_vector(ue, Eps)

        # stress & tangent
        mat.setStrain(Eps)
        Sig = mat.Stress()  # todo, replace with : mat.stress(Sig)
        C = mat.Tangent()  # todo, replace with : mat.tangent(C)

        # internal force
        elem.int_gradN_dot_tensor2_dV(Sig, fe)
        vector.assembleNode(fe, fint)

        # stiffness matrix
        elem.int_gradN_dot_tensor4_dot_gradNT_dV(C, Ke)
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
            # - print progress to screen
            print("inc = ", inc, "iiter = ", iiter, "res = ", res)
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
        du = Solver.Solve(K, fres, du)

        # add displacement update
        disp += du


vector.asElement(disp, ue)
elem.symGradN_vector(ue, Eps)
mat.setStrain(Eps)
Sig = mat.Stress()  # todo, replace with mat.stress(Sig)
sigeq = GMat.Sigeq(np.mean(Sig, axis=1))

fig, ax = plt.subplots()
gplt.patch(coor=coor + disp, conn=conn, cindex=sigeq, cmap="jet")
plt.show()
# fig.savefig("")
plt.close(fig)

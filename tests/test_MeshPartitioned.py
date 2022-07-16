import unittest

import GooseFEM
import numpy as np


class Test_MeshQuad4(unittest.TestCase):
    """ """

    def test_basic(self):

        mesh = GooseFEM.Mesh.Quad4.Regular(2, 2)

        nne = mesh.nne
        ndim = mesh.ndim
        nelem = mesh.nelem
        nnode = mesh.nnode
        dofs = mesh.dofs()
        iip = [0, 5, 7, 13]

        a = np.empty([nelem, nne * ndim, nne * ndim])
        x = np.random.random([nnode * ndim])

        for e in range(nelem):
            ae = np.random.random([nne * ndim, nne * ndim])
            a[e, ...] = 0.5 * (ae + ae.T)

        A = GooseFEM.MatrixPartitioned(mesh.conn(), dofs, iip)
        B = GooseFEM.Matrix(mesh.conn(), dofs)

        A.assemble(a)
        B.assemble(a)

        self.assertTrue(np.allclose(A.Todense(), B.Todense()))

        Solver = GooseFEM.MatrixPartitionedSolver()
        xc = Solver.Solve(A, A.Dot(x), x)
        self.assertTrue(np.allclose(x, xc))


if __name__ == "__main__":

    unittest.main()

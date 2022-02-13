import unittest

import GooseFEM
import numpy as np


class Test_Mesh(unittest.TestCase):
    """
    Simple mesh operations.
    """

    def test_center_of_gravity(self):
        """
        Get the center of gravity of a mesh.
        """

        mesh = GooseFEM.Mesh.Quad4.Regular(3, 3)
        coor = mesh.coor()
        conn = mesh.conn()

        self.assertTrue(np.allclose([1.5, 1.5], GooseFEM.Mesh.center_of_gravity(coor, conn)))

    def test_coordination_dofs(self):

        mesh = GooseFEM.Mesh.Quad4.Regular(4, 4)
        dofs = mesh.dofsPeriodic()
        real = np.array(
            [
                [4, 4],
                [2, 2],
                [2, 2],
                [2, 2],
                [2, 2],
                [1, 1],
                [1, 1],
                [1, 1],
                [2, 2],
                [1, 1],
                [1, 1],
                [1, 1],
                [2, 2],
                [1, 1],
                [1, 1],
                [1, 1],
            ]
        ).ravel()

        ret = np.array(GooseFEM.Mesh.coordination(dofs))

        self.assertTrue(np.all(real == ret))

    def test_coordination_node2dof(self):

        mesh = GooseFEM.Mesh.Quad4.Regular(4, 4)
        dofs = mesh.dofsPeriodic()
        real = [
            [0, 4, 20, 24],
            [0, 4, 20, 24],
            [1, 21],
            [1, 21],
            [2, 22],
            [2, 22],
            [3, 23],
            [3, 23],
            [5, 9],
            [5, 9],
            [6],
            [6],
            [7],
            [7],
            [8],
            [8],
            [10, 14],
            [10, 14],
            [11],
            [11],
            [12],
            [12],
            [13],
            [13],
            [15, 19],
            [15, 19],
            [16],
            [16],
            [17],
            [17],
            [18],
            [18],
        ]

        ret = GooseFEM.Mesh.node2dof(dofs)

        self.assertEqual(len(real), len(ret))

        for i in range(len(real)):
            self.assertTrue(np.all(np.array(real[i]) == np.array(ret[i])))

    def test_coordination_nodaltyings(self):

        mesh = GooseFEM.Mesh.Quad4.Regular(4, 4)
        dofs = mesh.dofsPeriodic()
        real = [
            [0, 4, 20, 24],
            [1, 21],
            [2, 22],
            [3, 23],
            [0, 4, 20, 24],
            [5, 9],
            [6],
            [7],
            [8],
            [5, 9],
            [10, 14],
            [11],
            [12],
            [13],
            [10, 14],
            [15, 19],
            [16],
            [17],
            [18],
            [15, 19],
            [0, 4, 20, 24],
            [1, 21],
            [2, 22],
            [3, 23],
            [0, 4, 20, 24],
        ]

        ret = GooseFEM.Mesh.nodaltyings(dofs)

        self.assertEqual(len(real), len(ret))

        for i in range(len(real)):
            self.assertTrue(np.all(np.array(real[i]) == np.array(ret[i])))

    def test_elem2node_periodic(self):

        mesh = GooseFEM.Mesh.Quad4.Regular(4, 4)
        conn = mesh.conn()
        dofs = mesh.dofsPeriodic()

        real = np.array(
            [
                [15, 12, 3, 0],
                [12, 13, 0, 1],
                [13, 14, 1, 2],
                [14, 15, 2, 3],
                [15, 12, 3, 0],
                [3, 0, 7, 4],
                [0, 1, 4, 5],
                [1, 2, 5, 6],
                [2, 3, 6, 7],
                [3, 0, 7, 4],
                [7, 4, 11, 8],
                [4, 5, 8, 9],
                [5, 6, 9, 10],
                [6, 7, 10, 11],
                [7, 4, 11, 8],
                [11, 8, 15, 12],
                [8, 9, 12, 13],
                [9, 10, 13, 14],
                [10, 11, 14, 15],
                [11, 8, 15, 12],
                [15, 12, 3, 0],
                [12, 13, 0, 1],
                [13, 14, 1, 2],
                [14, 15, 2, 3],
                [15, 12, 3, 0],
            ]
        )

        real = np.sort(real, axis=1)
        ret = np.array(GooseFEM.Mesh.elem2node(conn, dofs))

        self.assertTrue(np.all(real == ret))


if __name__ == "__main__":

    unittest.main()

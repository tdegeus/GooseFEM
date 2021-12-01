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


if __name__ == "__main__":

    unittest.main()

import unittest

import GooseFEM as gf
import numpy as np


class Test_MeshQuad4(unittest.TestCase):
    """
    Generation/manipulation of Quad4 meshes.
    """

    def test_FineLayer_replica(self):
        """
        Reconstruct Mesh.Quad4.FineLayer object from existing mesh.
        """

        mesh = gf.Mesh.Quad4.FineLayer(27, 27)
        replica = gf.Mesh.Quad4.FineLayer(mesh.coor(), mesh.conn())
        self.assertTrue(np.allclose(mesh.coor(), replica.coor()))
        self.assertTrue(np.all(np.equal(mesh.conn(), replica.conn())))


if __name__ == "__main__":

    unittest.main()

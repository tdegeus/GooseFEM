import unittest

import GooseFEM as gf
import numpy as np


class Test_MeshQuad4(unittest.TestCase):
    def test_FineLayer_replica(self):

        mesh = gf.Mesh.Quad4.FineLayer(27, 27)
        replica = gf.Mesh.Quad4.FineLayer(mesh.coor(), mesh.conn())
        self.assertTrue(np.allclose(mesh.coor(), replica.coor()))
        self.assertTrue(np.all(np.equal(mesh.conn(), replica.conn())))


if __name__ == "__main__":

    unittest.main()

import unittest

import GooseFEM
import numpy as np


class Test_MeshQuad4(unittest.TestCase):
    """
    Generation/manipulation of Quad4 meshes.
    """

    def test_Regular(self):
        coor = np.array(
            [
                [0, 0],
                [1, 0],
                [2, 0],
                [3, 0],
                [0, 1],
                [1, 1],
                [2, 1],
                [3, 1],
                [0, 2],
                [1, 2],
                [2, 2],
                [3, 2],
                [0, 3],
                [1, 3],
                [2, 3],
                [3, 3],
            ]
        )

        nodes = np.arange(4 * 4).reshape(4, 4)
        conn = np.hstack(
            (
                nodes[:-1, :-1].ravel().reshape(-1, 1),
                nodes[:-1, 1:].ravel().reshape(-1, 1),
                nodes[1:, 1:].ravel().reshape(-1, 1),
                nodes[1:, :-1].ravel().reshape(-1, 1),
            )
        )

        mesh = GooseFEM.Mesh.Quad4.Regular(3, 3)

        self.assertTrue(np.allclose(mesh.coor, coor))
        self.assertTrue(np.all(np.equal(mesh.conn, conn)))

        self.assertTrue(np.all(np.equal(mesh.nodesLeftEdge, nodes[:, 0])))
        self.assertTrue(np.all(np.equal(mesh.nodesRightEdge, nodes[:, -1])))
        self.assertTrue(np.all(np.equal(mesh.nodesBottomEdge, nodes[0, :])))
        self.assertTrue(np.all(np.equal(mesh.nodesTopEdge, nodes[-1, :])))

        self.assertTrue(np.all(np.equal(mesh.nodesLeftOpenEdge, mesh.nodesLeftEdge[1:-1])))
        self.assertTrue(np.all(np.equal(mesh.nodesRightOpenEdge, mesh.nodesRightEdge[1:-1])))
        self.assertTrue(np.all(np.equal(mesh.nodesBottomOpenEdge, mesh.nodesBottomEdge[1:-1])))
        self.assertTrue(np.all(np.equal(mesh.nodesTopOpenEdge, mesh.nodesTopEdge[1:-1])))

        self.assertTrue(np.all(np.equal(mesh.nodesBottomLeftCorner, mesh.nodesBottomEdge[0])))
        self.assertTrue(np.all(np.equal(mesh.nodesBottomRightCorner, mesh.nodesBottomEdge[-1])))
        self.assertTrue(np.all(np.equal(mesh.nodesTopLeftCorner, mesh.nodesTopEdge[0])))
        self.assertTrue(np.all(np.equal(mesh.nodesTopRightCorner, mesh.nodesTopEdge[-1])))

        # alias, needs checking only once
        self.assertTrue(np.all(np.equal(mesh.nodesLeftBottomCorner, mesh.nodesBottomLeftCorner)))
        self.assertTrue(np.all(np.equal(mesh.nodesRightBottomCorner, mesh.nodesBottomRightCorner)))
        self.assertTrue(np.all(np.equal(mesh.nodesLeftTopCorner, mesh.nodesTopLeftCorner)))
        self.assertTrue(np.all(np.equal(mesh.nodesRightTopCorner, mesh.nodesTopRightCorner)))

    def test_FineLayer_replica(self):
        """
        Reconstruct Mesh.Quad4.FineLayer object from existing mesh.
        """

        mesh = GooseFEM.Mesh.Quad4.FineLayer(27, 27)
        replica = GooseFEM.Mesh.Quad4.FineLayer(mesh.coor, mesh.conn)
        self.assertTrue(np.allclose(mesh.coor, replica.coor))
        self.assertTrue(np.all(np.equal(mesh.conn, replica.conn)))


if __name__ == "__main__":
    unittest.main()

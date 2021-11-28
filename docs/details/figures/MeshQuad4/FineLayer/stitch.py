import GooseFEM
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose"])

fig, ax = plt.subplots()

mesh = GooseFEM.Mesh.Quad4.FineLayer(6 * 9, 51)

coor_0 = mesh.coor()
coor_1 = mesh.coor()
coor_2 = mesh.coor()

H = np.max(coor_0[:, 1])

coor_1[:, 1] += H
coor_2[:, 1] += H * 2

stitch = GooseFEM.Mesh.Stitch()
stitch.push_back(coor_0, mesh.conn())
stitch.push_back(coor_1, mesh.conn())
stitch.push_back(coor_2, mesh.conn())

coor = stitch.coor()
conn = stitch.conn()

cindex = np.zeros(conn.shape[0])
cindex[stitch.elemset(np.arange(mesh.nelem()), 0)] = 1
cindex[stitch.elemset(np.arange(mesh.nelem()), 1)] = 2
cindex[stitch.elemset(np.arange(mesh.nelem()), 2)] = 3

left = stitch.nodeset([mesh.nodesLeftEdge(), mesh.nodesLeftEdge(), mesh.nodesLeftEdge()])
right = stitch.nodeset([mesh.nodesRightEdge(), mesh.nodesRightEdge(), mesh.nodesRightEdge()])
top = stitch.nodeset(mesh.nodesTopEdge(), 2)
bottom = stitch.nodeset(mesh.nodesBottomEdge(), 0)
left = np.setdiff1d(np.setdiff1d(left, top), bottom)
right = np.setdiff1d(np.setdiff1d(right, top), bottom)

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)
ax.plot(coor[left, 0], coor[left, 1], marker="o", c="b", ls="none")
ax.plot(coor[right, 0], coor[right, 1], marker="o", c="r", ls="none")
ax.plot(coor[bottom, 0], coor[bottom, 1], marker="o", c="g", ls="none")
ax.plot(coor[top, 0], coor[top, 1], marker="o", c="c", ls="none")

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.show()

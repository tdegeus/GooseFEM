import GooseFEM as gf
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose"])

# --------------------------------------------------------------------------------------------------

mesh = gf.Mesh.Quad4.FineLayer(6 * 9, 51)
coor = mesh.coor()
conn = mesh.conn()

cindex = np.arange(mesh.nelem())

Left = mesh.nodesLeftEdge()
Right = mesh.nodesRightEdge()
Bottom = mesh.nodesBottomEdge()
Top = mesh.nodesTopEdge()

# --------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 10))

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet")

ax.plot(coor[:, 0], coor[:, 1], marker="o", linestyle="none")

ax.plot(coor[Left, 0], coor[Left, 1], marker="o", linestyle="none", color="g")
ax.plot(coor[Right, 0], coor[Right, 1], marker="o", linestyle="none", color="b")
ax.plot(coor[Bottom, 0], coor[Bottom, 1], marker="o", linestyle="none", color="r")
ax.plot(coor[Top, 0], coor[Top, 1], marker="o", linestyle="none", color="y")

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.savefig("example.svg")
plt.close()

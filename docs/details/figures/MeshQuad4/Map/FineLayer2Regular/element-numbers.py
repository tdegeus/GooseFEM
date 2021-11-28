import GooseFEM as gf
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose"])

# --------------------------------------------------------------------------------------------------

fig, axes = plt.subplots(figsize=(10, 10), ncols=2)

# ---

ax = axes[0]

mesh = gf.Mesh.Quad4.FineLayer(6, 18)
coor = mesh.coor()
conn = mesh.conn()

cindex = np.arange(mesh.nelem())

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)

for e in range(mesh.nelem()):
    x = np.mean(coor[conn[e, :], 0])
    y = np.mean(coor[conn[e, :], 1])
    ax.text(x, y, str(e), ha="center", va="center")

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# ---

ax = axes[1]

new_mesh = gf.Mesh.Quad4.Regular(mesh.nelx(), mesh.nely())

coor = new_mesh.coor()
conn = new_mesh.conn()

cindex = np.arange(new_mesh.nelem())

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)

for e in range(new_mesh.nelem()):
    x = np.mean(coor[conn[e, :], 0])
    y = np.mean(coor[conn[e, :], 1])
    ax.text(x, y, str(e), ha="center", va="center")

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# ---

plt.savefig("element-numbers.svg")
plt.show()

import GooseFEM as gf
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose"])

# --------------------------------------------------------------------------------------------------

fig, axes = plt.subplots(figsize=(10, 10), ncols=2)

# ---

ax = axes[0]

mesh = gf.Mesh.Quad4.FineLayer(6 * 3, 18 * 2)
coor = mesh.coor
conn = mesh.conn

cindex = np.random.random(mesh.nelem())

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# ---

ax = axes[1]

mapping = gf.Mesh.Quad4.Map.FineLayer2Regular(mesh)

new_mesh = mapping.getRegularMesh()

coor = new_mesh.coor
conn = new_mesh.conn

c = mapping.mapToRegular(cindex)

im = gplt.patch(coor=coor, conn=conn, cindex=c, cmap="jet", axis=ax)

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# ---

plt.savefig("map.svg")
plt.show()

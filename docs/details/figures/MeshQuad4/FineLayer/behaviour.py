import GooseFEM as gf
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose"])

# --------------------------------------------------------------------------------------------------

fig, axes = plt.subplots(figsize=(18, 6), ncols=3)

# --------------------------------------------------------------------------------------------------

mesh = gf.Mesh.Quad4.FineLayer(6 * 9, 51)
coor = mesh.coor()
conn = mesh.conn()

cindex = np.arange(mesh.nelem())

# --------------------------------------------------------------------------------------------------

ax = axes[0]

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

ax.set_title(r"$n_x = 6 \times 9; n_{elem} = %d$" % mesh.nelem())

# --------------------------------------------------------------------------------------------------

mesh = gf.Mesh.Quad4.FineLayer(6 * 9 + 3, 51)
coor = mesh.coor()
conn = mesh.conn()

cindex = np.arange(mesh.nelem())

# --------------------------------------------------------------------------------------------------

ax = axes[1]

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

ax.set_title(r"$n_x = 6 \times 9 + 3; n_{elem} = %d$" % mesh.nelem())

# --------------------------------------------------------------------------------------------------

mesh = gf.Mesh.Quad4.FineLayer(6 * 9 + 1, 51)
coor = mesh.coor()
conn = mesh.conn()

cindex = np.arange(mesh.nelem())

# --------------------------------------------------------------------------------------------------

ax = axes[2]

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

ax.set_title(r"$n_x = 6 \times 9 + 1; n_{elem} = %d$" % mesh.nelem())

plt.savefig("behaviour.svg")
plt.show()

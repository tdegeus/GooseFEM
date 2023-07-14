import GooseFEM
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose"])

fig, ax = plt.subplots()

mesh = GooseFEM.Mesh.Quad4.FineLayer(6 * 9, 51)
coor = mesh.coor
conn = mesh.conn

cindex = np.zeros(mesh.nelem())

for i in range(mesh.elemrow_nelem().size):
    print(i, mesh.elementsLayer(i))
    cindex[mesh.elementsLayer(i)] = i

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.show()

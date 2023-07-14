import GooseFEM as gf
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose"])

mesh = gf.Mesh.Quad4.FineLayer(12, 38)
coor = mesh.coor
conn = mesh.conn
layer = mesh.elementsMiddleLayer
mid = layer[-1]
select = mesh.elementgrid_around_ravel(mid, 2)

iden = np.zeros(conn.shape[0])
iden[select] = 1
iden[mid] = 2

fig, ax = plt.subplots()

im = gplt.patch(coor=coor, conn=conn, cindex=iden, cmap="Reds", clim=[0, 2])

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

fig.savefig("elementgrid.svg")
plt.close(fig)

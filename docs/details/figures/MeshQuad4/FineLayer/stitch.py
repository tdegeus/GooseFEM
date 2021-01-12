
import GooseFEM
import matplotlib.pyplot as plt
import GooseMPL as gplt
import numpy as np

plt.style.use(['goose'])

fig, ax = plt.subplots()

mesh = GooseFEM.Mesh.Quad4.FineLayer(6 * 9, 51)

coor_a = mesh.coor()
conn_a = mesh.conn()

coor_b = mesh.coor()
conn_b = mesh.conn()

coor_c = mesh.coor()
conn_c = mesh.conn()

coor_b[:, 1] += np.max(coor_a[:, 1])
coor_c[:, 1] += np.max(coor_a[:, 1]) * 2

stitch = GooseFEM.Mesh.Stitch()
stitch.push_back(coor_a, conn_a)
stitch.push_back(coor_b, conn_b)
stitch.push_back(coor_c, conn_c)

coor = stitch.coor()
conn = stitch.conn()

cindex = np.zeros(conn.shape[0])
cindex[stitch.elementset(np.arange(mesh.nelem()), 0)] = 1
cindex[stitch.elementset(np.arange(mesh.nelem()), 1)] = 2
cindex[stitch.elementset(np.arange(mesh.nelem()), 2)] = 3

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap='jet', axis=ax)

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.show()

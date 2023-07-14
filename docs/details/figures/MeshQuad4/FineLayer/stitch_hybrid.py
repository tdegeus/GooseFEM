import GooseFEM
import GooseMPL as gplt
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["goose"])

fig, ax = plt.subplots()

finelayer = GooseFEM.Mesh.Quad4.FineLayer(6 * 9, 51)

coor_l0 = finelayer.coor
coor_l1 = finelayer.coor
coor_l2 = finelayer.coor

Hl = np.max(coor_l0[:, 1])

coor_l1[:, 1] += Hl
coor_l2[:, 1] += Hl * 2

h = finelayer.elemrow_nhx()[0]
nx = finelayer.elemrow_nelem()[0]
ny = np.ceil(((finelayer.nely() - 1) / 2) / h)

regular = GooseFEM.Mesh.Quad4.Regular(int(nx), int(ny), float(h))

coor_r0 = regular.coor
coor_r1 = regular.coor

Hr = np.max(coor_r0[:, 1])

coor_r0[:, 1] -= Hr
coor_r1[:, 1] += 3 * Hl

stitch = GooseFEM.Mesh.Stitch()
stitch.push_back(coor_r0, regular.conn)
stitch.push_back(coor_l0, finelayer.conn)
stitch.push_back(coor_l1, finelayer.conn)
stitch.push_back(coor_l2, finelayer.conn)
stitch.push_back(coor_r1, regular.conn)

coor = stitch.coor
conn = stitch.conn

cindex = np.zeros(conn.shape[0])
cindex[stitch.elemset(np.arange(regular.nelem()), 0)] = 1
cindex[stitch.elemset(np.arange(finelayer.nelem()), 1)] = 2
cindex[stitch.elemset(np.arange(finelayer.nelem()), 2)] = 3
cindex[stitch.elemset(np.arange(finelayer.nelem()), 3)] = 4
cindex[stitch.elemset(np.arange(regular.nelem()), 4)] = 5

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap="jet", axis=ax)

ax.set_aspect("equal")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.show()

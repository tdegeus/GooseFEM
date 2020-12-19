
import GooseFEM as gf
import matplotlib.pyplot as plt
import GooseMPL as gplt
import numpy as np

plt.style.use(['goose'])

mesh = gf.Mesh.Quad4.FineLayer(12, 38)
coor = mesh.coor()
conn = mesh.conn()
layer = mesh.elementsMiddleLayer()
mid = layer[-1]
select = mesh.elementgrid_around_ravel(mid, 2)

I = np.zeros((conn.shape[0]))
I[select] = 1
I[mid] = 2

fig, ax = plt.subplots()

im = gplt.patch(coor=coor, conn=conn, cindex=I, cmap='Reds', clim=[0, 2])

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

fig.savefig('elementgrid.svg')
plt.close(fig)


import GooseFEM as gf
import matplotlib.pyplot as plt
import GooseMPL as gplt
import numpy as np

plt.style.use(['goose'])

# --------------------------------------------------------------------------------------------------

mesh   = gf.Mesh.Quad4.FineLayer(6, 18)
coor   = mesh.coor()
conn   = mesh.conn()
cindex = np.arange(conn.shape[0])

# --------------------------------------------------------------------------------------------------

fig,ax = plt.subplots(figsize=(10,10))

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap='jet')

for e in range(conn.shape[0]):
  x = np.mean(coor[conn[e,:], 0])
  y = np.mean(coor[conn[e,:], 1])
  ax.text(x, y, str(e), ha='center', va='center')

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.savefig('element_numbers.svg')
plt.show()

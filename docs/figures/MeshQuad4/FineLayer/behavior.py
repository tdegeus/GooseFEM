
import GooseFEM as gf
import matplotlib.pyplot as plt
import goosempl as gplt
import numpy as np

plt.style.use(['goose'])

fig = plt.figure(figsize=(18,6))
fig.set_tight_layout(True)

# --------------------------------------------------------------------------------------------------

mesh   = gf.Mesh.Quad4.FineLayer(6*9,51)
coor   = mesh.coor()
conn   = mesh.conn()
cindex = np.arange(conn.shape[0])

print(mesh.nelem())

# --------------------------------------------------------------------------------------------------

ax = fig.add_subplot(1,3,1)

im = gplt.patch(coor=coor,conn=conn,cindex=cindex,cmap='jet')

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# --------------------------------------------------------------------------------------------------

mesh   = gf.Mesh.Quad4.FineLayer(6*9+3,51)
coor   = mesh.coor()
conn   = mesh.conn()
cindex = np.arange(conn.shape[0])

print(mesh.nelem())

# --------------------------------------------------------------------------------------------------

ax = fig.add_subplot(1,3,2)

im = gplt.patch(coor=coor,conn=conn,cindex=cindex,cmap='jet')

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

# --------------------------------------------------------------------------------------------------

mesh   = gf.Mesh.Quad4.FineLayer(6*9+1,51)
coor   = mesh.coor()
conn   = mesh.conn()
cindex = np.arange(conn.shape[0])

print(mesh.nelem())

# --------------------------------------------------------------------------------------------------

ax = fig.add_subplot(1,3,3)

im = gplt.patch(coor=coor,conn=conn,cindex=cindex,cmap='jet')

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.savefig('behavior.svg')
plt.show()

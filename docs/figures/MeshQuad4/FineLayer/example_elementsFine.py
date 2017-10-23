
import GooseFEM as gf
import matplotlib.pyplot as plt
import goosempl as gplt
import numpy as np

plt.style.use(['goose'])

# --------------------------------------------------------------------------------------------------

mesh   = gf.Mesh.Quad4.FineLayer(6*9,51)
coor   = mesh.coor()
conn   = mesh.conn()
cindex = np.zeros((conn.shape[0]))
cindex[mesh.elementsFine()] = np.linspace(1,2,len(mesh.elementsFine()))

# --------------------------------------------------------------------------------------------------

fig,ax = plt.subplots(figsize=(10,10))

im = gplt.patch(coor=coor,conn=conn,cindex=cindex,cmap='Reds')

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.savefig('example_elementsFine.svg')
plt.show()

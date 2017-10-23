
import GooseFEM as gf
import matplotlib.pyplot as plt
import goosempl as gplt
import numpy as np

plt.style.use(['goose'])

# --------------------------------------------------------------------------------------------------

mesh   = gf.Mesh.Quad4.FineLayer(6*9,51)
coor   = mesh.coor()
conn   = mesh.conn()
cindex = np.arange(conn.shape[0])

# --------------------------------------------------------------------------------------------------

fig,ax = plt.subplots(figsize=(10,10))

im = gplt.patch(coor=coor,conn=conn,cindex=cindex,cmap='jet')

ax.plot(coor[:                 ,0],coor[:                 ,1],marker='o',linestyle='none')
ax.plot(coor[mesh.nodesLeft  (),0],coor[mesh.nodesLeft  (),1],marker='o',linestyle='none',color='g')
ax.plot(coor[mesh.nodesRight (),0],coor[mesh.nodesRight (),1],marker='o',linestyle='none',color='b')
ax.plot(coor[mesh.nodesBottom(),0],coor[mesh.nodesBottom(),1],marker='o',linestyle='none',color='r')
ax.plot(coor[mesh.nodesTop   (),0],coor[mesh.nodesTop   (),1],marker='o',linestyle='none',color='y')

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

plt.savefig('example.svg')
plt.show()

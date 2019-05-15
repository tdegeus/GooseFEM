
import GooseFEM          as gf
import matplotlib.pyplot as plt
import GooseMPL          as gplt
import numpy             as np

plt.style.use(['goose'])

# --------------------------------------------------------------------------------------------------

mesh   = gf.Mesh.Quad4.FineLayer(6*9, 51)
coor   = mesh.coor()
conn   = mesh.conn()
cindex = np.arange(conn.shape[0])

# --------------------------------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10,10))

im = gplt.patch(coor=coor, conn=conn, cindex=cindex, cmap='jet')

ax.plot(coor[:,0], coor[:,1], marker='o', linestyle='none')

lft = mesh.nodesLeftOpenEdge()
rgt = mesh.nodesRightOpenEdge()
bot = mesh.nodesBottomEdge()
top = mesh.nodesTopEdge()

ax.plot(coor[lft,0], coor[lft,1], marker='o', linestyle='none', color='g', label='nodesLeftOpenEdge')
ax.plot(coor[rgt,0], coor[rgt,1], marker='o', linestyle='none', color='b', label='nodesRightOpenEdge')
ax.plot(coor[bot,0], coor[bot,1], marker='o', linestyle='none', color='r', label='nodesBottomEdge')
ax.plot(coor[top,0], coor[top,1], marker='o', linestyle='none', color='y', label='nodesTopEdge')

ax.set_aspect('equal')
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

ax.set_xlim([-1, 55])
ax.set_ylim([-1, 58])

ax.legend(ncol=2, loc='upper center')

plt.savefig('example.svg')
plt.show()

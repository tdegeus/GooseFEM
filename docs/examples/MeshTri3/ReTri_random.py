
import numpy as np
import GooseFEM

# define a mesh, and extract the relevant data
mesh  = GooseFEM.Mesh.Tri3.Regular(10,10)
coor  = mesh.coor()
conn  = mesh.conn()
nnode = coor.shape[0]
nelem = conn.shape[0]
per   = mesh.nodesPeriodic()

# define a random perturbation of the nodal coordinates
delta       = np.random.random(coor.shape)
delta[:,0] *= .1
delta[:,1] *= .1
# make periodic
delta[per[:,1],0] = delta[per[:,0],0]
delta[per[:,1],1] = delta[per[:,0],1]

# copy the coordinates and connectivity (for plotting)
newcoor = np.array(coor,copy=True)
newconn = np.array(conn,copy=True)

# add the perturbation incrementally, compute the new connectivity each time
for i in range(4):
  newcoor += delta/4.
  newconn  = GooseFEM.Mesh.Tri3.retriangulate(newcoor,newconn)

# ==================================================================================================

# define colors to clearly show the difference (and element overlaps)
orig = np.zeros(( nelem ))
diff = np.zeros(( nelem ))
# find the differences
for i in range(nelem):
  if not np.all( conn[i,:] == newconn[i,:] ):
    diff[i] = 1

# define shifts to show periodicity
DX      = np.zeros( coor.shape )
DY      = np.zeros( coor.shape )
DX[:,0] = +1.
DY[:,1] = +1.

# --------------------------------------------------------------------------------------------------

import matplotlib.pylab as plt
import goosempl as gplt

plt.style.use(['goose','goose-latex'])

fig = plt.figure(figsize=(18,6))
fig.set_tight_layout(True)

# --------------------------------------------------------------------------------------------------

ax  = fig.add_subplot(1,3,1)

gplt.patch(coor=coor,conn=conn,cindex=orig,clim=[-0.5,1.5],cmap='RdBu_r',alpha=.5)

gplt.patch(coor=coor-DX   ,conn=conn)
gplt.patch(coor=coor+DX   ,conn=conn)
gplt.patch(coor=coor+DX-DY,conn=conn)
gplt.patch(coor=coor+DX+DY,conn=conn)
gplt.patch(coor=coor-DX-DY,conn=conn)
gplt.patch(coor=coor-DX+DY,conn=conn)
gplt.patch(coor=coor   -DY,conn=conn)
gplt.patch(coor=coor   +DY,conn=conn)

plt.xlim([ -0.15 ,  1.2 ])
plt.ylim([ -0.15 ,  1.2 ])

ax.set_aspect('equal')

# --------------------------------------------------------------------------------------------------

ax  = fig.add_subplot(1,3,2)

gplt.patch(coor=newcoor,conn=conn,cindex=orig,clim=[-0.5,1.5],cmap='RdBu_r',alpha=.5)

gplt.patch(coor=newcoor-DX   ,conn=conn)
gplt.patch(coor=newcoor+DX   ,conn=conn)
gplt.patch(coor=newcoor+DX-DY,conn=conn)
gplt.patch(coor=newcoor+DX+DY,conn=conn)
gplt.patch(coor=newcoor-DX-DY,conn=conn)
gplt.patch(coor=newcoor-DX+DY,conn=conn)
gplt.patch(coor=newcoor   -DY,conn=conn)
gplt.patch(coor=newcoor   +DY,conn=conn)

plt.xlim([ -0.15 ,  1.2 ])
plt.ylim([ -0.15 ,  1.2 ])

ax.set_aspect('equal')

# --------------------------------------------------------------------------------------------------

ax  = fig.add_subplot(1,3,3)

gplt.patch(coor=newcoor,conn=newconn,cindex=diff,clim=[-0.5,1.5],cmap='RdBu_r',alpha=.5)

gplt.patch(coor=newcoor-DX   ,conn=newconn)
gplt.patch(coor=newcoor+DX   ,conn=newconn)
gplt.patch(coor=newcoor+DX-DY,conn=newconn)
gplt.patch(coor=newcoor+DX+DY,conn=newconn)
gplt.patch(coor=newcoor-DX-DY,conn=newconn)
gplt.patch(coor=newcoor-DX+DY,conn=newconn)
gplt.patch(coor=newcoor   -DY,conn=newconn)
gplt.patch(coor=newcoor   +DY,conn=newconn)

plt.xlim([ -0.15 ,  1.2 ])
plt.ylim([ -0.15 ,  1.2 ])

ax.set_aspect('equal')

# --------------------------------------------------------------------------------------------------

plt.savefig('ReTri_random.svg')
plt.show()

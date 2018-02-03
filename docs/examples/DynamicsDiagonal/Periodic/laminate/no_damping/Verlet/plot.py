
import sys,os,re
import h5py
import numpy as np
import matplotlib.pyplot as plt
import goosempl as gplt
import matplotlib as mpl

plt.style.use(['goose','goose-latex'])

# --------------------------------------------------------------------------------------------------

name = 'example.hdf5'
f    = h5py.File(name,'r')

fig = plt.figure(figsize=(20,10))
fig.set_tight_layout(True)

# --------------------------------------------------------------------------------------------------

ax = fig.add_subplot(1,2,1)

t = f['/global/t'   ][:]
E = f['/global/Ekin'][:] + f['/global/Epot'][:]; ax.plot(t,E,color='k',label=r'$E_\mathrm{pot} + E_\mathrm{kin}$')
E = f['/global/Epot'][:]                       ; ax.plot(t,E,color='r',label=r'$E_\mathrm{pot}$')
E = f['/global/Ekin'][:]                       ; ax.plot(t,E,color='b',label=r'$E_\mathrm{kin}$')

plt.legend(loc='lower right')

plt.xlabel(r'$t$')
plt.ylabel(r'$E$')

# --------------------------------------------------------------------------------------------------

ax = fig.add_subplot(1,2,2)

coor = f['/mesh/coor'][:]
conn = f['/mesh/conn'][:]
u    = f['/mesh/disp'][:]

im = gplt.patch(coor=coor+u,conn=conn)

Lx = np.max(coor[:,0])-np.min(coor[:,0])
Ly = np.max(coor[:,1])-np.min(coor[:,1])

plt.xlim([-0.1*Lx,1.3*Lx])
plt.ylim([-0.1*Ly,1.3*Ly])

ax.set_aspect('equal')

plt.savefig(re.sub(r'(.*)(\.hdf5)',r'\1.svg',name))
plt.show()

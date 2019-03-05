
import h5py

import matplotlib.pyplot as plt
import GooseMPL          as gplt
import numpy             as np

plt.style.use(['goose', 'goose-latex'])

# open file
file = h5py.File('main.h5', 'r')

# read fields
coor = file['/coor'][...]
conn = file['/conn'][...]
disp = file['/disp'][...]
Sig  = file['/Sig' ][...]

# extract dimension
nelem = conn.shape[0]

# tensor products
ddot22 = lambda A2,B2: np.einsum('eij  ,eji->e    ',A2,B2)
ddot42 = lambda A4,B2: np.einsum('eijkl,elk->eij  ',A4,B2)
dyad22 = lambda A2,B2: np.einsum('eij  ,ekl->eijkl',A2,B2)

# identity tensor (single tensor)
i    = np.eye(3)

# identity tensors (grid)
I    = np.einsum('ij  ,e'       ,                  i   ,np.ones([nelem]))
I4   = np.einsum('ijkl,e->eijkl',np.einsum('il,jk',i,i),np.ones([nelem]))
I4rt = np.einsum('ijkl,e->eijkl',np.einsum('ik,jl',i,i),np.ones([nelem]))
I4s  = (I4+I4rt)/2.
II   = dyad22(I,I)
I4d  = I4s-II/3.

# compute equivalent stress
Sigd  = ddot42(I4d, Sig)
sigeq = np.sqrt(3./2.*ddot22(Sigd,Sigd))

# plot

fig, ax = plt.subplots()

gplt.patch(coor=coor+disp, conn=conn, cindex=sigeq, cmap='jet', axis=ax, clim=(0,0.5))

gplt.patch(coor=coor, conn=conn, linestyle='--', axis=ax)

plt.show()


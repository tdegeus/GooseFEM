import h5py
import matplotlib.pyplot as plt
import GooseMPL as gplt
import numpy as np

plt.style.use(['goose', 'goose-latex'])

# load data

with h5py.File('output.h5', 'r') as data:
    coor = data['/coor'][...]
    conn = data['/conn'][...]
    disp = data['/disp'][...]
    Sig = data['/Sig'][...]
    nelem = conn.shape[0]

# tensor products

def ddot22(A2, B2):
    return np.einsum('eij, eji -> e', A2, B2)

def ddot42(A4, B2):
    return np.einsum('eijkl, elk -> eij', A4, B2)

def dyad22(A2, B2):
    return np.einsum('eij, ekl -> eijkl', A2, B2)

# identity tensors

i = np.eye(3)
I = np.einsum('ij, e', i, np.ones([nelem]))
I4 = np.einsum('ijkl, e -> eijkl', np.einsum('il, jk', i, i), np.ones([nelem]))
I4rt = np.einsum('ijkl, e -> eijkl', np.einsum('ik,jl', i, i), np.ones([nelem]))
I4s = 0.5 * (I4 + I4rt)
II = dyad22(I, I)
I4d = I4s - II / 3.0

# compute equivalent stress

Sigd = ddot42(I4d, Sig)
sigeq = np.sqrt(3.0 / 2.0 * ddot22(Sigd, Sigd))

# plot

fig, ax = plt.subplots()
gplt.patch(coor=coor + disp, conn=conn, cindex=sigeq, cmap='jet', axis=ax, clim=(0, 0.1))
gplt.patch(coor=coor, conn=conn, linestyle='--', axis=ax)
plt.savefig('plot.pdf')
plt.close()

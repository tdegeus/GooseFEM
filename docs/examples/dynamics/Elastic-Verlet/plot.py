
import h5py
import numpy as np
import matplotlib.pyplot as plt
import GooseMPL as gplt

plt.style.use(['goose', 'goose-latex'])

# --------------------------------------------------------------------------------------------------

with h5py.File('example.hdf5', 'r') as data:

    t = data['/global/t'][...]

    Epot = data['/global/Epot'][...]
    Ekin = data['/global/Ekin'][...]

    coor = data['/mesh/coor'][...]
    conn = data['/mesh/conn'][...]
    disp = data['/mesh/disp'][...]

# --------------------------------------------------------------------------------------------------

fig, axes = gplt.subplots(ncols=2)

ax = axes[0]

ax.plot(t, Ekin + Epot, color='k', label=r'$E_\mathrm{pot} + E_\mathrm{kin}$')
ax.plot(t, Epot, color='r', label=r'$E_\mathrm{pot}$')
ax.plot(t, Ekin, color='b', label=r'$E_\mathrm{kin}$')

ax.legend(loc='lower right')

ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$E$')

ax = axes[1]

im = gplt.patch(coor=coor + disp, conn=conn)

Lx = np.max(coor[:, 0]) - np.min(coor[:, 0])
Ly = np.max(coor[:, 1]) - np.min(coor[:, 1])

ax.set_xlim([-0.1 * Lx, 1.3 * Lx])
ax.set_ylim([-0.1 * Ly, 1.3 * Ly])

ax.set_aspect('equal')

plt.show()

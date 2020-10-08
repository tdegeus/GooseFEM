
import h5py

import matplotlib.pyplot as plt
import GooseMPL as gplt
import numpy as np

plt.style.use(['goose', 'goose-latex'])

# open file
file = h5py.File('main.h5', 'r')

# read stored increments
incs = file['/stored'][...]

# read fields
coor = file['/coor'][...]
conn = file['/conn'][...]
disp = file['/disp'][str(np.max(incs))][...][:, :2]
Sigeq = file['/sigeq'][str(np.max(incs))][...]
sigeq = file['/macroscopic/sigeq'][...]
epseq = file['/macroscopic/epseq'][...]

# plot

fig, axes = gplt.subplots(ncols=2)

axes[0].plot(epseq, sigeq)

gplt.patch(coor=coor + disp, conn=conn, cindex=Sigeq, cmap='jet', axis=axes[1], clim=(0.0, 0.1))

plt.show()

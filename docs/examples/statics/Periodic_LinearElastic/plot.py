
import h5py

import matplotlib.pyplot as plt
import GooseMPL as gplt
import numpy as np

plt.style.use(['goose', 'goose-latex'])

# open file
file = h5py.File('main.h5', 'r')

# read fields
coor = file['/coor'][...]
conn = file['/conn'][...]
disp = file['/disp'][...][:, :2]
Sigeq = file['/sigeq'][...]

# plot

fig, ax = plt.subplots()

gplt.patch(coor=coor + disp, conn=conn, cindex=Sigeq, cmap='jet', axis=ax, clim=(0, 0.5))

plt.show()

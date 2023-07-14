import os

import GooseFEM as gf
import GooseFEM.ParaView.HDF5 as pv
import h5py

# create mesh object
mesh = gf.Mesh.Quad4.Regular(9, 11)

# filename of the HDF5-file
fname = "paraview.hdf5"

# write HDF-file containing the data
with h5py.File(fname, "w") as data:
    data["coor"] = mesh.coor
    data["conn"] = mesh.conn

# write XDMF-file with metadata
pv.Mesh(
    pv.Connectivity(fname, "/conn", mesh.getElementType(), mesh.conn.shape),
    pv.Coordinates(fname, "/coor", mesh.coor.shape),
).write(os.path.splitext(fname)[0] + ".xdmf")

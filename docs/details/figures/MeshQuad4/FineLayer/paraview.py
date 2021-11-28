import os

import GooseFEM as gf
import GooseFEM.ParaView.HDF5 as pv
import h5py
import numpy as np

# create mesh object
mesh = gf.Mesh.Quad4.FineLayer(9, 17)

# filename of the HDF5-file
fname = "paraview.hdf5"

# define element set
elementsMiddleLayer = np.zeros((mesh.nelem()), dtype="int")
elementsMiddleLayer[mesh.elementsMiddleLayer()] = 1

# write HDF-file containing the data
with h5py.File(fname, "w") as data:
    data["coor"] = mesh.coor()
    data["conn"] = mesh.conn()
    data["elementsMiddleLayer"] = elementsMiddleLayer

# write XDMF-file with metadata
xdmf = pv.Mesh(
    pv.Connectivity(fname, "/conn", mesh.getElementType(), mesh.conn().shape),
    pv.Coordinates(fname, "/coor", mesh.coor().shape),
)

xdmf.push_back(
    pv.Attribute(
        fname,
        "/elementsMiddleLayer",
        "elementsMiddleLayer",
        pv.AttributeType.Cell,
        elementsMiddleLayer.shape,
    )
)

xdmf.write(os.path.splitext(fname)[0] + ".xdmf")

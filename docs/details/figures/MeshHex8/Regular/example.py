
import h5py, os
import numpy as np
import GooseFEM as gf
import GooseFEM.ParaView.HDF5 as pv

# create mesh object
mesh = gf.Mesh.Hex8.Regular(6,8,12)

# filename of the HDF5-file
fname = 'MeshHex8-Regular.hdf5'

# write HDF-file containing the data

with h5py.File(fname, 'w') as data:

  data['coor'] = mesh.coor()
  data['conn'] = mesh.conn()

# write XDMF-file with metadata

pv.Mesh(
  pv.Connectivity(fname, "/conn", pv.ElementType.Hexahedron, mesh.conn().shape),
  pv.Coordinates(fname, "/coor", mesh.coor().shape),
).write(os.path.splitext(fname)[0] + '.xdmf')

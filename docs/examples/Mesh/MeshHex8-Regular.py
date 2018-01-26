
import h5py
import numpy      as np
import GooseFEM   as gf
import lxml.etree as etree

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Hex8.Regular(6,8,12)

# open data file
name = 'MeshHex8-Regular'
file = h5py.File(name+'.hdf5','w')

# write nodal coordinates and connectivity
file['coor'] = mesh.coor()
file['conn'] = mesh.conn()

# ======================================== write XDMF-file =========================================

# get mesh dimensions
dims = dict(
  nnode = mesh.nnode(),
  ndim  = mesh.ndim(),
  nelem = mesh.nelem(),
  nne   = mesh.nne(),
)

# initialize file
root   = etree.fromstring('<Xdmf Version="2.0"></Xdmf>')
domain = etree.SubElement(root, "Domain")
grid   = etree.SubElement(domain, "Grid", Name="Mesh")

# add connectivity
conn = etree.SubElement(grid, "Topology", TopologyType="Hexahedron", NumberOfElements='{nelem:d}'.format(**dims))
data = etree.SubElement(conn, "DataItem", Dimensions='{nelem:d} {nne:d}'.format(**dims), Format="HDF")
data.text = "{0:s}.hdf5:/conn".format(name)

# add coordinates
coor = etree.SubElement(grid, "Geometry", GeometryType="XYZ")
data = etree.SubElement(coor, "DataItem", Dimensions='{nnode:d} {ndim:d}'.format(**dims), Format="HDF")
data.text = "{0:s}.hdf5:/coor".format(name)

# write to file
open(name+'.xdmf','wb').write(etree.tostring(root, pretty_print=True))


import h5py
import numpy      as np
import GooseFEM   as gf
import lxml.etree as etree

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Hex8.FineLayer(9,17,27)

# open data file
name = 'MeshHex8-FineLayer'
file = h5py.File(name+'.hdf5','w')

# element set
elementsMiddleLayer = np.zeros((mesh.nelem()),dtype='int')
elementsMiddleLayer[mesh.elementsMiddleLayer()] = 1

# write nodal coordinates and connectivity
file['coor'] = mesh.coor()
file['conn'] = mesh.conn()
file['elementsMiddleLayer'] = elementsMiddleLayer

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

# add attributes
attr = etree.SubElement(grid, "Attribute", Name="elementsMiddleLayer", AttributeType="Scalar", Center="Cell")
data = etree.SubElement(attr, "DataItem", Dimensions='{nelem:d}'.format(**dims), Format="HDF")
data.text = "{0:s}.hdf5:/elementsMiddleLayer".format(name)

# write to file
open(name+'.xdmf','wb').write(etree.tostring(root, pretty_print=True))

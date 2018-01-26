
import h5py
import numpy      as np
import GooseFEM   as gf
import lxml.etree as etree

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Quad4.Regular(9,11)

# initialize all node sets
nodesets = dict(
  nodesBottomEdge        = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopEdge           = np.zeros((mesh.nnode()),dtype='int'),
  nodesLeftEdge          = np.zeros((mesh.nnode()),dtype='int'),
  nodesRightEdge         = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottomOpenEdge    = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopOpenEdge       = np.zeros((mesh.nnode()),dtype='int'),
  nodesLeftOpenEdge      = np.zeros((mesh.nnode()),dtype='int'),
  nodesRightOpenEdge     = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottomLeftCorner  = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottomRightCorner = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopLeftCorner     = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopRightCorner    = np.zeros((mesh.nnode()),dtype='int'),
)

# define node-sets
nodesets['nodesBottomEdge'       ][mesh.nodesBottomEdge()       ] = 1
nodesets['nodesTopEdge'          ][mesh.nodesTopEdge()          ] = 1
nodesets['nodesLeftEdge'         ][mesh.nodesLeftEdge()         ] = 1
nodesets['nodesRightEdge'        ][mesh.nodesRightEdge()        ] = 1
nodesets['nodesBottomOpenEdge'   ][mesh.nodesBottomOpenEdge()   ] = 1
nodesets['nodesTopOpenEdge'      ][mesh.nodesTopOpenEdge()      ] = 1
nodesets['nodesLeftOpenEdge'     ][mesh.nodesLeftOpenEdge()     ] = 1
nodesets['nodesRightOpenEdge'    ][mesh.nodesRightOpenEdge()    ] = 1
nodesets['nodesBottomLeftCorner' ][mesh.nodesBottomLeftCorner() ] = 1
nodesets['nodesBottomRightCorner'][mesh.nodesBottomRightCorner()] = 1
nodesets['nodesTopLeftCorner'    ][mesh.nodesTopLeftCorner()    ] = 1
nodesets['nodesTopRightCorner'   ][mesh.nodesTopRightCorner()   ] = 1

# add DOF-numbers after eliminating periodicity
nodesets['dofsPeriodic'] = mesh.dofsPeriodic()[:,0]

# open data file
name = 'MeshQuad4-Regular-nodes'
file = h5py.File(name+'.hdf5','w')

# write nodal positions and a dummy connectivity
file['coor'] = mesh.coor()
file['conn'] = np.arange(mesh.nnode())

# write node-sets
for key,value in nodesets.items():
  file[key] = value

# ======================================== write XDMF-file =========================================

# get mesh dimensions
dims = dict(
  nnode = mesh.nnode(),
  ndim  = mesh.ndim(),
)

# initialize file
root   = etree.fromstring('<Xdmf Version="2.0"></Xdmf>')
domain = etree.SubElement(root, "Domain")
grid   = etree.SubElement(domain, "Grid", Name="Nodes")

# add connectivity
conn = etree.SubElement(grid, "Topology", TopologyType="Polyvertex", NumberOfElements='{nnode:d}'.format(**dims), NodesPerElement="1")
data = etree.SubElement(conn, "DataItem", Dimensions='{nnode:d} 1'.format(**dims), Format="HDF")
data.text = "{0:s}.hdf5:/conn".format(name)

# add coordinates
coor = etree.SubElement(grid, "Geometry", GeometryType="XY")
data = etree.SubElement(coor, "DataItem", Dimensions='{nnode:d} {ndim:d}'.format(**dims), Format="HDF")
data.text = "{0:s}.hdf5:/coor".format(name)

# add attributes
for key in nodesets:
  attr = etree.SubElement(grid, "Attribute", Name=key, AttributeType="Scalar", Center="Node")
  data = etree.SubElement(attr, "DataItem", Dimensions='{nnode:d}'.format(**dims), Format="HDF")
  data.text = "{name:s}.hdf5:/{key:s}".format(name=name,key=key)

# write to file
open(name+'.xdmf','wb').write(etree.tostring(root, pretty_print=True))

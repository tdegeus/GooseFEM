
import h5py
import numpy      as np
import GooseFEM   as gf
import lxml.etree as etree

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Hex8.FineLayer(9,17,27)

# initialize all node sets
nodesets = dict(
  nodesFront                  = np.zeros((mesh.nnode()),dtype='int'),
  nodesBack                   = np.zeros((mesh.nnode()),dtype='int'),
  nodesLeft                   = np.zeros((mesh.nnode()),dtype='int'),
  nodesRight                  = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottom                 = np.zeros((mesh.nnode()),dtype='int'),
  nodesTop                    = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontFace              = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackFace               = np.zeros((mesh.nnode()),dtype='int'),
  nodesLeftFace               = np.zeros((mesh.nnode()),dtype='int'),
  nodesRightFace              = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottomFace             = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopFace                = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontBottomEdge        = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontTopEdge           = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontLeftEdge          = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontRightEdge         = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackBottomEdge         = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackTopEdge            = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackLeftEdge           = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackRightEdge          = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottomLeftEdge         = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottomRightEdge        = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopLeftEdge            = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopRightEdge           = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontBottomOpenEdge    = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontTopOpenEdge       = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontLeftOpenEdge      = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontRightOpenEdge     = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackBottomOpenEdge     = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackTopOpenEdge        = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackLeftOpenEdge       = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackRightOpenEdge      = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottomLeftOpenEdge     = np.zeros((mesh.nnode()),dtype='int'),
  nodesBottomRightOpenEdge    = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopLeftOpenEdge        = np.zeros((mesh.nnode()),dtype='int'),
  nodesTopRightOpenEdge       = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontBottomLeftCorner  = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontBottomRightCorner = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontTopLeftCorner     = np.zeros((mesh.nnode()),dtype='int'),
  nodesFrontTopRightCorner    = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackBottomLeftCorner   = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackBottomRightCorner  = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackTopLeftCorner      = np.zeros((mesh.nnode()),dtype='int'),
  nodesBackTopRightCorner     = np.zeros((mesh.nnode()),dtype='int'),
)

# define node-sets
nodesets['nodesFront'                 ][mesh.nodesFront()                 ] = 1
nodesets['nodesBack'                  ][mesh.nodesBack()                  ] = 1
nodesets['nodesLeft'                  ][mesh.nodesLeft()                  ] = 1
nodesets['nodesRight'                 ][mesh.nodesRight()                 ] = 1
nodesets['nodesBottom'                ][mesh.nodesBottom()                ] = 1
nodesets['nodesTop'                   ][mesh.nodesTop()                   ] = 1
nodesets['nodesFrontFace'             ][mesh.nodesFrontFace()             ] = 1
nodesets['nodesBackFace'              ][mesh.nodesBackFace()              ] = 1
nodesets['nodesLeftFace'              ][mesh.nodesLeftFace()              ] = 1
nodesets['nodesRightFace'             ][mesh.nodesRightFace()             ] = 1
nodesets['nodesBottomFace'            ][mesh.nodesBottomFace()            ] = 1
nodesets['nodesTopFace'               ][mesh.nodesTopFace()               ] = 1
nodesets['nodesFrontBottomEdge'       ][mesh.nodesFrontBottomEdge()       ] = 1
nodesets['nodesFrontTopEdge'          ][mesh.nodesFrontTopEdge()          ] = 1
nodesets['nodesFrontLeftEdge'         ][mesh.nodesFrontLeftEdge()         ] = 1
nodesets['nodesFrontRightEdge'        ][mesh.nodesFrontRightEdge()        ] = 1
nodesets['nodesBackBottomEdge'        ][mesh.nodesBackBottomEdge()        ] = 1
nodesets['nodesBackTopEdge'           ][mesh.nodesBackTopEdge()           ] = 1
nodesets['nodesBackLeftEdge'          ][mesh.nodesBackLeftEdge()          ] = 1
nodesets['nodesBackRightEdge'         ][mesh.nodesBackRightEdge()         ] = 1
nodesets['nodesBottomLeftEdge'        ][mesh.nodesBottomLeftEdge()        ] = 1
nodesets['nodesBottomRightEdge'       ][mesh.nodesBottomRightEdge()       ] = 1
nodesets['nodesTopLeftEdge'           ][mesh.nodesTopLeftEdge()           ] = 1
nodesets['nodesTopRightEdge'          ][mesh.nodesTopRightEdge()          ] = 1
nodesets['nodesFrontBottomOpenEdge'   ][mesh.nodesFrontBottomOpenEdge()   ] = 1
nodesets['nodesFrontTopOpenEdge'      ][mesh.nodesFrontTopOpenEdge()      ] = 1
nodesets['nodesFrontLeftOpenEdge'     ][mesh.nodesFrontLeftOpenEdge()     ] = 1
nodesets['nodesFrontRightOpenEdge'    ][mesh.nodesFrontRightOpenEdge()    ] = 1
nodesets['nodesBackBottomOpenEdge'    ][mesh.nodesBackBottomOpenEdge()    ] = 1
nodesets['nodesBackTopOpenEdge'       ][mesh.nodesBackTopOpenEdge()       ] = 1
nodesets['nodesBackLeftOpenEdge'      ][mesh.nodesBackLeftOpenEdge()      ] = 1
nodesets['nodesBackRightOpenEdge'     ][mesh.nodesBackRightOpenEdge()     ] = 1
nodesets['nodesBottomLeftOpenEdge'    ][mesh.nodesBottomLeftOpenEdge()    ] = 1
nodesets['nodesBottomRightOpenEdge'   ][mesh.nodesBottomRightOpenEdge()   ] = 1
nodesets['nodesTopLeftOpenEdge'       ][mesh.nodesTopLeftOpenEdge()       ] = 1
nodesets['nodesTopRightOpenEdge'      ][mesh.nodesTopRightOpenEdge()      ] = 1
nodesets['nodesFrontBottomLeftCorner' ][mesh.nodesFrontBottomLeftCorner() ] = 1
nodesets['nodesFrontBottomRightCorner'][mesh.nodesFrontBottomRightCorner()] = 1
nodesets['nodesFrontTopLeftCorner'    ][mesh.nodesFrontTopLeftCorner()    ] = 1
nodesets['nodesFrontTopRightCorner'   ][mesh.nodesFrontTopRightCorner()   ] = 1
nodesets['nodesBackBottomLeftCorner'  ][mesh.nodesBackBottomLeftCorner()  ] = 1
nodesets['nodesBackBottomRightCorner' ][mesh.nodesBackBottomRightCorner() ] = 1
nodesets['nodesBackTopLeftCorner'     ][mesh.nodesBackTopLeftCorner()     ] = 1
nodesets['nodesBackTopRightCorner'    ][mesh.nodesBackTopRightCorner()    ] = 1

# add DOF-numbers after eliminating periodicity
nodesets['dofsPeriodic'] = mesh.dofsPeriodic()[:,0]

# open data file
name = 'MeshHex8-FineLayer-nodes'
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
coor = etree.SubElement(grid, "Geometry", GeometryType="XYZ")
data = etree.SubElement(coor, "DataItem", Dimensions='{nnode:d} {ndim:d}'.format(**dims), Format="HDF")
data.text = "{0:s}.hdf5:/coor".format(name)

# add attributes
for key in nodesets:
  attr = etree.SubElement(grid, "Attribute", Name=key, AttributeType="Scalar", Center="Node")
  data = etree.SubElement(attr, "DataItem", Dimensions='{nnode:d}'.format(**dims), Format="HDF")
  data.text = "{name:s}.hdf5:/{key:s}".format(name=name,key=key)

# write to file
open(name+'.xdmf','wb').write(etree.tostring(root, pretty_print=True))

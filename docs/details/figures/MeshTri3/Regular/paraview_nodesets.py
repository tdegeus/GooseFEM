
import h5py, os
import numpy as np
import GooseFEM as gf
import GooseFEM.ParaView.HDF5 as pv

# create mesh object
mesh = gf.Mesh.Tri3.Regular(9,11)

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

# filename of the HDF5-file
fname = 'paraview_nodesets.hdf5'

# write HDF-file containing the data
with h5py.File(fname, 'w') as data:
  data.file['coor'] = mesh.coor()
  data.file['conn'] = mesh.conn()

  for key, value in nodesets.items():
    data[key] = value

# write XDMF-file with metadata
xdmf = pv.Mesh(
  pv.Connectivity(fname, "/conn", mesh.getElementType(), mesh.conn().shape),
  pv.Coordinates(fname, "/coor", mesh.coor().shape),
)

for key, value in nodesets.items():
  xdmf.push_back(pv.Attribute(fname, "/"+key, key, pv.AttributeType.Node, value.shape))

xdmf.write(os.path.splitext(fname)[0] + '.xdmf')

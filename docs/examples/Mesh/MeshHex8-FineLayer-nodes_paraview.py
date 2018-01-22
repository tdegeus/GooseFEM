
import numpy as np
import h5py
import GooseFEM as gf

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Hex8.FineLayer(9,17,27)

# create node type, based on boundary sets
# - allocate all as background
nodes = np.ones((mesh.nnode()),dtype='float64')
# - faces
nodes[mesh.nodesFrontFace()             ] = 2.0
nodes[mesh.nodesBackFace()              ] = 2.0
nodes[mesh.nodesLeftFace()              ] = 3.0
nodes[mesh.nodesRightFace()             ] = 3.0
nodes[mesh.nodesBottomFace()            ] = 4.0
nodes[mesh.nodesTopFace()               ] = 4.0
# - edges (without corners)
nodes[mesh.nodesFrontBottomOpenEdge()   ] = 2.5
nodes[mesh.nodesFrontTopOpenEdge()      ] = 2.5
nodes[mesh.nodesFrontLeftOpenEdge()     ] = 2.5
nodes[mesh.nodesFrontRightOpenEdge()    ] = 2.5
nodes[mesh.nodesBackBottomOpenEdge()    ] = 2.5
nodes[mesh.nodesBackTopOpenEdge()       ] = 2.5
nodes[mesh.nodesBackLeftOpenEdge()      ] = 2.5
nodes[mesh.nodesBackRightOpenEdge()     ] = 2.5
nodes[mesh.nodesBottomLeftOpenEdge()    ] = 4.5
nodes[mesh.nodesBottomRightOpenEdge()   ] = 4.5
nodes[mesh.nodesTopLeftOpenEdge()       ] = 4.5
nodes[mesh.nodesTopRightOpenEdge()      ] = 4.5
# - corners
nodes[mesh.nodesFrontBottomLeftCorner() ] = 5.0
nodes[mesh.nodesFrontBottomRightCorner()] = 5.0
nodes[mesh.nodesFrontTopLeftCorner()    ] = 5.0
nodes[mesh.nodesFrontTopRightCorner()   ] = 5.0
nodes[mesh.nodesBackBottomLeftCorner()  ] = 5.0
nodes[mesh.nodesBackBottomRightCorner() ] = 5.0
nodes[mesh.nodesBackTopLeftCorner()     ] = 5.0
nodes[mesh.nodesBackTopRightCorner()    ] = 5.0

# DOFs after eliminating periodicity
dofsPeriodic = mesh.dofsPeriodic()[:,0]

# open data file
f = h5py.File('MeshHex8-FineLayer-nodes_paraview.hdf5','w')

# write particle positions, and a dummy connectivity
f.create_dataset('/coor'        ,data=mesh.coor()            )
f.create_dataset('/nodes'       ,data=nodes                  )
f.create_dataset('/dofsPeriodic',data=dofsPeriodic           )
f.create_dataset('/conn'        ,data=np.arange(mesh.nnode()))

# ======================================== write XDMF-file =========================================

xmf = '''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="Nodes">
      <Topology TopologyType="Polyvertex" NumberOfElements="{nnode:d}" NodesPerElement="1">
        <DataItem Dimensions="{nnode:d} 1" Format="HDF">
          MeshHex8-FineLayer-nodes_paraview.hdf5:/conn
        </DataItem>
      </Topology>
      <Geometry GeometryType="XYZ">
        <DataItem Dimensions="{nnode:d} 3" Format="HDF">
          MeshHex8-FineLayer-nodes_paraview.hdf5:/coor
        </DataItem>
      </Geometry>
      <Attribute Name="Node type" AttributeType="Scalar" Center="Node">
         <DataItem Dimensions="{nnode:d}" NumberType="Float" Precision="8" Format="HDF">
          MeshHex8-FineLayer-nodes_paraview.hdf5:/nodes
         </DataItem>
      </Attribute>
      <Attribute Name="dofsPeriodic" AttributeType="Scalar" Center="Node">
         <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-FineLayer-nodes_paraview.hdf5:/dofsPeriodic
         </DataItem>
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
'''

open('MeshHex8-FineLayer-nodes_paraview.xmf','w').write(xmf.format(
  nnode = mesh.nnode(),
))


import numpy as np
import h5py
import GooseFEM as gf

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Quad4.FineLayer(9,17)

# set node identifier for all boundary sets
Z                      = lambda mesh: np.zeros((mesh.nnode()),dtype='int')
nodesBottomEdge        = Z(mesh); nodesBottomEdge       [mesh.nodesBottomEdge()       ] = 1
nodesTopEdge           = Z(mesh); nodesTopEdge          [mesh.nodesTopEdge()          ] = 1
nodesLeftEdge          = Z(mesh); nodesLeftEdge         [mesh.nodesLeftEdge()         ] = 1
nodesRightEdge         = Z(mesh); nodesRightEdge        [mesh.nodesRightEdge()        ] = 1
nodesBottomOpenEdge    = Z(mesh); nodesBottomOpenEdge   [mesh.nodesBottomOpenEdge()   ] = 1
nodesTopOpenEdge       = Z(mesh); nodesTopOpenEdge      [mesh.nodesTopOpenEdge()      ] = 1
nodesLeftOpenEdge      = Z(mesh); nodesLeftOpenEdge     [mesh.nodesLeftOpenEdge()     ] = 1
nodesRightOpenEdge     = Z(mesh); nodesRightOpenEdge    [mesh.nodesRightOpenEdge()    ] = 1
nodesBottomLeftCorner  = Z(mesh); nodesBottomLeftCorner [mesh.nodesBottomLeftCorner() ] = 1
nodesBottomRightCorner = Z(mesh); nodesBottomRightCorner[mesh.nodesBottomRightCorner()] = 1
nodesTopLeftCorner     = Z(mesh); nodesTopLeftCorner    [mesh.nodesTopLeftCorner()    ] = 1
nodesTopRightCorner    = Z(mesh); nodesTopRightCorner   [mesh.nodesTopRightCorner()   ] = 1

# DOFs after eliminating periodicity
dofsPeriodic = mesh.dofsPeriodic()[:,0]

# open data file
f = h5py.File('MeshQuad4-FineLayer-nodes_paraview.hdf5','w')

# write nodal positions, and a dummy connectivity
f.create_dataset('/coor'                  ,data=mesh.coor()            )
f.create_dataset('/conn'                  ,data=np.arange(mesh.nnode()))
f.create_dataset('/nodesBottomEdge'       ,data=nodesBottomEdge        )
f.create_dataset('/nodesTopEdge'          ,data=nodesTopEdge           )
f.create_dataset('/nodesLeftEdge'         ,data=nodesLeftEdge          )
f.create_dataset('/nodesRightEdge'        ,data=nodesRightEdge         )
f.create_dataset('/nodesBottomOpenEdge'   ,data=nodesBottomOpenEdge    )
f.create_dataset('/nodesTopOpenEdge'      ,data=nodesTopOpenEdge       )
f.create_dataset('/nodesLeftOpenEdge'     ,data=nodesLeftOpenEdge      )
f.create_dataset('/nodesRightOpenEdge'    ,data=nodesRightOpenEdge     )
f.create_dataset('/nodesBottomLeftCorner' ,data=nodesBottomLeftCorner  )
f.create_dataset('/nodesBottomRightCorner',data=nodesBottomRightCorner )
f.create_dataset('/nodesTopLeftCorner'    ,data=nodesTopLeftCorner     )
f.create_dataset('/nodesTopRightCorner'   ,data=nodesTopRightCorner    )
f.create_dataset('/dofsPeriodic'          ,data=dofsPeriodic           )

# ======================================== write XDMF-file =========================================

# basic file-format
xdmf = '''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="Nodes">
      <Topology TopologyType="Polyvertex" NumberOfElements="{nnode:d}" NodesPerElement="1">
        <DataItem Dimensions="{nnode:d} 1" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/conn
        </DataItem>
      </Topology>
      <Geometry GeometryType="XY">
        <DataItem Dimensions="{nnode:d} 2" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/coor
        </DataItem>
      </Geometry>
      <Attribute Name="nodesBottomEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesBottomEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesTopEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesLeftEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesLeftEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesRightEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesRightEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottomOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesBottomOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesTopOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesLeftOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesLeftOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesRightOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesRightOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottomLeftCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesBottomLeftCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottomRightCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesBottomRightCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopLeftCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesTopLeftCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopRightCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/nodesTopRightCorner
        </DataItem>
      </Attribute>
      <Attribute Name="dofsPeriodic" AttributeType="Scalar" Center="Node">
         <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer-nodes_paraview.hdf5:/dofsPeriodic
         </DataItem>
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
'''

# write to file, fill mesh dimensions
open('MeshQuad4-FineLayer-nodes_paraview.xdmf','w').write(xdmf.format(
  nnode = mesh.nnode(),
))

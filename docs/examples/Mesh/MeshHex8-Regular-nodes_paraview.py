
import numpy as np
import h5py
import GooseFEM as gf

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Hex8.Regular(6,8,12)

# set node identifier for all boundary sets
Z                           = lambda mesh: np.zeros((mesh.nnode()),dtype='int')
nodesFront                  = Z(mesh); nodesFront                 [mesh.nodesFront()                 ] = 1
nodesBack                   = Z(mesh); nodesBack                  [mesh.nodesBack()                  ] = 1
nodesLeft                   = Z(mesh); nodesLeft                  [mesh.nodesLeft()                  ] = 1
nodesRight                  = Z(mesh); nodesRight                 [mesh.nodesRight()                 ] = 1
nodesBottom                 = Z(mesh); nodesBottom                [mesh.nodesBottom()                ] = 1
nodesTop                    = Z(mesh); nodesTop                   [mesh.nodesTop()                   ] = 1
nodesFrontFace              = Z(mesh); nodesFrontFace             [mesh.nodesFrontFace()             ] = 1
nodesBackFace               = Z(mesh); nodesBackFace              [mesh.nodesBackFace()              ] = 1
nodesLeftFace               = Z(mesh); nodesLeftFace              [mesh.nodesLeftFace()              ] = 1
nodesRightFace              = Z(mesh); nodesRightFace             [mesh.nodesRightFace()             ] = 1
nodesBottomFace             = Z(mesh); nodesBottomFace            [mesh.nodesBottomFace()            ] = 1
nodesTopFace                = Z(mesh); nodesTopFace               [mesh.nodesTopFace()               ] = 1
nodesFrontBottomEdge        = Z(mesh); nodesFrontBottomEdge       [mesh.nodesFrontBottomEdge()       ] = 1
nodesFrontTopEdge           = Z(mesh); nodesFrontTopEdge          [mesh.nodesFrontTopEdge()          ] = 1
nodesFrontLeftEdge          = Z(mesh); nodesFrontLeftEdge         [mesh.nodesFrontLeftEdge()         ] = 1
nodesFrontRightEdge         = Z(mesh); nodesFrontRightEdge        [mesh.nodesFrontRightEdge()        ] = 1
nodesBackBottomEdge         = Z(mesh); nodesBackBottomEdge        [mesh.nodesBackBottomEdge()        ] = 1
nodesBackTopEdge            = Z(mesh); nodesBackTopEdge           [mesh.nodesBackTopEdge()           ] = 1
nodesBackLeftEdge           = Z(mesh); nodesBackLeftEdge          [mesh.nodesBackLeftEdge()          ] = 1
nodesBackRightEdge          = Z(mesh); nodesBackRightEdge         [mesh.nodesBackRightEdge()         ] = 1
nodesBottomLeftEdge         = Z(mesh); nodesBottomLeftEdge        [mesh.nodesBottomLeftEdge()        ] = 1
nodesBottomRightEdge        = Z(mesh); nodesBottomRightEdge       [mesh.nodesBottomRightEdge()       ] = 1
nodesTopLeftEdge            = Z(mesh); nodesTopLeftEdge           [mesh.nodesTopLeftEdge()           ] = 1
nodesTopRightEdge           = Z(mesh); nodesTopRightEdge          [mesh.nodesTopRightEdge()          ] = 1
nodesFrontBottomOpenEdge    = Z(mesh); nodesFrontBottomOpenEdge   [mesh.nodesFrontBottomOpenEdge()   ] = 1
nodesFrontTopOpenEdge       = Z(mesh); nodesFrontTopOpenEdge      [mesh.nodesFrontTopOpenEdge()      ] = 1
nodesFrontLeftOpenEdge      = Z(mesh); nodesFrontLeftOpenEdge     [mesh.nodesFrontLeftOpenEdge()     ] = 1
nodesFrontRightOpenEdge     = Z(mesh); nodesFrontRightOpenEdge    [mesh.nodesFrontRightOpenEdge()    ] = 1
nodesBackBottomOpenEdge     = Z(mesh); nodesBackBottomOpenEdge    [mesh.nodesBackBottomOpenEdge()    ] = 1
nodesBackTopOpenEdge        = Z(mesh); nodesBackTopOpenEdge       [mesh.nodesBackTopOpenEdge()       ] = 1
nodesBackLeftOpenEdge       = Z(mesh); nodesBackLeftOpenEdge      [mesh.nodesBackLeftOpenEdge()      ] = 1
nodesBackRightOpenEdge      = Z(mesh); nodesBackRightOpenEdge     [mesh.nodesBackRightOpenEdge()     ] = 1
nodesBottomLeftOpenEdge     = Z(mesh); nodesBottomLeftOpenEdge    [mesh.nodesBottomLeftOpenEdge()    ] = 1
nodesBottomRightOpenEdge    = Z(mesh); nodesBottomRightOpenEdge   [mesh.nodesBottomRightOpenEdge()   ] = 1
nodesTopLeftOpenEdge        = Z(mesh); nodesTopLeftOpenEdge       [mesh.nodesTopLeftOpenEdge()       ] = 1
nodesTopRightOpenEdge       = Z(mesh); nodesTopRightOpenEdge      [mesh.nodesTopRightOpenEdge()      ] = 1
nodesFrontBottomLeftCorner  = Z(mesh); nodesFrontBottomLeftCorner [mesh.nodesFrontBottomLeftCorner() ] = 1
nodesFrontBottomRightCorner = Z(mesh); nodesFrontBottomRightCorner[mesh.nodesFrontBottomRightCorner()] = 1
nodesFrontTopLeftCorner     = Z(mesh); nodesFrontTopLeftCorner    [mesh.nodesFrontTopLeftCorner()    ] = 1
nodesFrontTopRightCorner    = Z(mesh); nodesFrontTopRightCorner   [mesh.nodesFrontTopRightCorner()   ] = 1
nodesBackBottomLeftCorner   = Z(mesh); nodesBackBottomLeftCorner  [mesh.nodesBackBottomLeftCorner()  ] = 1
nodesBackBottomRightCorner  = Z(mesh); nodesBackBottomRightCorner [mesh.nodesBackBottomRightCorner() ] = 1
nodesBackTopLeftCorner      = Z(mesh); nodesBackTopLeftCorner     [mesh.nodesBackTopLeftCorner()     ] = 1
nodesBackTopRightCorner     = Z(mesh); nodesBackTopRightCorner    [mesh.nodesBackTopRightCorner()    ] = 1

# DOFs after eliminating periodicity
dofsPeriodic = mesh.dofsPeriodic()[:,0]

# open data file
f = h5py.File('MeshHex8-Regular-nodes_paraview.hdf5','w')

# write particle positions, and a dummy connectivity
f.create_dataset('/coor'                       ,data=mesh.coor()                )
f.create_dataset('/conn'                       ,data=np.arange(mesh.nnode())    )
f.create_dataset('/nodesFront'                 ,data=nodesFront                 )
f.create_dataset('/nodesBack'                  ,data=nodesBack                  )
f.create_dataset('/nodesLeft'                  ,data=nodesLeft                  )
f.create_dataset('/nodesRight'                 ,data=nodesRight                 )
f.create_dataset('/nodesBottom'                ,data=nodesBottom                )
f.create_dataset('/nodesTop'                   ,data=nodesTop                   )
f.create_dataset('/nodesFrontFace'             ,data=nodesFrontFace             )
f.create_dataset('/nodesBackFace'              ,data=nodesBackFace              )
f.create_dataset('/nodesLeftFace'              ,data=nodesLeftFace              )
f.create_dataset('/nodesRightFace'             ,data=nodesRightFace             )
f.create_dataset('/nodesBottomFace'            ,data=nodesBottomFace            )
f.create_dataset('/nodesTopFace'               ,data=nodesTopFace               )
f.create_dataset('/nodesFrontBottomEdge'       ,data=nodesFrontBottomEdge       )
f.create_dataset('/nodesFrontTopEdge'          ,data=nodesFrontTopEdge          )
f.create_dataset('/nodesFrontLeftEdge'         ,data=nodesFrontLeftEdge         )
f.create_dataset('/nodesFrontRightEdge'        ,data=nodesFrontRightEdge        )
f.create_dataset('/nodesBackBottomEdge'        ,data=nodesBackBottomEdge        )
f.create_dataset('/nodesBackTopEdge'           ,data=nodesBackTopEdge           )
f.create_dataset('/nodesBackLeftEdge'          ,data=nodesBackLeftEdge          )
f.create_dataset('/nodesBackRightEdge'         ,data=nodesBackRightEdge         )
f.create_dataset('/nodesBottomLeftEdge'        ,data=nodesBottomLeftEdge        )
f.create_dataset('/nodesBottomRightEdge'       ,data=nodesBottomRightEdge       )
f.create_dataset('/nodesTopLeftEdge'           ,data=nodesTopLeftEdge           )
f.create_dataset('/nodesTopRightEdge'          ,data=nodesTopRightEdge          )
f.create_dataset('/nodesFrontBottomOpenEdge'   ,data=nodesFrontBottomOpenEdge   )
f.create_dataset('/nodesFrontTopOpenEdge'      ,data=nodesFrontTopOpenEdge      )
f.create_dataset('/nodesFrontLeftOpenEdge'     ,data=nodesFrontLeftOpenEdge     )
f.create_dataset('/nodesFrontRightOpenEdge'    ,data=nodesFrontRightOpenEdge    )
f.create_dataset('/nodesBackBottomOpenEdge'    ,data=nodesBackBottomOpenEdge    )
f.create_dataset('/nodesBackTopOpenEdge'       ,data=nodesBackTopOpenEdge       )
f.create_dataset('/nodesBackLeftOpenEdge'      ,data=nodesBackLeftOpenEdge      )
f.create_dataset('/nodesBackRightOpenEdge'     ,data=nodesBackRightOpenEdge     )
f.create_dataset('/nodesBottomLeftOpenEdge'    ,data=nodesBottomLeftOpenEdge    )
f.create_dataset('/nodesBottomRightOpenEdge'   ,data=nodesBottomRightOpenEdge   )
f.create_dataset('/nodesTopLeftOpenEdge'       ,data=nodesTopLeftOpenEdge       )
f.create_dataset('/nodesTopRightOpenEdge'      ,data=nodesTopRightOpenEdge      )
f.create_dataset('/nodesFrontBottomLeftCorner' ,data=nodesFrontBottomLeftCorner )
f.create_dataset('/nodesFrontBottomRightCorner',data=nodesFrontBottomRightCorner)
f.create_dataset('/nodesFrontTopLeftCorner'    ,data=nodesFrontTopLeftCorner    )
f.create_dataset('/nodesFrontTopRightCorner'   ,data=nodesFrontTopRightCorner   )
f.create_dataset('/nodesBackBottomLeftCorner'  ,data=nodesBackBottomLeftCorner  )
f.create_dataset('/nodesBackBottomRightCorner' ,data=nodesBackBottomRightCorner )
f.create_dataset('/nodesBackTopLeftCorner'     ,data=nodesBackTopLeftCorner     )
f.create_dataset('/nodesBackTopRightCorner'    ,data=nodesBackTopRightCorner    )
f.create_dataset('/dofsPeriodic'               ,data=dofsPeriodic               )

# ======================================== write XDMF-file =========================================

xmf = '''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="Nodes">
      <Topology TopologyType="Polyvertex" NumberOfElements="{nnode:d}" NodesPerElement="1">
        <DataItem Dimensions="{nnode:d} 1" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/conn
        </DataItem>
      </Topology>
      <Geometry GeometryType="XYZ">
        <DataItem Dimensions="{nnode:d} 3" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/coor
        </DataItem>
      </Geometry>
      <Attribute Name="nodesFront" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFront
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBack" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBack
        </DataItem>
      </Attribute>
      <Attribute Name="nodesLeft" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesLeft
        </DataItem>
      </Attribute>
      <Attribute Name="nodesRight" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesRight
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottom" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBottom
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTop" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesTop
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontFace" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontFace
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackFace" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackFace
        </DataItem>
      </Attribute>
      <Attribute Name="nodesLeftFace" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesLeftFace
        </DataItem>
      </Attribute>
      <Attribute Name="nodesRightFace" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesRightFace
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottomFace" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBottomFace
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopFace" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesTopFace
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontBottomEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontBottomEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontTopEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontTopEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontLeftEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontLeftEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontRightEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontRightEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackBottomEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackBottomEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackTopEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackTopEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackLeftEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackLeftEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackRightEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackRightEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottomLeftEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBottomLeftEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottomRightEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBottomRightEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopLeftEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesTopLeftEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopRightEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesTopRightEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontBottomOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontBottomOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontTopOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontTopOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontLeftOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontLeftOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontRightOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontRightOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackBottomOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackBottomOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackTopOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackTopOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackLeftOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackLeftOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackRightOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackRightOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottomLeftOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBottomLeftOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBottomRightOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBottomRightOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopLeftOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesTopLeftOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesTopRightOpenEdge" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesTopRightOpenEdge
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontBottomLeftCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontBottomLeftCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontBottomRightCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontBottomRightCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontTopLeftCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontTopLeftCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesFrontTopRightCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesFrontTopRightCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackBottomLeftCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackBottomLeftCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackBottomRightCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackBottomRightCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackTopLeftCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackTopLeftCorner
        </DataItem>
      </Attribute>
      <Attribute Name="nodesBackTopRightCorner" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/nodesBackTopRightCorner
        </DataItem>
      </Attribute>
      <Attribute Name="dofsPeriodic" AttributeType="Scalar" Center="Node">
         <DataItem Dimensions="{nnode:d}" NumberType="Int" Format="HDF">
          MeshHex8-Regular-nodes_paraview.hdf5:/dofsPeriodic
         </DataItem>
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
'''

open('MeshHex8-Regular-nodes_paraview.xmf','w').write(xmf.format(
  nnode = mesh.nnode(),
))

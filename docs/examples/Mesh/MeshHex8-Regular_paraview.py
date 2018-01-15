
# ========================================= load libraries =========================================

import GooseFEM as gf
import numpy as np
import h5py

# ========================================== define mesh ===========================================

mesh = gf.Mesh.Hex8.Regular(3,4,5)
coor = mesh.coor()
conn = mesh.conn()

# ======================================== write HDF5=file =========================================

# set file-name base
name = 'MeshHex8-Regular_paraview'

# open data file
f = h5py.File('{name:s}.hdf5'.format(name=name),'w')

# write particle positions, and a dummy connectivity
f.create_dataset('/coor',data=coor)
f.create_dataset('/conn',data=conn)

# create a sample deformation: simple shear
for inc,gamma in enumerate(np.linspace(0,1,100)):

  # - initialize displacement (should be 3-d for ParaView)
  U = np.zeros((mesh.nnode(),3),dtype='float64')
  # - set
  U[:,0] += gamma * coor[:,1]
  # - store
  f.create_dataset('/disp/%d'%inc,data=U)

# ======================================== write XDMF-file =========================================

# --------------------------------- format of the main structure ----------------------------------

xmf = '''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
{series:s}
    </Grid>
  </Domain>
</Xdmf>
'''

# ----------------------------- format of an increment in time-series ------------------------------

grid = '''<Grid Name="Increment = {inc:d}">
  <Time Value="{inc:d}"/>
  <Topology TopologyType="Hexahedron" NumberOfElements="{nelem:d}">
    <DataItem Dimensions="{nelem:d} {nne:d}" Format="HDF">
    {name:s}.hdf5:/conn
    </DataItem>
  </Topology>
  <Geometry GeometryType="XYZ">
    <DataItem Dimensions="{nnode:d} {ndim:d}" Format="HDF">
    {name:s}.hdf5:/coor
    </DataItem>
  </Geometry>
  <Attribute Name="Displacement" AttributeType="Vector" Center="Node">
     <DataItem Dimensions="{nnode:d} 3" NumberType="Float" Precision="8" Format="HDF">
      {name:s}.hdf5:/disp/{inc:d}
     </DataItem>
  </Attribute>
</Grid>
'''

# ------------------------------------------- write file -------------------------------------------

# initialize string that will contain the full time series
txt = ''

# loop over all increments, append the time series
for inc in range(100):
  txt += grid.format(
    name  = name,
    inc   = inc,
    nnode = mesh.nnode(),
    ndim  = mesh.ndim(),
    nelem = mesh.nelem(),
    nne   = mesh.nne()
  )

# write xmf-file, fix the indentation
fname = '{name:s}.xmf'.format(name=name)
open(fname,'w').write(xmf.format(series='      '+txt.replace('\n','\n      ')))

# ======================================== print node=sets =========================================

np.set_printoptions(linewidth=200)

print('nodesBottom                 =', mesh.nodesBottom())
print('nodesBottomFace             =', mesh.nodesBottomFace())
print('nodesBottomFrontEdge        =', mesh.nodesBottomFrontEdge())
print('nodesBottomBackEdge         =', mesh.nodesBottomBackEdge ())
print('nodesBottomLeftEdge         =', mesh.nodesBottomLeftEdge ())
print('nodesBottomRightEdge        =', mesh.nodesBottomRightEdge())
print('nodesTop                    =', mesh.nodesTop())
print('nodesTopFace                =', mesh.nodesTopFace())
print('nodesTopFrontEdge           =', mesh.nodesTopFrontEdge())
print('nodesTopBackEdge            =', mesh.nodesTopBackEdge ())
print('nodesTopLeftEdge            =', mesh.nodesTopLeftEdge ())
print('nodesTopRightEdge           =', mesh.nodesTopRightEdge())
print('nodesFront                  =', mesh.nodesFront())
print('nodesFrontFace              =', mesh.nodesFrontFace())
print('nodesFrontLeftEdge          =', mesh.nodesFrontLeftEdge())
print('nodesFrontRightEdge         =', mesh.nodesFrontRightEdge())
print('nodesFrontBottomEdge        =', mesh.nodesFrontBottomEdge())
print('nodesFrontTopEdge           =', mesh.nodesFrontTopEdge())
print('nodesBack                   =', mesh.nodesBack())
print('nodesBackFace               =', mesh.nodesBackFace())
print('nodesBackLeftEdge           =', mesh.nodesBackLeftEdge())
print('nodesBackRightEdge          =', mesh.nodesBackRightEdge())
print('nodesBackBottomEdge         =', mesh.nodesBackBottomEdge())
print('nodesBackTopEdge            =', mesh.nodesBackTopEdge())
print('nodesLeft                   =', mesh.nodesLeft())
print('nodesLeftFace               =', mesh.nodesLeftFace())
print('nodesLeftBottomEdge         =', mesh.nodesLeftBottomEdge())
print('nodesLeftTopEdge            =', mesh.nodesLeftTopEdge())
print('nodesLeftFrontEdge          =', mesh.nodesLeftFrontEdge())
print('nodesLeftBackEdge           =', mesh.nodesLeftBackEdge())
print('nodesRight                  =', mesh.nodesRight())
print('nodesRightFace              =', mesh.nodesRightFace())
print('nodesRightBottomEdge        =', mesh.nodesRightBottomEdge())
print('nodesRightTopEdge           =', mesh.nodesRightTopEdge())
print('nodesRightFrontEdge         =', mesh.nodesRightFrontEdge())
print('nodesRightBackEdge          =', mesh.nodesRightBackEdge())
print('nodesLeftBottomFrontCorner  =', mesh.nodesLeftBottomFrontCorner())
print('nodesRightBottomFrontCorner =', mesh.nodesRightBottomFrontCorner())
print('nodesRightBottomBackCorner  =', mesh.nodesRightBottomBackCorner())
print('nodesLeftBottomBackCorner   =', mesh.nodesLeftBottomBackCorner())
print('nodesLeftTopFrontCorner     =', mesh.nodesLeftTopFrontCorner())
print('nodesRightTopFrontCorner    =', mesh.nodesRightTopFrontCorner())
print('nodesRightTopBackCorner     =', mesh.nodesRightTopBackCorner())
print('nodesLeftTopBackCorner      =', mesh.nodesLeftTopBackCorner())

print('nodesPeriodic =')
print(mesh.nodesPeriodic())
print('nodesOrigin =')
print(mesh.nodesOrigin())

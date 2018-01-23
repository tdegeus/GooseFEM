
import numpy as np
import h5py
import GooseFEM as gf

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Quad4.FineLayer(9,17)

# open data file
f = h5py.File('MeshQuad4-FineLayer_paraview.hdf5','w')

# element set
elementsMiddleLayer = np.zeros((mesh.nelem()),dtype='int')
elementsMiddleLayer[mesh.elementsMiddleLayer()] = 1

# write nodal coordinates and connectivity
f.create_dataset('/coor'               ,data=mesh.coor()        )
f.create_dataset('/conn'               ,data=mesh.conn()        )
f.create_dataset('/elementsMiddleLayer',data=elementsMiddleLayer)

# ======================================== write XDMF-file =========================================

# basic file-format
xmf = '''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="Mesh">
      <Topology TopologyType="Quadrilateral" NumberOfElements="{nelem:d}">
        <DataItem Dimensions="{nelem:d} {nne:d}" Format="HDF">
          MeshQuad4-FineLayer_paraview.hdf5:/conn
        </DataItem>
      </Topology>
      <Geometry GeometryType="XY">
        <DataItem Dimensions="{nnode:d} 2" Format="HDF">
          MeshQuad4-FineLayer_paraview.hdf5:/coor
        </DataItem>
      </Geometry>
      <Attribute Name="elementsMiddleLayer" AttributeType="Scalar" Center="Cell">
        <DataItem Dimensions="{nelem:d}" NumberType="Int" Format="HDF">
          MeshQuad4-FineLayer_paraview.hdf5:/elementsMiddleLayer
        </DataItem>
      </Attribute>
      </Grid>
  </Domain>
</Xdmf>
'''

# write to file, fill mesh dimensions
open('MeshQuad4-FineLayer_paraview.xmf','w').write(xmf.format(
  nnode = mesh.nnode(),
  nelem = mesh.nelem(),
  nne   = mesh.nne(),
))

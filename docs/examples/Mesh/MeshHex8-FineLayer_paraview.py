
import numpy as np
import h5py
import GooseFEM as gf

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Hex8.FineLayer(9,17,27)

# open data file
f = h5py.File('MeshHex8-FineLayer_paraview.hdf5','w')

# element set
elementsMiddleLayer = np.zeros((mesh.nelem()),dtype='int')
elementsMiddleLayer[mesh.elementsMiddleLayer()] = 1

# write nodal coordinates and connectivity
f.create_dataset('/coor'               ,data=mesh.coor()        )
f.create_dataset('/conn'               ,data=mesh.conn()        )
f.create_dataset('/elementsMiddleLayer',data=elementsMiddleLayer)

# ======================================== write XDMF-file =========================================

# basic file-format
xdmf = '''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="Mesh">
      <Topology TopologyType="Hexahedron" NumberOfElements="{nelem:d}">
        <DataItem Dimensions="{nelem:d} {nne:d}" Format="HDF">
          MeshHex8-FineLayer_paraview.hdf5:/conn
        </DataItem>
      </Topology>
      <Geometry GeometryType="XYZ">
        <DataItem Dimensions="{nnode:d} 3" Format="HDF">
          MeshHex8-FineLayer_paraview.hdf5:/coor
        </DataItem>
      </Geometry>
      <Attribute Name="elementsMiddleLayer" AttributeType="Scalar" Center="Cell">
        <DataItem Dimensions="{nelem:d}" NumberType="Int" Format="HDF">
          MeshHex8-FineLayer_paraview.hdf5:/elementsMiddleLayer
        </DataItem>
      </Attribute>
      </Grid>
  </Domain>
</Xdmf>
'''

# write to file, fill mesh dimensions
open('MeshHex8-FineLayer_paraview.xdmf','w').write(xdmf.format(
  nnode = mesh.nnode(),
  nelem = mesh.nelem(),
  nne   = mesh.nne(),
))

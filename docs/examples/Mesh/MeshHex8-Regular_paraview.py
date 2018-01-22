
import numpy as np
import h5py
import GooseFEM as gf

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Hex8.Regular(6,8,12)

# open data file
f = h5py.File('MeshHex8-Regular_paraview.hdf5','w')

# write particle positions, and a dummy connectivity
f.create_dataset('/coor',data=mesh.coor())
f.create_dataset('/conn',data=mesh.conn())

# ======================================== write XDMF-file =========================================

xmf = '''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="Mesh">
      <Topology TopologyType="Hexahedron" NumberOfElements="{nelem:d}">
        <DataItem Dimensions="{nelem:d} {nne:d}" Format="HDF">
          MeshHex8-Regular_paraview.hdf5:/conn
        </DataItem>
      </Topology>
      <Geometry GeometryType="XYZ">
        <DataItem Dimensions="{nnode:d} 3" Format="HDF">
          MeshHex8-Regular_paraview.hdf5:/coor
        </DataItem>
      </Geometry>
      </Grid>
  </Domain>
</Xdmf>
'''

open('MeshHex8-Regular_paraview.xmf','w').write(xmf.format(
  nnode = mesh.nnode(),
  nelem = mesh.nelem(),
  nne   = mesh.nne(),
))

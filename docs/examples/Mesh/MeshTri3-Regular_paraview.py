
import numpy as np
import h5py
import GooseFEM as gf

# ====================== create fictitious configuration + store to HDF5-file ======================

# create mesh object
mesh = gf.Mesh.Tri3.Regular(9,11)

# open data file
f = h5py.File('MeshTri3-Regular_paraview.hdf5','w')

# write nodal coordinates and connectivity
f.create_dataset('/coor',data=mesh.coor())
f.create_dataset('/conn',data=mesh.conn())

# ======================================== write XDMF-file =========================================

# basic file-format
xmf = '''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
  <Domain>
    <Grid Name="Mesh">
      <Topology TopologyType="Triangle" NumberOfElements="{nelem:d}">
        <DataItem Dimensions="{nelem:d} {nne:d}" Format="HDF">
          MeshTri3-Regular_paraview.hdf5:/conn
        </DataItem>
      </Topology>
      <Geometry GeometryType="XY">
        <DataItem Dimensions="{nnode:d} 2" Format="HDF">
          MeshTri3-Regular_paraview.hdf5:/coor
        </DataItem>
      </Geometry>
      </Grid>
  </Domain>
</Xdmf>
'''

# write to file, fill mesh dimensions
open('MeshTri3-Regular_paraview.xmf','w').write(xmf.format(
  nnode = mesh.nnode(),
  nelem = mesh.nelem(),
  nne   = mesh.nne(),
))
